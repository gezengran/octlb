# OctLB 框架 PRD

> 版本：v0.1  
> 状态：草稿  
> 日期：2026-05-28

---

## 问题陈述

现有的大规模 CFD 并行求解框架面临两个核心矛盾：

1. **精度 vs 计算量**：均匀结构网格（如 OpenLB 的 CuboidDecomposition）在复杂几何体附近要达到足够分辨率，需要对全域加密，计算量激增；而自适应网格（AMR）可以将计算集中在需要的区域，但现有 LBM 框架缺乏健壮的 AMR 支持。

2. **物理丰富度 vs 拓扑灵活性**：OpenLB 拥有最完整的 LBM 物理算子库（BGK/MRT/LES/Smagorinsky/Bouzidi BC 等），但其拓扑层深度绑定均匀结构网格，无法直接承载八叉树 AMR；p4est 是工业级八叉树并行框架，提供精确的拓扑管理和 2:1 平衡保证，但自身不含任何物理求解能力。

前期的 octree-mesh (SAMR) 已验证了 p4est + 原生 LBM 的可行性，但存在以下关键缺口：邻居发现未使用 `p4est_iterate` face callback（O(n²) 空间搜索）、AMR 粗细界面插值不守恒、跨 rank 字段迁移仅标记 stale 未实际搬运数据、ghost exchange 重复造轮子（tag 7001/7002 自定义 MPI）。

**OctLB** 旨在通过将 p4est 的拓扑管理与 OpenLB 的物理算子深度融合，构建一个**正确、高效、可扩展的大规模并行 LBM-AMR 求解框架**。

---

## 解决方案

OctLB 采用**两模块架构**：

- **Mesh 模块**：以 p4est 为核心，负责八叉树拓扑、几何体素化、负载均衡分区，向 Solver 模块暴露纯拓扑接口（`FacePairList`、`BlockRegistry`、`MaterialField`），对物理量完全无感知。

- **Solver 模块**：内部分为两个子层——
  - `field/`：泛型字段容器与迭代器（`BlockCollection<T>`、`GhostSchedule<T>`、`FaceIterator`），未来可提取为独立的第三层供其他物理模型复用；
  - `lbm/`：从 OpenLB 按需复制并适配的 LBM 物理算子（`ConcreteBlockLattice`、Dynamics、Lagrava 耦合、边界条件）。

两个模块通过清晰定义的接口解耦，Solver 不直接调用任何 p4est API。

---

## 用户故事

1. 作为 CFD 仿真工程师，我希望用 STL 文件描述几何体并自动生成多层 AMR 网格，以便无需手工网格划分即可对复杂几何体进行高精度模拟。

2. 作为 CFD 仿真工程师，我希望框架在粗细网格界面处使用 Lagrava 守恒插值方案，以便 AMR 结果与均匀网格参考解的误差可控。

3. 作为 HPC 用户，我希望在 1000+ MPI rank 环境下运行仿真，每步的通信量与本地网格面积成正比而非体积，以便通信开销不成为瓶颈。

4. 作为 HPC 用户，我希望各 rank 的计算量自动均衡（细层块获得更高权重），以便不因 AMR 层级差异产生严重的负载不均。

5. 作为求解器开发者，我希望调用 OpenLB 中已有的 BGK/MRT/Smagorinsky 等 Dynamics 算子，以便不需要重新实现和验证 LBM 物理模型。

6. 作为求解器开发者，我希望 Mesh 模块完全不感知 LBM 物理量，以便将来在同一 Mesh 上叠加 FVM 或其他 PDE 求解器。

7. 作为框架维护者，我希望 `field/` 子层中的 `BlockCollection<T>` 和 `GhostSchedule<T>` 不依赖任何 LBM 头文件，以便将来将其提取为独立的字段容器中间层。

8. 作为仿真工程师，我希望以 VTK `.vts`/`.vtm`/`.pvd` 格式输出并行 AMR 仿真结果，以便在 ParaView 中直接可视化不同细化层级的流场。

9. 作为验证工程师，我希望用 cylinder3d 算例对比 OctLB 与 OpenLB 参考解的阻力系数（目标：相对误差 < 1%），以便确认 AMR 求解的正确性。

---

## 实现决策

### 架构总览

```
OctLB/
├── doc/                  # README; prd/, tasks/, dev/
├── src/
│   ├── mesh/             ← Mesh 模块（不 include OpenLB，不含物理量）
│   │   ├── forest/       # OctreeForest：p4est refine/balance/partition 封装
│   │   ├── topology/     # BlockRegistry, FacePairList
│   │   ├── geometry/     # GeometryEngine（CGAL 体素化）、MaterialField 产生
│   │   ├── io/           # Mesh 模块内共用 IO（不依赖 OpenLB）
│   │   │   └── stl_reader/   # 移植自 OpenLB，唯一消费方是 GeometryEngine
│   │   └── load_balance/ # WeightedLoadBalancer
│   │
│   └── solver/
│       ├── field/        ← FieldContainer 层（不 include LBM 头文件）
│       │   ├── block_collection.h
│       │   ├── block_iterator.h
│       │   ├── face_iterator.h
│       │   └── ghost_schedule.h
│       │
│       ├── io/           ← Solver 模块内共用 IO（不依赖 LBM 头文件）
│       │   ├── vtk_writer/   # 移植 blockVtkWriter3D，替换驱动层
│       │   └── gnuplot/      # 直接复制自 OpenLB
│       │
│       └── lbm/          ← LBM 专属（从 OpenLB 复制并适配）
│           ├── core/         # ConcreteBlockLattice, Cell, stages
│           ├── descriptor/   # D3Q19 等格子描述符
│           ├── dynamics/     # BGK, MRT, Smagorinsky 等
│           ├── boundary/     # Bouzidi 曲面 BC, 其他 BC
│           ├── refinement/   # Lagrava 算子（移植自 OpenLB src/refinement/）
│           ├── time_loop/    # 递归下降 TimeLoop
│           └── unit_converter/
│
├── examples/
│   ├── cavity3d/         # 初期验证：无 STL，无 AMR，Lid-driven cavity
│   ├── cylinder3d/       # 主验证：AMR + STL，对比 OpenLB/octree-mesh Cd
│   └── sphere3d/         # 大规模并行基准测试
│
└── tests/
    ├── unit/
    │   ├── mesh/
    │   └── solver/
    └── integration/
```

---

### 模块一：Mesh 模块

#### OctreeForest

- 封装 p4est 的 `p8est_t` forest，提供 refine / coarsen / balance / partition 操作。
- 对外暴露接口：`local_octants()`, `ghost_octants()`, `quadrant_bounds(octant_id)`, `quadrant_level(octant_id)`。
- 不暴露 p4est 内部类型给 Solver 模块。

#### FacePairList

- 通过 `p8est_iterate` + face callback 一次性遍历所有面，构建两类列表：
  - **SameLevelFaces**：`{local_block_id, face_dir, remote_block_id, remote_rank}`
  - **CoarseFineFaces**：`{coarse_block_id, fine_block_ids[4], face_normal, remote_ranks[4]}`
- 粗细面识别依据：face callback 中 `is_hanging == 1` 的一侧为 4 个细格，`is_hanging == 0` 的一侧为 1 个粗格。
- 静态分层下初始化一次，预留 `rebuild()` 接口供动态 AMR 扩展。

#### GeometryEngine / MaterialField

- 读入 STL，用 CGAL 体素化产生每个 octant 块内 N³ 格的 material 编号（fluid / solid / boundary）。
- 产生结果存储在 `MaterialField`（Mesh 模块持有），Solver 模块在初始化时**单次读取**，之后不再持有引用。

#### WeightedLoadBalancer

- 提供 `p8est_partition_ext` 的权重回调：weight(octant) = 2^level。
- 在 `OctreeForest::balance()` 完成后调用。
- 第一版静态分层只需初始化时调用一次。

---

### 模块二：Solver 模块

#### 核心索引约定

所有组件使用统一的格子地址 `(octant_id, i, j, k)`，无需翻译表：
- `octant_id` → `BlockStore` 中的 `ConcreteBlockLattice`
- `(i, j, k)` → 块内坐标，直接映射到 OpenLB `BlockLattice::cell(i, j, k)`

这是与 octree-mesh 中 `mb_to_mesh_` / `mesh_to_mb_` 翻译表的根本性改进。

#### field/ 子层（泛型，不依赖 LBM）

- **`BlockCollection<T>`**：管理本 rank 持有的所有块，每块存一个 `T` 类型对象，以 `octant_id` 为键。
- **`GhostSchedule<T>`**：从 `FacePairList::SameLevelFaces` 构建 MPI 通信计划。
  - 交换单位：N×N 面层（`N²×q` 个值），而非整块（N³×q）。
  - 与 N=8 相比，通信量是整块交换的 1/8；N=10 时为 1/10。
  - 提供 `pack_face(block, face_dir, buffer)` / `unpack_face(buffer, block, face_dir)` 回调接口，对数据格式无假设。
- **`BlockIterator`**：遍历本 rank 所有 `octant_id`。
- **`FaceIterator`**：遍历 `FacePairList` 中的面对，提供 `(coarse_block, fine_blocks[4], normal)` 视图。

#### lbm/ 子层（LBM 专属）

**BlockStore（LBM 具体化）**

- 每个 octant 持有一个 `ConcreteBlockLattice<T, D3Q19Descriptor, Platform::CPU_SISD>`（或 GPU 变体）。
- 块尺寸：`(N+2*overlap)³`，overlap = 1（D3Q19 streaming 只需 1 格 overlap）。
- 初始化时读取 `MaterialField`，为每格设置对应 `Dynamics`（fluid → BGK/MRT，solid → BounceBack，boundary → Bouzidi）。

**HaloExchange**

- 使用 `GhostSchedule<f_distribution>` 在每次 streaming 后交换 6 个方向的面层。
- 面层数据格式：19 个 double/float 的 N×N slab，对应 f 分布函数分量。

**LevelCoupler（Lagrava 方案）**

移植自 OpenLB `src/refinement/` 的 Lagrava 算子，适配到 `CoarseFineFaces` 拓扑：

- **耦合面来源**：`FacePairList::CoarseFineFaces`（由 `p8est_iterate` face callback 产生），替代 OpenLB 原版的 `SuperIndicatorDomainFrontierDistanceF`。
- **粗→细（prolongation）**：
  - `scalingFactor = (τ_coarse - 0.25) / τ_coarse`
  - `f_fine = f_eq(ρ, u) + scalingFactor × f_neq`
  - 在细层第一个子步前（half-time）和第二个子步前（full-time）各执行一次。
- **细→粗（restriction）**：
  - `scalingFactor = τ_coarse / (τ_coarse - 0.25)`
  - 对 4 个细格取均值后重建宏观量，再写回粗格分布函数。
  - 在两个细子步完成后执行。
- 跨 rank 的粗细耦合通过 `FaceIterator` 触发额外的 MPI 交换（coarse macro fields → fine side）。

**TimeLoop（递归下降）**

```
advance(level l):
  collide_and_stream(all blocks at level l)  // OpenLB ConcreteBlockLattice::collideAndStream()
  halo_exchange(level l)                      // GhostSchedule，面层 MPI
  coupler.apply_half_time(l → l+1)            // LevelCoupler prolongation half-time
  advance(level l+1)                          // 第一个细子步（递归）
  coupler.apply_full_time(l → l+1)            // LevelCoupler prolongation full-time
  advance(level l+1)                          // 第二个细子步（递归）
  coupler.restrict(l+1 → l)                   // LevelCoupler restriction
```

- 递归深度 = 细化层数（6-8，最大 10），不会溢出栈。
- 每层的 `collide_and_stream` 是 OpenLB `BlockLattice::collideAndStream()` 的并行循环。
- 第一版：静态分层，`advance(level 0)` 为入口，层级结构在初始化时确定。

**VTK Writer（IO）**

- 驱动层：`BlockStore` 中所有 `octant_id` 的遍历，替代 OpenLB 的 `CuboidDecomposition` 循环。
- 每个块输出一个 `.vts` 文件，空间范围和 `deltaX` 从 `OctreeForest::quadrant_bounds(octant_id)` 与 N 计算得到。
- Rank 0 汇总写 `.vtm`（MultiBlock 索引）和 `.pvd`（时间序列索引）。
- 支持输出字段：速度（u）、压力（p）、密度（ρ），可扩展。

---

### 负载均衡策略

- 权重：weight(octant) = 2^level（细层权重高）
- 接口：`p8est_partition_ext` 权重回调
- 父子 octant 允许跨 rank（不施加亲和性约束）
- 粗细耦合的跨 rank 额外通信量：N×N 宏观量（ρ, u, f_neq），数据量小于 halo 面层

---

### OpenLB 代码融入策略

- 从 OpenLB `src/` 按需复制并适配到以下位置：
  - `src/solver/lbm/core/`（BlockLattice、Cell、stages、platform）
  - `src/solver/lbm/descriptor/`（D3Q19、字段元描述）
  - `src/solver/lbm/dynamics/`（BGK、MRT、Smagorinsky 等）
  - `src/solver/lbm/boundary/`（Bouzidi、速度/压力 BC）
  - `src/solver/lbm/refinement/`（Lagrava/Rohde 算子，作为移植基础）
  - `src/solver/lbm/unit_converter/`（unitConverter 等）
  - `src/mesh/io/stl_reader/`（stlReader；不依赖 OpenLB 其他模块，唯一消费方是 GeometryEngine）
  - `src/solver/io/vtk_writer/`（blockVtkWriter3D；替换驱动层，不依赖 LBM 头文件）
  - `src/solver/io/gnuplot/`（gnuplotWriter；直接复制）

**IO 模块归属原则**：STL 读入是几何输入，属 Mesh 模块关切；VTK/gnuplot 是仿真输出，属 Solver 模块关切。两者均不依赖对方模块头文件，模块内共用，不在 lbm/ 子目录内重复放置。
- 不保留 OpenLB 的 build system（Makefile/config.mk）
- 版本基准：OpenLB 1.9.0
- 后续 OpenLB 升级：手动 diff 相关子目录，按需合并

---

## 测试决策

**好测试的定义**：只测对外行为（接口契约），不测实现细节（内部数据结构、MPI tag 值等）。

**开发策略：严格 TDD**——测试顺序即开发顺序，每个组件先写测试（红），再实现（绿），再重构。每条测试只因一个原因失败。

**单元测试 Fixture 原则**：凡依赖文件 I/O 的组件（如 GeometryEngine、VTK Writer），单元测试一律使用硬编码 fixture 数据（triangle soup、BlockCollection），不依赖真实文件。真实文件路径仅在集成测试中出现。

**测试顺序（= 开发顺序，各层独立可测）**：

| 顺序 | 测试名 | 验证内容 | 层级 |
|---|---|---|---|
| 0 | `test_octree_forest` | OctreeForest refine/balance/partition 在单位立方体上的基本正确性；验证 `local_octants()`、`quadrant_bounds()`、`quadrant_level()` | Mesh/forest |
| 1 | `test_stl_reader` | ASCII/binary STL 解析：三角面片数、法向量方向、包围盒正确性 | Mesh/io |
| 2 | `test_geometry_engine` | 硬编码 triangle soup（立方体）输入，CGAL 体素化产出正确 fluid/solid/boundary 标签；不依赖文件 I/O | Mesh/geometry |
| 3 | `test_face_pair_list` | p8est_iterate face callback 正确识别同级面和粗细 hanging 面 | Mesh/topology |
| 4 | `test_load_balancer` | 2^level 权重分区后各 rank 总权重差异 < 5% | Mesh |
| 5 | `test_ghost_topology` | ghost octant 列表与 FacePairList 一致（每条面对两端都认识对方） | Mesh |
| 6 | `test_block_collection` | `BlockCollection<T>` 增删查，以 `octant_id` 为键；`BlockIterator` 遍历正确 | Solver/field |
| 7 | `test_collision_bgk` | 单块 BGK 碰撞：质量守恒、动量守恒，收敛至 Maxwell 平衡态 | Solver/lbm |
| 8 | `test_ghost_schedule` | 2-rank，交换面层后两块 overlap 值一致（MPI 正确性） | Solver/field |
| 9 | `test_lagrava_coupler` | 2-block（1 粗 1 细），经 N 步后粗细界面处密度/速度连续，无虚假反射 | Solver/lbm |
| 10 | `test_time_loop_levels` | L=3 递归时间步，各层步数计数符合 1:2:4 关系 | Solver/lbm |
| 11 | `test_vtk_writer` | 硬编码 `BlockCollection` 输入，输出合法 VTK XML，空间范围与 `quadrant_bounds` 吻合 | Solver/io |
| 12 | `test_cavity3d_serial` | 单 rank，cavity3d，Re=100，中线速度剖面与 Ghia 1982 参考值误差 < 2% | Integration |
| 13 | `test_cylinder3d_parallel` | 多 rank（≥4），cylinder3d L=4，阻力系数 Cd 与 OpenLB 参考值误差 < 1% | Integration |
| 14 | `test_amr_convergence` | L=1→3 逐级细化，速度场 L2 误差收敛阶 ≈ 2 | Integration |

**依赖链**：0 → 3 → 5 → 8 → 9 → 10；1 → 2；6 → 8；7 → 9；0–11 全绿后进入 Integration。

**测试基础设施**：GTest + MPI（参考 octree-mesh 现有 `tests/` 结构）。

---

## 不在范围内（第一版）

- **动态 AMR**：运行时 refine/coarsen 及字段迁移（`p4est_transfer_fixed`）。预留接口，不实现。
- **GPU backend**：OpenLB 的 `Platform::GPU_CUDA` / `GPU_HIP`。预留 `Platform` 模板参数，不测试。
- **HDF5 checkpoint/restart**：预留 `io/checkpoint/` 目录，不实现。
- **FVM 或其他物理模型**：`field/` 子层设计为可复用，但第一版只跑 LBM。
- **自适应细化判据**：如 Q-criterion、速度梯度触发。静态分层由几何体素化决定。
- **多相流、热 LBM 等扩展物理模型**：第一版只验证单相不可压 NS。

---

## 其他说明

### 与 octree-mesh (SAMR) 的关系

OctLB 是对 octree-mesh 的**重新设计**，不是直接扩展。以下 octree-mesh 组件可作为参考或局部复用：

| octree-mesh 组件 | OctLB 对应 | 处理方式 |
|---|---|---|
| `OctreeMesh` + p4est 封装 | `OctreeForest` | 参考重写，正确使用 `p8est_iterate` face callback |
| `GeometryAdaptiveEngine` + CGAL | `GeometryEngine` | 参考复用 CGAL 体素化逻辑 |
| `MBArray` / `BlockField` | `BlockCollection<T>` + OpenLB `BlockLattice` | 废弃，换 OpenLB 存储 |
| `LbmSolver::ExchangeGhostFDistributions` | `GhostSchedule<f_distribution>` | 废弃，重写为面层 MPI |
| `AMRInterpolator` | `LevelCoupler`（Lagrava） | 废弃，换守恒插值 |
| `MigrationEngine` | 第一版不实现 | 预留接口 |
| Bouzidi 曲面 BC 逻辑 | OpenLB `boundary/` | 参考，以 OpenLB 为准 |
| `compare_cylinder3d.sh` | `test_cylinder3d_parallel` | 测试用例中复用数据文件和比对脚本 |

### p4est API 使用规范

OctLB 应完整使用 p4est 提供的通信接口，不重造轮子：
- ghost 拓扑：`p4est_ghost_new` / `p4est_ghost_destroy`
- 邻居遍历：`p8est_iterate` + face callback（含 `is_hanging` 检测）
- 重分区迁移：`p8est_transfer_fixed`（动态 AMR 阶段启用）

不应自行实现等价的 MPI 通信逻辑（octree-mesh 中 tag 7001/7002 的教训）。

---

## 需求完成情况

> 状态说明：`未开始` / `开发中` / `测试绿` / `完成`

### Mesh 模块

| 组件 | 所在路径 | 状态 | 对应测试 | 备注 |
|---|---|---|---|---|
| OctreeForest | `mesh/forest/` | 未开始 | `test_octree_forest` | 参考 octree-mesh 重写 |
| stl_reader | `mesh/io/stl_reader/` | 未开始 | `test_stl_reader` | 移植自 OpenLB |
| GeometryEngine + MaterialField | `mesh/geometry/` | 未开始 | `test_geometry_engine` | CGAL 体素化；参考 GeometryAdaptiveEngine |
| FacePairList | `mesh/topology/` | 未开始 | `test_face_pair_list` | 基于 `p8est_iterate` face callback |
| WeightedLoadBalancer | `mesh/load_balance/` | 未开始 | `test_load_balancer` | 权重 2^level |
| Ghost 拓扑一致性 | `mesh/topology/` | 未开始 | `test_ghost_topology` | 验证两端互认 |

### Solver/field 模块

| 组件 | 所在路径 | 状态 | 对应测试 | 备注 |
|---|---|---|---|---|
| BlockCollection\<T\> + BlockIterator | `solver/field/` | 未开始 | `test_block_collection` | 泛型，不依赖 LBM 头文件 |
| GhostSchedule\<T\> + HaloExchange | `solver/field/` | 未开始 | `test_ghost_schedule` | 面层 MPI，替代 tag 7001/7002 |

### Solver/lbm 模块

| 组件 | 所在路径 | 状态 | 对应测试 | 备注 |
|---|---|---|---|---|
| ConcreteBlockLattice + BGK | `solver/lbm/core/` + `dynamics/` | 未开始 | `test_collision_bgk` | 移植自 OpenLB |
| MRT / Smagorinsky Dynamics | `solver/lbm/dynamics/` | 未开始 | — | 移植自 OpenLB，集成测试间接覆盖 |
| Bouzidi BC | `solver/lbm/boundary/` | 未开始 | — | 移植自 OpenLB |
| LevelCoupler（Lagrava） | `solver/lbm/refinement/` | 未开始 | `test_lagrava_coupler` | 移植自 OpenLB src/refinement/ |
| TimeLoop（递归下降） | `solver/lbm/time_loop/` | 未开始 | `test_time_loop_levels` | |
| unit_converter | `solver/lbm/unit_converter/` | 未开始 | — | 移植自 OpenLB，LBM 专属 |

### Solver/io 模块

| 组件 | 所在路径 | 状态 | 对应测试 | 备注 |
|---|---|---|---|---|
| vtk_writer | `solver/io/vtk_writer/` | 未开始 | `test_vtk_writer` | 移植 blockVtkWriter3D，替换驱动层 |
| gnuplot | `solver/io/gnuplot/` | 未开始 | — | 直接复制自 OpenLB |

### 集成测试

| 测试 | 状态 | 前置条件 | 验收标准 |
|---|---|---|---|
| `test_cavity3d_serial` | 未开始 | 测试 0–11 全绿 | Re=100 中线速度剖面与 Ghia 1982 误差 < 2% |
| `test_cylinder3d_parallel` | 未开始 | `test_cavity3d_serial` 通过 | ≥4 rank，Cd 与 OpenLB 参考值误差 < 1% |
| `test_amr_convergence` | 未开始 | `test_cylinder3d_parallel` 通过 | L=1→3 速度场 L2 误差收敛阶 ≈ 2 |
