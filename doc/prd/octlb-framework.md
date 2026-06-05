# OctLB 框架 PRD

> 版本：v0.1  
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

7. 作为框架维护者，我希望 `field/` 子层中的 `BlockCollection<T>`、`GhostSchedule<T>`（及 `FacePackable`）不依赖任何 LBM 头文件；`GhostSchedule` 可通过 `octlb_field_schedule` 依赖 Mesh 拓扑，以便将来将无 Mesh 部分提取为独立字段容器中间层。

8. 作为仿真工程师，我希望以 VTK `.vts`/`.vtm`/`.pvd` 格式输出并行 AMR 仿真结果，以便在 ParaView 中直接可视化不同细化层级的流场。

9. 作为验证工程师，我希望用 cylinder3d 算例对比 OctLB 与 OpenLB 参考解的阻力系数（目标：相对误差 < 1%），以便确认 AMR 求解的正确性。

---

## 实现决策

### 架构总览

```
OctLB/
├── cmake/                # Find 模块：FindP4est.cmake 等
├── doc/                  # README; prd/, tasks/, dev/
├── src/
│   ├── common/           # 全项目共用基础类型（不依赖 p4est / MPI / OpenLB）
│   │   ├── types.h           # label, scalar, OctantId（OpenFOAM 风格）
│   │   └── bounding_box.h    # BoundingBox struct
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
│       │   ├── face_packable.h
│       │   ├── ghost_schedule.h      # T05；octlb_field_schedule
│       │   └── face_iterator.h       # T05；CoarseFineFaces 遍历
│       │
│       ├── io/           ← Solver 模块内共用 IO（不依赖 LBM 头文件）
│       │   ├── vtk_writer/   # 移植 blockVtkWriter3D，替换驱动层
│       │   └── gnuplot/      # 直接复制自 OpenLB
│       │
│       └── lbm/          ← LBM 专属；算法复用策略（见"OpenLB 代码融入策略"）
│           │             # ── OctLB 自建区（namespace octlb）──
│           ├── block_lattice.h/.cpp  # octlb::BlockLattice<T,D>，octant 感知格子数据结构
│           ├── level_coupler.h/.cpp  # Lagrava 粗细耦合（T06；算法复用 dynamics/）
│           ├── time_loop/            # 递归下降 TimeLoop（T06）
│           ├── unit_converter/       # LBM 单位换算（cavity3d 集成阶段）
│           │             # ── OpenLB 算法头文件（只读，namespace olb，无 .cpp）──
│           ├── descriptor/   # D3Q19 等格子常数（c/t/invCs2），~10 个头文件
│           ├── dynamics/     # BGK/MRT/Smagorinsky kernel、平衡态、矩计算，~12 个头文件
│           ├── core/         # MinimalCell concept、FieldD、Vector，~5 个头文件
│           ├── utilities/    # 数学工具（vectorHelpers 等），~5 个头文件
│           └── boundary/     # Bouzidi/bounce-back 算法头文件（T09 启用）
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

### 基础类型（src/common/）

全项目共用的基础类型，不依赖任何第三方库头文件：

```cpp
// types.h
namespace octlb {
using label    = int32_t;   // 整数索引，对应 OpenFOAM label
using scalar   = double;    // 浮点数，对应 OpenFOAM scalar
using OctantId = label;     // 本地 quadrant 线性序号（0…local_num_quadrants-1）

// 面方向编码（与 p4est 内部整数 0-5 一一对应，避免魔法数字）
enum class FaceDir : int {
  kXMin = 0, kXMax = 1,
  kYMin = 2, kYMax = 3,
  kZMin = 4, kZMax = 5,
};
}

// bounding_box.h
namespace octlb {
struct BoundingBox {
  scalar x_min, y_min, z_min;
  scalar x_max, y_max, z_max;
};
}
```

`OctantId` 是 p4est 本地 quadrant 的连续线性序号，与 `sc_array_t` 下标直接对应，`BlockCollection<T>` 可用 `std::vector<T>` 按下标 O(1) 访问。

---

### 模块一：Mesh 模块

#### OctreeForest

- 封装 p4est 的 `p8est_t` forest，提供 refine / balance / partition 操作。
- 构造函数：`OctreeForest(MPI_Comm comm, BoundingBox domain)`；内部固定使用 `p8est_connectivity_new_unitcube()`，物理坐标由 `domain` 缩放。
- `refine()` 签名：`void refine(std::function<bool(OctantId)> criterion, int max_level)`；内部通过 `user_pointer` 桥接 p4est C 回调。
- `partition()` 签名：`void partition(std::function<int(OctantId)> weight_fn = nullptr)`；`weight_fn == nullptr` 时均匀分区（`p8est_partition`），否则加权分区（`p8est_partition_ext`）；每次调用后自动重建 ghost 层。
- 对外暴露接口：`local_num_octants()`, `quadrant_bounds(OctantId)`, `quadrant_level(OctantId)`。
- `quadrant_bounds()` 返回 `BoundingBox`，将 p4est 整数坐标按 `domain` 缩放为物理坐标。
- 内部持有 `p8est_ghost_t*`，在 `partition()` 后重建、析构时销毁；仅供 Mesh 模块内部（`FacePairList`）通过模块内访问器使用。
- 不暴露任何 `p8est_*` 类型给 Solver 模块。

#### FacePairList

- 通过 `p8est_iterate` + face callback 一次性遍历所有面，构建两类列表：
  - **SameLevelFaces**：`{local_block_id, face_dir, remote_block_id, remote_rank, comm_tag}`
    - `comm_tag`：在 face callback 内由两侧 quadrant 几何信息生成的**对称** MPI tag（T05）；跨 rank send/recv 成对使用，同 rank 忽略。
    - 跨 rank 时 `remote_block_id` 为 p4est ghost 索引，**不是**远端 rank 的 `OctantId`。
  - **CoarseFineFaces**：`{coarse_block_id, fine_block_ids[4], face_normal, remote_ranks[4], comm_tags[4]}`
    - `comm_tags[i]`：在 face callback 内为 coarse 与 `fine_ids[i]` 之间生成的**对称** MPI tag（T06）；跨 rank 时 `LevelCoupler` 交换 coarse 宏观量使用，同 rank 忽略。
  - 粗细面识别依据：face callback 中 `is_hanging == 1` 的一侧为 4 个细格，`is_hanging == 0` 的一侧为 1 个粗格；**不**进入 `SameLevelFaces`。
- 静态分层下初始化一次，预留 `rebuild()` 接口供动态 AMR 扩展；`rebuild()` 后须重建 `GhostSchedule`、`LevelCoupler`（及 `TimeLoop` 层缓存）；通信计划与 tag 绑定拓扑。

#### GeometryEngine / MaterialField（T07）

**输入与多部件**

- `stl_reader`：ASCII/binary STL → `TriangleSoup`（移植 OpenLB；见 `mesh/io/stl_reader/`）。
- `GeometryAssembly`：一个或多个 `GeometryPart`（每个 part 一份 soup 或 STL 路径），**不**固定为两个文件。风洞典型配置：洞壁 STL → `kInternalChannel`；试验件 STL → `kExternalObstacle`（更高 `priority` 覆盖洞内流体区）。
- `GeometryPartRole`：
  - `kExternalObstacle`：实心闭体（如 cylinder、试验件）；STL **内** → solid，**外** → fluid。
  - `kInternalChannel`：第一版仅 **闭壳墙材**（watertight）；CGAL **inside** → solid（墙），**outside** → fluid（洞腔/管腔）。不支持「流体域 STL」（`kFluidLumenVolume`，P2）。

**几何自适应加密（静态 AMR，T07）**

- 参考 octree-mesh `GeometryAdaptiveEngine`：`ResolveBounding`（扩展包围盒 + wake）+ `ResolveSurface`（多轮表面相交 → `OctreeForest::refine`），每轮 `balance()`，表面附近密、远场疏；加密判据为 **全部 part 三角面并集**。
- 完成后 `partition(make_level_weight_fn(forest))`。
- `GeometryEngine::build(forest, assembly, config)` 修改 forest 并产出 `MaterialField`；**不**在 `build()` 内构造 `FacePairList`。调用方须在 `build()` 后重建 `FacePairList`、`GhostSchedule`、`LevelCoupler`（及 `TimeLoop` 层缓存，同 T06）。

**体素化与材料**

- CGAL 体素化（裁剪 octree-mesh `voxelization` 思路）：每 octant 块内 N³ 格 `MaterialKind`：`fluid` / `solid` / `boundary`。
- `boundary`：与三角网格 **表面相交** 的格（同 octree-mesh `CellVoxelizer`）；**不**预计算 Bouzidi 距离/法向（T09）。
- 分 part 标量场后按 `priority` 合并；冲突 **solid > boundary > fluid**。
- `MaterialField` 仅存本地 `OctantId` 对应块；Solver 初始化时**单次读取**，之后不再持有引用。

#### WeightedLoadBalancer

- 提供 free function `make_level_weight_fn(const OctreeForest&)`，返回 weight(octant) = 2^level 的 lambda。
- 调用方：`forest.partition(make_level_weight_fn(forest))`，在 `balance()` 完成后执行。
- 第一版静态分层只需初始化时调用一次。
- `WeightedLoadBalancer` 自身不持有任何 p4est 类型；分区操作由 `OctreeForest::partition()` 统一执行。

---

### 模块二：Solver 模块

#### 核心索引约定

所有组件使用统一的格子地址 `(octant_id, i, j, k)`，无需翻译表：
- `octant_id` → `BlockCollection` 中的 `octlb::BlockLattice`
- `(i, j, k)` → 块内坐标，直接映射到 `BlockLattice::get(i, j, k)`

这是与 octree-mesh 中 `mb_to_mesh_` / `mesh_to_mb_` 翻译表的根本性改进。

#### field/ 子层（泛型，不依赖 LBM）

- **`BlockCollection<T>`**：管理本 rank 持有的所有块，每块存一个 `T` 类型对象，以 `octant_id` 为键。
  - 内部使用 `std::vector<T>`，按 OctantId 下标 O(1) 访问。
  - 构造签名：`BlockCollection(label num_octants, std::function<T(OctantId)> factory)`；
    工厂函数在构造时逐 id 调用，不要求 T 默认可构造（兼容 `ConcreteBlockLattice` 等需要参数的类型）。
  - CMake target：`octlb_field` INTERFACE 库（`BlockCollection`、`BlockIterator`、`FacePackable` concept），不依赖 `octlb_mesh`、P4est、MPI。
- **`GhostSchedule<T>`**（T05）：从 `FacePairList::SameLevelFaces` **构造时**固化通信计划；`exchange()` 执行同级 halo 交换。
  - **范围**：仅同级、同尺寸块的 1:1 面；**不**处理 `CoarseFineFaces`（粗细 1:4 由 T06 `LevelCoupler` / Lagrava 负责）。
  - 交换单位：每面单层 N×N 字段（LBM 为 `N²×Q`），非整块 `N³×Q`。
  - **FacePackable**：`T` 提供 `pack_face(dir, buf, n)` / `unpack_face(dir, buf, n)`；pack 读 interior 最外层，unpack 写紧贴 interior 的 ghost 第一层（D3Q19 每步刷新 1 层，与 halo 总深度 `h` 无关）。
  - **同 rank**：`remote_rank == my_rank` 时 pack 后直接 `unpack` 到 `remote_id` 块的 `opposite(dir)` ghost，不走 MPI。
  - **跨 rank**：使用 `SameLevelFace::comm_tag` 匹配 `Isend`/`Irecv`；每步热路径无计划重建。
  - CMake target：`octlb_field_schedule` INTERFACE，依赖 `octlb_field`、`octlb_mesh`、MPI。
- **`BlockIterator`**：遍历本 rank 所有 `octant_id`（产出 OctantId 整数值，不是 T 引用）。
  - Level 过滤**不在** BlockIterator 内实现；TimeLoop 初始化时从 OctreeForest 缓存
    `level → [OctantId]` 映射，按层驱动 collide → halo → stream。
- **`FaceIterator`**（T05）：遍历 `CoarseFineFaces`，提供 `(coarse_block, fine_blocks[4], normal)` 视图，供 T06 `LevelCoupler` 使用；无 LBM 依赖，无耦合 MPI。

#### lbm/ 子层（LBM 专属）

**BlockStore（LBM 具体化）**

- 每个 octant 持有一个 `octlb::BlockLattice<T, D3Q19Descriptor>`——OctLB 自建的格子数据结构。
- 内存布局：`[Nx+2h][Ny+2h][Nz+2h][Q]` 平坦数组，`h=1` ghost halo，Q 在最内层。
- 块尺寸：`(N+2)³×Q`，halo 宽度 = 1（D3Q19 streaming 只需 1 格 overlap）。
- `BlockLattice::collide(omega)` 内部遍历所有非 ghost 格，调用 `olb::collision::BGK::type::apply(cell, params)`——OpenLB 的纯函数模板，只依赖 `concepts::MinimalCell`。
- `BlockLattice::pack_face` / `unpack_face`（T05）：满足 `FacePackable`，供 `GhostSchedule` 与测试使用。
- `BlockLattice::stream()` 使用 pull-scheme，读取已由 `GhostSchedule` 刷新的 ghost 层（`collide` 之后、`stream` 之前调用 `exchange()`）。
- 初始化时按 `MaterialField` 设置每格 omega（fluid）或标记为 solid/boundary（T09 实现 Bouzidi，依赖 T07 `MaterialField`）。
- `BlockCollection<BlockLattice>` 从 T03 泛型工厂直接复用，无需修改。

**HaloExchange**

- 每层时间步：`collide()` → `GhostSchedule::exchange()`（同级面）→ `stream()`。
- 面层数据格式：19 个 double/float 的 N×N slab（单层），对应 f 分布函数分量。
- 粗细界面不在此路径；由 `LevelCoupler` 在子步间对 `CoarseFineFaces` 做 Lagrava 插值/限制。

**LevelCoupler（Lagrava 方案，T06）**

**算法复用 + 自建数据结构**（与 T04 BGK 同路线）：按 Lagrava 公式在 `octlb::LevelCoupler` 内实现 prolongation/restriction；宏观量与平衡态复用 `olb::lbm::*` + `CellProxy`。**不**复制 OpenLB `refinement/` 框架头文件（`BlockRefinementContext` 等）。

- **耦合面来源**：`FacePairList::CoarseFineFaces` + `FaceIterator`（由 `p8est_iterate` face callback 产生），替代 OpenLB 原版的 `SuperIndicatorDomainFrontierDistanceF`。
- **连接点索引（Solver 侧）**：`LevelCoupler` 构造时从 `FaceIterator` + `OctreeForest::quadrant_bounds()` + 块尺寸 N 固化 **`coupling_plan_`**（cell 级 `(coarse_id,i,j,k) ↔ (fine_id,i,j,k)`）；Mesh 保持 block 级纯拓扑，不存 cell 级 coupling 列表。
- **粗→细（prolongation）**：
  - `scalingFactor = (τ - 0.25) / τ`（第一版全域统一 `τ = 1/omega`）
  - `f_fine = f_eq(ρ, u) + scalingFactor × f_neq`
  - half-time：粗侧宏观量取 prev 与 curr 平均；full-time：使用 curr
  - 在细层第一个子步前（half-time）和第二个子步前（full-time）各执行一次。
- **细→粗（restriction）**：
  - `scalingFactor = τ / (τ - 0.25)`
  - 对 2:1 子格（界面参与子集）取均值后重建宏观量，写回粗格分布函数。
  - 在两个细子步完成后执行。
- **跨 rank MPI**：使用 `CoarseFineFace::comm_tags[i]`；构造时固化 MPI 计划（coarse 宏观量 ρ、u、f_neq → fine 侧），模式对齐 `GhostSchedule`（Irecv → Isend/本地直写 → Waitall）。
- **拓扑 rebuild**：`FacePairList::rebuild()` 后须销毁并重建 `LevelCoupler`（与 `GhostSchedule` 同级）；第一版静态 AMR 不实现热路径。

**TimeLoop（递归下降，T06）**

```
advance(level l):
  collide(all blocks at level l)              // BlockLattice::collide()
  halo_exchange(level l)                      // GhostSchedule::exchange()，仅 SameLevelFaces
  stream(all blocks at level l)               // BlockLattice::stream()
  coupler.apply_half_time(l → l+1)            // LevelCoupler prolongation half-time
  advance(level l+1)                          // 第一个细子步（递归）
  coupler.apply_full_time(l → l+1)            // LevelCoupler prolongation full-time
  advance(level l+1)                          // 第二个细子步（递归）
  coupler.restrict(l+1 → l)                   // LevelCoupler restriction
```

- **`TimeLoop` 类 + 外部引用**：构造时注入 `BlockCollection`、`GhostSchedule`、`LevelCoupler`（调用方持有 Schedule/Coupler，便于 `rebuild()` 后整体替换）；`TimeLoop` 不创建或拥有上述对象。
- **对外 API**：`advance_one()` 从 level 0 递归执行一个粗层时间步；构造时缓存 `level → [OctantId]`。
- 递归深度 = 细化层数（6-8，最大 10），不会溢出栈。
- 每层块上顺序为 `collide()` → `GhostSchedule::exchange()` → `stream()` 的并行循环。
- 第一版：静态分层，层级结构在初始化时确定。
- **Level-to-OctantId 映射**：TimeLoop 在初始化时调用 `OctreeForest::quadrant_level(id)`
  遍历所有本地 octant，将 `level → [OctantId]` 映射缓存为 `std::vector<std::vector<OctantId>>`，
  之后每步直接按层索引，不再查询 OctreeForest。BlockCollection 和 BlockIterator 本身对 level 无感知。
- **测试 hook**：各层 `collide`/`stream` 调用计数（或等价可观测接口），供 `test_time_loop_levels` 验证 L=3 时步数比 1:2:4。

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

**融合方式：算法复用（Algorithm Reuse）**

OctLB 自建数据结构（`octlb::BlockLattice`），仅从 OpenLB 借用**无状态函数模板**（碰撞算子、
平衡态公式、矩计算），不引入 OpenLB 的格子管理框架（`ConcreteBlockLattice`、
`DynamicsPromise`、`SuperLattice`）。

**设计原因**

OpenLB 1.9.0 的 `ConcreteBlockLattice` 经 `blockLattice.h → analyticalF.h →
superGeometry.h` 的传递包含链拖入整个应用框架层；`DynamicsPromise` 在构造函数体内的
局部 lambda 类型强制要求所有模板实例化在同一翻译单元（GCC "local type, used but never
defined" 规则），唯一解是引用 `olb3D.h + olb3D.hh` umbrella，导致约 1477 个文件、16 MB
代码进入仓库。这与"按需融合"原则冲突。

**实际解**：OpenLB 的 BGK/MRT 碰撞 kernel 只依赖 `concepts::MinimalCell`——
即 `cell[iPop]` 操作符。OctLB 实现满足该 concept 的 `CellProxy`，直接调用
`olb::collision::BGK::type::apply(cell, params)`，彻底绕开框架层。

**只复制算法头文件**（约 25 个，无 `.cpp`，无框架层）：

| 目录 | 来源 | 内容 | 文件数 |
|---|---|---|---|
| `lbm/descriptor/` | OpenLB | D3Q19/D2Q9 等格子常数（c/t/invCs2），纯编译期 | ~10 |
| `lbm/dynamics/` | OpenLB | BGK/MRT/Smagorinsky kernel、平衡态、矩计算 | ~12 |
| `lbm/core/` | OpenLB | MinimalCell concept、FieldD、Vector 类型 | ~5 |
| `lbm/utilities/` | OpenLB | 数学工具（vectorHelpers、normSqr 等） | ~5 |
| `lbm/boundary/` | OpenLB | Bouzidi/bounce-back 算法头（T09 启用） | ~3 |

> **T06 Lagrava**：不复制 `refinement/`；`LevelCoupler` 按 PRD 公式实现，复用已移植的 `dynamics/lbm.h`（`computeRhoU`、`computeFneq`、`equilibrium`）。

**不复制**：`blockLattice.h`、`superLattice.h`、`case/`、`optimization/`、`particles/`、
`reaction/`、`uq/`、`functors/analytical/`、`geometry/superGeometry.h`、`communication/`（MPI
单例由 OctLB 自己管理或由 p4est 初始化）、`io/`（无 tinyxml2 依赖）。

**CMake target**（极简）：

```cmake
add_library(octlb_lbm STATIC
    block_lattice.cpp          # OctLB 自建格子，唯一 .cpp
)
target_include_directories(octlb_lbm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(octlb_lbm PUBLIC cxx_std_20)
target_link_libraries(octlb_lbm PUBLIC MPI::MPI_CXX)
# 不需要 PLATFORM_CPU_SISD / OLB_VERSION / tinyxml2
```

**关键约定**：

- **Namespace**：OpenLB 算法头文件保留 `namespace olb`；OctLB 自写代码使用 `namespace octlb`。
- **不修改 OpenLB 头文件**：只复制，不改写，保持与上游的 diff 可追踪性。
- **版本基准**：OpenLB 1.9.0 的算法头；后续升级只需 diff 对应的 ~25 个头文件。
- **GPU**：不定义 `PLATFORM_GPU_CUDA` / `PLATFORM_GPU_HIP`，预留 DESCRIPTOR 模板参数。
- **IO 模块归属**：`solver/io/vtk_writer/` 和 `solver/io/gnuplot/` 不依赖 LBM 头文件（T08 实现）。

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
| 2a | `test_geometry_adaptive_refine` | 硬编码 triangle soup + `OctreeForest`：表面附近叶层 `quadrant_level` 高于远场；可选 2-rank 集体 refine 不死锁 | Mesh/geometry |
| 2 | `test_geometry_engine` | 多场景 triangle soup（外流实心、内流方管、风洞 channel+obstacle 等）：`MaterialKind` 分布正确；非法洞壁 fixture 失败；不依赖大 STL 文件 | Mesh/geometry |
| 3 | `test_face_pair_list` | p8est_iterate face callback 正确识别同级面和粗细 hanging 面 | Mesh/topology |
| 4 | `test_load_balancer` | 2^level 权重分区后各 rank 总权重差异 < 5% | Mesh |
| 5 | `test_ghost_topology` | FacePairList 中每条跨 rank 面对的 remote_rank 与 p4est ghost 层记录一致（本地校验）；精心细化模式确保 hanging 面和跨 rank 面必然出现 | Mesh |
| 6 | `test_block_collection` | `BlockCollection<T>` 增删查，以 `octant_id` 为键；`BlockIterator` 遍历正确 | Solver/field |
| 7 | `test_collision_bgk` | 单块 BGK 碰撞：质量守恒、动量守恒，收敛至 Maxwell 平衡态 | Solver/lbm |
| 8 | `test_ghost_schedule` | L1：`DummyBlock` 覆盖同 rank / 跨 rank / 空计划 / 混合面 / 角点 / `comm_tag`；L2：2-rank `BlockLattice` 邻接面 populations 一致 | Solver/field |
| 9 | `test_lagrava_coupler` | L1：`coupling_plan_` 条数/索引/象限；L2：同 rank 1 粗+1 细界面 ρ/u 连续、质量守恒；L3：2-rank 跨 rank macro 交换 | Solver/lbm |
| 10 | `test_time_loop_levels` | L1：L=3 一次 `advance_one()` 后各层 collide 计数 1:2:4、coupler 调用顺序；P1：真实 BlockLattice 不崩溃 | Solver/lbm |
| 11 | `test_vtk_writer` | 硬编码 `BlockCollection` 输入，输出合法 VTK XML，空间范围与 `quadrant_bounds` 吻合 | Solver/io |
| 12 | `test_cavity3d_serial` | 单 rank，cavity3d，Re=100，中线速度剖面与 Ghia 1982 参考值误差 < 2% | Integration |
| 13 | `test_cylinder3d_parallel` | 多 rank（≥4），cylinder3d L=4，阻力系数 Cd 与 OpenLB 参考值误差 < 1% | Integration |
| 14 | `test_amr_convergence` | L=1→3 逐级细化，速度场 L2 误差收敛阶 ≈ 2 | Integration |

**依赖链**：0 → 3 → 5 → 8 → 9 → 10；0 → 2a；1 → 2a；2a → 2；1 → 2；6 → 8；7 → 9；0–11 全绿后进入 Integration。

**任务拆分（`doc/tasks/`，与测试序对照）**：T07 Mesh 几何链（#1、#2a、#2）→ T08 `vtk_writer`（#11）→ T09 Bouzidi + Lattice 材料初始化 → T10+ 集成算例（#12–#14）。T08 可与 T07 部分并行；T09 阻塞于 T07。

**测试基础设施**：GTest + MPI。
- `tests/mpi_main.cpp`：手写 `MPI_Init → RUN_ALL_TESTS → MPI_Finalize`，所有 MPI 测试 target 链接此 main 而非 `GTest::gtest_main`。
- CMake 用 `${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 --oversubscribe` 启动测试，跨平台（本地 `mpirun`、服务器 `srun` 均适用）。
- p4est 集成：`cmake/FindP4est.cmake` 暴露 `P4est::p4est` imported target（INTERFACE 依赖 `MPI::MPI_C`）；跨机器用 `-DP4EST_ROOT=...` 指定安装路径。

---

## 不在范围内（第一版）

- **动态 AMR**：运行时 refine/coarsen 及字段迁移（`p4est_transfer_fixed`）。预留 `FacePairList::rebuild()` 及 Solver 侧 `GhostSchedule` / `LevelCoupler` / `TimeLoop` 重建契约，第一版不实现热路径。
- **GPU backend**：OpenLB 的 `Platform::GPU_CUDA` / `GPU_HIP`。预留 `Platform` 模板参数，不测试。
- **HDF5 checkpoint/restart**：预留 `io/checkpoint/` 目录，不实现。
- **FVM 或其他物理模型**：`field/` 子层设计为可复用，但第一版只跑 LBM。
- **流场触发的自适应细化**：如 Q-criterion、速度梯度触发（第一版静态分层由 **几何表面** 驱动，见 T07 `GeometryEngine`）。
- **流体域 STL / 非闭壳薄壳洞壁**：T07 仅支持闭壳墙材 `kInternalChannel`；`kFluidLumenVolume`、通用 CSG 为 P2。
- **多相流、热 LBM 等扩展物理模型**：第一版只验证单相不可压 NS。

---

## 其他说明

### 与 octree-mesh (SAMR) 的关系

OctLB 是对 octree-mesh 的**重新设计**，不是直接扩展。以下 octree-mesh 组件可作为参考或局部复用：

| octree-mesh 组件 | OctLB 对应 | 处理方式 |
|---|---|---|
| `OctreeMesh` + p4est 封装 | `OctreeForest` | 参考重写，正确使用 `p8est_iterate` face callback |
| `GeometryAdaptiveEngine` + CGAL | `GeometryEngine`（T07） | 参考复用自适应加密 + CGAL 体素化；多 STL 用 `GeometryAssembly` Part 角色 |
| `MBArray` / `BlockField` | `BlockCollection<octlb::BlockLattice>` | 废弃，换 OctLB 自建格子存储 |
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
| OctreeForest | `mesh/forest/` | 测试绿 | `test_octree_forest` | T01；CI main 已通过 |
| stl_reader | `mesh/io/stl_reader/` | 未开始 | `test_stl_reader` | T07 W1；移植自 OpenLB |
| GeometryEngine + MaterialField | `mesh/geometry/` | 未开始 | `test_geometry_engine`、`test_geometry_adaptive_refine` | T07；自适应加密 + Part 合并 + CGAL 体素化 |
| FacePairList | `mesh/topology/` | 测试绿 | `test_face_pair_list` | T02；T05 `SameLevelFace::comm_tag`；T06 `CoarseFineFace::comm_tags[4]` + `coarse_remote_rank` |
| WeightedLoadBalancer | `mesh/load_balance/` | 测试绿 | `test_load_balancer` | T02；本地 4-rank 通过，待 CI |
| Ghost 拓扑一致性 | `mesh/topology/` | 测试绿 | `test_ghost_topology` | T02；本地 4-rank 通过，待 CI |

### Solver/field 模块

| 组件 | 所在路径 | 状态 | 对应测试 | 备注 |
|---|---|---|---|---|
| BlockCollection\<T\> + BlockIterator | `solver/field/` | 测试绿 | `test_block_collection` | 泛型，不依赖 LBM 头文件 |
| GhostSchedule\<T\> + FaceIterator | `solver/field/` | 测试绿 | `test_ghost_schedule` | T05；`octlb_field_schedule`；同级面层 MPI |
| FacePairList `comm_tag` | `mesh/topology/` | 测试绿 | `test_face_pair_list`（扩展） | T05 扩展 T02 |
| BlockLattice `pack_face` / `unpack_face` | `solver/lbm/` | 测试绿 | `test_ghost_schedule` L2 | T05 |

### Solver/lbm 模块

| 组件 | 所在路径 | 状态 | 对应测试 | 备注 |
|---|---|---|---|---|
| `octlb::BlockLattice` + BGK | `solver/lbm/block_lattice.*` | 测试绿 | `test_collision_bgk` | T04 已落地并通过本地 ctest 验证 |
| MRT / Smagorinsky Dynamics | `solver/lbm/dynamics/` | 未开始 | — | 复用 OpenLB 算法头文件，集成测试间接覆盖 |
| Bouzidi BC | `solver/lbm/boundary/` | 未开始 | — | T09；复用 OpenLB boundary/；依赖 T07 `MaterialField` |
| LevelCoupler（Lagrava） | `solver/lbm/level_coupler.*` | 测试绿 | `test_lagrava_coupler`、`test_lagrava_coupler_mpi2` | T06；算法复用 dynamics/，Solver 侧 coupling_plan_ |
| TimeLoop（递归下降） | `solver/lbm/time_loop/` | 测试绿 | `test_time_loop_levels` | T06；外部引用 GhostSchedule/LevelCoupler；v1 每步全局 `GhostSchedule::exchange()` |
| CoarseFineFace `comm_tags[4]` | `mesh/topology/` | 测试绿 | `test_face_pair_list`、`test_lagrava_coupler_mpi2` | T06 扩展 T02 |
| unit_converter | `solver/lbm/unit_converter/` | 未开始 | — | cavity3d 集成阶段；不在 T06 |

### Solver/io 模块

| 组件 | 所在路径 | 状态 | 对应测试 | 备注 |
|---|---|---|---|---|
| vtk_writer | `solver/io/vtk_writer/` | 未开始 | `test_vtk_writer` | T08；移植 blockVtkWriter3D，替换驱动层 |
| gnuplot | `solver/io/gnuplot/` | 未开始 | — | 直接复制自 OpenLB |

### 集成测试

| 测试 | 状态 | 前置条件 | 验收标准 |
|---|---|---|---|
| `test_cavity3d_serial` | 未开始 | 测试 0–11 全绿 | Re=100 中线速度剖面与 Ghia 1982 误差 < 2% |
| `test_cylinder3d_parallel` | 未开始 | `test_cavity3d_serial` 通过 | ≥4 rank，Cd 与 OpenLB 参考值误差 < 1% |
| `test_amr_convergence` | 未开始 | `test_cylinder3d_parallel` 通过 | L=1→3 速度场 L2 误差收敛阶 ≈ 2 |
