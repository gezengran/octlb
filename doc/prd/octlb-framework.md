# OctLB 框架 PRD

> 版本：v0.3  
> 日期：2026-05-28（进度更新：2026-06-24；T10 完成，见 `doc/tasks/T10-cavity3d.md`）  
> 推敲更新：2026-07-10——审计已实现栈发现缺陷 ①②③④⑤ + 根因，确立 T11 cylinder3d 治理路径（P0 + W1–W4，oracle 三阶 sanity→量级→Cd<1%），见 `doc/tasks/T11-cylinder3d.md` 与下文「已知缺陷与治理路径」节  
> 架构更新：2026-07-15——T11 W3 推进时发现 BC 调度架构缺陷 ⑥（per-face `DomainBcSpec` + 整 handler `boundary_lattice_mode_` 无法表达同 block 内混合 BC，每接入新算例需在唯一 FD 路径加分支）。决策：移植 OpenLB「per-cell Dynamics 多态 + material 号映射」**架构思想**为 OctLB「per-cell `BcKind` 枚举 + 中心化 dispatch」（**不引框架层**、无虚函数；枚举 switch 失效再切 OpenLB 验证过的多态兜底）。`MaterialKind` 保持几何 only，BC 角色求解器侧解析。先重构，压力出口（T11 W3 组件 2）作首个新客户。详见 `doc/tasks/T11-refactor-bc-dispatch.md` 与下文 BlockStore/DomainBoundaryHandler/缺陷 ⑥

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
│       │   ├── vtk_writer/   # StructuredGrid .vts + vtm/pvd（T08）；不复制 OpenLB ImageData 路径
│       │   └── gnuplot/      # 直接复制自 OpenLB
│       │
│       └── lbm/          ← LBM 专属；算法复用策略（见"OpenLB 代码融入策略"）
│           │             # ── OctLB 自建区（namespace octlb）──
│           ├── block_lattice.h/.cpp  # octlb::BlockLattice<T,D>，octant 感知格子数据结构
│           ├── cell_kind.h           # CellKind（T09-W1）
│           ├── domain_boundary_handler.*  # 域外 tree BC（T09-W1）
│           ├── bouzidi_link_data.*   # Bouzidi q_frac 缓存（T09-W2）
│           ├── lattice_material_init.*  # MaterialField → CellKind（T09-W2）
│           ├── level_coupler.h/.cpp  # Lagrava 粗细耦合（T06；算法复用 dynamics/）
│           ├── time_loop/            # 递归下降 TimeLoop（T06；T09-W1 注入 DomainBoundaryHandler）
│           ├── unit_converter/       # LBM 单位换算（T10 cavity3d）
│           │             # ── OpenLB 算法头文件（只读，namespace olb，无 .cpp）──
│           ├── descriptor/   # D3Q19 等格子常数（c/t/invCs2），~10 个头文件
│           ├── dynamics/     # BGK/MRT/Smagorinsky kernel、平衡态、矩计算，~12 个头文件
│           ├── core/         # MinimalCell concept、FieldD、Vector，~5 个头文件
│           ├── utilities/    # 数学工具（vectorHelpers 等），~5 个头文件
│           └── boundary/     # bounce-back / Zou-He（T09-W1）；Bouzidi pull（T09-W2）
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
  - **TreeBoundaryFaces**（T09-W1）：`{octant_id, face_dir}`；face callback 中 `tree_boundary != 0` 的域外树界面。**不**进入 `SameLevelFaces` / `CoarseFineFaces`；BC 类型由 example 配置，Mesh 只枚举拓扑。
  - 粗细面识别依据：face callback 中 `is_hanging == 1` 的一侧为 4 个细格，`is_hanging == 0` 的一侧为 1 个粗格；**不**进入 `SameLevelFaces`。
- 静态分层下初始化一次，预留 `rebuild()` 接口供动态 AMR 扩展；`rebuild()` 后须重建 `GhostSchedule`、`LevelCoupler`、`DomainBoundaryHandler`（T09-W1）、`BouzidiLinkData`（T09-W2）及 `TimeLoop` 层缓存；通信计划与 tag 绑定拓扑。

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
- 每 interior 格持有 **`CellKind`**（T09-W1）：`kFluid` / `kSolid` / `kBoundary`；默认 `kFluid`。`collide()` / `stream()` **跳过** `kSolid` 与 `kBoundary`。
  - **`BcKind` 超集重构（2026-07-15 立项，T11-refactor）**：`CellKind` 将被 `BcKind`（`kBulk`/`kBounceBack`/`kMovingBounceBack`/`kBouzidi`/`kVelocityDirichlet`/`kPressureDirichlet`/`kOutflow`/`kSolid`）取代——per-cell 标记 BC 类型，中心化 `BcDispatcher` 按 `BcKind` switch 分发 collide/computeRhoU/post_stream，**移除整 handler `boundary_lattice_mode_` 全局门控**，使同 block 内混合 BC（入口速度 + 出口压力 + 圆柱 Bouzidi + 壁 bounce-back）自由共存。`MaterialKind` 保持几何 only，BC 角色由求解器侧 `BcInstaller` 在 setup 时按 `(MaterialField + TreeBoundaryFace 面角色 + BouzidiLinkData)` stamp 为 per-cell `BcKind`。`BcDispatcher`/`BcInstaller` header-only（`lbm/` 唯一 `.cpp` 仍 `block_lattice.cpp`）。详见 `doc/tasks/T11-refactor-bc-dispatch.md`。
- `BlockLattice::collide(omega)` 在 `kFluid` 格上调用 `olb::collision::BGK::type::apply(cell, params)`——OpenLB 的纯函数模板，只依赖 `concepts::MinimalCell`。
- `BlockLattice::pack_face` / `unpack_face`（T05）：满足 `FacePackable`，供 `GhostSchedule` 与测试使用。
- `BlockLattice::stream()` 使用 pull-scheme，读取已由 `GhostSchedule` 刷新的 ghost 层及域边界 BC 写入的 tree ghost（见下）。**T09-W2**：pull 时若源邻居为 `kSolid` / `kBoundary`，用初始化缓存的 `q_frac` 在 `stream()` 内做 Bouzidi（**不**从 solid 格读 population）。
- **T09-W2**：`initialize_from_material(MaterialField)` 写入 `CellKind` 并为 fluid 格设统一 `omega`；**T09-W1 不**读取 `MaterialField`。
- `BlockCollection<BlockLattice>` 从 T03 泛型工厂直接复用，无需修改。
- **Overlap-padding（2026-07-15，T11 ④ HITL 决策落定）**：`BlockLattice` 在 core + `h=1` halo 之外，另带一层对齐 OpenLB `SuperLattice` padded-block 语义的 **overlap padding**（core + 2×OVERLAP，OVERLAP 默认 3）。相关机制：`collide_overlap_padding_bgk` / `stream_overlap_padding_shell` / `commit_overlap_padding_stream` / `OverlapPaddingCollideMode`（`kBgkOnMaterialNonZero` / `kNoDynamics`）/ `YminYmaxPaddingOutOfHaloMode`（`kOpenLbRotateWrap` / `kZero` / `kKeepSelf`）/ `fill_overlap_padding_from_core` / `fill_overlap_padding_bc_post_stream`（`solver/lbm/block_lattice.h`）。**保留原因**：cavity3d（T10）的 Cd ~2% 对齐 OpenLB 校准依赖此 padded 语义（stream 按列 rotate、PostStream `addPoints2CommBC` 与 padding 通信）；早期「halo `h=1`、无 padding」的最简假设不足以复现 OpenLB 结果。**设计原则不变**：padding 仍是 `octlb::BlockLattice` **自建存储**内部实现，**不**引入 OpenLB 框架层（`ConcreteBlockLattice`/`SuperLattice`/`DynamicsPromise`），不破坏 Mesh↔Solver seam、`field/` 可抽取性、`io/` 无 LBM 等硬不变量。即「最小自建」立论修正为「**自建存储 + 不引框架层**」，padding 属存储内部细节，见「OpenLB 代码融入策略」修订说明。

**DomainBoundaryHandler（T09-W1）**

- 构造时绑定 `FacePairList::TreeBoundaryFaces` + example 侧 BC 配置（如 `FaceDir → no-slip bounce-back` / `moving_lid` Zou-He）。
- `apply()` 在 post-collision populations 上写入 **tree 外侧 ghost 层**；**不**处理 STL 曲面（W2 Bouzidi）或同级 halo（T05）。
- **注入 `TimeLoop`**（外部引用、不持有）：`collide → domain_bc.apply()（pre-stream）→ GhostSchedule::exchange() → stream() → domain_bc.apply_post_stream()`；单元测试可用 `NoOpDomainBoundaryHandler`。`apply_post_stream()` 为 OpenLB PostStream 对齐段，承载 InterpolatedVelocity/PlaneFd 等**需 stream 后邻居**重建的 BC（见 TimeLoop 阶段顺序修订说明）。
- **per-cell dispatch 重构（2026-07-15 立项，T11-refactor）**：`apply()`/`apply_post_stream()` 将改为遍历边界格按 per-cell `BcKind` 中心化分发（`BcDispatcher::collide` / `post_stream`），**移除**整 handler `boundary_lattice_mode_` 与 `UsesInterpolatedVelocity()` 全局门控及 legacy `ApplyLegacyFaceBc` per-face ghost-fill 路径。prescribed 值仍由面 → `DomainBcSpec` 查表（`u_wall`/`inlet_field`/`rho_target`），与 `BcKind` dispatch 正交。详见 `doc/tasks/T11-refactor-bc-dispatch.md`。

**BouzidiLinkData（T09-W2，静态 AMR v1）**

- 初始化一次性构建：对每个 **相邻 solid/boundary 的 `kFluid` 格**、每个 cut link 缓存 `q_frac`（沿 lattice 方向到壁面的分数距离）。输入：`MaterialField` + `GeometryAssembly` + `OctreeForest::quadrant_bounds` + 格心物理坐标；**不**回写 `MaterialField`（T07 不预存距离/法向）。
- 运行时集成在 `BlockLattice::stream()` pull 阶段；**不**再注入 `TimeLoop`。算法核复用 OpenLB / octree-mesh 无状态 Bouzidi 公式。
  - **`kBouzidi` arm 归并（2026-07-15 立项，T11-refactor）**：Bouzidi pull 将由 per-cell `BcKind==kBouzidi` 触发（`BcInstaller` stamp），不再作为独立分类路径；`BouzidiLinkData` 数据结构与 pull 公式不变。

**HaloExchange**

- 每层时间步：`collide()` → `DomainBoundaryHandler::apply()`（域外 tree 面 pre-stream 段，T09-W1）→ `GhostSchedule::exchange()`（同级面）→ `stream()`（含 T09-W2 Bouzidi pull）→ `DomainBoundaryHandler::apply_post_stream()`（PostStream 段，InterpolatedVelocity 等）。
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
  collide(all blocks at level l)              // BlockLattice::collide()，跳过 kSolid/kBoundary
  domain_bc.apply()                           // 域外 tree ghost 预填 + 边界格 Dirichlet collide（T09-W1，pre-stream 段）
  halo_exchange(level l)                      // GhostSchedule::exchange()，仅 SameLevelFaces
  stream(all blocks at level l)               // BlockLattice::stream()，含 Bouzidi pull（T09-W2）
  domain_bc.apply_post_stream()              // InterpolatedVelocity（PlaneFd）等需 stream 后邻居重建的 BC（OpenLB PostStream 对齐）
  coupler.apply_half_time(l → l+1)            // LevelCoupler prolongation half-time
  advance(level l+1)                          // 第一个细子步（递归）
  coupler.apply_full_time(l → l+1)            // LevelCoupler prolongation full-time
  advance(level l+1)                          // 第二个细子步（递归）
  coupler.restrict(l+1 → l)                   // LevelCoupler restriction
```

> **阶段顺序修订（2026-07-15，T11 ③ HITL 决策落定）**：本文档早期版本为
> `collide → exchange → domain_bc.apply → stream`（单段 BC、BC 在 halo 交换之后）。
> 已验证实现为 `collide → domain_bc.apply → exchange → stream → domain_bc.apply_post_stream`
> （`solver/lbm/time_loop/time_loop.h:175-209`）。差异根因：`DomainBoundaryHandler` 拆为
> pre-stream `apply()`（域外 ghost 预填 / 边界格 Dirichlet collide）与 post-stream
> `apply_post_stream()`（InterpolatedVelocity / PlaneFd 有限差分需用 **stream 后** 的邻居分布
> 重建缺失 populations）两段。早期需求未预见到 InterpolatedVelocity BC 对 PostStream 阶段的依赖。
> **决策：以已验证实现为基准**——保留两段式 BC 与 PostStream 阶段，本文档同步更新为正式 spec。
> 此修订不改变设计原则（仍 `TimeLoop` 持外部引用、`DomainBoundaryHandler` 注入、seam 不破），
> 仅在 TimeLoop seam 内增加一个明确的 post-stream 阶段。

- **`TimeLoop` 类 + 外部引用**：构造时注入 `BlockCollection`、`GhostSchedule`、`LevelCoupler`、`DomainBoundaryHandler`（调用方持有，便于 `rebuild()` 后整体替换）；`TimeLoop` 不创建或拥有上述对象。
- **对外 API**：`advance_one()` 从 level 0 递归执行一个粗层时间步；构造时缓存 `level → [OctantId]`。
- 递归深度 = 细化层数（6-8，最大 10），不会溢出栈。
- 每层块上顺序为 `collide()` → `DomainBoundaryHandler::apply()`（pre-stream）→ `GhostSchedule::exchange()` → `stream()` → `DomainBoundaryHandler::apply_post_stream()` 的并行循环（详见上伪代码及阶段顺序修订说明）。
- 第一版：静态分层，层级结构在初始化时确定。
- **Level-to-OctantId 映射**：TimeLoop 在初始化时调用 `OctreeForest::quadrant_level(id)`
  遍历所有本地 octant，将 `level → [OctantId]` 映射缓存为 `std::vector<std::vector<OctantId>>`，
  之后每步直接按层索引，不再查询 OctreeForest。BlockCollection 和 BlockIterator 本身对 level 无感知。
- **测试 hook**：各层 `collide`/`stream` 调用计数（或等价可观测接口），供 `test_time_loop_levels` 验证 L=3 时步数比 1:2:4。

**VTK Writer（IO，T08）**

- **格式**：`VTKFile type="StructuredGrid"`，每 octant 一个 **`.vts`**；rank 0 写 **`.vtm`**（多块索引）与 **`.pvd`**（时间序列）。**不**使用 OpenLB `blockVtkWriter3D` 的 `ImageData`（`.vti`）路径。
- **驱动层**：`BlockIterator` 遍历本 rank 全部 `octant_id`，替代 OpenLB `CuboidDecomposition` / `SuperF3D` 循环。
- **几何**：块内 `N×N×N` LBM cell → `(N+1)³` 角点 `Points` + `N³` hexahedron；角点由 `OctreeForest::quadrant_bounds(octant_id)` 线性插值（`x(i)=x_min+(x_max-x_min)*i/N` 等）；`WholeExtent` 点索引为 `0 N  0 N  0 N`。
- **字段**：写在 **`CellData`**（cell-centered，仅 interior `0…N-1`，不含 ghost halo）。io 层定义 **`VtkCellField3D` concept**（`vtk_name`、`vtk_components`、`sample_cell`）；单块 `.vts` 可含多个数组（如 velocity、pressure、density）。`solver/io/vtk_writer/` **不** include LBM 头文件；`BlockLattice` 适配在 `solver/lbm/`（T08 W4，供 T10+）。
- **并行**：各 rank 写本地 `.vts`（命名 `{base}_r{rank}_oct{local_id}_T{iT:05d}.vts`，避免 partition 后本地 `OctantId` 重复）；rank 0 通过 MPI 收集文件名列表写当步 `.vtm` 并追加 `.pvd`（编排参考 OpenLB `SuperVTMwriter3D`，不引入其 cuboid 依赖）。

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

> **「最小自建」立论修订（2026-07-15，T11 ④ HITL 决策落定）**：项目首要目标是**对齐 OpenLB 结果**
> （cavity3d Cd ~2%、cylinder3d Cd<1%）。前期需求分析无法预见每个实现细节——验证发现
> `BlockLattice` 需在 core+halo 之外带一层 **overlap padding**（对齐 OpenLB `SuperLattice` padded-block
> 的 stream-rotate / PostStream `addPoints2CommBC` 语义）才能复现 OpenLB 校准结果，否则 halo `h=1`
> 无 padding 的最简假设不足以满足测试条件。**决策：保留 overlap-padding 机制，本文档同步更新**
> （见 BlockStore overlap-padding 说明）。**「最小自建」立论据此精确化为「自建存储 + 不引框架层」**：
> `BlockLattice` 仍是 OctLB 自建、内存自管、`CellProxy` 满足 `MinimalCell`、直接调 OpenLB 无状态 kernel，
> padding 仅是存储内部细节，**不**等价于引入 `ConcreteBlockLattice`/`SuperLattice` 框架层。
> 代码简洁原则与设计原则不因此改变——硬不变量（Mesh 不含 LBM、Solver 不调 p4est、`field/` 可抽取、
> `io/` 无 LBM、`lbm/` 唯一 `.cpp`）保持不变；padding 相关 API 须保持最小、可测、不外泄到 `field/`/`io/`。

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
| `lbm/boundary/` | OctLB 自建 + 参考 OpenLB / octree-mesh | bounce-back、Zou-He、Bouzidi pull 无状态核（T09） | ~3 |

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
- **IO 模块归属**：`solver/io/vtk_writer/` 不依赖 LBM 头文件（**T08**）；`solver/io/gnuplot/` 不依赖 LBM，单独后续任务（非 T08 范围）。

---

## 测试决策

**好测试的定义**：只测对外行为（接口契约），不测实现细节（内部数据结构）。**例外（2026-07-10 推敲纠正）**：跨 rank 契约（comm_tag 对称性——两侧 rank 为同一面分配同一 tag、边 ghost 填充覆盖 D3Q19 对角方向所需边线）虽涉及 tag/ghost 等"细节"，但属于 **send/recv 必须匹配的接口契约**，**必须测试**；此前笼统的"不测 MPI tag 值"导致缺陷 ①② 潜伏（见「已知缺陷与治理路径」）。

**开发策略：严格 TDD**——测试顺序即开发顺序，每个组件先写测试（红），再实现（绿），再重构。每条测试只因一个原因失败。

**单元测试 Fixture 原则**：凡依赖文件 I/O 的组件（如 GeometryEngine、VTK Writer），单元测试一律使用硬编码 fixture 数据（triangle soup、BlockCollection），不依赖真实文件。真实文件路径仅在集成测试中出现。

**测试顺序（= 开发顺序，各层独立可测）**：

| 顺序 | 测试名 | 验证内容 | 层级 |
|---|---|---|---|
| 0 | `test_octree_forest` | OctreeForest refine/balance/partition 在单位立方体上的基本正确性；验证 `local_octants()`、`quadrant_bounds()`、`quadrant_level()` | Mesh/forest |
| 1 | `test_stl_reader` | ASCII/binary STL 解析：三角面片数、法向量方向、包围盒正确性 | Mesh/io |
| 2a | `test_geometry_adaptive_refine` | 硬编码 triangle soup + `OctreeForest`：表面附近叶层 `quadrant_level` 高于远场；可选 2-rank 集体 refine 不死锁 | Mesh/geometry |
| 2 | `test_geometry_engine` | 多场景 triangle soup（外流实心、内流方管、风洞 channel+obstacle 等）：`MaterialKind` 分布正确；非法洞壁 fixture 失败；不依赖大 STL 文件 | Mesh/geometry |
| 3 | `test_face_pair_list` | p8est_iterate face callback 正确识别同级面、粗细 hanging 面；**T09** 域外 `TreeBoundaryFaces` 枚举且不进 `SameLevelFaces` | Mesh/topology |
| 4 | `test_load_balancer` | 2^level 权重分区后各 rank 总权重差异 < 5% | Mesh |
| 5 | `test_ghost_topology` | FacePairList 中每条跨 rank 面对的 remote_rank 与 p4est ghost 层记录一致（本地校验）；精心细化模式确保 hanging 面和跨 rank 面必然出现 | Mesh |
| 6 | `test_block_collection` | `BlockCollection<T>` 增删查，以 `octant_id` 为键；`BlockIterator` 遍历正确 | Solver/field |
| 7 | `test_collision_bgk` | 单块 BGK 碰撞：质量守恒、动量守恒，收敛至 Maxwell 平衡态 | Solver/lbm |
| 8 | `test_ghost_schedule` | L1：`DummyBlock` 覆盖同 rank / 跨 rank / 空计划 / 混合面 / 角点 / `comm_tag`；L2：2-rank `BlockLattice` 邻接面 populations 一致 | Solver/field |
| 9 | `test_lagrava_coupler` | L1：`coupling_plan_` 条数/索引/象限；L2：同 rank 1 粗+1 细界面 ρ/u 连续、质量守恒；L3：2-rank 跨 rank macro 交换 | Solver/lbm |
| 10 | `test_time_loop_levels` | L1：L=3 一次 `advance_one()` 后各层 collide 计数 1:2:4、coupler 调用顺序；P1：真实 BlockLattice 不崩溃 | Solver/lbm |
| 11 | `test_vtk_writer` | 硬编码 `BlockCollection` + `DummyCellField`：输出合法 StructuredGrid `.vts`（CellData、`(N+1)³` 角点），空间范围与 `quadrant_bounds` 吻合；W3 后补 2-rank `.vtm`/`.pvd` P1 | Solver/io |
| 11a | `test_lattice_cell_kind` | `CellKind` 默认全 fluid；`kSolid` 跳过 `collide`；与既有 BGK 行为一致 | Solver/lbm |
| 11b | `test_domain_boundary` | no-slip tree ghost bounce-back；moving lid Zou-He 宏观速度与配置一致 | Solver/lbm |
| 11c | `test_bouzidi_link` | `BouzidiLinkData::Build` 缓存 `q_frac`；`stream()` Bouzidi 替代 solid pull；`initialize_from_material` 映射 T07 fixture | Solver/lbm |
| 12 | `test_cavity3d_serial` | 单 rank，cavity3d（**单根 octant、L=0、无 AMR**），OpenLB 默认参数，**iT=5269**（OpenLB 收敛步）；`x=z=0.5` 竖直中线 17 点 vs OpenLB 参考剖面，`u_x/u_lid` **相对 L2** < 2% | Integration |
| 13 | `test_cylinder3d_parallel` | 多 rank（≥4），cylinder3d L=4，阻力系数 Cd 与 OpenLB 参考值误差 < 1% | Integration |
| 14 | `test_amr_convergence` | L=1→3 逐级细化，速度场 L2 误差收敛阶 ≈ 2 | Integration |

**依赖链**：0 → 3 → 5 → 8 → 9 → 10；0 → 2a；1 → 2a；2a → 2；1 → 2；6 → 8；7 → 9；0–11 + **T09（11a–11c）** 全绿后进入 Integration。

**任务拆分（`doc/tasks/`，与测试序对照）**：

| 任务 | 内容 | 阻塞于 | 状态 | 对应测试 / 算例 |
|------|------|--------|------|-----------------|
| T07 | Mesh 几何链 | T01、T02 | **测试绿** | #1、#2a、#2 |
| T08 | `vtk_writer` | T03、T01 | **测试绿** | #11（见 `T08-vtk-writer.md`） |
| **T09-W1** | 域边界 BC + `CellKind` + `TimeLoop` 注入 | T02、T04、T06 | **测试绿** | #11a、#11b；`test_face_pair_list` 扩展（见 `T09-boundary-bc.md`） |
| **T09-W2** | Bouzidi + `MaterialField` 初始化 | T09-W1、T07 | **测试绿** | #11c |
| **T10-W1** | `unit_converter` | T04 | **测试绿** | `test_unit_converter`（见 `T10-cavity3d.md`） |
| **T10-W2** | cavity 组装 + 冒烟 | T10-W1、T09-W1、T06 | **测试绿** | `test_cavity3d_serial`（smoke） |
| **T10-W3** | OpenLB #12 + `examples/cavity3d` + VTK P1 | T10-W2、T08 W4 | **测试绿** | #12、`examples/cavity3d` |
| **T11** | `cylinder3d` + 进出口 BC | T09-W2、#12 | 进行中（P0/W1/W2-sanity 绿，③④ 落定；W3 待开始） | #13 |
| **T11-refactor** | BC 调度架构重构（per-cell `BcKind` + 中心化 dispatch） | T09-W1/W2、T11-W2-sanity | 未开始（2026-07-15 立项） | 交付 T11 W3 组件 2；详见 `T11-refactor-bc-dispatch.md` |
| **T12** | AMR 收敛 | T11 | 未开始 | #14 |

T08 与 T07 无硬阻塞；T09-W2 阻塞于 T07；**当前里程碑**：T01–T10 已绿（含集成 #12 `test_cavity3d_serial`）；下一项 **T11 cylinder3d**（#13）已拆为 P0 + W1–W4 分波（见 `doc/tasks/T11-cylinder3d.md`），先修已知缺陷 ①–⑤ 再验收 Cd<1%。

### T10 · cavity3d 实现决策（2026-06-08，验收 2026-06-24）

| 项目 | 决策 |
|------|------|
| **算例对齐** | 离散参数与 OpenLB `examples/laminar/cavity3d` 默认一致：`N=30`，`τ=0.509`，`L=1`，`U_lid=1`，`ν=0.001`（**Re=1000**），`ρ=1`，`T_max=100` phys → `iT_max=getLatticeTime(100)` = **30000** |
| **网格** | `OctreeForest` **单根 octant**、`[0,1]³`、**不 refine**（`max_level=0`）；多级 cavity 不在 T10 |
| **时间推进** | **固定步数**；#12 与 example 默认 **iT=5269**（OpenLB ValueTracer 收敛步）；T10 **不**实现收敛检测 |
| **unit_converter** | `octlb` 轻量 struct，公式对齐 OpenLB `UnitConverterFromResolutionAndRelaxationTime`；不复制 OpenLB 完整 `unitConverter.h` |
| **#12 验收** | 17 点竖直中线（Ghia 高度采样 `x=z=0.5`）；`u_x/u_lid` 相对 **OpenLB 参考剖面** L2 < 2%；Ghia Re=100 表仅诊断 |
| **域 BC** | 六面 `kInterpolatedVelocity`（OpenLB `InterpolatedVelocity` / Skordos FD）；顶面 `u_wall=(u_lid_lattice,0,0)`，其余面 0 |
| **ConstRho** | `ConstRhoBGK` + `ConstRhoStatsScope::kFluidAndBoundary`（OpenLB 对齐）；不做 scope A/B 集成对比 |
| **example** | `examples/cavity3d/` 与集成测试共用 `cavity3d_case.h`；默认写 VTK + `centerline.csv` |
| **VTK** | T08 `AmrVtkWriter` + `vtk_lbm_fields.h`；P1 可视化（过域心 x–y 平面，ParaView Slice）；ctest 默认不写 VTK |
| **阻塞** | 仅 **T09-W1**；不依赖 T09-W2 / STL / Bouzidi |

**测试基础设施**：GTest + MPI。
- `tests/mpi_main.cpp`：手写 `MPI_Init → RUN_ALL_TESTS → MPI_Finalize`，所有 MPI 测试 target 链接此 main 而非 `GTest::gtest_main`。
- CMake 用 `${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 --oversubscribe` 启动测试，跨平台（本地 `mpirun`、服务器 `srun` 均适用）。
- p4est 集成：`cmake/FindP4est.cmake` 暴露 `P4est::p4est` imported target（INTERFACE 依赖 `MPI::MPI_C`）；跨机器用 `-DP4EST_ROOT=...` 指定安装路径。

---

## 已知缺陷与治理路径（2026-07-10 推敲）

T01–T10 已绿，但审计已实现栈发现若干潜伏缺陷——**根因是没有多块 / 多 rank / 多层端到端测试**：cavity3d（T10）是单根 octant、L=0、无 AMR（不触发跨 rank、粗细面、块边 ghost），单元测试用 `DummyBlock` 夹具。cylinder3d（T11，多 rank、L=4）是首个能同时暴露这些缺陷的真实算例。

| # | 缺陷 | 置信度 | 位置 |
|---|------|--------|------|
| ① | `comm_tag` 跨 rank 不对称（per-rank 冲突集不同→rehash 路径分叉→同面两 rank 不同 tag→`Isend`/`Irecv` 不匹配）；且 rehash mask 未对齐 `MPI_TAG_UB`，可超限 | CONFIRMED 读码 | `mesh/topology/face_pair_list.cpp:121,129,137` |
| ② | `GhostSchedule` 只刷面层 N×N，D3Q19 的 12 条对角方向在块边/角读 12 条边 ghost 线，未被任何面交换填充 | 合同层 CONFIRMED；运行时未确认 | `solver/lbm/block_lattice.cpp:15-29`（`face_buffer_count=ny*nz*kQ`）+ 本 PRD「每面单层 N×N」契约。**2026-07-11 更新**：`FacePackable` 概念已扩展 `pack_edge`/`unpack_edge`，`BlockLattice` 实现 12 条边 ghost 线 pack/unpack，`GhostSchedule` 加同 rank 边交换（face 邻居组合 `face_nbr[face_nbr[A][d1]][d2]`）——「Stage A 同 rank 边交换」已落。但 W2 multi-vs-single 探针**未能确认 ② 运行时显现**：探针被体素化混洞污染（per-octant vs whole-grid 在圆柱表面的 solid/fluid 分类差异；在 `f=f_eq−t` 约定下 solid 与 fluid-at-rest 都是 f=0→rho=1,u=0，step0 不可见，step1 才发散），已禁用该探针。② 需干净重探（无圆柱 multi-vs-single，或 toggle 有/无边交换）。跨 rank 边交换（Stage B，需 p4est corner callback 取对角邻居 rank）留 W3。 |
| ③ | TimeLoop 阶段顺序：早期本文档为 `collide→exchange→domain_bc.apply→stream`，实际为 `collide→domain_bc.apply→exchange→stream→apply_post_stream`（InterpolatedVelocity/PlaneFd 需 stream 后邻居） | **已决策落定（2026-07-15）**——以已验证实现为基准，本文档已更新为两段式 BC + PostStream 阶段为正式 spec；不改变 seam 设计原则 | `solver/lbm/time_loop/time_loop.h:175-209` |
| ④ | `BlockLattice` 带 overlap-padding/PostStream/rotate 机制（`collide_overlap_padding_bgk` 等），早期本文档「最小自建 BlockLattice」未预见 | **已决策落定（2026-07-15）**——保留该机制（cavity3d ~2% 校准所需），「最小自建」立论精确化为「自建存储+不引框架层」，padding 属存储内部细节；硬不变量不变 | `solver/lbm/block_lattice.h` |
| ⑤ | `LevelCoupler` 继承 ①：`MpiBatch` 按 `(peer,tag)` 去重——tag 碰撞把两个不同粗细面并成一 batch→错误宏量 | PROPAGATES | `solver/lbm/level_coupler.cpp:128,246` |
| ⑥ | BC 调度架构：per-face `DomainBcSpec` + 整 handler `boundary_lattice_mode_` 单一布尔 → 同 block 内混合 BC（入口速度+出口压力+圆柱 Bouzidi+壁 bounce-back）无法表达；`UsesInterpolatedVelocity()` 门控使 `apply()` 把所有 `kBoundary` 格统一当 Dirichlet 速度（含圆柱表面）。每接入 BC 组合不同的新算例需在唯一 FD 路径加 `spec.type` 分支 | CONFIRMED 读码（`examples/cylinder3d/cylinder3d_case.h:24` 注释自证） | `solver/lbm/domain_boundary_handler.h:42-49`、`domain_boundary_handler.cpp:164-186` |

**治理路径（T11 分波，oracle 三阶递进）**：

- **P0 前置**（修 ①）：`comm_tag` 改 per rank-pair 确定性枚举（按 canonical `face_key` 排序，两侧算出相同唯一 tag）根除碰撞类；rehash mask 对齐 `MPI_TAG_UB`；加跨 rank 对称性单测（强制碰撞→断言两侧同 tag）。同时消除 ⑤ 的传播根源。
- **W1**：阻力/Cd 后处理器（momentum-exchange，全新）+ 出口出流 BC（全新）组件单测；备圆柱 STL fixture。
- **W2 uniform cylinder3d**（多 rank，oracle=sanity：Cd 有限/符号正确/质量守恒/无 NaN）：**已完成**——`test_cylinder3d_uniform`（RunsNoCrash / MassBoundedDrift / DragFinitePositiveSign）全绿；legacy BC（Zou-He 入口 + kOutflow 出口 + kNoSlip 壁 + Bouzidi 圆柱）绕开 ④padding 与 InterpolatedVelocity/圆柱 kBoundary 冲突。NoEdgeArtifact（② 探针）因 multi-vs-single 体素化混洞**禁用**（见缺陷 ②）。**③④ HITL 已落定（2026-07-15）**：③ 以已验证实现（两段式 BC + PostStream）为基准更新 PRD；④ 保留 overlap-padding、PRD「最小自建」精确化为「自建存储+不引框架层」。W2 仍用 legacy BC 路径跑 sanity，InterpolatedVelocity/padding 路径在 W3/W4 切入。
- **W3 前置 · BC 调度架构重构（治 ⑥，T11-refactor）**：移植 OpenLB「per-cell + material 映射」架构思想为 OctLB「per-cell `BcKind` 枚举 + 中心化 `BcDispatcher`」（不引框架层、无虚函数；枚举失效再切多态兜底）；`MaterialKind` 保持几何 only，`BcInstaller` 求解器侧 stamp；移除 `boundary_lattice_mode_` 全局门控，混合 BC 自由共存。分 R0–R4 五波 TDD 推进，**每波保持 T01–T10 既有 ctest 全绿**（#12 L2<2% 不退化为硬门）。压力出口 `kPressureDirichlet`（T11 W3 组件 2）作为首个新客户在 R2 落位。详见 `doc/tasks/T11-refactor-bc-dispatch.md`。
- **W3 AMR L=4 + Lagrava**（oracle=量级：Cd 同量级/2x 内）：加几何 AMR + `LevelCoupler` 粗细耦合 + 跨 rank 粗细面（⑤ 路径首次真实流运行）；② Stage A 同 rank 边交换已落，**Stage B 跨 rank 边交换**（p4est corner callback 取对角邻居 rank + 跨 rank 边 MPI）在本波补；② 干净重探（无圆柱 multi-vs-single 或 toggle 对比）也在此波。
- **W4 Cd<1% 严格门**：参数对齐 OpenLB cylinder3d 参考、长时稳态、Cd<1%，验收 #13。

详见 `doc/tasks/T11-cylinder3d.md`。

## 不在范围内（第一版）

- **动态 AMR**：运行时 refine/coarsen 及字段迁移（`p4est_transfer_fixed`）。预留 `FacePairList::rebuild()` 及 Solver 侧 `GhostSchedule` / `LevelCoupler` / `DomainBoundaryHandler` / `BouzidiLinkData` / `TimeLoop` 重建契约，第一版不实现热路径。
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
| Bouzidi 曲面 BC 逻辑 | `solver/lbm/boundary/`、`bouzidi_link_data.*` | **T09-W2 已实现**（参考 OpenLB / octree-mesh） |
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
> **当前进度（2026-07-15）**：T01–T10 已绿（单元 + 集成 #12 `test_cavity3d_serial`）；`examples/cavity3d` 可跑；审计发现已实现栈缺陷 ①–⑤（见「已知缺陷与治理路径」），**T11** cylinder3d（#13）P0/W1/W2-sanity 已绿、③④ HITL 已落定、W3 待开始（见 `doc/tasks/T11-cylinder3d.md`）；W3 推进时发现 BC 调度架构缺陷 ⑥，立项 **T11-refactor**（per-cell `BcKind` + 中心化 dispatch，见 `doc/tasks/T11-refactor-bc-dispatch.md`）作为 W3 前置，压力出口作首个新客户；T12 AMR 收敛未开始（依赖 T11）。

### Mesh 模块

| 组件 | 所在路径 | 状态 | 对应测试 | 备注 |
|---|---|---|---|---|
| OctreeForest | `mesh/forest/` | 测试绿 | `test_octree_forest` | T01；CI main 已通过 |
| stl_reader | `mesh/io/stl_reader/` | 测试绿 | `test_stl_reader` | T07 W1；移植自 OpenLB |
| GeometryEngine + MaterialField | `mesh/geometry/` | 测试绿 | `test_geometry_engine`、`test_geometry_adaptive_refine` | T07；自适应加密 + Part 合并 + CGAL 体素化 |
| FacePairList | `mesh/topology/` | 测试绿 | `test_face_pair_list` | T02；T05 `comm_tag`；T06 `comm_tags[4]`；**T09-W1** `TreeBoundaryFaces` 枚举 |
| WeightedLoadBalancer | `mesh/load_balance/` | 测试绿 | `test_load_balancer` | T02；本地 4-rank 通过，待 CI |
| Ghost 拓扑一致性 | `mesh/topology/` | 测试绿 | `test_ghost_topology` | T02；本地 4-rank 通过，待 CI |

### Solver/field 模块

| 组件 | 所在路径 | 状态 | 对应测试 | 备注 |
|---|---|---|---|---|
| BlockCollection\<T\> + BlockIterator | `solver/field/` | 测试绿 | `test_block_collection` | 泛型，不依赖 LBM 头文件 |
| GhostSchedule\<T\> + FaceIterator | `solver/field/` | 测试绿 | `test_ghost_schedule` | T05；`octlb_field_schedule`；同级面层 MPI |
| FacePairList `comm_tag` | `mesh/topology/` | 测试绿 | `test_face_pair_list`（扩展） | T05 扩展 T02 |
| BlockLattice `pack_face` / `unpack_face` | `solver/lbm/` | 测试绿 | `test_ghost_schedule` L2 | T05 |
| `CellKind` + `populations_at_halo` | `solver/lbm/cell_kind.h`、`block_lattice.*` | 测试绿 | `test_lattice_cell_kind` | T09-W1；solid/boundary 跳过 collide/stream |

### Solver/lbm 模块

| 组件 | 所在路径 | 状态 | 对应测试 | 备注 |
|---|---|---|---|---|
| `octlb::BlockLattice` + BGK | `solver/lbm/block_lattice.*` | 测试绿 | `test_collision_bgk` | T04；T09-W2 扩展 Bouzidi pull |
| MRT / Smagorinsky Dynamics | `solver/lbm/dynamics/` | 未开始 | — | 复用 OpenLB 算法头文件，集成测试间接覆盖 |
| 域边界 BC（`DomainBoundaryHandler`） | `solver/lbm/domain_boundary_handler.*`、`boundary/bounce_back.h`、`boundary/zou_he_velocity.h`、`boundary/interpolated_velocity.h` | 测试绿 | `test_domain_boundary`、`test_interpolated_velocity_compute_rho` | T09-W1 no-slip + moving lid；**T10** 六面 `InterpolatedVelocity` + corner `ComputeRho` |
| Bouzidi BC + `MaterialField` 初始化 | `solver/lbm/boundary/bouzidi_pull.h`、`bouzidi_link_data.*`、`lattice_material_init.*` | 测试绿 | `test_bouzidi_link` | T09-W2；init 缓存 `q_frac` + `stream()` pull；静态 AMR v1 |
| LevelCoupler（Lagrava） | `solver/lbm/level_coupler.*` | 测试绿 | `test_lagrava_coupler`、`test_lagrava_coupler_mpi2` | T06；算法复用 dynamics/，Solver 侧 coupling_plan_ |
| TimeLoop（递归下降） | `solver/lbm/time_loop/` | 测试绿 | `test_time_loop_levels` | T06；T09-W1 注入 `DomainBoundaryHandler`（`NoOp` 回归）；v1 每步全局 `GhostSchedule::exchange()` |
| CoarseFineFace `comm_tags[4]` | `mesh/topology/` | 测试绿 | `test_face_pair_list`、`test_lagrava_coupler_mpi2` | T06 扩展 T02 |
| unit_converter | `solver/lbm/unit_converter/` | 测试绿 | `test_unit_converter`（T10-W1） | OpenLB `FromResolutionAndRelaxationTime` 公式；见 `T10-cavity3d.md` |
| cavity3d 集成 | `examples/cavity3d/`、`tests/integration/` | 测试绿 | `test_cavity3d_serial`（#12） | T10-W2/W3；单根 L=0、iT=5269、OpenLB 参考剖面相对 L2 |

### Solver/io 模块

| 组件 | 所在路径 | 状态 | 对应测试 | 备注 |
|---|---|---|---|---|
| vtk_writer | `solver/io/vtk_writer/` | 测试绿 | `test_vtk_writer`、`test_vtk_writer_two_rank` | T08；StructuredGrid `.vts` + CellData + `VtkCellField3D` + 并行 `.vtm`/`.pvd`；见 `doc/tasks/T08-vtk-writer.md` |
| gnuplot | `solver/io/gnuplot/` | 未开始 | — | 直接复制自 OpenLB |

### 集成测试

| 测试 | 状态 | 前置条件 | 验收标准 |
|---|---|---|---|
| `test_cavity3d_serial` | 测试绿 | 单元测试 0–11 + T09（11a–11c）全绿 | OpenLB 默认参数、**iT=5269**；17 点 `u_x/u_lid` 相对 OpenLB 参考 L2 < 2%（**T10-W3**）；见 `T10-cavity3d.md` |
| `test_cylinder3d_parallel` | 进行中（P0–W4 分波） | `test_cavity3d_serial` + **T09-W2** | ≥4 rank，Cd 与 OpenLB 参考值误差 < 1%（**T11**）；见 `doc/tasks/T11-cylinder3d.md` |
| `test_amr_convergence` | 未开始 | `test_cylinder3d_parallel` 通过 | L=1→3 速度场 L2 误差收敛阶 ≈ 2（**T12**） |
