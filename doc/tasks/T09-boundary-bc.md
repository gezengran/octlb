# T09 · 域边界 BC + Bouzidi 曲面 BC（W1 / W2）

> 类型：AFK  
> 阻塞于：T02（`FacePairList`）、T04（`BlockLattice`）、T06（`TimeLoop`）、T07（`MaterialField`，**仅 W2**）  
> 对应 PRD：`octlb/doc/prd/octlb-framework.md`（Solver/lbm + Mesh 拓扑扩展；用户故事 #5 边界算子、#9 cylinder 验证前置）  
> 状态：**W1+W2 已完成**（静态 AMR v1）

---

## 要做什么

在 T08 完成后，补齐 **LBM 边界条件** 的两条纵向切片。**第一版仅静态 AMR**（初始化一次；`FacePairList::rebuild()` 后由调用方整体替换 Handler / `TimeLoop`，不实现动态热路径）。

### T09-W1 · 域边界 + 格子材料骨架

1. **`TreeBoundaryFace`（Mesh 扩展）**  
   在 `FacePairList` face callback 中，对 `tree_boundary != 0` 的面收集 `{octant_id, face_dir}` 列表（**不**进入 `SameLevelFaces` / `CoarseFineFaces`）。BC **类型**由 example 配置，Mesh 只枚举拓扑。

2. **`CellKind` + `BlockLattice` 轻量预埋**  
   每 interior 格存 `kFluid` / `kSolid` / `kBoundary`（默认 `kFluid`）。`collide()` / `stream()` **跳过** `kSolid` 与 `kBoundary`。W1 **不**读取 `MaterialField`；W2 再灌 tag。

3. **`DomainBoundaryHandler`**  
   构造时绑定 `TreeBoundaryFaces` + BC 配置（如 `FaceDir → no-slip bounce-back` / `moving_lid` Zou-He）。`apply()` 在 post-collision  populations 上写入 **tree 外侧 ghost 层**（`collide` 之后、`stream` 之前）。

4. **`TimeLoop` 注入**  
   构造参数增加 `DomainBoundaryHandler&`（可 `NoOp`）；每层顺序：  
   `collide → GhostSchedule::exchange() → domain_bc.apply() → stream → … coupler …`  
   外部引用、不持有（对齐 T06 `GhostSchedule` / `LevelCoupler`）。

5. **OpenLB 算法头（W1 最小集）**  
   复制/适配 `bounceBack.h`、`zouHeVelocity.h`（或等价无状态模板）至 `solver/lbm/boundary/`；**不**引入 OpenLB 框架层。

W1 **不**实现：Bouzidi、`MaterialField` 绑定、cylinder 进出口 BC（**T11**）。

### T09-W2 · Bouzidi + MaterialField 初始化

1. **`initialize_from_material()`**  
   从 T07 `MaterialField` 写入每格 `CellKind`（`fluid` / `solid` / `boundary` 映射）；fluid 格设统一 `omega`。

2. **`BouzidiLinkData`（初始化缓存）**  
   Solver 侧一次性构建：对每个 **相邻 solid/boundary 的 `kFluid` 格**、每个 cut link 存 `q_frac`（沿 lattice 方向到壁面的分数距离）。输入：`MaterialField` + `GeometryAssembly`（三角面）+ `OctreeForest::quadrant_bounds` + 格心物理坐标。**不**回写 `MaterialField`（符合 T07「不预存 Bouzidi 距离/法向」）。参考 octree-mesh `ComputeWallDistances` + `BouzidiPostCollisionPull`。

3. **集成进 `BlockLattice::stream()`**  
   pull 阶段：若 stream 源邻居为 `kSolid` / `kBoundary`，用 post-collision `f` 与缓存 `q_frac` 调用 Bouzidi 公式，**不**从 solid 格读 population。域外 ghost 仍由 W1 `DomainBoundaryHandler` 负责；**不**再注入 `TimeLoop`。

4. **OpenLB 算法头（W2）**  
   补充 Bouzidi 无状态核（可复制 `bouzidiFields.h` 思路或自写 `BouzidiPostCollisionPull`，对齐 octree-mesh `bounce_back.h`）。

W2 **不**实现：动态 AMR 下 `BouzidiLinkData` 热重建、MRT/Smagorinsky、进出口 Poiseuille（**T11**）。

---

## 交付物

```
octlb/
├── src/
│   ├── mesh/topology/
│   │   ├── face_pair_list.h       # TreeBoundaryFace + tree_boundary_faces()
│   │   └── face_pair_list.cpp
│   └── solver/lbm/
│       ├── cell_kind.h            # CellKind enum
│       ├── block_lattice.h/.cpp   # CellKind 存储；solid/boundary 跳过；W2 stream Bouzidi
│       ├── boundary/
│       │   ├── bounce_back.h      # W1：no-slip
│       │   ├── zou_he_velocity.h  # W1：moving lid
│       │   └── bouzidi_pull.h     # W2：BouzidiPostCollisionPull
│       ├── domain_boundary_handler.h/.cpp   # W1
│       ├── bouzidi_link_data.h/.cpp         # W2：init 缓存 q_frac
│       ├── lattice_material_init.h/.cpp     # W2：MaterialField → CellKind
│       └── time_loop/time_loop.h  # W1：注入 DomainBoundaryHandler
└── tests/unit/
    ├── mesh/test_face_pair_list.cpp         # 扩展：TreeBoundary 枚举
    └── solver/
        ├── test_lattice_cell_kind.cpp       # W1：solid 跳过 collide
        ├── test_domain_boundary.cpp         # W1：tree ghost BC
        └── test_bouzidi_link.cpp            # W2：q_frac + stream 行为
```

---

## 关键设计决策（已对齐）

| 项目 | 决策 |
|------|------|
| **任务粒度** | 单一 T09 文档；**W1 → W2** 两波纵向切片，可分 PR merge |
| **域 BC 拓扑** | Mesh 枚举 `TreeBoundaryFaces`（**方案 B**）；example 配置 BC 类型 |
| **域 BC 调用** | `DomainBoundaryHandler` **注入 `TimeLoop`**（**方案 A**）；`rebuild()` 后整体替换 |
| **每步顺序（域 BC）** | `collide → exchange → domain_bc → stream` |
| **CellKind（W1）** | 轻量预埋，默认 `kFluid`；**不**读 `MaterialField`（**方案 B**） |
| **Bouzidi（W2）** | init 缓存 `BouzidiLinkData` + **集成 `stream()` pull**（**方案 A**）；静态 AMR v1 只算一次 |
| **Bouzidi 作用格** | **`kFluid` 格**指向 solid/boundary 邻居的 cut link；`kSolid`/`kBoundary` 不参与 collide/stream |
| **与 GhostSchedule 分工** | 同级 halo → `GhostSchedule`；域外 tree ghost → W1 handler；STL 曲面 → W2 stream 内 Bouzidi |
| **动态 AMR** | v1 不实现热路径；预留 `rebuild()` 后重建 Handler / `BouzidiLinkData` / `TimeLoop` 契约 |

---

## 波次与纵向切片

| 波次 | 标题 | 类型 | 阻塞于 | 用户故事 | 端到端交付 |
|------|------|------|--------|----------|------------|
| **W1** | 域边界 BC + CellKind + TimeLoop | AFK | T02、T04、T06 | —（cavity 前置） | `TreeBoundaryFaces` + moving lid / no-slip 单块测试绿 |
| **W2** | Bouzidi + MaterialField 初始化 | AFK | W1、T07 | #5、#9 前置 | T07 fixture + Bouzidi stream 单测；不挡 W1 merge |

---

## 测试决策

**原则**：只测对外行为；W1 不依赖 STL / `GeometryEngine`；W2 用 T07 hardcoded fixture + 极简 triangle soup。

### W1 · `test_face_pair_list`（扩展）

| 用例 | 断言 |
|------|------|
| `DomainTreeBoundary_NoSameLevelFace` | 域界面 **不**出现在 `SameLevelFaces` |
| `TreeBoundaryFaces_Enumerated` | 单位立方体单 octant 时 `tree_boundary_faces()` 含 6 条 `(id, dir)` |

### W1 · `test_lattice_cell_kind`

| 用例 | 断言 |
|------|------|
| `SolidCell_SkipsCollide` | 硬编码部分 `kSolid`：collide 前后 population 不变 |
| `DefaultAllFluid` | 未设 tag 时与现网 `collide` 行为一致 |

### W1 · `test_domain_boundary`

| 用例 | 设置 | 断言 |
|------|------|------|
| `NoSlip_FillsGhostAfterCollide` | 单块 + 一面 no-slip tree BC | ghost 层 post-collision 满足 bounce-back 关系 |
| `MovingLid_ZouHe` | 顶面 Zou-He，已知 `u_wall` | ghost 宏观速度与配置一致（容差 ε） |
| **P1** `TimeLoop_WithNoOpDomainBc` | `NoOpDomainBoundaryHandler` | `test_time_loop_levels` 计数与改前一致 |

### W2 · `test_bouzidi_link`

| 用例 | 设置 | 断言 |
|------|------|------|
| `LinkData_CutLinkQFrac` | 硬编码 fluid+solid 邻接 + 三角面 | 指定方向的 `q_frac ∈ (0,1)` |
| `Stream_BouzidiReplacesSolidPull` | 单块 stream 一步 | fluid 格从 solid 方向 pull 使用 Bouzidi 而非 solid population |
| **P1** `MaterialField_MapsCellKind` | T07 S1 fixture | `initialize_from_material` 后 `CellKind` 与 `MaterialKind` 一致 |

### 运行配置

- W1 默认 **1 rank**；扩展 face 测试可 1 rank  
- W2 P1 可 1 rank；跨 rank 曲面 BC 留 **T11** 集成验收  
- 复用 `tests/mpi_main.cpp`

---

## 验收标准

### W1

- [x] `cmake --build build -j4` 通过
- [x] `TreeBoundaryFaces` 与 T05「域面不进 `SameLevelFaces`」一致
- [x] `CellKind` 默认全 fluid；solid/boundary 跳过 collide/stream
- [x] `DomainBoundaryHandler` + `TimeLoop` 注入；`NoOp` 不破坏 T06 测试
- [x] `test_domain_boundary` / `test_lattice_cell_kind` / face 扩展 **全绿**
- [x] T01–T08 相关 ctest **仍绿**

### W2

- [x] `initialize_from_material()` 映射 T07 `MaterialField`
- [x] `BouzidiLinkData` 初始化缓存；**无**每步 runtime 三角查询
- [x] `stream()` 内 Bouzidi 替代 solid 方向 pull
- [x] `test_bouzidi_link` **全绿**
- [x] **不**修改 `TimeLoop` 签名（相对 W1 终态）

---

## 阻塞关系

```
T02（FacePairList）— 已完成
T04（BlockLattice）— 已完成
T06（TimeLoop）— 已完成
T07（MaterialField）— 已完成（W2 需要）
└── T09
    ├── W1 · 域边界 + CellKind + TimeLoop
    │       └── T10 · cavity3d + unit_converter（#12）
    └── W2 · Bouzidi + MaterialField 初始化
            └── T11 · cylinder3d（#13）→ T12 · AMR 收敛（#14）
```

---

## 构建与测试命令（给用户执行）

```bash
cd octlb
module load octlb
cmake -B build -DCMAKE_BUILD_TYPE=Debug \
      -DP4EST_ROOT=/home/dinglin/Public/p4est-2.8/opt
cmake --build build -j4

# W1
cd build && ctest --output-on-failure -R "test_domain_boundary|test_lattice_cell_kind|test_face_pair_list"

# W2
cd build && ctest --output-on-failure -R "test_bouzidi_link"
```
