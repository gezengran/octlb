# T06 · LevelCoupler（Lagrava）+ TimeLoop（递归下降）

> 类型：AFK  
> 阻塞于：T02（FacePairList）、T04（BlockLattice + BGK）、T05（GhostSchedule、FaceIterator）  
> 对应 PRD：`octlb/doc/prd/octlb-framework.md`（Solver/lbm 子层，测试顺序 #9、#10）  
> 状态：**已完成**（2026-06-04；本地 `ctest` 10/10 通过）

---

## 要做什么

在 `solver/lbm/` 实现 **粗细界面 Lagrava 守恒耦合** 与 **递归下降 TimeLoop**，并扩展 Mesh 侧 `CoarseFineFace` 的 MPI tag（与 T05 `SameLevelFace::comm_tag` 同模式）：

1. **`CoarseFineFace::comm_tags[4]`（T02 小扩展）**  
   在 `FacePairList` face callback 内，为每条粗细面、每个 `fine_ids[i]` 与 coarse 侧生成**对称** MPI tag（`remote_ranks[i]` 跨 rank 时使用；同 rank 忽略）。更新 `test_face_pair_list` / `test_ghost_topology` 中粗细面相关断言（若适用）。

2. **`LevelCoupler`**  
   - **算法复用**（同 T04）：按 PRD Lagrava 公式实现 `apply_half_time` / `apply_full_time` / `restrict`；宏观量与平衡态复用已有 `olb::lbm::*` + `CellProxy`。**不**复制 OpenLB `refinement/` 框架头文件。  
   - 构造时从 `FaceIterator` + `OctreeForest::quadrant_bounds()` + 块尺寸 `(N,N,N)` 固化 **`coupling_plan_`**（Solver 侧 cell 级连接：`coarse_id, ci,cj,ck` ↔ `fine_id, fi,fj,fk`，`normal`）。Mesh 保持 block 级纯拓扑。  
   - 跨 rank：构造时固化 MPI 计划（coarse 宏观量 ρ、u、f_neq → fine 侧）；`exchange()` 模式对齐 `GhostSchedule`（Irecv → Isend/本地直写 → Waitall）。

3. **`TimeLoop`**  
   - 类 + **外部引用**注入 `GhostSchedule<BlockLattice>`、`LevelCoupler`（调用方持有，便于 `rebuild()` 后整体替换）。  
   - 构造时缓存 `level → [OctantId]`；`advance_one()` 实现 PRD 递归下降（collide → halo → stream → half-time prolong → 细子步 ×2 → full-time prolong → restrict）。  
   - 测试用 collide/stream **调用计数 hook**（或等价可观测接口），供 `test_time_loop_levels` 验证 1:2:4 步数比。

本任务**不**实现 Bouzidi 域边界（T07）、VTK / STL（T07）、`unit_converter`（留到 `cavity3d` 集成）、动态 AMR 热路径（仅预留 rebuild 契约）。

---

## 交付物

```
octlb/
├── src/
│   ├── mesh/topology/
│   │   ├── face_pair_list.h      # CoarseFineFace 增加 comm_tags[4]
│   │   └── face_pair_list.cpp    # 对称 comm_tag 生成（per fine slot）
│   └── solver/lbm/
│       ├── level_coupler.h       # LevelCoupler + coupling plan 类型
│       ├── level_coupler.cpp     # Lagrava 公式 + MPI 交换
│       ├── time_loop/
│       │   └── time_loop.h       # TimeLoop（可 header-only 或 .cpp）
│       └── CMakeLists.txt        # octlb_lbm 追加源文件
└── tests/unit/solver/
    ├── test_lagrava_coupler.cpp  # L1 plan → L2 同 rank 物理 → L3 跨 rank
    ├── test_time_loop_levels.cpp # L1 步数 1:2:4；P1 真实 lattice
    └── CMakeLists.txt            # 2-rank MPI 用例
```

---

## 关键设计决策（已对齐）

| 项目 | 决策 |
|---|---|
| **任务粒度** | 单一 T06：Lagrava + LevelCoupler + TimeLoop 一并交付；测试做充分 |
| **Lagrava 落地** | **算法复用 + 自建数据结构**（T04 同路线）；`lagrava.h` 仅作公式对照，**不**引入 OpenLB `BlockRefinementContext` / `refinement/` 框架 |
| **连接点索引** | **Solver 侧 `coupling_plan_`**：构造时由 `FaceIterator` + `quadrant_bounds` + N 计算；Mesh 不存 cell 级 coupling |
| **动态 AMR 兼容** | `FacePairList::rebuild()` 后须**销毁并重建** `LevelCoupler`、`GhostSchedule`、`TimeLoop`（及 level 缓存）；第一版静态 AMR 不实现热路径，T06 写 P1 rebuild 测试 |
| **跨 rank MPI** | `CoarseFineFace::comm_tags[4]` 在 Mesh callback 对称生成；`LevelCoupler` 构造时固化 MPI 计划（对齐 T05 `GhostSchedule`） |
| **TimeLoop 所有权** | `TimeLoop` **不**拥有 `GhostSchedule` / `LevelCoupler`；构造参数为引用，与 Schedule 独立可测一致 |
| **τ / ω** | 第一版全域统一 `omega`（BGK）；`τ = 1/omega` 用于 Lagrava scalingFactor |
| **与 GhostSchedule 分工** | 同级 1:1 面仍仅 `GhostSchedule`；粗细 1:4 **不能**用 halo slab 替代 |
| **unit_converter** | 不在 T06；`cavity3d` 集成阶段再引入 |

---

## Lagrava 公式（LevelCoupler 行为契约）

**粗→细（prolongation，`apply_half_time` / `apply_full_time`）**

- `scalingFactor = (τ - 0.25) / τ`（`τ = 1/omega`）
- half-time：粗侧宏观量取 **prev 与 curr 的算术平均**（需在上层 coarse collide 前/后缓存 prev，或在 coupler 内按 OpenLB Lagrava 语义维护）
- full-time：使用 curr 宏观量
- `f_fine[iPop] = f_eq(ρ, u) + scalingFactor × f_neq[iPop]`

**细→粗（restriction，`restrict`）**

- `scalingFactor = τ / (τ - 0.25)`
- 对体积 2:1 子格（界面层参与子集）取均值得 `(ρ, u)`，再经 scaling 写回粗格 `f`

宏观量计算：`olb::lbm::computeRhoU` / `computeFneq` + `olb::equilibrium`（经 `CellProxy`）。

---

## coupling_plan_ 与 CoarseFineFace 扩展

```cpp
// face_pair_list.h — T06 扩展
struct CoarseFineFace {
  OctantId coarse_id;
  OctantId fine_ids[4];
  FaceDir normal;
  int remote_ranks[4];
  int comm_tags[4];   // Symmetric MPI tag per fine slot (T06).
};

// level_coupler.h — Solver 侧（示意）
struct CouplingPoint {
  OctantId coarse_id;
  int ci, cj, ck;
  OctantId fine_id;
  int fi, fj, fk;
  FaceDir normal;
  int remote_rank;    // coarse or fine owner for MPI leg
  int comm_tag;       // from CoarseFineFace::comm_tags[i] when cross-rank
};

class LevelCoupler {
 public:
  LevelCoupler(MPI_Comm comm,
               const FacePairList& faces,
               const OctreeForest& forest,
               BlockCollection<BlockLattice>& blocks,
               int nx, int ny, int nz,
               double omega);

  void apply_half_time(int coarse_level);
  void apply_full_time(int coarse_level);
  void restrict(int fine_level);

 private:
  // coupling_plan_、macro send/recv bufs、prev macro cache — 构造时分配
};
```

**Plan 构建**：对每个 `CoarseFineFace`，用 coarse/fine 的 `quadrant_bounds` 确定 4 个 fine block 在 coarse 面上的切向象限，映射到块内 `(i,j,k)` 界面层 cell；条数与拓扑 + N 确定性相关。

---

## TimeLoop

```cpp
class TimeLoop {
 public:
  TimeLoop(const OctreeForest& forest,
           BlockCollection<BlockLattice>& blocks,
           GhostSchedule<BlockLattice>& ghosts,
           LevelCoupler& coupler,
           double omega);

  void advance_one();

  int max_level() const;
  // 测试：各层 collide/stream 调用计数（或 injectable hook）
};
```

**每 coarse 步顺序**（与 PRD 一致）：

```
advance(level l):
  collide(blocks at l)
  ghosts.exchange()          // 仅同级 SameLevelFaces；按层块子集或全局 exchange（实现时二选一并文档化）
  stream(blocks at l)
  coupler.apply_half_time(l)
  advance(l+1)               // 第一个细子步
  coupler.apply_full_time(l)
  advance(l+1)               // 第二个细子步
  coupler.restrict(l+1)
```

第一版入口：`advance_one()` 从 level 0 递归；静态分层，初始化时确定 `max_level`。

---

## 拓扑 rebuild 契约（P1 测试 + 文档）

| 事件 | 动作 |
|---|---|
| `FacePairList::rebuild(forest)` | 重建 `GhostSchedule`、`LevelCoupler`；重建或刷新 `TimeLoop` 的 level 缓存 |
| 第一版运行期 | 不调用；初始化建一次 |
| 未来动态 AMR | `p8est_transfer_fixed` 完成后、下一步 `advance_one()` 前执行上表 |

---

## CMake

- `octlb_lbm` STATIC 追加 `level_coupler.cpp`（及 `time_loop` 若有 `.cpp`）
- `test_lagrava_coupler` / `test_time_loop_levels` → `octlb_lbm` + `octlb_field_schedule` + `octlb_mesh` + MPI
- **不**新增 `refinement/` 目录；**不**复制 OpenLB `src/refinement/` 头文件

---

## 测试决策

**原则**：只测对外行为；`test_lagrava_coupler` 与 `test_time_loop_levels` 分 target，职责分离。

### `test_lagrava_coupler`（PRD #9）

| 层级 | 用例 | 设置 | 断言 |
|---|---|---|---|
| **L1** | `CouplingPlan_MatchesTopology` | T02 中心细化 fixture，1 rank | plan 条数 > 0；无重复 `(coarse_id,ci,cj,ck,fine_id,fi,fj,fk)`；象限与 `fine_ids[4]` 一致 |
| **L1** | `CouplingPlan_InterfaceIndicesInRange` | 同上 | 所有界面索引 ∈ `[0,N)`，法向侧为界面第一层 |
| **L2** | `OneRank_CoarseFine_InterfaceContinuous` | 1 粗 + 1 细（或 1 粗 + 4 细）同 rank，均匀初场 + 若干手动 coupler 步 | 界面邻域 ρ、\|u\| 连续（阈值文档化，如相对误差 < 1e-6 均匀场） |
| **L2** | `OneRank_ProlongRestrict_MassConserved` | 同上，仅 fluid，无 BC | 单步 prolong+restrict 后总质量变化 < 1e-10 |
| **L3** | `TwoRank_CoarseFine_MacroExchange` | 2 rank，coarse/fine 分 rank，可识别 macro 模式 | 跨 rank prolong 后 fine 侧界面 macro 与 coarse 发送一致 |

### `test_time_loop_levels`（PRD #10）

| 层级 | 用例 | 设置 | 断言 |
|---|---|---|---|
| **L1** | `ThreeLevels_StepCountRatio_1_2_4` | L=0,1,2 静态 forest（如逐层中心细化），hook 计数 | 一次 `advance_one()` 后 L0:L1:L2 的 `collide`（或等价层步）次数 = 1:2:4 |
| **L1** | `ThreeLevels_CouplerCallOrder` | 同上，记录 coupler 调用序列 | half → 细子步 → full → 细子步 → restrict 顺序正确 |
| **P1** | `AdvanceOne_WithBlockLattice_NoThrow` | 真实 `BlockLattice`，小 N，少量 `advance_one()` | 不崩溃；质量无 NaN |

### P1（不挡 T06 合并）

- `Rebuild_RecreateCoupler`：`pairs.rebuild()` 后新建 `LevelCoupler`，plan 条数与拓扑一致
- `CoarseFine_CommTagsUnique`：多 rank 跨 rank 粗细面 `comm_tags[i]` 在本 rank 有效 send 条目内两两不同

### 运行配置

- L1 plan / 同 rank 物理：1 rank
- L3 / comm_tag：2 rank `--oversubscribe`
- 复用 `tests/mpi_main.cpp`

---

## 验收标准

- [x] `cmake --build build -j4` 通过
- [x] T02 mesh 测试仍绿（`comm_tags[4]` 扩展后）
- [x] `test_lagrava_coupler` L1 **全部**通过
- [x] `test_lagrava_coupler` L2 **至少 1 例**通过
- [x] `test_lagrava_coupler` L3 **至少 1 例**通过（`test_lagrava_coupler_mpi2`，2 rank）
- [x] `test_time_loop_levels` L1 **全部**通过
- [x] `octlb_lbm` **不**链接 OpenLB `refinement/` 或新增 `refinement/` 目录
- [x] `LevelCoupler` / `TimeLoop` **不** include p4est 头文件（经 `OctreeForest` 公共 API 取 bounds/level）
- [x] `field/` 子层仍无 LBM 依赖（粗细 MPI 在 `lbm/`）

---

## 阻塞关系

```
T02（FacePairList）— 本任务扩展 comm_tags[4]
T04（BlockLattice + BGK）— 已完成
T05（GhostSchedule + FaceIterator）— 已完成
└── T06（本任务）
    └── T07 · VTK Writer + STL Reader + Bouzidi BC（可与 T06 部分并行，集成依赖 T06）
```

---

## 构建与测试命令（给用户执行）

```bash
cd octlb
module load octlb   # 若需要
cmake -B build -DCMAKE_BUILD_TYPE=Debug \
      -DP4EST_ROOT=/home/dinglin/Public/p4est-2.8/opt
cmake --build build -j4

cd build && ctest --output-on-failure -R "test_face_pair_list|test_ghost_topology|test_lagrava_coupler|test_time_loop_levels"
```
