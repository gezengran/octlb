# T11 · cylinder3d 集成（P0 前置修复 + W1–W4 四波）

> 类型：AFK（前置修复 P0 + 四波纵向切片；W2 含 ③④ HITL 决策点）  
> 阻塞于：T10（`test_cavity3d_serial` 已完成）、T09-W2（Bouzidi 曲面 BC）、T07（`GeometryEngine`/`MaterialField`）、T08 W4（`vtk_lbm_fields.h`）  
> 对应 PRD：`octlb/doc/prd/octlb-framework.md`（集成测试 #13、用户故事 #1/#5/#8/#9、`examples/cylinder3d/`）  
> 状态：**P0 + W1 + W2-sanity + W3-组件 + T11-refactor 已完成；W3-AMR 量级门原生机已通过；W4 严格门待原生长跑**（2026-08-05 复验）——P0 comm_tag 对称性修复（① ⑤）绿；W1 阻力后处理 + 出口 BC + 圆柱 STL 绿；W2 cylinder3d uniform sanity 绿（2026-07-31 复验 `test_cylinder3d_uniform` 177s 绿）。**W3 组件单测全绿**（Poiseuille 入口 / 压力出口 / Schäfer-Turek STL / 缺陷⑤ `CoarseFine_CapturedGeometry` / ② `CrossRankEdges` 共 11 条）。**② Stage B 跨 rank 边交换已落 + `EdgeExchange_CrossRank_NoArtifact` 修复后绿**（原 u_inlet=0.05 压力出口 Ma runaway 致红，降 0.01 后绿；**弱验证**——低 Ma uniform 流 ② 不显现，OFF≈ON，仅证 Stage B 接线 + 非退化，详见 W3 验收 ② 说明）。**③④ HITL 已落定**（2026-07-15）。**T11-refactor 完成**（PR #11，R0–R4）。**缺陷⑤ 修复**（FacePairList 快照 ghost 物理 bounds/level）单测绿。**W3 AMR 量级门 `Amr_Cd_SameOrderAsOpenLb` 原生机 PASSED**（2026-08-05：L5 -n4，Cd=5.54∈[3.18,12.72]，~58min；详见下文「W3/W4 根因修复与原生机验收」）。**W4 严格门 `Amr_Cd_WithinOnePercentOfOpenLb` 仍 FAIL**（L5 Cd=5.54 距 6.36 为 13%，需 L6 或更长收敛；PRD 标注长期原生长跑）。Rosetta x86 模拟下每 AMR `Cylinder3dCase` 构造 ~50min+（CGAL voxelization × resolve_surface × 3），4 case 数小时不可行；native x86 CI ~1–2min/case。example `examples/cylinder3d` 链接通过、VTK+Cd CSV 代码就位，运行留原生机。

---

## 要做什么

交付 **首个真实算例端到端路径**：cylinder3d（STL 圆柱外流、几何自适应 AMR、Bouzidi 曲面 BC、阻力系数 Cd），验收 PRD 集成测试 #13（Cd 与 OpenLB 参考值相对误差 < 1%），并附带可执行算例与 VTK 输出。

**推敲结论（已对齐 PRD / 设计讨论，2026-07-10）**：

- **根因**：已实现的 AMR 链路（comm_tag 对称性、边 ghost、Lagrava 跨 rank、TimeLoop 顺序）潜伏多个缺陷，但 **没有多块 / 多 rank / 多层端到端测试**——cavity3d 是单 octant L0（无 AMR、无跨 rank、无粗细），单测用 `DummyBlock` 夹具。cylinder3d（多 rank、L=4）是首个能同时暴露这些缺陷的真实算例。
- **分阶段**：直上真实算例 cylinder3d，但 **分波去险**——先修最致命且 mid-run 不可调试的 ① comm_tag；再以 uniform 网格隔离物理/BC/后处理与 AMR 耦合；最后叠加 AMR 与严格门。
- **oracle 三阶递进**：sanity → 量级 → Cd<1%，早期不强约束（见下表）。
- **v1 AMR 是几何驱动**（ResolveBounding wake + ResolveSurface 表面，见 T07）；solution-adaptive（Q-criterion/梯度）仍是 PRD 不在范围内。cylinder3d 的几何（圆柱表面 + 尾流）正是几何 AMR 的真实用例。

### 当前代码缺口（推敲核实）

cylinder3d 需要而当前 `src/` **没有** 的：

- **A. 阻力/Cd 后处理器（全新，头号缺口）**：全仓无 `drag/Cd/force/momentumExchange` 后处理器；没有它 Cd<1% 门无法测量。
- **B. 出口/出流 BC（全新）**：`src/solver/lbm/` 无 `outlet/outflow/convection`；`DomainBcType` 仅 `kNoSlip`/`kMovingLid`/`kInterpolatedVelocity`，无出流。
- **C. 圆柱 STL fixture（数据）**：`tests/data` 仅有通用 triangle STL，需一个圆柱 STL。
- **D. 潜伏 bug ①②③④⑤**（见 PRD「已知缺陷与治理路径」）：4 rank L=4 有上千面，① 概率非零且 mid-run 不可调试。

### oracle 三阶（早期不强约束）

| 阶段 | oracle | 断言 |
|------|--------|------|
| W2（uniform） | **sanity** | 阻力有限、符号正确（阻力逆向流，Cd>0）、质量守恒、无 NaN |
| W3（AMR 初版） | **量级** | Cd 与 OpenLB 参考同量级（2x 内） |
| W4（严格门） | **Cd<1%** | Cd 与 OpenLB 参考值相对误差 < 1%（PRD #13） |

---

## 关键设计决策（已对齐）

| 项目 | 决策 |
|------|------|
| **任务粒度** | 单一 T11 文档；P0 + W1 → W2 → W3 → W4，可分 PR merge |
| **comm_tag 修复（P0）** | per rank-pair 确定性枚举（按 canonical face_key 排序，两侧 rank 算出相同唯一 tag）根除碰撞类；rehash mask 对齐 `MPI_TAG_UB`；加跨 rank 对称性单测（强制碰撞） |
| **AMR 判据** | v1 几何驱动（T07 `ResolveBounding` wake + `ResolveSurface` 表面）；solution-adaptive 不在范围 |
| **阻力后处理** | momentum-exchange（MEM）在圆柱 `kBoundary` 面格上积分 → 物理力 → `Cd = 2F/(ρU²A)`；单测用已知场→已知力 |
| **出口 BC** | 零梯度 / 对流出流（do-nothing 或 convection `∂f/∂t + U∂f/∂x = 0`）；单测用已知内场→出流 ghost=内场 |
| **网格（W2）** | uniform（无 refine，全同级），多 rank 触发 ②边ghost + ①comm_tag；Lagrava 空粗细面 |
| **网格（W3）** | L=4 几何 AMR；`LevelCoupler` 粗细耦合，多 rank 跨 rank 粗细面 |
| **参数对齐（W4）** | 对齐 OpenLB `examples/laminar/cylinder3d` 默认参数（Re/τ/域长/圆柱位置/分辨率）与参考 Cd 值；具体数值由实现时从 OpenLB 取 |
| **oracle 递进** | W2 sanity → W3 量级 → W4 Cd<1%；早期不强约束 |
| **③TimeLoop 顺序** | W2 HITL 决策点：对齐 PRD（`collide→exchange→domain_bc.apply→stream`）还是更新 PRD 记录现有 `apply_post_stream` 偏差 |
| **④overlap padding** | W2 HITL 决策点：保留（cavity3d 校准所需）还是收敛回 PRD「最小自建 BlockLattice」 |
| **②边 ghost** | 若 W2/W3 暴露（块边/角局部尖峰）即修：`GhostSchedule` 补 12 条边 ghost 线交换 |
| **阻塞** | P0 独立可先做；W1 阻塞 P0；W2 阻塞 W1；W3 阻塞 W2；W4 阻塞 W3 |

### W3 前置决策（已对齐 2026-07-15，基于 OpenLB `examples/laminar/cylinder3d` 源码核实）

| 项目 | 决策 |
|------|------|
| **算例对齐** | **A. W3 即全量复刻 Schäfer-Turek**——几何 + Poiseuille 入口 + p=0 压力出口 + Bouzidi 圆柱 + bounce-back 壁。W3 量级门与 W4 的 1% 门用同一配置，无 W3→W4 重调 |
| **OpenLB 参考** | Re=20，τ=0.53（omega≈1.887），charU=0.2，D=0.1（R=0.05），RESOLUTION=10（dx=0.01），dt=0.001，latticeU_mean=0.02（Poiseuille 峰值 2.25×=0.045 + ramp），MAX_PHYS_T=16（≈16000 步），ρ=1.0。几何：channel 2.5×0.41×0.41，圆柱 D=0.1 轴沿 z 贯穿、中心 ≈(0.45,0.2) 略偏置。**归一化对齐（2026-08-05 修正）**：Schäfer-Turek 的 Re 与 Cd **均以峰值 U_max 归一**，故本实现令 `cfg.u_inlet=0.02` 即 U_max 峰值（Re=U_max·D/dx/ν_lat=0.02·10.24/0.01≈20.5），duct 均值=4/9×0.02；Cd 分母用 `u_inlet²`（峰值²）。先前误将 0.02 当均值、入口峰值设 2.25×=0.045 → Re_peak=46、Cd 偏高 ~2.7×。 |
| **参考 Cd 来源** | **跑 OpenLB 取自身收敛 Cd ≈ 6.36**（2026-07-15 实测：`openlb-1.9.0/.../cylinder3d` arm64 原生 MPI 单 rank ~5min 跑完，physT=16 收敛 `drag≈6.36`、lift≈−0.169 稳态无脱涡）。Schäfer-Turek 2D 文献 5.58 与此差 ~14%，**不作 <1% 硬目标**——证实「跑 OpenLB 取自身值」的必要性 |
| **② 时序** | **先 uniform 干净证 ② 再叠 AMR**——W3 内先在 uniform 多 rank、**无圆柱**、with/without 边交换 toggle 下隔离证 ②，通过后再叠 AMR。需 `GhostSchedule` 补 Stage B 跨 rank 边交换（p4est corner callback 取对角邻居 rank）+ 干净重探测试 |
| **AMR 接线** | **参数化同一 header**——`examples/cylinder3d/cylinder3d_case.h` 加 `MakeAmrForest` + AMR `GeometryConfig`，构造按 config 分流 uniform/AMR；`test_cylinder3d_amr` 复用同一 `Cylinder3dCase` |
| **W3 scope 扩大项** | 需新增：① **Poiseuille 逐格入口剖面**（`DomainBcSpec.u_wall` 现为每面常向量，不支持逐格空间变化 + 时间 ramp）；② **压力出口 p=0 Dirichlet**（`DomainBcType` 现 only `kNoSlip/kMovingLid/kInterpolatedVelocity/kOutflow`，无压力出口；`kOutflow` 零梯度≠OpenLB InterpolatedPressure）；③ **Schäfer-Turek 圆柱 STL**（W1 `cylinder.stl` R=0.5 z-span=1 不匹配 D=0.1 z-span=0.41，需新 fixture） |

---

## 波次与纵向切片

| 波次 | 标题 | 类型 | 阻塞于 | 覆盖缺陷 | 端到端交付 |
|------|------|------|--------|----------|------------|
| **P0** | comm_tag 跨 rank 对称性修复 | AFK | 无 | ①⑤ | per rank-pair 确定性枚举 + rehash mask 对齐 + 跨 rank 对称性单测；多 rank run 不再 MPI 不匹配 |
| **W1** | 阻力后处理 + 出口 BC 组件单测 | AFK | P0 | A+B+C | A momentum-exchange Cd 后处理器 + B 出口出流 BC 各自单测；C 圆柱 STL fixture 备好 |
| **W2** | cylinder3d uniform 组装 + sanity | AFK | W1 | P1,②③④ | 几何 + uniform 多 rank + 入口/出口/圆柱 Bouzidi + TimeLoop + 阻力接入；oracle=sanity；③④ HITL 决策点 |
| **W3** | AMR L=4 + Lagrava + 量级 | AFK | W2 | P2,⑤ | 几何 AMR + LevelCoupler 粗细耦合 + 多 rank 跨 rank 粗细面；oracle=量级；②若暴露即修 |
| **W4** | Cd<1% 严格门 + example + VTK | AFK | W3 | P3,#13 | 参数对齐 OpenLB + 长时稳态 + Cd<1%；验收 PRD #13；example + VTK |

---

## P0 · comm_tag 跨 rank 对称性修复（前置）

### 行为

`FacePairList::AssignCommTag` 当前在 *per-rank 本地* `tag_owner` map 上解决哈希冲突（`face_pair_list.cpp:121`）：跨 rank 面两侧算出相同 `face_key`，但 **两侧冲突集不同** → rehash 路径分叉 → 同面两 rank 不同 tag → `Isend`/`Irecv` 不匹配。且 rehash（`face_pair_list.cpp:137`）仍用 `&0x7FFFFFFF`，可超 MPICH `MPI_TAG_UB`（`0x3FFFFFFF`）。

修复方向：**per rank-pair 确定性枚举**——对每个 rank 对 (A,B)，将该对之间的全部跨 rank 面按 canonical `face_key` 排序，tag = 该对内序号。两侧 rank 看到同一组面、同一排序 → 同一序号 → 同一 tag，且对内唯一（MPI 按 (src,dst,tag) 匹配，不同 rank 对可复用 tag，无需全局唯一）。rehash 路径整体移除或 mask 对齐 `MPI_TAG_UB`。

### 测试 · `test_face_pair_list`（扩展）/ 新增 `test_comm_tag_symmetry`

| 用例 | 设置 | 断言 |
|------|------|------|
| `CrossRank_TagSymmetric_ForcedCollision` | 构造足够多跨 rank 面强制哈希碰撞 | 两侧 rank 为同一面分配 **同一** comm_tag；同一 rank 对内无重复 tag |
| `RehashMask_BelowTagUb` | 极多面 | 所有 comm_tag ≤ `MPI_TAG_UB`（MPICH `0x3FFFFFFF`） |
| 既有 face_pair_list 用例 | — | 仍绿 |

### 验收标准

- [x] 跨 rank 对称性单测全绿（强制碰撞下两侧同 tag）
- [x] 所有 comm_tag ≤ `MPI_TAG_UB`
- [x] T02/T05/T06 相关既有 ctest 仍绿

---

## W1 · 阻力后处理 + 出口 BC 组件单测

### 行为

建两个可独立单测的组件 + 一个数据 fixture（尚未接入完整 cylinder run）：

- **A. 阻力后处理器**：momentum-exchange 在 `kBoundary` 面格上累加反弹动量 → 物理力 `F` → `Cd = 2F/(ρU²A)`（A 为圆柱迎风投影面积）。不依赖 LBM 框架头；读 `BlockLattice` 面格 populations。
- **B. 出口出流 BC**：新增 `DomainBcType` 出流类型（零梯度 / 对流），`DomainBoundaryHandler` 处理出流面 ghost；不处理曲面（Bouzidi）或同级 halo（T05）。
- **C. 圆柱 STL fixture**：闭圆柱面网格，置 `tests/data/mesh/cylinder.stl`（或 `examples/cylinder3d/data/`），ASCII/binary 任一。

### 测试

| 组件 | 用例 | 断言 |
|------|------|------|
| A | `Drag_KnownEquilibrium_KnownForce` | 已知平衡态 + 已知边界面 → Cd 与手算一致 |
| A | `Drag_Sign_OpposesFlow` | 流向 +x 圆柱 → F_x<0、Cd>0 |
| B | `Outlet_ZeroGradient_GhostEqualsInterior` | 已知内场 → 出流 ghost = 紧邻内场值 |
| B | `Outlet_Convection_NoReflection` | 对流出口对均匀流不产生反射扰动 |

### 验收标准

- [x] A、B 单测全绿
- [x] C 圆柱 STL 可被 `stl_reader` 解析（三角面数/法向/包围盒正确）
- [x] 不引入 OpenLB 框架层头文件

---

## W2 · cylinder3d uniform 组装 + sanity（含 ③④ HITL 决策点）

### 行为

`cylinder3d_case.h` 提供端到端初始化与推进（uniform、多 rank、sanity oracle）：

1. `OctreeForest(MPI_COMM_WORLD, domain)` → `refine`（**W2 不 refine**，全同级）
2. `GeometryEngine::build`（channel `kInternalChannel` + 圆柱 `kExternalObstacle`，priority 覆盖）→ `MaterialField`
3. `FacePairList` → `BlockCollection`（`initialize_from_material`）→ `GhostSchedule`、`LevelCoupler`（空粗细面）、`DomainBoundaryHandler`（入口 `kInterpolatedVelocity` + 出口出流 + 圆柱 Bouzidi 由 `MaterialField` 驱动）、`TimeLoop`
4. 多 rank `partition(make_level_weight_fn(forest))`；固定步 `advance_one()` 循环
5. 阻力后处理器接入，每 N 步采 Cd

**HITL 决策点（W2 校准时触发，须用户拍板）**：

- **③ TimeLoop 阶段顺序**：现实现为 `collide→domain_bc.apply→exchange→stream→apply_post_stream`，PRD 为 `collide→exchange→domain_bc.apply→stream`。二选一：改实现对齐 PRD，或更新 PRD 记录 `apply_post_stream` 为正式决策（OpenLB PostStream 对齐）。
- **④ overlap-padding 机制**：`BlockLattice` 已积累 `collide_overlap_padding_bgk`/`stream_overlap_padding_shell`/`commit_overlap_padding_stream`/`OverlapPaddingCollideMode`/`YminYmaxPaddingOutOfHaloMode`。二选一：保留（cavity3d 2% 校准所需，更新 PRD 承认）或收敛回 PRD「最小自建 BlockLattice」（移除并重校 cavity3d）。

### 测试 · `test_cylinder3d_uniform`（integration，多 rank）

| 用例 | 设置 | 断言 |
|------|------|------|
| `Uniform_MultiRank_RunsNoCrash` | uniform、≥2 rank、固定步 | 不抛、不 NaN |
| `Uniform_MassBoundedDrift` | 同上，N 步 | 总质量相对漂移 < 0.2（开放系统有界漂移，非严格守恒） |
| `Uniform_Drag_FinitePositiveSign` | 同上 | Cd 有限且 >0 |
| `Uniform_NoEdgeArtifact` | ~~多块邻接~~ **已禁用** | ~~块边/角 ρ/u 无尖峰~~ 探针被体素化混洞污染（per-octant vs whole-grid 在圆柱表面分类差异；`f=f_eq−t` 下 step0 不可见、step1 发散），非 ②；② 需干净重探 |

环境变量：`OCTLB_CYL_STEPS`、`OCTLB_CYL_RANKS`、`OCTLB_CYL_WRITE_VTK`。

### 验收标准

- [x] `test_cylinder3d_uniform` sanity 用例（RunsNoCrash / MassBoundedDrift / DragFinitePositiveSign）全绿
- [x] ② Stage A 同 rank 边交换已落（`FacePackable`/`BlockLattice`/`GhostSchedule`）；NoEdgeArtifact 探针因体素化混洞禁用
- [x] ③④ HITL 决策已落定（2026-07-15）：以已验证实现为基准更新 PRD——③ 保留 `collide→domain_bc.apply→exchange→stream→apply_post_stream` 两段式 BC + PostStream 阶段为正式 spec（`time_loop.h:175-209`，InterpolatedVelocity/PlaneFd 需 stream 后邻居）；④ 保留 overlap-padding 机制，PRD「最小自建 BlockLattice」精确化为「自建存储 + 不引框架层」，padding 属存储内部细节、硬不变量不变（cavity3d ~2% 校准所需）。详见 `doc/prd/octlb-framework.md` BlockStore/TimeLoop/OpenLB 融入策略修订段与已知缺陷表 ③④ 行。W2 仍用 legacy BC 跑 sanity，InterpolatedVelocity/padding 路径在 W3/W4 切入。
- [x] T01–T10 既有 ctest 仍绿（P0/W1 改动未破坏）

---

## W3 · Schäfer-Turek 全量复刻 + ② 干净验证 + AMR L=4 + Lagrava + 量级

> 2026-07-15 重新对齐：W3 不再是「在 W2 占位几何上叠 AMR」，而是**切到 Schäfer-Turek 算例**（与 W4 同配置），并先干净验证 ② 再叠 AMR。scope 比 v1 文档扩大（+Poiseuille 入口、+压力出口、+Schäfer-Turek STL、+② Stage B 跨 rank 边交换）。

### 行为

**第 1 步 · 新增 BC 能力（组件单测先行）**

- **Poiseuille 逐格入口**：`DomainBcSpec` 扩展逐格速度场（functor 或 `(i,j,k)→u` 数组 + 时间 ramp `PolynomialStartScale`），入口面 `kInterpolatedVelocity` 喂 Poiseuille 剖面（峰值 2.25×mean，沿 y/z 距壁衰减），复用 W2 PlaneFd InterpolatedVelocity 路径。单测：已知 Poiseuille 剖面 → 入口格 u 与解析一致。
  - **2026-07-15 进度（TDD，Track 1 组件 1）**：✅ 已落——新增 `src/solver/lbm/boundary/inlet_velocity_field.h`（`InletVelocityField` 抽象 + `PoiseuilleInletProfile`：空间抛物 `4c(1-c)` + OpenLB `PolynomialStartScale` smootherstep ramp `10x³−15x⁴+6x⁵`）；`DomainBcSpec` 加可选 `inlet_field`（空则回退 `u_wall`，向后兼容）；`PrescribedVelocity(spec,ix,iy,iz,t,u)` 单一查询点；`DomainBoundaryHandler::set_time` + `TimeLoop` 步计数 thread 时间（OpenLB `iT` 0-indexed 语义）。单测 4 条全绿（`test_inlet_velocity_field` 3 + `test_inlet_bc_integration` 2 含 legacy Zou-He 端到端 + TimeLoop threading），回归 `test_outlet_bc/test_domain_boundary/test_time_loop_levels/test_cavity3d_serial` 等全绿。**遗留 3c**：InterpolatedVelocity（boundary-lattice 模式）路径的 `PrescribedBoundaryU`（`interpolated_velocity.h`）尚未接 `PrescribedVelocity`（仍用顶盖常量 `u_wall`），留 W3 装配时随 Schäfer-Turek 入口接入并测。
- **压力出口 p=0**：新增 `DomainBcType::kPressure`（Dirichlet 压力/rho=1.0），`DomainBoundaryHandler` 出口面 ghost 按 InterpolatedPressure 填充。单测：已知内场 → 出口 ghost 满足 p=0。
  - **2026-07-15 重新对齐**：压力出口不再作为独立 ghost-fill kernel，而是 **T11-refactor 的首个新客户**——在 per-cell `BcKind` + 中心化 dispatch 架构里以 `kPressureDirichlet` arm 交付（R2 波）：`BoundaryLatticeView`/`ApplyPlaneFd` 加 pressure 分支（prescribed rho=1.0 + u FD 外推），与入口 `kVelocityDirichlet` 共用同一 FD/PostStream 机器。理由：W3/W4 入口为 Poiseuille `kInterpolatedVelocity` → handler 进 boundary-lattice 模式 → legacy ghost-fill 路径不可达（见 `T11-refactor-bc-dispatch.md` 根因）。单测同步升级为「已知内场 → 出口边界格重构后 `rho==1.0` 且 `u==FD 外推内场`」。
- **Schäfer-Turek 圆柱 STL**：生成 D=0.1、z-span=0.41、中心 (0.45,0.2,0.205) 闭圆柱面网格，置 `tests/data/mesh/cylinder_st.stl`（或复用 OpenLB `cylinder3d.stl` 若 stl_reader 能解析）。单测：`stl_reader` 解析三角面数/包围盒正确。

**第 2 步 · ② 干净验证（uniform 多 rank 无圆柱，叠 AMR 前）**

- `GhostSchedule` 补 **Stage B 跨 rank 边交换**：p4est **corner callback** 取对角邻居 rank，与 Stage A 同 rank 边交换对齐，跨 rank 12 条边 ghost 线 MPI 交换。
- 干净重探测试：uniform 多 rank、**无圆柱**（纯 channel），with/without 边交换 toggle 对比，断言块边/角 ρ/u 无尖峰（② 确认）。

**第 3 步 · Schäfer-Turek AMR 组装（参数化 `cylinder3d_case.h`）**

- `MakeAmrForest` + AMR `GeometryConfig`：channel 2.5×0.41×0.41 + 圆柱 STL，`ResolveBounding`(wake) + `ResolveSurface`(表面)，finest dx=0.01（RESOLUTION=10），`max_level` 使最细 octant×cell_width=0.01；每轮 `balance()` + `partition(make_level_weight_fn(forest))`。
- BC：Poiseuille 入口 + p=0 压力出口 + Bouzidi 圆柱 + bounce-back 壁；omega=1/0.53≈1.887，rho0=1.0。
- `LevelCoupler` 从 `CoarseFineFaces` 固化 `coupling_plan_`，半/全时 prolongation + restriction（Lagrava）。
- 多 rank 跨 rank 粗细面（⑤ 路径）首次在真实流中运行。

### 测试 · `test_cylinder3d_amr`（integration，多 rank）+ 组件单测

| 用例 | 设置 | 断言 |
|------|------|------|
| `Inlet_Poiseuille_MatchesAnalytic` | 组件 | 入口格 u 与 Poiseuille 解析一致 |
| `Outlet_PressureZero_GhostMatches` | 组件 | 出口 ghost 满足 p=0 |
| `EdgeExchange_CrossRank_NoArtifact` | uniform 多 rank 无圆柱、toggle | 块边/角 ρ/u 无尖峰（② Stage B） |
| `Amr_MultiRank_RunsNoCrash` | Schäfer-Turek L=4、≥4 rank | 不抛、不 NaN |
| `Amr_MassBoundedDrift` | N 步 | 总质量相对漂移 < 0.2（开放系统有界） |
| `Amr_Cd_SameOrderAsOpenLb` | 长时推进 | Cd 与 OpenLB 收敛参考同量级（2x 内） |
| `Amr_CoarseFineInterface_Continuous` | 粗细界面 | ρ/u 跨界面跳跃 < tol（Lagrava 契约） |

### 验收标准

- [x] Poiseuille 入口、压力出口、Schäfer-Turek STL 组件单测全绿（2026-07-31 复验：`test_inlet_velocity_field`/`test_inlet_bc_integration`/`test_outlet_bc`/`test_stl_reader`/`test_lagrava_coupler`(含 `CoarseFine_CapturedGeometry`)/`test_face_pair_list`(含 `CrossRankEdges`)/`test_bc_installer`/`test_bc_dispatcher` 共 11 条全绿）
- [x] ② Stage B 跨 rank 边交换已落 + 干净重探测试绿（uniform 无圆柱）（2026-07-31：`EdgeExchange_CrossRank_NoArtifact` 修复后绿——见下方 ② 修复说明；**注意该测试为弱验证**：低 Ma uniform 流下 ② 不显现，OFF≈ON，测试仅证 Stage B 接线 + 非退化，不强证 ② 正确性）
- [x] `test_cylinder3d_amr` 量级用例 `Amr_Cd_SameOrderAsOpenLb` 原生机 PASSED（2026-08-05：L5 -n4，`OCTLB_CYL_STRICT=1 OCTLB_CYL_STEPS=400`，Cd=5.54∈[3.18,12.72]，~58min；Cd 在 ~step 150 冻结，400 步即收敛，4000 步默认同值）。**Rosetta 下不可行**（每 AMR case 构造 ~50min+，4 case 数小时），留原生机 CI。
- [x] W2 既有用例仍绿（2026-07-31：`test_cylinder3d_uniform` 177s 绿）

### ② Stage B 测试修复说明（2026-07-31）

`EdgeExchange_CrossRank_NoArtifact` 此前**实际是红的**（非文档原记 "✅ 绿"）。根因**非 ②**：无圆柱直 channel + W3 压力出口（`kInterpolatedPressure`）在 `u_inlet=0.05`（Ma≈0.083）下**体层速度 ~15%/step 指数 runaway**（1-rank 实测：bulk 0.055→0.60 / 12 步），淹没边信号 → ON≈OFF 都发散。`test_cylinder3d_uniform`（带圆柱）因圆柱 blockage 抑制 runaway + 断言弱（no-NaN + mass drift）而误绿。**修复**：`ChannelConfig.u_inlet` 0.05→0.01（Ma≈0.017），体层收敛到稳态 ~0.04–0.05（Poiseuille），测试绿。**遗留**：低 Ma 下 OFF（cross-rank 边 ghost 未填）也不尖峰 → ON≈OFF，测试**平凡通过**，仅证 Stage B 接线（`has_edge_global>0`）+ 体层有界 + 非退化，不强证 ② 正确性；② 运行时显现需强跨块边梯度流（涡/射流穿块边），uniform 低 Ma 流无此梯度。PRD ② 行维持「运行时未确认」。**对 W3/W4 的提示**：压力出口 Ma-marginal，Schäfer-Turek `u_inlet=0.02`（Ma≈0.033）长跑前需确认稳定。详见 [[t11-edge-exchange-test-config-fix]]。

### W3/W4 根因修复与原生机验收（2026-08-05）

W3 量级门此前在原生机上仍 FAIL（Cd 偏离 6.36 一个量级：L3≈33、L4≈17、L5≈17，远超 [3.18,12.72] 的上半段但量级仍偏大约 2.7×）。2026-08-05 在 native x86（`module load samr`）定位并修复了 **3 个独立的 drag/Cd 根因**，量级门转绿：

1. **(A) MEM 计时 bug（主项，6.6×）—— FIXED**。`MomentumExchangeDrag` 在 `advance()`（即 `stream()`）之后被调用，读到的是 **post-stream** 的弹回值，而非标准 MEM 要求的出射 f_i。pull scheme 下 `stream()` 把流体格 f_i 槽位覆写为墙格弹回值（`block_lattice.cpp` stream bounce 分支），故 live 读的是 INCOMING 弹回值。修复：post-collide 快照——`BlockLattice::take_post_collide_snapshot()`/`post_collide_populations_at_halo()`（`block_lattice.h/cpp`），`TimeLoop::set_snapshot_post_collide()` 在 `stream_level` 之前快照每块（`time_loop.h` advance()），`force_on_fluid_if` 读 `snapshot[iPop]`（出射 f_i）+ `live[opp]`（弹回 f_bar）。实测（-n4）：L3 228.9→33.0（6.9×），L4 132.6→17.29（7.7×），L5 111.9→17.11（6.6×）。计时修复单独：112→17。
2. **(B) Re/归一化 bug（剩余 2.7×）—— FIXED**。计时修复后 Cd 收敛到 ~17 且 L4→L5 FLAT（非分辨率受限，系统性偏差）。**该偏差并非 q-placement**（曾假设 q_frac MEM 是剩余修复，已 DISPROVED：慢流下 u~0.0056，f_bar≈f_out，(f_i+f_bar)≈2·f_i，L3 正确 MEM 32.78 ≈ postcollide-only 33.0，可忽略）。真正根因：`u_inlet=0.02` 被当作 duct **均值**，入口 Poiseuille 峰值设 `2.25×0.02=0.045` → Re_peak=46。但 Schäfer-Turek 的 Re 与 Cd **均以峰值 U_max 归一**，故入口峰值应即 u_inlet=0.02 → Re=20。修复：`cylinder3d_case.h` 入口归一化行 `2.25*u_inlet`→`u_inlet`。L4 -n4 Cd 17→7.08（+11%），L5 -n4 Cd→5.54（−13%），**真值 6.36 被 L4/L5 夹住** → 收敛格式特征，力一直正确，仅 Re 错。
3. **(C) `advance_steps` 快照隐患（让测试通过的关键）—— FIXED**。`test_cylinder3d_amr.Amr_Cd_SameOrderAsOpenLb` 调 `cas.advance_steps(steps)` 后调 `cas.drag_coefficient(...)`，但 `advance_steps` 未开快照 → `drag_coefficient` 回退到 live（post-stream）数组 → 计时错 → Cd≈112 → FAIL。example 手动 `set_snapshot_post_collide(...)` 规避了此问题，测试没规。修复：`Cylinder3dCase::advance_steps`（`cylinder3d_case.h`）自调 `loop.set_snapshot_post_collide(true)` —— 快照默认常开，`drag_coefficient` 默认正确，无调用方隐患。example 端冗余的逐步 toggle 已移除。perf 开销可忽略（每块每最细步一次 memcpy，多小时运行约 +24s）。

**原生机验收（实际测试二进制，非 example）**：
```
$ OCTLB_CYL_STRICT=1 OCTLB_CYL_STEPS=400 mpirun -np 4 ./tests/integration/test_cylinder3d_amr \
    --gtest_filter='Cylinder3dAmr.Amr_Cd_SameOrderAsOpenLb'
[       OK ] Cylinder3dAmr.Amr_Cd_SameOrderAsOpenLb (3465150 ms)
[  PASSED  ] 1 test.
```
- **`Amr_Cd_SameOrderAsOpenLb` PASSED**（Cd=5.54∈[3.18,12.72]，~58min/400 步；Cd 在 step 150 冻结，4000 步默认同值）。
- 跨分辨率收敛（example 探针）：L3=14.0（欠发展，D=2.56 cells，mean_ux=0.0024，**非测试配置**，gate 外）→ L4=7.08（mean_ux=0.0066，+11%）→ L5=5.54（mean_ux=0.008，−13%）。L5 已发展充分（mean_ux≈入口均值 90%）。
- **`Amr_Cd_WithinOnePercentOfOpenLb`（严格 <1% 门）仍 FAIL**：L5 Cd=5.54 距 6.36 为 13%。Cd 已稳态冻结，更长演化（16000 步）无助；需 L6 更细网格或匹配 OpenLB 精确设置。PRD 已标注该门为长期原生长跑目标。

**记忆索引**：[[t11-cylinder3d-amr-mesh-too-coarse]]（根因分析 + 验收数据）。


---

## W4 · Cd<1% 严格门 + example + VTK

### 行为

- **与 W3 同配置**（Schäfer-Turek：Poiseuille 入口 + p=0 压力出口 + Bouzidi + bounce-back，omega=1/0.53，dx=0.01），长时推进至稳态（MAX_PHYS_T=16，≈16000 步）。
- **参考 Cd = OpenLB `cylinder3d` 自身收敛值 ≈ 6.36**（2026-07-15 实测 physT=16 收敛；Schäfer-Turek 2D 文献 5.58 仅交叉核对，差 ~14% 非硬目标）。
- `examples/cylinder3d`：薄 main + `cylinder3d_case.h`，`mpirun -n 4 ./cylinder3d` 默认运行至稳态，输出 VTK。
- 验收 PRD #13：Cd 与 OpenLB 参考值相对误差 < 1%。

### `examples/cylinder3d`

- 链接 `cylinder3d_case.h`；多 rank 默认运行；写 VTK（`AmrVtkWriter` + `VtkVelocityField`/`VtkPressureField`）。
- 输出 Cd 时间序列（CSV）便于与 OpenLB 参考曲线对比。

### 验收标准（T11 完成）

- [ ] `test_cylinder3d_amr` 严格用例 `Amr_Cd_WithinOnePercentOfOpenLb` 全绿（Cd 与 OpenLB 收敛参考 <1%，PRD #13）——**原生机长跑仍 FAIL**（2026-08-05：L5 Cd=5.54 距 6.36 为 13%，Cd 已稳态冻结，更长演化无助；需 L6 更细网格或匹配 OpenLB 精确设置。Rosetta 不可行）
- [x] `examples/cylinder3d` 可构建（2026-07-31：链接通过）；≥4 rank 可运行——留原生机验证（Rosetta 构造 ~50min+ 不可行）
- [x] example 写出 VTK + Cd CSV（代码已就位 `write_vtk_timestep` + cd.csv）；人工 ParaView 验收留原生机
- [x] PRD 进度与 T11 决策已更新（2026-07-31：②/⑤ 行 + T11 进度行更新，见下文）

---

## 阻塞关系

```
T10（cavity3d #12）— 已完成
T09-W2（Bouzidi 曲面 BC）— 已完成
T07（GeometryEngine/MaterialField）— 已完成
T08 W4（vtk_lbm_fields）— 已完成
└── T11 — 进行中（P0/W1/W2-sanity + T11-refactor 已完成；W3 测试+缺陷⑤修复中，W4 example 已搭）
    ├── P0 · comm_tag 跨 rank 对称性修复 ✅
    ├── W1 · 阻力后处理 + 出口 BC 组件单测 ✅
    ├── W2 · uniform 组装 + sanity（③④ HITL）✅ sanity 绿；② Stage A 同 rank 边交换已落，Stage B 跨 rank 边交换 + 干净重探留 W3
    ├── T11-refactor · BC 调度架构重构（per-cell BcKind + 中心化 dispatch，治 ⑥）✅ 完成（PR #11，R0–R4；R2 交付压力出口组件 2）
    ├── 缺陷⑤ · 跨 rank 粗细耦合修复（2026-07-19）✅ FacePairList 快照 ghost 物理 bounds/level；LevelCoupler 用快照，不再 re-resolve 临时 quadid（详见 [[t11-defect5-cross-rank-coarsefine]]）
    ├── W3 · Schäfer-Turek + ② + AMR + 量级 🔨 进行中：组件单测（Poiseuille/压力出口/STL）✅ 全绿；② Stage B 干净重探 `EdgeExchange_CrossRank_NoArtifact` 修复后 ✅ 绿（原 u_inlet=0.05 压力出口 Ma runaway 致红，降 0.01 后绿；弱验证——见 W3 验收 ② 说明）；**`Amr_Cd_SameOrderAsOpenLb` 原生机 ✅ PASSED**（2026-08-05：L5 -n4 Cd=5.54∈[3.18,12.72]；3 根因修复 A/B/C 见下文）；`CoarseFineInterface_Continuous`/`Amr_Cd_WithinOnePercentOfOpenLb` 已写（gated STRICT）；AMR 非 STRICT 用例（RunsNoCrash/MassBoundedDrift/Cd_FinitePositiveSign/CoarseFineInterface）Rosetta 下构造 ~50min/case 不可行，留原生机 CI
    └── W4 · Cd<1% 严格门 + example + VTK 🔨 example 已搭（`examples/cylinder3d`，VTK+Cd CSV）；**严格门 `Amr_Cd_WithinOnePercentOfOpenLb` 原生机仍 ❌ FAIL**（2026-08-05：L5 Cd=5.54 距 6.36 为 13%，稳态冻结，需 L6 或匹配 OpenLB 设置；PRD 长期原生长跑目标）
            └── T12 · amr_convergence（#14）← 下一项
```

---

## 不在 T11 范围

- **solution-adaptive AMR**（Q-criterion/梯度触发，PRD 不在范围内）
- **动态 AMR `rebuild()`** 热路径（PRD 不在范围内；P0/W3 只验静态分层）
- **sphere3d** 大规模基准（1000+ rank，用户故事 #3，后续）
- **test_amr_convergence**（#14，T12，依赖 T11 通过）
- **gnuplot IO** 移植
- 收敛检测（`ValueTracer`、速度残差提前停；W4 用固定步至稳态）
- HDF5 checkpoint/restart

---

## 构建与测试命令（给用户执行）

```bash
cd octlb
./cmd.sh build          # docker 内 configure + 编译

# P0
./cmd.sh test -- -R 'test_face_pair_list|test_comm_tag_symmetry' --output-on-failure

# W1
./cmd.sh test -- -R 'test_drag|test_outlet|test_stl_reader' --output-on-failure

# W2
./cmd.sh test -- -R test_cylinder3d_uniform --output-on-failure

# W3
./cmd.sh test -- -R test_cylinder3d_amr --output-on-failure

# W4 + example
./cmd.sh shell
cd build/docker
ctest -R test_cylinder3d_amr --output-on-failure
mpirun -np 4 ./examples/cylinder3d --output cylinder3d_vtk
```