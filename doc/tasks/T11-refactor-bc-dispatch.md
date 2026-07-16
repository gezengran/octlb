# T11-refactor · BC 调度架构重构（per-cell BcKind + 中心化 dispatch）

> 类型：AFK 架构重构（跨模块；TDD 纵向切片，每波保持既有测试绿）
> 阻塞于：T09-W1/W2（现有 BC 模型：`CellKind` + `DomainBcSpec` + `DomainBoundaryHandler`）、T11-W2-sanity（legacy BC 已在 uniform sanity 验证）
> 阻塞：T11-W3 组件 2（压力出口 p=0，作为首个新客户）、T11-W3/W4、后续新算例（sphere3d 等）
> 对应 PRD：`doc/prd/octlb-framework.md`（BlockStore / DomainBoundaryHandler / BouzidiLinkData / OpenLB 融入策略节，本次同步更新）
> 状态：**已完成（2026-07-16，R0–R4 TDD 落地）**——决策已对齐（见下「关键设计决策」），R0–R4 已按 TDD 纵向切片实现并通过验收。

---

## 要做什么

把 OctLB 的域边界 BC 模型从 **「每面 `DomainBcSpec` + 整 handler `boundary_lattice_mode_` 单一布尔」** 重构为 **「per-cell `BcKind` + 中心化 dispatch」**，移植 OpenLB「per-cell Dynamics 多态 + material 号映射」的**架构思想**，但**不引入其框架层**（`DynamicsVector`/`MixinsDynamics`/`ExternalField`/虚函数分发），保持 OctLB「自建存储 + 不引框架层」立论与所有硬不变量。

### 根因（为什么要重构）

现状无法表达「一个 block 内混合多种 BC」。证据：`examples/cylinder3d/cylinder3d_case.h:24` 注释——`boundary_lattice_mode=true` 会让 `DomainBoundaryHandler::apply()` 把**所有** `kBoundary` 格（含圆柱表面）统一按树面 Dirichlet 速度碰撞。Schäfer-Turek 需要「入口=速度 + 出口=压力 + 圆柱=Bouzidi + 壁=bounce-back」共存于同一 block，现状表达不了。每接入一个 BC 组合不同的新算例，都得在唯一一条 FD 路径里加 `spec.type` 分支——重蹈覆辙。

OpenLB 用 per-cell Dynamics 多态消解该多样性；OctLB 因「不引框架层」约束没有 per-cell `Dynamics*` 机器，需用 **`BcKind` 枚举 switch + 中心化 dispatch** 达到同样的解耦：**BC 类型成为 per-cell 属性（setup 时 stamp），核心循环按 cell 分发，无全局 mode 标志；新增算例 = 新的 setup stamping 规则，不改核心循环。**

---

## 关键设计决策（已对齐 2026-07-15）

| 项目 | 决策 |
|------|------|
| **分发机制** | **`BcKind` 枚举 switch**（OctLB-flavored：per-cell uint8、无虚函数、热路径缓存友好）。若后续发现枚举不能覆盖某 BC 或出现分发错误，再切换到 OpenLB 已验证的 per-cell `Dynamics*` 多态（虚函数）作为兜底 |
| **MaterialKind seam** | `MaterialKind` **保持几何 only** `{kFluid, kSolid, kBoundary}`（mesh 不含物理，硬不变量不破）。BC 角色由**求解器侧**在 setup 时解析：`(MaterialField 几何 + TreeBoundaryFace 面角色 + BouzidiLinkData 标记) → per-cell BcKind`。面角色来自 `TreeBoundaryFace`（mesh 拓扑，非物理），seam 不破 |
| **顺序** | **先重构，压力出口（T11 W3 组件 2）作为首个新客户**——一次连贯推进，不在旧架构上写压力出口再返工 |
| **prescribed 值存储** | `BcKind` 只管 dispatch 类别；per-cell prescribed 值（`u_wall` / 入口 `inlet_field` / `rho_target`）仍由**面 → `DomainBcSpec`** 查表得到（Bouzidi 值来自 `BouzidiLinkData`）。dispatch 按 cell 的 `BcKind`，值按 cell 所属面的 spec——二者正交，混合 BC 自由共存 |
| **CMake / .cpp 不变量** | `BcDispatcher` / `BcInstaller` **header-only inline**（`solver/lbm/` 唯一 `.cpp` 仍为 `block_lattice.cpp`，PRD 硬不变量不破） |
| **复用而非重写** | 现有 kernel（`ApplyNoSlipGhost` / `ApplyZouHeVelocityGhost` / `ApplyOutflowGhost` / `ApplyPlaneFdBoundary` / Bouzidi pull）作为 switch 各 arm 的实现，**算法不重写**，仅重组 dispatch |
| **`cell_kind` 兼容** | `CellKind`（T09-W1）在 R0 被 `BcKind` 超集取代；提供 `CellKind→BcKind` 默认映射使既有代码/测试在重构期逐步迁移、每波保持绿 |

---

## 目标架构

### `BcKind`（`solver/lbm/bc_kind.h`，物理侧）

```cpp
enum class BcKind : std::uint8_t {
  kBulk,                // 流体 BGK
  kBounceBack,          // 静止壁（full-way bounce-back）
  kMovingBounceBack,    // 动壁（bounce-back + 壁速，cavity3d 顶盖）
  kBouzidi,             // 曲面（圆柱），stream-time pull（T09-W2）
  kVelocityDirichlet,   // 入口 InterpolatedVelocity FD（prescribed u）
  kPressureDirichlet,   // 出口 InterpolatedPressure FD（prescribed rho）  ← 首个新客户
  kOutflow,             // 零梯度 do-nothing（W1，保留为一种 kind）
  kSolid,               // 惰性固体
};
```

### `BcDispatcher`（`solver/lbm/bc_dispatcher.h`，header-only）

三个入口，按 `lat.bc_kind(ix,iy,iz)` switch，各 arm 委托现有 kernel：
- `collide(lat, ix,iy,iz, omega, specs, rho_stats, ...)`：`kBulk`→BGK；`kBounceBack/kMovingBounceBack`→反弹；`kVelocityDirichlet/kPressureDirichlet`→Dirichlet collide（prescribed u / prescribed rho）；其余跳过。
- `compute_rho_u(lat, ix,iy,iz, rho, u, specs)`：velocity 面 prescribed u + rho-from-pop；pressure 面 prescribed rho + u FD 外推；Bouzidi/bounce-back 按 OpenLB `computeU` 语义。**PostStream FD 邻域靠此多态读邻居**（OpenLB 邻居 `computeRhoU` 分发到该邻居自己的 dynamics 的 OctLB 等价）。
- `post_stream(lat, ix,iy,iz, omega, specs, padding_mode)`：仅 `kVelocityDirichlet`/`kPressureDirichlet` 走 `ApplyPlaneFdBoundary`（Edge/Corner）；其余 no-op。

### `BcInstaller`（`solver/lbm/bc_installer.h`，header-only，OpenLB `setVelocityBoundary(mat)` 等价）

setup 时一次性 stamp per-cell `BcKind`：
- 流体格 → `kBulk`；solid → `kSolid`。
- Bouzidi 标记格（`BouzidiLinkData` 覆盖的 cut-link 邻接 fluid 格）→ `kBouzidi`。
- 域外 tree 面 `kBoundary` 格：按该格所处面查 `DomainBcSpec.type` → `kVelocityDirichlet`/`kPressureDirichlet`/`kOutflow`/`kBounceBack`/`kMovingBounceBack`。
- 替换 `lattice_material_init.*` 中的 `MaterialKind→CellKind` 直映射为 `BcInstaller`（保留 `initialize_from_material` 期间兼容）。

### `DomainBoundaryHandler` 重构

- **移除**整 handler `boundary_lattice_mode_` 与 `UsesInterpolatedVelocity()` 全局门控。
- `apply()`：遍历 `TreeBoundaryFaces` 的边界格 → `BcDispatcher::collide`（per-cell）。混合 BC 自然共存。
- `apply_post_stream()`：遍历需 FD 重建的格（`kVelocityDirichlet`/`kPressureDirichlet`）→ `BcDispatcher::post_stream`。
- legacy `ApplyLegacyFaceBc`（per-face ghost-fill）在 R4 随迁移完成移除；迁移期保留以保持既有测试绿。

### `interpolated_velocity.h` 泛化

`BoundaryLatticeView` / `ApplyPlaneFdBoundary`（及 Edge/Corner）按 cell 的 `BcKind` 分两支：
- velocity：prescribed u（`inlet_field`/`u_wall`）+ rho from pop（现状）。
- pressure：prescribed `rho_target`（默认 1.0，p=0）+ u 从内场 FD 外推。
FD pi 计算两支共用（邻域 u 梯度）；仅 `WriteFromPi` 的 `rho`/`u` 来源不同。**这是 OpenLB `setVelocityBoundary`/`setPressureBoundary` 工厂对称性的 OctLB 等价。**

### Bouzidi 归并

`kBouzidi` arm 收进 `BlockLattice::stream()` pull 阶段（现状已是 stream-time），由 `BcKind==kBouzidi` 触发，不再作为独立路径。

---

## 波次（TDD 纵向切片，每波保持既有 ctest 全绿）

| 波 | 标题 | 红→绿切片 | 既有测试 |
|----|------|----------|----------|
| **R0** | `BcKind` 引入 + BlockLattice 存取 + 默认 `CellKind→BcKind` 映射 | 一条：stamp `kBounceBack` 的格 collide 后等于反弹结果 | 全绿（默认映射无行为变更） |
| **R1** | `BcDispatcher::collide` 中心化；`handler.apply` 接 per-cell dispatch | 一条：混合 `kBounceBack`+`kVelocityDirichlet` 格在同 block 各按其 kind collide | 全绿（`kBoundary` 默认走 velocity-Dirichlet，cavity3d 行为不变） |
| **R2** | `BoundaryLatticeView`/`ApplyPlaneFd` 加 pressure 分支 + `kPressureDirichlet` arm | **T11 W3 组件 2 单测**：已知内场 → 出口边界格重构后 `rho==1.0`(p=0) 且 `u==FD 外推内场` | 全绿 |
| **R3** | `BcInstaller` + `cylinder3d_case` 重排（去 `boundary_lattice_mode`）+ Bouzidi 归 `kBouzidi` arm | 一条：`BcInstaller` stamping（Schäfer-Turek 面→BcKind 分布正确）；`cylinder3d_uniform` sanity 仍绿 | `test_cylinder3d_uniform` 全绿 |
| **R4** | cavity3d + 既有 BC 测试迁移到 `BcInstaller`；移除 legacy `ApplyLegacyFaceBc` + `boundary_lattice_mode_` | 回归：`test_cavity3d_serial`(#12)、`test_domain_boundary`、`test_time_loop_levels`、`test_bouzidi_link`、`test_outlet_bc`、`test_inlet_*` 全绿 | 全绿 |

> 组件 3（Schäfer-Turek 圆柱 STL）**独立于本重构**（纯数据 + `stl_reader` 解析单测），可与本重构并行以 TDD 完成，不阻塞。

### 每轮 TDD 自检（按 `/skills/tdd`）

```
[ ] 测行为不测实现（dispatch 结果可观测，不查内部 switch 结构）
[ ] 只用公共接口（BlockLattice::bc_kind / DomainBoundaryHandler::apply 等）
[ ] 既有测试不因内部 dispatch 重组而挂
[ ] 代码仅为当前波服务，不预判后续波
```

---

## 测试

| 用例 | 设置 | 断言 |
|------|------|------|
| `BcDispatcher.BounceBack_CollidesAsBounceBack` (R0) | stamp `kBounceBack` 格 + 已知入射 pop | collide 后等于反弹结果 |
| `BcDispatcher.MixedKinds_PerCellDispatch` (R1) | 同 block 内 `kBounceBack` 与 `kVelocityDirichlet` 共存 | 各格按各自 kind collide，互不串扰 |
| `Outlet_PressureZero_GhostMatches` (R2) | 已知内场 + `kPressureDirichlet` 出口格 | 重构后边界格 `rho==1.0` 且 `u==FD 外推内场`（T11 W3 组件 2） |
| `BcInstaller.SchaeferTurek_StampingCorrect` (R3) | Schäfer-Turek 几何 + 面→BcType map + Bouzidi | 入口格 `kVelocityDirichlet`、出口格 `kPressureDirichlet`、圆柱格 `kBouzidi`、壁格 `kBounceBack` |
| `BcInstaller.Cavity_LegacyEquivalent` (R4) | cavity3d 六面 `kInterpolatedVelocity` | stamping 与旧 `MarkDomainBoundaryCellKinds` 等价；`test_cavity3d_serial` #12 L2<2% 不退化 |

### 验收标准

- [x] R0–R4 每波新增单测红→绿
- [x] **每波结束 T01–T10 既有 ctest 全绿**（#12 `test_cavity3d_serial` L2<2% 不退化为硬门；R3/R4 复跑全量 ctest + L2 通过）
- [x] `boundary_lattice_mode_` 与 `UsesInterpolatedVelocity()` 全局门控移除；legacy `ApplyLegacyFaceBc` 移除
- [x] 压力出口 `kPressureDirichlet` 作为首个新客户落位（T11 W3 组件 2 交付）
- [x] 硬不变量不破：`MaterialKind` 几何 only；`field/` 无 LBM；`io/` 无 LBM；`lbm/` 唯一 `.cpp`=`block_lattice.cpp`（`BcDispatcher`/`BcInstaller` header-only）；Solver 不调 p4est
- [x] PRD BlockStore / DomainBoundaryHandler / Bouzidi / OpenLB 融入策略节已同步

---

## 阻塞关系

```
T09-W1/W2（现 BC 模型）— 已完成
T11-W2-sanity（legacy BC uniform sanity）— 已完成
└── T11-refactor · BC 调度架构重构 ✅ 完成
    ├── R0 BcKind 引入 ✅
    ├── R1 dispatch 中心化 ✅
    ├── R2 pressure 分支 ← 交付 T11 W3 组件 2（压力出口）✅
    ├── R3 BcInstaller + cylinder3d 重排 ✅
    └── R4 迁移 + 移除 legacy ✅
        └── T11-W3 继续（AMR L=4 + Lagrava + 量级；组件 3 STL 可并行）
            └── T11-W4（Cd<1%）
                └── T12 · amr_convergence（#14）
```

---

## 不在范围内

- **per-cell `Dynamics*` 虚函数多态**：仅作枚举 switch 失效时的兜底，本重构不实现
- **动态 AMR `rebuild()` 热路径**：`BcInstaller` 静态一次 stamp；动态 AMR 的 re-stamp 留动态 AMR 任务
- **新物理**：MRT/Smagorinsky、多相、热 LBM（本重构只重组 dispatch，不加算子）
- **solution-adaptive AMR**
- **`field/` / `io/` 改动**：本重构限于 `solver/lbm/` + example/case 头；`field/`、`io/` 不动