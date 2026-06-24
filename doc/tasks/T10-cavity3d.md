# T10 · unit_converter + cavity3d 集成（W1 / W2 / W3）

> 类型：AFK（三波纵向切片）  
> 阻塞于：T09-W1（`DomainBoundaryHandler`、`TimeLoop`）；T08 W4（`vtk_lbm_fields.h`，W3 可视化 P1）  
> 对应 PRD：`octlb/doc/prd/octlb-framework.md`（集成测试 #12、`examples/cavity3d/`）  
> 状态：**完成**（2026-06-24 本地 ctest 验收：`test_unit_converter`、`test_cavity3d_serial`、`test_interpolated_velocity_compute_rho` 全绿）

---

## 要做什么

在 T09-W1 单元测试全绿后，交付 **首个端到端 LBM 集成路径**：lid-driven cavity3d（**无 STL、无 AMR**），验收 PRD 集成测试 #12，并附带可执行算例与可选 VTK 流场切片输出。

**第一版约定（已对齐 PRD / 设计讨论）**：

- **网格**：`OctreeForest` 单根 octant、`[0,1]³`、**不 refine**（`max_level=0`）；多级 cavity 留 **T12** 或后续。
- **参数**：与 OpenLB `examples/laminar/cavity3d` 默认一致（见下表），保证可比性。
- **时间推进**：**固定步数** `iT_max = getLatticeTime(MAX_PHYS_T)`；**不**实现收敛检测（`ValueTracer` / 速度残差留 example 后续或 T11+）。
- **验收**：OpenLB `cavity3d` 同默认参数、**iT=5269**（ValueTracer 收敛步）；竖直中线 17 点 `u_x/u_lid` 相对 OpenLB 参考剖面 **L2 < 2%**（采样 `x=0.5`、`z=0.5`，含 ghost 边界层插值）。Ghia Re=100 表仅作诊断，**不作**硬验收（默认参数 Re=1000）。
- **可视化（P1）**：`examples/cavity3d` 默认写 VTK（过域心 x–y 平面切片，ParaView 对比）；集成测试默认**不写** VTK，仅断言 L2。

### OpenLB 对齐参数（默认值）

| 符号 | OpenLB 参数 | 默认值 | OctLB 用途 |
|------|-------------|--------|------------|
| N | `RESOLUTION` | 30 | 块内 `BlockLattice` 尺寸 |
| τ | `LATTICE_RELAXATION_TIME` | 0.509 | `omega = 1/tau` |
| L | `PHYS_CHAR_LENGTH` | 1.0 | `BoundingBox` 边长 |
| U | `PHYS_CHAR_VELOCITY` | 1.0 | 顶盖物理速度 → `u_lid_lattice` |
| ν | `PHYS_CHAR_VISCOSITY` | 0.001 | Re = U·L/ν = **1000**（Ghia Re=100 表不匹配） |
| ρ | `PHYS_CHAR_DENSITY` | 1.0 | 初始密度 |
| `iT_max = getLatticeTime(100)` | `get_lattice_time(100)` 默认 **30000** 步（非文档误写的 ~3333；见 OpenLB 公式） |

域 BC：六面 `kInterpolatedVelocity`（OpenLB `InterpolatedVelocity` / Skordos FD）；顶面 `u_wall = (u_lid_lattice, 0, 0)`，其余面 0。

---

## 交付物

```
octlb/
├── src/solver/lbm/unit_converter/
│   ├── unit_converter.h          # W1：FromResolutionAndRelaxationTime 公式
│   └── unit_converter.cpp        # 可选：print / 自洽性辅助
├── examples/cavity3d/
│   ├── cavity3d.cpp                # W3：薄 main
│   └── cavity3d_case.h           # W2/W3：森林/格子/BC/TimeLoop 组装 + 运行 + 可选 IO
├── tests/
│   ├── unit/solver/
│   │   └── test_unit_converter.cpp   # W1
│   └── integration/
│       └── test_cavity3d_serial.cpp    # W2 smoke + W3 OpenLB L2
└── doc/tasks/T10-cavity3d.md
```

共享逻辑放在 `cavity3d_case.h`（或 `examples/cavity3d/` 下头文件，由 test 与 example 共同 include），避免重复。

---

## 关键设计决策（已对齐）

| 项目 | 决策 |
|------|------|
| **任务粒度** | 单一 T10 文档；**W1 → W2 → W3** 三波，可分 PR merge |
| **unit_converter** | OpenLB `UnitConverterFromResolutionAndRelaxationTime` 公式（`octlb` 轻量 struct）；**不**复制 OpenLB 完整 `unitConverter.h`；octree-mesh `LbmUnitConverter` 仅作输出换算参考，**不**单独作为仿真设定来源 |
| **收敛** | T10 **不**实现；固定 `iT_max`；#12 以 OctLB vs OpenLB 同参数 L2 为准 |
| **网格** | 单根 octant、L=0；`TimeLoop` 在 `max_level=0` 不进入 `LevelCoupler` |
| **Ghia 比对** | 17 点表仅作诊断；**硬验收**为 OpenLB 参考剖面 `u_x/u_lid` 相对 L2 < 2% |
| **VTK** | W3 P1；example 默认开启；ctest 默认关闭（可用环境变量 opt-in） |
| **gnuplot** | 不在 T10；可选 CSV 中线剖面供手动画图 |
| **ConstRho** | `TimeLoop` 默认 `ConstRhoStatsScope::kFluidAndBoundary`（OpenLB 对齐）；**不**做 A/B 集成对比 |
| **阻塞** | 仅 T09-W1；**不**依赖 T09-W2 / T07 |

---

## 波次与纵向切片

| 波次 | 标题 | 类型 | 阻塞于 | 用户故事 | 端到端交付 |
|------|------|------|--------|----------|------------|
| **W1** | `unit_converter` + 单测 | AFK | T04（D3Q19 `invCs2`） | — | Re/τ/N/ν 自洽；`u_lid_lattice`、`iT_max` 可算 |
| **W2** | cavity 组装 + 冒烟 | AFK | W1、T09-W1、T06 | — | 单 octant 固定步推进不崩溃；顶盖速度非零 |
| **W3** | OpenLB L2 #12 + example + VTK P1 | AFK | W2、T08 W4 | —（正确性基线） | `test_cavity3d_serial` 绿；`examples/cavity3d` 可跑 |

---

## W1 · unit_converter

### 行为

实现 `octlb::UnitConverter`（命名可微调），构造输入：

`(resolution, lattice_relaxation_time, char_phys_length, char_phys_velocity, phys_viscosity, phys_density)`

核心公式（与 OpenLB 一致）：

- `phys_delta_x = L / N`
- `phys_delta_t = (τ - 0.5) / invCs2 × (L/N)² / ν`
- `conversion_velocity = phys_delta_x / phys_delta_t`
- `char_lattice_velocity = U / conversion_velocity`
- `omega = 1 / τ`；`reynolds = U × L / ν`
- `get_lattice_time(phys_t) = round(phys_t / phys_delta_t)`

**不**引入 OpenLB 框架头（`xmlReader`、`ostreamManager` 等）。

### 测试 · `test_unit_converter`

| 用例 | 断言 |
|------|------|
| `OpenLbDefaults_Re1000` | 默认参数 Re = 1000（容差 ε） |
| `TauOmegaConsistent` | `omega = 1/tau` |
| `LatticeLidVelocity` | `char_lattice_velocity` 与手算 OpenLB 默认一致 |
| `LatticeTime_MaxPhysT` | `get_lattice_time(100)` = 30000（OpenLB 公式） |

---

## W2 · cavity 组装 + 冒烟

### 行为

`cavity3d_case.h` 提供端到端初始化与推进（无 Ghia 断言）：

1. `OctreeForest(MPI_COMM_WORLD, domain)` → `partition()`（不 refine）
2. `FacePairList` → `BlockCollection`（N³，`initialize(ρ, u=0)`）
3. `GhostSchedule`、`LevelCoupler`（空粗细面）、`ConcreteDomainBoundaryHandler`、`TimeLoop`
4. 固定 `iT_smoke`（如 50 步）`advance_one()` 循环

### 测试 · `test_cavity3d_serial`（integration，单文件）

| 用例 | 设置 | 断言 |
|------|------|------|
| `Smoke_RunsWithoutCrash` | OpenLB 默认参数，50 步 | 不抛、不 NaN |
| `Smoke_MovingLid_NonZeroUx` | 同上 | 顶盖附近 `u_x > 0` |
| `OpenLbCenterline_RelativeL2_Below2Pct` | OpenLB 默认参数，`iT=5269` | 17 点 Ghia 高度 trilinear 采样 vs OpenLB 参考，相对 L2 < 0.02 |

环境变量：

- `OCTLB_CAVITY_STEPS`：覆盖默认 `iT=5269`
- `OCTLB_CAVITY_WRITE_VTK=1`：验收用例结束后写一步 VTK

（调试期 ymax/genesis 逐步审计已随 T10 收敛从 `cavity3d_case.h` 移除；OpenLB genesis dump 仍保留在 `tests/data/` 供手动画图参考。）

---

## W3 · OpenLB 验收 + example + 可视化

### `examples/cavity3d`

- 链接共享 `cavity3d_case.h`；`mpirun -n 1 ./cavity3d` 默认 **iT=5269**（OpenLB 收敛步；可用 `--steps` / `OCTLB_CAVITY_STEPS` 覆盖；全量 `iT_max=30000` 无收敛检测会发散）
- 结束后写 VTK：`AmrVtkWriter` + `VtkVelocityField`（对齐 T08 W4 集成示例）
- 可选：输出 `centerline.csv`（`y/H`, `u/u_lid`, `ghia`）便于与参考曲线对比
- 切片可视化：ParaView 对 velocity 在 **z = 0.5** 处 Slice（等价 OpenLB 过域心 x–y 平面）

---

## 验收标准

### W1

- [x] `cmake --build build -j4` 通过
- [x] `test_unit_converter` 上表 **全绿**
- [x] 不引入 OpenLB 框架层头文件

### W2

- [x] `cavity3d_case.h` 单 octant 路径可编译
- [x] `test_cavity3d_serial` 冒烟用例 **全绿**
- [x] T01–T09 相关 ctest **仍绿**

### W3（T10 完成）

- [x] `test_cavity3d_serial` 验收用例 **全绿**（#12，OctLB vs OpenLB L2 < 2%）
- [x] `examples/cavity3d` 可构建、单 rank 可运行
- [x] example 写出至少一步 VTK（P1；人工 ParaView 验收）
- [x] PRD `octlb-framework.md` 进度与 T10 决策已更新

---

## 阻塞关系

```
T09-W1（域 BC + TimeLoop）— 已完成
T08 W4（vtk_lbm_fields）— 已完成
└── T10 — 已完成（2026-06-24）
    ├── W1 · unit_converter
    ├── W2 · cavity 组装 + 冒烟
    └── W3 · OpenLB #12 + examples/cavity3d
            └── T11 · cylinder3d（#13）← 下一项
```

---

## 不在 T10 范围

- Bouzidi / `MaterialField`（T09-W2、T11）
- 多级 AMR cavity（T12）
- `gnuplot` 移植
- 收敛检测（`ValueTracer`、速度残差提前停）
- 多 rank cavity 并行（cylinder 在 T11）

---

## 构建与测试命令（给用户执行）

```bash
cd octlb
module load octlb
cmake -B build -DCMAKE_BUILD_TYPE=Debug \
      -DP4EST_ROOT=/home/dinglin/Public/p4est-2.8/opt
cmake --build build -j4

# W1
cd build && ctest --output-on-failure -R test_unit_converter

# W2/W3（单文件：冒烟 + OpenLB L2 验收）
cd build && ctest --output-on-failure -R test_cavity3d_serial

# Example + VTK
mpirun -n 1 ./build/examples/cavity3d --output cavity3d_vtk
```
