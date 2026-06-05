# T04 · OctLB BlockLattice + BGK（算法复用，自建数据结构）

> 类型：AFK  
> 阻塞于：T03（BlockCollection\<T\> 已完成）  
> 对应 PRD：`octlb/doc/prd/octlb-framework.md`（OpenLB 代码融入策略，测试顺序 #7）  
> 状态：**完成（测试绿）**（已按“算法复用 + 自建数据结构”落地，见下方验收）

---

## 要做什么

实现 `octlb::BlockLattice<T, DESCRIPTOR>`——OctLB 自己的**面向 octant 的格子数据结构**。

设计原则：
- **数据结构自建**：`BlockLattice` 管理内存布局（每 octant 一个实例，带 ghost halo），以 `OctantId` 为键整合进 `BlockCollection`。
- **算法复用**：直接调用 OpenLB `dynamics/collision.h` 里的 BGK/MRT kernel（纯函数模板，只依赖 `MinimalCell` concept），以及 `dynamics/lbm.h` 里的平衡态公式，**不**引入 OpenLB 的 `ConcreteBlockLattice`、`DynamicsPromise`、`SuperLattice` 等框架层。
- **依赖最小化**：只复制 OpenLB 的算法头文件（descriptor/、dynamics/算法层、core/concepts.h、utilities/数学工具），约 25 个头文件，无 `.cpp`，无 tinyxml2，无 MPI 单例。

本任务**不**实现边界条件（Bouzidi / bounce-back），留到 T09；  
本任务**不**对接 `BlockCollection<BlockLattice>` 集成，留到 T06（TimeLoop）。

---

## 交付物

```
octlb/
├── src/
│   └── solver/
│       └── lbm/
│           ├── block_lattice.h          # octlb::BlockLattice<T, DESCRIPTOR>（新建）
│           ├── block_lattice.cpp        # 非模板部分（构造/析构/stream 主循环）
│           │
│           │   ── OpenLB 算法头文件（只读，不修改，namespace olb 保持）──
│           ├── descriptor/              # D3Q19 等格子常数（c/t/invCs2），~10 个头文件
│           ├── dynamics/                # BGK/MRT kernel、平衡态公式、矩计算，~12 个头文件
│           ├── core/                    # MinimalCell concept、FieldD、Vector 等，~5 个头文件
│           ├── utilities/               # 数学工具（vectorHelpers 等），~5 个头文件
│           │
│           └── CMakeLists.txt           # 定义 octlb_lbm STATIC target
└── tests/
    └── unit/
        └── solver/
            ├── CMakeLists.txt           # 已有，追加 test_collision_bgk
            └── test_collision_bgk.cpp   # 使用 octlb::BlockLattice 的行为测试
```

> `solver/io/`、`solver/lbm/time_loop/`、`solver/lbm/unit_converter/` 均不在本任务中。

---

## 关键设计决策

| 项目 | 决策 |
|---|---|
| **融合策略** | 算法复用（Algorithm Reuse）：OpenLB 只提供无状态函数模板；数据结构、内存布局、streaming 全由 OctLB 自实现 |
| **OpenLB 依赖范围** | 仅算法头文件（descriptor/、dynamics/算法层、core/concepts.h、utilities/）；不复制 blockLattice.h / superLattice.h / case/ / optimization/ 等框架层 |
| **复制 vs 外部引用** | 约 25 个头文件复制进 `lbm/` tree，彻底切断对 `../openlb/` 外部路径的依赖 |
| **Namespace** | OpenLB 算法头文件保留 `namespace olb`；OctLB 自写代码用 `namespace octlb` |
| **DynamicsPromise / lambda 问题** | 完全绕开：不使用 `ConcreteBlockLattice` 及其 `DynamicsPromise`，不需要 `olb3D.h + olb3D.hh` umbrella |
| **CMake target** | `octlb_lbm` STATIC，只编译 `block_lattice.cpp`，无 tinyxml2，无 MPI 单例依赖 |
| **CPU 平台** | BGK kernel 是纯 C++ 函数模板，不需要 `PLATFORM_CPU_SISD` 宏 |
| **Streaming** | Pull-scheme（拉取），OctLB 自实现；ghost halo 由 `GhostSchedule<BlockLattice>`（T05）负责在 streaming 前刷新 |
| **边界条件** | T04 不实现；周期性测试直接用 halo 参与 stream，T09 加 Bouzidi |
| **与 BlockCollection 的关系** | `BlockCollection<BlockLattice>` 自然适配（T03 泛型工厂构造），T06 集成 |

---

## octlb::BlockLattice 接口设计

```cpp
// src/solver/lbm/block_lattice.h
namespace octlb {

template <typename T, typename DESCRIPTOR>
class BlockLattice {
 public:
  // 构造一个 Nx × Ny × Nz 的格子，halo_width 层 ghost cell（streaming 用）
  BlockLattice(int nx, int ny, int nz, int halo_width = 1);

  // 返回内部 Cell 代理——满足 olb::concepts::MinimalCell
  CellProxy<T, DESCRIPTOR> get(int i, int j, int k);

  // 全块碰撞（遍历内部格，调用 olb::collision::BGK::type::apply）
  void collide(T omega);

  // 全块 streaming（pull-scheme，需 ghost 已刷新）
  void stream();

  // 初始化：均匀密度 rho0 + 速度 u0
  void initialize(T rho0, const T* u0);

  // ghost halo 访问（供 GhostSchedule 读写面层）
  T* face_data(int face_dir);              // 返回 face 方向的 halo 起始指针
  const T* face_data(int face_dir) const;

  int nx() const; int ny() const; int nz() const;

 private:
  // 内存布局：[Nx+2h][Ny+2h][Nz+2h][Q]，row-major，Q 在最内层
  std::vector<T> populations_;
  int nx_, ny_, nz_, h_;
};

// Cell 代理：满足 olb::concepts::MinimalCell 和 olb::concepts::Cell
template <typename T, typename DESCRIPTOR>
struct CellProxy {
  using value_t      = T;
  using descriptor_t = DESCRIPTOR;

  T& operator[](unsigned iPop) { return data_[iPop]; }
  const T& operator[](unsigned iPop) const { return data_[iPop]; }

  // FieldD / getField 接口（供 concepts::Cell 使用）
  template <typename FIELD>
  auto getField() const;

  template <typename FIELD>
  auto* getFieldPointer();

  void computeRhoU(T& rho, T* u) const;
  void defineRhoU(T rho, const T* u);
  void iniEquilibrium(T rho, const T* u);

 private:
  T* data_;   // 指向 populations_ 中该 cell 的起始位置
};

}  // namespace octlb
```

---

## Params 设计（满足 olb::concepts::Parameters）

BGK kernel 通过 `parameters.template get<descriptors::OMEGA>()` 获取 omega。
`OctParams` 是一个轻量包装，满足 `olb::concepts::Parameters`：

```cpp
// src/solver/lbm/block_lattice.h（内嵌）
template <typename T, typename DESCRIPTOR>
struct OctParams {
  using value_t      = T;
  using descriptor_t = DESCRIPTOR;

  template <typename FIELD> bool provides() const;
  template <typename FIELD> auto get() const;  // OMEGA → T
  template <typename FIELD> void set(auto val);
};
```

---

## 测试接口速览

```cpp
// tests/unit/solver/test_collision_bgk.cpp
#include "block_lattice.h"   // octlb::BlockLattice
// 无需 olb3D.h / olb3D.hh / ConcreteBlockLattice

using T = double;
using D = olb::descriptors::D3Q19<>;
using Lattice = octlb::BlockLattice<T, D>;

constexpr int N = 8;
Lattice lattice(N, N, N, /*halo=*/1);
lattice.initialize(/*rho0=*/1.0, /*u0=*/{0, 0, 0});

for (int t = 0; t < steps; ++t) {
    lattice.collide(/*omega=*/1.4);
    lattice.stream();    // 周期性 halo（测试自己填 ghost）
}

// 遍历格点，积分质量、动量，断言守恒
```

---

## CMake 配置

```cmake
# src/solver/lbm/CMakeLists.txt
add_library(octlb_lbm STATIC
    block_lattice.cpp
)

target_include_directories(octlb_lbm
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(octlb_lbm PUBLIC cxx_std_20)

# 不需要 PLATFORM_CPU_SISD（不用 OpenLB 的 dispatch 框架）
# 不需要 OLB_VERSION（不调 olbInit）
# 不需要 tinyxml2

target_link_libraries(octlb_lbm
    PUBLIC MPI::MPI_CXX
)
```

`src/solver/CMakeLists.txt` 中 `add_subdirectory(lbm)` 保持不变。

---

## 验收标准

- [x] `cmake --build build -j4` 编译通过，无 warning
- [x] `test_collision_bgk`（1 rank）通过：
  - [x] N=8 均匀初始条件，运行 100 步后总质量变化 < 1e-10
  - [x] 同等条件下 x/y/z 方向总动量守恒误差 < 1e-10
  - [x] 单格初速扰动，经足够步数后速度场趋向均匀（Maxwell 平衡态）
- [x] `octlb_lbm` 只编译 `block_lattice.cpp`，不编译任何 OpenLB `.cpp` 文件
- [x] `lbm/` 下的 OpenLB 算法头文件均在 `namespace olb`，OctLB 代码均在 `namespace octlb`
- [x] `lbm/CMakeLists.txt` 不引用 `../openlb/` 外部路径
- [x] `solver/io/` 目录本任务不触碰

---

## 阻塞关系

```
T02（FacePairList，已完成）
T03（BlockCollection<T>，已完成）
└── T04（本任务）← 可立即开始
    ├── T05 · GhostSchedule<T>（依赖 T02 + T03，与 T04 并行）
    └── T06 · LevelCoupler + TimeLoop（依赖 T04 + T05）
            └── T07 · Mesh 几何链 → T08 VTK → T09 Bouzidi（见 `T07-geometry-mesh.md`）
```

---

## 构建与测试命令（给用户执行）

```bash
cd octlb
cmake -B build -DCMAKE_BUILD_TYPE=Debug \
      -DP4EST_ROOT=/home/dinglin/Public/p4est-2.8/opt
cmake --build build -j4

cd build && ctest --output-on-failure -R "test_collision_bgk"
```

---

## 历史记录：原"照搬"方案与废弃原因

> 保留此节作为技术备忘，供后续遇到类似问题时参考。

### 原方案概述

T04 第一版试图将 `openlb/src/` **整体复制**进 `lbm/`，使用 `olb::ConcreteBlockLattice`
作为格子存储，通过 `olb3D.h + olb3D.hh` umbrella 头文件触发整体编译。

### 发现的五个问题（均已在第一版中修复，但暴露了方向性问题）

| 问题 | 根因 | 第一版修复 |
|---|---|---|
| P1 `PLATFORM_CPU_SISD` 未定义 | OpenLB Makefile 自动注入，CMake 须显式定义 | 加 `target_compile_definitions` |
| P2 `tinyxml2` 引用外部路径 | `xmlReader.h` 是传递依赖 | 复制到 `lbm/external/tinyxml2/` |
| P3 `core.h` 拉入 SuperLattice | `core.h` 是全量入口，无单块轻量路径 | 改用 `olb3D.h + olb3D.hh` umbrella |
| P4 复制范围超出预期（1477 文件） | `DynamicsPromise` lambda 问题强制使用 umbrella | 接受全量复制 |
| P5 测试退化为 smoke test | P1+P3 综合后果 | 修复后恢复三条行为断言 |

### 废弃根因

`DynamicsPromise` 构造函数体内的局部 lambda 类型必须与 `std::function<>` 实例化在同一翻译单元；
OpenLB 的解法是 `olb3D.h + olb3D.hh`，这意味着**用一个库就必须拖入整个应用框架**
（`optimization/`、`particles/`、`reaction/`、`uq/` 等与 BGK 无关的 1477 个文件）。
这与"按需融合"的设计原则冲突，因此废弃，转向方案 C（算法复用，自建数据结构）。

### 第一版验收状态（历史）

- [x] 编译通过（使用 olb3D.h + olb3D.hh umbrella）
- [x] `ConservesMassForUniformState` 通过（误差 < 1e-10）
- [x] `ConservesMomentumXYZForUniformState` 通过（误差 < 1e-10）
- [x] `RelaxesSingleCellPerturbationToEquilibrium` 通过（400 步后 spread < 1e-5）
