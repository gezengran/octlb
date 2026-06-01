# T01 · OctreeForest TDD

> 类型：AFK  
> 阻塞于：T00（已完成）  
> 对应 PRD：`octlb/doc/prd/octlb-framework.md`（Mesh 模块 / OctreeForest 节，测试顺序 #0）

---

## 要做什么

实现 `OctreeForest`：封装 p4est `p8est_t` forest，提供 refine / balance / partition
操作，向上层暴露纯拓扑接口（不暴露任何 `p8est_*` 类型）。同步完成三件配套工作：

1. **CMake 集成**：`cmake/FindP4est.cmake` + 根 `CMakeLists.txt` 引入 MPI
2. **基础类型**：`src/common/types.h`（label、scalar）和 `src/common/bounding_box.h`
3. **MPI 测试基础设施**：`tests/mpi_main.cpp`，使后续所有 MPI 测试复用同一入口

严格 TDD：先写红测试，再实现，再重构。

---

## 交付物

```
octlb/
├── cmake/
│   └── FindP4est.cmake              # 暴露 P4est::p4est imported target
├── src/
│   ├── CMakeLists.txt               # 定义 octlb_mesh STATIC target
│   ├── common/
│   │   ├── types.h                  # label, scalar, OctantId（OpenFOAM 风格）
│   │   └── bounding_box.h           # BoundingBox struct
│   └── mesh/
│       ├── CMakeLists.txt           # 将 forest/ 源文件追加进 octlb_mesh
│       └── forest/
│           ├── octree_forest.h
│           └── octree_forest.cpp
├── tests/
│   ├── mpi_main.cpp                  # MPI_Init → RUN_ALL_TESTS → MPI_Finalize
│   └── unit/
│       └── mesh/
│           ├── CMakeLists.txt       # 更新：新增 test_octree_forest target
│           └── test_octree_forest.cpp
└── CMakeLists.txt                   # 更新：find_package(MPI/P4est)、add_subdirectory(src)
```

`src/` 目录结构遵照 PRD，不新建未经批准的顶层目录。

---

## 关键设计决策（已对齐）

| 项目 | 决策 |
|---|---|
| p4est CMake 集成 | `cmake/FindP4est.cmake`，暴露 `P4est::p4est` imported target；跨机器用 `-DP4EST_ROOT=...` |
| MPI | 根 `CMakeLists.txt` 显式 `find_package(MPI REQUIRED)`；`P4est::p4est` INTERFACE 依赖 `MPI::MPI_C` |
| 生产库粒度 | 单一 `octlb_mesh` STATIC；各子目录 CMakeLists 向其追加源文件 |
| 基础类型 | `src/common/types.h`：`using label = int32_t; using scalar = double;`（OpenFOAM 风格） |
| `OctantId` | `using OctantId = label`（本地 quadrant 线性序号，0…local_num_quadrants-1） |
| `quadrant_bounds` 返回值 | `BoundingBox`（`src/common/bounding_box.h`）；内部把 p4est 整数坐标缩放到物理域 |
| `OctreeForest` 构造函数 | `OctreeForest(MPI_Comm comm, BoundingBox domain)` |
| `refine()` 签名 | `void refine(std::function<bool(OctantId)> criterion, int max_level)` |
| MPI 测试运行 | `tests/mpi_main.cpp` 手写 init/finalize；CMake 用 `MPIEXEC_EXECUTABLE --oversubscribe -n 4` |

---

## 接口速览

```cpp
// src/common/types.h
namespace octlb {
using label  = int32_t;
using scalar = double;
using OctantId = label;
}  // namespace octlb

// src/common/bounding_box.h
namespace octlb {
struct BoundingBox {
  scalar x_min, y_min, z_min;
  scalar x_max, y_max, z_max;
};
}  // namespace octlb

// src/mesh/forest/octree_forest.h
namespace octlb {
class OctreeForest {
 public:
  OctreeForest(MPI_Comm comm, BoundingBox domain);
  ~OctreeForest();

  // 不可复制；可移动
  OctreeForest(const OctreeForest&) = delete;
  OctreeForest& operator=(const OctreeForest&) = delete;

  void refine(std::function<bool(OctantId)> criterion, int max_level);
  void balance();
  void partition();

  label      local_num_octants() const;
  BoundingBox quadrant_bounds(OctantId id) const;
  int        quadrant_level(OctantId id) const;
};
}  // namespace octlb
```

---

## 验收标准

- [ ] `cmake -B build -DP4EST_ROOT=<path> && cmake --build build -j4` 编译通过，无 warning
- [ ] `ctest --output-on-failure` 中 `test_octree_forest` 通过（4 rank，`--oversubscribe`）
- [ ] 单位立方体均匀细化 2 层后，所有 rank 的 `local_num_octants()` 之和 = 64
- [ ] `quadrant_bounds()` 对所有 octant 返回不重叠、拼合覆盖整个 domain 的包围盒
- [ ] `quadrant_level()` 在 level-2 均匀细化后返回 2
- [ ] `balance()` 后相邻 octant 层级差 ≤ 1（p4est 2:1 balance 保证）
- [ ] `partition()` 后各 rank 的 `local_num_octants()` 差异 ≤ 1（均匀分区）
- [ ] `src/common/types.h` 和 `bounding_box.h` 不 include 任何 p4est / MPI 头文件
- [ ] `octree_forest.h` 的 public API 中不出现任何 `p8est_*` 类型

---

## 阻塞关系

```
T00（已完成）
└── T01（本任务）
    ├── T03 · FacePairList（依赖 OctreeForest）
    └── T05 · GhostTopology（依赖 OctreeForest）
```

---

## 构建与测试命令（给用户执行）

```bash
cd octlb
cmake -B build -DCMAKE_BUILD_TYPE=Debug \
      -DP4EST_ROOT=/home/dinglin/Public/p4est-2.8/opt
cmake --build build -j4
cd build && ctest --output-on-failure
```

服务器上替换 `-DP4EST_ROOT` 为实际安装路径即可。
