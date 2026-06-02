# T03 · BlockCollection\<T\> TDD

> 类型：AFK  
> 阻塞于：无（T02 完成后可立即开始）  
> 对应 PRD：`octlb/doc/prd/octlb-framework.md`（Solver/field 子层，测试顺序 #6）

---

## 要做什么

实现 `solver/field/` 子层的基础组件，为后续 GhostSchedule（T05）和 TimeLoop（T06）
提供块存储与遍历原语：

1. **`BlockCollection<T>`**：泛型块容器，以 `OctantId` 为键，内部 `std::vector<T>` 保证
   O(1) 随机访问。T 通过工厂函数在构造时一次性生成，支持非默认可构造的 T（如 `ConcreteBlockLattice`）。
2. **`BlockIterator`**：遍历所有 `OctantId`（0…size-1）的轻量范围类，对 level 完全
   无感知——level 过滤是 TimeLoop 的职责，不在此层实现。

严格 TDD：先写红测试，再实现，再重构。

---

## 交付物

```
octlb/
├── src/
│   ├── CMakeLists.txt               # 更新：追加 add_subdirectory(solver)
│   └── solver/
│       ├── CMakeLists.txt           # 新建：定义 octlb_field INTERFACE target
│       └── field/
│           ├── CMakeLists.txt       # 新建：挂 INTERFACE include 目录
│           ├── block_collection.h   # 新建：BlockCollection<T> 模板定义
│           └── block_iterator.h     # 新建：BlockIterator 定义
└── tests/
    └── unit/
        ├── CMakeLists.txt           # 更新：追加 add_subdirectory(solver)
        └── solver/
            ├── CMakeLists.txt       # 新建：test_block_collection target（1 rank）
            └── test_block_collection.cpp  # 新建：TDD 测试
```

---

## 关键设计决策（已对齐）

| 项目 | 决策 |
|---|---|
| 存储结构 | `std::vector<T>`，OctantId 即下标，O(1) 随机访问 |
| 构造方式 | `BlockCollection(label n, std::function<T(OctantId)> factory)`；工厂在构造时逐 id 调用，不要求 T 默认可构造 |
| Level 过滤 | 不在 BlockCollection 内；TimeLoop 初始化时从 OctreeForest 缓存 `level→[OctantId]` 映射 |
| BlockIterator 遍历单元 | 产出 `OctantId` 值（整数），**不是** T 引用；调用方再用 `collection[id]` 取块 |
| Mesh 依赖 | BlockCollection 和 BlockIterator 均不 include 任何 p4est / Mesh 模块头文件 |
| CMake target | `octlb_field` INTERFACE 库（全为头文件模板），不依赖 `octlb_mesh`、P4est、MPI |
| 测试运行 | 单 rank（`mpirun -n 1 --oversubscribe`），复用 `tests/mpi_main.cpp` |

---

## 接口速览

```cpp
// src/solver/field/block_collection.h
namespace octlb {

/** Stores one T per local octant, keyed by OctantId (0…size-1).
 *
 *  T need not be default-constructible; factory(id) is called once
 *  per octant at construction time.
 */
template <typename T>
class BlockCollection {
 public:
  BlockCollection(label num_octants, std::function<T(OctantId)> factory);

  T& operator[](OctantId id);
  const T& operator[](OctantId id) const;

  label size() const;
};

}  // namespace octlb

// src/solver/field/block_iterator.h
namespace octlb {

/** Iterates OctantId values 0…num_octants-1.
 *
 *  Level-agnostic: does not know which octants belong to which level.
 *  TimeLoop is responsible for filtering by level using OctreeForest.
 */
class BlockIterator {
 public:
  explicit BlockIterator(label num_octants);

  // Supports range-for: yields OctantId
  class Iterator {
   public:
    explicit Iterator(OctantId current) : current_(current) {}
    OctantId operator*() const { return current_; }
    Iterator& operator++() { ++current_; return *this; }
    bool operator!=(const Iterator& other) const {
      return current_ != other.current_;
    }
   private:
    OctantId current_;
  };

  Iterator begin() const { return Iterator(0); }
  Iterator end()   const { return Iterator(num_octants_); }

 private:
  label num_octants_;
};

}  // namespace octlb
```

### 典型用法

```cpp
// 构造（lbm/ 层会传入真正的 ConcreteBlockLattice 工厂）
BlockCollection<int> col(forest.local_num_octants(),
                         [](OctantId id) { return id * 10; });

// 全量遍历
for (OctantId id : BlockIterator(col.size())) {
    process(col[id]);
}

// TimeLoop 的 level 过滤（在 lbm/ 层实现，不在 field/ 层）
// std::vector<std::vector<OctantId>> level_octants_;  // 初始化时从 OctreeForest 缓存
// for (OctantId id : level_octants_[l]) { col[id].collideAndStream(); }
```

---

## 验收标准

- [x] `cmake --build build -j4` 编译通过，无 warning
- [x] `test_block_collection`（1 rank）通过：
  - [x] `BlockCollection<int>` 以工厂 `[](OctantId id) { return id * 10; }` 构造后，`collection[id] == id * 10` 对所有 id 成立
  - [x] `BlockCollection<std::string>` 以工厂 `[](OctantId id) { return std::to_string(id); }` 构造后，`collection[id] == std::to_string(id)` 对所有 id 成立（验证非 trivial T 可用）
  - [x] `collection.size()` 返回构造时指定的 `num_octants`
  - [x] `BlockIterator(n)` 遍历产出恰好 n 个 OctantId，值为 0, 1, …, n-1，顺序正确
  - [x] `BlockIterator(0)` 立即结束，无迭代
  - [x] `block_collection.h` 和 `block_iterator.h` 中不出现任何 p4est / Mesh 模块类型
- [x] `octlb_field` target 不依赖 `octlb_mesh`、P4est、MPI（CMake link 无这些依赖）

---

## 阻塞关系

```
T02（已完成）
└── T03（本任务）← 无前置依赖，可立即开始
    ├── T04 · OpenLB 最小移植（ConcreteBlockLattice + BGK，依赖 T03）
    └── T05 · GhostSchedule<T>（依赖 T02 FacePairList + T03 BlockCollection）
        └── T06 · LevelCoupler + TimeLoop（依赖 T04 + T05）
```

T07（stl_reader + GeometryEngine + VTK Writer）与 T03–T06 可并行推进，
仅需在集成测试前完成。

---

## 构建与测试命令（给用户执行）

```bash
cd octlb
cmake -B build -DCMAKE_BUILD_TYPE=Debug \
      -DP4EST_ROOT=/home/dinglin/Public/p4est-2.8/opt
cmake --build build -j4
cd build && ctest --output-on-failure -R "test_block_collection"
```
