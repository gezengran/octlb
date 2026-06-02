# T02 · Mesh 拓扑层 TDD

> 类型：AFK  
> 阻塞于：T01（已完成）  
> 对应 PRD：`octlb/doc/prd/octlb-framework.md`（Mesh 模块 / FacePairList、WeightedLoadBalancer 节，测试顺序 #3 #4 #5）

---

## 要做什么

完成 Mesh 模块的拓扑层，覆盖三件事：

1. **FacePairList**：通过 `p8est_iterate` face callback 一次性遍历所有面，构建同级面列表（SameLevelFaces）和粗细 hanging 面列表（CoarseFineFaces）。这是下游 GhostSchedule 和 LevelCoupler 的共同前置。
2. **WeightedLoadBalancer**：为 `OctreeForest::partition()` 提供 weight(octant)=2^level 的权重 lambda，让细层 octant 获得更高分区权重。
3. **GhostTopology 一致性验证**：验证 FacePairList 中每条跨 rank 面对引用的 ghost octant 与 p4est ghost 层记录的 owner rank 一致。

同步完成两件配套工作：
- `src/common/types.h` 追加 `enum class FaceDir`
- `OctreeForest` 追加 ghost 层管理和 `partition(weight_fn)` 重载

严格 TDD：先写红测试，再实现，再重构。

---

## 交付物

```
octlb/
├── src/
│   ├── common/
│   │   └── types.h                  # 更新：追加 FaceDir enum class
│   └── mesh/
│       ├── forest/
│       │   ├── octree_forest.h      # 更新：partition(weight_fn)、内部访问器
│       │   └── octree_forest.cpp    # 更新：ghost 层生命周期管理
│       ├── topology/
│       │   ├── CMakeLists.txt       # 源文件追加进 octlb_mesh
│       │   ├── face_pair_list.h
│       │   └── face_pair_list.cpp
│       └── load_balance/
│           ├── CMakeLists.txt       # 源文件追加进 octlb_mesh
│           └── weighted_load_balancer.h   # 仅 header，提供 weight lambda
└── tests/
    └── unit/
        └── mesh/
            ├── CMakeLists.txt       # 更新：新增三个 test target
            ├── test_face_pair_list.cpp
            ├── test_load_balancer.cpp
            └── test_ghost_topology.cpp
```

---

## 关键设计决策（已对齐）

| 项目 | 决策 |
|---|---|
| 面方向编码 | `enum class FaceDir : int { kXMin=0, kXMax=1, kYMin=2, kYMax=3, kZMin=4, kZMax=5 }`，加入 `src/common/types.h` |
| Ghost 层生命周期 | `OctreeForest` 持有 `p8est_ghost_t*`；`partition()` 完成后自动销毁旧 ghost、重建新 ghost；析构时销毁 |
| OctreeForest 内部访问器 | ⚠ **请审查**：`OctreeForest::Impl` 暴露 `p8est_t*` 和 `p8est_ghost_t*`，仅限 Mesh 模块内部使用（注释标明，不出现在 `mesh/` 对外头文件中） |
| FacePairList 构造 | `FacePairList(const OctreeForest& forest)`；通过模块内访问器拿到两个裸指针后调用 `p8est_iterate` |
| partition 签名 | `void partition(std::function<int(OctantId)> weight_fn = nullptr)`；有权重时走 `p8est_partition_ext`，无权重时走 `p8est_partition` |
| WeightedLoadBalancer | 头文件内 free function：`std::function<int(OctantId)> make_level_weight_fn(const OctreeForest&)`；不直接碰 p4est |
| is_hanging 检测 | ⚠ **请审查**：face callback 中 `side[i].is_hanging == 1` 的一侧为 4 个细格（fine），`is_hanging == 0` 的一侧为 1 个粗格（coarse）；需正确处理两侧都非 hanging 和至少一侧为 hanging 两种情况 |
| test_ghost_topology 策略 | 本地对 ghost 层交叉校验（不需要额外 MPI）；精心设计细化模式保证 hanging 面和跨 rank 面必然出现；4 rank 运行 |

---

## 接口速览

```cpp
// src/common/types.h（追加）
namespace octlb {
enum class FaceDir : int {
  kXMin = 0, kXMax = 1,
  kYMin = 2, kYMax = 3,
  kZMin = 4, kZMax = 5,
};
}  // namespace octlb

// src/mesh/topology/face_pair_list.h
namespace octlb {

struct SameLevelFace {
  OctantId local_id;
  FaceDir  dir;
  OctantId remote_id;   // ghost octant 下标（跨 rank 时）或本地 octant 下标
  int      remote_rank; // 同 rank 时等于本 rank
};

struct CoarseFineFace {
  OctantId coarse_id;
  OctantId fine_ids[4];
  FaceDir  normal;           // 从粗侧看向细侧的方向
  int      remote_ranks[4];  // fine 侧各格的 owner rank
};

class FacePairList {
 public:
  explicit FacePairList(const OctreeForest& forest);

  const std::vector<SameLevelFace>& same_level_faces() const;
  const std::vector<CoarseFineFace>& coarse_fine_faces() const;

  // 预留接口，动态 AMR 阶段实现
  void rebuild(const OctreeForest& forest);
};

}  // namespace octlb

// src/mesh/load_balance/weighted_load_balancer.h
namespace octlb {

// 返回 weight(octant) = 2^level 的权重 lambda，
// 传入 OctreeForest::partition(weight_fn)
std::function<int(OctantId)> make_level_weight_fn(const OctreeForest& forest);

}  // namespace octlb

// src/mesh/forest/octree_forest.h（变更部分）
namespace octlb {
class OctreeForest {
 public:
  // ...（原有接口不变）

  // weight_fn == nullptr 时均匀分区（p8est_partition）
  // weight_fn != nullptr 时加权分区（p8est_partition_ext）
  void partition(std::function<int(OctantId)> weight_fn = nullptr);
};
}  // namespace octlb
```

---

## ⚠ 关键分支：p8est_iterate face callback 实现

`face_pair_list.cpp` 中的 callback 是本任务最容易出错的地方，实现前请先阅读 p4est 文档中 `p8est_iter_face_info_t` 的字段说明。

关键点：
- `info->sides.elem_count` == 2（所有内部面有两侧；边界面不进入 face callback）
- `side->is_hanging == 0`：该侧是完整的单个 quadrant（`side->full`）
- `side->is_hanging == 1`：该侧是 4 个 hanging quadrant（`side->hanging`）
- 跨 rank 时 `side->full.is_ghost == 1`（或 `side->hanging.is_ghost[i] == 1`），ghost octant 下标从 ghost 层映射到 remote rank
- 粗细 hanging 面：一定是一侧 `is_hanging==0`（粗）、另一侧 `is_hanging==1`（细）；其余为同级面

---

## 验收标准

- [x] `cmake --build build -j4` 编译通过，无 warning
- [x] `test_face_pair_list`（4 rank）通过：
  - [x] 均匀细化 2 层的单位立方体中，所有面对均为 SameLevelFaces，CoarseFineFaces 为空
  - [x] 中心 octant 额外细化 1 级后，CoarseFineFaces 数量与预期 hanging 面数一致（`balance` 后全局 9 条：6 个块面 + 3 个角 hanging）
  - [x] 每条 CoarseFineFace 的 `fine_ids` 数组恰好有 4 个条目
  - [x] `FacePairList` 头文件中不出现任何 `p8est_*` 类型
- [x] `test_load_balancer`（4 rank）通过：
  - [x] 使用 `make_level_weight_fn` 后调用 `partition(weight_fn)`，各 rank 总权重差异 < 5%
  - [x] 加权分区相对均匀分区的权重不平衡度更低（细层权重 2^level 被纳入分区）
- [x] `test_ghost_topology`（4 rank）通过：
  - [x] 精心细化模式确保跨 rank 面和 hanging 面必然出现（断言 `same_level_faces().size() > 0` 且 `coarse_fine_faces().size() > 0`）
  - [x] 所有 SameLevelFaces 中 `remote_rank` 与 p4est ghost 层记录的 owner rank 一致
  - [x] 所有 CoarseFineFaces 中每个 `remote_ranks[i]` 与 p4est ghost 层记录一致
- [ ] `OctreeForest` 析构时 ghost 层正确销毁（Valgrind 无泄漏）
- [x] `OctreeForest` 内部访问器仅在 `mesh/` 模块内头文件中使用，不出现在 `solver/` 下任何文件

---

## 阻塞关系

```
T01（已完成）
└── T02（本任务）
    ├── T03 · BlockCollection<T>（solver/field，无 Mesh 依赖，可并行启动）
    └── （T04 · OpenLB 移植依赖 T03）
        └── T05 · GhostSchedule（依赖 T02 FacePairList + T03 BlockCollection）
```

---

## 构建与测试命令（给用户执行）

```bash
cd octlb
cmake -B build -DCMAKE_BUILD_TYPE=Debug \
      -DP4EST_ROOT=/home/dinglin/Public/p4est-2.8/opt
cmake --build build -j4
cd build && ctest --output-on-failure -R "test_face_pair_list|test_load_balancer|test_ghost_topology"
```
