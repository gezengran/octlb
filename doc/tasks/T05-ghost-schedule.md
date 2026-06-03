# T05 · GhostSchedule\<T\> + FaceIterator + BlockLattice 面层

> 类型：AFK  
> 阻塞于：T02（FacePairList）、T03（BlockCollection\<T\>）、T04（BlockLattice，L2 测试）  
> 对应 PRD：`octlb/doc/prd/octlb-framework.md`（Solver/field 子层，测试顺序 #8）  
> 状态：**完成（测试绿）**

---

## 要做什么

在 `solver/field/` 实现同级块 **ghost halo 交换** 与 **粗细面遍历原语**，并扩展 Mesh / LBM 的最小接口：

1. **`FacePackable` concept + `GhostSchedule<T>`**  
   从 `FacePairList::SameLevelFaces` 在**构造时**固化通信计划；`exchange()` 执行 pack →（同 rank 直写 | 跨 rank MPI）→ unpack。  
   **不**处理 `CoarseFineFaces`（粗细搭接由 T06 `LevelCoupler` / Lagrava 负责）。

2. **`FaceIterator`**  
   遍历 `CoarseFineFaces`，供 T06 `LevelCoupler` 消费；T05 仅实现迭代器与条数/字段合法性测试。

3. **`SameLevelFace::comm_tag`（T02 小扩展）**  
   在 `FacePairList` face callback 内为每条同级面生成**对称** MPI tag，供 `GhostSchedule` 直接使用（避免 Schedule 构造期 Allgatherv 配对）。

4. **`BlockLattice::pack_face` / `unpack_face`**  
   满足 `FacePackable`；语义：pack = interior 最外一层，unpack = 紧贴 interior 的 ghost 第一层（D3Q19 单步 stream 每层 exchange 刷新 1 层）。

本任务**不**实现 `LevelCoupler`、TimeLoop、域边界 Bouzidi（T06/T07）。

---

## 交付物

```
octlb/
├── src/
│   ├── mesh/topology/
│   │   ├── face_pair_list.h      # SameLevelFace 增加 comm_tag
│   │   └── face_pair_list.cpp    # 对称 comm_tag 生成
│   └── solver/field/
│       ├── face_packable.h       # FacePackable concept
│       ├── ghost_schedule.h      # GhostSchedule<T>
│       ├── face_iterator.h       # FaceIterator
│       └── CMakeLists.txt        # octlb_field + octlb_field_schedule
└── tests/unit/solver/
    ├── test_ghost_schedule.cpp   # L1 DummyBlock 全集 + L2 BlockLattice×1
    └── CMakeLists.txt            # 2-rank / 4-rank MPI 用例
```

---

## 关键设计决策（已对齐）

| 项目 | 决策 |
|---|---|
| **与粗细搭接分工** | `GhostSchedule` 仅 `SameLevelFaces`（1:1 同级块）；1 粗面 : 4 细面走 `CoarseFineFaces` + T06 Lagrava，**不能**用多层 halo slab 替代 |
| **面层语义** | `pack_face(dir)` 读 interior 最外层；`unpack_face(dir, buf)` 写 `dir` 侧紧贴 interior 的 **1 层** ghost；buffer = N×N×Q，与 `h` 总深度无关 |
| **T 的接口** | **方案 B**：`FacePackable` concept 要求 `pack_face` / `unpack_face`（及 `face_buffer_count`）；`field/` 不 include LBM 头文件 |
| **同 rank** | `remote_rank == my_rank`：`pack` 后 `blocks[remote_id].unpack_face(opposite(dir), buf)`，不走 MPI |
| **跨 rank** | 使用 `SameLevelFace::comm_tag`；send 与 recv 成对；`remote_id` 为 ghost 索引，**不**用于本地块下标 |
| **通信计划** | **8a**：构造时固化 `entries_` 与 buffer；`FacePairList::rebuild()` 后须**重建** `GhostSchedule` |
| **MPI tag** | **9a**：`comm_tag` 在 Mesh 侧 face callback 中对称生成；无 Schedule 构造期配对 Allgatherv |
| **CMake** | **7b**：`octlb_field` 保持无 mesh/MPI；新建 `octlb_field_schedule` INTERFACE → `octlb_field` + `octlb_mesh` + `MPI` |
| **每步顺序** | `collide()` → `GhostSchedule::exchange()` → `stream()`（与 T04 `block_lattice` 注释一致） |
| **域边界** | `tree_boundary` 面不进入 `SameLevelFaces`；外侧 ghost 不由本任务填充（T07 BC） |

---

## FacePackable 与 GhostSchedule

```cpp
// face_packable.h — 仅 concept，无 mesh / LBM
template <typename T>
concept FacePackable = requires(const T& ct, T& t, FaceDir dir, T* buf, int n) {
  { ct.pack_face(dir, buf, n) } -> std::same_as<void>;
  { t.unpack_face(dir, buf, n) } -> std::same_as<void>;
  { T::face_buffer_count(nx, ny, nz) } -> std::same_as<int>;
};

// ghost_schedule.h
template <FacePackable T>
class GhostSchedule {
 public:
  GhostSchedule(MPI_Comm comm,
                const FacePairList& faces,
                BlockCollection<T>& blocks,
                int nx, int ny, int nz);

  void exchange();

 private:
  struct Entry {
    OctantId local_id;
    FaceDir dir;
    OctantId remote_id;   // 同 rank：本地 OctantId；跨 rank：不用于 unpack 目标
    int remote_rank;
    int comm_tag;         // 跨 rank；同 rank 忽略
  };
  // entries_、send_bufs_、MPI 请求等构造时分配
};
```

`exchange()` 推荐顺序：先 posting 全部跨 rank `Irecv`，再 `Isend` + 同 rank 直写，最后 `Waitall` 并 `unpack` 本 rank ghost。

---

## SameLevelFace 扩展（comm_tag）

在 `SameLevelFace` 增加 `int comm_tag`，于 `FaceCallback` 推送同级面时由 **p8est face 几何信息** 生成对称 tag（同一物理面两侧 rank 得到相同值）。  
更新 `test_face_pair_list` / `test_ghost_topology`：跨 rank 面 `comm_tag` 一致、本地条目互异（若适用）。

---

## BlockLattice 面层（L2）

在 `block_lattice.h/.cpp` 增加：

- `static int face_buffer_count(int nx, int ny, int nz);` — 按 `FaceDir` 返回 N×N×Q 等
- `void pack_face(FaceDir dir, T* buffer, int count) const;`
- `void unpack_face(FaceDir dir, const T* buffer, int count);`

索引与 `fill_periodic_halo()` 一致，但 **pack 源 = interior 边界**，**unpack 目标 = ghost 第一层**（`h=1` 时与 periodic 拷贝方向相同、角色相反）。  
`h>1` 时仍只 exchange **1 层**；`fill_periodic_halo()` 若需支持 `h>1` 可单独修，不挡 T05。

---

## FaceIterator

```cpp
class FaceIterator {
 public:
  explicit FaceIterator(const FacePairList& faces);
  // 迭代 CoarseFineFace：coarse_id, fine_ids[4], normal, remote_ranks[4]
};
```

无 MPI、无 LBM；T05 只测迭代条数与 `fine_ids` 合法性（复用 T02 中心细化 fixture）。

---

## CMake

```cmake
# solver/field/CMakeLists.txt
# octlb_field — 不变（block_collection, block_iterator, face_packable）

add_library(octlb_field_schedule INTERFACE)
target_link_libraries(octlb_field_schedule INTERFACE
  octlb_field
  octlb_mesh
  MPI::MPI_CXX
)
```

- `test_block_collection` → 仅 `octlb_field`
- `test_ghost_schedule` → `octlb_field_schedule`（+ L2 时 `octlb_lbm`）

---

## 测试决策

**原则**：只测对外行为；L1 用 `DummyBlock` 隔离 MPI 与拓扑；L2 验证真实 `pack_face` / `unpack_face`。

### L1 · 硬性验收（`DummyBlock`，链 `octlb_field_schedule`）

| 用例 | 设置 | 断言 |
|---|---|---|
| `OneRank_TwoAdjacentBlocks` | 1 rank，2 块，1 条同 rank 面 | `exchange()` 后双方 ghost 与对方 pack 模式一致 |
| `TwoRank_OneBlockEach` | 2 rank，各 1 块，1 条跨 rank 面 | 两侧 ghost 与对方 interior pack 一致 |
| `EmptySchedule_NoOp` | 1 rank，1 块，无同级邻面 | `exchange()` 不崩溃、数据不变 |
| `MixedLocalAndRemoteFaces` | ≥2 rank，单 rank 上既有同 rank 又有跨 rank 邻块 | 一次 `exchange()` 后两类面均正确 |
| `CornerGhostConsistent` | 1 rank，2×2（或 2×2×1）块网格 | 角点 ghost 与多面 exchange 后的期望值一致 |
| `CommTagUnique` | 多 rank + 多跨 rank 面（如 4 rank 细化） | 本 rank 跨 rank 条目 `comm_tag` 两两不同 |
| `OppositeDirPairing` | 2 rank 一对邻块 | A 的 `kXMax` 与 B 的 `kXMin` 使用相同 `comm_tag` |
| `FaceIterator_MatchesCoarseFineCount` | T02 中心细化 forest | 迭代数 = `coarse_fine_faces().size()`，`fine_ids[4]` 合法 |

### L2 · 硬性验收（1 个，`BlockLattice`，链 `octlb_lbm`）

| 用例 | 设置 | 断言 |
|---|---|---|
| `TwoRank_BlockLattice_FaceValuesMatch` | 2 rank，各 1 块 `BlockLattice`，初始化可识别 f 模式 → `collide`（可选）→ `exchange()` | 邻接面 ghost 与对方 interior 边界 populations 一致（逐分量或 checksum） |

### P1（文档记录，不挡 T05 合并）

- `DomainTreeBoundary_NoSameLevelFace`：域外树边界无同级条目
- `Rebuild_RecreateSchedule`：`pairs.rebuild()` 后新建 `GhostSchedule` 仍正确

### 运行配置

- L1 多数用例：2 rank `--oversubscribe`；`CommTagUnique` / `MixedLocalAndRemoteFaces` 可用 4 rank
- `FaceIterator` / `EmptySchedule`：1 rank 即可
- 复用 `tests/mpi_main.cpp`

---

## 验收标准

- [x] `cmake --build build -j4` 通过
- [x] T02 现有 mesh 测试仍绿（`comm_tag` 扩展后）
- [x] `test_ghost_schedule` L1 上表 **全部**通过
- [x] `test_ghost_schedule` L2 **1 例**通过
- [x] `octlb_field` 仍不链接 `octlb_mesh` / MPI
- [x] `ghost_schedule.h` / `face_iterator.h` 不 include LBM 头文件
- [x] `FaceIterator` 不实现粗细耦合 MPI

---

## 阻塞关系

```
T02（FacePairList）— 本任务扩展 comm_tag
T03（BlockCollection）— 已完成
T04（BlockLattice）— 已完成；本任务扩展 pack/unpack
└── T05（本任务）
    └── T06 · LevelCoupler + TimeLoop（依赖 T04 + T05）
```

---

## 构建与测试命令（给用户执行）

```bash
cd octlb
module load octlb   # 若需要
cmake -B build -DCMAKE_BUILD_TYPE=Debug \
      -DP4EST_ROOT=/home/dinglin/Public/p4est-2.8/opt
cmake --build build -j4

cd build && ctest --output-on-failure -R "test_face_pair_list|test_ghost_topology|test_ghost_schedule"
```
