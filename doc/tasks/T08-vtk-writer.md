# T08 · VTK StructuredGrid Writer（`.vts` / `.vtm` / `.pvd`）

> 类型：AFK  
> 阻塞于：T03（`BlockCollection<T>`）、T01（`OctreeForest::quadrant_bounds`，W2 正式驱动）  
> 对应 PRD：`octlb/doc/prd/octlb-framework.md`（Solver/io 子层，测试顺序 #11；用户故事 #8）  
> 状态：**已完成**（W1–W4；`test_vtk_writer` / `test_vtk_writer_two_rank` 本地 ctest 绿，2026-06-07 复审）

---

## 要做什么

在 `solver/io/vtk_writer/` 实现 **AMR 并行 VTK 输出**，对齐 PRD 与用户故事 #8：每个 octant 一块 **StructuredGrid**（`.vts`），rank 0 写 **`.vtm`** 多块索引与 **`.pvd`** 时间序列。

**不**沿用 OpenLB `blockVtkWriter3D` 的 `ImageData`（`.vti`）路径；仅参考 `SuperVTMwriter3D` 的 `.vtm`/`.pvd` 编排思路，驱动层改为 `BlockCollection` + `OctreeForest::quadrant_bounds`。

内部分 **W1–W4** 波次（可分 PR merge）：

1. **W1 · 单块 `.vts` XML 内核**  
   给定块元数据 + 一组 `VtkCellField3D` 字段，写出合法 `VTKFile type="StructuredGrid"`；无 Mesh/LBM 依赖。

2. **W2 · AMR 块驱动 + `test_vtk_writer`（PRD #11）**  
   遍历本 rank 全部 `octant_id`，按 `quadrant_bounds` 与块内 `N` 调用 W1；每块一个 `.vts`。**W2 完成即 #11 可绿。**

3. **W3 · 并行 `.vtm` + `.pvd`**  
   各 rank 写本地 `.vts`；rank 0 汇总当步 `.vtm` 并维护 `.pvd` 时间序列（用户故事 #8 补齐）。

4. **W4 · LBM 场适配（不挡 T08 合并）**  
   在 `solver/lbm/` 提供 `velocity` / `pressure` / `density` 等 `VtkCellField3D` 薄适配，供 T10+ 集成算例调用；**不**将 `octlb_lbm` 链入 `octlb_io`。

本任务**不**实现：`gnuplot`（单独后续任务）、HDF5 checkpoint、域边界/solid 遮罩（T09 后由调用方决定采样）、OpenLB `SuperF3D` / `CuboidDecomposition` 驱动。

---

## 交付物

```
octlb/
├── src/
│   └── solver/
│       ├── io/
│       │   ├── CMakeLists.txt              # octlb_io
│       │   └── vtk_writer/
│       │       ├── vtk_cell_field.h        # VtkCellField3D concept
│       │       ├── vtk_block_meta.h        # bounds + nx,ny,nz + octant_id
│       │       ├── structured_grid_writer.h
│       │       ├── structured_grid_writer.cpp
│       │       ├── amr_vtk_writer.h        # W2：forest + BlockCollection 驱动
│       │       ├── amr_vtk_writer.cpp
│       │       ├── parallel_vtk_index.h    # W3：CollectVts / WriteVtm / AppendPvd
│       │       ├── parallel_vtk_index.cpp
│       │       ├── vtk_binary_codec.h        # base64 Float64 编解码
│       │       └── vtk_binary_codec.cpp
│       └── lbm/
│           └── vtk_lbm_fields.h            # W4：BlockLattice 适配（可选 .cpp）
└── tests/
    └── unit/solver/
        ├── test_vtk_writer.cpp             # PRD #11（W2 主验收）
        └── vtk_writer_fixtures.h           # DummyCellField + 假 bounds
```

CMake：`octlb_io` 依赖 `octlb_field`、`octlb_mesh`（W2 起）、`MPI::MPI_CXX`；**不**依赖 `octlb_lbm`。`test_vtk_writer` 链 `octlb_io`；W4 适配器仅在集成目标或独立小测试中链 `octlb_lbm`。

---

## 关键设计决策（已对齐）

| 项目 | 决策 |
|------|------|
| **VTK 类型** | `StructuredGrid`（`.vts`），**非** OpenLB 默认的 `ImageData`（`.vti`） |
| **多块索引** | 每时间步 `.vtm` + 时间序列 `.pvd`（rank 0） |
| **字段位置** | **CellData**（cell-centered，与 LBM interior 一致） |
| **几何** | 块内 `N×N×N` LBM cell → `(N+1)³` **角点** `Points` + `N³` hexahedron cells；角点由 `quadrant_bounds` **线性插值**：`x(i)=x_min+(x_max-x_min)*i/N`（y、z 同理） |
| **Extent** | `WholeExtent` / `Piece Extent` 为点索引：`0 N  0 N  0 N` |
| **采样范围** | 仅 interior `(i,j,k)∈[0,N)`，**不含** ghost halo |
| **字段抽象** | **`VtkCellField3D` concept**（对齐 T05 `FacePackable`）；单块 `.vts` 内可写**多个** `CellData` 数组 |
| **io 与 LBM 解耦** | `vtk_writer/` 不 include LBM 头文件；W4 适配在 `solver/lbm/` |
| **编码** | 第一版 **binary**（base64，`byte_order="LittleEndian"`）；ASCII 为 P2 |
| **OpenLB 移植** | **不**复制 `blockVtkWriter3D`；W3 参考 `SuperVTMwriter3D` 的 vtm/pvd 链接逻辑，去掉 `CuboidDecomposition` / `SuperF3D` |

---

## 类型与 API（示意）

```cpp
// vtk_cell_field.h — 无 mesh / LBM
template <typename T>
concept VtkCellField3D = requires(const T& f, int i, int j, int k, double* out) {
  { f.vtk_name() } -> std::convertible_to<std::string_view>;
  { f.vtk_components() } -> std::same_as<int>;  // 1=scalar, 3=vector
  { f.sample_cell(i, j, k, out) } -> std::same_as<void>;  // i,j,k ∈ [0, N)
};

struct VtkBlockMeta {
  OctantId id;
  int nx, ny, nz;
  BoundingBox bounds;  // OctreeForest::quadrant_bounds(id)
};

// W1：单块写出
void WriteStructuredGridVts(
    const std::string& filename,
    const VtkBlockMeta& meta,
    std::span<const VtkCellField3D auto* const> fields);

// W2：本 rank 全部本地块
class AmrVtkWriter {
 public:
  AmrVtkWriter(MPI_Comm comm, const OctreeForest& forest, int nx, int ny, int nz,
               std::string output_dir, std::string base_name);

  template <typename BlockT>
  std::vector<std::string> WriteTimestep(
      int iT, const BlockCollection<BlockT>& blocks,
      std::span<const VtkCellFieldView> fields);

  /** Gathers paths from all ranks; rank 0 writes `.vtm` and appends `.pvd`. */
  void WriteVtmAndPvd(int iT, const std::vector<std::string>& local_vts_paths);
};

// W3：底层索引 API（由 AmrVtkWriter::WriteVtmAndPvd 调用）
std::vector<std::string> CollectVtsRelativePaths(
    MPI_Comm comm, const std::vector<std::string>& local_absolute_paths,
    const std::string& output_dir);
void WriteVtmIndex(MPI_Comm comm, int iT, const std::string& output_dir,
                   const std::string& base_name,
                   const std::vector<std::string>& vts_relative_paths);
void AppendPvdTimestep(int iT, const std::string& pvd_path,
                       const std::string& vtm_relative_path);
```

**Cell 索引顺序**：VTK StructuredGrid 约定 x 最快，再 y，再 z；`sample_cell(i,j,k)` 与 LBM interior 下标一致。

**字段传入**：实现使用 type-erased `VtkCellFieldView::From(field)` 传入 `std::span<const VtkCellFieldView>`，支持单步写出多个 `VtkCellField3D` 场。

**文件命名**：`{base_name}_r{rank}_oct{local_octant_id}_T{iT:05d}.vts`（`OctantId` 为 rank 内本地序号，并行时必须带 `rank` 前缀以免冲突）；`.vtm` 为 `{base_name}_T{iT:05d}.vtm`，`.pvd` 为 `{base_name}.pvd`，均置于 `output_dir` 下并以相对路径互相引用。

---

## 波次与纵向切片

| 波次 | 标题 | 类型 | 阻塞于 | 用户故事 | 端到端交付 |
|------|------|------|--------|----------|------------|
| **W1** | 单块 StructuredGrid `.vts` 内核 | AFK | T03（concept 可独立） | — | 常量 `DummyCellField` → 合法 `.vts` XML |
| **W2** | AMR 驱动 + `test_vtk_writer` | AFK | W1、T01、T03 | #8（几何部分） | 多 octant + 假 bounds → 每块 `.vts` 与 bounds 一致 |
| **W3** | 并行 `.vtm` / `.pvd` | AFK | W2 | #8（完整） | 2 rank 各 1 块 → rank0 `.vtm` 含 2 条 `DataSet` |
| **W4** | LBM 场适配 u/p/ρ | AFK | W2、T04 | #8（流场） | `BlockLattice` 驱动写出（供 T10+；不挡 #11） |

---

## 测试决策

**原则**：只测对外行为；#11 用硬编码 `BlockCollection` + 极简 forest fixture（`OctreeForest(domain, 2,1,1)` brick → 1 rank 恰好 2 octant；`partition()` → 2 rank 各 1 octant），**不**依赖真实 STL 或大算例。

### W1 · 内核（可合入 `test_vtk_writer` 或独立 `test_structured_grid_writer`）

| 用例 | 设置 | 断言 |
|------|------|------|
| `SingleBlock_ConstantScalar` | `N=2`，1 个标量 `DummyCellField` | XML 含 `StructuredGrid`；`NumberOfPoints=(N+1)³`；`NumberOfCells=N³` |
| `SingleBlock_CornerBoundsMatch` | 已知 `bounds`，`N=4` | 解析 `Points`：min/max 与 `x_min/x_max` 等一致（容差 ε） |
| `SingleBlock_MultiField` | 标量 + 3 分量矢量 | 两个 `CellData` 数组，`Name` 与 `vtk_components` 一致 |

### W2 · `test_vtk_writer`（PRD #11，硬性验收）

| 用例 | 设置 | 断言 |
|------|------|------|
| `TwoOctants_DistinctBounds` | 1 rank，2 个 octant，不同 `quadrant_bounds`，共享 `N` | 输出 2 个 `.vts`；各自角点范围与对应 bounds 一致 |
| `CellData_CountMatchesN` | `N=8` | 每数组 `NumberOfTuples=N³` |
| `InteriorOnly_NoGhost` | 场值含 halo 标记 | 输出仅反映 `sample_cell(0…N-1)`，不测 ghost |

### W3 · P1（不挡 W2 / #11 合并）

| 用例 | 设置 | 断言 |
|------|------|------|
| `TwoRank_VtmListsAllBlocks` | 2 rank，各 1 octant | rank0 `.vtm` 引用 2 个 `.vts`；`.pvd` 至少 1 个 timestep 条目 |

### W4 · P1（T10+ 前文档记录）

- `BlockLattice_VelocityField_ParaViewLoads`：单块写出后 ParaView 可开（人工或集成阶段验收）

### 运行配置

- #11 默认 **1 rank**
- W3 P1：2 rank `--oversubscribe`
- 复用 `tests/mpi_main.cpp`

---

## 验收标准

### W1 + W2（T08 核心，对应 PRD #11）

- [x] `cmake --build build -j4` 通过
- [x] `octlb_io` **不**链接 `octlb_lbm`
- [x] `vtk_writer/` 头文件 **不** include LBM / OpenLB 算法头
- [x] W1 单块 `.vts`：StructuredGrid + CellData + 角点 bounds 断言通过
- [x] `test_vtk_writer` W2 上表 **全部**通过（PRD #11）
- [x] T01–T07 相关 ctest **仍绿**（T08 合入前基线；合入后由用户本地验证）

### W3（用户故事 #8 补齐）

- [x] 2-rank P1：`TwoRank_VtmListsAllBlocks` 通过（`test_vtk_writer_two_rank`）
- [x] rank 0 可生成当步 `.vtm` 与可追加的 `.pvd`

### W4（可选，不挡 T08 合并）

- [x] `vtk_lbm_fields.h` 提供 velocity / pressure / density 适配
- [x] 文档说明 T10+ 集成调用方式（见下「W4 集成示例」）

---

## 复审结论（2026-06-07）

对照 PRD 用户故事 #8 与测试 #11，实现与验收**已对齐**；复审中发现并修复一项并行文件名冲突（见下）。

| 项 | 结论 |
|----|------|
| W1–W3 功能与 ctest | 与任务表一致；用户本地 `ctest -R test_vtk_writer` 全绿 |
| io / LBM 解耦 | `octlb_io` 仅链 `octlb_field`、`octlb_mesh`、MPI；`vtk_lbm_fields.h` 在 `solver/lbm/` |
| 并行 `.vts` 命名 | **已修复**：原 `{base}_oct{local_id}` 在 partition 后各 rank 均为 `oct0`，会覆盖；现改为 `{base}_r{rank}_oct{local_id}`；W3 测试增加「2 个不同 `.vts` 路径且文件存在」断言 |
| 多 timestep `.pvd` | `AppendPvdTimestep` 已实现，**未**单测多步追加；T10+ 集成时补测或人工 ParaView 验收 |
| W4 ParaView | `BlockLattice_VelocityField_ParaViewLoads` 仍留 T10+ 集成阶段人工/集成验收 |
| ASCII 编码 | P2，未实现（与决策表一致） |
| `BlockIterator` | PRD 驱动层描述为 `BlockIterator`；实现用 `local_num_octants` 下标循环，语义等价 |

### W4 集成示例（T10+）

```cpp
#include "src/solver/lbm/vtk_lbm_fields.h"
#include "src/solver/io/vtk_writer/amr_vtk_writer.h"

// 每步 collide/stream 之后：
AmrVtkWriter writer(comm, forest, N, N, N, output_dir, "cylinder");
const auto& lattice = blocks[octant_id];  // 在循环外按块构造 field 亦可
VtkVelocityField vel(lattice);
VtkPressureField pres(lattice);
const VtkCellFieldView fields[] = {VtkCellFieldView::From(vel),
                                   VtkCellFieldView::From(pres)};
const auto paths = writer.WriteTimestep(iT, blocks, fields);
writer.WriteVtmAndPvd(iT, paths);
```

链入：`examples/*` 或集成测试 target 同时链 `octlb_io` 与 `octlb_lbm`；**不**把 LBM 链入 `octlb_io`。

---

## 阻塞关系

```
T03（BlockCollection）— 已完成
T01（OctreeForest::quadrant_bounds）— 已完成
T07（Mesh 几何链）— 测试绿；与 T08 无硬阻塞，可并行
└── T08（本任务）
    ├── W1 → W2（#11 绿）
    ├── W2 → W3（#8 完整并行输出）
    ├── W2 → W4（LBM 适配，T10+ 消费）
    ├── T09 · Bouzidi + Lattice 材料初始化（阻塞于 T07）
    └── T10+ · cavity/cylinder 集成（#12–#14；阻塞于 #0–#11 全绿）
```

---

## 构建与测试命令（给用户执行）

```bash
cd octlb
module load octlb
cmake -B build -DCMAKE_BUILD_TYPE=Debug \
      -DP4EST_ROOT=/home/dinglin/Public/p4est-2.8/opt
cmake --build build -j4

cd build && ctest --output-on-failure -R "test_vtk_writer"
```

W3 / 2-rank（`test_vtk_writer_two_rank` 仅跑 `TwoRank_VtmListsAllBlocks`）：

```bash
cd build && ctest --output-on-failure -R "test_vtk_writer_two_rank"
```

复审后若修改了 `amr_vtk_writer.cpp` 或 W3 测试，请重新 `cmake --build build -j4` 再跑上述两条 ctest。
