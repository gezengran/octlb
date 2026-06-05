# T07 · STL Reader + GeometryEngine + MaterialField（几何自适应加密）

> 类型：AFK  
> 阻塞于：T01（OctreeForest）、T02（FacePairList，仅 P1 冒烟）  
> 对应 PRD：`octlb/doc/prd/octlb-framework.md`（Mesh 模块，测试顺序 #1、#2a、#2；用户故事 #1）  
> 状态：**开发中**（W1–W4 代码与 #1/#2/#2a 单测已落地；待本地全量 ctest 确认）

---

## 要做什么

在 Mesh 模块实现 **STL 输入 → 几何自适应 AMR 加密 → 块内体素化 → `MaterialField`** 的端到端路径，对齐 octree-mesh 的 `GeometryAdaptiveEngine` + CGAL 体素化思路，并支持 **多部件（Part）** 组合（如风洞洞壁 + 试验件）。

1. **`stl_reader`**（W1）  
   移植 OpenLB STL 读取逻辑（ASCII/binary），产出 `TriangleSoup` 与包围盒/法向；`tests/data/mesh/` 小 STL 仅用于 #1。

2. **`GeometryAssembly` + 自适应加密**（W2）  
   - 多 `GeometryPart`（`role` + `priority`），加密判据为 **全部 part 三角面并集** 的表面附近 refine。  
   - `ResolveBounding` + `ResolveSurface`（多轮 `FastVoxelize` 式相交检测）驱动 `OctreeForest::refine` / `balance` / `partition`（`make_level_weight_fn`）。  
   - **不**在 `build()` 内构造 `FacePairList`（见下方契约）。

3. **体素化 + `MaterialField`**（W3）  
   - 按 part 角色映射 CGAL inside/outside/相交 → `fluid` / `solid` / `boundary`；**boundary** = 与三角相交的格（同 octree-mesh `CellVoxelizer`），**不**预计算 Bouzidi 距离/法向。  
   - 按 `priority` 合并：**solid > boundary > fluid**。  
   - `kInternalChannel` 第一版仅 **闭壳墙材**（CGAL outside → 洞腔 fluid）；加载时 watertight 检查。

4. **`GeometryEngine::build`**（W4）  
   编排 W2+W3；返回（或填充）`MaterialField`；API/文档强制 **调用方** 在拓扑变更后重建 `FacePairList` 及下游 Schedule/Coupler。

本任务**不**实现：VTK、`gnuplot`（**T08**）、Bouzidi / `BlockLattice` 材料写入（**T09**）、`unit_converter` 与 cavity/cylinder 集成（**T10+**）、`kFluidLumenVolume` STL、曲率 `ResolveSharp`、通用 CSG、动态 AMR 热路径、`FacePairList` 实现。

---

## 交付物

```
octlb/
├── src/
│   └── mesh/
│       ├── io/stl_reader/
│       │   ├── stl_reader.h
│       │   └── stl_reader.cpp
│       └── geometry/
│           ├── triangle_soup.h
│           ├── geometry_assembly.h
│           ├── material_field.h
│           ├── geometry_config.h
│           ├── geometry_engine.h
│           ├── geometry_engine.cpp
│           ├── geometry_adaptive_refinement.cpp   # ResolveBounding/Surface
│           └── voxelization/                    # 裁剪自 octree-mesh + CGAL
│               └── ...
├── cmake/                                       # 若需要 FindCGAL.cmake
└── tests/
    ├── data/mesh/                               # 极小 ASCII/binary STL（#1）
    └── unit/mesh/
        ├── test_stl_reader.cpp                  # PRD #1
        ├── test_geometry_adaptive_refine.cpp    # PRD #2a
        └── test_geometry_engine.cpp             # PRD #2（多场景 S1–S7）
```

CMake：几何体素化链 CGAL；`test_stl_reader` **不**链接 CGAL；MPI 加密路径与 octree-mesh 一致（`MPI_Allreduce` 保证各 rank 同步 `refine`/`balance`）。

---

## 关键设计决策（已对齐）

| 项目 | 决策 |
|------|------|
| **任务边界** | 单一 T07；内部分 W1–W4 波次，可分 PR merge |
| **多几何** | `GeometryAssembly` = `vector<GeometryPart>`；部件数不固定（1/2/N STL） |
| **风洞** | 洞壁 → `kInternalChannel`（低 `priority`）；试验件 → `kExternalObstacle`（高 `priority`） |
| **洞壁 STL** | 仅 **闭壳墙材** `kSolidWallWatertight`；腔体 = CGAL outside → fluid |
| **试验件 STL** | 实心闭体；inside → solid，outside → fluid |
| **材料合并** | 按 priority 从低到高；冲突 **solid > boundary > fluid** |
| **boundary** | 与三角 **表面相交** 的格；无 Bouzidi 距离/法向（**T09**） |
| **加密** | 全 part 表面并集；参考 octree-mesh `GeometryAdaptiveEngine` |
| **`build()` 与拓扑** | **A**：`build()` 只改 `OctreeForest` + `MaterialField`；**不**创建 `FacePairList`；调用方 rebuild（同 T06） |
| **P1 冒烟** | S5 风洞 fixture：`build` 后 `FacePairList(forest)` 同级/粗细列表非空 |
| **体素化来源** | 裁剪 octree-mesh `voxelization` + CGAL；`stl_reader` 以 OpenLB 移植为主 |
| **单测几何** | #2 / #2a 用 **硬编码 triangle soup**；不依赖大 STL 文件 |

---

## 类型与 API（示意）

```cpp
enum class GeometryPartRole : std::uint8_t {
  kInternalChannel,   // 洞壁/管道：闭壳墙材，腔体 fluid
  kExternalObstacle,  // 试验件/钝体：体内 solid，体外 fluid
};

enum class MaterialKind : std::uint8_t {
  kFluid,
  kSolid,
  kBoundary,
};

struct GeometryPart {
  TriangleSoup soup;
  GeometryPartRole role = GeometryPartRole::kExternalObstacle;
  int priority = 0;
  std::string name;  // e.g. "tunnel", "model"
};

struct GeometryAssembly {
  std::vector<GeometryPart> parts;
};

struct GeometryConfig {
  int max_level = 6;
  int resolve_surface_times = 3;
  double bound_width = 0.0;
  double wake_length = 0.0;
  int wake_direction = 0;  // 0=x, 1=y, 2=z
  int cell_width = 4;    // 块内 Nx=Ny=Nz
};

class MaterialField {
 public:
  MaterialKind at(OctantId id, int i, int j, int k) const;
  // 构造时绑定 local_num_octants 与 (nx,ny,nz)
};

class GeometryEngine {
 public:
  // 修改 forest（refine/balance/partition），产出 MaterialField。
  // 调用方须在返回后：FacePairList pairs(forest); 并重建 GhostSchedule/LevelCoupler/TimeLoop。
  MaterialField build(OctreeForest& forest,
                      const GeometryAssembly& assembly,
                      const GeometryConfig& config) const;
};
```

**Part 体素映射（共用 CGAL 分类，按 role 翻表）**

| CGAL / 体素中间态 | `kExternalObstacle` | `kInternalChannel` |
|-------------------|---------------------|---------------------|
| outside STL 体 | fluid | fluid（洞腔） |
| inside STL 体 | solid | solid（墙材） |
| 与表面相交 | boundary | boundary |

---

## 几何自适应加密（W2 行为契约）

参考 octree-mesh `GeometryAdaptiveEngine`：

1. **ResolveBounding**：扩展几何包围盒（`bound_width`、可选 `wake_length` 沿 `wake_direction`），从粗到细标记与 bbox 相交的 octant 并 `refine`，每轮后 `balance()`，全局 `MPI_Allreduce` 同步是否 refine。  
2. **ResolveSurface**：对当前最细叶层做 fast 相交检测，相交则 refine；迭代 `resolve_surface_times` 次；首轮后 `partition(make_level_weight_fn(forest))`。  
3. 加密判据几何 = **assembly 内所有 part 的三角并集**（实现可合并 soup 或逐 part 相交 OR）。

---

## `build()` 后拓扑契约

| 事件 | 调用方动作 |
|------|------------|
| `GeometryEngine::build()` 返回 | 新建 `FacePairList pairs(forest)` |
| 同上 | 重建 `GhostSchedule`、`LevelCoupler`；刷新或重建 `TimeLoop` 层缓存（T06） |
| 第一版运行期 | 初始化调用一次；无动态 AMR 热路径 |

T07 **不**实现 `FacePairList::rebuild` 新逻辑，仅消费 T02 已有 API。

---

## 测试决策

**原则**：只测对外行为；几何/scene 用 triangle soup fixture；STL 文件仅 #1。

### `test_stl_reader`（PRD #1）

| 用例 | 断言 |
|------|------|
| ASCII 小模型 | 三角面片数、法向方向、包围盒与已知值一致 |
| binary 小模型 | 同上 |

### `test_geometry_adaptive_refine`（PRD #2a）

| 用例 | 设置 | 断言 |
|------|------|------|
| `SurfaceRefinement_IncreasesLevelNearGeometry` | 单位域 + 中心小立方 soup，config 使 `max_level≥2` | 贴表面叶 octant 的 `quadrant_level` 最大值 > 远场角点叶 octant |
| `ResolveBounding_RefinesInsideExtendedBBox` | 仅 bbox 步骤或完整 build 前半 | bbox 内层数高于域角（阈值文档化） |
| **P1** `TwoRank_RefineCollective` | 2 rank，几何跨 rank | 无死锁；各 rank `local_num_octants` 合理 |

### `test_geometry_engine`（PRD #2）

| ID | 场景 | parts | 断言要点 |
|----|------|-------|----------|
| **S1** | 外流实心块 | 1× `kExternalObstacle` 立方 | 内 solid、外 fluid、表面 boundary |
| **S3** | 内流方管 | 1× `kInternalChannel` 薄壳方管 soup | 腔 fluid、壁 solid、壁面 boundary |
| **S5** | 风洞组合 | 方管 channel + 管内球 obstacle，model priority 更高 | 管腔 fluid；球内 solid；球外管腔仍 fluid |
| **S6** | 覆盖顺序 | 交换 priority 或仅 obstacle | 与 S5 对比，锁定 priority 语义 |
| **S7** | 非法洞壁 | 实心盒标为 `kInternalChannel` | `build` 失败或明确错误（watertight/腔体体积诊断） |

**P1（不挡 T07 合并）**

- `Build_ThenFacePairList_NonEmpty`：S5 上 `build` 后 `FacePairList`，`SameLevelFaces` 与 `CoarseFineFaces` 条目数 > 0（1 rank 即可）。

### 运行配置

- #1、#2、#2a 默认 1 rank；P1 两 rank 用 `--oversubscribe`  
- 复用 `tests/mpi_main.cpp`

---

## 验收标准

- [ ] `cmake --build build -j4` 通过（`module load octlb` 含 CGAL）
- [ ] `test_stl_reader` **全部**通过（PRD #1）
- [ ] `test_geometry_adaptive_refine` **至少**表面加密用例通过（PRD #2a）
- [ ] `test_geometry_engine` S1、S3、S5 **通过**；S7 非法洞壁 **失败符合预期**
- [ ] T01–T06 相关 ctest **仍绿**
- [ ] `GeometryEngine::build` **不**链接/构造 `FacePairList`；文档与注释含 rebuild 契约
- [ ] Mesh 几何代码 **不** include OpenLB LBM / `BlockLattice` 头文件
- [ ] **不**在 T07 写入 Bouzidi 距离/法向

---

## 阻塞关系

```
T01（OctreeForest）— 已完成
T02（FacePairList）— 已完成；T07 P1 仅冒烟
└── T07（本任务）· Mesh 几何链
    ├── T08 · vtk_writer（#11；可与 T07 W1 并行）
    ├── T09 · Bouzidi + BlockLattice 材料初始化（阻塞于 T07）
    └── T10+ · unit_converter + cavity/cylinder 集成（#12–#14；阻塞于 #0–#11 全绿）
```

---

## 构建与测试命令（给用户执行）

```bash
cd octlb
module load octlb
cmake -B build -DCMAKE_BUILD_TYPE=Debug \
      -DP4EST_ROOT=/home/dinglin/Public/p4est-2.8/opt
cmake --build build -j4

cd build && ctest --output-on-failure -R "test_stl_reader|test_geometry_adaptive_refine|test_geometry_engine"
```
