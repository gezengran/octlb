# OctLB 文档

`doc/` 采用**产品 / 工程分离方案**，顶层仅三个目录（**英文通用缩写**）：**做什么**、**细化怎么做**、**具体实现细节**。

```
doc/
├── README.md
├── prd/               ← 产品需求（Product Requirements Document）
├── tasks/             ← 拆分的任务
└── dev/               ← 架构、接口、ADR、构建说明
```

---

## prd/

**回答「做什么、做到什么算完成」。** 读者：产品负责人、仿真工程师、评审人。

| 放这里 | 不放这里 |
|--------|----------|
| PRD、问题陈述、用户故事 | 类图、头文件清单、CMake 细节 |
| 验收标准（如 Cd 相对误差 &lt; 1%） | 具体实现步骤、Issue 拆解 |
| 范围与非目标（out of scope） | 已拍板的技术方案长文（见 `dev/`） |

**现有文件**

| 文件 | 说明 |
|------|------|
| [prd/octlb-framework.md](prd/octlb-framework.md) | OctLB 框架总 PRD（v0.2；含 T09/T10 进度） |

**tasks/ 任务文档**

| 文件 | 说明 |
|------|------|
| [tasks/T10-cavity3d.md](tasks/T10-cavity3d.md) | T10：unit_converter + cavity3d（W1/W2/W3，已完成） |
| [tasks/T11-cylinder3d.md](tasks/T11-cylinder3d.md) | T11：cylinder3d（P0 修 ①–⑤ + W1–W4，oracle 三阶 sanity→量级→Cd<1%，进行中） |

**后续可增示例**

- `prd/acceptance.md` — 跨算例的统一验收表
- `prd/glossary.md` — 术语（octant、FacePair、Lagrava 等）

---

## tasks/

**回答「按什么顺序做、做到哪一步算阶段性完成」。** 读者：实现者、自己排期时用。

| 放这里 | 不放这里 |
|--------|----------|
| 里程碑 / 阶段划分（M0 拓扑、M1 ghost…） | PRD 级「为什么要做」 |
| 可领取的工作项说明（对应 Issue 或 checklist） | 架构决策理由（见 `dev/adr/`） |
| 验证**执行**清单（跑哪个 example、对比哪条曲线） | 模块接口与数据结构设计（见 `dev/design/`） |

**与 `prd/` 的边界**：`prd/` 写「cylinder3d 的 Cd 误差 &lt; 1%」；`tasks/` 写「第 3 周跑 cylinder3d，对比 OpenLB 参考值，记录于 xxx」。

**后续可增示例**

- `tasks/milestone-m0-mesh.md`
- `tasks/verify-cylinder3d.md`

---

## dev/

**回答「系统长什么样、为何这样设计、本地如何构建」。** 读者：框架开发者、维护者。

| 放这里 | 不放这里 |
|--------|----------|
| 架构总览、模块边界、数据流 | 用户故事、商业优先级 |
| 接口约定（`FacePairList`、索引 `(octant_id,i,j,k)`） | 迭代排期与工单状态 |
| ADR（Architecture Decision Records） | 验收指标定义 |
| 构建与环境说明（`module load octlb` 等） | |

建议在 `dev/` **内部**用子目录组织（仍只占一个顶层文件夹）：

```
dev/
├── design/       # 设计说明（可从 PRD「实现决策」拆出）
├── adr/          # 一条决策一个文件，如 0001-ghost-face-exchange.md
└── build.md      # 编译、依赖、并行度约定（-j4）
```

**后续可增示例**

- `dev/design/mesh.md`、`dev/design/solver-field.md`、`dev/design/solver-lbm.md`
- `dev/adr/0001-face-pair-via-p8est-iterate.md`

---

## 三目录对照（缩写含义）

| 目录 | 含义 | 典型内容 |
|------|------|----------|
| `prd/` | Product Requirements Document | 用户故事、验收标准、范围 |
| `tasks/` | 执行与排期 | 里程碑、验证 checklist、Issue 对应说明 |
| `dev/` | Development / engineering | `design/`、`adr/`、`build.md` |

---

## 维护约定

1. **新增文档前先定类**：`prd/` / `tasks/` / `dev/` 三者选一；拿不准时先问再落盘。
2. **Agent 默认以** `prd/octlb-framework.md` **为需求基线**；设计与构建以 `dev/` 为准。
