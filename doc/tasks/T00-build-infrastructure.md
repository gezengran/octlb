# T00 · Build Infrastructure

> 类型：AFK  
> 阻塞于：无，可立即开始  
> 对应 PRD：`octlb/doc/prd/octlb-framework.md`（测试决策节，测试顺序 #0 前置）

---

## 要做什么

搭建 OctLB 项目的 CMake 骨架，使一条空测试可以通过编译、链接并被 `ctest` 执行。
本任务**不引入任何生产源文件**，只验证工具链通路。

---

## 交付物

```
octlb/
├── CMakeLists.txt                        # 根：项目声明、C++17、find_package(GTest)
├── tests/
│   ├── CMakeLists.txt                    # enable_testing() + add_subdirectory
│   └── unit/
│       ├── CMakeLists.txt                # add_subdirectory(mesh)
│       └── mesh/
│           ├── CMakeLists.txt            # 空测试可执行文件
│           └── test_placeholder.cc      # 含一条 EXPECT_TRUE(true)
```

`src/` 目录保持空，等 T01 写第一个生产组件时自然生长。

---

## 关键设计决策（已对齐）

| 项目 | 决策 |
|---|---|
| GTest 引入方式 | `find_package(GTest REQUIRED)`（系统已安装） |
| CMakeLists 层级 | 分层：根 → tests/ → unit/ → unit/mesh/ |
| C++ 标准 | C++17（与 OpenLB 1.9.0 保持一致） |
| src/ 占位文件 | 无 |
| MPI | 本任务不引入，T01 引入 p4est 时再添加 |

---

## 验收标准

- [ ] `cmake -B build && cmake --build build -j4` 成功，无 warning
- [ ] `cd build && ctest --output-on-failure` 输出 `1/1 Test #1 ... Passed`
- [ ] `src/` 目录下无任何文件
- [ ] CMakeLists.txt 中 C++ 标准锁定为 17，通过 `target_compile_features` 或 `CMAKE_CXX_STANDARD` 设置

---

## 阻塞关系

```
T00（本任务）
└── T01 · OctreeForest（test_octree_forest，依赖本任务工具链）
```

---

## 构建与测试命令（给用户执行）

```bash
cd octlb
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j4
cd build && ctest --output-on-failure
```

若 GTest 安装在非标准路径，可追加：

```bash
cmake -B build -DGTest_DIR=/path/to/gtest/lib/cmake/GTest
```
