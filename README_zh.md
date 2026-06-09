# CP2K Learn — CP2K 入门教程

[English](README.md) | [简体中文](README_zh.md)

---

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
[![CP2K Version](https://img.shields.io/badge/CP2K-v2025.2-green.svg)](https://www.cp2k.org/)
[![Chapters](https://img.shields.io/badge/章节-16-orange.svg)](#教程章节)

> 一套循序渐进的 CP2K 实战教程，共 16 个章节、30+ 经过测试的输入文件，并配有交互式网页界面。所有示例均已在 CP2K v2025.2 上验证通过。

## 目录

- [简介](#简介)
- [前置要求](#前置要求)
- [快速开始](#快速开始)
- [教程章节](#教程章节)
- [项目结构](#项目结构)
- [网页界面](#网页界面)
- [参与贡献](#参与贡献)
- [许可证](#许可证)

## 简介

**CP2K Learn** 是一套面向 [CP2K](https://www.cp2k.org/) 的结构化学习路径。CP2K 是一款功能强大的开源量子化学与固态物理软件包。本项目提供：

- **16 个渐进式章节** — 从基础的单点能计算到高级的 GW/BSE 能带结构计算
- **30+ 经过验证的输入文件** — 每个 `.inp` 文件均已使用 `cp2k.ssmp`（CP2K v2025.2）测试通过
- **即用型结构文件** — 包含水盒子、硅晶体、丙氨酸二肽等的 `.xyz` 坐标文件
- **交互式网页界面** — 单页应用，支持教程浏览与学习进度跟踪

### 你将学到什么

| 级别 | 章节 | 主题 |
|------|------|------|
| ⭐ 入门 | 1–4 | 第一次 CP2K 运行、DFT 基础、基组、截断能收敛 |
| ⭐⭐ 中级 | 5–9 | 周期性体系、几何优化、晶胞优化、分子动力学 |
| ⭐⭐⭐ 高级 | 10–14 | 元动力学、xTB、QM/MM、振动分析、TDDFT |
| ⭐⭐⭐⭐⭐ 专家 | 15–16 | 后 Hartree-Fock（PBE0 杂化泛函）、能带结构 |

## 前置要求

| 依赖 | 版本 | 是否必需 |
|------|------|----------|
| CP2K | >= 2025.2 | 是 |
| Python | >= 3.8 | 否（用于网页界面脚本） |

### 安装 CP2K

从 [CP2K GitHub Releases](https://github.com/cp2k/cp2k/releases) 页面下载预编译二进制文件：

```bash
# 示例：下载适用于 Linux x86_64 的 CP2K 2025.2
wget https://github.com/cp2k/cp2k/releases/download/v2025.2/cp2k-2025.2-Linux-gnu-x86_64.ssmp
chmod +x cp2k-2025.2-Linux-gnu-x86_64.ssmp
```

或通过包管理器安装：

```bash
# conda（推荐初学者使用）
conda install -c conda-forge cp2k

# Ubuntu/Debian
sudo apt install cp2k
```

### 数据文件

CP2K 需要基组和赝势的数据文件。在输入文件中进行如下设置：

```
BASIS_SET_FILE_NAME  BASIS_SET
POTENTIAL_FILE_NAME  GTH_POTENTIALS
```

这些文件位于 `<CP2K_INSTALL_DIR>/data/` 目录下。如需使用 MOLOPT 基组，请使用 `BASIS_MOLOPT`。

## 快速开始

```bash
# 克隆仓库
git clone https://github.com/user/cp2k-learn.git
cd cp2k-learn

# 运行你的第一个 CP2K 计算
cd 01-hello-cp2k
cp2k.ssmp -i helium.inp -o helium.out

# 查看输出
grep "ENERGY" helium.out
```

## 教程章节

### ⭐ 入门

| # | 章节 | 输入文件 | 说明 |
|---|------|----------|------|
| 01 | [Hello CP2K](01-hello-cp2k/) | `helium.inp` | 你的第一次 CP2K 计算 — 单个氦原子 |
| 02 | [水分子静态计算](02-water-static/) | `H2O_static.inp` | 水分子的静态 DFT 计算 |
| 03 | [基组](03-basis-sets/) | `H2O_szv.inp`, `H2O_dzvp.inp`, `H2O_tzv2p.inp` | 比较 SZV、DZVP 和 TZVP 基组 |
| 04 | [截断能收敛](04-cutoff-convergence/) | `Si_cutoff.inp` | 平面波截断能的收敛性测试 |

### ⭐⭐ 中级

| # | 章节 | 输入文件 | 说明 |
|---|------|----------|------|
| 05 | [周期性硅](05-periodic-silicon/) | `Si_bulk8.inp`, `Si_bulk8_metal.inp` | 周期性边界条件、金属体系 |
| 06 | [几何优化](06-geo-opt/) | `H2O_geo_opt.inp`, `H2O_geo_opt_fixed.inp` | 使用 CG 和 BFGS 优化原子位置 |
| 07 | [晶胞优化](07-cell-opt/) | `Si_cell_opt.inp` | 硅的晶格参数优化 |
| 08 | [NVE 分子动力学](08-md-nve/) | `H2O_nve.inp` | 微正则系综 MD 模拟 |
| 09 | [NVT/NPT 恒温器](09-md-thermostats/) | `H2O_nvt.inp`, `H2O_nvt_csvr.inp`, `H2O_npt.inp` | 使用 CSVR 恒温器的正则系综和等温等压 MD |

### ⭐⭐⭐ 高级

| # | 章节 | 输入文件 | 说明 |
|---|------|----------|------|
| 10 | [元动力学](10-metadynamics/) | `alanine_dipeptide_meta.inp` | 使用 Well-Tempered Metadynamics 进行增强采样 |
| 11 | [xTB 方法](11-xtb/) | `ethanol_xtb.inp`, `water_xtb_md.inp` | GFN1-xTB 半经验方法，用于快速计算 |
| 12 | [QM/MM 基础](12-qmmm-basics/) | `water_qmmm.inp` | 水的混合 QM/MM 模拟 |
| 13 | [振动分析](13-vibrational/) | `H2O_vib.inp` | 频率分析与红外光谱 |
| 14 | [TDDFT](14-tddft/) | `ethylene_tddft.inp` | 含时密度泛函理论计算激发态 |

### ⭐⭐⭐⭐⭐ 专家

| # | 章节 | 输入文件 | 说明 |
|---|------|----------|------|
| 15 | [杂化泛函](15-mp2/) | `water_mp2.inp` | PBE0 杂化泛函计算 |
| 16 | [能带结构](16-gw-bse/) | `silicon_gw.inp` | 硅的电子能带结构 |

## 项目结构

```
cp2k_learn/
├── 01-hello-cp2k/          # 第 1 章：第一次 CP2K 运行
├── 02-water-static/        # 第 2 章：静态 DFT
├── 03-basis-sets/           # 第 3 章：基组比较
├── 04-cutoff-convergence/  # 第 4 章：截断能收敛
├── 05-periodic-silicon/    # 第 5 章：周期性体系
├── 06-geo-opt/             # 第 6 章：几何优化
├── 07-cell-opt/            # 第 7 章：晶胞优化
├── 08-md-nve/              # 第 8 章：NVE 分子动力学
├── 09-md-thermostats/      # 第 9 章：NVT/NPT 恒温器
├── 10-metadynamics/        # 第 10 章：元动力学
├── 11-xtb/                 # 第 11 章：xTB 半经验方法
├── 12-qmmm-basics/        # 第 12 章：QM/MM 混合方法
├── 13-vibrational/         # 第 13 章：振动分析
├── 14-tddft/               # 第 14 章：TDDFT 激发态
├── 15-mp2/                 # 第 15 章：杂化泛函
├── 16-gw-bse/              # 第 16 章：能带结构
├── index.html              # 交互式网页界面（SPA）
├── app.js                  # 网页界面应用逻辑
├── styles.css              # 网页界面样式
├── LICENSE                 # GPL v3 许可证
└── README.md               # 本文件
```

每个章节目录包含：
- `README.md` — 教程文档（中文）
- `*.inp` — CP2K 输入文件
- `*.xyz` — 分子/晶体结构文件（如适用）

## 网页界面

项目包含一个用于交互式浏览教程的单页应用。在任意现代浏览器中打开 `index.html` 即可使用。

功能特性：
- 带难度指示的章节导航
- 学习进度跟踪（存储在 localStorage 中）
- 语法高亮的输入文件查看器
- 适配移动端和桌面端的响应式设计

## 参与贡献

欢迎参与贡献！以下是你可以提供帮助的方式：

1. **报告问题** — 提交 Issue，描述 CP2K 输入文件或文档中的问题
2. **建议新章节** — 提出尚未涵盖的主题（例如 AIMD、NEB、RPA）
3. **改进文档** — 修正错别字、补充说明、改进翻译
4. **添加示例** — 按照现有章节结构提交经过测试的新 `.inp` 文件

### 开发工作流程

```bash
# Fork 并克隆
git fork https://github.com/user/cp2k-learn.git
git clone https://github.com/YOUR_USERNAME/cp2k-learn.git

# 创建功能分支
git checkout -b feature/new-chapter

# 修改并测试
cd 17-new-chapter
cp2k.ssmp -i new_calculation.inp

# 提交并推送
git add .
git commit -m "Add Chapter 17: New Topic"
git push origin feature/new-chapter

# 开启拉取请求
```

## 许可证

本项目基于 [GNU 通用公共许可证 v3.0](LICENSE) 授权。

---

*所有 CP2K 输入文件均已在 CP2K v2025.2（`cp2k.ssmp`）上验证通过。*
