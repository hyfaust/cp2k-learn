# 第五章：周期性体系 — 硅晶体

**难度：⭐⭐**

在前面的章节中，我们学习了气相分子的计算（如水分子）。然而，真实世界中的许多材料（金属、半导体、矿物等）都是以晶体形式存在的——原子在三维空间中周期性排列。本章将介绍如何在 CP2K 中处理**周期性体系**，以硅（Si）晶体作为示例。

---

## 目录

1. [周期性边界条件 (PBC)](#1-周期性边界条件-pbc)
2. [晶格与坐标](#2-晶格与坐标)
3. [K点采样](#3-k点采样)
4. [能带理论基础](#4-能带理论基础)
5. [展宽 (Smearing)](#5-展宽-smearing)
6. [Si 晶体结构](#6-si-晶体结构)
7. [输入文件详解](#7-输入文件详解)
8. [输出解读](#8-输出解读)
9. [练习](#9-练习)

---

## 1. 周期性边界条件 (PBC)

### 1.1 什么是 PBC？

在凝聚态物理和材料科学中，晶体的原子排列具有**平移对称性**——整个晶体可以看作是一个称为**原胞 (unit cell)** 的基本单元在三维空间中无限重复。

**周期性边界条件 (Periodic Boundary Conditions, PBC)** 的核心思想是：

- 我们只需计算**一个原胞**内的电子结构
- 该原胞在三个方向上无限重复
- 一个粒子从原胞一侧"出去"，会从对面"回来"
- 这消除了表面效应，模拟"无限大"体相材料

```
┌─────────┐┌─────────┐┌─────────┐
│  Unit   ││  Unit   ││  Unit   │
│  Cell   ││  Cell   ││  Cell   ││ → 无限重复...
│         ││         ││         │
└─────────┘└─────────┘└─────────┘
```

### 1.2 CP2K 中的 PBC

在 CP2K 中，通过 `&CELL` 部分定义周期性。每个方向可以独立设置为周期性或非周期性：

```
&CELL
  ABC  5.43 5.43 5.43   ! 立方晶胞，边长 5.43 Å
  PERIODIC XYZ           ! 三个方向均周期性
&END CELL
```

`PERIODIC` 关键字的选项：
| 值 | 含义 |
|------|------|
| `XYZ` | 三维周期性（体相材料）|
| `XY` | 二维周期性（表面/薄膜）|
| `X` | 一维周期性（纳米线）|
| `NONE` | 无周期性（气相分子）|

> **注意：** 对于气相分子（如前面章节的水分子），我们使用了较大的盒子但设置了 `PERIODIC NONE`。对于晶体，我们必须使用 `PERIODIC XYZ`。

### 1.3 为什么 PBC 对计算很重要？

PBC 从根本上改变了哈密顿量的形式：

- 气相分子：只需考虑分子内原子间的相互作用
- 周期性体系：需要考虑**Ewald 求和**来处理长程静电相互作用，以及对**布里渊区**的积分

CP2K 自动处理这些技术细节，但理解其背后的物理概念有助于正确设置参数。

---

## 2. 晶格与坐标

### 2.1 &CELL 部分

CP2K 中的 `&CELL` 部分定义晶胞的形状和大小。有两种指定方式：

#### 方式一：ABC 格式（边长 + 角度）

```
&CELL
  ABC  5.43  5.43  5.43    ! a, b, c (Å)
  ALPHA_BETA_GAMMA  90.0 90.0 90.0  ! α, β, γ (度)
&END CELL
```

- `ABC`：三条晶格矢量的长度（单位 Å）
- `ALPHA_BETA_GAMMA`：三条晶格矢量之间的夹角
  - `ALPHA`：b 与 c 之间的夹角
  - `BETA`：a 与 c 之间的夹角
  - `GAMMA`：a 与 b 之间的夹角
- 如果省略 `ALPHA_BETA_GAMMA`，默认为 90° 90° 90°（正交晶胞）

#### 方式二：A-B-C 格式（直接指定晶格矢量）

```
&CELL
  A  5.43  0.00  0.00
  B  0.00  5.43  0.00
  C  0.00  0.00  5.43
&END CELL
```

这种方式更加灵活，可以定义任意形状的晶胞（包括非正交的）。

### 2.2 晶格常数

**晶格常数 (lattice constant)** 是晶胞边长的物理量，是晶体最基本的参数之一。硅的晶格常数为：

$$a_{\text{Si}} = 5.4306975 \text{ Å}$$

在几何优化或晶胞优化中，我们可以让 CP2K 自动寻找最优晶格常数（见第七章）。

### 2.3 坐标格式

CP2K 支持多种坐标格式，在 `&COORD` 部分中指定：

#### Cartesian 坐标（笛卡尔坐标）

直接给出原子在三维空间中的 (x, y, z) 位置（单位 Å）：

```
&COORD
  Si  0.000  0.000  0.000
  Si  1.358  1.358  1.358
  ...
&END COORD
```

#### Scaled/Fractional 坐标（分数坐标）

给出原子在晶格矢量中的比例位置（0 到 1 之间）：

```
&COORD
  SCALED
  Si  0.000  0.000  0.000
  Si  0.250  0.250  0.250
  ...
&END COORD
```

使用 `SCALED` 关键字切换到分数坐标模式。

> **何时用哪种？** 分数坐标与晶胞大小无关——改变晶格常数时，原子的相对位置自动保持正确。在晶胞优化中推荐使用分数坐标。Cartesian 坐标更直观，适合对具体位置进行调整。

---

## 3. K点采样

### 3.1 布里渊区与 Bloch 定理

根据 **Bloch 定理**，周期性势场中的电子波函数可以写为：

$$\psi_{n\mathbf{k}}(\mathbf{r}) = e^{i\mathbf{k}\cdot\mathbf{r}} u_{n\mathbf{k}}(\mathbf{r})$$

其中 $\mathbf{k}$ 是**波矢**，$u_{n\mathbf{k}}$ 具有晶格周期性。

不同 $\mathbf{k}$ 对应不同的电子态。$\mathbf{k}$ 的取值构成**布里渊区 (Brillouin Zone, BZ)**——这是倒空间中的一个基本区域。

### 3.2 Monkhorst-Pack 网格

对布里渊区的积分需要在一系列 k 点上采样。**Monkhorst-Pack** 方法在布里渊区内均匀生成 k 点网格：

- N×N×N 网格：N 越大，采样越密，计算越精确但越耗时
- 1×1×1 网格 = 仅 Gamma 点（Γ-point only）

### 3.3 何时仅用 Gamma 点？

**仅 Gamma 点** 的情况：
- 原胞足够大（通常 > 10 Å 每个方向）
- 绝缘体/宽带隙半导体
- 分子在大盒子中

**需要 k 点采样** 的情况：
- 原胞较小（如 Si 的 8 原子胞 = 5.43 Å）
- 金属或窄带隙半导体
- 需要精确的能量和力

> **重要提示：** CP2K 中 k 点采样通过 `&KPOINTS` 部分实现，但目前**仅适用于部分方法**（如 GPW 方案配合某些设置）。对于本教程中的常规计算，我们主要使用 Gamma-point。如需精确的周期性计算，应使用较大的超胞或确认 CP2K 版本支持所需功能。

### 3.4 CP2K 中的 &KPOINTS

```
&KPOINTS
  SCHEME MONKHORST-PACK  2 2 2   ! 2×2×2 Monkhorst-Pack 网格
  SYMMETRY  .TRUE.                ! 利用对称性减少 k 点数
  EPS_GEO  1.0E-6                 ! k 点几何精度
&END KPOINTS
```

> **注意：** 在本章的示例中，我们使用 Gamma-point only 计算，这对教学目的足够。对于研究级 Si 计算，应进行 k 点收敛测试。

---

## 4. 能带理论基础

### 4.1 能带结构

在孤立原子中，电子占据离散的能级。当原子周期性排列形成晶体时，这些离散能级展宽为**能带 (energy bands)**：

```
孤立原子              晶体
  ── E₃           ═══════════  能带3
                   
  ── E₂           ═══════════  能带2
                   
  ── E₁           ═══════════  能带1
```

- 每条能带对应一组电子态，能量随 k 点变化
- 能带之间可能存在**带隙 (band gap)**——不允许电子占据的能量范围

### 4.2 金属、绝缘体与半导体

根据电子填充和带隙大小，材料分为：

| 类型 | 带隙 | 特征 |
|------|------|------|
| **金属** | 0 eV | 最高占据带未填满，或占据带与空带重叠 |
| **半导体** | 0.1 ~ 4 eV | 有较小带隙，热激发可产生载流子 |
| **绝缘体** | > 4 eV | 大带隙，电子难以跃迁 |

硅是**半导体**，实验带隙约 1.12 eV（间接带隙）。注意 DFT 通常会低估带隙（LDA/GGA 给出约 0.5 eV）。

### 4.3 费米能级

**费米能级 (Fermi level)** 是电子化学势，在 T=0 时是最高占据态的能量：

- 金属：费米能级位于能带内部
- 半导体/绝缘体：费米能级位于带隙中

---

## 5. 展宽 (Smearing)

### 5.1 为什么需要展宽？

对于**金属**或**窄带隙半导体**，HOMO 和 LUMO 能级非常接近（甚至简并）。这会导致 SCF 收敛困难——电子在占据和非占据态之间"振荡"。

**展宽 (Smearing)** 的解决方案：不将电子严格占据到费米能级以下，而是用一个连续的**费米函数**来平滑占据数：

$$f_i = \frac{1}{1 + e^{(\varepsilon_i - \mu) / k_B T_{\text{elec}}}}$$

其中 $T_{\text{elec}}$ 是**电子温度**，$\mu$ 是化学势（费米能级）。

### 5.2 Fermi-Dirac 展宽

CP2K 默认且最常用的展宽方案是 **Fermi-Dirac (FD)** 展宽：

```
&SMEAR
  METHOD  FERMI_DIRAC
  ELECTRONIC_TEMPERATURE  300   ! 电子温度 (K)
&END SMEAR
```

- 电子温度 $T_{\text{elec}}$ 控制展宽宽度
- $T_{\text{elec}} = 0$ K → 阶梯函数（无展宽）
- $T_{\text{elec}} = 300$ K → $k_BT \approx 0.026$ eV，较温和的展宽
- $T_{\text{elec}} = 1000$ K → 更大的展宽，更容易收敛但物理意义改变

### 5.3 ADDED_MOS

使用展宽时，必须提供**额外的分子轨道 (MOS)** 来允许电子占据费米能级以上：

```
&DFT
  ADDED_MOS  10     ! 添加 10 个额外的空轨道
  &SMEAR
    METHOD  FERMI_DIRAC
    ELECTRONIC_TEMPERATURE  300
  &END SMEAR
  &MIXING
    METHOD  BROYDEN_MIXING
    ALPHA  0.4
    NBUFFER  7
  &END MIXING
&END DFT
```

- `ADDED_MOS` 指定额外轨道数
- 对于金属，通常需要较多的额外轨道（如 10~50）
- 对于半导体，少量即可（如 5~10）

### 5.4 展宽对能量的影响

使用展宽后，CP2K 输出中的总能量包含**熵贡献**：

$$E_{\text{total}} = E_{\text{band-structure}} + E_{\text{Hartree}} + E_{\text{xc}} + ... - T_{\text{elec}} \cdot S$$

CP2K 同时输出 `ENERGY| Total` 和 `ENTROPIC| Entropic contribution`。对于物理量的比较，应使用**无熵修正的能量**（即 `Total` 减去 `Entropic contribution`），或者在 $T_{\text{elec}} \to 0$ 极限下取值。

---

## 6. Si 晶体结构

### 6.1 金刚石立方结构

硅具有**金刚石立方 (diamond cubic)** 结构，是面心立方 (FCC) 的变体：

- 基本框架：FCC 格子
- 在每个 FCC 格点上有两个原子（间距为体对角线的 1/4）
- 因此**常规 FCC 晶胞含有 8 个 Si 原子**

### 6.2 8 原子 FCC 晶胞的坐标

晶格常数 $a = 5.4306975$ Å，8 个 Si 原子的 Cartesian 坐标为：

| 原子 | x (Å) | y (Å) | z (Å) | 分数坐标 |
|------|--------|--------|--------|----------|
| Si 1 | 0.0000 | 0.0000 | 0.0000 | (0, 0, 0) |
| Si 2 | 1.3577 | 1.3577 | 1.3577 | (1/4, 1/4, 1/4) |
| Si 3 | 2.7153 | 2.7153 | 0.0000 | (1/2, 1/2, 0) |
| Si 4 | 4.0730 | 4.0730 | 1.3577 | (3/4, 3/4, 1/4) |
| Si 5 | 2.7153 | 0.0000 | 2.7153 | (1/2, 0, 1/2) |
| Si 6 | 4.0730 | 1.3577 | 4.0730 | (3/4, 1/4, 3/4) |
| Si 7 | 0.0000 | 2.7153 | 2.7153 | (0, 1/2, 1/2) |
| Si 8 | 1.3577 | 4.0730 | 4.0730 | (1/4, 3/4, 3/4) |

其中 $1.3577 \approx 5.4306975 / 4$。

### 6.3 键长与配位

- 每个 Si 原子与 4 个最近邻 Si 原子形成共价键
- Si-Si 键长 = $\frac{\sqrt{3}}{4} a \approx 2.352$ Å
- 配位数 = 4（sp³ 杂化）
- 每个原子贡献 4 个价电子，形成完全填满的价带

---

## 7. 输入文件详解

### 7.1 Si_bulk8.inp — 能量/力计算

这是最基本的 Si 晶体计算，计算总能量和原子力。

```bash
# 运行命令
cp2k.psmp -i Si_bulk8.inp -o Si_bulk8.out
```

关键部分逐行解释：

#### GLOBAL 部分

```
&GLOBAL
  PROJECT  Si_bulk8        ! 项目名称，输出文件以此为前缀
  RUN_TYPE ENERGY_FORCE     ! 计算类型：单点能量和力
  PRINT_LEVEL LOW           ! 输出详细程度
&END GLOBAL
```

- `RUN_TYPE ENERGY_FORCE`：这是最基本的计算类型，只计算一次电子结构，输出总能量和原子力
- `PROJECT Si_bulk8`：所有输出文件将命名为 `Si_bulk8-*`

#### FORCE_EVAL 部分

```
&FORCE_EVAL
  METHOD  QS                ! 使用 Quickstep (DFT) 方法
  &DFT
    BASIS_SET_FILE_NAME  BASIS_SET       ! 基组文件
    POTENTIAL_FILE_NAME  GTH_POTENTIALS  ! 赝势文件
    
    &MGRID
      CUTOFF  300           ! 平面波截断能 (Ry)
      REL_CUTOFF  60        ! 相对截断能 (Ry)
    &END MGRID
    
    &QS
      EPS_DEFAULT  1.0E-10  ! 各种默认精度
    &END QS
    
    &MIXING
      METHOD  BROYDEN_MIXING  ! Broyden 密度混合
      ALPHA  0.4              ! 混合参数
      NBUFFER  7              ! 存储的历史步数
    &END MIXING
    
    &SCF
      SCF_GUESS  ATOMIC      ! 初始猜测：原子密度叠加
      EPS_SCF  1.0E-6        ! SCF 收敛判据
      MAX_SCF  50            ! 最大 SCF 迭代次数
    &END SCF
    
    &XC
      &XC_FUNCTIONAL PADE    ! LDA (PADE 近似)
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  
  &SUBSYS
    &CELL
      ABC  5.4306975  5.4306975  5.4306975  ! 立方晶胞 (Å)
      PERIODIC  XYZ                          ! 三维周期性
    &END CELL
    
    &TOPOLOGY
      COORD_FILE_FORMAT  XYZ   ! 坐标从下面的 &COORD 读取
    &END TOPOLOGY
    
    &KIND Si
      BASIS_SET  DZVP-GTH-PADE   ! 双ζ价极化基组
      POTENTIAL  GTH-PADE-q4     ! GTH 赝势，Si 有 4 个价电子
    &END KIND
    
    &COORD
      Si   0.000000   0.000000   0.000000
      Si   1.357674   1.357674   1.357674
      Si   2.715349   2.715349   0.000000
      Si   4.073023   4.073023   1.357674
      Si   2.715349   0.000000   2.715349
      Si   4.073023   1.357674   4.073023
      Si   0.000000   2.715349   2.715349
      Si   1.357674   4.073023   4.073023
    &END COORD
  &END SUBSYS
&END FORCE_EVAL
```

**关键参数说明：**

| 参数 | 值 | 说明 |
|------|-----|------|
| `CUTOFF` | 300 Ry | 平面波截断能，控制实空间网格精度 |
| `REL_CUTOFF` | 60 Ry | 相对截断能，控制最精细网格 |
| `BROYDEN_MIXING` | — | 适合固体的密度混合方法 |
| `SCF_GUESS` | ATOMIC | 用原子密度叠加作为 SCF 初始猜测 |
| `EPS_SCF` | 1E-6 | SCF 收敛判据（能量变化） |
| `DZVP-GTH-PADE` | — | 双ζ价极化基组 |
| `GTH-PADE-q4` | — | GTH 赝势，Si 有 4 个价电子 (3s²3p²) |

### 7.2 Si_bulk8_metal.inp — 带展宽的计算

此文件在基本计算的基础上增加了 **Fermi-Dirac 展宽**，适合处理金属性或窄带隙体系。

**与 Si_bulk8.inp 的主要区别：**

```
&GLOBAL
  PROJECT  Si_bulk8_metal      ! 不同的项目名
  RUN_TYPE ENERGY_FORCE
  PRINT_LEVEL LOW
&END GLOBAL

&FORCE_EVAL
  ...
  &DFT
    ...
    ADDED_MOS  10              ! ★ 新增：10 个额外分子轨道
    
    &SMEAR                     ! ★ 新增：展宽部分
      METHOD  FERMI_DIRAC      ! Fermi-Dirac 展宽
      ELECTRONIC_TEMPERATURE  300  ! 电子温度 300 K
    &END SMEAR
    
    &MIXING
      METHOD  BROYDEN_MIXING
      ALPHA  0.4
      NBUFFER  7
    &END MIXING
    
    &PRINT                     ! ★ 新增：输出部分
      &MO                      ! 输出分子轨道信息
        OCCUPATION_NUMBERS  .TRUE.  ! 打印占据数
      &END MO
    &END PRINT
    ...
  &END DFT
  ...
&END FORCE_EVAL
```

**新增参数说明：**

| 参数 | 值 | 说明 |
|------|-----|------|
| `ADDED_MOS` | 10 | 额外分子轨道数，展宽需要这些轨道容纳被"激发"的电子 |
| `METHOD FERMI_DIRAC` | — | Fermi-Dirac 展宽方案 |
| `ELECTRONIC_TEMPERATURE` | 300 K | 电子温度，控制展宽宽度 |
| `OCCUPATION_NUMBERS` | .TRUE. | 打印每个轨道的占据数，可以看到非整数占据 |

---

## 8. 输出解读

### 8.1 Si_bulk8.out 的关键输出

#### 总能量

```
ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):         -31.36850849423
```

这是最终的总能量（Hartree 单位）。对于 `ENERGY_FORCE` 计算，这就是单点能量。

#### 原子力

```
ATOMIC FORCES in [a.u.]

#   Atom   Kind   Element          X              Y              Z
      1      1      Si          0.00000000     0.00000000     0.00000000
      2      1      Si          0.00000000     0.00000000     0.00000000
      ...
```

对于完美的晶体结构（且没有几何优化），力应非常接近零。如果力不为零，说明原子位置可能需要优化。

#### SCF 迭代信息

```
  Step     Update method      Time    Convergence         Total energy   Change
  ------------------------------------------------------------------------------
     1 OT DIIS     0.15E+00    3.2     0.00014587      -31.36848954   -3.14E+01
     2 OT DIIS     0.10E+00    2.8     0.00000892      -31.36850849   -1.90E-05
     3 OT DIIS     0.85E-01    2.6     0.00000012      -31.36850849   -2.15E-08
```

观察 SCF 收敛过程——能量和密度残差应单调递减。

### 8.2 Si_bulk8_metal.out 的关键输出

#### 总能量（含熵贡献）

```
ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):         -31.36882546891

ENTROPIC| Entropic contribution (a.u.) :                  -0.00012345678
```

**解读：**
- `Total`：包含熵贡献的总能量
- `Entropic contribution`：$-T_{\text{elec}} \cdot S$ 项
- 对于能量比较，如果两个计算使用相同的电子温度，直接比较 `Total` 即可
- 如果要与零温度计算比较，应使用 `Total - Entropic`

#### 占据数

使用展宽后，某些轨道的占据数不再是严格的 0 或 1：

```
MO| Molecular Orbital     5:                       Occupancy: 1.99998
MO| Molecular Orbital     6:                       Occupancy: 1.99996
MO| Molecular Orbital     7:                       Occupancy: 0.50012
MO| Molecular Orbital     8:                       Occupancy: 0.00003
```

注意轨道 7 的占据数约为 0.5——它位于费米能级附近，被部分占据。这正是展宽的效果。

---

## 9. 练习

### 练习 1：截断能收敛测试

对 `Si_bulk8.inp`，分别设置 `CUTOFF` 为 200, 300, 400, 500, 600 Ry，运行计算并绘制总能量 vs 截断能曲线。确定收敛所需的最小截断能。

### 练习 2：比较有无展宽的计算

运行 `Si_bulk8.inp` 和 `Si_bulk8_metal.inp`，比较：
- 总能量差异
- SCF 迭代次数
- `Si_bulk8_metal.out` 中的占据数

### 练习 3：改变电子温度

修改 `Si_bulk8_metal.inp` 中的 `ELECTRONIC_TEMPERATURE`，分别设为 100, 300, 1000, 3000 K。观察：
- 熵贡献如何变化
- 占据数如何变化
- SCF 收敛行为有何不同

### 练习 4：超胞计算

将 8 原子晶胞扩展为 64 原子（2×2×2 超胞），仅用 Gamma 点计算。与 8 原子结果比较（注意：总能量不是可直接比较的量，应比较**每原子能量**）。

---

## 小结

本章学习了 CP2K 中处理周期性体系的关键概念：

| 概念 | 关键词/参数 |
|------|------------|
| 周期性边界条件 | `PERIODIC XYZ` |
| 晶胞定义 | `&CELL` ABC / A-B-C |
| 坐标格式 | `&COORD` Cartesian / SCALED |
| K 点采样 | `&KPOINTS` MONKHORST-PACK |
| 展宽 | `&SMEAR` FERMI_DIRAC, `ADDED_MOS` |
| 赝势与基组 | `GTH-PADE-q4`, `DZVP-GTH-PADE` |

下一章将介绍**几何优化**——如何找到原子在势能面上的最低能量构型。

---

**参考资源：**
- CP2K 手册 - &CELL：https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/SUBSYS/CELL.html
- CP2K 手册 - &SMEAR：https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/SCF/SMEAR.html
- CP2K 手册 - &KPOINTS：https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/KPOINTS.html
