# 第六章：几何优化

**难度：⭐⭐**

前几章中，我们计算了给定结构的单点能量和力。但在实际应用中，我们通常不知道原子的"最优"位置——我们需要通过**几何优化 (Geometry Optimization)** 来找到体系在势能面上的最低能量构型。

---

## 目录

1. [几何优化原理](#1-几何优化原理)
2. [RUN_TYPE GEO_OPT](#2-run_type-geo_opt)
3. [优化算法](#3-优化算法)
4. [收敛判据](#4-收敛判据)
5. [约束优化](#5-约束优化)
6. [输出文件](#6-输出文件)
7. [SCF 混合策略](#7-scf-混合策略)
8. [输入文件详解](#8-输入文件详解)
9. [练习](#9-练习)

---

## 1. 几何优化原理

### 1.1 势能面 (PES)

**势能面 (Potential Energy Surface, PES)** 是体系能量作为原子坐标的函数：

$$E = E(\mathbf{R}_1, \mathbf{R}_2, ..., \mathbf{R}_N)$$

其中 $\mathbf{R}_i$ 是第 $i$ 个原子的位置矢量。对于 $N$ 个原子的体系，PES 是 $3N$ 维的超曲面。

```
能量 E
  │    ╱╲
  │   ╱  ╲        ╱╲
  │  ╱    ╲      ╱  ╲
  │ ╱      ╲    ╱    ╲
  │╱        ╲──╱      ╲──→ 原子坐标 R
  │     *              *
  │   最小值         鞍点
```

### 1.2 Born-Oppenheimer 近似

在 Born-Oppenheimer 近似下，原子核在由电子提供的势能面上运动：

1. 给定原子核位置 → 求解电子薛定谔方程 → 得到电子基态能量 $E_{\text{elec}}$
2. $E_{\text{elec}} + V_{\text{nuc-nuc}} = E_{\text{total}}(\mathbf{R})$ → 这就是 PES
3. 原子核在 PES 上运动

几何优化就是在这个 PES 上寻找**局部极小值 (local minimum)** 的过程。

### 1.3 梯度与力

PES 的**梯度 (gradient)** 等于原子力的负值：

$$\mathbf{F}_i = -\frac{\partial E}{\partial \mathbf{R}_i}$$

在能量极小值处，所有原子力为零（或足够小）。几何优化的目标就是找到力为零的构型。

### 1.4 优化的基本流程

```
初始结构 → 计算能量和力 → 是否收敛？→ 是 → 输出最优结构
                ↑                              ↓
                │                              否
                └── 根据力和算法移动原子 ←──────┘
```

每一步的循环称为一个**优化步 (optimization step)**，每步中包含一次完整的电子结构计算（SCF 循环）。

---

## 2. RUN_TYPE GEO_OPT

在 CP2K 中执行几何优化，需要将 `RUN_TYPE` 设为 `GEO_OPT`：

```
&GLOBAL
  PROJECT  H2O_geo_opt
  RUN_TYPE GEO_OPT       ! ★ 几何优化
&END GLOBAL
```

CP2K 会自动：
1. 执行 SCF 计算得到能量和力
2. 调用优化算法更新原子位置
3. 重复直到满足收敛判据或达到最大步数

优化参数在 `&GEO_OPT` 部分设置（位于 `&MOTION` 下）：

```
&MOTION
  &GEO_OPT
    OPTIMIZER  CG          ! 优化算法
    MAX_ITER  200           ! 最大优化步数
    MAX_DR    1.0E-3       ! 最大位移收敛判据 (Bohr)
    RMS_DR    1.0E-3       ! 均方位移收敛判据 (Bohr)
    MAX_FORCE  1.0E-3      ! 最大力收敛判据 (Bohr⁻¹·Ha)
    RMS_FORCE  1.0E-3      ! 均方力收敛判据 (Bohr⁻¹·Ha)
  &END GEO_OPT
&END MOTION
```

---

## 3. 优化算法

### 3.1 最速下降法 (Steepest Descent, SD)

最简单的方法：沿力（梯度负方向）移动原子。

$$\mathbf{R}_{n+1} = \mathbf{R}_n + \alpha \cdot \mathbf{F}_n$$

其中 $\alpha$ 是步长。

**优点：**
- 实现简单，总是下降方向
- 适合远离极小值的初始阶段

**缺点：**
- 收敛速度慢（线性收敛）
- 会在峡谷形 PES 上"之字形"振荡

```
SD 收敛路径示意：
    ╲  │  ╱
     ╲ │ ╱
      ╲│╱
       │  ← 之字形振荡
      ╱│╲
     ╱ │ ╲
```

### 3.2 共轭梯度法 (Conjugate Gradients, CG)

改进最速下降法：当前搜索方向与前一步方向**共轭**，避免重复搜索同一方向。

$$\mathbf{d}_{n+1} = \mathbf{F}_{n+1} + \beta_n \cdot \mathbf{d}_n$$

其中 $\beta_n$ 有不同的选择公式（Fletcher-Reeves, Polak-Ribiere 等）。

**优点：**
- 比 SD 快得多（超线性收敛）
- 不需要存储 Hessian 矩阵
- 内存需求低

**缺点：**
- 需要精确的线搜索
- 可能需要周期性重启

### 3.3 BFGS (Broyden-Fletcher-Goldfarb-Shanno)

**准牛顿方法**：逐步构建 Hessian 矩阵（二阶导数矩阵）的近似：

$$\mathbf{R}_{n+1} = \mathbf{R}_n - \mathbf{H}_n^{-1} \cdot \mathbf{g}_n$$

其中 $\mathbf{g}_n$ 是梯度，$\mathbf{H}_n$ 是近似 Hessian。

**优点：**
- 二次收敛速度（非常快）
- 适合接近极小值的精细优化

**缺点：**
- 内存需求较高（需存储 Hessian 矩阵）
- 对大体系可能不适用

### 3.4 L-BFGS (Limited-memory BFGS)

BFGS 的内存优化版本：只存储最近几步的梯度和位移信息，而不是完整的 Hessian 矩阵。

**优点：**
- 保持 BFGS 的快速收敛
- 内存需求与体系大小成线性关系
- **CP2K 中的默认推荐算法**

### 3.5 算法选择指南

| 场景 | 推荐算法 |
|------|---------|
| 初始结构远离极小值 | SD 或 CG |
| 一般几何优化 | CG 或 BFGS |
| 大体系（>几百原子）| L-BFGS |
| 快速收敛 | BFGS |

在 CP2K 中通过 `OPTIMIZER` 关键字指定：

```
&GEO_OPT
  OPTIMIZER  CG     ! 或 SD, BFGS, L-BFGS
&END GEO_OPT
```

---

## 4. 收敛判据

### 4.1 四个收敛标准

CP2K 的几何优化使用**四个收敛判据**，必须**全部满足**才认为优化收敛：

| 判据 | 关键字 | 含义 | 默认值 |
|------|--------|------|--------|
| 最大位移 | `MAX_DR` | 任何原子的最大位移 | 3.0E-3 Bohr |
| 均方位移 | `RMS_DR` | 所有原子位移的均方根 | 1.5E-3 Bohr |
| 最大力 | `MAX_FORCE` | 任何原子上力分量的最大值 | 4.5E-4 Ha/Bohr |
| 均方力 | `RMS_FORCE` | 所有原子力分量的均方根 | 3.0E-4 Ha/Bohr |

### 4.2 如何设置收敛判据

```
&GEO_OPT
  OPTIMIZER  CG
  MAX_ITER  200
  
  ! 收敛判据（全部需满足）
  MAX_DR    1.0E-3     ! 最大位移 (Bohr)
  RMS_DR    1.0E-3     ! RMS 位移 (Bohr)
  MAX_FORCE  1.0E-3     ! 最大力 (Ha/Bohr)
  RMS_FORCE  1.0E-3     ! RMS 力 (Ha/Bohr)
&END GEO_OPT
```

> **实用建议：**
> - 初步优化（探索构型）：使用较松的标准，如 1E-3
> - 精确优化（用于频率计算等）：使用较严的标准，如 1E-4 或 1E-5
> - 力的判据通常比位移判据更重要

### 4.3 单位说明

CP2K 内部使用原子单位 (a.u.)：

- 长度：1 Bohr = 0.529177 Å
- 力：1 Ha/Bohr = 51.422 eV/Å = 51.422 nN
- 能量：1 Ha = 27.2114 eV

---

## 5. 约束优化

### 5.1 什么是约束优化？

有时我们希望固定某些原子的位置，只优化其余原子。常见场景：

- 研究表面吸附：固定表面原子，优化吸附分子
- 模拟缺陷：固定远离缺陷的原子
- 计算特定坐标的势能曲线
- 降低计算成本

### 5.2 &FIXED_ATOMS 部分

在 `&MOTION` 中添加 `&FIXED_ATOMS` 部分：

```
&MOTION
  &GEO_OPT
    OPTIMIZER  CG
    MAX_ITER  200
    MAX_FORCE  1.0E-3
    RMS_FORCE  1.0E-3
  &END GEO_OPT
  
  &CONSTRAINT
    &FIXED_ATOMS
      LIST  1            ! 固定第 1 号原子
      COMPONENTS_TO_FIX  XYZ  ! 固定所有三个分量
    &END FIXED_ATOMS
  &END CONSTRAINT
&END MOTION
```

### 5.3 FIXED_ATOMS 关键字详解

| 关键字 | 说明 |
|--------|------|
| `LIST` | 要固定的原子序号（从 1 开始），空格分隔多个原子 |
| `COMPONENTS_TO_FIX` | 要固定的分量：`XYZ`（全部）、`XY`、`XZ`、`YZ`、`X`、`Y`、`Z` |

**示例：**

```
! 固定原子 1 和 2
&FIXED_ATOMS
  LIST  1  2
  COMPONENTS_TO_FIX  XYZ
&END FIXED_ATOMS

! 只固定原子 1 的 z 坐标（x 和 y 自由移动）
&FIXED_ATOMS
  LIST  1
  COMPONENTS_TO_FIX  Z
&END FIXED_ATOMS
```

---

## 6. 输出文件

### 6.1 主输出文件

`H2O_geo_opt.out`：包含每次优化步的详细信息：
- SCF 收敛过程
- 总能量
- 原子力
- 优化算法的位移向量
- 收敛判据的当前值

### 6.2 轨迹文件

`H2O_geo_opt-pos-1.xyz`：记录每一步优化的原子坐标。

可以用 VMD、OVITO 或其他可视化软件打开，观察优化过程：

```bash
vmd H2O_geo_opt-pos-1.xyz
```

### 6.3 重启文件

`H2O_geo_opt-1.restart`：包含最后一优化步的所有信息，可用于：
- 继续未收敛的优化
- 作为后续计算（如频率分析）的输入结构

使用重启文件：

```
&EXT_RESTART
  RESTART_FILE_NAME  H2O_geo_opt-1.restart
&END EXT_RESTART
```

---

## 7. SCF 混合策略

### 7.1 DIIS/Pulay 混合

在几何优化过程中，原子位置每步都在变化。高效的 SCF 混合策略可以利用前几步的密度信息加速收敛。

```
&MIXING
  METHOD  BROYDEN_MIXING    ! 混合方法
  ALPHA  0.4                ! 新旧密度的混合比例
  NBUFFER  7                ! 存储的历史密度数
&END MIXING
```

### 7.2 参数说明

| 参数 | 说明 | 建议值 |
|------|------|--------|
| `ALPHA` | 线性混合参数：$\rho_{\text{new}} = \alpha \cdot \rho_{\text{output}} + (1-\alpha) \cdot \rho_{\text{input}}$ | 0.2 ~ 0.5 |
| `NBUFFER` | 存储的历史密度/DIIS 向量数 | 4 ~ 8 |

**调节建议：**
- SCF 不收敛 → 减小 `ALPHA`（如 0.1）
- SCF 收敛慢 → 增大 `ALPHA` 或增大 `NBUFFER`
- Broyden 混合通常比简单的 Pulay 混合更稳健

### 7.3 混合方法选择

| 方法 | 关键字 | 适用场景 |
|------|--------|---------|
| Broyden | `BROYDEN_MIXING` | 通用，推荐默认 |
| Pulay | `PULAY_MIXING` | SCF 收敛良好的体系 |
| Kerker | `KERKER_MIXING` | 金属体系 |
| 直接混合 | `DIRECT_P_MIXING` | 简单体系，调试用 |

---

## 8. 输入文件详解

### 8.1 H2O_geo_opt.inp — 水分子几何优化

这是标准的水分子几何优化：从初始猜测结构出发，优化到最近的能量极小值。

```bash
# 运行命令
cp2k.psmp -i H2O_geo_opt.inp -o H2O_geo_opt.out
```

```
&GLOBAL
  PROJECT  H2O_geo_opt        ! 项目名称
  RUN_TYPE GEO_OPT             ! ★ 几何优化
  PRINT_LEVEL LOW
&END GLOBAL

&FORCE_EVAL
  METHOD  QS

  &DFT
    BASIS_SET_FILE_NAME  BASIS_SET
    POTENTIAL_FILE_NAME  GTH_POTENTIALS

    &MGRID
      CUTOFF  300
      REL_CUTOFF  60
    &END MGRID

    &QS
      EPS_DEFAULT  1.0E-10
    &END QS

    &MIXING
      METHOD  BROYDEN_MIXING    ! Broyden 混合
      ALPHA  0.4
      NBUFFER  7
    &END MIXING

    &SCF
      SCF_GUESS  ATOMIC
      EPS_SCF  1.0E-6
      MAX_SCF  50
    &END SCF

    &XC
      &XC_FUNCTIONAL PADE       ! LDA (PADE)
      &END XC_FUNCTIONAL
    &END XC
  &END DFT

  &SUBSYS
    &CELL
      ABC  12.4 12.4 12.4       ! 大盒子 (Å)，足够大以避免自相互作用
      PERIODIC  NONE             ! 非周期性（气相分子）
    &END CELL

    &TOPOLOGY
      COORD_FILE_FORMAT  XYZ
    &END TOPOLOGY

    &KIND O
      BASIS_SET  DZVP-GTH-PADE
      POTENTIAL  GTH-PADE-q6
    &END KIND

    &KIND H
      BASIS_SET  DZVP-GTH-PADE
      POTENTIAL  GTH-PADE-q1
    &END KIND

    &COORD
      O   0.000   0.000   0.117   ! 氧原子
      H   0.000   0.757  -0.469   ! 氢原子 1
      H   0.000  -0.757  -0.469   ! 氢原子 2
    &END COORD
  &END SUBSYS
&END FORCE_EVAL

&MOTION
  &GEO_OPT
    OPTIMIZER  CG               ! 共轭梯度法
    MAX_ITER  200               ! 最大 200 步
    MAX_DR    1.0E-3            ! 最大位移收敛 (Bohr)
    RMS_DR    1.0E-3            ! RMS 位移收敛 (Bohr)
    MAX_FORCE  1.0E-3           ! 最大力收敛 (Ha/Bohr)
    RMS_FORCE  1.0E-3           ! RMS 力收敛 (Ha/Bohr)
  &END GEO_OPT
&END MOTION
```

**关键要点：**
- 初始结构的 O-H 键长和 H-O-H 键角只是近似值
- 优化后应得到约 0.97 Å 的 O-H 键长和约 104.5° 的键角
- 使用大盒子 (12.4 Å) 和 `PERIODIC NONE` 避免周期性镜像相互作用

### 8.2 H2O_geo_opt_fixed.inp — 约束几何优化

此文件演示如何固定氧原子，只优化两个氢原子的位置。

```
&GLOBAL
  PROJECT  H2O_geo_opt_fixed
  RUN_TYPE GEO_OPT
  PRINT_LEVEL LOW
&END GLOBAL

&FORCE_EVAL
  ! ... 与 H2O_geo_opt.inp 相同 ...
&END FORCE_EVAL

&MOTION
  &GEO_OPT
    OPTIMIZER  CG
    MAX_ITER  200
    MAX_DR    1.0E-3
    RMS_DR    1.0E-3
    MAX_FORCE  1.0E-3
    RMS_FORCE  1.0E-3
  &END GEO_OPT

  &CONSTRAINT
    &FIXED_ATOMS
      LIST  1                   ! ★ 固定第 1 号原子（氧原子）
      COMPONENTS_TO_FIX  XYZ   ! 固定 x, y, z 三个方向
    &END FIXED_ATOMS
  &END CONSTRAINT
&END MOTION
```

**与 H2O_geo_opt.inp 的区别：**

仅增加了 `&CONSTRAINT` / `&FIXED_ATOMS` 部分。

**物理解释：**
- 固定 O 原子后，两个 H 原子围绕固定的 O 原子重新排列
- 这模拟了在"表面"上（O 代表表面位点）H 原子的弛豫
- 优化后的 O-H 键长应与无约束情况相同（如果初始结构接近平衡）
- 但由于对称性破缺，两个 O-H 键长可能不完全相等

---

## 9. 练习

### 练习 1：比较优化算法

对 `H2O_geo_opt.inp`，分别使用 SD、CG、BFGS 优化，比较：
- 收敛所需的步数
- 每步的计算时间
- 最终能量是否一致

### 练习 2：初始结构的影响

修改 `H2O_geo_opt.inp` 中的初始坐标，使 O-H 键长偏离更多（如 1.5 Å），观察：
- 优化步数是否增加
- 最终结构是否相同

### 练习 3：收敛标准的影响

分别使用 1E-2, 1E-3, 1E-4, 1E-5 作为力收敛标准，比较：
- 收敛步数
- 最终能量差异
- 最终结构差异

### 练习 4：约束优化的物理

运行 `H2O_geo_opt_fixed.inp`，然后：
1. 比较约束和无约束优化的最终结构
2. 交换固定原子（固定 H，优化 O 和另一个 H），结果如何？
3. 只固定 O 的 z 坐标（`COMPONENTS_TO_FIX Z`），观察约束的效果

---

## 小结

本章学习了几何优化的核心概念和在 CP2K 中的实现：

| 概念 | CP2K 实现 |
|------|----------|
| 几何优化 | `RUN_TYPE GEO_OPT` |
| 优化算法 | `&GEO_OPT` OPTIMIZER (SD/CG/BFGS/L-BFGS) |
| 收敛判据 | MAX_DR, RMS_DR, MAX_FORCE, RMS_FORCE |
| 约束优化 | `&FIXED_ATOMS` LIST + COMPONENTS_TO_FIX |
| SCF 混合 | `&MIXING` METHOD, ALPHA, NBUFFER |
| 轨迹输出 | `project-pos-1.xyz` |

下一章将介绍**晶胞优化**——同时优化原子位置和晶胞参数。

---

**参考资源：**
- CP2K 手册 - &GEO_OPT：https://manual.cp2k.org/trunk/CP2K_INPUT/MOTION/GEO_OPT.html
- CP2K 手册 - &FIXED_ATOMS：https://manual.cp2k.org/trunk/CP2K_INPUT/MOTION/CONSTRAINT/FIXED_ATOMS.html
