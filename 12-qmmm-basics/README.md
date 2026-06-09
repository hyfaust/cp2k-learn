# 第12章：QM/MM混合方法基础

> **难度：** ⭐⭐⭐⭐  
> **前置知识：** 第2章（DFT基础）、第8-9章（分子动力学）、第11章（半经验方法）  
> **预计学习时间：** 3-4小时

## 目录

- [QM/MM的基本思想](#qmmm的基本思想)
- [嵌入方案](#嵌入方案)
- [力场基础](#力场基础)
- [连接原子方法](#连接原子方法)
- [CP2K中的QM/MM设置](#cp2k中的qmmm设置)
- [示例：水溶液中的QM/MM计算](#示例水溶液中的qmmm计算)
- [QM/MM的应用场景](#qmmm的应用场景)
- [练习](#练习)

---

## QM/MM的基本思想

### 为什么需要QM/MM

在实际的化学和生物学体系中，我们经常遇到这样的问题：

- **酶催化反应** — 反应中心需要量子力学描述，但整个酶（数千原子）太昂贵
- **溶液中的化学反应** — 溶质需要QM，溶剂分子用MM足够
- **表面吸附** — 吸附物需要QM，大块材料用MM

QM/MM（Quantum Mechanics/Molecular Mechanics）混合方法的核心思想是：

> **将体系分为两个区域：QM区域用量子力学精确处理，MM区域用经典力场近似处理。**

```
┌──────────────────────────────────┐
│         MM 区域（力场）          │
│    ┌──────────────────────┐      │
│    │    QM 区域（DFT）    │      │
│    │    ●───● 反应中心    │      │
│    │    ●                 │      │
│    └──────────────────────┘      │
│  ○ ○ ○ ○ ○ ○ ○ ○ ○ ○ ○ ○ ○ ○  │
└──────────────────────────────────┘
```

### 多尺度建模的哲学

QM/MM是多尺度建模的一个范例：

| 尺度 | 方法 | 体系大小 | 时间尺度 |
|------|------|----------|----------|
| 电子 | DFT, MP2 | ~100原子 | ~fs |
| 原子 | 力场 (AMBER, CHARMM) | ~100,000原子 | ~ns-μs |
| 粗粒化 | MARTINI等 | ~1,000,000原子 | ~μs-ms |

QM/MM桥接了电子尺度和原子尺度。

---

## 嵌入方案

### 机械嵌入（Mechanical Embedding）

最简单的方案：QM和MM区域只通过非键力场相互作用，没有静电耦合。

$$E_{\text{total}} = E_{\text{QM}} + E_{\text{MM}} + E_{\text{QM-MM}}^{\text{vdW}}$$

- QM区域的电子结构不受MM环境影响
- 只考虑QM-MM之间的范德华相互作用
- 精度较低，但实现简单

**CP2K设置：** `E_COUPL NONE`

### 静电嵌入（Electrostatic Embedding）

更精确的方案：MM原子的电荷被包含在QM哈密顿量中。

$$E_{\text{total}} = E_{\text{QM+MM电荷}} + E_{\text{MM}} + E_{\text{QM-MM}}^{\text{vdW}}$$

- QM区域的电子密度被MM点电荷极化
- 更准确地描述环境对反应中心的影响
- 是最常用的QM/MM方案

**CP2K设置：** `E_COUPL COULOMB`

CP2K使用GEEP（Gaussian Expansion of the Electrostatic Potential）方法将MM点电荷展宽为高斯函数，以避免QM区域的电子密度坍缩到MM点电荷上。

### 极化嵌入（Polarizable Embedding）

最先进的方案：MM区域也可以被QM区域极化。

- MM原子带有诱导偶极
- 双向极化耦合
- 计算成本最高

CP2K支持极化力场，通过 `POLAR` 关键字设置。

---

## 力场基础

### 什么是力场

力场是原子间相互作用的参数化数学模型。典型的力场势能函数为：

$$E_{\text{FF}} = \underbrace{\sum_{\text{bonds}} k_r(r - r_0)^2}_{\text{键伸缩}} + \underbrace{\sum_{\text{angles}} k_\theta(\theta - \theta_0)^2}_{\text{键角弯曲}} + \underbrace{\sum_{\text{dihedrals}} V_n[1 + \cos(n\phi - \gamma)]}_{\text{二面角扭转}} + \underbrace{\sum_{i<j}\left[\frac{A_{ij}}{r_{ij}^{12}} - \frac{B_{ij}}{r_{ij}^6} + \frac{q_iq_j}{r_{ij}}\right]}_{\text{非键相互作用}}$$

### 常用力场

| 力场 | 主要应用 | 来源 |
|------|----------|------|
| AMBER | 蛋白质、核酸、药物 | AMBER社区 |
| CHARMM | 蛋白质、脂质、糖类 | CHARMM社区 |
| OPLS-AA | 有机液体、蛋白质 | Jorgensen组 |
| GROMOS | 生物分子 | GROMOS社区 |

### CP2K中使用力场

```
&MM
  &FORCEFIELD
    PARM_FILE_NAME amber.prmtop
    PARM_TYPE AMBER
    &SPLINE
      EMAX_SPLINE 1.0
      RCUT_NB 12.0
    &END SPLINE
  &END FORCEFIELD
  &POISSON
    &EWALD
      EWALD_TYPE SPME
      ALPHA 0.35
      GMAX 80
    &END EWALD
  &END POISSON
&END MM
```

---

## 连接原子方法

### QM/MM边界问题

当QM和MM区域之间有共价键被切断时，需要处理"悬挂键"问题：

```
QM区域          MM区域
  ...─C─C─|─C─C─...
         ↑
       切断位置
```

解决方案：在切断位置添加**连接原子（Link Atom）**，通常是氢原子：

```
QM区域              MM区域
  ...─C─C─H  C─C─...
         ↑
      连接原子(H)
```

### CP2K中的连接原子设置

```
&QMMM
  &LINK
    QM_INDEX 15       ! QM侧原子索引
    MM_INDEX 16       ! MM侧原子索引
    LINK_TYPE IMOMM   ! 连接原子类型
    ALPHA 0.77        ! 连接原子位置缩放因子
  &END LINK
&END QMMM
```

**LINK_TYPE选择：**
- **IMOMM** — 最常用，连接原子放在QM原子和MM原子的连线上
- **MOMM** — 使用分子力学修正
- **PSEUDO** — 使用赝势替代

**边界切割原则：**
- 尽量在C-C、C-N等非极性键处切割
- 避免在C=O、N-H等极性键处切割
- 连接原子应放在远离反应中心的位置

---

## CP2K中的QM/MM设置

### 完整的QM/MM输入结构

```
&GLOBAL
  PROJECT qmmm_system
  RUN_TYPE MD
&END GLOBAL

&FORCE_EVAL
  METHOD QMMM          ! 使用QM/MM方法

  &DFT                 ! QM部分的DFT设置
    BASIS_SET_FILE_NAME BASIS_SET
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &QS
      EPS_DEFAULT 1.0E-10
    &END QS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-6
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT

  &MM                  ! MM部分的力场设置
    &FORCEFIELD
      PARM_FILE_NAME system.prmtop
    &END FORCEFIELD
    &POISSON
      &EWALD
        EWALD_TYPE SPME
      &END EWALD
    &END POISSON
  &END MM

  &QMMM               ! QM/MM耦合设置
    E_COUPL COULOMB    ! 静电嵌入
    &CELL
      ABC 15.0 15.0 15.0
    &END CELL
    &QM_KIND C         ! 哪些C原子在QM区域
      MM_INDEX 1 2 3 4
    &END QM_KIND
    &QM_KIND O
      MM_INDEX 5
    &END QM_KIND
    &QM_KIND H
      MM_INDEX 6 7 8 9 10
    &END QM_KIND
    &LINK
      QM_INDEX 4
      MM_INDEX 11
      LINK_TYPE IMOMM
    &END LINK
  &END QMMM

  &SUBSYS
    &CELL
      ABC 50.0 50.0 50.0    ! 整个体系的胞
    &END CELL
    &TOPOLOGY
      PARM_FILE_NAME system.prmtop
    &END TOPOLOGY
  &END SUBSYS
&END FORCE_EVAL
```

### 关键设置说明

**METHOD QMMM：** 告诉CP2K使用QM/MM模式而不是纯QM或纯MM。

**E_COUPL：** 
- `COULOMB` — 静电嵌入（推荐）
- `NONE` — 机械嵌入

**&QM_KIND节：** 定义哪些MM原子被提升到QM级别。每个&QM_KIND节对应一种元素，`MM_INDEX` 列出该元素在QM区域的所有原子的MM拓扑索引。

**&CELL节（在QMMM下）：** QM区域的胞大小，通常比MM胞小得多。

---

## 示例：水溶液中的QM/MM计算

### 系统描述

将一个水分子放在QM级别（使用DFT），周围3 Å范围内的水分子使用MM力场描述：

### 输入文件：`water_qmmm.inp`

```
&GLOBAL
  PROJECT water_qmmm
  RUN_TYPE ENERGY
  PRINT_LEVEL LOW
&END GLOBAL

&FORCE_EVAL
  METHOD QMMM

  &DFT
    BASIS_SET_FILE_NAME BASIS_SET
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &QS
      EPS_DEFAULT 1.0E-10
    &END QS
    &MGRID
      CUTOFF 300
      NGRIDS 4
      REL_CUTOFF 60
    &END MGRID
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-6
      MAX_SCF 200
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.4
        NBROYDEN 8
      &END MIXING
    &END SCF
    &XC
      &XC_FUNCTIONAL PADE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT

  &MM
    &FORCEFIELD
      PARM_FILE_NAME water_ff.pot
      PARM_TYPE USER  ! 或AMBER
      &SPLINE
        EMAX_SPLINE 1.0
        RCUT_NB 12.0
      &END SPLINE
    &END FORCEFIELD
    &POISSON
      &EWALD
        EWALD_TYPE SPME
        ALPHA 0.35
        GMAX 80
      &END EWALD
    &END POISSON
  &END MM

  &QMMM
    E_COUPL COULOMB
    USE_GEEP_LIB 7
    &CELL
      ABC 10.0 10.0 10.0
    &END CELL
    &QM_KIND O
      MM_INDEX 1
    &END QM_KIND
    &QM_KIND H
      MM_INDEX 2 3
    &END QM_KIND
  &END QMMM

  &SUBSYS
    &CELL
      ABC 20.0 20.0 20.0
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME water_box.xyz
      COORD_FILE_FORMAT XYZ
    &END TOPOLOGY
  &END SUBSYS
&END FORCE_EVAL
```

### 输出解读

QM/MM计算的输出包含：
- **QM能量** — DFT计算的QM区域能量
- **MM能量** — 力场计算的MM区域能量
- **QM-MM耦合能量** — 两个区域之间的静电和vdW相互作用
- **总能量** — 所有贡献之和

---

## QM/MM的应用场景

### 酶催化反应

最经典的应用：
1. 用AMBER力场处理整个酶（~10,000原子）
2. 反应底物和关键残基用DFT处理（~50-100原子）
3. 使用元动力学计算反应自由能面

### 溶液化学

研究溶液中的化学反应：
1. 溶质分子用DFT
2. 溶剂分子用MM力场
3. 需要足够的溶剂层来屏蔽边界效应

### 材料表面

研究表面吸附和反应：
1. 吸附物和表面活性位用DFT
2. 大块材料用MM或简单模型

### 蛋白质-配体相互作用

药物设计中：
1. 药物分子和蛋白质活性位用DFT
2. 蛋白质其余部分和溶剂用MM

---

## 练习

### 练习1：QM/MM基础

构建一个包含5个水分子的系统，将中心水分子设为QM（DFT/PADE），其余4个设为MM。计算总能量并与纯QM和纯MM结果对比。

### 练习2：嵌入方案对比

对同一系统分别使用机械嵌入（E_COUPL NONE）和静电嵌入（E_COUPL COULOMB），比较QM区域水分子的O-H键长和偶极矩差异。

### 练习3：QM区域大小的影响

逐步增大QM区域（1个、2个、3个水分子），观察QM区域水分子几何构型的变化趋势。

### 练习4：QM/MM MD

对QM/MM水系统进行10ps的NVT-MD，分析QM水分子与MM水分子之间的径向分布函数。

---

## 参考文献

1. Warshel, A., & Levitt, M. (1976). Theoretical studies of enzymic reactions. *J. Mol. Biol.*, 103, 227-249.
2. Lonsdale, R., Harvey, J. N., & Mulholland, A. J. (2012). A practical guide to modelling enzyme-catalysed reactions. *Chem. Soc. Rev.*, 41, 3025-3038.
3. CP2K Documentation — QM/MM: https://www.cp2k.org/docs
