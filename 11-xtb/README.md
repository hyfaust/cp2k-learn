# 第11章：GFN1-xTB半经验方法

> **难度：** ⭐⭐⭐  
> **前置知识：** 第1-2章（CP2K基础）、第6章（几何优化）  
> **预计学习时间：** 2-3小时

## 目录

- [什么是半经验方法](#什么是半经验方法)
- [GFN1-xTB理论基础](#gfn1-xtb理论基础)
- [CP2K中使用xTB](#cp2k中使用xtb)
- [示例1：乙醇分子优化](#示例1乙醇分子优化)
- [示例2：水盒子MD模拟](#示例2水盒子md模拟)
- [xTB vs DFT对比](#xtb-vs-dft对比)
- [自定义非键相互作用](#自定义非键相互作用)
- [性能优化技巧](#性能优化技巧)
- [练习](#练习)

---

## 什么是半经验方法

### 为什么需要半经验方法

在前面的章节中，我们使用了密度泛函理论（DFT）进行各种计算。DFT虽然比波函数方法便宜得多，但对于以下场景仍然可能太慢：

- **大体系**（数千原子的蛋白质、材料表面）
- **长时间MD**（纳秒级模拟需要数百万步）
- **高通量筛选**（快速评估大量候选结构）
- **预平衡**（在昂贵的DFT MD之前快速弛豫）

半经验量子化学方法通过以下策略大幅降低计算成本：

1. **简化电子积分** — 使用参数化的解析公式替代精确的多中心积分
2. **冻结内层电子** — 只处理价电子
3. **经验参数** — 用实验数据或高级别计算拟合参数
4. **简化哈密顿量** — 使用零阶或低阶展开

### 半经验方法的层次

```
精度 ↑
  │  CCSD(T) ── 波函数方法
  │  MP2      ── 后HF方法
  │  DFT      ── 密度泛函理论
  │  xTB/DFTB ── 紧束缚方法
  │  AM1/PM3  ── NDDO方法
  │  Force Field ── 经典力场
  └───────────────────────→ 速度
```

CP2K支持的半经验方法包括：
- **GFN1-xTB** — Grimme组开发的通用紧束缚方法
- **DFTB** — 密度泛函紧束缚方法（需要专门的参数文件）
- **AM1, PM3, PM6, RM1, MNDO** — NDDO族方法

本章重点介绍 **GFN1-xTB**，因为它无需额外参数文件，开箱即用。

---

## GFN1-xTB理论基础

### GFN1-xTB能量表达式

GFN1-xTB（Geometry, Frequency, Non-covalent interactions, extended Tight Binding）的总能量为：

$$E_{\text{total}} = E_{\text{el}} + E_{\text{rep}} + E_{\text{disp}} + E_{\text{XB}}$$

#### 1. 电子能量 $E_{\text{el}}$

电子能量使用零阶哈密顿量 $H^0$ 的期望值加上二阶和三阶修正：

$$E_{\text{el}} = \sum_i n_i \langle \psi_i | H^0 | \psi_i \rangle + E^{(2)} + E^{(3)}$$

其中：
- $n_i$ 是轨道占据数（通过Fermi展宽确定）
- $H^0$ 包含动能、核-电子吸引、有效电子-电子排斥
- $E^{(2)}$ 是二阶电子排斥修正，使用原子对相关的 $\gamma_{AB}$ 函数
- $E^{(3)}$ 是三阶Mulliken电荷修正

**关键区别于DFT：** xTB不计算四中心电子积分，而是使用参数化的解析公式近似所有电子相互作用。

#### 2. 排斥能 $E_{\text{rep}}$

原子对间的短程排斥源于核心-核心排斥和Pauli排斥：

$$E_{\text{rep}} = \sum_{A<B} Z_A^{\text{eff}} Z_B^{\text{eff}} \cdot \alpha_{AB} \cdot \exp(-\beta_{AB} \cdot r_{AB})$$

这是一个指数衰减的原子对势，参数 $Z^{\text{eff}}$ 是有效核电荷，$\alpha$ 和 $\beta$ 是拟合参数。

#### 3. 色散能 $E_{\text{disp}}$

使用D3方法（Grimme色散校正）的BJ阻尼版本：

$$E_{\text{disp}} = -\sum_{n=6,8} s_n \sum_{A<B} \frac{C_n^{AB}}{r_{AB}^n + (R_0^{AB})^n}$$

其中 $C_6^{AB}$ 和 $C_8^{AB}$ 是色散系数，$R_0$ 是阻尼距离。

#### 4. 卤键校正 $E_{\text{XB}}$

针对卤素键（Cl, Br与O, N之间的相互作用）的专门校正：

$$E_{\text{XB}} = -\sum_{A \in \text{halogen}} \sum_{B} \frac{C_6^{AB}}{r_{AB}^6} \cdot f_{\text{damp}}(r_{AB})$$

通过 `USE_HALOGEN_CORRECTION` 关键字控制。

### 与DFT的关键区别

| 特性 | DFT | GFN1-xTB |
|------|-----|----------|
| 基组 | 需要BASIS_SET文件 | 内置，无需外部文件 |
| 赝势 | 需要GTH_POTENTIALS文件 | 不需要 |
| 电子积分 | 精确四中心积分 | 参数化近似 |
| XC泛函 | 需要选择 | 内置参数化 |
| 速度 | 中等 | 快10-100倍 |
| 精度 | 较高 | 中等（依赖参数化） |
| 适用体系 | 几乎所有 | 主族元素（H-Rn） |

---

## CP2K中使用xTB

### 基本设置

在CP2K中启用xTB非常简单，不需要任何外部参数文件：

```
&QS
  METHOD XTB
&END QS
```

仅此一行！CP2K会自动使用内置的GFN1-xTB参数。

### 关键输入节

#### `&QS` 节 — 方法选择

```
&QS
  METHOD XTB
  &XTB
    DO_EWALD .TRUE.             ! 长程静电使用Ewald求和
    USE_HALOGEN_CORRECTION .TRUE. ! 卤键校正
    CHECK_ATOMIC_CHARGES .TRUE.  ! 检查Mulliken电荷合理性
  &END XTB
&END QS
```

#### `&OT` 节 — 轨道变换

xTB使用OT（Orbital Transformation）方法进行SCF求解，而不是传统的对角化：

```
&SCF
  &OT
    MINIMIZER DIIS              ! DIIS最小化器
    PRECONDITIONER FULL_ALL     ! 预条件器
  &END OT
  MAX_SCF 50
  EPS_SCF 1.0E-6
&END SCF
```

#### 拓扑设置

xTB不需要键连信息，因此应关闭连通性检测：

```
&TOPOLOGY
  CONNECTIVITY OFF
&END TOPOLOGY
```

### xTB的SCF收敛

xTB的SCF收敛通常比DFT快得多，因为：
1. 哈密顿矩阵更简单
2. 没有复杂的XC泛函
3. 初始猜测更准确

但偶尔也会遇到收敛困难，此时可以使用 `&OUTER_SCF` 节控制：

```
&SCF
  &OUTER_SCF
    MAX_SCF 20
    EPS_SCF 1.0E-6
    TYPE NONE                   ! 或 MULTI_P_GEM
  &END OUTER_SCF
&END SCF
```

---

## 示例1：乙醇分子优化

### 输入文件：`ethanol_xtb.inp`

```
&GLOBAL
  PROJECT ethanol_xtb
  RUN_TYPE GEO_OPT
  PRINT_LEVEL LOW
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep

  &DFT
    &QS
      METHOD XTB
      &XTB
        DO_EWALD .FALSE.
        USE_HALOGEN_CORRECTION .TRUE.
        CHECK_ATOMIC_CHARGES .TRUE.
      &END XTB
    &END QS

    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-6
      MAX_SCF 50
      &OT
        MINIMIZER DIIS
        PRECONDITIONER FULL_ALL
      &END OT
    &END SCF
  &END DFT

  &SUBSYS
    &CELL
      ABC 15.0 15.0 15.0
    &END CELL
    &COORD
      C    0.000   0.000   0.000
      C    1.520   0.000   0.000
      O    2.020   1.340   0.000
      H    2.920   1.340   0.000
      H   -0.510   0.510   0.890
      H   -0.510   0.510  -0.890
      H   -0.510  -1.020   0.000
      H    1.920  -0.530   0.890
      H    1.920  -0.530  -0.890
    &END COORD
  &END SUBSYS
&END FORCE_EVAL

&MOTION
  &GEO_OPT
    TYPE MINIMIZATION
    OPTIMIZER BFGS
    MAX_ITER 200
    MAX_DR    1.0E-03
    MAX_FORCE 1.0E-03
    RMS_DR    1.0E-03
    RMS_FORCE 1.0E-03
  &END GEO_OPT
&END MOTION
```

### 逐行解析

**GLOBAL节：**
- `RUN_TYPE GEO_OPT` — 几何优化
- 不需要BASIS_SET_FILE_NAME和POTENTIAL_FILE_NAME

**DFT/QS节：**
- `METHOD XTB` — 使用GFN1-xTB方法
- `DO_EWALD .FALSE.` — 对于孤立分子关闭Ewald求和
- `CHECK_ATOMIC_CHARGES .TRUE.` — 检查Mulliken电荷是否合理

**SCF节：**
- 使用OT方法而非对角化
- `PRECONDITIONER FULL_ALL` — 使用全矩阵预条件器

**SUBSYS节：**
- 只需要坐标和胞参数，**不需要指定基组和赝势**
- 拓扑中应设置 `CONNECTIVITY OFF`

### 运行和输出

```bash
cp2k.popt -o ethanol_xtb.out ethanol_xtb.inp
```

输出中需要关注：
- **总能量** — 以Hartree为单位，但绝对值不如DFT精确
- **优化构型** — 键长、键角应与实验值合理一致
- **Mulliken电荷** — 检查电荷分布是否合理
- **收敛步数** — xTB通常比DFT收敛更快

---

## 示例2：水盒子MD模拟

### 输入文件：`water_xtb_md.inp`

xTB的一个主要应用场景是大体系的长时间分子动力学。以下示例模拟64个水分子：

```
&GLOBAL
  PROJECT water_xtb_md
  RUN_TYPE MD
  PRINT_LEVEL LOW
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep

  &DFT
    &QS
      METHOD XTB
      &XTB
        DO_EWALD .TRUE.
      &END XTB
      EXTRAPOLATION ASPC
      EXTRAPOLATION_ORDER 3
    &END QS

    &SCF
      SCF_GUESS RESTART
      EPS_SCF 1.0E-4
      MAX_SCF 20
      &OT
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
    &END SCF
  &END DFT

  &SUBSYS
    &CELL
      ABC 12.4138 12.4138 12.4138
    &END CELL
    &TOPOLOGY
      CONNECTIVITY OFF
      COORD_FILE_FORMAT XYZ
      COORD_FILE_NAME water64.xyz
    &END TOPOLOGY
  &END SUBSYS
&END FORCE_EVAL

&MOTION
  &MD
    ENSEMBLE NVT
    STEPS 10000
    TIMESTEP 1.0
    TEMPERATURE 300.0
    &THERMOSTAT
      TYPE CSVR
      &CSVR
        TIMECON 100.0
      &END CSVR
    &END THERMOSTAT
  &END MD

  &PRINT
    &TRAJECTORY
      FORMAT XYZ
      &EACH
        MD 10
      &END EACH
    &END TRAJECTORY
    &VELOCITIES OFF
    &FORCES OFF
    &RESTART
      &EACH
        MD 1000
      &END EACH
    &END RESTART
  &END PRINT
&END MOTION
```

### 关键设置说明

**EPS_SCF放宽：** 对于MD，SCF收敛精度可以适当放宽（1.0E-4 vs DFT的1.0E-7），因为每一步的力只要求统计正确，不需要精确到机器精度。

**ASPC外推：** 对于NVT/NPT MD，ASPC外推方法配合放宽的SCF收敛可以显著加速计算，同时保持足够的精度。

**CONNECTIVITY OFF：** 关闭键连检测，xTB不需要此信息。

### xTB MD的优势

| 指标 | DFT (PBE/DZVP) | GFN1-xTB |
|------|----------------|----------|
| 单步耗时 | ~10-30秒 | ~0.5-2秒 |
| 10ps模拟 | ~10-30小时 | ~1-4小时 |
| 最大体系 | ~100-500原子 | ~5000-10000原子 |
| 精度 | 较高 | 中等 |

---

## xTB vs DFT对比

### 精度对比

GFN1-xTB在以下性质上表现良好：

- **几何构型** — 键长误差通常 < 0.03 Å，键角误差 < 3°
- **相对能量** — 异构体能量排序通常正确，但绝对值偏差较大
- **振动频率** — 通常偏高约5-10%，需要缩放因子
- **反应能垒** — 定性正确，定量可能偏差较大

### 适用场景

**适合使用xTB的情况：**
- 大体系（>500原子）的初步研究
- 长时间MD的预平衡
- 构象搜索和构象采样
- 高通量筛选
- 热浴化预处理

**不适合使用xTB的情况：**
- 需要精确电子结构（带隙、激发能）
- 含过渡金属的精确计算
- 需要精确热化学数据
- 弱相互作用精确定量

---

## 自定义非键相互作用

CP2K的xTB实现支持通过 `&GENPOT` 添加用户自定义的非键相互作用势：

```
&QS
  METHOD XTB
  &XTB
    DO_NONBONDED .TRUE.
    &NONBOND
      &GENPOT
        ATOMS O Zn
        FUNCTION A*exp(-B*r)
        PARAMETERS A B
        VALUES 0.5 2.0
        VARIABLES r
      &END GENPOT
    &END NONBOND
  &END XTB
&END QS
```

这允许您在xTB的基础上添加系统特定的修正势，例如：
- 特定原子对的额外排斥势
- 表面吸附的修正势
- 自定义约束势

---

## 性能优化技巧

### 1. SCF收敛加速

```
&SCF
  &OT
    MINIMIZER DIIS
    PRECONDITIONER FULL_SINGLE_INVERSE  ! 对大体系更省内存
    ENERGY_GAP 0.1                      ! xTB能隙估计
  &END OT
&END SCF
```

### 2. 并行化

xTB的计算量主要在矩阵构建和OT最小化：
- MPI并行：`mpirun -n N cp2k.popt ...`
- xTB的并行效率通常比DFT好（因为计算更简单）

### 3. MD中的SCF步数

```
&SCF
  MAX_SCF 10           ! MD中通常10步内收敛
  EPS_SCF 1.0E-4       ! 适当放宽收敛标准
  &OT
    MINIMIZER DIIS
  &END OT
&END SCF
```

### 4. 输出控制

```
&MOTION
  &PRINT
    &TRAJECTORY OFF       ! 减少I/O
    &ENERGY
      &EACH
        MD 100            ! 每100步输出一次能量
      &END EACH
    &END ENERGY
  &END PRINT
&END MOTION
```

---

## 练习

### 练习1：xTB几何优化

对以下分子进行xTB几何优化，比较优化后的键长与实验值：
- 甲烷 CH₄
- 乙烷 C₂H₆
- 苯 C₆H₆

### 练习2：xTB vs DFT能量对比

对水二聚体分别进行xTB和DFT(PBE/DZVP)计算，比较：
- 氢键长度
- 结合能
- 计算时间

### 练习3：大体系MD

构建一个包含256个水分子的盒子（使用Packmol或其他工具），运行100ps的NVT-MD，分析：
- 径向分布函数 g(r)
- 扩散系数

### 练习4：自定义势

为一个蛋白质-配体体系添加自定义约束势，限制配体在活性位点附近。

---

## 参考文献

1. Grimme, S., Bannwarth, C., & Shushkov, P. (2017). A Robust and Accurate Tight-Binding Quantum Chemical Method for Structures, Vibrational Frequencies, and Noncovalent Interactions of Large Molecular Systems. *J. Chem. Theory Comput.*, 13, 1989-2009.
2. Bannwarth, C., Ehlert, S., & Grimme, S. (2019). GFN2-xTB—An Accurate and Broadly Parametrized Self-Consistent Tight-Binding Quantum Chemical Method with Multipole Electrostatics and Density-Dependent Dispersion Contributions. *J. Chem. Theory Comput.*, 15, 1652-1671.
3. CP2K Documentation — Semi-Empirical Methods: https://www.cp2k.org/docs
