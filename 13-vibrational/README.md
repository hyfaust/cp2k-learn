# 第13章：振动分析与红外光谱

> **难度：** ⭐⭐⭐  
> **关键词：** 振动分析、简正模式、红外光谱、Hessian矩阵、频率计算  
> **RUN_TYPE：** `VIBRATIONAL_ANALYSIS`

## 目录

1. [振动分析原理](#1-振动分析原理)
2. [有限差分法](#2-有限差分法)
3. [RUN_TYPE VIBRATIONAL_ANALYSIS](#3-run_type-vibrational_analysis)
4. [简正模式](#4-简正模式)
5. [对称性](#5-对称性)
6. [红外光谱](#6-红外光谱)
7. [结果解读](#7-结果解读)
8. [注意事项](#8-注意事项)
9. [输入文件详解](#9-输入文件详解)

---

## 1. 振动分析原理

### 1.1 Born-Oppenheimer势能面

在 Born-Oppenheimer 近似下，电子和原子核的运动被分离处理。电子在原子核构型固定的条件下求解薛定谔方程，得到势能面 (Potential Energy Surface, PES)：

$$E(\mathbf{R}) = E_{\text{elec}}(\mathbf{R}) + V_{\text{nn}}(\mathbf{R})$$

其中 $\mathbf{R}$ 代表所有原子核的坐标集合。势能面 $E(\mathbf{R})$ 包含了分子的全部静态信息——平衡构型对应势能面上的极小值点，振动性质则由极小值附近的曲率决定。

### 1.2 谐振近似

在平衡构型 $\mathbf{R}_0$ 附近，我们可以对势能面做 Taylor 展开：

$$E(\mathbf{R}) = E(\mathbf{R}_0) + \sum_i \left.\frac{\partial E}{\partial R_i}\right|_0 \Delta R_i + \frac{1}{2}\sum_{i,j} \left.\frac{\partial^2 E}{\partial R_i \partial R_j}\right|_0 \Delta R_i \Delta R_j + \cdots$$

在平衡构型处，一阶项为零（梯度为零），截断到二阶项即为**谐振近似 (Harmonic Approximation)**：

$$E(\mathbf{R}) \approx E(\mathbf{R}_0) + \frac{1}{2}\sum_{i,j} H_{ij} \Delta R_i \Delta R_j$$

其中 **Hessian 矩阵** 定义为：

$$H_{ij} = \frac{\partial^2 E}{\partial R_i \partial R_j}$$

这是一个 $3N \times 3N$ 的对称矩阵（$N$ 为原子数），矩阵元素是能量对原子坐标的二阶偏导数。

### 1.3 简正模式求解

为了求解振动频率，需要求解广义本征值方程：

$$\mathbf{H} \mathbf{L} = \mathbf{M} \mathbf{L} \boldsymbol{\Lambda}$$

其中：
- $\mathbf{H}$ 是质量加权 Hessian 矩阵
- $\mathbf{M}$ 是质量矩阵（对角阵，元素为原子质量）
- $\mathbf{L}$ 是本征向量矩阵（简正模式位移方向）
- $\boldsymbol{\Lambda}$ 是本征值对角阵

经过质量加权变换 $\tilde{H}_{ij} = H_{ij} / \sqrt{m_i m_j}$，本征值 $\lambda_k$ 与振动频率的关系为：

$$\nu_k = \frac{1}{2\pi}\sqrt{\lambda_k}$$

或换算为波数（cm⁻¹）：

$$\tilde{\nu}_k = \frac{1}{2\pi c}\sqrt{\lambda_k}$$

其中 $c$ 为光速。

---

## 2. 有限差分法

### 2.1 为什么使用有限差分？

CP2K 中的振动分析默认通过**有限差分法 (Finite Difference Method)** 计算 Hessian 矩阵。虽然理论上可以用解析二阶导数，但对于 GPW/GAPW 方法，有限差分法实现更简单且足够精确。

### 2.2 正向差分方案

CP2K 采用**正向差分 (Forward Difference)** 方案来计算 Hessian 矩阵元素：

$$H_{i\alpha, j\beta} = \frac{\partial^2 E}{\partial R_{i\alpha} \partial R_{j\beta}} \approx \frac{F_{i\alpha}(\mathbf{R}_0 - \Delta R_{j\beta} \hat{e}_{j\beta}) - F_{i\alpha}(\mathbf{R}_0)}{-\Delta R}$$

其中：
- $R_{i\alpha}$ 是第 $i$ 个原子在 $\alpha$ 方向（$x, y, z$）的坐标
- $\Delta R$ 是位移量（由 `DX` 参数控制）
- $F_{i\alpha}$ 是对应方向上的力分量

### 2.3 位移策略

对每个原子的每个方向（$x$, $y$, $z$）进行位移，总共有 $3N$ 个位移方向。对于非周期体系，还需扣除平动和转动自由度（6个或5个），但位移数仍为 $3N$。

如果利用对称性，实际需要的位移数可以显著减少。CP2K 会在计算前自动分析分子的点群对称性，只进行独立的位移操作，然后通过对称操作恢复完整的 Hessian 矩阵。

### 2.4 位移量 DX 的选择

`DX` 参数控制有限差分的位移步长，典型值为 0.01 bohr。选择时需平衡：

- **太小**：数值精度不足，力的微小差异被噪声淹没
- **太大**：高阶项贡献过大，偏离谐振近似

经验上，0.001–0.02 bohr 范围内结果通常对 DX 不敏感。

---

## 3. RUN_TYPE VIBRATIONAL_ANALYSIS

### 3.1 基本设置

在 CP2K 中执行振动分析，需要将 `RUN_TYPE` 设为 `VIBRATIONAL_ANALYSIS`：

```fortran
&GLOBAL
  PROJECT H2O_vib
  RUN_TYPE VIBRATIONAL_ANALYSIS
&END GLOBAL
```

### 3.2 &VIBRATIONAL_ANALYSIS 部分

振动分析的核心参数在 `&VIBRATIONAL_ANALYSIS` 部分设置：

```fortran
&VIBRATIONAL_ANALYSIS
  DX 0.01                          # 位移步长 (bohr)
  NPROC_REP 1                      # 每个位移使用的 MPI 进程数
  &PRINT
    &MOMENTS                       # 输出偶极矩（用于IR强度）
      PERIODIC .FALSE.             # 非周期体系
    &END MOMENTS
    &PROGRAM_RUN_INFO ON           # 输出详细运行信息
    &END
  &END PRINT
&END VIBRATIONAL_ANALYSIS
```

### 3.3 关键参数说明

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `DX` | 0.01 | 有限差分位移步长，单位 bohr |
| `INTENSITIES` | .FALSE. | 是否计算 IR 强度（需要偶极矩导数） |
| `TC` | .FALSE. | 是否计算热力学量（焓、熵、热容） |
| `THERMOCHEMISTRY` | .FALSE. | 同上 |
| `ROTATION` | .TRUE. | 是否扣除转动自由度 |
| `PRINT_LEVEL` | MEDIUM | 输出详细程度 |

### 3.4 计算流程

1. 首先对输入结构做几何优化（推荐在单独计算中完成，或确保结构已在极小值）
2. CP2K 分析分子对称性
3. 对每个独立位移方向，执行 SCF 计算得到能量/力
4. 组装 Hessian 矩阵
5. 质量加权对角化
6. 输出振动频率和简正模式

> **重要提示：** 振动分析要求输入结构是势能面上的驻点（至少是极小值）。如果结构不在极小值处，会出现虚频（负频率），这不代表真实振动模式。

---

## 4. 简正模式

### 4.1 自由度分析

一个由 $N$ 个原子组成的分子有 $3N$ 个振动自由度：

- **非线性分子：** $3N - 6$ 个振动模式（扣除 3 个平动 + 3 个转动）
- **线性分子：** $3N - 5$ 个振动模式（扣除 3 个平动 + 2 个转动）

例如，水分子（非线性，$N=3$）：$3 \times 3 - 6 = 3$ 个振动模式。

### 4.2 水分子的三个简正模式

| 模式 | 对称性 | 频率 (cm⁻¹) 实验值 | 描述 |
|------|--------|---------------------|------|
| $\nu_1$ | $A_1$ | 3657 | 对称伸缩 |
| $\nu_2$ | $A_1$ | 1595 | 弯曲（剪式） |
| $\nu_3$ | $B_2$ | 3756 | 反对称伸缩 |

### 4.3 简正模式的物理意义

每个简正模式是一个集体振动，所有原子以相同频率、同相位振动。模式的本征向量给出了每个原子的振动方向和相对振幅。

简正模式的特征：
- 不同简正模式之间完全独立（正交）
- 每个模式可以看作一个量子谐振子
- 零点能：$E_0 = \sum_k \frac{1}{2}\hbar\omega_k$

### 4.4 虚频

如果计算中出现**虚频**（频率为负值，CP2K 输出中显示为负数），通常意味着：

1. 结构不在势能面的极小值（最常见原因）
2. 位移步长太大，数值问题
3. 对称性约束导致的鞍点

对于过渡态搜索，有且仅有一个虚频对应反应坐标方向，是正常的。

---

## 5. 对称性

### 5.1 点群对称性

CP2K 可以利用分子的**点群对称性 (Point Group Symmetry)** 来减少振动分析的计算量。原理如下：

如果分子具有对称操作 $\{E, C_2, \sigma_v, \cdots\}$，则 Hessian 矩阵也必须满足这些对称性约束。利用对称性可以：

1. **减少独立位移数：** 对称等价方向只需计算一个
2. **提高精度：** 对称化后的 Hessian 矩阵数值误差更小

### 5.2 CP2K 中的对称性利用

CP2K 会自动检测输入结构的点群（最高到 $D_{2h}$ 子群）。对称性检测基于坐标的容差（`SYMMETRY_EPS` 参数）。

如果不想使用对称性（例如结构接近但不完全对称），可以关闭：

```fortran
&VIBRATIONAL_ANALYSIS
  DX 0.01
  &PRINT
    &MOMENTS
      PERIODIC .FALSE.
    &END MOMENTS
  &END PRINT
&END VIBRATIONAL_ANALYSIS
```

CP2K 的对称性利用是在内部自动完成的，无需用户额外设置。

### 5.3 对称性对 IR 选择定则的影响

对称性不仅节省计算，还决定了哪些振动模式是**红外活性 (IR active)** 的。一个振动模式具有 IR 活性的条件是：该模式的对称表示与偶极矩分量的对称表示有非零交集。

对于水分子（$C_{2v}$ 点群），所有三个模式（$2A_1 + B_2$）都是 IR 活性的。

---

## 6. 红外光谱

### 6.1 IR 强度的物理基础

红外光谱的强度正比于振动跃迁偶极矩的平方。在谐振近似下，第 $k$ 个模式的 IR 强度为：

$$I_k \propto \left|\frac{\partial \boldsymbol{\mu}}{\partial Q_k}\right|^2$$

其中 $\boldsymbol{\mu}$ 是分子偶极矩，$Q_k$ 是第 $k$ 个简正坐标。

$\partial \boldsymbol{\mu} / \partial Q_k$ 称为**偶极矩导数 (Dipole Derivative)**，反映了分子振动时电荷分布的变化。

### 6.2 CP2K 中计算 IR 强度

CP2K 通过对偶极矩做有限差分来计算 IR 强度。需要在 `&VIBRATIONAL_ANALYSIS` 的 `&PRINT` 部分启用偶极矩输出：

```fortran
&VIBRATIONAL_ANALYSIS
  DX 0.01
  INTENSITIES .TRUE.               # 启用 IR 强度计算
  &PRINT
    &MOMENTS
      PERIODIC .FALSE.             # 非周期体系使用 Berry phase 或直接求和
    &END MOMENTS
  &END PRINT
&END VIBRATIONAL_ANALYSIS
```

### 6.3 偶极矩的计算方式

对于**非周期体系**（分子），偶极矩通过对原子电荷和位置直接求和：

$$\boldsymbol{\mu} = \sum_i q_i \mathbf{r}_i$$

或使用密度矩阵积分：

$$\mu_\alpha = -\int r_\alpha \rho(\mathbf{r}) d\mathbf{r} + \sum_I Z_I R_{I,\alpha}$$

设置 `PERIODIC .FALSE.` 时，CP2K 使用直接求和方法。

对于**周期体系**，需要使用 Berry phase 方法（`PERIODIC .TRUE.`），但这通常用于固态 IR 光谱，不在本章讨论范围。

### 6.4 从偶极矩导数到 IR 强度

CP2K 通过数值微分获得偶极矩对简正坐标的导数：

$$\frac{\partial \mu_\alpha}{\partial Q_k} \approx \sum_i \frac{\partial \mu_\alpha}{\partial R_i} \cdot L_{ik}$$

其中 $L_{ik}$ 是简正模式位移向量。最终的 IR 强度（km/mol 单位）：

$$I_k = \frac{N_A \pi}{3c^2} \left|\frac{\partial \boldsymbol{\mu}}{\partial Q_k}\right|^2$$

---

## 7. 结果解读

### 7.1 频率输出

振动分析完成后，CP2K 会在输出文件中打印所有振动频率。典型的输出格式如下：

```
VIB|                        NORMAL MODES - CARTESIAN DISPLACEMENTS
VIB|
VIB|                     1         2         3         4         5         6
VIB|Frequency (cm^-1)   0.0000    0.0000    0.0000    0.0000    0.0000 1594.79
VIB|Intensities           0.000     0.000     0.000     0.000     0.000   68.42
VIB|Red. masses          ...
VIB|
VIB|                     7         8         9
VIB|Frequency (cm^-1) 3657.23  3756.12     ...
VIB|Intensities          5.21    42.87     ...
```

前 6 个模式（非线性分子）为零频率或近零频率，对应三个平动和三个转动自由度。

### 7.2 虚频的判读

- **频率 = 0 或非常接近 0**：平动或转动自由度（数值不精确可能导致微小的非零值）
- **频率 < 0（虚频）**：CP2K 通常用负数表示，有些程序用虚数表示。如果出现在优化后的结构中，说明优化不充分
- **频率 > 0**：真实的振动模式

### 7.3 模式可视化

简正模式可以可视化来理解振动特征。CP2K 输出的本征向量（位移向量）可以导入可视化软件（如 VMD、Jmol、GaussView）：

1. CP2K 输出的 `.mol` 文件（如果设置了 `&MOL_SET`）可直接读取
2. 或者将频率和位移向量转换为 Gaussian 格式的输出，用 VMD 等工具读取

### 7.4 与实验 IR 光谱的比较

比较计算与实验 IR 光谱时，需注意：

1. **频率标度因子：** DFT 计算的频率通常系统性偏高，需乘以标度因子（见注意事项章节）
2. **峰型：** 计算只给出离散的频率和强度，实验光谱有展宽。可以用 Lorentz 或 Gaussian 展宽函数模拟
3. **基质效应：** 实验可能在溶液或固态基质中进行，与气相计算有差异

### 7.5 热力学量

如果设置了 `TC .TRUE.`，CP2K 还会基于振动频率计算热力学量：

- **零点能 (ZPE)：** $E_0 = \sum_k \frac{1}{2}\hbar\omega_k$
- **热容 (Heat Capacity)：** $C_V = \sum_k R \left(\frac{\hbar\omega_k}{k_BT}\right)^2 \frac{e^{\hbar\omega_k/k_BT}}{(e^{\hbar\omega_k/k_BT}-1)^2}$
- **熵 (Entropy)：** 包含平动、转动、振动贡献
- **焓 (Enthalpy)：** $H = E + k_BT + \text{振动贡献}$

这些热力学量在计算反应能垒和自由能时非常有用。

---

## 8. 注意事项

### 8.1 非谐性 (Anharmonicity)

谐振近似是振动分析的首要近似。实际分子的势能面包含高阶项，导致：

- **振动频率偏移：** 真实频率通常低于谐振频率（对于大多数振动模式）
- **倍频和组合频：** 出现在非倍频位置（如 $2\nu_1$, $\nu_1 + \nu_2$ 等）
- **费米共振：** 某些频率接近的倍频与基频混合

DFT 谐振频率的典型误差（相对于实验值）：
- HF：偏高 ~12%
- B3LYP：偏高 ~5%
- PBE/PADE：偏高 ~3-5%

### 8.2 基组对频率的影响

基组的选择对振动频率有显著影响：

| 基组 | 特点 | 频率精度 |
|------|------|----------|
| SZV | 最小基组 | 较差，频率偏高 |
| DZVP | 双ζ+极化 | 合理，常用于初步计算 |
| TZV2P | 三ζ+双极化 | 较好 |
| cc-pVDZ | Dunning 相关一致 | 好 |
| cc-pVTZ | 更大 | 很好 |

经验表明，DZVP 级别的基组配合 GGA 泛函，经过标度因子校正后，可以得到与实验值偏差在 50 cm⁻¹ 以内的结果。

### 8.3 频率标度因子

由于谐振近似和方法的系统误差，计算频率需要乘以**标度因子 (Scaling Factor)**：

| 方法 | 标度因子 |
|------|----------|
| HF/STO-3G | 0.8929 |
| HF/6-31G(d) | 0.8905 |
| BLYP/DZVP | 0.9940 |
| B3LYP/6-31G(d) | 0.9614 |
| PBE/DZVP | ~0.97-0.98 |
| PADE/DZVP | ~0.97 |

标度因子来自 NIST 数据库的系统研究。使用时，将计算频率乘以标度因子后再与实验值比较。

### 8.4 几何优化的重要性

**振动分析前必须做充分的几何优化！** 推荐：

1. 先用较低精度快速优化
2. 再用与振动分析相同的方法和基组做高精度优化
3. 确认优化收敛后，再做振动分析
4. 检查梯度是否足够小（最大力 < 10⁻⁴ Hartree/bohr）

### 8.5 内存和计算量

振动分析的计算量约为单点计算的 $3N$ 倍（不利用对称性）或更少（利用对称性）。对于大分子：

- 计算时间随原子数 $N$ 线性增长（每个位移的 SCF 独立）
- 可以通过 `NPROC_REP` 控制并行度
- 考虑使用低精度方法（如 DFTB 或半经验方法）做初步振动分析

### 8.6 周期性体系的振动分析

对于周期性体系（晶体），振动分析需要特殊考虑：

- 没有真正的平动和转动自由度
- 使用 Gamma 点振动，对应声子的 $\Gamma$ 点
- 如果需要完整的声子色散关系，应使用 Phonopy + CP2K 的方案
- IR 强度的计算需要 Berry phase 方法求偶极矩

---

## 9. 输入文件详解

本章的示例文件 `H2O_vib.inp` 对水分子进行振动分析。以下逐节解释：

### 9.1 &GLOBAL 部分

```fortran
&GLOBAL
  PROJECT H2O_vib
  RUN_TYPE VIBRATIONAL_ANALYSIS
&END GLOBAL
```

- `PROJECT`：项目名称，所有输出文件以此为前缀
- `RUN_TYPE VIBRATIONAL_ANALYSIS`：指定计算类型为振动分析

### 9.2 &FORCE_EVAL 部分

```fortran
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME  BASIS_SET
    POTENTIAL_FILE_NAME  GTH_POTENTIALS
    &MGRIT ... &END        # 多网格设置
    &QS
      EPS_DEFAULT 1.0E-12  # 高精度积分阈值
    &END QS
    &XC
      &XC_FUNCTIONAL PADE   # LDA (PADE) 泛函
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 10.0 10.0 10.0   # 大盒子，避免镜像相互作用
    &END CELL
    &COORD
      O   0.000   0.000   0.117
      H   0.000   0.757  -0.469
      H   0.000  -0.757  -0.469
    &END COORD
    ...
  &END SUBSYS
&END FORCE_EVAL
```

关键点：
- 使用 Quickstep (GPW) 方法
- PADE (LDA) 泛函，计算速度快，适合演示
- DZVP-GTH 基组，双ζ+极化级别
- 盒子大小 10 Å³，足够隔离水分子

### 9.3 &VIBRATIONAL_ANALYSIS 部分

```fortran
&VIBRATIONAL_ANALYSIS
  DX 0.01                   # 位移步长 0.01 bohr
  &PRINT
    &MOMENTS
      PERIODIC .FALSE.      # 非周期体系
    &END MOMENTS
  &END PRINT
&END VIBRATIONAL_ANALYSIS
```

- `DX 0.01`：有限差分位移步长。0.01 bohr 是常用的默认值，在精度和数值稳定性之间取得平衡
- `&MOMENTS` + `PERIODIC .FALSE.`：输出偶极矩信息，用于计算 IR 强度。`PERIODIC .FALSE.` 指定使用直接求和方法计算偶极矩

### 9.4 运行与输出

```bash
# 运行计算
mpirun -np 4 cp2k.psmp -i H2O_vib.inp -o H2O_vib.out

# 查看振动频率
grep "Frequency" H2O_vib.out

# 查看 IR 强度
grep "Intensit" H2O_vib.out
```

预期输出：
- 9 个模式（$3N = 9$），前 6 个为平动/转动（近零频率）
- 3 个振动模式：
  - ~1600 cm⁻¹：弯曲模式
  - ~3650 cm⁻¹：对称伸缩
  - ~3750 cm⁻¹：反对称伸缩

由于使用 LDA/DZVP，频率会与实验值有偏差，但模式的排序和相对强度应与实验一致。

---

## 参考文献

1. P. Pulay, *Ab initio calculation of force constants and equilibrium geometries in polyatomic molecules*, Mol. Phys. **17**, 197 (1969).
2. J. Hutter et al., *Ab initio molecular dynamics and vibrational spectroscopy*, in *Reviews in Computational Chemistry*, Vol. 24.
3. The CP2K manual: https://manual.cp2k.org/
4. NIST CCCBDB: Computational Chemistry Comparison and Benchmark Database — 频率标度因子来源.
5. D. I. Bichoutskaya et al., *Vibrational frequencies scaling factors*, J. Phys. Chem. Ref. Data.

---

## 练习

1. **基础练习：** 运行 `H2O_vib.inp`，确认得到 3 个振动模式，比较频率与实验值的偏差。

2. **基组效应：** 将基组改为 TZV2P，比较频率变化。

3. **对称伸缩 vs 反对称伸缩：** 观察模式的位移向量，理解两种 O-H 伸缩模式的区别。

4. **同位素效应：** 将 H 替换为 D（氘），重新计算频率。验证频率比是否接近 $\sqrt{m_H / m_D} \approx 0.707$。

5. **热力学量：** 添加 `TC .TRUE.`，计算 298.15 K 下水分子的零点能和热容。
