# 第14章：TDDFT 激发态计算

> **难度：** ⭐⭐⭐⭐  
> **关键词：** TDDFT、激发态、Casida方程、吸收光谱、振子强度  
> **RUN_TYPE：** `ENERGY`（TDDFT 在 `&TDDFT` 部分设置）

## 目录

1. [激发态理论基础](#1-激发态理论基础)
2. [Casida 方程](#2-casida-方程)
3. [交换关联核 f_xc](#3-交换关联核-f_xc)
4. [CP2K 中的 TDDFT](#4-cp2k-中的-tddft)
5. [吸收光谱](#5-吸收光谱)
6. [适用范围](#6-适用范围)
7. [自旋翻转 TDDFT](#7-自旋翻转-tddft)
8. [输入文件详解](#8-输入文件详解)

---

## 1. 激发态理论基础

### 1.1 基态与激发态

Kohn-Sham DFT 是一种基态理论——它精确地（原则上）给出系统的基态能量和电子密度，但 KS 轨道能量并不直接对应激发能。理解这一点非常重要：

**KS 物理 vs 真实物理：**

| 特性 | Kohn-Sham 体系 | 真实体系 |
|------|----------------|----------|
| 基态能量 | 精确（原则上） | 精确 |
| 基态密度 | 精确（原则上） | 精确 |
| 轨道能量 $\varepsilon_i$ | Lagrange 乘子 | 无直接物理意义 |
| HOMO 能量 | $-I$（电离能） | 仅在精确 XC 泛函下成立 |
| LUMO-HOMO 差 | **不等于** 带隙/激发能 | — |

KS 带隙（$\varepsilon_{\text{LUMO}} - \varepsilon_{\text{HOMO}}$）通常**严重低估**真实带隙，原因在于 XC 泛函的**导数不连续性 (Derivative Discontinuity)**：

$$E_{\text{gap}}^{\text{true}} = E_{\text{gap}}^{\text{KS}} + \Delta_{xc}$$

其中 $\Delta_{xc}$ 是交换关联势的导数不连续项，LDA/GGA 中缺失。

### 1.2 Time-Dependent DFT (TDDFT)

TDDFT 是计算激发态的标准第一性原理方法，基于 **Runge-Gross 定理**（1984）：

> 对于给定初始态的多体系统，含时外势 $v_{\text{ext}}(\mathbf{r}, t)$ 与含时电子密度 $\rho(\mathbf{r}, t)$ 之间存在一一对应关系。

这意味着，原则上所有含时性质都可以从密度获得。

**绝热近似 (Adiabatic Approximation)：** TDDFT 的实际应用中最关键的近似。它假设含时 XC 势只依赖于**当前时刻**的密度（而非历史）：

$$v_{xc}^{\text{adiabatic}}[\rho](\mathbf{r}, t) = \left.\frac{\delta E_{xc}[\rho]}{\delta \rho(\mathbf{r})}\right|_{\rho = \rho(t)}$$

### 1.3 线性响应 TDDFT

在实际计算中，我们通常不需要真正求解含时 KS 方程。对于弱外场（如光吸收），**线性响应理论**已经足够：

1. 系统处于基态 $\rho_0$
2. 受到微扰 $\delta v_{\text{ext}}(\mathbf{r}, \omega)$（振荡电场）
3. 密度响应 $\delta\rho(\mathbf{r}, \omega) = \int \chi(\mathbf{r}, \mathbf{r}', \omega) \delta v_{\text{ext}}(\mathbf{r}', \omega) d\mathbf{r}'$

其中 $\chi$ 是**密度-密度响应函数**。通过 Dyson 方程将 $\chi$ 与非相互作用 KS 体系的响应函数 $\chi_0$ 关联：

$$\chi = \chi_0 + \chi_0 (v_c + f_{xc}) \chi$$

这里 $v_c$ 是库仑核，$f_{xc} = \delta v_{xc} / \delta \rho$ 是交换关联核。

---

## 2. Casida 方程

### 2.1 Casida 本征值问题

线性响应 TDDFT 在频域中等价于求解 **Casida 方程**（Casida, 1995），这是一个矩阵本征值问题：

$$\begin{pmatrix} \mathbf{A} & \mathbf{B} \\ \mathbf{B}^* & \mathbf{A}^* \end{pmatrix} \begin{pmatrix} \mathbf{X}_n \\ \mathbf{Y}_n \end{pmatrix} = \omega_n \begin{pmatrix} \mathbf{1} & \mathbf{0} \\ \mathbf{0} & -\mathbf{1} \end{pmatrix} \begin{pmatrix} \mathbf{X}_n \\ \mathbf{Y}_n \end{pmatrix}$$

其中 $\omega_n$ 是第 $n$ 个激发态的激发能，$(\mathbf{X}_n, \mathbf{Y}_n)$ 是跃迁振幅向量。

### 2.2 A 和 B 矩阵的物理意义

**A 矩阵**（占据→虚轨道跃迁之间的耦合）：

$$A_{ia, jb} = \delta_{ij}\delta_{ab}(\varepsilon_a - \varepsilon_i) + (ia|jb) + (ia|f_{xc}|jb)$$

各项含义：
- $\varepsilon_a - \varepsilon_i$：单粒子激发能（KS 轨道能量差）
- $(ia|jb)$：库仑耦合项（Hartree response）
- $(ia|f_{xc}|jb)$：交换关联耦合项

**B 矩阵**（描述去激发耦合）：

$$B_{ia, jb} = (ia|bj) + (ia|f_{xc}|bj)$$

### 2.3 Tamm-Dancoff 近似 (TDA)

如果忽略 B 矩阵（即令 $\mathbf{B} = 0$），Casida 方程简化为：

$$\mathbf{A} \mathbf{F}_n = \omega_n \mathbf{F}_n$$

这就是 **Tamm-Dancoff 近似 (TDA)**。优势：
- 方程变为 Hermitian 本征值问题，求解更稳定
- 计算量减半
- 对大多数分子，精度损失很小
- 可以避免一些非物理的不稳定态

### 2.4 振子强度

从 Casida 方程的解可以获得**振子强度 (Oscillator Strength)**，它决定了光吸收的强度：

$$f_n = \frac{2}{3}\omega_n \left|\langle 0 | \hat{\boldsymbol{\mu}} | n \rangle\right|^2 = \frac{2\omega_n}{3} \left|\sum_{ia} \left(X_{ia}^n + Y_{ia}^n\right) \langle i | \hat{\mathbf{r}} | a \rangle \right|^2$$

振子强度的物理意义：
- $f_n \propto$ 吸收截面
- 满足 Thomas-Reiche-Kuhn 求和规则：$\sum_n f_n = N_e$（电子总数）
- 只有 $f_n > 0$ 的态才能通过电偶极跃迁被激发

---

## 3. 交换关联核 f_xc

### 3.1 核函数的定义

交换关联核 $f_{xc}$ 是 TDDFT 的核心量，定义为：

$$f_{xc}(\mathbf{r}, \mathbf{r}', \omega) = \frac{\delta v_{xc}[\rho](\mathbf{r}, \omega)}{\delta \rho(\mathbf{r}', \omega)}$$

它描述了密度扰动如何影响交换关联势。

### 3.2 ALDA (Adiabatic LDA)

最简单也是最常用的近似是 **ALDA (Adiabatic LDA)**：

$$f_{xc}^{\text{ALDA}}(\mathbf{r}, \mathbf{r}') = \frac{d^2 e_{xc}^{\text{LDA}}}{d\rho^2} \delta(\mathbf{r} - \mathbf{r}')$$

特点：
- 局部的、不含频的、无记忆的
- 对**局域激发**（如 $\pi \to \pi^*$）表现良好
- 对**Rydberg 态**和**电荷转移态**表现差

### 3.3 杂化泛函

使用杂化泛函（如 B3LYP、PBE0）时，$f_{xc}$ 包含精确交换的贡献：

$$f_{xc}^{\text{hybrid}} = (1-\alpha) f_{xc}^{\text{GGA}} + \alpha f_x^{\text{exact}}$$

杂化泛函通常改善了激发能，但精确交换的比例 $\alpha$ 是经验参数。

### 3.4 范围分离杂化泛函

对于**电荷转移 (Charge Transfer, CT)** 激发，标准泛函（包括杂化泛函）严重低估激发能。原因是 CT 涉及电子从给体到受体的长程移动，而标准泛函不能正确描述长程交换作用。

**范围分离杂化泛函 (Range-Separated Hybrid, RSH)** 通过将库仑核分为短程和长程两部分来解决：

$$\frac{1}{r_{12}} = \frac{\text{erfc}(\omega r_{12})}{r_{12}} + \frac{\text{erf}(\omega r_{12})}{r_{12}}$$

其中 $\omega$ 是范围分离参数。短程部分使用 GGA 交换，长程部分使用精确交换：

$$E_x^{\text{RSH}} = E_x^{\text{SR-GGA}} + E_x^{\text{LR-HF}}$$

常用 RSH 泛函：
- **CAM-B3LYP**：Coulomb Attenuating Method
- **ωB97X-D**：含色散校正
- **LC-ωPBE**：长程校正 PBE

在 CP2K 中可以使用 `&RANGE_SEPARATED` 来实现范围分离杂化。

---

## 4. CP2K 中的 TDDFT

### 4.1 输入结构

CP2K 的 TDDFT 模块在 `&DFT` 部分下通过 `&TDDFT` 关键字启用：

```fortran
&DFT
  ...  # 基态 SCF 设置
  &TDDFT
    NSTATES 10                    # 计算 10 个激发态
    KERNEL ALDA                   # 使用 ALDA 核
    &DIAG
      ALGORITHM STANDARD          # 对角化算法
    &END DIAG
  &END TDDFT
&END DFT
```

### 4.2 关键参数

| 参数 | 说明 | 典型值 |
|------|------|--------|
| `NSTATES` | 激发态数目 | 5–50（根据需要） |
| `KERNEL` | XC 核类型：`ALDA`、`NONE`（HF） | `ALDA` |
| `CONVERGENCE` | 本征值收敛阈值 | 1.0E-6 |
| `NEV` | Davidson 算法的试探向量数 | 略大于 NSTATES |

### 4.3 对角化算法

CP2K 提供两种对角化方案：

1. **`STANDARD`**：直接对角化完整的 A 矩阵。适用于小体系（少量基函数和激发态）。
2. **`DAVIDSON`**：迭代 Davidson 算法，只需要少数几个本征值/本征向量。适用于大体系。

```fortran
&TDDFT
  NSTATES 10
  KERNEL ALDA
  &DIAG
    ALGORITHM DAVIDSON
  &END DIAG
&END TDDFT
```

### 4.4 RS (Reference State) 部分

在某些情况下，可以指定 TDDFT 的参考态：

```fortran
&RS
  # 非占据轨道数（如果少于全部虚轨道，加速计算）
  # 通常不需要手动设置，CP2K 会自动选择
&END RS
```

### 4.5 基组和泛函的选择

- **基组：** DZVP 足够用于定性分析；定量计算推荐 TZVP 或更大基组
- **泛函：** GGA (PBE) + ALDA 核用于初步计算；杂化泛函 (PBE0, B3LYP) 改善精度；RSH 泛函用于 CT 激发

### 4.6 输出信息

TDDFT 计算完成后，CP2K 输出包含：
- 每个激发态的能量（eV 和 nm）
- 振子强度
- 主要跃迁贡献（占据→虚轨道）
- 总激发能

典型输出：

```
TDDFT| Excited State  1:     7.81 eV    158.7 nm    f = 0.3512
TDDFT|   Occupied  -> Virtual    Coefficient
TDDFT|     4 (HOMO)  ->   5 (LUMO)    0.9876
TDDFT|
TDDFT| Excited State  2:     9.24 eV    134.2 nm    f = 0.0000
TDDFT|   Occupied  -> Virtual    Coefficient
TDDFT|     3         ->   5 (LUMO)    0.9921
```

---

## 5. 吸收光谱

### 5.1 从振子强度到吸收截面

光吸收的宏观可观测量是**吸收截面** $\sigma(\omega)$ 或**消光系数** $\varepsilon(\omega)$。在偶极近似下：

$$\sigma(\omega) = \frac{\pi e^2}{2 m_e c \varepsilon_0} \sum_n f_n \, g(\omega - \omega_n)$$

其中 $g(\omega - \omega_n)$ 是展宽函数。

### 5.2 光谱展宽

计算得到的是离散的 $(\omega_n, f_n)$ 数据点，需要展宽才能与实验光谱比较。常用展宽函数：

**Lorentz 展宽：**

$$g_L(\omega - \omega_n) = \frac{\Gamma / (2\pi)}{(\omega - \omega_n)^2 + (\Gamma/2)^2}$$

其中 $\Gamma$ 是半高全宽 (FWHM)，典型值 0.1–0.5 eV。

**Gaussian 展宽：**

$$g_G(\omega - \omega_n) = \frac{1}{\sigma\sqrt{2\pi}} \exp\left(-\frac{(\omega - \omega_n)^2}{2\sigma^2}\right)$$

Lorentz 展宽反映自然线宽和寿命展宽，Gaussian 展宽反映热运动和仪器展宽。实验光谱通常是两者的卷积（Voigt 线型）。

### 5.3 CP2K 输出的跃迁偶极矩

CP2K 的 TDDFT 输出包含了跃迁偶极矩的分量：

$$\langle 0 | \hat{\mu}_\alpha | n \rangle = \sum_{ia} (X_{ia}^n + Y_{ia}^n) \langle i | r_\alpha | a \rangle$$

$\alpha = x, y, z$。振子强度是三个分量的平方和乘以 $2\omega_n/3$。

### 5.4 光谱生成脚本

可以使用如下 Python 脚本从 CP2K TDDFT 输出中生成吸收光谱：

```python
import numpy as np
import matplotlib.pyplot as plt

# 从 CP2K 输出中提取的能量和振子强度
energies = [7.81, 9.24, 10.56]  # eV
osc_strengths = [0.3512, 0.0000, 0.1245]

# Lorentz 展宽
omega = np.linspace(4, 15, 1000)  # eV
gamma = 0.3  # FWHM in eV

spectrum = np.zeros_like(omega)
for e, f in zip(energies, osc_strengths):
    spectrum += f * (gamma/2) / ((omega - e)**2 + (gamma/2)**2)

plt.plot(omega, spectrum)
plt.xlabel('Energy (eV)')
plt.ylabel('Absorption (arb. units)')
plt.title('TDDFT Absorption Spectrum')
plt.savefig('absorption.png', dpi=150)
```

---

## 6. 适用范围

### 6.1 TDDFT 的优势

- **计算效率：** 标度为 $O(N^3)$ 到 $O(N^4)$，远低于波函数方法（CISD, CASPT2, etc.）
- **精度：** 对于局域激发（valence transitions），误差通常在 0.2–0.5 eV
- **可扩展性：** 可以处理数百个原子的体系

### 6.2 适合的体系和激发类型

| 激发类型 | 推荐泛函 | TDDFT 表现 |
|----------|----------|------------|
| 局域 valence 激发 | PBE0, B3LYP | 好 (±0.3 eV) |
| π→π* 跃迁 | ALDA / 杂化 | 很好 |
| n→π* 跃迁 | 杂化泛函 | 好 |
| Rydberg 态 | RSH (CAM-B3LYP) | 中等 |
| 电荷转移 CT | RSH 必须 | ALDA 差，RSH 好 |
| 双电子激发 | 所有 | 定性失败 |
| 暗态（振子强度为零） | — | 能量可用，强度为零 |

### 6.3 不适合的情况

**半导体和绝缘体的光学性质：**

TDDFT 在分子体系中表现良好，但对于周期性半导体/绝缘体，**不应使用 TDDFT**。原因：

1. 固态中的激子效应（electron-hole interaction）不能被 ALDA 描述
2. GW/BSE 方法是计算固体光学性质的标准方法（见第16章）
3. CP2K 的 TDDFT 模块目前只支持分子（非周期）体系

**多参考态问题：**

如果基态本身就是多参考态特征（如双自由基、过渡金属配合物），单参考 TDDFT 可能给出错误结果。应考虑 CASPT2 或 MRCI 方法。

### 6.4 精度预期

| 泛函 | 典型误差 (eV) | 适用范围 |
|------|---------------|----------|
| PBE/ALDA | 0.5–1.0 | 定性分析 |
| PBE0 | 0.2–0.5 | 定量 |
| B3LYP | 0.2–0.5 | 定量（有机分子） |
| CAM-B3LYP | 0.2–0.4 | CT 态，Rydberg |
| ωB97X-D | 0.1–0.3 | 最佳精度 |

---

## 7. 自旋翻转 TDDFT

### 7.1 基本思想

**自旋翻转 TDDFT (Spin-Flip TDDFT, SF-TDDFT)** 是一种计算开壳层激发态（特别是三重态）的方法：

- 参考态：高自旋（如三重态或开壳层单重态）
- 激发：α 电子翻转为 β 电子（$\Delta S = -1$）
- 结果：从高自旋参考态出发，可以同时获得不同自旋多重度的态

### 7.2 三重态的获取

标准 TDDFT 从闭壳层单重态参考态出发，不能直接给出三重态。而 SF-TDDFT：

1. 以三重态为参考态
2. 自旋翻转产生 $\Delta M_S = -1$ 的态
3. 结果包括：单重态（从三重态翻转得到）和三重态（保持）

这使得可以用一次计算同时获得单重态和三重态激发能。

### 7.3 在 CP2K 中的实现

CP2K 目前的 TDDFT 模块对自旋翻转的支持有限。如果需要 SF-TDDFT，可能需要使用其他程序（如 Q-Chem、ORCA）。

### 7.4 应用场景

- 有机发光材料中的磷光发射（需要单重态-三重态能隙）
- 光催化反应中的系间窜越
- 双自由基体系

---

## 8. 输入文件详解

本章示例文件 `ethylene_tddft.inp` 对乙烯分子进行 TDDFT 吸收光谱计算。以下逐节解释：

### 8.1 &GLOBAL 部分

```fortran
&GLOBAL
  PROJECT ethylene_tddft
  RUN_TYPE ENERGY
&END GLOBAL
```

TDDFT 不需要特殊的 `RUN_TYPE`，使用 `ENERGY` 即可。TDDFT 模块在基态 SCF 收敛后自动执行。

### 8.2 基态 SCF 设置

```fortran
&DFT
  BASIS_SET_FILE_NAME  BASIS_SET
  POTENTIAL_FILE_NAME  GTH_POTENTIALS
  ...
  &XC
    &XC_FUNCTIONAL PBE       # GGA 泛函
    &END XC_FUNCTIONAL
  &END XC
&END DFT
```

- 使用 PBE 泛函计算基态
- 基态 SCF 必须完全收敛，否则 TDDFT 结果不可靠
- `EPS_SCF` 设置较紧的收敛阈值（1.0E-8）

### 8.3 &TDDFT 部分

```fortran
&TDDFT
  NSTATES 10                   # 计算 10 个激发态
  KERNEL ALDA                  # ALDA 交换关联核
  CONVERGENCE 1.0E-6           # 激发能收敛阈值
  &DIAG
    ALGORITHM STANDARD          # 直接对角化（小体系）
  &END DIAG
&END TDDFT
```

- `NSTATES 10`：计算前 10 个激发态。对于乙烯，这足以覆盖低能吸收带
- `KERNEL ALDA`：使用绝热 LDA 核。因为基态用 PBE（GGA），所以核函数为 GGA 级别的 ALDA
- `CONVERGENCE`：激发能的收敛阈值，单位 Hartree
- `ALGORITHM STANDARD`：直接对角化。对于乙烯这样的小分子，标准算法即可

### 8.4 分子结构

```fortran
&COORD
  C   0.000   0.000   0.6695
  C   0.000   0.000  -0.6695
  H   0.000   0.928   1.2320
  H   0.000  -0.928   1.2320
  H   0.000   0.928  -1.2320
  H   0.000  -0.928  -1.2320
&END COORD
```

乙烯分子（C₂H₄）的平衡构型，z 轴方向为 C=C 键轴。

### 8.5 基组和赝势

```fortran
&KIND C
  BASIS_SET DZVP-GTH-PBE
  POTENTIAL GTH-PBE-q4
&END KIND

&KIND H
  BASIS_SET DZVP-GTH-PBE
  POTENTIAL GTH-PBE-q1
&END KIND
```

- DZVP-GTH-PBE：双ζ+极化基组，与 PBE 赝势配套
- C 使用 GTH-PBE-q4（4 个价电子）
- H 使用 GTH-PBE-q1（1 个价电子）

### 8.6 运行与结果分析

```bash
# 运行
mpirun -np 4 cp2k.psmp -i ethylene_tddft.inp -o ethylene_tddft.out

# 提取激发态信息
grep -A 5 "Excited State" ethylene_tddft.out

# 查看振子强度
grep "f =" ethylene_tddft.out
```

**预期结果：**

乙烯的最低激发态是 π→π* 跃迁：
- 垂直激发能：约 7.8 eV（实验值 ~7.6 eV）
- 强度：振子强度大（允许跃迁）
- 对称性：B₁ᵤ

使用 PBE/ALDA，典型误差在 0.5 eV 左右。使用杂化泛函（PBE0）可改善至 0.2 eV。

### 8.7 进阶建议

1. **增加基组：** 使用 TZVP 或 cc-pVTZ 改善精度
2. **杂化泛函：** 将 PBE 替换为 PBE0 或 B3LYP
3. **更多激发态：** 增加 `NSTATES` 以获得更宽的光谱范围
4. **溶剂效应：** 使用 COSMO 隐式溶剂模型（CP2K 支持）

---

## 参考文献

1. M. E. Casida, *Time-dependent density functional response theory for molecules*, in *Recent Advances in Density Functional Methods*, Part I, p. 155 (1995).
2. E. Runge and E. K. U. Gross, *Density-functional theory for time-dependent systems*, Phys. Rev. Lett. **52**, 997 (1984).
3. A. Dreuw and M. Head-Gordon, *Single-reference ab initio methods for the calculation of excited states of large molecules*, Chem. Rev. **105**, 4009 (2005).
4. R. Bauernschmitt and R. Ahlrichs, *Treatment of electronic excitations within the adiabatic approximation of time dependent density functional theory*, Chem. Phys. Lett. **256**, 454 (1996).
5. T. Yanai, D. P. Tew, and N. C. Handy, *A new hybrid exchange-correlation functional using the Coulomb-attenuating method (CAM-B3LYP)*, Chem. Phys. Lett. **393**, 51 (2004).
6. The CP2K manual: https://manual.cp2k.org/

---

## 练习

1. **基础练习：** 运行 `ethylene_tddft.inp`，记录前 5 个激发态的能量和振子强度。

2. **泛函比较：** 将 PBE 替换为 PBE0（需要在 `&XC_FUNCTIONAL` 中使用 `&PBE0`），比较激发能的变化。

3. **基组效应：** 分别使用 SZV 和 TZVP 基组计算，观察基组大小对激发能的影响。

4. **光谱绘制：** 使用 Python 脚本将 TDDFT 结果绘制成吸收光谱，设置不同的展宽参数。

5. **电荷转移态：** 构建一个简单的给体-受体分子（如甲醛 + 氨的复合物），比较 ALDA 和 CAM-B3LYP 对 CT 激发的描述差异。
