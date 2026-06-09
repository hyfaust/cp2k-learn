# 第15章：后 Hartree-Fock 方法 —— RI-MP2

> **难度：** ⭐⭐⭐⭐⭐  
> **关键词：** Hartree-Fock、MP2、RI近似、电子关联、相关能  
> **RUN_TYPE：** `ENERGY`（MP2 在 `&WF_CORRELATION` 部分设置）

## 目录

1. [Hartree-Fock 方法回顾](#1-hartree-fock-方法回顾)
2. [Møller-Plesset 微扰理论](#2-møller-plesset-微扰理论)
3. [RI 近似 (Resolution of Identity)](#3-ri-近似-resolution-of-identity)
4. [CP2K 的 MP2 实现](#4-cp2k-的-mp2-实现)
5. [基组要求](#5-基组要求)
6. [双杂化泛函](#6-双杂化泛函)
7. [RI-MP2 梯度](#7-ri-mp2-梯度)
8. [BSSE 校正](#8-bsse-校正)
9. [输入文件详解](#9-输入文件详解)

---

## 1. Hartree-Fock 方法回顾

### 1.1 平均场近似

Hartree-Fock (HF) 方法是量子化学的基础。其核心思想是**平均场近似 (Mean-Field Approximation)**：每个电子在其他所有电子产生的平均场中运动。

HF 能量可以写为：

$$E_{\text{HF}} = \sum_i h_{ii} + \frac{1}{2}\sum_{ij} (J_{ij} - K_{ij})$$

其中：
- $h_{ii}$：单电子积分（动能 + 核-电子吸引）
- $J_{ij} = (ii|jj)$：库仑积分（经典电子-电子排斥）
- $K_{ij} = (ij|ji)$：交换积分（量子力学效应，仅存在于自旋平行电子间）

### 1.2 HF 的精确交换

HF 方法的一个重要优势是它包含了**精确交换 (Exact Exchange)**。DFT 中 LDA/GGA 的交换泛函是近似的，而 HF 交换是精确的——它是基于 Slater 行列式波函数计算得到的。

### 1.3 关联能的缺失

HF 方法的主要缺陷是完全忽略了**电子关联 (Electron Correlation)**。关联能定义为：

$$E_{\text{corr}} = E_{\text{exact}} - E_{\text{HF}}$$

关联能来源于电子的瞬时相互作用（超越平均场的部分）。对于化学键的断裂、范德华相互作用等，电子关联至关重要。

典型关联能大小：
- He 原子：~0.042 hartree (~1.1 eV)
- H₂O 分子：~0.3 hartree (~8 eV)
- 每个电子对约 0.04 hartree

### 1.4 HF 轨道和轨道能量

HF 迭代求解得到的轨道 $\{\phi_i\}$ 和轨道能量 $\{\varepsilon_i\}$ 满足：

$$\hat{f} \phi_i = \varepsilon_i \phi_i$$

其中 $\hat{f}$ 是 Fock 算符。这些轨道和轨道能量是后 HF 方法的出发点。特别地：

- **Koopmans 定理：** $-\varepsilon_{\text{HOMO}} \approx I$（电离能），$-\varepsilon_{\text{LUMO}} \approx A$（电子亲和能）
- 轨道能量差 $\varepsilon_a - \varepsilon_i$ 构成 MP2 公式中的分母

---

## 2. Møller-Plesset 微扰理论

### 2.1 基本框架

**Møller-Plesset (MP) 微扰理论**是基于 HF 参考态的 Rayleigh-Schrödinger 微扰理论。将 Hamiltonian 分解为：

$$\hat{H} = \hat{H}_0 + \hat{V}$$

其中 $\hat{H}_0 = \sum_i \hat{f}_i$ 是 Fock 算符之和（零级 Hamiltonian），$\hat{V} = \hat{H} - \hat{H}_0$ 是微扰。

- 零级能量：$E^{(0)} = \sum_i \varepsilon_i$
- 一级校正：$E^{(1)} = -\frac{1}{2}\sum_{ij}(J_{ij} - K_{ij})$
- HF 能量：$E_{\text{HF}} = E^{(0)} + E^{(1)}$

因此，**关联能从二级校正 $E^{(2)}$ 开始**。

### 2.2 MP2 能量公式

MP2 是最简单的后 HF 方法，其能量校正为：

$$E^{(2)}_{\text{MP2}} = -\sum_{i<j}\sum_{a<b} \frac{|(ia|jb) - (ib|ja)|^2}{\varepsilon_a + \varepsilon_b - \varepsilon_i - \varepsilon_j}$$

等价地（对于闭壳层体系）：

$$E^{(2)}_{\text{MP2}} = -\sum_{ia}\sum_{jb} \frac{(ia|jb)[2(ia|jb) - (ib|ja)]}{\varepsilon_a + \varepsilon_b - \varepsilon_i - \varepsilon_j}$$

其中：
- $i, j$：占据轨道指标
- $a, b$：虚轨道指标
- $(ia|jb)$：双电子积分 $\int \phi_i^*(\mathbf{r}_1) \phi_a(\mathbf{r}_1) \frac{1}{r_{12}} \phi_j^*(\mathbf{r}_2) \phi_b(\mathbf{r}_2) d\mathbf{r}_1 d\mathbf{r}_2$
- $\varepsilon_a + \varepsilon_b - \varepsilon_i - \varepsilon_j$：能量分母（始终为正）

### 2.3 物理理解

MP2 能量的物理图像清晰：

- **分子：** 同时激发两个电子从占据轨道 $i, j$ 到虚轨道 $a, b$
- **贡献：** 所有双激发组态对关联能的贡献（一阶波函数系数）
- **分母：** 能量守恒——激发能越高，贡献越小
- **分子：** 电子间库仑和交换相互作用

### 2.4 MP2 的标度

MP2 的计算量主要在于四指标双电子积分 $(ia|jb)$ 的求和：

- **计算标度：** $O(N^5)$（$N$ 为基函数数目）
- **存储标度：** $O(N^4)$（存储所有双电子积分）

这比 DFT ($O(N^3)$) 昂贵得多，但比 CCSD ($O(N^6)$) 或 CCSD(T) ($O(N^7)$) 便宜。

### 2.5 MP2 的优缺点

**优点：**
- 系统性改进 HF 结果
- 无经验参数
- 大小一致性 (Size-Consistent)
- 为更高级别方法（CCSD, CCSD(T)）提供参考

**缺点：**
- 对于强关联体系（多参考态特征）可能发散或给出错误结果
- $O(N^5)$ 标度限制了可处理的体系大小
- 色散相互作用描述虽好但需要大基组
- 基组收敛慢（需基组外推）

---

## 3. RI 近似 (Resolution of Identity)

### 3.1 动机

MP2 的计算瓶颈是四指标双电子积分：

$$(ia|jb) = \int \phi_i(\mathbf{r}_1) \phi_a(\mathbf{r}_1) \frac{1}{r_{12}} \phi_j(\mathbf{r}_2) \phi_b(\mathbf{r}_2) d\mathbf{r}_1 d\mathbf{r}_2$$

对于 $N$ 个基函数，这样的积分有 $O(N^4)$ 个。**Resolution of Identity (RI)** 近似（也称为密度拟合 Density Fitting）是降低计算量的关键技术。

### 3.2 RI 的基本思想

RI 的核心思想是引入一组辅助基函数 $\{P\}$（auxiliary basis set），将两中心乘积 $\phi_i(\mathbf{r})\phi_a(\mathbf{r})$ 展开：

$$\phi_i(\mathbf{r}_1)\phi_a(\mathbf{r}_1) \approx \sum_P c_{ia}^P \chi_P(\mathbf{r}_1)$$

其中展开系数为：

$$c_{ia}^P = \sum_Q (ia|Q) (P|Q)^{-1}$$

代入四指标积分：

$$(ia|jb) \approx \sum_{PQ} (ia|P) (P|Q)^{-1} (Q|jb)$$

### 3.3 RI 的计算优势

| 量 | 原始 | RI 近似后 |
|----|------|-----------|
| 积分数量 | $O(N^4)$ | $O(N^3)$（三指标积分） |
| RI-MP2 标度 | $O(N^5)$ | 保持 $O(N^5)$，但前因子大幅减小 |
| 存储 | $O(N^4)$ | $O(N^3)$ |

RI-MP2 的总标度仍为 $O(N^5)$，因为积分变换从原子轨道 (AO) 到分子轨道 (MO) 的标度为 $O(N^5)$。但前因子减小了约 $10-100$ 倍，使得实际计算快了很多。

### 3.4 RI 基组

为了使用 RI 近似，需要专用的辅助基组。在 CP2K 中：

- 主基组：`BASIS_SET cc-DZ` 或 `BASIS_SET cc-TZ`
- RI 基组：`BASIS_RI_SET RI_cc-DZ` 或 `BASIS_RI_SET RI_cc-TZ`

RI 基组通常是为特定主基组优化的。CP2K 提供了预定义的 RI 基组文件 `BASIS_RI`。

### 3.5 RI 误差

RI 近似引入的误差通常很小：
- 对于总能量：误差 < 0.1 kcal/mol
- 对于能量差：误差可以忽略
- RI 误差通常远小于基组不完备误差

因此，RI 近似在实际计算中几乎是"免费的"——精度损失可忽略，但计算效率大幅提升。

---

## 4. CP2K 的 MP2 实现

### 4.1 三种实现方式

CP2K 提供了三种 MP2 计算方案：

1. **Canonical MP2：** 传统 MP2，直接计算四指标积分。适合小分子和小基组，但标度为 $O(N^5)$。
2. **GPW-based MP2：** 利用 CP2K 的 Gaussian-Plane Wave 框架。通过辅助平面波展开库仑核。
3. **RI-MP2 (RI_MP2_GPW)：** 最常用的方案。结合 RI 近似和 GPW 方法，是 CP2K 中最经济的 MP2 实现。

### 4.2 &WF_CORRELATION 和 &MP2 部分

MP2 计算通过 `&WF_CORRELATION` 部分触发：

```fortran
&DFT
  ...
  &WF_CORRELATION
    METHOD RI_MP2_GPW           # RI-MP2 with GPW
    &MP2
      METHOD RI_MP2_GPW
      RI_METRIC RI              # 使用 RI 库仑度规
      &RI
        RI_SIGMA 1.0E-6         # RI 收敛参数
      &END RI
    &END MP2
    &INTERACTION_POTENTIAL
      POTENTIAL_TYPE TRUNCATED  # 截断库仑势（周期体系）
      CUTOFF_RADIUS 6.0         # 截断半径 (Angstrom)
    &END INTERACTION_POTENTIAL
    &MEMORY
      MAX_MEMORY 2000           # 每个 MPI 进程的最大内存 (MB)
    &END MEMORY
    EPS_SCF 1.0E-6              # SCF 收敛阈值
  &END WF_CORRELATION
&END DFT
```

### 4.3 关键参数详解

| 参数 | 说明 | 建议值 |
|------|------|--------|
| `METHOD` | MP2 方法 | `RI_MP2_GPW`（推荐） |
| `EPS_SCF` | 求解 CP-HF 方程的收敛阈值 | 1.0E-6 |
| `MAX_MEMORY` | 内存限制 (MB/进程) | 根据可用内存设置 |
| `CUTOFF_RADIUS` | 截断库仑势半径 | 根据盒子大小设置 |
| `NUM_INTEG_GROUPS` | 积分组数（并行） | 根据进程数优化 |

### 4.4 前置 HF 计算

RI-MP2 需要先获得 HF 轨道。在 CP2K 中，这通过设置 `&XC` 部分为 HF 交换来实现：

```fortran
&XC
  &XC_FUNCTIONAL NONE
  &END XC_FUNCTIONAL
  &HF
    FRACTION 1.0               # 100% 精确交换 = HF
    &SCREENING
      EPS_SCHWARZ 1.0E-7       # Schwarz 筛选阈值
    &END SCREENING
    &INTERACTION_POTENTIAL
      POTENTIAL_TYPE TRUNCATED  # 截断库仑势
      CUTOFF_RADIUS 6.0
    &END INTERACTION_POTENTIAL
  &END HF
&END XC
```

### 4.5 计算流程

1. SCF 迭代获得 HF 轨道（`&XC_FUNCTIONAL NONE` + `&HF`）
2. 轨道从 AO 基变换到 MO 基
3. 三指标积分计算 $(ia|P)$
4. 求解 CP-HF 方程（获得弛豫轨道）
5. 计算 MP2 相关能

---

## 5. 基组要求

### 5.1 后 HF 方法对基组的敏感性

与 DFT 不同，后 HF 方法的精度**强烈依赖基组大小**。这是因为：

1. **电子关联发生在电子之间距离很近时**：需要描述电子云重叠区域的精细结构
2. **基组完备性极限 (CBS)：** 只有在完备基组极限下，MP2 能量才是"正确的" MP2 能量
3. **基组收敛慢：** 关联能的基组收敛行为为 $E \propto X^{-3}$（$X$ 为基组级别）

### 5.2 推荐基组

| 精度级别 | 基组 | 典型应用 |
|----------|------|----------|
| 定性 | cc-DZ | 初步筛选、大体系 |
| 定量 | cc-TZ | 标准计算 |
| 高精度 | cc-QZ | 基组外推、小分子 |
| 带弥散 | aug-cc-DZ | 阴离子、弱相互作用 |

### 5.3 基组外推

为获得 CBS 极限的关联能，使用两点外推公式：

$$E_{\text{corr}}(X) = E_{\text{corr}}^{\text{CBS}} + \frac{A}{X^3}$$

其中 $X$ 是基组的 zeta 数（2, 3, 4, ...）。使用两个不同大小的基组（如 cc-TZ 和 cc-QZ），可以拟合得到 $E_{\text{corr}}^{\text{CBS}}$。

HF 能量的基组收敛更快（指数级），通常 cc-TZ 已足够。

### 5.4 弥散函数

对于以下情况，**必须**使用弥散基函数（aug-前缀）：
- 阴离子
- 弱相互作用（氢键、范德华力）
- 激发态
- 极化率计算

弥散函数描述了电子密度在分子远端的分布，对这些性质至关重要。

### 5.5 CP2K 中的基组设置

CP2K 使用 GTH 赝势基组，但也可以使用全电子基组。对于 MP2，推荐：

```fortran
&BASIS_SET cc-DZ                  # 或 cc-TZ
&BASIS_RI_SET RI_cc-DZ            # RI 辅助基组，与主基组配套
```

如果使用赝势基组，需要配套的 RI 基组。CP2K 的 `BASIS_RI` 文件中包含了预优化的 RI 基组。

---

## 6. 双杂化泛函

### 6.1 DFT + MP2 = 双杂化泛函

**双杂化泛函 (Double-Hybrid Functional)** 将 DFT 和 MP2 结合：

$$E_{xc}^{\text{DH}} = a_x E_x^{\text{exact}} + (1-a_x) E_x^{\text{DFT}} + a_c E_c^{\text{MP2}} + (1-a_c) E_c^{\text{DFT}}$$

其中 $a_x$ 和 $a_c$ 是混合系数。

### 6.2 典型双杂化泛函

| 泛函 | $a_x$ | $a_c$ | 特点 |
|------|-------|-------|------|
| B2PLYP | 0.53 | 0.27 | 最早的双杂化泛函 |
| B2PLYP-D3 | 0.53 | 0.27 | 含色散校正 |
| PWPB95 | 0.50 | 0.22 | mGGA + MP2 |
| DSD-PBEP86 | 0.69 | 0.22 | 色散校正，高精度 |

### 6.3 CP2K 中的实现

在 CP2K 中，双杂化泛函通过调节 HF 交换比例和 MP2 关联比例实现：

```fortran
&XC
  &XC_FUNCTIONAL
    &LYP
      SCALE_C 0.73
    &END LYP
    &BECKE88
      SCALE_X 0.47
    &END BECKE88
  &END XC_FUNCTIONAL
  &HF
    FRACTION 0.53
    &SCREENING
      EPS_SCHWARZ 1.0E-7
    &END SCREENING
  &END HF
&END XC

&WF_CORRELATION
  &MP2
    SCALE_S 0.0                 # 闭壳层单重态 MP2 标度
    SCALE_T 0.27                # 三重态 MP2 标度（对应 a_c）
  &END MP2
&END WF_CORRELATION
```

注意：`SCALE_S` 和 `SCALE_T` 分别调节闭壳层自旋对和三重态自旋对的 MP2 贡献。对于 B2PLYP，$a_c = 0.27$。

### 6.4 双杂化泛函的优势

- **精度高：** 通常比纯 DFT 和纯 MP2 都好
- **误差抵消：** DFT 和 MP2 的误差有时可以抵消
- **基组要求低：** 由于 DFT 部分对基组不敏感，可以使用较小基组
- **适合热化学：** 反应能、活化能等

### 6.5 局限性

- 计算量介于 DFT 和 MP2 之间（仍需 $O(N^5)$ 的 MP2 步骤）
- 经验参数多（$a_x$, $a_c$ 拟合得到）
- 对强关联体系无改善

---

## 7. RI-MP2 梯度

### 7.1 解析梯度的必要性

如果需要在 MP2 级别做几何优化或分子动力学，就需要 MP2 能量对原子坐标的梯度（解析梯度）。

**数值梯度**（有限差分）虽然可行，但对于 $3N$ 个坐标需要 $2 \times 3N$ 次能量计算，计算量巨大。**解析梯度**只需一次计算。

### 7.2 Z-向量方法

MP2 梯度的核心是 **CPHF (Coupled Perturbed Hartree-Fock)** 方程，即求解轨道对核坐标的响应。通过 **Z-向量方法 (Z-vector Method)**：

$$E^{(1)} = \sum_\mu \frac{\partial E}{\partial \mathbf{R}_\mu} = \sum_\mu \left(\frac{\partial E}{\partial \mathbf{R}_\mu}\right)_{\text{Hellmann-Feynman}} + \sum_{pq} Z_{pq} \frac{\partial F_{pq}}{\partial \mathbf{R}_\mu}$$

其中 $Z_{pq}$ 是通过求解 CPHF 方程得到的辅助量。

### 7.3 CP2K 中的梯度计算

RI-MP2 梯度在 CP2K 中可以通过设置 `&RESPONSE` 部分来计算：

```fortran
&WF_CORRELATION
  METHOD RI_MP2_GPW
  &MP2
    METHOD RI_MP2_GPW
  &END MP2
  &RESPONSE
    MAX_ITER 50                  # CPHF 迭代最大次数
    EPS 1.0E-6                   # 收敛阈值
  &END RESPONSE
&END WF_CORRELATION
```

### 7.4 应用

RI-MP2 梯度可以用于：
- MP2 级别的几何优化（`RUN_TYPE GEO_OPT`）
- MP2 级别的频率计算（`RUN_TYPE VIBRATIONAL_ANALYSIS`）
- MP2 级别的分子动力学（`RUN_TYPE MD`）

这些计算比 DFT 贵得多，但对于需要高精度的体系非常有价值。

---

## 8. BSSE 校正

### 8.1 什么是 BSSE？

**基组重叠误差 (Basis Set Superposition Error, BSSE)** 是有限基组计算中一个重要的系统误差。当两个分子接近时，每个分子可以"借用"对方的基函数来改善自身的描述，人为地降低了相互作用能。

$$E_{\text{int}}^{\text{BSSE}} = E_{AB}^{AB} - E_A^{AB} - E_B^{AB}$$

其中上标表示使用的基组，下标表示分子。$E_A^{AB}$ 是在 AB 复合物的基组中计算的单体 A 的能量。

### 8.2 Counterpoise 校正

**Counterpoise (CP) 校正**是消除 BSSE 的标准方法：

$$\Delta E_{\text{int}}^{\text{CP}} = E_{AB}^{AB}(\mathbf{R}_{AB}) - E_A^{AB}(\mathbf{R}_{AB}) - E_B^{AB}(\mathbf{R}_{AB})$$

具体步骤：
1. 在**复合物几何构型**下，用**完整基组**计算复合物能量 $E_{AB}^{AB}$
2. 在**相同几何构型**下，用**完整基组**（含 ghost 原子基函数）计算单体 A 的能量 $E_A^{AB}$
3. 同理计算 $E_B^{AB}$
4. 相互作用能 = $E_{AB}^{AB} - E_A^{AB} - E_B^{AB}$

### 8.3 Ghost 原子

Ghost 原子是只提供基函数但不贡献核-电子势和核-核排斥的原子。在 CP2K 中，可以通过设置 `CHARGE 0` 和去除赝势来近似实现。

### 8.4 CP2K 中的 BSSE 校正

CP2K 没有内置的 counterpoise 校正功能，需要用户手动进行三个计算：

```fortran
! 计算 1: 复合物 AB（完整系统）
&FORCE_EVAL
  ...
  &SUBSYS
    &KIND A ... &END
    &KIND B ... &END
  &END SUBSYS
&END FORCE_EVAL

! 计算 2: 单体 A + ghost B（手动将 B 的赝势设为零，保留基函数）
! 计算 3: 单体 B + ghost A
```

### 8.5 BSSE 的大小

BSSE 的大小取决于：
- **基组大小：** 越小的基组 BSSE 越大。cc-DZ 时 BSSE 可能占相互作用能的 10-30%
- **体系类型：** 氢键、范德华复合物受 BSSE 影响最大
- **计算方法：** MP2 的 BSSE 通常比 DFT 大

经验法则：使用 cc-TZ 或更大的基组时，BSSE 通常可以忽略。

---

## 9. 输入文件详解

本章示例文件 `water_mp2.inp` 对水二聚体进行 RI-MP2 能量计算。以下逐节解释：

### 9.1 整体结构

```fortran
&GLOBAL
  PROJECT water_mp2
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    ...HF 设置...
    &WF_CORRELATION
      ...MP2 设置...
    &END WF_CORRELATION
  &END DFT
  &SUBSYS
    ...结构和基组...
  &END SUBSYS
&END FORCE_EVAL
```

关键点：MP2 计算嵌套在 `&DFT` 部分的 `&WF_CORRELATION` 中。

### 9.2 HF 前置计算

```fortran
&XC
  &XC_FUNCTIONAL NONE           # 不使用 DFT XC 泛函
  &END XC_FUNCTIONAL
  &HF
    FRACTION 1.0                 # 100% 精确交换 = Hartree-Fock
    &SCREENING
      EPS_SCHWARZ 1.0E-7         # Schwarz 积分筛选
    &END SCREENING
    &INTERACTION_POTENTIAL
      POTENTIAL_TYPE TRUNCATED   # 截断库仑势
      CUTOFF_RADIUS 6.0          # 截断半径
    &END INTERACTION_POTENTIAL
  &END HF
&END XC
```

- `XC_FUNCTIONAL NONE`：关闭 DFT 交换关联
- `HF FRACTION 1.0`：纯 HF 交换
- `POTENTIAL_TYPE TRUNCATED`：对于周期体系，截断库仑势避免镜像相互作用

### 9.3 RI-MP2 设置

```fortran
&WF_CORRELATION
  METHOD RI_MP2_GPW
  &MP2
    METHOD RI_MP2_GPW
    RI_METRIC RI
  &END MP2
  &INTERACTION_POTENTIAL
    POTENTIAL_TYPE TRUNCATED
    CUTOFF_RADIUS 6.0
  &END INTERACTION_POTENTIAL
  &MEMORY
    MAX_MEMORY 2000
  &END MEMORY
&END WF_CORRELATION
```

- `METHOD RI_MP2_GPW`：使用 RI 近似的 MP2，GPW 积分方案
- `MAX_MEMORY 2000`：每个进程 2 GB 内存限制
- `CUTOFF_RADIUS 6.0`：截断半径（应小于盒子半径）

### 9.4 基组设置

```fortran
&KIND O
  BASIS_SET cc-DZ               # 主基组：相关一致双ζ
  BASIS_RI_SET RI_cc-DZ         # RI 辅助基组
  POTENTIAL GTH-HF-q6           # HF 赝势
&END KIND

&KIND H
  BASIS_SET cc-DZ
  BASIS_RI_SET RI_cc-DZ
  POTENTIAL GTH-HF-q1
&END KIND
```

- `cc-DZ`：Dunning 相关一致双ζ基组（配合 GTH 赝势）
- `RI_cc-DZ`：配套的 RI 辅助基组
- `GTH-HF-q6`：GTH 赝势，HF 参数化

### 9.5 盒子大小

```fortran
&CELL
  ABC 12.0 12.0 12.0
&END CELL
```

MP2 计算需要比 DFT 更大的盒子，因为关联效应涉及更远的电子-电子相互作用。对于水二聚体，12 Å 应该足够。但需要测试盒子大小的收敛性。

### 9.6 运行与分析

```bash
# 运行（需要较多内存和时间）
mpirun -np 8 cp2k.psmp -i water_mp2.inp -o water_mp2.out

# 查看 HF 能量
grep "ENERGY| Total" water_mp2.out | head -1

# 查看 MP2 相关能
grep "MP2" water_mp2.out

# 总能量 = HF 能量 + MP2 相关能
```

### 9.7 预期结果

对于水二聚体：
- HF 相互作用能：~5-8 kcal/mol（包含 BSSE）
- MP2 相互作用能：~4-6 kcal/mol（包含 BSSE 和色散贡献）
- CCSD(T)/CBS 参考值：~5.0 kcal/mol

MP2 相比 HF 增加了色散相互作用的描述，这对水二聚体很重要。

---

## 参考文献

1. C. Møller and M. S. Plesset, *Note on an approximation treatment for many-electron systems*, Phys. Rev. **46**, 618 (1934).
2. M. Feyereisen, G. Fitzgerald, and A. Komornicki, *Use of approximate integrals in ab initio theory. An application in MP2 energy calculations*, Chem. Phys. Lett. **208**, 359 (1993).
3. S. Grimme, *Semiempirical hybrid density functional with perturbative second-order correlation*, J. Chem. Phys. **124**, 034108 (2006).
4. M. Del Ben, J. Hutter, and J. VandeVondele, *Second-order Møller-Plesset perturbation theory in the condensed phase: An efficient and massively parallel Gaussian and plane waves approach*, J. Chem. Theory Comput. **8**, 4177 (2012).
5. S. Boys and F. Bernardi, *The calculation of small molecular interactions by the differences of separate total energies*, Mol. Phys. **19**, 553 (1970).
6. The CP2K manual: https://manual.cp2k.org/

---

## 练习

1. **基础练习：** 运行 `water_mp2.inp`，分别提取 HF 能量和 MP2 相关能。

2. **基组收敛：** 将基组从 cc-DZ 改为 cc-TZ，观察 MP2 相关能的变化。使用外推公式估算 CBS 极限。

3. **BSSE 估计：** 手动进行 counterpoise 计算，估计水二聚体的 BSSE 大小。

4. **泛函比较：** 将 HF 交换比例改为 0.53，MP2 标度改为 0.27，模拟 B2PLYP 双杂化泛函。

5. **盒子大小收敛：** 将盒子从 10 Å 增加到 15 Å，检查总能量的变化。
