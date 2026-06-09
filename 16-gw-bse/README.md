# 第16章：GW/BSE 能带结构计算

> **难度：** ⭐⭐⭐⭐⭐  
> **关键词：** GW 近似、BSE、准粒子能、带隙、激子效应  
> **RUN_TYPE：** `ENERGY`（GW/BSE 在 `&WF_CORRELATION` / `&RI_RPA` 部分设置）

## 目录

1. [DFT 带隙问题](#1-dft-带隙问题)
2. [GW 近似](#2-gw-近似)
3. [自洽级别](#3-自洽级别)
4. [BSE (Bethe-Salpeter 方程)](#4-bse-bethe-salpeter-方程)
5. [CP2K 的 GW 实现](#5-cp2k-的-gw-实现)
6. [关键参数](#6-关键参数)
7. [基组要求](#7-基组要求)
8. [能带结构计算](#8-能带结构计算)
9. [输入文件详解](#9-输入文件详解)

---

## 1. DFT 带隙问题

### 1.1 Kohn-Sham 本征值的物理意义

KS-DFT 中的本征值方程：

$$\hat{H}_{\text{KS}} \phi_i(\mathbf{r}) = \varepsilon_i \phi_i(\mathbf{r})$$

本征值 $\varepsilon_i$ 是 Lagrange 乘子，它们确保轨道正交归一化，但**并非**系统的真实单粒子激发能。仅在精确 XC 泛函下，HOMO 能量才等于（负的）第一电离能：

$$-\varepsilon_{\text{HOMO}} = I \quad \text{(仅精确泛函)}$$

而带隙：

$$E_{\text{gap}}^{\text{exact}} = I - A = (\varepsilon_{\text{LUMO}} - \varepsilon_{\text{HOMO}}) + \Delta_{xc}$$

其中 **导数不连续性 (Derivative Discontinuity)**：

$$\Delta_{xc} = \left.\frac{\delta E_{xc}[\rho]}{\delta \rho(\mathbf{r})}\right|_{\rho+\delta} - \left.\frac{\delta E_{xc}[\rho]}{\delta \rho(\mathbf{r})}\right|_{\rho}$$

### 1.2 LDA/GGA 的系统性低估

LDA 和 GGA 泛函在描述半导体和绝缘体的带隙时，存在严重的系统性低估：

| 材料 | 实验带隙 (eV) | PBE 带隙 (eV) | 低估幅度 |
|------|---------------|---------------|----------|
| Si | 1.17 | 0.61 | ~48% |
| GaAs | 1.52 | 0.53 | ~65% |
| LiF | 14.2 | 8.9 | ~37% |
| 金刚石 | 5.48 | 3.92 | ~28% |

低估的原因有两个方面：
1. **自相互作用误差 (Self-Interaction Error, SIE)：** LDA/GGA 中电子与自身的库仑排斥未完全消除
2. **导数不连续性缺失：** LDA/GGA 的 XC 势是连续函数，$\Delta_{xc} = 0$

### 1.3 为什么带隙重要？

带隙决定了材料的电子和光学性质：
- **半导体器件：** 带隙决定导电类型、阈值电压
- **光吸收：** 带隙决定了光吸收的起始波长
- **催化：** 带边位置决定氧化还原电位

因此，准确预测带隙对材料科学至关重要。

### 1.4 改善带隙的替代方法

| 方法 | 带隙精度 | 计算量 | 适用范围 |
|------|----------|--------|----------|
| DFT+U | 中等 | 低 | 强关联体系 |
| 杂化泛函 (HSE06) | 好 | 中等 | 固体 |
| GW | 很好 | 高 | 通用 |
| DMFT | 很好 | 很高 | 强关联 |

---

## 2. GW 近似

### 2.1 多体微扰理论框架

GW 近似基于 Hedin 方程的循环自洽方案。在多体格林函数理论中，**自能 (Self-Energy)** $\Sigma$ 包含了所有交换-关联效应。

Hedin 方程的循环：

1. 格林函数 $G$
2. 屏蔽库仑相互作用 $W$
3. 自能 $\Sigma = iGW$
4. 顶角函数 $\Gamma$
5. 极化率 $P = -iG\Gamma G$

### 2.2 GW 近似

**GW 近似**保留了自能的最低阶项：

$$\Sigma(\mathbf{r}, \mathbf{r}', \omega) = \frac{i}{2\pi} \int G(\mathbf{r}, \mathbf{r}', \omega + \omega') W(\mathbf{r}, \mathbf{r}', \omega') d\omega'$$

其中：
- $G$ 是单粒子格林函数
- $W = \varepsilon^{-1} v$ 是屏蔽库仑相互作用（$v$ 为裸库仑势，$\varepsilon$ 为介电函数）

物理理解：电子在传播过程中与极化介质相互作用——$G$ 代表电子传播，$W$ 代表屏蔽的库仑相互作用。

### 2.3 G0W0（单次修正）

最常用的 GW 方案是 **G0W0**（one-shot GW），也称为单次 GW：

1. 先做一次 DFT 或 HF 计算，得到 $G_0$（零级格林函数，由 KS 轨道构建）
2. 用 $G_0$ 计算极化率和屏蔽库仑相互作用 $W_0$
3. 计算自能修正：$\Sigma^{G_0W_0} = iG_0W_0$
4. 准粒子能量：

$$\varepsilon_n^{G_0W_0} = \varepsilon_n^{\text{DFT}} + Z_n \langle \phi_n | \text{Re}[\Sigma(\varepsilon_n^{\text{DFT}})] - v_{xc} | \phi_n \rangle$$

其中 $Z_n$ 是重整化因子：

$$Z_n = \left[1 - \left.\frac{\partial \text{Re}\Sigma}{\partial \omega}\right|_{\varepsilon_n}\right]^{-1}$$

### 2.4 物理直觉

GW 修正的物理图像：

$$\varepsilon_n^{GW} = \varepsilon_n^{\text{KS}} + \underbrace{\langle \Sigma - v_{xc} \rangle}_{\text{自能修正}}$$

- $\Sigma$ 包含了精确的交换-关联效应（在 GW 近似下）
- $v_{xc}$ 是 DFT 的近似交换-关联势
- 修正项 $\langle \Sigma - v_{xc} \rangle$ 弥补了 DFT 的不足

对于带隙：
- 价带顶上移（减少束缚能）
- 导带底下移（增加束缚能，即增加电子亲和能）
- 带隙增大，接近实验值

### 2.5 屏蔽的计算

$W$ 中的介电函数 $\varepsilon(\omega)$ 通常在 **RPA (Random Phase Approximation)** 下计算：

$$\varepsilon^{-1}(\mathbf{G}, \mathbf{G}', \omega) = \delta_{\mathbf{G}\mathbf{G}'} + v(\mathbf{G}) \sum_{ia} \frac{(\phi_i | e^{i(\mathbf{G}-\mathbf{G}')\mathbf{r}} | \phi_a)(\phi_a | e^{-i(\mathbf{G}-\mathbf{G}')\mathbf{r}} | \phi_i)}{\varepsilon_a - \varepsilon_i - \omega - i\eta}$$

RPA 意味着极化率仅由独立粒子激发贡献（无局部场效应的交换部分）。

---

## 3. 自洽级别

### 3.1 G0W0（单次修正）

如上所述，G0W0 是从 DFT 出发的一次性修正。它的优势和局限：

**优势：**
- 计算量相对较小
- 如果起点选择得当（如 PBE 或 HSE06），结果很好

**局限：**
- 结果依赖于起始 DFT 交换关联泛函
- 对于某些体系（如过渡金属氧化物），不同起点导致的结果差异较大
- 不满足自洽条件

### 3.2 evGW0（部分自洽）

**特征值自洽 GW (Eigenvalue Self-Consistent GW, evGW0)** 只更新格林函数 $G$ 中的特征值，保持 W 固定：

1. DFT → 得到初始 $\{\varepsilon_n^0\}$
2. G0W0 → 得到 $\{\varepsilon_n^1\}$
3. 用 $\{\varepsilon_n^1\}$ 重建 $G_1$，保持 $W_0$ 不变
4. G1W0 → 得到 $\{\varepsilon_n^2\}$
5. 迭代直到收敛

**优势：**
- 消除了对起点泛函的依赖
- 计算量比 evGW 小（只更新 G，不更新 W）

### 3.3 evGW（完全自洽）

**完全自洽 evGW** 同时更新 G 和 W 中的特征值：

1. DFT → 初始
2. G0W0 → 更新 G 和 W
3. G1W1 → 更新
4. 直到收敛

**优势：**
- 完全自洽，无起点依赖
- 对于孤立分子，结果可靠

**局限：**
- 计算量大
- 对于固体，自洽有时会过估带隙

### 3.4 推荐策略

| 体系 | 推荐方案 | 原因 |
|------|----------|------|
| 分子 | G0W0@PBE 或 evGW0@PBE | 起点依赖性小 |
| 半导体 | G0W0@PBE (保守) 或 evGW0@PBE | 平衡精度和计算量 |
| 绝缘体 | evGW0@PBE 或 G0W0@HSE06 | 大带隙需更好起点 |
| 强关联 | 需要特殊处理 | GW 可能不够 |

---

## 4. BSE (Bethe-Salpeter 方程)

### 4.1 激子效应

**激子 (Exciton)** 是固体中的束缚电子-空穴对。在光吸收过程中，入射光子产生一个电子-空穴对。如果电子和空穴之间存在吸引相互作用，它们可以形成束缚态（激子）。

**为什么 GW 不够？**

GW 计算的是**准粒子能量**（添加/移除一个粒子的能量），但光学吸收涉及的是**电子-空穴对**的产生。电子和空穴之间的相互作用（特别是吸引）降低了激发能，导致：

- 光学带隙 < 准粒子带隙
- 在吸收光谱中，激子峰出现在带边以下

对于半导体，激子束缚能通常为几十 meV；对于绝缘体和分子，可以达到 1 eV 或更多。

### 4.2 Bethe-Salpeter 方程

**BSE (Bethe-Salpeter Equation)** 是描述电子-空穴对激发的标准方程：

$$\sum_{a'j'} H_{aj, a'j'}^{\text{BSE}} A_{a'j'}^S = \Omega_S A_{aj}^S$$

其中 BSE Hamiltonian：

$$H_{aj, a'j'}^{\text{BSE}} = (\varepsilon_a^{\text{QP}} - \varepsilon_j^{\text{QP}}) \delta_{aa'}\delta_{jj'} + K_{aj, a'j'}^{\text{dir}} - K_{aj, a'j'}^{\text{x}}$$

各项含义：
- $\varepsilon_a^{\text{QP}} - \varepsilon_j^{\text{QP}}$：对角项，准粒子能量差（来自 GW 计算）
- $K^{\text{dir}}$：直接（屏蔽）相互作用项——电子-空穴吸引，**降低**激发能
- $K^{\text{x}}$：交换（非屏蔽）相互作用项——仅对明态有贡献

### 4.3 光学性质

从 BSE 解可以获得**宏观介电函数**：

$$\varepsilon_2(\omega) = \frac{16\pi^2 e^2}{\omega^2} \sum_S |\langle 0 | \hat{\mathbf{e}} \cdot \mathbf{v} | S \rangle|^2 \delta(\omega - \Omega_S)$$

其中 $|S\rangle$ 是 BSE 激发态，$\hat{\mathbf{e}}$ 是光的偏振方向。

### 4.4 GW+BSE 的精度

GW+BSE 是目前计算固体光学性质的**黄金标准方法**：

| 材料 | 实验光学带隙 | GW+BSE 预测 | 误差 |
|------|-------------|-------------|------|
| Si | 3.4 eV (direct) | 3.2-3.5 eV | ~5% |
| LiF | 12.8 eV | 12.5-13.0 eV | ~2% |
| GaAs | 1.85 eV | 1.7-1.9 eV | ~5% |

### 4.5 TDDFT vs BSE

| 特性 | TDDFT | BSE |
|------|-------|-----|
| 理论框架 | 含时 KS 方程 | 多体格林函数 |
| 激子效应 | ALDA 无法描述 | 自然包含 |
| 计算量 | $O(N^3)$ | $O(N^4)$-$O(N^6)$ |
| 适用体系 | 分子 | 分子和固体 |
| 交换关联核 | $f_{xc}$ (近似) | 屏蔽库仑 $W$ (更准确) |
| 电荷转移 | 需要 RSH | 正确描述 |

---

## 5. CP2K 的 GW 实现

### 5.1 CP2K 中 GW 的三种方案

CP2K 提供了三种 GW 计算方案，适用于不同场景：

#### 方案 1：分子体系 (Molecular GW)

- 适用于孤立分子
- 标度为 $O(N^4)$
- 不涉及 k 点和周期性

#### 方案 2：周期体系 + k 点 (Periodic with k-points)

- 适用于固体
- 使用 k 点采样布里渊区
- 计算量大，但可以做能带结构

#### 方案 3：Gamma-only 大单胞 (Gamma-only Large Cell)

- 使用 Gamma 点近似
- 需要大单胞来消除镜像相互作用
- 适用于缺陷、表面等非完美周期体系

### 5.2 &RI_RPA 部分

GW 计算通过 `&WF_CORRELATION` 的 `&RI_RPA` 部分触发：

```fortran
&WF_CORRELATION
  METHOD RI_RPA                    # RI-RPA 方法（GW 的基础）
  &RI_RPA
    QUADRATURE_POINTS 20           # Minimax 积分点数
    &GW
      CORR_MOS_OCC   10            # 修正的占据轨道数
      CORR_MOS_VIRT  10            # 修正的虚轨道数
      ANALYTIC_CONTINUATION        # 解析延拓方案
    &END GW
    &RI
      RI_SIGMA 1.0E-6
    &END RI
  &END RI_RPA
  &INTERACTION_POTENTIAL
    POTENTIAL_TYPE TRUNCATED
    CUTOFF_RADIUS 5.0
  &END INTERACTION_POTENTIAL
  &MEMORY
    MAX_MEMORY 4000
  &END MEMORY
&END WF_CORRELATION
```

### 5.3 计算流程

1. **SCF 计算：** 获取 DFT 基态轨道和能量
2. **RI-RPA 能量：** 在 RI-RPA 框架下计算极化率
3. **GW 自能：** 从极化率计算屏蔽库仑相互作用 $W$ 和自能 $\Sigma$
4. **准粒子能量：** 对角修正：$\varepsilon^{GW} = \varepsilon^{DFT} + \langle \Sigma - v_{xc} \rangle$

### 5.4 HF 轨道预收敛

对于 GW 计算，有时推荐先用 HF 轨道作为起点（而非 PBE）。这可以：

- 改善 G0W0 的起点依赖性
- 特别是对于带隙较大的绝缘体

```fortran
! 步骤 1: 先做 HF 计算
&XC
  &XC_FUNCTIONAL NONE
  &END XC_FUNCTIONAL
  &HF
    FRACTION 1.0
    ...
  &END HF
&END XC

! 步骤 2: 用 HF 轨道做 GW
&WF_CORRELATION
  METHOD RI_RPA
  &RI_RPA
    &GW
      ...
    &END GW
  &END RI_RPA
&END WF_CORRELATION
```

---

## 6. 关键参数

### 6.1 QUADRATURE_POINTS

**Minimax 积分点数**控制频率积分（自能 $\Sigma(\omega)$ 中的频率卷积）的精度。

```fortran
&RI_RPA
  QUADRATURE_POINTS 20             ! 10-50 个点
&END RI_RPA
```

- 更多点数 → 更精确的频率积分 → 更长的计算时间
- 通常 10-20 个点足够
- 对于金属或小带隙体系，可能需要更多点

### 6.2 CORR_MOS_OCC 和 CORR_MOS_VIRT

指定需要修正的占据和虚轨道数：

```fortran
&GW
  CORR_MOS_OCC   10                ! 修正前 10 个占据轨道
  CORR_MOS_VIRT  10                ! 修正前 10 个虚轨道
&END GW
```

- 这决定了哪些轨道获得 GW 修正
- 对于带隙，至少需要修正 HOMO 和 LUMO
- 如果需要能带结构，需要覆盖所有感兴趣的能带

### 6.3 EPS_FILTER

积分筛选阈值：

```fortran
&RI_RPA
  EPS_FILTER 1.0E-5                ! 通常 1.0E-5 到 1.0E-8
&END RI_RPA
```

- 控制矩阵元素的截断阈值
- 较小值更精确但更慢
- 影响内存使用

### 6.4 MEMORY_PER_PROC

每个 MPI 进程的内存限制：

```fortran
&MEMORY
  MAX_MEMORY 4000                  ! 4 GB per process
&END MEMORY
```

GW 计算非常消耗内存，特别是极化率矩阵的存储。建议每进程至少 4-8 GB。

### 6.5 解析延拓方案

从虚频域到实频域的解析延拓：

```fortran
&GW
  ANALYTIC_CONTINUATION             ! 使用 Padé 解析延拓
  NUM_TIME_FREQ_POINTS 20           ! 时间-频率网格点数
&END GW
```

CP2K 的 GW 实现使用 minimax 方法在虚时/虚频域中进行积分，然后通过解析延拓得到实频自能。

---

## 7. 基组要求

### 7.1 后 HF 方法的基组需求

GW 和 BSE 属于后 HF 方法（多体微扰理论），因此对基组的要求与 MP2 类似：

- **最小要求：** aug-cc-pVDZ（含弥散函数的双ζ基组）
- **推荐：** cc-pVTZ 或 aug-cc-pVDZ
- **高精度：** aug-cc-pVTZ + 基组外推

### 7.2 为什么需要弥散函数？

GW 涉及准粒子波函数的尾部（远离原子核的区域），弥散函数对于正确描述这些区域至关重要。没有弥散函数：
- 介电函数被低估
- 屏蔽效应被低估
- 带隙被高估

### 7.3 CP2K 中的 GTH 赝势基组

CP2K 使用 GTH 赝势，配套基组为：

```fortran
&BASIS_SET cc-DZ                    ! 基础相关一致双ζ
&BASIS_RI_SET RI_cc-DZ              ! RI 辅助基组
```

CP2K 也提供了带弥散函数的基组（如 `cc-DZ-MOLOPT-SR` 和 `cc-TZ-MOLOPT-SR`）。

### 7.4 基组外推

与 MP2 类似，GW 准粒子能量也需要基组外推。常用的外推方案：

$$E_{\text{gap}}(X) = E_{\text{gap}}^{\text{CBS}} + \frac{A}{X^3}$$

使用 cc-DZ 和 cc-TZ 的结果进行两点外推。

### 7.5 平面波截断

在 GPW 框架下，平面波截断能 (`CUTOFF`) 也需要收敛：

```fortran
&MGRID
  CUTOFF 800                       ! GW 需要更高的截断能
  REL_CUTOFF 80
  NGRIDS 5
&END MGRID
```

通常 GW 计算需要比 DFT 更高的截断能（600-1000 Ry）。

---

## 8. 能带结构计算

### 8.1 k 点采样

对于周期性固体，能带结构沿布里渊区的高对称路径计算：

```fortran
&KPOINTS
  SCHEME MONKHORST-Pack 4 4 4      ! SCF 的 k 点网格
  &BAND_STRUCTURE
    NPOINTS 20                      ! 每段路径的点数
    &SPECIAL_POINT
      GAMMA  0.0  0.0  0.0         ! Γ 点
    &END SPECIAL_POINT
    &SPECIAL_POINT
      X     0.5  0.0  0.5          ! X 点
    &END SPECIAL_POINT
    &SPECIAL_POINT
      W     0.5  0.25 0.75         ! W 点
    &END SPECIAL_POINT
    &SPECIAL_POINT
      L     0.5  0.5  0.5          ! L 点
    &END SPECIAL_POINT
  &END BAND_STRUCTURE
&END KPOINTS
```

### 8.2 GW 能带结构

GW 能带结构计算分两步：

1. **SCF：** 用 DFT + 稠密 k 网格收敛基态
2. **GW 修正：** 在高对称路径上的 k 点，计算 GW 准粒子能量

这比 DFT 能带结构计算昂贵得多，因为每个 k 点都需要完整的 GW 计算。

### 8.3 布里渊区路径

常见晶体结构的高对称路径：

**FCC（如 Si, GaAs）：** $\Gamma - X - W - L - \Gamma - K$

**BCC：** $\Gamma - H - N - \Gamma - P - H$

**简单立方：** $\Gamma - X - M - \Gamma - R - X$

### 8.4 能带结构的解读

从 GW 能带结构可以提取：
- **带隙：** 价带顶到导带底的能量差
  - **直接带隙：** 价带顶和导带底在同一点（如 GaAs 在 Γ 点）
  - **间接带隙：** 不在同一点（如 Si 的 VBM 在 Γ，CBM 在 X 附近）
- **能带色散：** 有效质量 $m^* = \hbar^2 / (\partial^2 E / \partial k^2)$
- **能带宽度：** 决定了电子态密度

---

## 9. 输入文件详解

本章示例文件 `silicon_gw.inp` 对体硅进行 G0W0 带隙计算。以下逐节解释：

### 9.1 整体结构

```fortran
&GLOBAL
  PROJECT silicon_gw
  RUN_TYPE ENERGY
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    ...DFT 基态设置...
    &XC
      ...HF 或 PBE...
    &END XC
    &WF_CORRELATION
      &RI_RPA
        &GW ...&END GW
      &END RI_RPA
    &END WF_CORRELATION
  &END DFT
  &SUBSYS
    ...Si 结构和基组...
  &END SUBSYS
&END FORCE_EVAL
```

### 9.2 硅结构

```fortran
&CELL
  ABC 5.43 5.43 5.43              ! Si 晶格常数 (Angstrom)
  ALPHA 90.0
  BETA 90.0
  GAMMA 90.0
&END CELL

&COORD
  SCALED                          ! 使用分数坐标
  Si   0.00  0.00  0.00
  Si   0.50  0.50  0.00
  Si   0.50  0.00  0.50
  Si   0.00  0.50  0.50
  Si   0.25  0.25  0.25
  Si   0.75  0.75  0.25
  Si   0.75  0.25  0.75
  Si   0.25  0.75  0.75
&END COORD
```

硅的金刚石结构，8 个原子/单胞，晶格常数 5.43 Å。

### 9.3 HF 预收敛

```fortran
&XC
  &XC_FUNCTIONAL
    &NONE
    &END NONE
  &END XC_FUNCTIONAL
  &HF
    FRACTION 1.0
    &SCREENING
      EPS_SCHWARZ 1.0E-7
    &END SCREENING
    &INTERACTION_POTENTIAL
      POTENTIAL_TYPE TRUNCATED
      CUTOFF_RADIUS 3.0
    &END INTERACTION_POTENTIAL
  &END HF
&END XC
```

对于硅，推荐先做 HF 计算作为 GW 的起点。HF 轨道比 PBE 轨道更"局域化"，产生的 G0W0 带隙通常更接近实验值。

### 9.4 GW 设置

```fortran
&WF_CORRELATION
  METHOD RI_RPA
  &RI_RPA
    QUADRATURE_POINTS 16          ! Minimax 积分点数
    &GW
      CORR_MOS_OCC   8            ! 修正 8 个占据轨道
      CORR_MOS_VIRT  8            ! 修正 8 个虚轨道
    &END GW
    &RI
      RI_SIGMA 1.0E-6
    &END RI
  &END RI_RPA
  &INTERACTION_POTENTIAL
    POTENTIAL_TYPE TRUNCATED
    CUTOFF_RADIUS 3.0
  &END INTERACTION_POTENTIAL
  &MEMORY
    MAX_MEMORY 4000
  &END MEMORY
&END WF_CORRELATION
```

关键参数：
- `QUADRATURE_POINTS 16`：minimax 积分精度
- `CORR_MOS_OCC 8`：修正前 8 个占据轨道（Si 的 4 个价带 × 2 自旋）
- `CORR_MOS_VIRT 8`：修正前 8 个虚轨道
- `MAX_MEMORY 4000`：4 GB/进程

### 9.5 基组设置

```fortran
&KIND Si
  BASIS_SET cc-DZ                ! 相关一致双ζ基组
  BASIS_RI_SET RI_cc-DZ          ! RI 辅助基组
  POTENTIAL GTH-HF-q4            ! HF 赝势 (4 个价电子)
&END KIND
```

### 9.6 运行与分析

```bash
# 运行（需要大量内存和时间）
mpirun -np 16 cp2k.psmp -i silicon_gw.inp -o silicon_gw.out

# 查看 DFT 能量
grep "ENERGY| Total" silicon_gw.out | head -1

# 查看 GW 修正
grep "GW" silicon_gw.out
grep "QP" silicon_gw.out              ! 准粒子能量
grep "band gap" silicon_gw.out
```

### 9.7 预期结果

对于体硅（8 原子单胞，Gamma 点）：

| 方法 | 带隙 (eV) |
|------|-----------|
| PBE (DFT) | ~0.6 |
| HF | ~5-6 |
| G0W0@PBE | ~1.1-1.3 |
| G0W0@HF | ~1.2-1.4 |
| 实验值 | 1.17 |

G0W0 修正将 PBE 带隙从 ~0.6 eV 提升到 ~1.2 eV，与实验值 1.17 eV 非常接近。

### 9.8 进阶建议

1. **k 点收敛：** 使用更稠密的 k 网格（2x2x2, 4x4x4）检查带隙收敛
2. **基组外推：** 使用 cc-TZ 基组，通过两点外推得到 CBS 极限
3. **evGW0：** 进行特征值自洽计算，消除起点依赖
4. **BSE：** 在 GW 基础上做 BSE 计算，得到光学吸收谱
5. **缺陷态：** 使用大单胞 + Gamma 点近似计算点缺陷的能级

---

## 参考文献

1. L. Hedin, *New method for calculating the one-particle Green's function with application to the electron-gas problem*, Phys. Rev. **139**, A796 (1965).
2. M. S. Hybertsen and S. G. Louie, *Electron correlation in semiconductors and insulators: Band gaps and quasiparticle energies*, Phys. Rev. B **34**, 5390 (1986).
3. G. Onida, L. Reining, and A. Rubio, *Electronic excitations: density-functional versus many-body Green's-function approaches*, Rev. Mod. Phys. **74**, 601 (2002).
4. M. Del Ben et al., *Linear scaling DFT-based GW approximation for large systems*, Phys. Rev. B **93**, 035101 (2016).
5. M. Shishkin and G. Kresse, *Implementation and performance of the frequency-dependent GW method within the PAW framework*, Phys. Rev. B **74**, 035101 (2006).
6. The CP2K manual: https://manual.cp2k.org/

---

## 练习

1. **基础练习：** 运行 `silicon_gw.inp`，比较 DFT 带隙和 GW 带隙的差异。

2. **起点依赖性：** 分别以 PBE 和 HF 为起点做 G0W0，比较准粒子能量的差异。

3. **基组效应：** 将基组从 cc-DZ 改为 cc-TZ，观察带隙变化。

4. **截断能收敛：** 将 CUTOFF 从 500 改为 800，检查总能量和带隙的收敛。

5. **更多材料：** 将 Si 替换为 GaAs 或 Ge，计算它们的 GW 带隙。
