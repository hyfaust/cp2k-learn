# 第9章：NVT/NPT分子动力学与热浴

> **难度：** ⭐⭐⭐  
> **前置知识：** 第8章 NVE分子动力学基础  
> **关键词：** 系综、Nose-Hoover、CSVR、Langevin、压浴、平衡化

---

## 目录

1. [系综理论](#1-系综理论)
2. [Nosé-Hoover热浴](#2-nosé-hoover热浴)
3. [CSVR热浴](#3-csvr热浴)
4. [速度重标定](#4-速度重标定)
5. [Langevin动力学](#5-langevin动力学)
6. [NPT与压浴](#6-npt与压浴)
7. [平衡化策略](#7-平衡化策略)
8. [三个输入文件详解](#8-三个输入文件详解)
9. [常见问题与调试](#9-常见问题与调试)
10. [进阶阅读](#10-进阶阅读)

---

## 1. 系综理论

### 1.1 什么是系综？

在统计力学中，**系综（Ensemble）** 是大量相同宏观条件下系统的集合。分子动力学模拟中，我们通过选择不同的系综来控制系统的宏观热力学条件：

| 系综 | 固定量 | 控制量 | 符号 | 典型应用 |
|------|--------|--------|------|----------|
| 微正则 | N, V, E | T 波动 | NVE | 验证能量守恒 |
| 正则 | N, V, T | E 波动 | NVT | 恒温平衡化 |
| 等温等压 | N, P, T | V, E 波动 | NPT | 液体/溶液模拟 |
| 等压等焓 | N, P, H | V, T 波动 | NPH | 相变研究 |

其中 N 是粒子数、V 是体积、E 是能量、T 是温度、P 是压力、H 是焓。

### 1.2 为什么需要热浴？

在NVE系综中，总能量守恒，但温度会随动能自然波动。对于实际应用，我们通常需要在特定温度下进行模拟：

- **实验条件：** 大多数实验在恒温条件下进行（如室温300K）
- **平衡化：** 需要将系统从初始构型弛豫到目标温度
- **采样：** 在特定温度下收集统计样本
- **相变研究：** 需要在不同温度下研究系统行为

**热浴（Thermostat）** 就是用来维持系统温度恒定的算法。它通过适当修改粒子速度来调节动能，从而控制温度。

### 1.3 温度的定义

瞬时温度通过动能定义：

$$T = \frac{2 \cdot E_{\text{kin}}}{k_B \cdot N_{\text{df}}}$$

其中：
- $E_{\text{kin}}$ 是总动能
- $k_B$ 是Boltzmann常数
- $N_{\text{df}}$ 是系统自由度数（通常为 $3N - 3$，减去质心平动）

在CP2K中，`TEMPERATURE` 关键词设置目标温度，`PRINT` 部分可以输出瞬时温度。

### 1.4 系综的选择指南

```
初始构型 → GEO_OPT (消除不合理接触)
         → NVE (短程，检查能量守恒)
         → NVT (平衡温度，5000-50000步)
         → NPT (平衡压力，如果需要的话)
         → Production (NVT或NPT，收集数据)
```

---

## 2. Nosé-Hoover热浴

### 2.1 理论基础

Nosé-Hoover热浴是最经典的确定性热浴方法之一，由Shuichi Nosé（1984）和William Hoover（1985）独立提出。

**核心思想：** 引入一个额外的自由度（热浴变量 $\eta$），将系统和热浴作为一个整体来处理。这个扩展系统的运动方程可以严格产生正则系综分布。

**扩展拉格朗日量：**

$$\mathcal{L} = \sum_i \frac{m_i s^2 \dot{\mathbf{r}}_i^2}{2} - V(\mathbf{r}) + \frac{Q \dot{s}^2}{2} - g k_B T \ln s$$

其中 $s$ 是热浴变量，$Q$ 是热浴"质量"，$g$ 是自由度数。

**等价的Nosé-Hoover运动方程：**

$$\dot{\mathbf{r}}_i = \frac{\mathbf{p}_i}{m_i}$$

$$\dot{\mathbf{p}}_i = \mathbf{F}_i - \eta \mathbf{p}_i$$

$$\dot{\eta} = \frac{1}{Q}\left(\sum_i \frac{\mathbf{p}_i^2}{m_i} - g k_B T\right)$$

这里 $\eta$ 是与热浴相关的摩擦系数。当系统动能高于目标值时，$\eta$ 增大，产生摩擦使动能降低；反之亦然。

### 2.2 热浴质量与TIMECON

热浴"质量" $Q$ 是一个关键参数，它决定了温度弛豫的时间尺度。在CP2K中，通过 `TIMECON` 参数间接控制：

$$Q = g k_B T \cdot \tau^2$$

其中 $\tau$ = `TIMECON`，称为**热浴振荡周期**。

**TIMECON的选择：**
- **太小（< 50 fs）：** 热浴响应太快，对系统动力学干扰过大，温度剧烈振荡
- **太大（> 1000 fs）：** 热浴响应太慢，温度平衡化需要很长时间
- **推荐值：** 100 fs（对于水等常见液体），通常在50-200 fs之间

### 2.3 链式热浴（Nose-Hoover Chain）

单个Nosé-Hoover热浴在小系统中可能无法产生遍历性（ergodicity）。解决方案是使用**链式热浴**：将多个热浴变量串联起来。

在CP2K中，可以通过 `NCOLS` 关键词控制链的长度：

```
&NOSE
  TIMECON 100.0    ! 热浴振荡周期 (fs)
  NCOLS 3          ! 链长度（默认1）
&END NOSE
```

### 2.4 CP2K中的Nosé-Hoover设置

在MD部分使用 `NOSE` 热浴类型：

```
&MOTION
  &MD
    ENSEMBLE NVT
    TEMPERATURE 300.0       ! 目标温度 (K)
    TIMESTEP 1.0            ! 时间步长 (fs)
    STEPS 5000              ! 总步数
    &THERMOSTAT
      REGION GLOBAL         ! 全局热浴
      &NOSE
        TIMECON 100.0       ! 振荡周期 (fs)
      &END NOSE
    &END THERMOSTAT
  &END MD
&END MOTION
```

### 2.5 全局 vs 局部热浴

CP2K支持两种热浴区域：

- **GLOBAL（全局热浴）：** 整个系统使用同一个热浴。适合较小的系统。
- **MASSIVE（质量热浴）：** 每个原子（或每组原子）使用独立的热浴。在NVE测试后、正式模拟中更推荐使用。

```
&THERMOSTAT
  REGION MASSIVE     ! 每个原子独立热浴
  &NOSE
    TIMECON 100.0
  &END NOSE
&END THERMOSTAT
```

---

## 3. CSVR热浴

### 3.1 理论基础

CSVR（Canonical Sampling through Velocity Rescaling）由Bussi、Donadio和Parrinello于2007年提出。它是一种**随机性**热浴，通过精确的速度重标定来产生正确的正则系综分布。

**核心思想：** 在每个时间步，对速度进行一次微小的随机重标定：

$$\mathbf{v}_{\text{new}} = \alpha \cdot \mathbf{v}_{\text{old}}$$

其中重标定因子 $\alpha$ 被精心设计，使得：
1. 正确采样正则系综的动能分布
2. 温度弛豫速率可控
3. 不会破坏系统的动力学行为

### 3.2 CSVR vs Nosé-Hoover

| 特性 | Nosé-Hoover | CSVR |
|------|-------------|------|
| 类型 | 确定性 | 随机性 |
| 系综采样 | 需要遍历性 | 严格正则 |
| 温度波动 | 系统性振荡 | 瞬时正确 |
| 平衡化速度 | 中等 | 快 |
| 动力学保真度 | 高 | 高 |
| 推荐场景 | 长时间生产模拟 | 平衡化和生产模拟 |

### 3.3 CP2K中的CSVR设置

```
&MOTION
  &MD
    ENSEMBLE NVT
    TEMPERATURE 300.0
    TIMESTEP 1.0
    STEPS 5000
    &THERMOSTAT
      REGION GLOBAL
      &CSVR
        TIMECON 100.0       ! 热浴弛豫时间 (fs)
      &END CSVR
    &END THERMOSTAT
  &END MD
&END MOTION
```

**TIMECON的物理意义：** 在CSVR中，TIMECON定义温度弛豫的特征时间。温度偏差以指数衰减：

$$\Delta T(t) \approx \Delta T(0) \cdot e^{-t/\tau}$$

- TIMECON = 100 fs 意味着温度偏差在约100 fs内衰减到 $1/e$
- 更小的TIMECON → 更快的温度控制 → 更大的动力学扰动

### 3.4 推荐使用场景

- **平衡化阶段：** CSVR特别适合快速将系统从非平衡态拉到目标温度
- **生产模拟：** 由于其随机性，CSVR能正确采样正则系综
- **初学者首选：** CSVR比Nosé-Hoover更容易设置和使用，不容易出错

---

## 4. 速度重标定

### 4.1 TEMP_TOL方法

在某些情况下（特别是NVE模拟中），我们需要一种简单的方法来防止温度漂移。`TEMP_TOL` 提供了一种"粗糙"的速度重标定方法：

```
&MOTION
  &MD
    ENSEMBLE NVE
    TEMPERATURE 300.0
    &THERMOSTAT
      REGION GLOBAL
      &NOSE
        TIMECON 100.0
      &END NOSE
    &END THERMOSTAT
    TEMP_TOL 50.0       ! 允许温度偏离50K
  &END MD
&END MOTION
```

### 4.2 TEMP_TOL的工作原理

当系统温度偏离目标值超过 `TEMP_TOL` 设定的范围时，所有粒子的速度会被均匀重标定，使温度回到目标值。

**适用场景：**
- 初始构型能量很高，需要逐步降温
- 防止NVE模拟中的温度漂移
- 简单的温度控制（不推荐用于正式模拟）

**注意事项：**
- TEMP_TOL不是严格的系综控制方法
- 它会产生非物理的温度跳跃
- 仅用于辅助目的，正式模拟应使用正统的热浴方法

---

## 5. Langevin动力学

### 5.1 理论基础

Langevin动力学通过在运动方程中添加摩擦力和随机力来实现温度控制：

$$m_i \ddot{\mathbf{r}}_i = \mathbf{F}_i - \gamma_i m_i \dot{\mathbf{r}}_i + \mathbf{R}_i(t)$$

其中：
- $\gamma_i$ 是摩擦系数（CP2K中通过 `GAMMA` 参数设置）
- $\mathbf{R}_i(t)$ 是随机力，满足涨落-耗散定理：$\langle \mathbf{R}_i(t) \cdot \mathbf{R}_j(t') \rangle = 6 k_B T m_i \gamma_i \delta_{ij} \delta(t-t')$

### 5.2 CP2K中的Langevin设置

```
&MOTION
  &MD
    ENSEMBLE LANGEVIN
    TEMPERATURE 300.0
    TIMESTEP 1.0
    STEPS 5000
    &LANGEVIN
      GAMMA 0.001        ! 摩擦系数 (fs^-1)
    &END LANGEVIN
  &END MD
&END MOTION
```

### 5.3 GAMMA参数

GAMMA控制摩擦力的强度：

- **GAMMA大（~0.01 fs⁻¹）：** 强摩擦，快速平衡化，但严重扰乱动力学
  - 适用：快速平衡化阶段
  - 类似于"在糖浆中游泳"

- **GAMMA小（~0.0001 fs⁻¹）：** 弱摩擦，对动力学扰动小，但平衡化慢
  - 适用：长时间生产模拟
  - 类似于"在稀薄气体中运动"

- **推荐值：** 平衡化用0.001-0.01 fs⁻¹，生产模拟用0.0001 fs⁻¹

### 5.4 Langevin动力学的适用场景

**优势：**
- 实现简单，参数少
- 快速达到目标温度
- 严格产生正则系综分布
- 不需要扩展拉格朗日量

**劣势：**
- 摩擦力改变粒子动力学
- 扩散系数、粘度等动力学性质需要修正
- 不适合研究真实的动力学过程

**推荐场景：**
- 快速平衡化（大GAMMA）
- 隐式溶剂模型中的采样
- 粗粒化模型
- 不关心真实动力学的热力学性质计算

---

## 6. NPT与压浴

### 6.1 等温等压系综

NPT系综同时控制温度和压力，是最接近实验条件的系综。实验中，样品通常暴露在恒定的大气压下，温度由恒温器控制。

在NPT中，盒子体积可以变化以维持目标压力：

$$\frac{dV}{dt} \propto (P_{\text{inst}} - P_{\text{target}})$$

### 6.2 压浴方法

CP2K中常用的压浴方法：

#### 6.2.1 Martyna-Tobias-Klein (MTK) 方法

MTK方法是Nosé-Hoover方法在NPT中的推广，产生严格的等温等压系综：

```
&MOTION
  &MD
    ENSEMBLE NPT_I    ! 各向同性压浴
    TEMPERATURE 300.0
    &BAROSTAT
      PRESSURE 1.01325    ! 目标压力 (bar), 1 atm
      TIMESCALE 1000.0    ! 压浴时间尺度 (fs)
    &END BAROSTAT
    &THERMOSTAT
      REGION GLOBAL
      &CSVR
        TIMECON 100.0
      &END CSVR
    &END THERMOSTAT
  &END MD
&END MOTION
```

#### 6.2.2 各向异性压浴

对于需要独立控制不同方向压力的情况：

```
&MD
  ENSEMBLE NPT_F    ! 各向异性压浴（自由形状）
  &BAROSTAT
    PRESSURE 1.01325
    TIMESCALE 1000.0
    PRESSURE_TENSOR 1.01325 1.01325 1.01325 0.0 0.0 0.0
  &END BAROSTAT
&END MD
```

**NPT_I vs NPT_F：**
- `NPT_I`：各向同性，盒子在三个方向均匀缩放。适合液体。
- `NPT_F`：各向异性，盒子形状可以改变。适合晶体、表面等。

### 6.3 TIMESCALE参数

`TIMESCALE` 控制压力弛豫的特征时间：

- **太小（< 100 fs）：** 盒子尺寸剧烈振荡，可能不稳定
- **太大（> 10000 fs）：** 压力平衡化极慢
- **推荐值：** 1000 fs 作为起点，根据系统调整

**经验值：**
- 液体：TIMESCALE = 500-2000 fs
- 晶体：TIMESCALE = 500-1000 fs
- 气体：通常不需要NPT

### 6.4 压力的计算

CP2K中压力（应力张量）的计算包括：
- 动能贡献（理想气体部分）
- virial（力与距离的乘积，来自原子间相互作用）

输出压力时，CP2K会在 `.ener` 文件中报告：
- `PRESSURE`：瞬时压力
- `PRESSURE_TENSOR`：完整的应力张量（3x3矩阵）

---

## 7. 平衡化策略

### 7.1 推荐工作流程

一个稳健的分子动力学模拟通常需要以下步骤：

```
Step 1: 能量最小化 (GEO_OPT)
        目的：消除初始构型中的不合理接触和过大的力
        关键：检查力的范数是否合理

Step 2: NVT平衡化
        目的：将系统从初始温度加热/冷却到目标温度
        热浴：CSVR，TIMECON 100 fs
        步数：5000-20000步
        监控：温度是否收敛到目标值

Step 3: NPT平衡化（如需要）
        目的：让盒子尺寸达到目标压力下的平衡值
        热浴：CSVR + 压浴
        步数：10000-50000步
        监控：压力和体积是否稳定

Step 4: 生产模拟
        目的：收集统计数据
        系综：NVT或NPT，取决于研究目标
        步数：取决于时间尺度
        注意：不要在生产阶段改变模拟参数！
```

### 7.2 平衡化的监控指标

| 指标 | 平衡标志 | 文件来源 |
|------|----------|----------|
| 温度 | 在目标值附近平稳波动 | .ener |
| 压力 | 波动减小，均值接近目标 | .ener |
| 体积 | 收敛到稳定值 | .ener |
| 势能 | 没有系统性漂移 | .ener |
| 密度 | 与实验值一致 | 计算得出 |
| RMSD | 构型稳定 | .xyz 分析 |

### 7.3 常见的平衡化错误

**错误1：跳过能量最小化**
```
# 不推荐：直接从初始构型开始MD
后果：大接触力导致粒子高速飞出，模拟崩溃
```

**错误2：TIMECON设置过小**
```
# 不推荐
&CSVR
  TIMECON 10.0     ! 太短，温度剧烈振荡
&END CSVR
```

**错误3：NPT中盒子太小**
```
# 水盒子中水分子太少，表面效应显著
# 建议至少使用 64 个水分子（~12.4 Å 盒子）
```

**错误4：平衡化时间不足**
```
# 只跑了1000步NVT就开始收集数据
# 应该至少观察温度和能量的时间序列，确认已平衡
```

### 7.4 温度与压力的典型波动

对于 $N$ 个原子的系统，温度的相对波动为：

$$\frac{\sigma_T}{\langle T \rangle} \approx \sqrt{\frac{2}{3N}}$$

- 72个原子（24个水）：$\sigma_T / T \approx 9.6\%$（300K时波动约29K）
- 192个原子（64个水）：$\sigma_T / T \approx 5.9\%$（300K时波动约18K）

压力的波动通常非常大（数百bar），这是正常现象。小系统中压力波动更大。

---

## 8. 三个输入文件详解

### 8.1 H2O_nvt.inp — Nosé-Hoover NVT

这个输入文件使用Nosé-Hoover热浴对24个水分子进行NVT分子动力学模拟。

**关键设置：**
- `ENSEMBLE NVT`：正则系综
- `NOSE` 热浴，TIMECON = 100 fs
- 目标温度300K，5000步，1 fs时间步长
- 使用DFTB（近似DFT方法）计算力
- GTH-DZVP基组和势

**注意事项：**
- 这是一个教学示例，模拟时间很短（5 ps）
- 实际研究中通常需要更长的模拟时间
- Nosé-Hoover是确定性的，温度会有系统性振荡

### 8.2 H2O_nvt_csvr.inp — CSVR NVT

与Nosé-Hoover版本类似，但使用CSVR热浴。

**区别：**
- CSVR是随机性热浴
- 温度分布更快收敛到正确的正则分布
- 没有温度振荡问题
- 更适合平衡化和一般用途

### 8.3 H2O_npt.inp — NPT模拟

在NVT基础上增加压浴控制，实现等温等压系综。

**额外设置：**
- `ENSEMBLE NPT_I`：各向同性等温等压系综
- `BAROSTAT` 部分：目标压力1.01325 bar（1 atm）
- 压浴时间尺度1000 fs
- CSVR热浴保持温度

**输出分析：**
- `.ener` 文件包含温度、压力、体积时间序列
- 可以计算平衡密度并与实验值比较

---

## 9. 常见问题与调试

### 9.1 温度不稳定

**症状：** 温度远离目标值或剧烈振荡

**可能原因和解决方法：**
1. 能量最小化不充分 → 重新进行更严格的GEO_OPT
2. 时间步长太大 → 减小TIMESTEP（对于含氢系统，建议1 fs或更小）
3. 热浴TIMECON太小 → 增大TIMECON到100 fs以上
4. 系统太小 → 使用更大的系统（>64个原子）

### 9.2 模拟崩溃（原子重叠）

**症状：** CP2K报错"SCF did not converge"或力的范数异常大

**解决方法：**
1. 检查初始构型是否合理
2. 先做NVT平衡化，再做NPT
3. 使用更小的时间步长
4. 添加 `MAX_FORCE 0.01` 限制最大力

### 9.3 压力不收敛（NPT）

**症状：** 体积持续增大或减小

**可能原因和解决方法：**
1. TIMESCALE太小 → 增大到1000 fs
2. 系统太小 → 小系统的压力波动很大是正常的
3. 使用NPT_I而非NPT_A（各向异性可能在某些方向不稳定）

### 9.4 如何选择热浴？

```
需要严格正则系综？
├── 是 → CSVR
└── 否 → 需要确定性动力学？
    ├── 是 → Nosé-Hoover
    └── 否 → 需要快速平衡化？
        ├── 是 → Langevin (大GAMMA) 或 CSVR (小TIMECON)
        └── 否 → CSVR (标准TIMECON)
```

---

## 10. 进阶阅读

### 10.1 参考文献

1. **Nosé, S.** (1984). A unified formulation of the constant temperature molecular dynamics methods. *J. Chem. Phys.*, 81(1), 511-519.
2. **Hoover, W. G.** (1985). Canonical dynamics: Equilibrium phase-space distributions. *Phys. Rev. A*, 31(3), 1695.
3. **Bussi, G., Donadio, D., & Parrinello, M.** (2007). Canonical sampling through velocity rescaling. *J. Chem. Phys.*, 126(1), 014101.
4. **Martyna, G. J., Klein, M. L., & Tuckerman, M.** (1992). Nosé-Hoover chains: The canonical constant-temperature molecular dynamics method. *J. Chem. Phys.*, 97(4), 2635-2643.

### 10.2 CP2K相关文档

- CP2K Manual: MOTION > MD > ENSEMBLE
- CP2K Manual: MOTION > MD > THERMOSTAT
- CP2K Manual: MOTION > MD > BAROSTAT
- CP2K HOWTO: Molecular Dynamics

### 10.3 后续章节预告

- **第10章：** 元动力学增强采样 — 使用元动力学方法克服自由能势垒
- **第11章：** xTB半经验方法 — 快速进行大尺度分子动力学
- **第12章：** QM/MM混合方法 — 处理溶剂化、酶催化等复杂体系

---

## 附录：关键词速查表

| 关键词 | 位置 | 说明 |
|--------|------|------|
| `ENSEMBLE` | `MOTION > MD` | 系综类型：NVE, NVT, NPT_I, NPT_F, LANGEVIN |
| `TEMPERATURE` | `MOTION > MD` | 目标温度 (K) |
| `TIMESTEP` | `MOTION > MD` | 时间步长 (fs) |
| `STEPS` | `MOTION > MD` | 总步数 |
| `THERMOSTAT` | `MOTION > MD` | 热浴设置块 |
| `REGION` | `THERMOSTAT` | 热浴区域：GLOBAL, MASSIVE |
| `NOSE` | `THERMOSTAT` | Nosé-Hoover热浴 |
| `CSVR` | `THERMOSTAT` | CSVR热浴 |
| `TIMECON` | `NOSE` / `CSVR` | 热浴振荡周期 (fs) |
| `NCOLS` | `NOSE` | 链式热浴长度 |
| `BAROSTAT` | `MOTION > MD` | 压浴设置块 |
| `PRESSURE` | `BAROSTAT` | 目标压力 (bar) |
| `TIMESCALE` | `BAROSTAT` | 压浴时间尺度 (fs) |
| `LANGEVIN` | `MOTION > MD` | Langevin热浴设置 |
| `GAMMA` | `LANGEVIN` | 摩擦系数 (fs⁻¹) |
| `TEMP_TOL` | `MOTION > MD` | 温度容差 (K) |
