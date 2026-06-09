# 第八章：NVE 分子动力学

**难度：⭐⭐⭐**

前面几章我们计算了静态的电子结构和优化的几何构型。然而，真实世界中的原子总是在运动——它们在有限温度下振动、扩散、反应。**分子动力学 (Molecular Dynamics, MD)** 模拟原子在势能面上随时间的运动轨迹，是研究材料动态性质的强大工具。

---

## 目录

1. [分子动力学基础](#1-分子动力学基础)
2. [NVE 系综](#2-nve-系综)
3. [CP2K 的 MD 设置](#3-cp2k-的-md-设置)
4. [密度外推](#4-密度外推)
5. [初速度](#5-初速度)
6. [轨迹与能量文件](#6-轨迹与能量文件)
7. [NVE 能量守恒检验](#7-nve-能量守恒检验)
8. [改善能量守恒](#8-改善能量守恒)
9. [输入文件详解](#9-输入文件详解)
10. [练习](#10-练习)

---

## 1. 分子动力学基础

### 1.1 核心思想

分子动力学的基本思路非常直观：

1. 给定所有原子的初始位置和速度
2. 计算每个原子受到的力（从 DFT 电子结构计算得到）
3. 用牛顿运动方程更新原子位置和速度
4. 回到第 2 步，重复

$$\mathbf{F}_i = m_i \mathbf{a}_i = m_i \frac{d^2 \mathbf{R}_i}{dt^2}$$

这就是 **Born-Oppenheimer MD (BOMD)**——在每个时间步，电子都瞬间弛豫到基态，原子核在 Born-Oppenheimer 势能面上运动。

### 1.2 Velocity Verlet 算法

CP2K 使用 **Velocity Verlet** 积分算法，它是经典 MD 中最常用的算法之一：

$$\mathbf{R}(t + \Delta t) = \mathbf{R}(t) + \mathbf{v}(t)\Delta t + \frac{1}{2}\mathbf{a}(t)\Delta t^2$$

$$\mathbf{v}(t + \Delta t) = \mathbf{v}(t) + \frac{1}{2}[\mathbf{a}(t) + \mathbf{a}(t + \Delta t)]\Delta t$$

**算法流程：**

```
时间步 n:
  1. 已知 R(n), v(n), a(n)
  2. 计算新位置: R(n+1) = R(n) + v(n)·Δt + ½·a(n)·Δt²
  3. 用 DFT 计算 R(n+1) 处的力 → 得到 a(n+1)
  4. 计算新速度: v(n+1) = v(n) + ½·[a(n) + a(n+1)]·Δt
  5. 记录能量、速度等
  6. 进入时间步 n+1
```

### 1.3 时间步长 (Timestep)

时间步长 $\Delta t$ 是 MD 中最关键的参数之一：

| $\Delta t$ | 适用场景 |
|-------------|---------|
| 0.5 fs | 含氢体系（O-H 振动周期 ~10 fs）|
| 1.0 fs | 不含氢的轻元素体系 |
| 2.0 fs | 使用 SHAKE/RATTLE 约束的含氢体系 |
| 0.25 fs | 高频振动或高精度要求 |

**选择原则：** 时间步长应远小于体系中最快的振动周期。O-H 伸缩振动周期约 10 fs，因此 $\Delta t \leq 0.5$ fs 是安全的。

> **注意：** 1 fs = $10^{-15}$ s。在 CP2K 中，时间步长的单位是 **fs**。

---

## 2. NVE 系综

### 2.1 统计力学中的系综

在 MD 模拟中，**系综 (ensemble)** 定义了哪些物理量保持恒定：

| 系综 | 恒定量 | 名称 | 实验对应 |
|------|--------|------|---------|
| **NVE** | 粒子数 N, 体积 V, 能量 E | 微正则系综 | 孤立体系 |
| NVT | N, 体积 V, 温度 T | 正则系综 | 恒温实验 |
| NPT | N, 压力 P, 温度 T | 等温等压系综 | 常温常压实验 |
| NPH | N, 压力 P, 焓 H | 等焓等压系综 | 绝热加压 |

### 2.2 NVE 的特殊性

NVE（微正则）系综是最"纯粹"的 MD 模拟：

- **没有任何外部干预**——不使用恒温器 (thermostat) 或恒压器 (barostat)
- 总能量（动能 + 势能）理论上严格守恒
- 这使得 NVE 成为检验 MD 质量的**基准测试**

$$E_{\text{total}} = E_{\text{kinetic}} + E_{\text{potential}} = \text{const}$$

### 2.3 NVE 中的温度

虽然 NVE 中没有外部控制温度，但温度仍可通过**动能**定义：

$$E_{\text{kin}} = \frac{1}{2}\sum_i m_i |\mathbf{v}_i|^2 = \frac{3}{2}N k_B T_{\text{kin}}$$

在 NVE 中，温度会**自然波动**——这不是错误，而是统计力学的正常行为。温度的相对波动量级约为：

$$\frac{\delta T}{T} \sim \frac{1}{\sqrt{N}}$$

对于 3 个原子的水分子（N=3），温度波动很大；对于 1000 个原子，波动很小。

### 2.4 为什么从 NVE 开始？

NVE 是学习 MD 的理想起点：

1. **最简单的设置**：不需要配置恒温器/恒压器
2. **质量检验**：能量守恒程度直接反映模拟质量
3. **理解基础**：NVT 和 NPT 都是在 NVE 基础上添加约束
4. **物理透明**：结果直接反映势能面上的动力学

---

## 3. CP2K 的 MD 设置

### 3.1 基本结构

```
&GLOBAL
  PROJECT  H2O_nve
  RUN_TYPE MD              ! ★ 分子动力学
&END GLOBAL
```

### 3.2 &MD 部分

MD 参数在 `&MOTION` / `&MD` 中设置：

```
&MOTION
  &MD
    ENSEMBLE  NVE            ! 系综：NVE（微正则）
    STEPS  1000              ! 总步数
    TIMESTEP  0.5            ! 时间步长 (fs)
    TEMPERATURE  300.0       ! 初始温度 (K)，仅用于生成初速度
  &END MD
&END MOTION
```

### 3.3 关键参数

| 参数 | 说明 | 建议值 |
|------|------|--------|
| `ENSEMBLE` | 系综选择 | NVE, NVT, NPT |
| `STEPS` | MD 总步数 | 100~1000000 |
| `TIMESTEP` | 时间步长 (fs) | 0.5 fs (含氢) |
| `TEMPERATURE` | 初始温度 (K) | 300 K |

### 3.4 模拟总时间

$$t_{\text{total}} = \text{STEPS} \times \text{TIMESTEP}$$

例如：1000 步 × 0.5 fs = 500 fs = 0.5 ps

> **典型模拟时间尺度：**
> - 分子振动：~100 fs
> - 液体扩散：~10 ps
> - 蛋白质折叠：~1 μs
> - Ab initio MD 通常限于 ~100 ps（因为每个时间步都需要一次 DFT 计算）

---

## 4. 密度外推

### 4.1 问题：每个时间步都要 SCF

在 BOMD 中，每个时间步都需要进行完整的 SCF 计算来得到电子密度和能量。这对于 ab initio MD 来说非常昂贵。

### 4.2 解决方案：密度外推 (Density Extrapolation)

相邻时间步的原子位置变化很小，因此电子密度也变化很小。我们可以用前几步的密度来**外推**当前步的初始猜测，大幅减少 SCF 迭代次数。

CP2K 支持多种外推方案：

| 方案 | 关键字 | 说明 |
|------|--------|------|
| 线性外推 | `LINEAR_PS` | 用前 2 步密度线性外推 |
| 多项式外推 | `PS` | 用前 N 步密度多项式外推（推荐用于 MD）|
| ASPC | `ASPC` | Always Stable Predictor-Corrector |

### 4.3 设置密度外推

```
&QS
  EXTRAPOLATION  PS         ! 多项式外推
  EXTRAPOLATION_ORDER  3    ! 使用前 3 步的信息
&END QS
```

**参数选择：**
- `EXTRAPOLATION_ORDER 2`：使用前 2 步（线性外推）
- `EXTRAPOLATION_ORDER 3`：使用前 3 步（推荐，平衡精度和成本）
- `EXTRAPOLATION_ORDER 4`：更精确，但需要更多内存

### 4.4 外推的效果

使用密度外推后：

| | 无外推 | 有外推 (PS, order 3) |
|---|--------|---------------------|
| SCF 迭代步数 | 8~15 | 2~4 |
| 每个 MD 步时间 | 100% | ~30% |
| 能量守恒 | 更好 | 略差（但通常可接受）|

> **注意：** MD 第一步（或重启后的第一步）没有历史密度可用，SCF 会需要更多迭代。

---

## 5. 初速度

### 5.1 Maxwell-Boltzmann 分布

MD 模拟需要给每个原子赋予初始速度。在平衡态下，原子速度遵循 **Maxwell-Boltzmann 分布**：

$$P(v_\alpha) = \sqrt{\frac{m}{2\pi k_B T}} \exp\left(-\frac{m v_\alpha^2}{2 k_B T}\right)$$

CP2K 会根据 `TEMPERATURE` 关键字自动从 Maxwell-Boltzmann 分布中采样初始速度。

### 5.2 设置初始温度

```
&MD
  TEMPERATURE  300.0    ! 初始温度 (K)
&END MD
```

> **注意：** `TEMPERATURE` 仅用于**生成初速度**。在 NVE 系综中，之后温度会根据能量守恒自然变化。

### 5.3 从重启文件继续

如果 MD 模拟需要继续（例如 1000 步不够，再跑 1000 步），使用重启文件：

```
&EXT_RESTART
  RESTART_FILE_NAME  H2O_nve-1.restart    ! 重启文件
&END EXT_RESTART
```

重启文件包含：
- 最后一步的原子位置
- 最后一步的原子速度
- SCF 状态信息（用于密度外推）

---

## 6. 轨迹与能量文件

### 6.1 轨迹文件

`H2O_nve-pos-1.xyz`：记录每个时间步的原子坐标。

```
3                              ! 原子数
 i =       1, time =       0.250, E =      -17.1635744058
 O         0.000000    0.000000    0.117000
 H         0.000000    0.757000   -0.469000
 H         0.000000   -0.757000   -0.469000
3
 i =       2, time =       0.750, E =      -17.1635744012
 O         0.001234    0.000567    0.118234
 H        -0.000890    0.759123   -0.467890
 H         0.000345   -0.755234   -0.468123
...
```

每帧包含：
- 步号 `i`
- 时间 (fs)
- 总能量 (Ha)
- 原子坐标 (Å)

**可视化：**

```bash
# 使用 VMD 查看轨迹
vmd H2O_nve-pos-1.xyz

# 使用 OVITO
ovito H2O_nve-pos-1.xyz
```

### 6.2 能量文件

`H2O_nve-1.ener`：记录每个时间步的热力学量。

```
#     Step   Time[fs]       Kin.[a.u.]      Temp[K]           Pot.[a.u.]        Cons Qty[a.u.]    CPU time[s]
        0       0.000       0.005432123     300.000        -17.1635744058     -17.1581422828         1.234
        1       0.500       0.005398765     298.150        -17.1635410483     -17.1581422818         0.987
        2       1.000       0.005465432     301.325        -17.1636077714     -17.1581422822         0.956
        3       1.500       0.005512345     303.875        -17.1636546847     -17.1581422847         0.945
...
```

**各列含义：**

| 列 | 含义 |
|----|------|
| Step | 时间步编号 |
| Time [fs] | 模拟时间 |
| Kin. [a.u.] | 动能 (Hartree) |
| Temp [K] | 瞬时温度 |
| Pot. [a.u.] | 势能 (Hartree) |
| Cons Qty [a.u.] | 守恒量（NVE 中 = 总能量）|
| CPU time [s] | 每步计算时间 |

> **Cons Qty（守恒量）** 是 NVE 中最重要的检验指标——它应保持恒定。

---

## 7. NVE 能量守恒检验

### 7.1 为什么检验能量守恒？

在 NVE 系综中，总能量理论上严格守恒。数值误差（截断、迭代收敛不完全等）会导致能量漂移。能量守恒的程度直接反映模拟的**数值质量**。

### 7.2 定量标准

| 能量守恒程度 | 评价 | 典型情况 |
|-------------|------|---------|
| $\Delta E / E_{\text{kin}} < 10^{-4}$ | 优秀 | 精确设置 |
| $\Delta E / E_{\text{kin}} < 10^{-3}$ | 良好 | 推荐标准 |
| $\Delta E / E_{\text{kin}} < 10^{-2}$ | 可接受 | 可用于初步探索 |
| $\Delta E / E_{\text{kin}} > 10^{-2}$ | 差 | 需要改善参数 |

其中 $\Delta E$ 是总能量的波动范围，$E_{\text{kin}}$ 是平均动能。

### 7.3 分析方法

用 Python 分析 `.ener` 文件：

```python
import numpy as np
import matplotlib.pyplot as plt

# 读取 .ener 文件（跳过注释行）
data = np.loadtxt('H2O_nve-1.ener', comments='#')
step    = data[:, 0]
time    = data[:, 1]
ekin    = data[:, 2]
temp    = data[:, 3]
epot    = data[:, 4]
etot    = data[:, 5]

# 总能量守恒
e_avg = np.mean(etot)
e_drift = (etot[-1] - etot[0]) / np.mean(ekin)
e_fluct = (np.max(etot) - np.min(etot)) / np.mean(ekin)

print(f"平均总能量: {e_avg:.10f} Ha")
print(f"能量漂移 (相对): {e_drift:.2e}")
print(f"能量波动 (相对): {e_fluct:.2e}")

# 绘图
fig, axes = plt.subplots(3, 1, figsize=(10, 8))

axes[0].plot(time, etot)
axes[0].set_ylabel('Total Energy (Ha)')
axes[0].set_title('NVE Energy Conservation')

axes[1].plot(time, ekin, label='Kinetic')
axes[1].plot(time, epot - np.mean(epot), label='Potential (shifted)')
axes[1].set_ylabel('Energy (Ha)')
axes[1].legend()

axes[2].plot(time, temp)
axes[2].set_ylabel('Temperature (K)')
axes[2].set_xlabel('Time (fs)')

plt.tight_layout()
plt.savefig('nve_analysis.png', dpi=150)
plt.show()
```

### 7.4 能量守恒图解读

```
Total Energy (Ha)
│  ~~~~~~~~~~~~~~~~~~~~  ← 应近似水平直线
│ ~  ~  ~  ~  ~  ~  ~
│~  ~  ~  ~  ~  ~  ~  ~
└────────────────────────→ Time

Kinetic / Potential
│  ╱╲    ╱╲    ╱╲    ╱╲   ← 动能（与势能反相位振荡）
│ ╱  ╲  ╱  ╲  ╱  ╲  ╱  ╲
│╱    ╲╱    ╲╱    ╲╱    ╲
│╲    ╱╲    ╱╲    ╱╲    ╱
│ ╲  ╱  ╲  ╱  ╲  ╱  ╲  ╱  ← 势能
│  ╲╱    ╲╱    ╲╱    ╲╱
└────────────────────────→ Time
```

- **总能量**：应为近似水平线，微小振荡是正常的
- **动能与势能**：应反相位振荡（一个升高，另一个降低）
- **漂移**：如果总能量单调增加或减少，说明数值精度不足

---

## 8. 改善能量守恒

### 8.1 SCF 收敛

SCF 收敛不完全是能量不守恒的最常见原因：

```
&SCF
  EPS_SCF  1.0E-7       ! 更严格的 SCF 收敛（默认 1E-5 ~ 1E-6）
  MAX_SCF  50
&END SCF
```

### 8.2 网格精度

提高平面波截断能和网格精度：

```
&MGRID
  CUTOFF  400            ! 增大截断能（从 300 → 400 Ry）
  REL_CUTOFF  80         ! 增大相对截断能
  NGRIDS  5              ! 使用更多网格层次
&END MGRID
```

### 8.3 默认精度

提高全局精度参数：

```
&QS
  EPS_DEFAULT  1.0E-12   ! 更严格的默认精度
  EXTRAPOLATION  PS
  EXTRAPOLATION_ORDER  3
&END QS
```

### 8.4 时间步长

减小时间步长：

```
&MD
  TIMESTEP  0.25         ! 从 0.5 fs 减小到 0.25 fs
&END MD
```

更小的时间步长意味着：
- 每步原子位移更小
- 密度外推更准确
- 能量守恒更好
- 但模拟总时间相同时需要更多步

### 8.5 改善优先级

| 优先级 | 参数 | 效果 |
|--------|------|------|
| 1 | `EPS_SCF` | 最显著——SCF 收敛不好会导致能量跳跃 |
| 2 | `CUTOFF` / `REL_CUTOFF` | 影响力的计算精度 |
| 3 | `EPS_DEFAULT` | 影响所有数值积分精度 |
| 4 | `TIMESTEP` | 影响积分精度和外推准确性 |
| 5 | `EXTRAPOLATION_ORDER` | 影响 SCF 初始猜测质量 |

> **实用经验：** 如果能量守恒很差（$>10^{-2}$），先检查 `EPS_SCF`。如果已经不错但想更好，再调整其他参数。

---

## 9. 输入文件详解

### 9.1 H2O_nve.inp — 水分子 NVE 分子动力学

```bash
# 运行命令
cp2k.psmp -i H2O_nve.inp -o H2O_nve.out
```

完整输入文件及逐行注释：

```
&GLOBAL
  PROJECT  H2O_nve            ! 项目名称
  RUN_TYPE MD                  ! ★ 分子动力学
  PRINT_LEVEL LOW
&END GLOBAL

&FORCE_EVAL
  METHOD  QS                   ! DFT (Quickstep)

  &DFT
    BASIS_SET_FILE_NAME  BASIS_SET
    POTENTIAL_FILE_NAME  GTH_POTENTIALS

    &MGRID
      CUTOFF  300               ! 平面波截断能 (Ry)
      REL_CUTOFF  60
    &END MGRID

    &QS
      EPS_DEFAULT  1.0E-10      ! 全局精度
      EXTRAPOLATION  PS         ! ★ 多项式密度外推（MD 推荐）
      EXTRAPOLATION_ORDER  3    ! ★ 使用前 3 步信息外推
    &END QS

    &MIXING
      METHOD  BROYDEN_MIXING
      ALPHA  0.4
      NBUFFER  7
    &END MIXING

    &SCF
      SCF_GUESS  ATOMIC         ! 初始猜测（第一步用原子密度）
      EPS_SCF  1.0E-6           ! SCF 收敛判据
      MAX_SCF  50
    &END SCF

    &XC
      &XC_FUNCTIONAL PADE       ! LDA (PADE)
      &END XC_FUNCTIONAL
    &END XC
  &END DFT

  &SUBSYS
    &CELL
      ABC  12.4 12.4 12.4       ! 大盒子
      PERIODIC  NONE             ! 非周期性（气相水分子）
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
      O   0.000   0.000   0.117
      H   0.000   0.757  -0.469
      H   0.000  -0.757  -0.469
    &END COORD
  &END SUBSYS
&END FORCE_EVAL

&MOTION
  &MD
    ENSEMBLE  NVE               ! ★ NVE 系综（微正则）
    STEPS  1000                 ! 1000 个时间步
    TIMESTEP  0.5               ! 时间步长 0.5 fs
    TEMPERATURE  300.0           ! 初始温度 300 K（仅用于生成初速度）
  &END MD

  &PRINT
    &TRAJECTORY
      FORMAT  XYZ               ! 轨迹格式
      EVERY  1                  ! 每步输出
    &END TRAJECTORY

    &VELOCITIES
      FORMAT  XYZ               ! 速度轨迹
      EVERY  10                 ! 每 10 步输出
    &END VELOCITIES

    &FORCES
      FORMAT  XYZ
      EVERY  10
    &END FORCES
  &END PRINT
&END MOTION
```

### 9.2 关键参数总结

| 参数 | 值 | 作用 |
|------|-----|------|
| `RUN_TYPE` | MD | 启动分子动力学 |
| `ENSEMBLE` | NVE | 微正则系综 |
| `STEPS` | 1000 | 模拟 1000 步 |
| `TIMESTEP` | 0.5 fs | 每步 0.5 飞秒 |
| `TEMPERATURE` | 300 K | 初始温度 |
| `EXTRAPOLATION` | PS | 密度外推方案 |
| `EXTRAPOLATION_ORDER` | 3 | 外推阶数 |
| `EPS_SCF` | 1E-6 | SCF 收敛精度 |

### 9.3 输出文件列表

运行完成后，将生成以下文件：

| 文件 | 内容 |
|------|------|
| `H2O_nve.out` | 主输出文件（详细日志）|
| `H2O_nve-pos-1.xyz` | 原子位置轨迹 |
| `H2O_nve-vel-1.xyz` | 原子速度轨迹 |
| `H2O_nve-frc-1.xyz` | 原子力轨迹 |
| `H2O_nve-1.ener` | 热力学量（能量、温度）|
| `H2O_nve-1.restart` | 重启文件 |

---

## 10. 练习

### 练习 1：能量守恒分析

运行 `H2O_nve.inp`，分析 `.ener` 文件：
1. 绘制总能量 vs 时间
2. 计算能量漂移和波动
3. 绘制动能和势能的反相位振荡

### 练习 2：温度分析

从 `.ener` 文件中：
1. 计算平均温度
2. 计算温度的标准差
3. 与理论预期 $\delta T / T \sim 1/\sqrt{N}$ 比较

### 练习 3：时间步长的影响

分别用 $\Delta t = 0.25, 0.5, 1.0$ fs 运行 MD，比较：
- 能量守恒质量
- 计算时间
- 物理结果（温度、振动频率）差异

### 练习 4：改善能量守恒

从默认设置出发，逐步改善参数：
1. 只改 `EPS_SCF` (1E-6 → 1E-8)
2. 只改 `CUTOFF` (300 → 500)
3. 只改 `TIMESTEP` (0.5 → 0.25)
4. 全部改善

记录每种情况的能量守恒程度。

### 练习 5：长时间模拟

将 MD 步数增加到 5000 步（2.5 ps），观察：
- 能量漂移是否随时间累积
- 温度是否保持大致稳定
- 是否出现任何非物理行为

### 练习 6：可视化

使用 VMD 或 OVITO 打开轨迹文件：
1. 观察水分子的振动
2. 测量 O-H 键长的振动幅度
3. 估算 O-H 振动频率

---

## 小结

本章学习了 NVE 分子动力学的基础知识和 CP2K 实现：

| 概念 | CP2K 实现 |
|------|----------|
| MD 运行 | `RUN_TYPE MD` |
| NVE 系综 | `ENSEMBLE NVE` |
| 时间步长 | `TIMESTEP 0.5` (fs) |
| 总步数 | `STEPS 1000` |
| 初始温度 | `TEMPERATURE 300` (K) |
| 密度外推 | `EXTRAPOLATION PS`, `ORDER 3` |
| 能量守恒 | 检查 `.ener` 文件的 Cons Qty 列 |
| 轨迹输出 | `project-pos-1.xyz` |

**关键要点：**
1. NVE 是最纯粹的 MD，能量应严格守恒
2. 密度外推（PS, order 3）可大幅加速 MD
3. 能量守恒是判断 MD 质量的金标准
4. EPS_SCF 是影响能量守恒的最关键参数
5. 时间步长应小于体系最快振动周期的 1/10

下一章将介绍 **NVT 系综**和恒温器——如何在 MD 中控制温度。

---

**参考资源：**
- CP2K 手册 - &MD：https://manual.cp2k.org/trunk/CP2K_INPUT/MOTION/MD.html
- CP2K 手册 - Velocity Verlet：https://manual.cp2k.org/trunk/CP2K_INPUT/MOTION/MD.html
- Allen & Tildesley, "Computer Simulation of Liquids" — 经典 MD 教材
