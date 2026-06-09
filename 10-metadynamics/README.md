# 第10章：元动力学增强采样

> **难度：** ⭐⭐⭐⭐  
> **前置知识：** 第8-9章 分子动力学基础与热浴  
> **关键词：** 元动力学、集合变量、自由能面、高斯沉积、HILLS文件

---

## 目录

1. [增强采样问题](#1-增强采样问题)
2. [元动力学原理](#2-元动力学原理)
3. [集合变量（Collective Variables）](#3-集合变量collective-variables)
4. [CP2K中的元动力学设置](#4-cp2k中的元动力学设置)
5. [高斯参数选择](#5-高斯参数选择)
6. [反射壁](#6-反射壁)
7. [预运行](#7-预运行)
8. [FES重建](#8-fes重建)
9. [丙氨酸二肽示例](#9-丙氨酸二肽示例)
10. [进阶阅读与常见问题](#10-进阶阅读与常见问题)

---

## 1. 增强采样问题

### 1.1 稀有事件问题

分子动力学（MD）模拟的一个核心局限性是**时间尺度问题**。许多重要的生物学和化学过程涉及跨越自由能势垒，其特征时间远超MD可直接模拟的范围：

| 过程 | 典型时间尺度 | MD可直接模拟？ |
|------|-------------|---------------|
| 化学键振动 | ~10 fs | 是 |
| 溶剂化壳层重排 | ~10 ps | 是 |
| 蛋白质侧链旋转 | ~1 ns | 是 |
| 蛋白质折叠 | ~μs - ms | 勉强 |
| 化学反应 | ~μs - s | 困难 |
| 药物解离 | ~ms - s | 不可能 |

标准MD按照Boltzmann分布采样构型空间。如果自由能面（Free Energy Surface, FES）存在高势垒，系统将被困在局部最小值附近，无法在合理时间内访问所有相关构型。

### 1.2 自由能面

自由能面 $F(\xi)$ 是集合变量 $\xi$ 的函数：

$$F(\xi) = -k_B T \ln P(\xi)$$

其中 $P(\xi)$ 是在集合变量 $\xi$ 处的概率分布。自由能面的地形决定了：
- **最小值：** 稳定构象或反应物/产物
- **鞍点/过渡态：** 连接不同最小值的路径上的最高点
- **势垒高度：** 决定转换速率，通过Arrhenius方程 $k \propto e^{-\Delta F / k_B T}$

### 1.3 增强采样方法概览

为克服时间尺度限制，发展了多种增强采样方法：

| 方法 | 原理 | 优点 | 缺点 |
|------|------|------|------|
| **元动力学** | 添加历史依赖偏差势 | 可重建完整FES | 需要选择CV |
| 平行回火 | 不同温度间交换构型 | 无需CV | 需要大量副本 |
| 伞形采样 | 沿CV施加偏置势 | 精确 | 需要预定义路径 |
| 加速MD | 降低势垒 | 简单 | 难以定量重建FES |
| 适应性偏置力 | 基于力的偏置 | 精确 | 实现复杂 |

元动力学是其中最流行的方法之一，因为它能系统地重建完整的自由能面。

---

## 2. 元动力学原理

### 2.1 基本思想

元动力学（Metadynamics）由Laio和Parrinello于2002年提出。其核心思想是：

> 在MD模拟过程中，定期在集合变量空间中沉积小的高斯偏差势。这些高斯势逐步"填满"自由能最小值，迫使系统探索新的构型区域。当模拟结束时，所有最小值都被填平，累积的偏差势就是自由能面的负像。

### 2.2 数学表达

偏差势的时间演化：

$$V(\xi, t) = \sum_{t' = \tau_G, 2\tau_G, \ldots}^{t} W \cdot \exp\left(-\frac{[\xi(\mathbf{r}) - \xi(t')]^2}{2\sigma^2}\right)$$

其中：
- $W$ 是高斯高度（`WW`），单位为能量（通常为Hartree）
- $\sigma$ 是高斯宽度（`SCALE`），单位与CV相同
- $\tau_G$ 是沉积间隔（`NT_HILLS`步，每步一个时间步长）
- $\xi(\mathbf{r})$ 是当前构型的CV值

### 2.3 收敛性

当偏差势足以填平所有最小值时，模拟收敛：

$$V(\xi, t \to \infty) \approx -F(\xi) + C$$

其中 $C$ 是一个常数。因此，自由能面可以通过偏差势的负值来估算：

$$F(\xi) \approx -V(\xi, t_{\text{final}}) + C$$

### 2.4 Well-Tempered Metadynamics

标准元动力学存在过度填充（overfilling）的风险。**Well-Tempered Metadynamics**通过逐渐减小高斯高度来解决这个问题：

$$W(t) = W_0 \cdot \exp\left(-\frac{V(\xi, t)}{k_B \Delta T}\right)$$

其中 $\Delta T$ 是一个"虚假温度"参数，$\gamma = (T + \Delta T)/T$ 是偏差因子。这保证了收敛到缩放的自由能面：

$$V(\xi, t \to \infty) = -\frac{\Delta T}{T + \Delta T} F(\xi) + C$$

在CP2K中，通过设置 `DO_WELLTEMPERED .TRUE.` 和 `BIAS_TEMPERATURE` 来启用。

### 2.5 元动力学的优缺点

**优点：**
- 能够系统地探索高维集合变量空间
- 可以从偏差势直接重建自由能面
- 算法简单，实现直观
- 与其他方法（如REMD）可以组合使用

**缺点：**
- 需要选择合适的集合变量（最关键也最困难）
- 高维CV空间（>3D）收敛困难
- 参数选择需要经验
- 收敛判据不明确（需要后验分析）

---

## 3. 集合变量（Collective Variables）

### 3.1 CV选择的原则

集合变量（Collective Variable, CV）的选择是元动力学模拟中最关键的决定。好的CV应该满足：

1. **区分性：** 能区分反应物态、产物态和过渡态
2. **描述性：** 能捕捉系统的关键自由度
3. **低维性：** 维度尽可能低（1-3维），便于收敛
4. **连续性：** 关于原子位置是连续可微的
5. **相关性：** 与感兴趣的物理/化学过程密切相关

### 3.2 常用CV类型

#### 3.2.1 距离（Distance）

两个原子或原子组之间的距离：

$$d = |\mathbf{r}_i - \mathbf{r}_j|$$

适用于：键的形成/断裂、配位过程、分子间相互作用

在CP2K中通过 `&DISTANCE` 定义：
```
&COLVAR
  &DISTANCE
    ATOMS 1 2
  &END DISTANCE
&END COLVAR
```

#### 3.2.2 配位数（Coordination Number）

原子 $i$ 周围配位原子的数目，通过切换函数平滑化：

$$CN_i = \sum_{j \neq i} \frac{1}{1 + \exp[n(r_{ij} - r_0)]}$$

其中 $r_0$ 是特征距离，$n$ 控制切换函数的锐度。

适用于：溶剂化结构、氢键网络、金属配位

在CP2K中通过 `&COORDINATION` 定义：
```
&COLVAR
  &COORDINATION
    KIND_A O
    KIND_B H
    R_0 2.5       ! 特征距离 (Bohr)
    NN 8           ! 切换函数指数
    ND 14
  &END COORDINATION
&END COLVAR
```

#### 3.2.3 二面角/扭转角（Torsion / Dihedral）

四个原子定义的二面角：

$$\phi = \arctan2\left(\frac{(\mathbf{r}_{12} \times \mathbf{r}_{23}) \times (\mathbf{r}_{23} \times \mathbf{r}_{34})}{|\mathbf{r}_{23}|}, (\mathbf{r}_{12} \times \mathbf{r}_{23}) \cdot (\mathbf{r}_{23} \times \mathbf{r}_{34})\right)$$

适用于：构象变化（phi/psi角）、手性中心翻转、旋转异构体

在CP2K中通过 `&TORSION` 定义：
```
&COLVAR
  &TORSION
    ATOMS 5 7 9 15
  &END TORSION
&END COLVAR
```

#### 3.2.4 角度（Angle）

三个原子定义的夹角：

$$\theta = \arccos\left(\frac{\mathbf{r}_{12} \cdot \mathbf{r}_{32}}{|\mathbf{r}_{12}||\mathbf{r}_{32}|}\right)$$

适用于：弯曲运动、反应坐标

```
&COLVAR
  &ANGLE
    ATOMS 1 2 3
  &END ANGLE
&END COLVAR
```

#### 3.2.5 其他CV

- **RMSD：** 与参考结构的均方根偏差
- **GYRATION RADIUS：** 回转半径，描述蛋白质紧密程度
- **PATH：** 沿预定义路径的进展参数
- **DRMSD：** 距离矩阵的RMSD

### 3.3 CV选择的实用建议

**初学者指南：**

1. **从简单开始：** 先用1-2个CV进行探索
2. **物理直觉：** 基于对过程的理解选择CV
3. **文献参考：** 查看类似系统的已有研究
4. **预运行分析：** 先做NVT模拟，分析CV的波动范围
5. **逐步增加：** 如果2D不够，逐步添加CV，但避免超过3D

**常见错误：**
- CV太多导致维度灾难（curse of dimensionality）
- CV与反应坐标无关，无法区分产物和反应物
- CV变化不连续（如跨越周期性边界）

---

## 4. CP2K中的元动力学设置

### 4.1 主要输入结构

元动力学的设置位于 `MOTION > FREE_ENERGY > METADYN` 部分：

```
&MOTION
  &FREE_ENERGY
    &METADYN
      DO_HILLS .TRUE.              ! 启用高斯沉积
      NT_HILLS 50                  ! 每50步沉积一个高斯
      WW 0.001                     ! 高斯高度 (Hartree)
      LAGRANGE                   .TRUE.      ! 使用拉格朗日方法约束CV
      &METAVAR
        COLVAR 1                   ! 引用第1个COLVAR
        SCALE 0.1                  ! 高斯宽度 (CV单位)
      &END METAVAR
      &METAVAR
        COLVAR 2                   ! 引用第2个COLVAR
        SCALE 0.1                  ! 高斯宽度 (CV单位)
      &END METAVAR
    &END METADYN
  &END FREE_ENERGY
&END MOTION
```

### 4.2 关键参数说明

#### DO_HILLS
- `.TRUE.`：启用高斯沉积（正式元动力学）
- `.FALSE.`：仅监控CV，不添加偏差（预运行模式）

#### NT_HILLS
- 高斯沉积的频率（每多少个MD步沉积一个）
- 太频繁：偏差势增长太快，可能导致系统被强制推向不合理构型
- 太稀少：模拟时间大大增加
- 推荐值：50-500步（取决于时间步长和系统弛豫时间）

#### WW（Well Width / 高斯高度）
- 每个高斯势的能量高度，单位Hartree
- 太大：偏差势增长过快，可能跳过重要特征
- 太小：需要极长模拟时间才能覆盖FES
- 推荐值：0.0001 - 0.005 Hartree（约0.06 - 3 kcal/mol）

#### SCALE（高斯宽度）
- 每个高斯在CV空间中的宽度，单位与CV相同
- 太大：平滑掉FES的细节特征
- 太小：FES表面噪声大，需要更多高斯
- 推荐值：约为CV在最小值附近波动幅度的1/3到1/2

### 4.3 CV定义块

在输入文件的 `FORCE_EVAL > SUBSYS > COLVAR` 部分定义集合变量：

```
&SUBSYS
  ...
  &COLVAR
    ! 第1个CV: 二面角 (phi)
    &TORSION
      ATOMS 5 7 9 15
    &END TORSION
  &END COLVAR
  &COLVAR
    ! 第2个CV: 二面角 (psi)
    &TORSION
      ATOMS 7 9 15 17
    &END TORSION
  &END COLVAR
&END SUBSYS
```

### 4.4 METADYN中的METAVAR与COLVAR的关系

`METAVAR` 块中通过 `COLVAR` 索引引用 `SUBSYS` 中定义的集合变量：
- `COLVAR 1` 引用第一个 `&COLVAR` 块
- `COLVAR 2` 引用第二个 `&COLVAR` 块
- 以此类推

### 4.5 Well-Tempered Metadynamics设置

```
&METADYN
  DO_HILLS .TRUE.
  NT_HILLS 50
  WW 0.001
  DO_WELLTEMPERED .TRUE.
  BIAS_TEMPERATURE 1500.0     ! 虚假温度 (K), 通常是实际温度的3-5倍
  &METAVAR
    COLVAR 1
    SCALE 0.1
  &END METAVAR
  &METAVAR
    COLVAR 2
    SCALE 0.1
  &END METAVAR
&END METADYN
```

**BIAS_TEMPERATURE的选择：**
- 越大 → 更慢的填充，更好的收敛性，FES更平滑
- 越小 → 更快的探索，但可能过度填充
- 推荐：实际温度的3-10倍

---

## 5. 高斯参数选择

### 5.1 参数选择的权衡

元动力学参数的选择存在内在矛盾：

```
探索速度 ←——————————————→ 采样精度

  大WW, 小NT_HILLS         小WW, 大NT_HILLS
  快速覆盖FES               精细采样
  可能跳过细节               需要很长模拟
  适合初步探索               适合精确FES
```

### 5.2 推荐的参数选择策略

**Step 1: 估计CV的波动范围**

先做标准MD（不开元动力学），分析CV的时间序列：
- 二面角：波动范围通常 ±30°
- 距离：波动范围取决于键的刚度
- 配位数：波动范围 ±0.5

**Step 2: 设定SCALE**

```
SCALE ≈ CV波动范围 × 0.5
```

例如，phi角在最小值附近波动约 ±20°（0.35 rad），则 SCALE ≈ 0.17

**Step 3: 设定WW**

初始选择：
```
WW ≈ kBT × (CV波动范围 / SCALE)^2 × 0.01
```

在300K时，$k_BT$ ≈ 0.00094 Hartree（约0.6 kcal/mol）

对于典型的二面角CV：WW ≈ 0.001 Hartree 是一个合理的起点。

**Step 4: 设定NT_HILLS**

```
NT_HILLS × TIMESTEP ≈ CV弛豫时间 × 5-10
```

例如，对于二面角弛豫时间约 1 ps，TIMESTEP = 1 fs：
```
NT_HILLS ≈ 1000 / 1 × 5 = 500
```

但这通常太慢。实践中，NT_HILLS = 50-100 是常见的起点。

### 5.3 参数敏感性测试

建议用不同参数组合进行短时间测试：

```
测试1: WW=0.0005, SCALE=0.05, NT_HILLS=100
测试2: WW=0.001,  SCALE=0.1,  NT_HILLS=50
测试3: WW=0.002,  SCALE=0.2,  NT_HILLS=100
```

比较FES结果，选择能平衡分辨率和收敛性的参数。

---

## 6. 反射壁

### 6.1 为什么需要反射壁？

在元动力学模拟中，CV可能演化到物理上不合理的区域。例如：
- 二面角接近 ±180°（跨越周期性边界）
- CV远离反应区域，浪费计算时间
- CV进入高能区域，可能导致系统崩溃

**反射壁（Reflecting Walls）** 在CV达到边界时施加额外的排斥势，防止其继续前进。

### 6.2 CP2K中的反射壁设置

```
&METADYN
  ...
  &WALL
    TYPE QUADRATIC        ! 二次势壁
    POSITION 3.14         ! 壁的位置 (CV单位)
    K 0.1                 ! 力常数 (Hartree/CV单位^2)
    &WALL_ENV
      COLVAR 1            ! 对第1个CV施加壁
      DIRECTION POSITIVE  ! 只在正方向施加
    &END WALL_ENV
  &END WALL
  &WALL
    TYPE QUADRATIC
    POSITION -3.14
    K 0.1
    &WALL_ENV
      COLVAR 1
      DIRECTION NEGATIVE  ! 只在负方向施加
    &END WALL_ENV
  &END WALL
&END METADYN
```

### 6.3 反射壁的类型

- **QUADRATIC：** 二次排斥势，$V_{\text{wall}} = K \cdot (s - s_{\text{wall}})^2$，当 $s > s_{\text{wall}}$
- **HARMONIC：** 简谐势形式
- **WELLTEMP：** Well-tempered形式的壁

### 6.4 实用建议

- 壁的位置应该在CV感兴趣的范围之外
- 力常数 $K$ 要足够大使CV不会明显穿透壁，但不能太大导致数值不稳定
- 对于二面角CV（范围 -180° 到 +180°），由于周期性，通常不需要壁
- 对于距离CV，需要设置合理的上下限

---

## 7. 预运行

### 7.1 预运行的重要性

在正式启动元动力学之前，建议进行**预运行（Pre-production run）**：

1. **确认CV波动范围：** 确保CV在没有偏差的情况下在合理范围内波动
2. **确认系统稳定性：** 检查MD是否稳定运行
3. **设定参数基准：** 根据CV波动幅度设定SCALE
4. **检查初始构型：** 确保系统已经充分平衡化

### 7.2 预运行设置

在CP2K中，设置 `DO_HILLS .FALSE.` 即为预运行模式：

```
&METADYN
  DO_HILLS .FALSE.          ! 关闭高斯沉积
  NT_HILLS 50
  WW 0.001
  &METAVAR
    COLVAR 1
    SCALE 0.1
  &END METAVAR
  &METAVAR
    COLVAR 2
    SCALE 0.1
  &END METAVAR
&END METADYN
```

### 7.3 预运行输出分析

预运行的输出文件中应该检查：

**HILLS文件（当DO_HILLS=.FALSE.时只记录CV值）：**
- CV的平均值和波动范围
- CV是否在预期范围内
- CV之间是否有强相关性

**ener文件：**
- 温度是否稳定
- 势能没有系统性漂移
- SCF收敛正常

**建议：** 预运行至少进行1000步（1 ps），观察CV是否稳定。

### 7.4 从预运行到正式运行

预运行完成后：
1. 记录CV的波动范围（如 ±0.3 rad）
2. 设定 SCALE = 波动范围 × 0.5（如 0.15）
3. 根据CV波动幅度调整WW（通常 0.0005 - 0.002 Hartree）
4. 确认NT_HILLS（每步沉积频率）
5. 设置 DO_HILLS .TRUE.，正式启动元动力学

---

## 8. FES重建

### 8.1 HILLS文件格式

元动力学模拟会生成 `HILLS` 文件，包含每个沉积的高斯的信息：

```
# step  phi[rad]  psi[rad]  sigma_phi  sigma_psi  height  biasf
100     -1.234     0.567     0.1        0.1        0.001   1.0
150     -1.156     0.612     0.1        0.1        0.001   1.0
200     -1.089     0.589     0.1        0.1        0.001   1.0
...
```

各列含义：
1. 沉积时的MD步数
2. 第1个CV的值
3. 第2个CV的值
4. 第1个CV的高斯宽度
5. 第2个CV的高斯宽度
6. 高斯高度（对于Well-Tempered逐渐减小）
7. 偏差因子（对于Well-Tempered > 1）

### 8.2 使用graph工具重建FES

CP2K提供了 `graph` 工具（或 `graph.psmp`）来从HILLS文件重建自由能面：

```bash
graph.psmp -i HILLS -stride 10 -ndim 2 -ndw 1 2 -cp2k -integrated_fes
```

参数说明：
- `-i HILLS`：输入HILLS文件
- `-stride 10`：每10个高斯取一个采样点
- `-ndim 2`：CV空间维度（2D自由能面）
- `-ndw 1 2`：使用的CV编号（使用第1和第2个CV）
- `-cp2k`：指定CP2K格式的HILLS文件
- `-integrated_fes`：输出积分的自由能面

### 8.3 其他分析工具

也可以使用其他工具分析HILLS文件：

**Python + FESIL/plumed：**
```python
import numpy as np
# 读取HILLS文件
data = np.loadtxt('HILLS', comments='#')
# 手动重建FES（简化版）
```

**PLUMED sum_hills工具：**
如果将CP2K的HILLS文件转换为PLUMED格式，可以使用 `plumed sum_hills` 工具。

### 8.4 收敛性判断

**方法1：比较不同时刻的FES**
```bash
# 从模拟中段重建
graph.psmp -i HILLS -stride 10 -ndim 2 -ndw 1 2 -cp2k -integrated_fes \
  -skip 5000   # 跳过前5000个高斯

# 从模拟后段重建
graph.psmp -i HILLS -stride 10 -ndim 2 -ndw 1 2 -cp2k -integrated_fes
```

如果两次重建的FES差异小于 $k_BT$（约0.6 kcal/mol at 300K），则认为收敛。

**方法2：监控CV遍历**
- CV应该多次访问所有感兴趣的区域
- 不应持续停留在某个区域

**方法3：偏差势平坦度**
- 收敛后，偏差势应该覆盖所有最小值
- 不应有"孔洞"（未被填满的区域）

### 8.5 可视化FES

使用Python绘制自由能面等高线图：

```python
import numpy as np
import matplotlib.pyplot as plt

# 读取FES文件（由graph工具生成）
fes = np.loadtxt('fes.dat')
phi = fes[:, 0].reshape(nphi, npsi)    # 或根据输出格式调整
psi = fes[:, 1].reshape(nphi, npsi)
energy = fes[:, 2].reshape(nphi, npsi)

# 绘制等高线图
plt.figure(figsize=(8, 6))
levels = np.arange(0, 15, 0.5)   # 0-15 kcal/mol, 间隔0.5
cs = plt.contourf(phi*180/np.pi, psi*180/np.pi, energy*627.5, levels=levels, cmap='RdYlBu_r')
plt.colorbar(cs, label='Free Energy (kcal/mol)')
plt.xlabel(r'$\phi$ (degrees)')
plt.ylabel(r'$\psi$ (degrees)')
plt.title('Free Energy Surface of Alanine Dipeptide')
plt.savefig('fes_contour.png', dpi=300)
plt.show()
```

---

## 9. 丙氨酸二肽示例

### 9.1 为什么选择丙氨酸二肽？

丙氨酸二肽（Ala dipeptide，或称N-acetyl-N'-methylalaninamide，简称NANMA）是蛋白质构象研究的"氢原子"：

- **结构简单：** CH3-CO-NH-CH(CH3)-CO-NH-CH3
- **生物相关：** 蛋白质骨架的最小模型
- **构象丰富：** phi/psi Ramachandran空间的典型构象（αR, β, αL等）
- **文献丰富：** 大量参考数据可用于验证

### 9.2 丙氨酸二肽的构象

丙氨酸二肽的构象主要由两个骨架二面角决定：

**Phi角（φ）：** C(i-1)-N-Cα-C(i) 的二面角
**Psi角（ψ）：** N-Cα-C(i)-N(i+1) 的二面角

主要构象：
| 构象 | φ (degrees) | ψ (degrees) | 描述 |
|------|-------------|-------------|------|
| αR (alpha-right) | -60 | -40 | 右手α螺旋构象 |
| β (beta) | -120 | 130 | 伸展/β折叠构象 |
| αL (alpha-left) | 60 | 40 | 左手α螺旋构象 |
| γ | -80 | 80 | γ-turn构象 |

### 9.3 原子编号与CV定义

在CP2K输入文件中，丙氨酸二肽的原子编号（注意：CP2K从1开始计数）：

```
原子序列（典型XYZ文件）：
1:  C  (CH3)
2:  C  (羰基C)
3:  O  (羰基O)
4:  N  (肽键N)
5:  Cα (中心Cα)
6:  Hα (Cα上的H)
7:  Cβ (侧链CH3)
8:  ...
9:  C  (第二个羰基C)
...
15: N  (第二个肽键N)
...
17: C  (NHCH3的C)
```

Phi角的定义：原子 2-4-5-9（C-N-Cα-C）
Psi角的定义：原子 4-5-9-15（N-Cα-C-N）

**注意：** 原子编号取决于具体的XYZ文件顺序，实际使用时需要根据您的输入结构确认。

### 9.4 使用xTB还是DFT？

对于教学示例，使用xTB方法（GFN1-xTB）是最实际的选择：
- 计算速度快（比DFT快100-1000倍）
- 对有机分子的构象能有合理精度
- 足以展示元动力学的工作原理

对于正式研究：
- 推荐使用PBE/DZVP或B3LYP/TZVP等DFT方法
- 考虑溶剂效应（显式或隐式）

### 9.5 输入文件详解

输入文件 `alanine_dipeptide_meta.inp` 的关键部分将在下文详细解读。

---

## 10. 进阶阅读与常见问题

### 10.1 常见问题

**Q: 元动力学跑了很长时间，FES还是不平怎么办？**
A: 可能原因：
1. 模拟时间不够 → 增加步数
2. NT_HILLS太大 → 减小以加快偏差势积累
3. WW太小 → 增大高斯高度
4. CV选择不当 → 重新评估CV是否充分描述了自由度

**Q: 高斯势沉积到边界后怎么办？**
A: 添加反射壁，或者使用周期性CV（如二面角天然具有周期性）。

**Q: 可以同时运行多个元动力学副本吗？**
A: 可以。使用不同的初始速度或起始构型，然后平均多个FES以提高统计精度。

**Q: Well-Tempered和标准元动力学该选哪个？**
A: 优先选择Well-Tempered，因为它有更好的收敛性。标准元动力学主要用于初步探索。

**Q: CV之间需要正交吗？**
A: 不需要正交，但高相关性的CV会导致采样效率降低。选择物理上独立的CV更有效。

### 10.2 进阶主题

- **多Walker Metadynamics：** 多个独立的模拟共享同一个偏差势，加速探索
- **Bias-Exchange Metadynamics：** 不同的模拟使用不同的CV，定期交换构型
- **Metadynamics + REMD：** 结合平行回火和元动力学
- **Funnel Metadynamics：** 专门用于配体-蛋白质结合自由能计算
- **On-the-fly Probability Enhanced Sampling (OPES)：** 元动力学的新发展

### 10.3 参考文献

1. **Laio, A. & Parrinello, M.** (2002). Escaping free-energy minima. *PNAS*, 99(20), 12562-12566.
2. **Barducci, A., Bussi, G., & Parrinello, M.** (2008). Well-tempered metadynamics: A smoothly converging and tunable free-energy method. *Phys. Rev. Lett.*, 100(2), 020603.
3. **Laio, A. & Gervasio, F. L.** (2008). Metadynamics: A method to simulate rare events and reconstruct the free energy in biophysics, chemistry and material science. *Rep. Prog. Phys.*, 71(12), 126601.
4. **Valsson, O., Tiwary, P., & Parrinello, M.** (2016). Enhancing important fluctuations: Rare events and metadynamics from a conceptual viewpoint. *Annu. Rev. Phys. Chem.*, 67, 159-184.

### 10.4 后续章节预告

- **第11章：** GFN1-xTB半经验方法 — 更快速的力场计算，适合大系统和长时间MD
- **第12章：** QM/MM混合方法 — 将高精度量子化学与力场结合，处理复杂生物体系

---

## 附录A：关键词速查表

| 关键词 | 位置 | 说明 |
|--------|------|------|
| `DO_HILLS` | `METADYN` | 启用/禁用高斯沉积 |
| `NT_HILLS` | `METADYN` | 高斯沉积频率（步数） |
| `WW` | `METADYN` | 高斯高度（Hartree） |
| `DO_WELLTEMPERED` | `METADYN` | 启用Well-Tempered元动力学 |
| `BIAS_TEMPERATURE` | `METADYN` | 虚假温度（K） |
| `METAVAR` | `METADYN` | 集合变量的偏差设置 |
| `COLVAR` | `METAVAR` | 引用的COLVAR编号 |
| `SCALE` | `METAVAR` | 高斯宽度 |
| `&TORSION` | `COLVAR` | 二面角CV定义 |
| `&DISTANCE` | `COLVAR` | 距离CV定义 |
| `&COORDINATION` | `COLVAR` | 配位数CV定义 |
| `&ANGLE` | `COLVAR` | 角度CV定义 |
| `&WALL` | `METADYN` | 反射壁定义 |

## 附录B：推荐工作流程

```
1. 构建初始结构
   └── 使用GaussView, Avogadro, 或手动构建

2. 能量最小化
   └── GEO_OPT, 确保结构合理

3. NVT平衡化
   └── 5000-10000步, 热浴平衡

4. 元动力学预运行
   └── DO_HILLS .FALSE., 1000-5000步
   └── 分析CV波动范围

5. 设定元动力学参数
   └── 根据预运行设定WW, SCALE, NT_HILLS

6. 元动力学模拟
   └── DO_HILLS .TRUE., 足够长的模拟

7. FES重建
   └── graph.psmp工具, 检查收敛性

8. 后处理与分析
   └── 可视化FES, 提取构象, 计算势垒高度
```
