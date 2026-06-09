# 第二章：水分子的静态DFT计算

> 难度: ⭐ | 预计用时: 45 分钟

---

## 目录

1. [DFT基础回顾](#1-dft基础回顾)
2. [CP2K的DFT实现：QUICKSTEP](#2-cp2k的dft实现quickstep)
3. [关键概念](#3-关键概念)
4. [水分子输入文件详解](#4-水分子输入文件详解)
5. [输出解读](#5-输出解读)
6. [练习](#6-练习)

---

## 1. DFT基础回顾

### 1.1 什么是密度泛函理论

密度泛函理论 (Density Functional Theory, DFT) 是目前量子化学和凝聚态物理中应用最广泛的电子结构计算方法。它的核心思想是：**用电子密度 ρ(r) 代替波函数 Ψ 作为基本变量来描述多电子体系的基态性质**。

这一思想基于两个 Hohenberg-Kohn 定理（1964年）：

1. **定理一**：外势 v(r) 是基态电子密度 ρ(r) 的唯一泛函（即知道了 ρ(r)，就确定了体系的所有基态性质）
2. **定理二**：存在一个能量泛函 E[ρ]，其在真实基态密度处取极小值

### 1.2 Kohn-Sham 方程

1965年，Kohn 和 Sham 提出了一个实用的方案：将多电子问题映射到一组无相互作用的单电子方程——**Kohn-Sham 方程**：

```
[-½∇² + v_eff(r)] φᵢ(r) = εᵢ φᵢ(r)
```

其中有效势 v_eff 包含三部分：

```
v_eff(r) = v_ext(r) + v_H(r) + v_xc(r)
```

- **v_ext(r)**：外部势（通常是原子核产生的库仑势）
- **v_H(r)**：Hartree 势，即经典电子-电子库仑排斥
- **v_xc(r)**：交换关联势，包含所有量子力学的多体效应

电子密度由 Kohn-Sham 轨道构建：

```
ρ(r) = Σᵢ |φᵢ(r)|²     （对所有占据态求和）
```

### 1.3 交换关联泛函

DFT 的精度完全取决于**交换关联泛函** E_xc[ρ] 的近似形式。常用的泛函按精度递增分为几个层级（"Jacob's Ladder"）：

| 层级 | 类型 | 代表泛函 | 特点 |
|------|------|----------|------|
| 1 | LDA (局域密度近似) | **PADE**, VWN, PZ | 只依赖 ρ(r)，最简单 |
| 2 | GGA (广义梯度近似) | **PBE**, BLYP, PW91 | 依赖 ρ 和 ∇ρ |
| 3 | meta-GGA | TPSS, SCAN | 依赖 ρ、∇ρ 和动能密度 |
| 4 | 杂化泛函 (Hybrid) | **B3LYP**, PBE0, HSE06 | 混入部分 Hartree-Fock 交换 |
| 5 | 双杂化泛函 | B2PLYP | 混入 HF 交换和 MP2 关联 |

- **PADE**：LDA 泛函的一种参数化形式，计算速度快但精度有限，适合教学和初步测试
- **PBE**：GGA 泛函，是凝聚态物理中最常用的泛函，精度和效率的平衡较好
- **B3LYP**：杂化泛函，化学中最常用的泛函之一，但计算成本较高（需要计算 HF 交换）

> **CP2K 中的泛函命名**：在 CP2K 中，PADE 对应 LDA 的 Perdew-Zunger 参数化（也常简写为 LDA）。如果看到 `&XC_FUNCTIONAL PADE`，就相当于其他软件中的 LDA/VWN 或 LDA/PZ。

### 1.4 自洽场 (SCF) 迭代

由于 Kohn-Sham 方程中的有效势依赖于电子密度，而电子密度又由 Kohn-Sham 轨道决定，因此必须进行**自洽迭代 (SCF, Self-Consistent Field)**：

```
初始猜测 ρ₀(r)
    ↓
计算 v_eff[ρₙ](r)
    ↓
求解 Kohn-Sham 方程 → {φᵢ}
    ↓
计算新密度 ρₙ₊₁(r) = Σ|φᵢ|²
    ↓
判断收敛？ → 是 → 输出结果
    ↓ 否
混合 ρₙ 和 ρₙ₊₁ → 更新 ρₙ
    ↓
回到计算 v_eff
```

SCF 收敛的快慢取决于混合策略。常见方法：

- **线性混合**：ρ_new = (1-α)ρ_old + α·ρ_new，简单但收敛慢
- **Broyden 混合**：利用前几步的迭代信息，自适应地估计最优混合，收敛较快
- **Pulay (DIIS) 混合**：直接反演迭代子空间方法，广泛使用

---

## 2. CP2K的DFT实现：QUICKSTEP

### 2.1 GPW 方法概述

CP2K 的 DFT 模块名为 **QUICKSTEP**，它实现了一种称为 **GPW (Gaussian and Plane Waves)** 的混合基组方法。该方法由 Lippert, Hutter 和 Parrinello 于 1999 年提出。

GPW 方法的核心思想：

- **高斯型基函数 (GTO)**：用于展开 Kohn-Sham 轨道（分子轨道）。使用高斯型轨道 (Gaussian Type Orbitals) 的线性组合 (LCAO) 来表示轨道
- **平面波/实空间网格**：用于表示电子密度和进行各种实空间积分操作

这种混合方案的优势：

1. 高斯基组对原子核附近电子描述效率高
2. 实空间网格对 Hartree 势和交换关联势的计算方便高效
3. 可以使用高效的 FFT 和多重网格技术
4. 计算标度可以做得很好，适合大体系

### 2.2 计算流程

QUICKSTEP 的一次 SCF 迭代流程如下：

```
1. 电子密度 ρ(r) 在高斯基组中展开（系数矩阵 C）
       ↓
2. 将高斯密度映射 (collocation) 到实空间网格
       ↓
3. 在实空间网格上计算：
   - Hartree 势（通过求解 Poisson 方程，使用多重网格法）
   - XC 势（通过数值积分）
       ↓
4. 将网格上的势能映射回高斯基组 (integration)
       ↓
5. 构建 Kohn-Sham 矩阵 H
       ↓
6. 求解广义特征值问题 HC = SCE → 获得新的轨道系数 C
       ↓
7. 用新的 C 计算新的密度 → 检查 SCF 收敛
```

### 2.3 GPW 与 GAPW

CP2K 提供两种 Quickstep 方法：

| 方法 | 全称 | 特点 |
|------|------|------|
| **GPW** | Gaussian and Plane Waves | 赝势 + 平面波，效率高，适合大多数情况 |
| **GAPW** | Gaussian Augmented Plane Waves | 全电子计算，区分核心/价电子区域 |

本教程所有示例均使用 GPW 方法（默认）。

---

## 3. 关键概念

### 3.1 基组 (Basis Sets)

在 DFT 计算中，分子轨道 φᵢ 用一组预定义的**基函数**的线性组合来展开：

```
φᵢ(r) = Σ_μ c_μᵢ χ_μ(r)
```

其中 χ_μ(r) 是基函数，c_μᵢ 是展开系数。

CP2K 使用**高斯型轨道 (GTO)** 作为基函数。GTO 的优点是多中心积分可以解析计算，速度快。

#### 基组命名约定

CP2K 的基组名称遵循以下格式：

```
精度等级-GTH-泛函名
```

| 前缀 | 含义 | 说明 |
|------|------|------|
| SZV | Single Zeta Valence | 最小基组，每个价轨道一个基函数 |
| DZVP | Double Zeta Valence Polarized | 双 zeta + 极化函数，性价比高 |
| TZV2P | Triple Zeta Valence 2-Polarized | 三 zeta + 两个极化函数，精度高 |
| TZVP | Triple Zeta Valence Polarized | 三 zeta + 一个极化函数 |
| MOLOPT | 优化的分子基组 | 专门为分子计算优化 |

**"zeta"的含义**：类比于量子力学中氢原子轨道的描述。单 zeta = 每个价轨道用一个基函数描述；双 zeta = 用两个不同宽度的基函数描述（允许轨道"呼吸"），以此类推。

**极化函数**：比价轨道更高角动量的基函数。例如，对氧原子，价轨道是 2s 和 2p，极化函数就是 d 轨道。极化函数允许轨道变形，对化学键和分子极性的描述很重要。

#### 选择基组的经验法则

| 用途 | 推荐基组 |
|------|----------|
| 初步测试/教学 | SZV 或 DZVP |
| 日常计算/分子动力学 | DZVP |
| 高精度单点能 | TZV2P 或更大 |
| 气相分子 | MOLOPT 基组 |
| 周期性体系 | DZVP-MOLOPT-SR（短程优化） |

### 3.2 赝势 (Pseudopotentials)

#### 什么是赝势

赝势方法将原子中**核心电子**的效应用一个有效势（赝势）替代，只显式计算**价电子**。

为什么要这样做？

1. **减少电子数目**：例如硅 (Si, Z=14) 有 14 个电子，但核心电子有 10 个（1s²2s²2p⁶），只需显式计算 4 个价电子
2. **消除核心态的振荡**：核心电子的波函数在原子核附近快速振荡，需要非常细的网格才能描述，赝势使波函数变得光滑
3. **提高计算效率**：电子数减少 → 基组更小 → 计算更快

#### GTH 赝势

CP2K 使用 **GTH (Goedecker-Teter-Hutter) 赝势**，这是一种解析的模守恒赝势 (norm-conserving pseudopotential)。

GTH 赝势的关键特点：

- 参数化形式，有解析公式
- 与 LDA/GGA 等特定泛函匹配
- 在 CP2K 数据库中有预定义的完整参数集

#### 赝势命名

```
GTH-泛函名-qN
```

- `GTH`：Goedecker-Teter-Hutter 赝势
- `泛函名`：如 `PADE`（对应 LDA）、`PBE`（对应 GGA-PBE）
- `qN`：赝势处理的价电子数

**常见元素的赝势选择**：

| 元素 | 原子序数 | 电子构型 | GTH-PADE-qN | 说明 |
|------|---------|----------|-------------|------|
| H | 1 | 1s¹ | q1 | 1个价电子 |
| C | 6 | [He]2s²2p² | q4 | 4个价电子 |
| N | 7 | [He]2s²2p³ | q5 | 5个价电子 |
| O | 8 | [He]2s²2p⁴ | q6 | 6个价电子 |
| Si | 14 | [Ne]3s²3p² | q4 | 4个价电子 |
| He | 2 | 1s² | q2 | 2个价电子 |

### 3.3 交换关联泛函的选择

在 CP2K 中指定泛函：

```fortran
&XC
  &XC_FUNCTIONAL PADE       ! LDA 泛函
  &END XC_FUNCTIONAL
&END XC
```

常用泛函及其在 CP2K 中的名称：

| 泛函 | CP2K 关键字 | 类型 | 说明 |
|------|-------------|------|------|
| LDA (Perdew-Zunger) | PADE | LDA | 最简单的泛函 |
| PBE | PBE | GGA | 凝聚态常用 |
| BLYP | BLYP | GGA | 量子化学常用 |
| PBE0 | PBE0 | Hybrid | 含 25% HF 交换 |
| B3LYP | B3LYP | Hybrid | 化学最常用杂化泛函 |

> **重要**：泛函选择必须与赝势和基组匹配！使用 PADE 泛函时用 `GTH-PADE` 系列赝势，使用 PBE 泛函时用 `GTH-PBE` 系列赝势。

### 3.4 SCF 收敛控制

```fortran
&SCF
  SCF_GUESS ATOMIC          ! 初始猜测：原子密度叠加
  EPS_SCF 1.0E-7            ! 收敛判据（能量变化）
  MAX_SCF 50                ! 最大迭代步数
  &DIAGONALIZATION
    ALGORITHM STANDARD       ! 标准对角化
  &END DIAGONALIZATION
  &MIXING
    METHOD BROYDEN_MIXING    ! Broyden 混合
    ALPHA 0.4                ! 混合参数
    NBROYDEN 8               ! 历史步数
  &END MIXING
&END SCF
```

**SCF 收敛失败的常见原因及解决方法**：

| 问题 | 解决方法 |
|------|----------|
| 混合参数太大 | 减小 `ALPHA`（如 0.1~0.2） |
| 体系有金属特性 | 使用更保守的混合或 smearing |
| 初始猜测不好 | 尝试 `SCF_GUESS RESTART` |
| 基组/赝势不匹配 | 检查基组和赝势是否与泛函对应 |
| 网格太粗 | 增加 `CUTOFF` |

---

## 4. 水分子输入文件详解

### 4.1 完整输入文件

以下是 `H2O_static.inp` 的完整内容：

```fortran
&GLOBAL
  PROJECT H2O_static
  RUN_TYPE ENERGY_FORCE
  PRINT_LEVEL LOW
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &SUBSYS
    &KIND O
      ELEMENT O
      BASIS_SET DZVP-GTH-PADE
      POTENTIAL GTH-PADE-q6
    &END KIND
    &KIND H
      ELEMENT H
      BASIS_SET DZVP-GTH-PADE
      POTENTIAL GTH-PADE-q1
    &END KIND
    &CELL
      ABC 12.4 12.4 12.4
    &END CELL
    &COORD
      O   0.000000   0.000000   0.117310
      H   0.000000   0.756950  -0.469241
      H   0.000000  -0.756950  -0.469241
    &END COORD
  &END SUBSYS
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
      EPS_SCF 1.0E-7
      MAX_SCF 50
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
    &PRINT
      &FORCES ON
      &END FORCES
    &END PRINT
  &END DFT
&END FORCE_EVAL
```

### 4.2 逐段解析

#### 全局设置

```fortran
&GLOBAL
  PROJECT H2O_static
  RUN_TYPE ENERGY_FORCE
  PRINT_LEVEL LOW
&END GLOBAL
```

- `PROJECT H2O_static`：输出文件前缀
- **`RUN_TYPE ENERGY_FORCE`**：计算单点能**和**原子受力。这与第一章的 `ENERGY` 不同——`ENERGY_FORCE` 会额外输出每个原子上的力。对于几何优化和分子动力学，力是必需的。
- `PRINT_LEVEL LOW`：简洁输出

#### 体系定义

```fortran
  &SUBSYS
    &KIND O
      ELEMENT O
      BASIS_SET DZVP-GTH-PADE
      POTENTIAL GTH-PADE-q6
    &END KIND
    &KIND H
      ELEMENT H
      BASIS_SET DZVP-GTH-PADE
      POTENTIAL GTH-PADE-q1
    &END KIND
```

水分子中有两种元素（O 和 H），因此需要**两个 `&KIND` 节段**：

- **氧 (O)**：价电子数 6（2s²2p⁴），所以用 `GTH-PADE-q6`
- **氢 (H)**：价电子数 1（1s¹），所以用 `GTH-PADE-q1`
- 两者都使用 `DZVP-GTH-PADE` 基组

```fortran
    &CELL
      ABC 12.4 12.4 12.4
    &END CELL
```

立方盒子，边长 12.4 Å。对于孤立水分子，这个尺寸足够大。

> **经验法则**：盒子边长应至少比分子最大尺寸大 6 Å 以上，以避免周期性镜像效应。水分子的 O-H 键长约 0.96 Å，分子直径约 2.5 Å，12.4 Å 的盒子留出了约 10 Å 的真空空间。

```fortran
    &COORD
      O   0.000000   0.000000   0.117310
      H   0.000000   0.756950  -0.469241
      H   0.000000  -0.756950  -0.469241
    &END COORD
  &END SUBSYS
```

水分子的实验几何构型（实验平衡键长 O-H = 0.9572 Å，键角 H-O-H = 104.52°）。

从坐标可以计算：
- O-H 距离 = sqrt(0.756950² + (0.117310 + 0.469241)²) ≈ 0.9572 Å
- H-O-H 角度 ≈ 104.52°

#### DFT 参数

与第一章类似，但有几处新增：

```fortran
    &XC
      &XC_FUNCTIONAL PADE
      &END XC_FUNCTIONAL
    &END XC
    &PRINT
      &FORCES ON
      &END FORCES
    &END PRINT
```

`&PRINT` 子节段控制额外的输出。`&FORCES ON` 表示在输出中打印原子受力信息。这对于 `ENERGY_FORCE` 运行类型来说很有用，可以确认力的数值是否合理（例如，接近平衡几何时力应接近零）。

---

## 5. 输出解读

### 5.1 运行计算

```bash
cp2k.popt -o H2O_static.out H2O_static.inp
```

### 5.2 关键输出内容

#### 总能量

```
 ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):           -17.18927716540000
```

水分子在 PADE/LDA 泛函 + DZVP 基组下的总能量约为 **-17.189 Ha**。

> **注意**：总能量的绝对值没有直接的物理意义（它取决于参考能量的定义）。有意义的是能量差，如反应能、结合能等。

#### 电子计数

```
 Number of electrons:                                                  8
 Number of occupied Orbitals:                                          4
```

水分子有 8 个价电子（O: 6 个，2×H: 各 1 个），占据 4 个分子轨道。这是因为赝势只处理价电子（O 的核心电子已被赝势替代）。

#### 力的输出

由于我们使用了 `ENERGY_FORCE` 和 `&FORCES ON`，输出中会包含每个原子上的力：

```
 ATOMIC FORCES in [a.u.]

 # Atom   Kind   Element          X              Y              Z
      1      1     O           0.00000000     0.00000000    -0.00000000
      2      2     H          -0.00000000     0.00000000     0.00000000
      3      2     H           0.00000000    -0.00000000     0.00000000

 # Total force:     0.00000000
```

因为我们使用的是实验平衡几何构型，所以力应该非常小（接近零）。如果使用的是非平衡几何，力会较大，可以作为几何优化的驱动力。

#### SCF 收敛信息

观察 SCF 迭代过程，确认收敛情况。典型的水分子 SCF 迭代步数为 3-10 步。

### 5.3 计算结果摘要

| 物理量 | 值 |
|--------|-----|
| 总能量 | ≈ -17.189 Ha |
| 价电子数 | 8 |
| 占据轨道数 | 4 |
| SCF 收敛步数 | 约 3-5 步 |
| 原子受力 | ≈ 0（平衡几何）|

---

## 6. 练习

### 练习 1：更换交换关联泛函

将 `PADE` 替换为 `PBE`，同时将赝势和基组也换成 PBE 系列：

```fortran
&KIND O
  BASIS_SET DZVP-GTH-PBE
  POTENTIAL GTH-PBE-q6
&END KIND
&KIND H
  BASIS_SET DZVP-GTH-PBE
  POTENTIAL GTH-PBE-q1
&END KIND
```

以及：

```fortran
&XC_FUNCTIONAL PBE
```

比较 PADE 和 PBE 的总能量。哪一个能量更低？

### 练习 2：添加打印输出

在 `&DFT` 节段内添加更多打印选项：

```fortran
&PRINT
  &FORCES ON
  &END FORCES
  &DENSITY_CUBE
    STRIDE 1 1 1
  &END DENSITY_CUBE
  &MO_CUBES
    NHOMO 4
    NLUMO 4
  &END MO_CUBES
&END PRINT
```

- `DENSITY_CUBE` 会输出电子密度的 Cube 文件，可用 VMD 或 Avogadro 可视化
- `MO_CUBES` 会输出分子轨道的 Cube 文件（4个最高占据和 4个最低未占据轨道）

> **提示**：Cube 文件可以用 VMD 软件打开并可视化：`vmd *.cube`

### 练习 3：改变盒子大小

将晶胞大小从 12.4 Å 改为 8.0 Å 和 20.0 Å，分别计算能量。观察盒子大小对总能量的影响有多大？

预期结果：盒子足够大时（> 10 Å），能量变化应该很小（< 1 meV）。

### 练习 4：尝试不同基组

分别使用 `SZV-GTH-PADE` 和 `TZV2P-GTH-PADE` 替代 `DZVP-GTH-PADE`，计算水分子能量。观察基组大小对总能量和 SCF 收敛行为的影响。

### 思考题

1. 水分子有 10 个全电子（O: 8, 2×H: 各 1），但计算中只处理了 8 个价电子。缺失的 2 个电子在哪里？
2. 为什么 CP2K 需要同时指定基组文件 (`BASIS_SET`) 和赝势文件 (`GTH_POTENTIALS`)？它们各自的作用是什么？
3. 如果 SCF 不收敛，可以尝试调整哪些参数？

---

## 参考资源

- Kohn, W.; Sham, L. J. *Phys. Rev.* **140**, A1133 (1965) — Kohn-Sham 方程原始论文
- Hohenberg, P.; Kohn, W. *Phys. Rev.* **136**, B864 (1964) — Hohenberg-Kohn 定理
- Perdew, J. P.; Zunger, A. *Phys. Rev. B* **23**, 5048 (1981) — LDA 参数化
- Perdew, J. P.; Burke, K.; Ernzerhof, M. *Phys. Rev. Lett.* **77**, 3865 (1996) — PBE 泛函
- Goedecker, S.; Teter, M.; Hutter, J. *Phys. Rev. B* **54**, 1703 (1996) — GTH 赝势
