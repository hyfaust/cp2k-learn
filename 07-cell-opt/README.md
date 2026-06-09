# 第七章：晶胞优化

**难度：⭐⭐⭐**

在第六章中，我们学习了几何优化——在固定的晶胞内优化原子位置。但对于周期性体系，晶胞参数本身也是变量。**晶胞优化 (Cell Optimization)** 同时优化原子位置和晶胞参数，找到能量最低的晶体结构。

---

## 目录

1. [晶胞优化 vs 几何优化](#1-晶胞优化-vs-几何优化)
2. [应力张量](#2-应力张量)
3. [CELL_OPT 设置](#3-cell_opt-设置)
4. [优化策略](#4-优化策略)
5. [Si 晶体示例](#5-si-晶体示例)
6. [收敛判据](#6-收敛判据)
7. [输入文件详解](#7-输入文件详解)
8. [练习](#8-练习)

---

## 1. 晶胞优化 vs 几何优化

### 1.1 关键区别

| 特征 | GEO_OPT | CELL_OPT |
|------|---------|----------|
| 优化变量 | 仅原子位置 $\{\mathbf{R}_i\}$ | 原子位置 $\{\mathbf{R}_i\}$ + 晶胞参数 $\{a, b, c, \alpha, \beta, \gamma\}$ |
| 晶胞 | 固定 | 可变 |
| 输出 | 优化后的坐标 | 优化后的坐标 + 晶胞参数 |
| 适用场景 | 分子、固定晶胞内的弛豫 | 确定最优晶格常数、结构预测 |

### 1.2 何时需要晶胞优化？

- **确定材料的平衡晶格常数**：实验值可能存在温度效应，计算给出 0 K 理论值
- **从近似结构出发**：初始晶格常数不精确，需要自动优化
- **研究压力效应**：在不同外部压力下优化结构
- **结构预测**：探索不同组分的合金结构
- **为后续计算准备结构**：MD、声子计算等需要精确的平衡结构

### 1.3 计算成本

晶胞优化比几何优化**昂贵得多**，因为：

1. 每步需要计算**应力张量**（比只算力更耗时）
2. 优化变量更多（6 个晶胞参数 + 3N 个原子坐标）
3. 晶胞变化导致基组和网格变化

---

## 2. 应力张量

### 2.1 什么是应力张量？

**应力张量 (stress tensor)** $\sigma_{\alpha\beta}$ 描述材料内部的应力状态：

$$\sigma = \begin{pmatrix} \sigma_{xx} & \sigma_{xy} & \sigma_{xz} \\ \sigma_{yx} & \sigma_{yy} & \sigma_{yz} \\ \sigma_{zx} & \sigma_{zy} & \sigma_{zz} \end{pmatrix}$$

- 对角元素 $\sigma_{xx}, \sigma_{yy}, \sigma_{zz}$：法向应力（压/拉）
- 非对角元素 $\sigma_{xy}$ 等：剪切应力
- 应力张量是对称的（$\sigma_{\alpha\beta} = \sigma_{\beta\alpha}$），因此有 6 个独立分量

### 2.2 Virial 定理

在 DFT 中，应力张量通过 **virial 定理** 计算：

$$\sigma_{\alpha\beta} = \frac{1}{V}\left(\sum_i m_i v_{i\alpha} v_{i\beta} + \sum_i R_{i\alpha} F_{i\beta}\right)$$

对于静态计算（$T=0$），动能项为零：

$$\sigma_{\alpha\beta} = \frac{1}{V}\sum_i R_{i\alpha} F_{i\beta}$$

在 CP2K 中，还需要加上电子部分的贡献。

### 2.3 压力

**静水压力 (hydrostatic pressure)** 是应力张量的迹：

$$P = -\frac{1}{3}\text{Tr}(\sigma) = -\frac{1}{3}(\sigma_{xx} + \sigma_{yy} + \sigma_{zz})$$

平衡态下，内部压力应为零（或等于外部施加的压力）。

### 2.4 压力单位

| 单位 | 换算关系 |
|------|---------|
| 1 GPa | = 10,000 bar |
| 1 GPa | = 10⁹ Pa |
| 1 Ha/Bohr³ | ≈ 29421.03 GPa |

CP2K 输出中通常使用 GPa 或 bar。

---

## 3. CELL_OPT 设置

### 3.1 基本设置

```
&GLOBAL
  PROJECT  Si_cell_opt
  RUN_TYPE CELL_OPT       ! ★ 晶胞优化
&END GLOBAL
```

### 3.2 &CELL_OPT 部分

在 `&MOTION` 中设置晶胞优化参数：

```
&MOTION
  &CELL_OPT
    OPTIMIZER  BFGS            ! 优化算法
    TYPE  DIRECT_CELL_OPT      ! 直接优化晶胞参数
    KEEP_ANGLES  .FALSE.       ! 是否保持晶胞角度不变
    PRESSURE_TOLERANCE  0.1    ! 压力收敛判据 (GPa)
    MAX_ITER  200              ! 最大步数
  &END CELL_OPT
&END MOTION
```

### 3.3 TYPE 关键字

| 值 | 说明 |
|-----|------|
| `DIRECT_CELL_OPT` | 直接优化晶胞参数（推荐）|
| `CG_CELL_OPT` | 使用 CG 算法优化晶胞 |
| `LBFGS_CELL_OPT` | 使用 L-BFGS 算法优化晶胞 |

> **注意：** `OPTIMIZER` 和 `TYPE` 是配合使用的。推荐使用 `DIRECT_CELL_OPT`。

### 3.4 KEEP_ANGLES

```
KEEP_ANGLES  .TRUE.    ! 保持 α=β=γ=90° 不变（仅优化 a, b, c）
KEEP_ANGLES  .FALSE.   ! 允许角度变化（优化全部 6 个参数）
```

- 对于立方晶体，设为 `.TRUE.` 可以节省计算
- 对于低对称性晶体，应设为 `.FALSE.`

### 3.5 CONSTRAINT（体积/形状约束）

可以施加晶胞约束：

```
&CELL_OPT
  CONSTRAINT  VOLUME    ! 固定体积，只优化形状
  ! CONSTRAINT  NONE    ! 无约束（默认）
  ! CONSTRAINT  CELL_SHAPE  固定形状，只优化体积
&END CELL_OPT
```

---

## 4. 优化策略

### 4.1 交替优化策略

一种高效的做法是交替进行几何优化和晶胞优化：

```
初始结构
  │
  ├→ GEO_OPT（固定晶胞，优化原子位置）
  │
  ├→ CELL_OPT（优化晶胞 + 原子位置）
  │
  └→ 最终平衡结构
```

但 CP2K 的 CELL_OPT 已经包含了原子位置的优化，通常不需要手动交替。

### 4.2 各向同性 vs 各向异性

**各向同性 (Isotropic) 优化：**
- 晶胞保持立方形状（a = b = c, α = β = γ = 90°）
- 只有一个自由度（晶格常数 a）
- 适用于立方晶体的简单优化

**各向异性 (Anisotropic) 优化：**
- 晶胞的三个边长和三个角度都可以独立变化
- 适用于低对称性晶体

在 CP2K 中通过 `KEEP_ANGLES` 和外部压力的对称性来控制。

### 4.3 外部压力

可以在晶胞优化中施加外部压力：

```
&CELL_OPT
  EXTERNAL_PRESSURE  1000.0  ! 外部压力 (bar)
  PRESSURE_TOLERANCE  0.1    ! 压力收敛判据 (GPa)
&END CELL_OPT
```

- `EXTERNAL_PRESSURE` = 0：零压平衡（默认）
- `EXTERNAL_PRESSURE` > 0：模拟高压条件
- 也可以用 3×3 矩阵指定各向异性压力

---

## 5. Si 晶体示例

### 5.1 问题设置

我们将从一个**略微偏离平衡值**的晶格常数出发，让 CP2K 优化到正确的晶格常数。

- 实验值：$a_{\text{Si}} = 5.4306975$ Å
- 初始猜测：$a_0 = 5.30$ Å（偏小约 2.4%）

优化过程将：
1. 调整晶格常数（增大，以降低应变能）
2. 在每个晶格常数下优化原子位置
3. 直到应力和力都满足收敛判据

### 5.2 预期结果

- 最终晶格常数应接近 5.43 Å（与实验值的偏差取决于 DFT 方法）
- LDA (PADE) 通常**低估**晶格常数约 1-2%
- GGA (PBE) 通常**高估**晶格常数约 1-2%

### 5.3 验证结果

优化完成后，检查输出中的：

1. **最终晶格参数**：
   ```
   CELL| Volume [angstrom^3]:  xxx.xx
   CELL| Vector a [angstrom]:   5.xxx  0.000  0.000
   CELL| Vector b [angstrom]:   0.000  5.xxx  0.000
   CELL| Vector c [angstrom]:   0.000  0.000  5.xxx
   ```

2. **压力/应力**：
   ```
   PRESSURE| Pressure [GPa]:      0.000xxx
   ```
   压力应接近零（在收敛判据范围内）。

3. **总能量**：应低于初始结构的能量。

---

## 6. 收敛判据

### 6.1 压力收敛

晶胞优化的主要收敛判据是**内部压力与外部压力之差**：

```
PRESSURE_TOLERANCE  0.1    ! 压力差 < 0.1 GPa 时收敛
```

### 6.2 原子位置收敛

同时，原子力也必须满足与几何优化相同的判据：

```
MAX_DR    3.0E-3     ! 最大位移 (Bohr)
RMS_DR    1.5E-3     ! RMS 位移 (Bohr)
MAX_FORCE  4.5E-4    ! 最大力 (Ha/Bohr)
RMS_FORCE  3.0E-4    ! RMS 力 (Ha/Bohr)
```

### 6.3 晶胞变化收敛

晶胞参数的变化也需收敛：

```
MAX_DR    3.0E-3     ! 晶胞矢量的最大变化
RMS_DR    1.5E-3     ! 晶胞矢量变化的 RMS
```

### 6.4 实用建议

| 计算目的 | PRESSURE_TOLERANCE | MAX_FORCE |
|---------|-------------------|-----------|
| 初步优化 | 0.5 GPa | 1E-3 Ha/Bohr |
| 标准优化 | 0.1 GPa | 4.5E-4 Ha/Bohr |
| 高精度优化 | 0.01 GPa | 1E-4 Ha/Bohr |

---

## 7. 输入文件详解

### 7.1 Si_cell_opt.inp — Si 晶体晶胞优化

完整输入文件分析：

```bash
# 运行命令
cp2k.psmp -i Si_cell_opt.inp -o Si_cell_opt.out
```

#### GLOBAL 部分

```
&GLOBAL
  PROJECT  Si_cell_opt
  RUN_TYPE CELL_OPT       ! 晶胞优化
  PRINT_LEVEL LOW
&END GLOBAL
```

- `RUN_TYPE CELL_OPT`：CP2K 将同时优化晶胞参数和原子位置

#### FORCE_EVAL 部分

```
&FORCE_EVAL
  METHOD  QS
  STRESS_TENSOR  ANALYTICAL   ! ★ 使用解析应力张量（更精确）

  &DFT
    BASIS_SET_FILE_NAME  BASIS_SET
    POTENTIAL_FILE_NAME  GTH_POTENTIALS

    &MGRID
      CUTOFF  300
      REL_CUTOFF  60
    &END MGRID

    &QS
      EPS_DEFAULT  1.0E-10
    &END QS

    &MIXING
      METHOD  BROYDEN_MIXING
      ALPHA  0.4
      NBUFFER  7
    &END MIXING

    &SCF
      SCF_GUESS  ATOMIC
      EPS_SCF  1.0E-6
      MAX_SCF  50
    &END SCF

    &XC
      &XC_FUNCTIONAL PADE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT

  &SUBSYS
    &CELL
      ABC  5.30  5.30  5.30     ! ★ 故意偏小的初始晶格常数
      PERIODIC  XYZ
    &END CELL

    &TOPOLOGY
      COORD_FILE_FORMAT  XYZ
    &END TOPOLOGY

    &KIND Si
      BASIS_SET  DZVP-GTH-PADE
      POTENTIAL  GTH-PADE-q4
    &END KIND

    &COORD
      ! 使用分数坐标 —— 更适合晶胞优化
      SCALED
      Si   0.000   0.000   0.000
      Si   0.250   0.250   0.250
      Si   0.500   0.500   0.000
      Si   0.750   0.750   0.250
      Si   0.500   0.000   0.500
      Si   0.750   0.250   0.750
      Si   0.000   0.500   0.500
      Si   0.250   0.750   0.750
    &END COORD
  &END SUBSYS
&END FORCE_EVAL
```

**关键点：**

1. **`STRESS_TENSOR ANALYTICAL`**：晶胞优化**必须**计算应力张量。`ANALYTICAL` 比数值方法更精确。

2. **初始晶格常数 5.30 Å**：故意偏离平衡值 5.43 Å，让优化过程将其调整到正确值。

3. **分数坐标 (`SCALED`)**：在晶胞优化中使用分数坐标非常重要！
   - 使用 Cartesian 坐标时，晶胞变化会导致原子位置相对改变
   - 使用分数坐标时，原子相对于晶胞的位置保持正确

#### MOTION 部分

```
&MOTION
  &CELL_OPT
    OPTIMIZER  BFGS              ! BFGS 优化器
    TYPE  DIRECT_CELL_OPT        ! 直接优化晶胞
    KEEP_ANGLES  .TRUE.          ! 保持立方对称性（α=β=γ=90°）
    PRESSURE_TOLERANCE  0.1      ! 压力收敛 0.1 GPa
    MAX_ITER  200                ! 最大步数
  &END CELL_OPT

  &GEO_OPT                        ! 同时设置几何优化参数
    OPTIMIZER  BFGS
    MAX_DR    3.0E-3
    RMS_DR    1.5E-3
    MAX_FORCE  4.5E-4
    RMS_FORCE  3.0E-4
  &END GEO_OPT
&END MOTION
```

### 7.2 输出解读

#### 每步输出

```
-------------------------------- Cell Optimization --------------------------------
  Step     Update method     Time        Change           Pressure      Volume
-----------------------------------------------------------------------------------
     1 BFGS                  2.5        0.12345         -15.234      148.87
     2 BFGS                  2.3        0.04567          -3.456      155.23
     3 BFGS                  2.4        0.01234          -0.567      158.91
     4 BFGS                  2.2        0.00345          -0.034      160.12
```

- `Change`：晶胞参数的变化量
- `Pressure`：内部压力（GPa），应趋向零
- `Volume`：晶胞体积（Å³）

#### 最终结果

```
CELL| Volume [angstrom^3]:     160.159
CELL| Vector a [angstrom    5.4307     0.0000     0.0000
CELL| Vector b [angstrom    0.0000     5.4307     0.0000
CELL| Vector c [angstrom    0.0000     0.0000     5.4307

ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):       -31.36850849423

PRESSURE| Pressure [GPa]:     0.000123
```

---

## 8. 练习

### 练习 1：从不同初始值出发

分别从 $a_0 = 5.0, 5.2, 5.4, 5.6, 5.8$ Å 出发进行晶胞优化，比较：
- 最终晶格常数是否一致
- 收敛步数
- 最终能量

### 练习 2：各向异性优化

修改输入文件，设置 `KEEP_ANGLES .FALSE.`，使用非立方初始晶胞（如 ABC = 5.3 5.4 5.5）。观察：
- 晶胞最终是否恢复立方形状
- 如果不恢复，说明什么？

### 练习 3：压力的影响

在不同外部压力（0, 5, 10, 20 GPa）下优化 Si 晶体：
- 晶格常数如何随压力变化
- 绘制 $a(P)$ 曲线

### 练习 4：体积与能量的关系

从一系列不同的晶格常数出发，只做单点能量计算（`ENERGY_FORCE`），绘制 $E(a)$ 曲线。与晶胞优化的结果对比。

### 练习 5：分数坐标 vs 笛卡尔坐标

分别使用 SCALED 和 Cartesian 坐标进行晶胞优化，观察结果是否一致。

---

## 小结

本章学习了晶胞优化的核心概念和 CP2K 实现：

| 概念 | CP2K 实现 |
|------|----------|
| 晶胞优化 | `RUN_TYPE CELL_OPT` |
| 应力张量 | `STRESS_TENSOR ANALYTICAL` |
| 优化算法 | `&CELL_OPT` OPTIMIZER BFGS |
| 晶胞约束 | `KEEP_ANGLES`, `CONSTRAINT` |
| 压力收敛 | `PRESSURE_TOLERANCE` |
| 坐标选择 | 推荐使用 SCALED 分数坐标 |
| 外部压力 | `EXTERNAL_PRESSURE` |

下一章将介绍**分子动力学**——如何模拟原子在有限温度下的运动。

---

**参考资源：**
- CP2K 手册 - &CELL_OPT：https://manual.cp2k.org/trunk/CP2K_INPUT/MOTION/CELL_OPT.html
- CP2K 手册 - STRESS_TENSOR：https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/STRESS_TENSOR.html
