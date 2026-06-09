# 第一章：你的第一个CP2K计算

> 难度: ⭐ | 预计用时: 30 分钟

---

## 目录

1. [什么是CP2K](#1-什么是cp2k)
2. [CP2K的输入文件结构](#2-cp2k的输入文件结构)
3. [关键概念解析](#3-关键概念解析)
4. [示例详解：helium.inp](#4-示例详解heliuminp)
5. [如何运行CP2K](#5-如何运行cp2k)
6. [输出文件解读](#6-输出文件解读)
7. [练习](#7-练习)

---

## 1. 什么是CP2K

### 1.1 简介

CP2K 是一款功能强大的开源量子化学与固态物理模拟软件包，专门用于原子尺度的模拟计算。它最初由苏黎世联邦理工学院 (ETH Zürich) 的 Jürg Hutter 教授课题组开发，目前由一个国际化的开发者社区共同维护。

CP2K 的核心能力包括：

- **密度泛函理论 (DFT)**：通过 QUICKSTEP 模块实现高效的 GPW/GAPW 方法
- **Møller-Plesset 微扰理论 (MP2)** 及其变体 (RI-MP2, RPA 等)
- **半经验量子力学方法**：如 DFTB（密度泛函紧束缚）
- **经典力场 (Force Field)**：通过 FIST 模块实现分子动力学
- **QM/MM 混合方法**：量子力学/分子力学耦合
- **从头算分子动力学 (AIMD)**：Car-Parrinello 和 Born-Oppenheimer 方案
- **元动力学 (Metadynamics)**：增强采样方法

### 1.2 技术特点

CP2K 的技术架构有几个重要特点：

| 特性 | 说明 |
|------|------|
| 编程语言 | Fortran 2008 |
| 并行方式 | MPI + OpenMP 混合并行，支持 CUDA/GPU 加速 |
| 周期性体系 | 支持三维周期性（晶体）、二维（表面）、一维（链）、零维（分子） |
| 电子结构方法 | GPW（高斯和平面波）混合基组方法 |
| 可执行文件 | `cp2k.ssmp`（单节点OpenMP）、`cp2k.psmp`（MPI+OpenMP）、`cp2k.popt`（纯MPI）等 |

### 1.3 适用场景

CP2K 特别适合以下类型的计算：

- **大体系的从头算分子动力学**：得益于 GPW 方法的线性标度特性
- **凝聚态体系**：晶体、液体、界面、表面等周期性体系
- **生物分子模拟**：蛋白质、DNA 等大分子的 QM/MM 模拟
- **材料科学**：催化剂、电池材料、半导体等

---

## 2. CP2K的输入文件结构

CP2K 使用一种特定的文本格式作为输入文件（通常以 `.inp` 或 `.cp2k` 为后缀）。理解这种格式是使用 CP2K 的第一步。

### 2.1 基本语法

CP2K 的输入文件使用 **节段（Section）** 结构，由 `&` 关键字标识：

```
&SECTION_NAME
   keyword1  value1
   keyword2  value2

   &SUBSECTION
      keyword3  value3
   &END SUBSECTION

&END SECTION_NAME
```

**语法规则**：

1. **节段开始**：以 `&` 开头，后跟节段名称（如 `&GLOBAL`）
2. **节段结束**：以 `&END` 开头，后跟节段名称（如 `&END GLOBAL`）
3. **关键字赋值**：在节段内部，使用 `关键字  值` 的格式（等号可选，空格分隔即可）
4. **嵌套结构**：节段可以嵌套，形成层级结构
5. **注释**：以 `!` 或 `#` 开头的行为注释

### 2.2 单位系统

CP2K 使用以下默认单位（重要！）：

| 物理量 | 默认单位 | 转换关系 |
|--------|----------|----------|
| 长度 | Angstrom (Å) | 1 Å = 0.529177 Bohr |
| 能量 | Hartree (Ha) | 1 Ha = 27.211386 eV = 627.509 kcal/mol |
| 力 | Hartree/Bohr | |
| 截断能 | Rydberg (Ry) | 1 Ry = 13.605693 eV |
| 温度 | Kelvin | |
| 时间 | 飞秒 (fs) | |

> **注意**：CP2K 中的 CUTOFF（截断能）使用 **Rydberg** 单位，而非大多数其他 DFT 软件使用的 Hartree 或 eV。这是初学者常见的混淆点。

如果需要指定其他单位，可以在关键字后添加单位标记：

```
ABC  5.43  5.43  5.43                     ! 默认单位 Angstrom
ABC  10.258  10.258  10.258  [bohr]       ! 使用 Bohr 单位
CUTOFF  300                               ! 默认单位 Rydberg
CUTOFF  40.817  [eV]                      ! 使用 eV 单位
```

### 2.3 输入文件的整体结构

一个典型的 CP2K 输入文件由以下主要节段组成：

```
&GLOBAL
   ...全局设置...
&END GLOBAL

&FORCE_EVAL
   ...力评估设置（方法、DFT参数等）...
   &SUBSYS
      ...体系定义（原子坐标、晶胞等）...
   &END SUBSYS
&END FORCE_EVAL
```

后续章节中我们将逐步介绍更多节段，如 `&MOTION`（用于几何优化和分子动力学）等。

---

## 3. 关键概念解析

### 3.1 `&GLOBAL` 节段

`&GLOBAL` 节段包含控制整个计算流程的全局参数。它**必须**出现在输入文件的最外层。

**常用关键字**：

| 关键字 | 说明 | 常用值 |
|--------|------|--------|
| `PROJECT` | 项目名称，决定输出文件的前缀 | 任意字符串，如 `He_atom` |
| `RUN_TYPE` | 计算类型 | `ENERGY`, `ENERGY_FORCE`, `GEO_OPT`, `CELL_OPT`, `MD`, `VIBRATIONAL_ANALYSIS` 等 |
| `PRINT_LEVEL` | 输出详细程度 | `LOW`, `MEDIUM`, `HIGH`, `DEBUG` |

**RUN_TYPE 常用选项详解**：

- **`ENERGY`**：仅计算单点能（本章使用）
- **`ENERGY_FORCE`**：计算能量和原子受力
- **`GEO_OPT`**：几何优化（寻找能量最低的原子构型）
- **`CELL_OPT`**：晶胞优化（同时优化晶胞参数和原子位置）
- **`MD`**：分子动力学模拟
- **`VIBRATIONAL_ANALYSIS`**：振动频率分析（简正模式）
- **`BAND`**：NEB（Nudged Elastic Band）方法计算反应路径

**示例**：

```fortran
&GLOBAL
  PROJECT  He_atom       ! 输出文件将命名为 He_atom-*.out, He_atom-*.ener 等
  RUN_TYPE ENERGY         ! 进行单点能计算
  PRINT_LEVEL LOW         ! 输出信息较少，适合初学时使用
&END GLOBAL
```

### 3.2 `&FORCE_EVAL` 节段

`&FORCE_EVAL`（Force Evaluation）节段是输入文件的核心，定义了如何计算原子间的相互作用力和能量。

**关键关键字**：

| 关键字 | 说明 | 常用值 |
|--------|------|--------|
| `METHOD` | 力评估方法 | `Quickstep`（DFT）, `Fist`（经典力场）, `QMMM`（QM/MM混合）|

- **`Quickstep`**：CP2K 的 DFT 引擎，使用 GPW（Gaussian and Plane Waves）方法。这是做量子力学计算时最常用的选项。
- **`Fist`**：经典分子力学模块，使用力场势函数，适用于大体系的经典模拟。
- **`QMMM`**：量子力学/分子力学耦合方法，体系的一部分用 QM 处理，其余用 MM 处理。

在 `&FORCE_EVAL` 节段内部，通常包含：

- `&DFT`：DFT 相关参数（基组文件、截断能、SCF 参数、泛函等）
- `&SUBSYS`：体系结构定义

### 3.3 `&SUBSYS` 节段

`&SUBSYS`（Subsystem）节段定义了模拟体系的结构信息，包括原子种类、原子坐标和模拟晶胞。

**常用子节段**：

#### `&CELL` — 晶胞定义

定义模拟晶胞的形状和大小。对于周期性 DFT 计算，即使计算一个孤立分子，也需要定义一个足够大的晶胞以避免周期性镜像之间的相互作用。

```fortran
&CELL
  ABC  10.0  10.0  10.0    ! 晶胞的三条边长（Å），默认正交晶胞
&END CELL
```

对于非正交晶胞，可以指定晶格矢量：

```fortran
&CELL
  A  5.0  0.0  0.0
  B  2.5  4.33  0.0
  C  0.0  0.0  7.0
&END CELL
```

#### `&COORD` — 原子坐标

定义体系中所有原子的元素种类和笛卡尔坐标（默认单位：Å）。

```fortran
&COORD
  He  0.0  0.0  0.0      ! 元素符号  X  Y  Z
&END COORD
```

也可以使用 `SCALED` 关键字输入分数坐标：

```fortran
&COORD
  Si  0.0  0.0  0.0
  Si  0.25  0.25  0.25
&END COORD
  SCALED .TRUE.           ! 使用分数坐标
```

### 3.4 `&KIND` 节段

`&KIND` 节段定义了每种元素的计算参数。它是 `&SUBSYS` 的子节段，**每种出现的元素都需要一个对应的 `&KIND` 定义**。

```fortran
&KIND He
  ELEMENT   He                     ! 元素符号
  BASIS_SET DZVP-GTH-PADE          ! 基组名称
  POTENTIAL GTH-PADE-q2            ! 赝势名称
&END KIND
```

**关键说明**：

- **`ELEMENT`**：指定元素符号
- **`BASIS_SET`**：指定基组，命名格式为 `类型-赝势类型-泛函`，如 `DZVP-GTH-PADE`
- **`POTENTIAL`**：指定赝势，命名格式为 `GTH-泛函-价电子数`，如 `GTH-PADE-q2` 表示使用 PADE 泛函对应的 GTH 赝势，价电子数为 2
- **赝势中 `qN` 的含义**：`qN` 表示赝势处理的价电子数。例如 He 的电子构型为 1s²，全部 2 个电子都是价电子，所以用 `q2`

> **重要提示**：基组、赝势和交换关联泛函三者必须匹配！例如使用 PADE (LDA) 泛函时，应选择 `GTH-PADE` 系列的赝势和基组。

---

## 4. 示例详解：helium.inp

下面我们逐行分析本章的示例文件 `helium.inp`——一个氦原子的单点能计算。

### 4.1 完整文件

```fortran
&GLOBAL
  PROJECT He_atom
  RUN_TYPE ENERGY
  PRINT_LEVEL LOW
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &SUBSYS
    &KIND He
      ELEMENT He
      BASIS_SET DZVP-GTH-PADE
      POTENTIAL GTH-PADE-q2
    &END KIND
    &CELL
      ABC 10.0 10.0 10.0
    &END CELL
    &COORD
      He 0.0 0.0 0.0
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
  &END DFT
&END FORCE_EVAL
```

### 4.2 逐段解析

#### 全局设置

```fortran
&GLOBAL
  PROJECT He_atom
  RUN_TYPE ENERGY
  PRINT_LEVEL LOW
&END GLOBAL
```

- `PROJECT He_atom`：项目名称设为 `He_atom`，输出文件将以 `He_atom` 为前缀
- `RUN_TYPE ENERGY`：只计算单点能（不做几何优化、动力学等）
- `PRINT_LEVEL LOW`：输出级别设为低，只打印必要信息

#### 力评估方法

```fortran
&FORCE_EVAL
  METHOD Quickstep
```

使用 `Quickstep`（即 GPW 方法的 DFT 引擎）来计算力和能量。

#### 体系定义

```fortran
  &SUBSYS
    &KIND He
      ELEMENT He
      BASIS_SET DZVP-GTH-PADE
      POTENTIAL GTH-PADE-q2
    &END KIND
```

定义氦原子的计算参数：
- 基组：`DZVP-GTH-PADE`（Double-Zeta Valence with Polarization，针对 PADE 泛函优化的 GTH 赝势配套基组）
- 赝势：`GTH-PADE-q2`（PADE 泛函对应的 GTH 赝势，价电子数 2）

```fortran
    &CELL
      ABC 10.0 10.0 10.0
    &END CELL
```

定义一个 10x10x10 Å 的立方模拟盒子。对于孤立原子，这个盒子足够大，可以忽略周期性镜像的相互作用（经验法则：盒子边长应至少比原子间距大 5-6 Å）。

```fortran
    &COORD
      He 0.0 0.0 0.0
    &END COORD
  &END SUBSYS
```

将氦原子放在盒子的原点位置。

#### DFT 参数

```fortran
  &DFT
    BASIS_SET_FILE_NAME BASIS_SET
    POTENTIAL_FILE_NAME GTH_POTENTIALS
```

指定基组和赝势的查找文件。`BASIS_SET` 和 `GTH_POTENTIALS` 是 CP2K 自带的数据文件，通常位于 CP2K 的 `data/` 目录下或通过环境变量 `CP2K_DATA_DIR` 指定。

```fortran
    &QS
      EPS_DEFAULT 1.0E-10
    &END QS
```

`&QS` 节段控制 Quickstep 方法的通用参数。`EPS_DEFAULT 1.0E-10` 设置各种数值积分的默认精度阈值（较高精度，但计算稍慢）。

```fortran
    &MGRID
      CUTOFF 300
      NGRIDS 4
      REL_CUTOFF 60
    &END MGRID
```

`&MGRID` 节段定义实空间网格（Multi-Grid）参数：

- **`CUTOFF 300`**：平面波截断能为 300 Rydberg。该值控制最精细网格的分辨率。越大越精确，但计算量也越大。后面第 4 章将详细讲解。
- **`NGRIDS 4`**：使用 4 层网格。CP2K 的多网格方法将不同宽度的高斯函数分配到不同粗细的网格上。
- **`REL_CUTOFF 60`**：参考截断能为 60 Ry。决定高斯函数如何分配到各层网格。后面第 4 章将详细讲解。

```fortran
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-7
      MAX_SCF 50
```

`&SCF` 节段控制自洽场 (SCF) 迭代过程：

- **`SCF_GUESS ATOMIC`**：初始电子密度猜测方式为原子密度叠加。其他选项还有 `RESTART`（从上次计算重启）、`RANDOM`（随机）等。
- **`EPS_SCF 1.0E-7`**：SCF 收敛判据——当总能量变化小于此值时，认为 SCF 收敛。
- **`MAX_SCF 50`**：最大 SCF 迭代次数。如果 50 步内未收敛，CP2K 将报错退出。

```fortran
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
```

选择标准对角化算法求解 Kohn-Sham 方程。标准对角化的计算标度为 O(N³)，适合小体系。大体系可以使用 OT（Orbital Transformation）方法。

```fortran
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.4
        NBROYDEN 8
      &END MIXING
    &END SCF
```

密度混合方法：

- **`METHOD BROYDEN_MIXING`**：使用 Broyden 混合方案加速 SCF 收敛。其他常用选项：`PULAY_MIXING`、`DIRECT_P_MIXING`
- **`ALPHA 0.4`**：混合参数（0 到 1 之间）。值越小越稳定但收敛越慢。
- **`NBROYDEN 8`**：Broyden 方法使用的前几步信息数量。

#### 交换关联泛函

```fortran
    &XC
      &XC_FUNCTIONAL PADE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
&END FORCE_EVAL
```

指定交换关联泛函为 **PADE**。PADE 是 LDA（Local Density Approximation，局域密度近似）泛函的一种参数化形式，是最简单但也是最基础的泛函。后续章节会介绍更多泛函选择。

---

## 5. 如何运行CP2K

### 5.1 单进程运行

最简单的运行方式：

```bash
cp2k.popt -o output.out input.inp
```

- `cp2k.popt`：CP2K 的可执行文件（纯 MPI 版本）
- `-o output.out`：指定输出文件名
- `input.inp`：输入文件

对于本章示例：

```bash
cp2k.popt -o He_atom.out helium.inp
```

### 5.2 多进程并行运行

CP2K 支持 MPI 并行，可以显著加速计算：

```bash
mpirun -n 4 cp2k.popt -o He_atom.out helium.inp
```

将 `-n` 后的数字替换为您想使用的 MPI 进程数。

### 5.3 混合并行（MPI + OpenMP）

对于支持 OpenMP 的编译版本（`cp2k.psmp`）：

```bash
export OMP_NUM_THREADS=2
mpirun -n 4 cp2k.psmp -o He_atom.out helium.inp
```

总核心数 = MPI 进程数 × OpenMP 线程数 = 4 × 2 = 8。

### 5.4 可执行文件说明

| 可执行文件 | 说明 |
|-----------|------|
| `cp2k.popt` | 纯 MPI 并行（优化版） |
| `cp2k.psmp` | MPI + OpenMP 混合并行 |
| `cp2k.ssmp` | 纯 OpenMP（单节点） |
| `cp2k.sopt` | 串行（无并行） |

> **提示**：不同系统上的可执行文件名称可能不同。使用 `which cp2k.popt` 或 `which cp2k.psmp` 查找安装位置。

### 5.5 运行后生成的文件

运行成功后，会生成以下文件：

| 文件 | 说明 |
|------|------|
| `He_atom.out` | 主输出文件，包含计算结果和详细日志 |
| `He_atom-RESTART.wfn` | 波函数重启文件（可用于续算） |
| `He_atom-1_0.restart` | 重启文件 |

---

## 6. 输出文件解读

打开 `He_atom.out` 文件，以下是需要关注的关键输出部分：

### 6.1 版本和编译信息

```
 DBCSR| CPU multiprocessor for this message...
 **** **** ******  **  PROGRAM STARTED AT               2024-...
 CP2K| version string:                                  CP2K 2024.1
```

确认 CP2K 版本信息。

### 6.2 输入参数回显

CP2K 会在输出文件开头回显完整的输入文件内容（经过解析后），这是检查输入是否正确的好机会。

### 6.3 SCF 迭代过程

```
 SCF WAVEFUNCTION OPTIMIZATION

  Step     Update method      Time    Convergence         Total energy    Change
  ------------------------------------------------------------------------------
    1 Broy./Diag. 0.40E+00    0.5     0.00000045        -2.9077946612     -2.91E+00
    2 Broy./Diag. 0.40E+00    0.3     0.00000003        -2.9077954801     -8.2E-07
    3 Broy./Diag. 0.40E+00    0.3     0.00000000        -2.9077954823     -2.2E-09

  *** SCF run converged in     3 steps ***
```

- **Convergence**：SCF 残差，越小越好
- **Total energy**：当前步骤的总能量（Hartree）
- **Change**：与上一步的能量差
- 当 Change 小于 `EPS_SCF`（1.0E-7）时，SCF 收敛

### 6.4 最终能量

```
 ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):           -2.90779548230000
```

这是最终的总能量，以 **Hartree** 为单位。对于氦原子，PADE/LDA 泛函下的参考值约为 -2.9078 Ha。

### 6.5 计算时间统计

输出文件末尾包含详细的计时信息：

```
 -------------------------------------------------------------------------------
 -                                T I M I N G                                -
 -------------------------------------------------------------------------------
 ...
 TOTAL                                                  ...
```

---

## 7. 练习

完成以下练习来巩固本章所学内容：

### 练习 1：计算氖原子 (Ne) 的能量

创建一个新文件 `neon.inp`，将氦改为氖：
- `&KIND` 中将 `He` 改为 `Ne`
- 赝势改为 `GTH-PADE-q8`（氖有 8 个价电子：2s²2p⁶）
- 基组保持 `DZVP-GTH-PADE`
- `PROJECT` 改为 `Ne_atom`

预期总能量约为 -128.5 Ha。

### 练习 2：更换基组

在同一目录下，将 `helium.inp` 复制一份，把基组从 `DZVP-GTH-PADE` 改为 `SZV-GTH-PADE`（单 zeta 基组）和 `TZV2P-GTH-PADE`（三 zeta 基组），分别计算能量。

观察不同基组对总能量的影响。预期结果：基组越大，能量越低（变分原理）。

### 练习 3：更改输出级别

将 `PRINT_LEVEL` 从 `LOW` 改为 `MEDIUM` 或 `HIGH`，重新运行并观察输出文件的变化。输出中增加了哪些信息？

### 练习 4：尝试不同的 RUN_TYPE

将 `RUN_TYPE` 改为 `ENERGY_FORCE`，观察输出文件中是否多出了力的信息。对于处在原点的单个原子，总力应该是多少？（提示：对称性）

### 思考题

1. 为什么计算一个孤立原子需要定义晶胞？如果晶胞太小会怎样？
2. PADE 是什么类型的泛函？它和 PBE 有什么区别？
3. `GTH-PADE-q2` 中的 `q2` 代表什么含义？

---

## 参考资源

- CP2K 官方文档：https://manual.cp2k.org/
- CP2K 基组和赝势数据：https://github.com/cp2k/cp2k/tree/master/data
- Quickstep (GPW) 方法论文：J. VandeVondele, M. Krack, F. Mohamed, M. Parrinello, T. Chassaing, J. Hutter, *Comput. Phys. Commun.* **167**, 103 (2005)
