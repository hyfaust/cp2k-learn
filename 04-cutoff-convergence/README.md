# 第四章：收敛CUTOFF与REL_CUTOFF

> 难度: ⭐⭐ | 预计用时: 60-90 分钟

---

## 目录

1. [多网格 (Multi-Grid) 系统](#1-多网格-multi-grid-系统)
2. [CUTOFF 详解](#2-cutoff-详解)
3. [REL_CUTOFF 详解](#3-rel_cutoff-详解)
4. [收敛测试流程](#4-收敛测试流程)
5. [硅晶体示例](#5-硅晶体示例)
6. [能量收敛判据](#6-能量收敛判据)
7. [脚本详解](#7-脚本详解)
8. [结果分析](#8-结果分析)

---

## 1. 多网格 (Multi-Grid) 系统

### 1.1 为什么需要多网格

CP2K 的 QUICKSTEP 模块使用 GPW（Gaussian and Plane Waves）方法，其中电子密度和势函数在**实空间网格**上表示。不同的物理量具有不同的空间频率特征：

- **Hartree 势**：变化平缓，可以用较粗的网格
- **交换关联势**：变化较快，需要较细的网格
- **核心区域的密度**：变化剧烈，需要非常细的网格

如果只用一层最细的网格，虽然精度高，但计算和存储成本巨大。CP2K 的解决方案是使用**多网格 (Multi-Grid)** 技术：同时使用多层不同分辨率的网格，将不同宽度的高斯函数分配到合适的网格上。

### 1.2 多网格结构

CP2K 的多网格系统由参数 `NGRIDS` 控制（默认值为 4），形成一个网格层级：

```
网格层级 1 (最粗)  ←──  宽高斯函数（变化缓慢）
网格层级 2
网格层级 3
网格层级 4 (最细)  ←──  窄高斯函数（变化快速）
```

每个网格层级的间距 (h) 是前一层的 1/α 倍，其中 α 是**递进因子 (progression factor)**，默认值为 3.0：

$$h_i = h_{\text{finest}} \times \alpha^{\text{NGRIDS} - i}$$

其中 $i = 1, 2, \ldots, \text{NGRIDS}$，$\alpha \approx 2$ 是递进因子（refinement factor）。

### 1.3 高斯函数到网格的映射

每个高斯基函数根据其宽度（指数参数 α）被分配到**合适的网格层级**上：

- **窄高斯（大 α）**：变化剧烈 → 放在最细的网格
- **宽高斯（小 α）**：变化平缓 → 放在较粗的网格

具体的分配规则：一个高斯函数被放在满足以下条件的**最粗网格**上：

$$E_{\text{cut}}(\text{grid}) \geq \text{REL\_CUTOFF} \times \frac{4\alpha_{\text{Gaussian}}}{\pi}$$

即网格能够分辨该高斯函数的特征。这就是 `REL_CUTOFF` 的作用——它决定了高斯函数在不同网格间的分配方式。

### 1.4 三层网格示意图

假设使用 NGRIDS=4 和 CUTOFF=300 Ry：

```
网格 1 (最粗):  CUTOFF/α³ = 300/27 ≈ 11.1 Ry    ←  最宽的高斯函数
网格 2:         CUTOFF/α² = 300/9  ≈ 33.3 Ry    ←  
网格 3:         CUTOFF/α   = 300/3  = 100.0 Ry   ←  
网格 4 (最细):  CUTOFF = 300 Ry                   ←  最窄的高斯函数
```

---

## 2. CUTOFF 详解

### 2.1 什么是 CUTOFF

`CUTOFF` 是最精细网格的**平面波截断能**，以 **Rydberg (Ry)** 为单位。它决定了实空间网格的分辨率：

$$h = \frac{\pi}{\sqrt{2 \times \text{CUTOFF}}}$$

（在原子单位中）

或者用更直观的形式：

```
h (Å) ≈ 3.14159 / √(2 × CUTOFF × 0.367493)   (Ry 转 Ha 再转 Å)
```

对于 CUTOFF = 300 Ry：
```
h ≈ π / √(2 × 300 × 0.3675) ≈ π / √(220.5) ≈ π / 14.85 ≈ 0.2114 Bohr ≈ 0.112 Å
```

### 2.2 CUTOFF 的物理意义

在平面波基组中，截断能 $E_{\text{cut}}$ 限定了平面波展开的最大波矢：

$$G_{\max} = \sqrt{2 \times E_{\text{cut}}}$$

只有波矢 $|\mathbf{G}| \leq G_{\max}$ 的平面波分量才被包含。CUTOFF 越大：
- 包含越多的高频分量
- 网格越细（格点越密）
- 描述电子密度的精度越高
- 计算量越大（格点数正比于 $\text{CUTOFF}^{3/2}$）

### 2.3 CUTOFF 与各层网格截断能的关系

各层网格的截断能由以下公式决定：

$$E_{\text{cut}}^{(i)} = \frac{E_{\text{cut}}^{(1)}}{\alpha^{i-1}}$$

其中：
- $i = 1$ 为最粗网格，$i = \text{NGRIDS}$ 为最细网格
- $\alpha$ 是递进因子（默认 3.0）
- $E_{\text{cut}}^{(\text{NGRIDS})} = \text{CUTOFF}$（即输入中的 CUTOFF 值）

因此，当 α = 3.0，NGRIDS = 4，CUTOFF = 300 Ry 时：

| 网格层级 i | 截断能 (Ry) | 截断能 (Ha) | 网格间距 (Bohr) |
|-----------|------------|------------|----------------|
| 1 (最粗) | 300/27 ≈ 11.1 | 5.56 | 0.750 |
| 2 | 300/9 ≈ 33.3 | 16.67 | 0.433 |
| 3 | 300/3 = 100.0 | 50.00 | 0.250 |
| 4 (最细) | 300 | 150.00 | 0.144 |

### 2.4 CUTOFF 选择的经验法则

| 体系类型 | 推荐 CUTOFF (Ry) | 说明 |
|---------|------------------|------|
| 仅 s, p 元素 (H-F, Na-Cl) | 250-400 | 价电子波函数较光滑 |
| 包含 d 过渡金属 | 400-600 | d 轨道需要更细的网格 |
| 包含 f 镧系/锕系 | 600-1000 | f 轨道更局域 |
| 使用 MOLOPT 基组 | 300-600 | MOLOPT 基组可能需要更高的 CUTOFF |

> **永远不要只用一个 CUTOFF 值**！必须进行收敛测试，确认结果不随 CUTOFF 变化。

---

## 3. REL_CUTOFF 详解

### 3.1 什么是 REL_CUTOFF

`REL_CUTOFF` (Reference Cutoff) 是一个参考截断能（单位：Ry），控制高斯函数如何分配到不同粗细的网格层级上。

### 3.2 工作原理

对于每个高斯基函数，CP2K 根据其"宽度"（与指数参数 α 相关）确定其应该放在哪一层网格上。

判断规则：一个高斯函数被分配到满足以下条件的最粗网格：

$$E_{\text{cut}}^{\text{eff}} \geq \text{REL\_CUTOFF} \times \Delta G_{\text{Gaussian}}$$

其中 ΔG_Gaussian 与高斯函数的频率宽度相关。

直觉上理解：
- **REL_CUTOFF 较高** → 更多高斯函数被分配到细网格 → 更精确但更慢
- **REL_CUTOFF 较低** → 更多高斯函数被分配到粗网格 → 较快但可能不够精确

### 3.3 REL_CUTOFF 的影响

REL_CUTOFF 影响的是**精度**，而非**分辨率**。它决定了我们是否正确地将每个高斯函数放在了合适的网格上：

- 如果 REL_CUTOFF 太小：一些高斯函数被放在过粗的网格上，积分精度下降
- 如果 REL_CUTOFF 适中：每个高斯函数都在足够细的网格上计算
- 增大 REL_CUTOFF 超过某一阈值后：所有高斯函数都已在正确位置，继续增大无意义

### 3.4 默认值和典型选择

REL_CUTOFF 的**默认值为 50 Ry**。通常推荐使用 **60 Ry**，这是一个比较安全的选择。

在进行收敛测试时，通常扫描范围为：

```
REL_CUTOFF: 20, 30, 40, 50, 60, 70, 80, 90, 100, 120 Ry
```

大多数体系在 REL_CUTOFF = 60 Ry 时已经充分收敛。

---

## 4. 收敛测试流程

### 4.1 两步收敛策略

CP2K 中的网格收敛测试应分两步进行：

#### 第一步：收敛 CUTOFF

1. **固定 REL_CUTOFF** 为一个较高的值（如 60 Ry）
2. **扫描 CUTOFF** 从低到高（如 50 → 500 Ry，步长 50 Ry）
3. 绘制总能量 vs CUTOFF 曲线
4. 确定能量变化 < 阈值的 CUTOFF 值

#### 第二步：确认 REL_CUTOFF

1. **固定 CUTOFF** 为第一步确定的收敛值
2. **扫描 REL_CUTOFF**（如 20 → 120 Ry）
3. 确认能量不随 REL_CUTOFF 变化（或变化在阈值内）

### 4.2 为什么要先收敛 CUTOFF

CUTOFF 直接决定了网格分辨率，对总能量的影响更大。如果 CUTOFF 不够，增大 REL_CUTOFF 也没有意义——网格本身就太粗了。

REL_CUTOFF 的影响通常较小，且在 REL_CUTOFF ≥ 50-60 Ry 时很快饱和。

### 4.3 两步流程的代码逻辑

```
步骤 1:
  对于 CUTOFF = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]:
    REL_CUTOFF = 60  (固定)
    运行 CP2K 计算
    记录总能量
  绘图，确定收敛的 CUTOFF

步骤 2:
    CUTOFF = 已收敛值  (固定)
  对于 REL_CUTOFF = [20, 30, 40, 50, 60, 70, 80, 100, 120]:
    运行 CP2K 计算
    记录总能量
  确认收敛
```

### 4.4 收敛判据选择

常用的能量收敛判据：

| 判据 | 含义 | 适用场景 |
|------|------|----------|
| $\Delta E < 1$ meV/atom | 每原子能量变化 < 1 meV | 日常计算 |
| $\Delta E < 0.1$ meV/atom | 每原子能量变化 < 0.1 meV | 高精度计算 |
| $\Delta E < 10^{-6}$ Ha | 总能量变化 < 10⁻⁶ Hartree | 严格基准 |
| $\Delta F < 1$ meV/Å | 力的变化很小 | 几何优化/MD |

对于初步测试，$\Delta E < 1$ meV 通常已经足够。对于严格的基准计算或振动频率分析，需要更严格的判据。

---

## 5. 硅晶体示例

### 5.1 体系描述

本章使用**硅 (Si) 的 FCC 晶体结构**作为示例体系。选择硅的原因：

1. 硅是半导体材料的代表，具有重要应用价值
2. FCC 结构简单，容易理解
3. 硅的 d 电子效应适中，不会出现极端情况
4. 8 原子超胞大小适中，计算快速

### 5.2 硅的晶体结构

硅为金刚石结构（两个相互穿透的 FCC 子晶格），晶格常数 a = 5.4306975 Å。

8 原子超胞中 Si 原子的坐标（分数坐标和笛卡尔坐标）：

| 原子 | 分数坐标 | 笛卡尔坐标 (Å) |
|------|---------|----------------|
| Si₁ | (0, 0, 0) | (0.000, 0.000, 0.000) |
| Si₂ | (0, 0.5, 0.5) | (0.000, 2.715, 2.715) |
| Si₃ | (0.5, 0, 0.5) | (2.715, 0.000, 2.715) |
| Si₄ | (0.5, 0.5, 0) | (2.715, 2.715, 0.000) |
| Si₅ | (0.25, 0.25, 0.25) | (1.358, 1.358, 1.358) |
| Si₆ | (0.25, 0.75, 0.75) | (1.358, 4.073, 4.073) |
| Si₇ | (0.75, 0.25, 0.75) | (4.073, 1.358, 4.073) |
| Si₈ | (0.75, 0.75, 0.25) | (4.073, 4.073, 1.358) |

### 5.3 输入文件模板 (Si_cutoff.inp)

```fortran
&GLOBAL
  PROJECT Si_cutoff_CUTOFFVAL
  RUN_TYPE ENERGY
  PRINT_LEVEL LOW
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &SUBSYS
    &KIND Si
      ELEMENT Si
      BASIS_SET DZVP-GTH-PADE
      POTENTIAL GTH-PADE-q4
    &END KIND
    &CELL
      ABC 5.4306975 5.4306975 5.4306975
    &END CELL
    &COORD
      Si   0.000000   0.000000   0.000000
      Si   0.000000   2.715349   2.715349
      Si   2.715349   0.000000   2.715349
      Si   2.715349   2.715349   0.000000
      Si   1.357674   1.357674   1.357674
      Si   1.357674   4.073023   4.073023
      Si   4.073023   1.357674   4.073023
      Si   4.073023   4.073023   1.357674
    &END COORD
  &END SUBSYS
  &DFT
    ...
    &MGRID
      CUTOFF CUTOFFVAL        ! ← 占位符，由扫描脚本替换
      NGRIDS 4
      REL_CUTOFF 60
    &END MGRID
    ...
  &END DFT
&END FORCE_EVAL
```

**注意** `CUTOFF CUTOFFVAL` 中的 `CUTOFFVAL` 是一个占位符。扫描脚本 `scan_cutoff.sh` 会用 `sed` 将其替换为实际的 CUTOFF 值。

### 5.4 计算参数说明

| 参数 | 值 | 说明 |
|------|-----|------|
| 基组 | DZVP-GTH-PADE | 双 zeta + 极化 |
| 赝势 | GTH-PADE-q4 | Si 的 4 个价电子 |
| 泛函 | PADE | LDA |
| 晶胞 | 5.431×5.431×5.431 Å | 实验晶格常数 |
| 周期性 | XYZ（默认） | 三维周期性 |
| 原子数 | 8 | FCC 超胞 |

---

## 6. 能量收敛判据

### 6.1 绝对能量差

最直接的判据是相邻 CUTOFF 值之间的总能量差：

$$\Delta E = E(\text{CUTOFF}_n) - E(\text{CUTOFF}_{n-1})$$

当 $|\Delta E| <$ 阈值时，认为 CUTOFF 已收敛。

### 6.2 每原子能量差

对于不同大小的体系，使用每原子能量差更为公平：

$$\Delta E_{\text{per atom}} = \frac{\Delta E}{N_{\text{atoms}}}$$

### 6.3 常用阈值

| 阈值 | Hartree | eV | meV |
|------|---------|----|-----|
| 宽松 | 10⁻⁴ | 2.72 | 2721 |
| 日常 | 10⁻⁵ | 0.272 | 272 |
| 严格 | 10⁻⁶ | 0.0272 | 27.2 |
| 极严格 | 10⁻⁷ | 0.00272 | 2.72 |

**推荐**：对于 CUTOFF 收敛测试，使用 $\Delta E < 10^{-5}$ Ha（约 0.27 meV）作为收敛判据。这意味着进一步增大 CUTOFF 不会改变总能量超过 0.27 meV。

### 6.4 力的收敛

除了能量，还应检查力的收敛：

$$\Delta F_{\max} = \max_i |F_i(\text{CUTOFF}_n) - F_i(\text{CUTOFF}_{n-1})|$$

对于几何优化和分子动力学，力的收敛可能比能量收敛需要更高的 CUTOFF。

---

## 7. 脚本详解

### 7.1 scan_cutoff.sh

这个脚本自动扫描一系列 CUTOFF 值：

#### 核心逻辑

```bash
# 定义要扫描的 CUTOFF 值
CUTOFF_VALUES=(50 100 150 200 250 300 350 400 450 500)

for CUTOFF in "${CUTOFF_VALUES[@]}"; do
    # 创建独立的运行目录
    DIR="run_cutoff_${CUTOFF}"
    mkdir -p "$DIR"

    # 从模板生成输入文件（替换占位符）
    sed "s/CUTOFFVAL/${CUTOFF}/g" "$TEMPLATE" > "$DIR/$INPUT"

    # 运行 CP2K
    cd "$DIR"
    cp2k.popt -o "$OUTPUT" "$INPUT"
    cd ..

    # 提取能量
    ENERGY=$(grep "ENERGY| Total FORCE_EVAL" "$DIR/$OUTPUT" | tail -1 | awk '{print $NF}')
    
    # 记录结果
    echo "$CUTOFF $ENERGY" >> "$RESULTS_FILE"
done
```

#### 关键设计

1. **每个 CUTOFF 值独立目录**：避免输出文件互相覆盖，便于回溯检查
2. **sed 替换**：从模板文件动态生成输入，保证除 CUTOFF 外一切一致
3. **结果记录到文件**：`cutoff_results.dat` 格式为 `CUTOFF 值  能量`，方便后续绘图
4. **实时输出**：脚本运行时实时显示当前进度和能量值

#### 使用方法

```bash
# 确保在 04-cutoff-convergence/ 目录下
chmod +x scan_cutoff.sh

# 如需要，修改 CP2K_EXE 和 NP 变量
vim scan_cutoff.sh

# 运行
./scan_cutoff.sh
```

运行完成后，会生成：
- `run_cutoff_*/` 目录：每个 CUTOFF 值的计算结果
- `cutoff_results.dat`：汇总的 CUTOFF-能量数据

### 7.2 scan_rel_cutoff.sh

结构与 `scan_cutoff.sh` 类似，但扫描 REL_CUTOFF：

```bash
REL_CUTOFF_VALUES=(20 30 40 50 60 70 80 90 100 120)

for REL_CUTOFF in "${REL_CUTOFF_VALUES[@]}"; do
    # 从模板生成输入：CUTOFF 固定为 300，替换 REL_CUTOFF
    sed -e "s/CUTOFFVAL/300/g" \
        -e "s/REL_CUTOFF 60/REL_CUTOFF ${REL_CUTOFF}/g" \
        "$TEMPLATE" > "$DIR/$INPUT"
    ...
done
```

注意这里使用了两条 sed 替换规则：
1. 将 CUTOFF 占位符替换为固定值 300
2. 将默认的 REL_CUTOFF 60 替换为扫描值

### 7.3 plot_convergence.py

Python 脚本，用于可视化收敛数据：

#### 功能

1. 读取 `.dat` 数据文件
2. 绘制总能量 vs 参数曲线
3. 绘制能量差 $|\Delta E|$ vs 参数曲线（对数坐标）
4. 标注收敛阈值线（1 meV 和 $10^{-6}$ Ha）
5. 纯文本输出（无 matplotlib 时的回退方案）

#### 使用方法

```bash
# 绘制 CUTOFF 收敛曲线
python plot_convergence.py --data cutoff_results.dat

# 绘制 REL_CUTOFF 收敛曲线
python plot_convergence.py --data rel_cutoff_results.dat --rel-cutoff --title "REL_CUTOFF 收敛"

# 自定义输出文件名
python plot_convergence.py --data cutoff_results.dat --output my_convergence
```

#### 输出

- 终端文本表格（始终输出，包含能量和 $\Delta E$ 列表）
- `cutoff_convergence.png`（或自定义名称）：包含两个子图的收敛曲线
  - 上图：总能量 vs 参数
  - 下图：$|\Delta E|$ vs 参数（对数坐标，含阈值线）

#### 依赖

```bash
# 必需（Python 标准库即可运行文本模式）
python3 --version

# 可选（用于生成图片）
pip install matplotlib numpy
```

---

## 8. 结果分析

### 8.1 预期的 CUTOFF 收敛行为

对于 Si + DZVP-GTH-PADE + PADE 体系，预期的收敛行为如下：

```
CUTOFF (Ry)  | 能量 (Ha)          | ΔE (meV)
-------------|--------------------|----------
50           | 较高的能量（不精确）|  ---
100          | 显著下降           |  数百 meV
150          | 继续下降           |  数十 meV
200          | 趋于稳定           |  ~10 meV
250          | 基本收敛           |  < 5 meV
300          | 已收敛             |  < 1 meV
350-500      | 无显著变化         |  < 0.1 meV
```

### 8.2 收敛曲线的典型形状

CUTOFF 收敛曲线通常呈现以下特征：

```
能量
  ↑
  |  *
  |   *
  |     *
  |       *
  |         *  *  *  *  *    ← 已收敛区域（平台）
  |
  +---+---+---+---+---+---→ CUTOFF (Ry)
     50 100 150 200 250 300
```

- **低 CUTOFF 区域**：能量随 CUTOFF 增大快速下降（不精确区域）
- **过渡区域**：下降速率减缓
- **收敛区域**：能量基本不再变化（平台）
- **目标**：找到平台开始的 CUTOFF 值

### 8.3 常见问题及诊断

#### 问题一：能量振荡不收敛

```
能量
  ↑
  |  *     *
  |   *   * *   *
  |    * *   * *
  |     *      *
  +---+---+---+---→ CUTOFF
```

**可能原因**：CUTOFF 值的间距太大，没有捕捉到收敛的细节。

**解决**：缩小扫描步长（如 25 Ry 间距），特别是收敛过渡区域。

#### 问题二：能量持续下降不收敛

```
能量
  ↑
  |  *
  |   *
  |    *
  |     *
  |      *
  |       *
  +---+---+---+---→ CUTOFF
```

**可能原因**：
1. CUTOFF 范围不够大 → 增大最大 CUTOFF
2. 基组和赝势不匹配 → 检查 `&KIND` 中的设置
3. REL_CUTOFF 太小 → 试试 REL_CUTOFF = 60 或更高

#### 问题三：REL_CUTOFF 依赖性强

如果能量强烈依赖 REL_CUTOFF（即使在大值时），这可能意味着：
1. CUTOFF 不够大
2. 使用了非常弥散的基函数（如 aug-cc-pVDZ 等），需要更高的网格精度

### 8.4 最终推荐参数

基于收敛测试结果，选择满足以下条件的最小参数组合：

```
CUTOFF:     使得 ΔE < 1 meV 的最小值（通常为 300-400 Ry 对于 s/p 元素）
REL_CUTOFF: 60 Ry（通常足够，除非有特殊需求）
NGRIDS:     4（默认值，通常不需要修改）
```

对于后续所有计算，使用这些收敛的参数值。

### 8.5 不同体系的收敛行为对比

不同元素组成的体系，收敛速度不同：

| 体系 | 推荐 CUTOFF | 原因 |
|------|------------|------|
| 有机分子 (H, C, N, O) | 280-350 Ry | 仅 s, p 价电子 |
| 过渡金属配合物 | 400-600 Ry | d 轨道局域性强 |
| 含重元素 | 500-800 Ry | f 轨道更局域 |
| 使用 MOLOPT 基组 | 300-600 Ry | MOLOPT 基组的高斯可能较窄 |

### 8.6 将收敛参数应用到实际计算

确定了收敛的 CUTOFF 和 REL_CUTOFF 后，在后续所有计算中使用这些值：

```fortran
&MGRID
  CUTOFF 350        ! 收敛测试确定的值
  NGRIDS 4
  REL_CUTOFF 60     ! 通常 60 Ry 足够
&END MGRID
```

---

## 总结

| 概念 | 关键要点 |
|------|---------|
| 多网格 | 不同分辨率的网格，用于高效处理不同宽度的高斯函数 |
| CUTOFF | 最细网格的截断能，决定网格分辨率 |
| REL_CUTOFF | 参考截断能，控制高斯函数到网格的分配 |
| NGRIDS | 网格层数，默认 4 |
| $\alpha$ 因子 | 相邻网格的截断能比，默认 3.0 |
| 收敛策略 | 先固定 REL_CUTOFF 扫 CUTOFF，再确认 REL_CUTOFF |
| 收敛判据 | $\Delta E < 10^{-5}$ Ha（约 0.27 meV）或 $\Delta E < 1$ meV/atom |
| 典型值 | s/p 元素用 300 Ry，d 金属用 500 Ry |

---

## 参考资源

- Lippert, G.; Hutter, J.; Parrinello, M. *Mol. Phys.* **92**, 477 (1997) — GPW 方法
- Lippert, G.; Hutter, J.; Parrinello, M. *Theor. Chem. Acc.* **103**, 124 (1999) — GAPW 方法
- CP2K 手册 - MGRID 部分：https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/MGRID.html
- CP2K 手册 - QS 部分：https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/DFT/QS.html
