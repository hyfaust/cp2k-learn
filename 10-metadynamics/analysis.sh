#!/bin/bash
# ============================================================
# analysis.sh — 元动力学HILLS文件分析脚本
# ============================================================
# 用法: bash analysis.sh [HILLS_FILE]
# 默认读取当前目录下的HILLS文件
# ============================================================

set -e

# ---- 配置 ----
HILLS_FILE="${1:-HILLS}"
OUTPUT_FES="fes.dat"
OUTPUT_CVS="cv_timeseries.dat"
OUTPUT_HILLS_STATS="hills_statistics.dat"

echo "============================================"
echo " 元动力学分析脚本"
echo "============================================"
echo "HILLS文件: ${HILLS_FILE}"
echo ""

# ---- 检查文件存在 ----
if [ ! -f "${HILLS_FILE}" ]; then
    echo "错误: 找不到HILLS文件 '${HILLS_FILE}'"
    echo "请确认元动力学模拟已经运行，并生成了HILLS文件。"
    echo ""
    echo "用法: bash analysis.sh [HILLS_FILE]"
    exit 1
fi

# ---- 基本统计 ----
echo "=== HILLS文件基本信息 ==="
N_LINES=$(wc -l < "${HILLS_FILE}")
N_COMMENTS=$(grep -c '^#' "${HILLS_FILE}" || true)
N_HILLS=$((N_LINES - N_COMMENTS))
echo "总行数: ${N_LINES}"
echo "注释行数: ${N_COMMENTS}"
echo "高斯数目: ${N_HILLS}"
echo ""

# ---- 提取CV时间序列 ----
echo "=== 提取CV时间序列 ==="
grep -v '^#' "${HILLS_FILE}" | awk '{
    print $1, $2, $3
}' > "${OUTPUT_CVS}"
echo "CV时间序列已保存到: ${OUTPUT_CVS}"

# ---- CV统计 ----
echo ""
echo "=== CV统计信息 ==="
echo "--- CV1 (第1个集合变量) ---"
grep -v '^#' "${HILLS_FILE}" | awk '{
    sum += $2; sumsq += $2*$2; n++;
    if(n==1 || $2 < min) min=$2;
    if(n==1 || $2 > max) max=$2;
} END {
    mean = sum/n;
    std = sqrt(sumsq/n - mean*mean);
    printf "  最小值: %.4f\n", min;
    printf "  最大值: %.4f\n", max;
    printf "  平均值: %.4f\n", mean;
    printf "  标准差: %.4f\n", std;
    printf "  波动范围: %.4f\n", max-min;
}'

echo "--- CV2 (第2个集合变量) ---"
grep -v '^#' "${HILLS_FILE}" | awk '{
    sum += $3; sumsq += $3*$3; n++;
    if(n==1 || $3 < min) min=$3;
    if(n==1 || $3 > max) max=$3;
} END {
    mean = sum/n;
    std = sqrt(sumsq/n - mean*mean);
    printf "  最小值: %.4f\n", min;
    printf "  最大值: %.4f\n", max;
    printf "  平均值: %.4f\n", mean;
    printf "  标准差: %.4f\n", std;
    printf "  波动范围: %.4f\n", max-min;
}'
echo ""

# ---- 高斯高度统计 ----
echo "=== 高斯高度统计 ==="
echo "(用于Well-Tempered Metadynamics时，高度会逐渐减小)"
grep -v '^#' "${HILLS_FILE}" | awk '{
    if(NR==1) first_h=$4;
    last_h=$4;
    sum += $4; n++;
    if(n==1 || $4 < min) min=$4;
    if(n==1 || $4 > max) max=$4;
} END {
    printf "  初始高度: %.6e\n", first_h;
    printf "  最终高度: %.6e\n", last_h;
    printf "  最小高度: %.6e\n", min;
    printf "  最大高度: %.6e\n", max;
    printf "  平均高度: %.6e\n", sum/n;
    if(first_h > 0) {
        printf "  最终/初始比: %.4f\n", last_h/first_h;
    }
}'
echo ""

# ---- 尝试使用CP2K graph工具重建FES ----
echo "=== FES重建 ==="
GRAPH_TOOL=""
# 搜索graph工具
for tool in graph graph.psmp graph.ssmp; do
    if command -v "${tool}" &> /dev/null; then
        GRAPH_TOOL="${tool}"
        break
    fi
done

if [ -n "${GRAPH_TOOL}" ]; then
    echo "找到graph工具: ${GRAPH_TOOL}"
    echo "正在重建自由能面..."

    # 检测维度
    NCOLS=$(grep -v '^#' "${HILLS_FILE}" | head -1 | awk '{print NF}')
    # 通常: step, cv1, cv2, sigma1, sigma2, height, biasf = 7列 (2D)
    #        step, cv1, sigma1, height, biasf = 5列 (1D)

    if [ "${NCOLS}" -ge 7 ]; then
        NDIM=2
        echo "检测到2D集合变量空间"
        ${GRAPH_TOOL} -i "${HILLS_FILE}" -stride 10 -ndim 2 -ndw 1 2 -cp2k -integrated_fes 2>/dev/null && \
            echo "FES已保存到: fes.dat" || \
            echo "警告: graph工具执行失败，请手动运行。"
    elif [ "${NCOLS}" -ge 5 ]; then
        NDIM=1
        echo "检测到1D集合变量空间"
        ${GRAPH_TOOL} -i "${HILLS_FILE}" -stride 10 -ndim 1 -ndw 1 -cp2k -integrated_fes 2>/dev/null && \
            echo "FES已保存到: fes.dat" || \
            echo "警告: graph工具执行失败，请手动运行。"
    else
        echo "警告: 无法自动检测CV维度 (${NCOLS}列)"
    fi
else
    echo "未找到graph工具。"
    echo "请手动运行以下命令重建FES："
    echo ""
    echo "  # 2D FES (例如 phi/psi):"
    echo "  graph.psmp -i HILLS -stride 10 -ndim 2 -ndw 1 2 -cp2k -integrated_fes"
    echo ""
    echo "  # 1D FES:"
    echo "  graph.psmp -i HILLS -stride 10 -ndim 1 -ndw 1 -cp2k -integrated_fes"
    echo ""
    echo "  # 注意: graph工具是CP2K安装的一部分，"
    echo "  # 通常位于CP2K的tools/目录或已加入PATH。"
fi
echo ""

# ---- 收敛性检查建议 ----
echo "=== 收敛性检查建议 ==="
echo "1. 比较不同时刻的FES："
HALF_HILLS=$((N_HILLS / 2))
echo "   - 前半段: graph.psmp -i HILLS -stride 10 -ndim 2 -ndw 1 2 -cp2k -integrated_fes -skip ${HALF_HILLS}"
echo "   - 全部:   graph.psmp -i HILLS -stride 10 -ndim 2 -ndw 1 2 -cp2k -integrated_fes"
echo ""
echo "2. 检查CV是否多次访问所有区域（查看${OUTPUT_CVS}）"
echo ""
echo "3. 对于Well-Tempered Metadynamics，检查高度是否趋近于零"
echo ""

# ---- 生成gnuplot脚本（可选） ----
echo "=== 生成可视化脚本 ==="
cat > plot_cv.gp << 'GNUPLOT_EOF'
# CV时间序列绘图脚本 (gnuplot)
# 用法: gnuplot plot_cv.gp

set terminal pngcairo size 1000,800 enhanced font 'Arial,14'
set output 'cv_timeseries.png'

set title 'Collective Variable Time Evolution' font 'Arial,16'
set xlabel 'Step'
set ylabel 'CV Value (rad)'
set grid

set key top right

plot 'cv_timeseries.dat' using 1:2 with lines title 'CV1 (phi)' lw 1, \
     'cv_timeseries.dat' using 1:3 with lines title 'CV2 (psi)' lw 1
GNUPLOT_EOF
echo "gnuplot脚本已生成: plot_cv.gp"

cat > plot_fes.gp << 'GNUPLOT_EOF'
# FES等高线绘图脚本 (gnuplot)
# 用法: gnuplot plot_fes.gp
# 前提: 已运行graph工具生成fes.dat

set terminal pngcairo size 900,800 enhanced font 'Arial,14'
set output 'fes_contour.png'

set title 'Free Energy Surface (kcal/mol)' font 'Arial,16'
set xlabel '{/Symbol f} (rad)'
set ylabel '{/Symbol y} (rad)'
set colorbox
set palette defined (0 "blue", 1 "cyan", 2 "green", 3 "yellow", 4 "red")

# 等高线图
set view map
set contour base
set cntrparam levels incremental 0, 0.5, 15

# 注意: fes.dat的格式可能需要根据实际输出调整
# splot 'fes.dat' using 1:2:3 with pm3d notitle
GNUPLOT_EOF
echo "gnuplot FES脚本已生成: plot_fes.gp"
echo ""

# ---- Python分析脚本 ----
cat > analyze_fes.py << 'PYTHON_EOF'
#!/usr/bin/env python3
"""
元动力学FES分析脚本
用法: python3 analyze_fes.py [fes.dat]
"""
import sys
import numpy as np

fes_file = sys.argv[1] if len(sys.argv) > 1 else 'fes.dat'

try:
    data = np.loadtxt(fes_file, comments='#')
except FileNotFoundError:
    print(f"错误: 找不到文件 {fes_file}")
    print("请先运行graph工具生成FES文件。")
    sys.exit(1)

print(f"FES文件: {fes_file}")
print(f"数据点数: {len(data)}")
print(f"列数: {data.shape[1]}")

if data.shape[1] >= 3:
    energy = data[:, -1]  # 最后一列是能量
    # 转换为kcal/mol
    energy_kcal = energy * 627.509

    print(f"\n自由能统计 (kcal/mol):")
    print(f"  最小值: {np.min(energy_kcal):.2f}")
    print(f"  最大值: {np.max(energy_kcal):.2f}")
    print(f"  势垒高度: {np.max(energy_kcal) - np.min(energy_kcal):.2f}")

    # 找到全局最小值
    min_idx = np.argmin(energy_kcal)
    print(f"\n全局最小值位置:")
    for i in range(data.shape[1] - 1):
        print(f"  CV{i+1}: {data[min_idx, i]:.4f}")
    print(f"  能量: {energy_kcal[min_idx]:.2f} kcal/mol")

print("\n提示: 使用matplotlib可视化FES:")
print("  import matplotlib.pyplot as plt")
print("  import numpy as np")
print(f"  data = np.loadtxt('{fes_file}', comments='#')")
print("  # 根据数据格式进行reshape和绘图")
PYTHON_EOF
echo "Python分析脚本已生成: analyze_fes.py"
echo ""

echo "============================================"
echo " 分析完成!"
echo "============================================"
echo ""
echo "生成的文件:"
echo "  ${OUTPUT_CVS}       — CV时间序列"
echo "  plot_cv.gp          — CV绘图gnuplot脚本"
echo "  plot_fes.gp         — FES绘图gnuplot脚本"
echo "  analyze_fes.py      — FES分析Python脚本"
echo ""
echo "后续步骤:"
echo "  1. gnuplot plot_cv.gp        — 查看CV演化"
echo "  2. graph.psmp -i HILLS ...   — 重建FES"
echo "  3. gnuplot plot_fes.gp       — 可视化FES"
echo "  4. python3 analyze_fes.py    — 定量分析"
