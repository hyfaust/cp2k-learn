#!/bin/bash
# =============================================================================
# scan_cutoff.sh — CUTOFF 收敛性扫描脚本
# 固定 REL_CUTOFF = 60 Ry，扫描 CUTOFF 从 50 到 500 Ry
# =============================================================================

set -e

# CP2K 可执行文件路径（请根据实际安装路径修改）
CP2K_EXE="cp2k.popt"
NP=1  # MPI 进程数

# 模板输入文件
TEMPLATE="Si_cutoff.inp"

# CUTOFF 值列表 (Ry)
CUTOFF_VALUES=(50 100 150 200 250 300 350 400 450 500)

# 结果输出文件
RESULTS_FILE="cutoff_results.dat"

echo "============================================================"
echo "  CUTOFF 收敛性扫描"
echo "  固定 REL_CUTOFF = 60 Ry"
echo "============================================================"
echo ""

# 检查模板文件是否存在
if [ ! -f "$TEMPLATE" ]; then
    echo "[错误] 模板文件 $TEMPLATE 不存在！"
    exit 1
fi

# 写入结果文件头
echo "# CUTOFF(Ry)  Energy(Hartree)" > "$RESULTS_FILE"

for CUTOFF in "${CUTOFF_VALUES[@]}"; do
    INPUT="Si_cutoff_${CUTOFF}Ry.inp"
    OUTPUT="Si_cutoff_${CUTOFF}Ry.out"
    DIR="run_cutoff_${CUTOFF}"

    echo ">>> CUTOFF = ${CUTOFF} Ry"

    # 创建运行目录
    mkdir -p "$DIR"

    # 用 sed 替换模板中的 CUTOFF 占位符
    sed "s/CUTOFFVAL/${CUTOFF}/g" "$TEMPLATE" > "$DIR/$INPUT"

    # 运行 CP2K
    cd "$DIR"
    if [ "$NP" -gt 1 ]; then
        mpirun -n "$NP" "../$CP2K_EXE" -o "$OUTPUT" "$INPUT"
    else
        "../$CP2K_EXE" -o "$OUTPUT" "$INPUT"
    fi
    cd ..

    # 提取总能量
    ENERGY=$(grep "ENERGY| Total FORCE_EVAL" "$DIR/$OUTPUT" | tail -1 | awk '{print $NF}')

    if [ -z "$ENERGY" ]; then
        echo "    [警告] 未能提取能量，检查 $DIR/$OUTPUT"
        continue
    fi

    echo "    总能量: $ENERGY Hartree"

    # 记录结果
    echo "$CUTOFF $ENERGY" >> "$RESULTS_FILE"
done

echo ""
echo "============================================================"
echo "  扫描完成！结果保存在: $RESULTS_FILE"
echo "============================================================"
echo ""

# 打印结果表格
echo "结果汇总："
echo "------------------------------------------------------------"
printf "%-10s %20s %20s\n" "CUTOFF(Ry)" "能量(Hartree)" "ΔE(与上一个)"
echo "------------------------------------------------------------"

PREV_ENERGY=""
while IFS=' ' read -r cut en; do
    # 跳过注释行
    [[ "$cut" == \#* ]] && continue

    if [ -n "$PREV_ENERGY" ]; then
        DIFF=$(awk "BEGIN {printf \"%.10f\", $en - $PREV_ENERGY}")
        DIFF_EV=$(awk "BEGIN {printf \"%.6f\", ($en - $PREV_ENERGY) * 27.211386}")
        printf "%-10s %20s %12s Ha (%s eV)\n" "$cut" "$en" "$DIFF" "$DIFF_EV"
    else
        printf "%-10s %20s\n" "$cut" "$en"
    fi
    PREV_ENERGY="$en"
done < "$RESULTS_FILE"

echo ""
echo "提示：使用 plot_convergence.py 绘制收敛曲线。"
