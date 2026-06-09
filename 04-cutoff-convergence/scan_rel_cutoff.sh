#!/bin/bash
# =============================================================================
# scan_rel_cutoff.sh — REL_CUTOFF 收敛性扫描脚本
# 固定 CUTOFF = 300 Ry，扫描 REL_CUTOFF 从 20 到 120 Ry
# =============================================================================

set -e

# CP2K 可执行文件路径（请根据实际安装路径修改）
CP2K_EXE="cp2k.popt"
NP=1  # MPI 进程数

# 模板输入文件 — 使用 Si_cutoff.inp 作为基础
TEMPLATE="Si_cutoff.inp"

# REL_CUTOFF 值列表 (Ry)
REL_CUTOFF_VALUES=(20 30 40 50 60 70 80 90 100 120)

# 结果输出文件
RESULTS_FILE="rel_cutoff_results.dat"

echo "============================================================"
echo "  REL_CUTOFF 收敛性扫描"
echo "  固定 CUTOFF = 300 Ry"
echo "============================================================"
echo ""

# 检查模板文件是否存在
if [ ! -f "$TEMPLATE" ]; then
    echo "[错误] 模板文件 $TEMPLATE 不存在！"
    exit 1
fi

# 写入结果文件头
echo "# REL_CUTOFF(Ry)  Energy(Hartree)" > "$RESULTS_FILE"

for REL_CUTOFF in "${REL_CUTOFF_VALUES[@]}"; do
    # 创建一个临时输入文件，同时替换 CUTOFF 和 REL_CUTOFF
    # 首先将模板的 CUTOFF 固定为 300，REL_CUTOFF 替换为当前值
    INPUT="Si_rel_cutoff_${REL_CUTOFF}Ry.inp"
    OUTPUT="Si_rel_cutoff_${REL_CUTOFF}Ry.out"
    DIR="run_rel_cutoff_${REL_CUTOFF}"

    echo ">>> REL_CUTOFF = ${REL_CUTOFF} Ry"

    # 创建运行目录
    mkdir -p "$DIR"

    # 用 sed 替换模板中的占位符
    # CUTOFFVAL -> 300 (固定)，REL_CUTOFFVAL -> 当前值
    sed -e "s/CUTOFFVAL/300/g" \
        -e "s/REL_CUTOFF 60/REL_CUTOFF ${REL_CUTOFF}/g" \
        "$TEMPLATE" > "$DIR/$INPUT"

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
    echo "$REL_CUTOFF $ENERGY" >> "$RESULTS_FILE"
done

echo ""
echo "============================================================"
echo "  扫描完成！结果保存在: $RESULTS_FILE"
echo "============================================================"
echo ""

# 打印结果表格
echo "结果汇总："
echo "------------------------------------------------------------"
printf "%-14s %20s %20s\n" "REL_CUTOFF(Ry)" "能量(Hartree)" "ΔE(与上一个)"
echo "------------------------------------------------------------"

PREV_ENERGY=""
while IFS=' ' read -r relcut en; do
    # 跳过注释行
    [[ "$relcut" == \#* ]] && continue

    if [ -n "$PREV_ENERGY" ]; then
        DIFF=$(awk "BEGIN {printf \"%.10f\", $en - $PREV_ENERGY}")
        DIFF_EV=$(awk "BEGIN {printf \"%.6f\", ($en - $PREV_ENERGY) * 27.211386}")
        printf "%-14s %20s %12s Ha (%s eV)\n" "$relcut" "$en" "$DIFF" "$DIFF_EV"
    else
        printf "%-14s %20s\n" "$relcut" "$en"
    fi
    PREV_ENERGY="$en"
done < "$RESULTS_FILE"

echo ""
echo "提示：使用 plot_convergence.py 绘制收敛曲线（传入 --rel-cutoff 选项）。"
