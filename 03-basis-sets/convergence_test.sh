#!/bin/bash
# =============================================================================
# convergence_test.sh — 基组收敛性测试脚本
# 使用不同基组（SZV, DZVP, TZV2P）计算水分子能量并比较
# =============================================================================

set -e  # 遇到错误立即退出

# CP2K 可执行文件路径（请根据实际安装路径修改）
CP2K_EXE="cp2k.popt"
NP=1  # MPI 进程数，根据您的机器调整

echo "============================================================"
echo "  基组收敛性测试 — 水分子 (H2O)"
echo "============================================================"
echo ""

# 定义基组列表
BASIS_SETS=("SZV" "DZVP" "TZV2P")
INPUT_FILES=("H2O_szv.inp" "H2O_dzvp.inp" "H2O_tzv2p.inp")

# 存储结果
declare -A ENERGIES

for i in "${!BASIS_SETS[@]}"; do
    BASIS="${BASIS_SETS[$i]}"
    INPUT="${INPUT_FILES[$i]}"
    OUTPUT="H2O_${BASIS,,}.out"  # 转小写

    echo ">>> 运行基组: $BASIS"
    echo "    输入文件: $INPUT"
    echo "    输出文件: $OUTPUT"

    if [ ! -f "$INPUT" ]; then
        echo "    [错误] 输入文件 $INPUT 不存在，跳过。"
        continue
    fi

    # 运行 CP2K
    if [ "$NP" -gt 1 ]; then
        mpirun -n "$NP" "$CP2K_EXE" -o "$OUTPUT" "$INPUT"
    else
        "$CP2K_EXE" -o "$OUTPUT" "$INPUT"
    fi

    # 提取总能量
    ENERGY=$(grep "ENERGY| Total FORCE_EVAL" "$OUTPUT" | tail -1 | awk '{print $NF}')
    ENERGIES["$BASIS"]="$ENERGY"
    echo "    总能量: $ENERGY Hartree"
    echo ""
done

# 打印汇总表
echo "============================================================"
echo "  结果汇总"
echo "============================================================"
printf "%-12s %20s\n" "基组" "总能量 (Hartree)"
echo "------------------------------------------------------------"

# 用于计算能量差
PREV_ENERGY=""
for BASIS in "${BASIS_SETS[@]}"; do
    ENERGY="${ENERGIES[$BASIS]:-N/A}"
    printf "%-12s %20s\n" "$BASIS" "$ENERGY"

    # 计算与上一个基组的能量差（使用 awk）
    if [ -n "$PREV_ENERGY" ] && [ "$ENERGY" != "N/A" ]; then
        DIFF=$(awk "BEGIN {printf \"%.10f\", $ENERGY - $PREV_ENERGY}")
        DELTA_EV=$(awk "BEGIN {printf \"%.6f\", ($ENERGY - $PREV_ENERGY) * 27.211386}")
        printf "%-12s %20s %15s eV\n" "" "ΔE = $DIFF Ha" "($DELTA_EV eV)"
    fi
    PREV_ENERGY="$ENERGY"
done

echo "============================================================"
echo ""
echo "预期结果："
echo "  - SZV 总能量最高（精度最低）"
echo "  - DZVP 能量显著低于 SZV"
echo "  - TZV2P 能量低于 DZVP，但与 DZVP 的差距远小于 SZV→DZVP"
echo "  - 如果 DZVP 和 TZV2P 能量差 < 0.001 Hartree (~0.027 eV)，"
echo "    则 DZVP 已经足够收敛。"
echo ""
echo "注意：能量本身不重要，重要的是能量差（如反应能、结合能）。"
echo "基组收敛性应当在能量差的层面上进行评估。"
