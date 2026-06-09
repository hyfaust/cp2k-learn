#!/bin/bash
# ============================================================
# compare_speed.sh — xTB vs DFT 速度比较脚本
# ============================================================
# 用法: bash compare_speed.sh
# 前提: 已准备好xTB和DFT输入文件，以及坐标文件
# ============================================================

set -e

# ---- 配置 ----
CP2K_BIN="cp2k.psmp"               # CP2K可执行文件
CP2K_XTB_BIN="cp2k.psmp"           # xTB也使用相同的CP2K
WATER_XYZ="water64.xyz"            # 水分子坐标文件
ETHANOL_XYZ="ethanol.xyz"          # 乙醇坐标文件

# 测试的系统大小列表
SYSTEM_SIZES=("1" "8" "27" "64")   # 水分子数目

echo "============================================"
echo " xTB vs DFT 速度比较"
echo "============================================"
echo ""

# ---- 检查CP2K可执行文件 ----
if ! command -v "${CP2K_BIN}" &> /dev/null; then
    echo "警告: 未找到CP2K可执行文件 '${CP2K_BIN}'"
    echo "请确保CP2K已安装并加入PATH。"
    echo "或者修改脚本中的CP2K_BIN变量。"
    echo ""
    echo "继续生成测试输入文件..."
fi

# ---- 生成DFT输入文件模板 ----
generate_dft_input() {
    local n_mol=$1
    local n_atoms=$((n_mol * 3))
    local box_size=$(echo "scale=2; ${n_mol}^(1.0/3.0) * 3.1" | bc -l 2>/dev/null || echo "12.4")
    local coord_file="water${n_mol}.xyz"
    local input_file="water${n_mol}_dft.inp"

    cat > "${input_file}" << EOF
@SET SYSTEM water${n_mol}_dft

&GLOBAL
  PROJECT \${SYSTEM}
  RUN_TYPE MD
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep

  &DFT
    BASIS_SET_FILE_NAME  BASIS_MOLOPT
    POTENTIAL_FILE_NAME  GTH_POTENTIALS

    &MGRID
      CUTOFF 400
      NGRIDS 4
    &END MGRID

    &QS
      EPS_DEFAULT 1.0E-12
      EXTRAPOLATION ASPC
      EXTRAPOLATION_ORDER 3
    &END QS

    &SCF
      MAX_SCF 50
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-6
      &OT
        PRECONDITIONER FULL_ALL
        MINIMIZER DIIS
      &END OT
      &OUTER_SCF
        MAX_SCF 10
        EPS_SCF 1.0E-6
      &END OUTER_SCF
    &END SCF

    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
      &VDW_POTENTIAL
        POTENTIAL_TYPE PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE DFTD3
          PARAMETER_FILE_NAME dftd3.dat
          REFERENCE_FUNCTIONAL PBE
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
    &END XC

    &POISSON
      PERIODIC XYZ
    &END POISSON

  &END DFT

  &SUBSYS
    &CELL
      ABC ${box_size} ${box_size} ${box_size}
      PERIODIC XYZ
    &END CELL

    &TOPOLOGY
      COORD_FILE_NAME ${coord_file}
      COORD_FILE_FORMAT XYZ
      CONNECTIVITY OFF
    &END TOPOLOGY

    &KIND H
      BASIS_SET DZVP-MOLOPT-GTH
      POTENTIAL GTH-PBE-q1
    &END KIND
    &KIND O
      BASIS_SET DZVP-MOLOPT-GTH
      POTENTIAL GTH-PBE-q6
    &END KIND

  &END SUBSYS

&END FORCE_EVAL

&MOTION
  &MD
    ENSEMBLE NVT
    TEMPERATURE 300.0
    TIMESTEP 1.0
    STEPS 5
    &THERMOSTAT
      REGION GLOBAL
      &CSVR
        TIMECON 100.0
      &END CSVR
    &END THERMOSTAT
  &END MD
&END MOTION
EOF
    echo "${input_file}"
}

# ---- 生成xTB输入文件模板 ----
generate_xtb_input() {
    local n_mol=$1
    local n_atoms=$((n_mol * 3))
    local box_size=$(echo "scale=2; ${n_mol}^(1.0/3.0) * 3.1" | bc -l 2>/dev/null || echo "12.4")
    local coord_file="water${n_mol}.xyz"
    local input_file="water${n_mol}_xtb.inp"

    cat > "${input_file}" << EOF
@SET SYSTEM water${n_mol}_xtb

&GLOBAL
  PROJECT \${SYSTEM}
  RUN_TYPE MD
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep

  &DFT
    &QS
      METHOD xTB
      &xTB
        CHECK_ATOMIC_CHARGES .FALSE.
        DO_EWALD .TRUE.
      &END xTB
    &END QS

    &SCF
      MAX_SCF 30
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-6
      &OT
        PRECONDITIONER FULL_ALL
        MINIMIZER DIIS
      &END OT
      &OUTER_SCF
        MAX_SCF 10
        EPS_SCF 1.0E-6
      &END OUTER_SCF
    &END SCF

    &POISSON
      PERIODIC XYZ
    &END POISSON

  &END DFT

  &SUBSYS
    &CELL
      ABC ${box_size} ${box_size} ${box_size}
      PERIODIC XYZ
    &END CELL

    &TOPOLOGY
      COORD_FILE_NAME ${coord_file}
      COORD_FILE_FORMAT XYZ
      CONNECTIVITY OFF
    &END TOPOLOGY

  &END SUBSYS

&END FORCE_EVAL

&MOTION
  &MD
    ENSEMBLE NVT
    TEMPERATURE 300.0
    TIMESTEP 1.0
    STEPS 5
    &THERMOSTAT
      REGION GLOBAL
      &CSVR
        TIMECON 100.0
      &END CSVR
    &END THERMOSTAT
  &END MD
&END MOTION
EOF
    echo "${input_file}"
}

# ---- 生成水分子XYZ文件 ----
generate_water_xyz() {
    local n_mol=$1
    local output="water${n_mol}.xyz"
    local n_atoms=$((n_mol * 3))

    echo "${n_atoms}" > "${output}"
    echo "Water box with ${n_mol} molecules" >> "${output}"

    # 生成简单的水分子坐标（简单立方排列）
    local idx=0
    local spacing=3.1   # Angstrom between molecules
    local n_side=$(echo "scale=0; l(${n_mol})/l(8)*3" | bc -l 2>/dev/null || echo "4")

    for i in $(seq 0 $((n_mol - 1))); do
        # 简单立方网格排列
        local ix=$((i % 4))
        local iy=$(((i / 4) % 4))
        local iz=$((i / 16))

        local x=$(echo "scale=3; ${ix} * ${spacing}" | bc -l)
        local y=$(echo "scale=3; ${iy} * ${spacing}" | bc -l)
        local z=$(echo "scale=3; ${iz} * ${spacing}" | bc -l)

        # 水分子几何 (O在中心, H在两侧)
        local h1x=$(echo "scale=3; ${x} + 0.757" | bc -l)
        local h1y=$(echo "scale=3; ${y} + 0.586" | bc -l)
        local h2x=$(echo "scale=3; ${x} - 0.757" | bc -l)
        local h2y=$(echo "scale=3; ${y} + 0.586" | bc -l)

        printf "O    %12.6f  %12.6f  %12.6f\n" "${x}" "${y}" "${z}" >> "${output}"
        printf "H    %12.6f  %12.6f  %12.6f\n" "${h1x}" "${h1y}" "${z}" >> "${output}"
        printf "H    %12.6f  %12.6f  %12.6f\n" "${h2x}" "${h2y}" "${z}" >> "${output}"
    done

    echo "${output}"
}

# ---- 主测试流程 ----
echo "=== 生成测试输入文件 ==="
echo ""

RESULTS_FILE="speed_comparison.dat"
echo "# N_molecules  N_atoms  DFT_time(s)  xTB_time(s)  Speedup" > "${RESULTS_FILE}"

for n_mol in "${SYSTEM_SIZES[@]}"; do
    n_atoms=$((n_mol * 3))
    echo "--- 测试 ${n_mol} 个水分子 (${n_atoms} 原子) ---"

    # 生成坐标文件（如果不存在）
    if [ ! -f "water${n_mol}.xyz" ]; then
        echo "  生成坐标文件: water${n_mol}.xyz"
        generate_water_xyz "${n_mol}" > /dev/null
    fi

    # 生成输入文件
    dft_input=$(generate_dft_input "${n_mol}")
    xtb_input=$(generate_xtb_input "${n_mol}")
    echo "  DFT输入: ${dft_input}"
    echo "  xTB输入: ${xtb_input}"

    # 运行DFT测试（如果CP2K可用）
    dft_time="N/A"
    xtb_time="N/A"
    speedup="N/A"

    if command -v "${CP2K_BIN}" &> /dev/null; then
        echo "  运行DFT测试..."
        start_time=$(date +%s.%N)
        ${CP2K_BIN} -i "${dft_input}" -o "${dft_input%.inp}.out" 2>/dev/null || true
        end_time=$(date +%s.%N)
        dft_time=$(echo "${end_time} - ${start_time}" | bc -l)
        echo "  DFT时间: ${dft_time} 秒"

        echo "  运行xTB测试..."
        start_time=$(date +%s.%N)
        ${CP2K_BIN} -i "${xtb_input}" -o "${xtb_input%.inp}.out" 2>/dev/null || true
        end_time=$(date +%s.%N)
        xtb_time=$(echo "${end_time} - ${start_time}" | bc -l)
        echo "  xTB时间: ${xtb_time} 秒"

        if [ "${dft_time}" != "N/A" ] && [ "${xtb_time}" != "N/A" ] && [ "$(echo "${xtb_time} > 0" | bc -l)" -eq 1 ]; then
            speedup=$(echo "scale=1; ${dft_time} / ${xtb_time}" | bc -l)
            echo "  加速比: ${speedup}x"
        fi
    else
        echo "  (跳过实际运行 - CP2K未找到)"
    fi

    echo "${n_mol}  ${n_atoms}  ${dft_time}  ${xtb_time}  ${speedup}" >> "${RESULTS_FILE}"
    echo ""
done

echo "============================================"
echo " 比较结果摘要"
echo "============================================"
echo ""
cat "${RESULTS_FILE}"
echo ""

# ---- 生成绘图脚本 ----
cat > plot_speed.gp << 'GNUPLOT_EOF'
# 速度比较绘图脚本 (gnuplot)
# 用法: gnuplot plot_speed.gp

set terminal pngcairo size 1000,600 enhanced font 'Arial,14'
set output 'speed_comparison.png'

set title 'xTB vs DFT Computational Time' font 'Arial,16'
set xlabel 'Number of Water Molecules'
set ylabel 'Time (seconds)'
set grid
set key top left

set logscale xy

plot 'speed_comparison.dat' using 1:3 with linespoints title 'DFT (PBE/DZVP)' lw 2 pt 7 ps 1.5, \
     'speed_comparison.dat' using 1:4 with linespoints title 'GFN1-xTB' lw 2 pt 5 ps 1.5

set output 'speedup.png'
set title 'xTB Speedup over DFT' font 'Arial,16'
set ylabel 'Speedup Factor'
unset logscale y

plot 'speed_comparison.dat' using 1:5 with linespoints title 'Speedup (DFT/xTB)' lw 2 pt 9 ps 1.5 lc rgb "red"
GNUPLOT_EOF
echo "绘图脚本已生成: plot_speed.gp"
echo ""

# ---- 使用说明 ----
echo "============================================"
echo " 使用说明"
echo "============================================"
echo ""
echo "本脚本的工作流程："
echo "  1. 生成不同大小的水盒子坐标文件"
echo "  2. 生成对应的DFT和xTB输入文件"
echo "  3. 如果CP2K可用，运行并计时"
echo "  4. 计算加速比并保存结果"
echo ""
echo "手动运行方法："
echo "  # DFT计算"
echo "  cp2k.psmp -i water64_dft.inp -o water64_dft.out"
echo ""
echo "  # xTB计算"
echo "  cp2k.psmp -i water64_xtb.inp -o water64_xtb.out"
echo ""
echo "  # 比较输出文件中的'CP2K| Total elapsed time'"
echo "  grep 'Total elapsed time' water64_dft.out water64_xtb.out"
echo ""
echo "注意事项："
echo "  - 确保BASIS_MOLOPT和GTH_POTENTIALS文件在CP2K路径中"
echo "  - 测试时只运行少量步数（5步MD）以节省时间"
echo "  - 实际加速比取决于系统大小和硬件"
echo "  - 对于小系统（<50原子），DFT已经很快，xTB优势不明显"
echo "  - 对于大系统（>500原子），xTB的优势最为显著"
