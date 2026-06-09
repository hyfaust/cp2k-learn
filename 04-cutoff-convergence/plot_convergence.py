#!/usr/bin/env python3
"""
plot_convergence.py — 绘制 CUTOFF / REL_CUTOFF 收敛曲线

用法:
    python plot_convergence.py
    python plot_convergence.py --data cutoff_results.dat
    python plot_convergence.py --data rel_cutoff_results.dat --title "REL_CUTOFF 收敛" --rel-cutoff

功能:
    - 读取 scan_cutoff.sh 或 scan_rel_cutoff.sh 生成的 .dat 文件
    - 绘制能量随参数变化的收敛曲线
    - 绘制能量差 (ΔE) 曲线，标注收敛阈值线
    - 保存为 PNG 图片
"""

import argparse
import sys
import os

try:
    import matplotlib
    matplotlib.use('Agg')  # 无显示器环境下使用非交互后端
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


def read_data(filename):
    """读取 .dat 结果文件，返回 (参数列表, 能量列表)"""
    params = []
    energies = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                params.append(float(parts[0]))
                energies.append(float(parts[1]))
    return params, energies


def compute_energy_deltas(params, energies):
    """计算相邻能量差 ΔE (Hartree 和 eV)"""
    deltas_ha = [0.0]  # 第一个值没有差值
    deltas_ev = [0.0]
    for i in range(1, len(energies)):
        d = energies[i] - energies[i - 1]
        deltas_ha.append(d)
        deltas_ev.append(d * 27.211386)  # Hartree -> eV
    return deltas_ha, deltas_ev


def plot_with_matplotlib(params, energies, deltas_ha, deltas_ev, 
                         param_label, title, output_prefix):
    """使用 matplotlib 绘图"""
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle(title, fontsize=14, fontweight='bold')

    # --- 子图 1: 总能量 ---
    ax1 = axes[0]
    ax1.plot(params, energies, 'bo-', markersize=8, linewidth=2, label='Total Energy')
    ax1.set_ylabel('Energy (Hartree)', fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10)

    # 在每个点上标注能量值
    for x, y in zip(params, energies):
        ax1.annotate(f'{y:.6f}', (x, y), textcoords="offset points",
                     xytext=(0, 12), ha='center', fontsize=8, color='blue')

    # --- 子图 2: 能量差 ΔE ---
    ax2 = axes[1]
    # 跳过第一个点（ΔE = 0）
    ax2.semilogy(params[1:], [abs(d) for d in deltas_ev[1:]], 
                 'rs-', markersize=8, linewidth=2, label='|ΔE| (eV)')
    
    # 收敛阈值线: 1 meV = 0.001 eV
    threshold = 0.001
    ax2.axhline(y=threshold, color='green', linestyle='--', linewidth=1.5, 
                label=f'Threshold = {threshold} eV')
    # 收敛阈值线: 10^-6 Hartree ≈ 0.0272 meV
    threshold_ha = 1e-6 * 27.211386  # 转换为 eV
    ax2.axhline(y=threshold_ha, color='orange', linestyle=':', linewidth=1.5,
                label=f'Threshold = 10⁻⁶ Ha ({threshold_ha:.4f} eV)')

    ax2.set_xlabel(param_label, fontsize=12)
    ax2.set_ylabel('|ΔE| (eV)', fontsize=12)
    ax2.grid(True, alpha=0.3, which='both')
    ax2.legend(fontsize=10)

    plt.tight_layout()
    
    filename = f"{output_prefix}.png"
    plt.savefig(filename, dpi=150, bbox_inches='tight')
    print(f"图片已保存: {filename}")


def plot_text(params, energies, deltas_ha, deltas_ev, param_label):
    """无 matplotlib 时的纯文本输出"""
    print("")
    print(f"{'='*70}")
    print(f"  收敛性分析 ({param_label})")
    print(f"{'='*70}")
    print(f"  {param_label:<14} {'Energy (Ha)':>18} {'ΔE (Ha)':>14} {'ΔE (eV)':>14} {'ΔE (meV)':>14}")
    print(f"  {'-'*70}")

    for i in range(len(params)):
        d_ha = deltas_ha[i]
        d_ev = deltas_ev[i]
        d_mev = d_ev * 1000
        delta_str = f"{d_ha:>14.10f}" if i > 0 else f"{'---':>14}"
        delta_ev_str = f"{d_ev:>14.8f}" if i > 0 else f"{'---':>14}"
        delta_mev_str = f"{d_mev:>14.6f}" if i > 0 else f"{'---':>14}"
        print(f"  {params[i]:<14.1f} {energies[i]:>18.10f} {delta_str} {delta_ev_str} {delta_mev_str}")

    print(f"  {'-'*70}")
    
    # 判断收敛点
    converged_idx = None
    for i in range(len(deltas_ev) - 1, 0, -1):
        if abs(deltas_ev[i]) > 0.001:  # 1 meV threshold
            converged_idx = i + 1
            break
    
    if converged_idx is not None and converged_idx < len(params):
        print(f"\n  收敛判断 (阈值 < 1 meV):")
        print(f"  参数 {param_label} = {params[converged_idx]:.1f} 已收敛")
    else:
        print(f"\n  在当前范围内未找到明确收敛点，请扩大扫描范围。")


def main():
    parser = argparse.ArgumentParser(
        description='绘制 CP2K CUTOFF / REL_CUTOFF 收敛曲线')
    parser.add_argument('--data', '-d', default='cutoff_results.dat',
                        help='结果数据文件 (默认: cutoff_results.dat)')
    parser.add_argument('--title', '-t', default='CUTOFF 收敛性测试',
                        help='图表标题')
    parser.add_argument('--rel-cutoff', action='store_true',
                        help='标记为 REL_CUTOFF 扫描（修改 x 轴标签）')
    parser.add_argument('--output', '-o', default=None,
                        help='输出文件前缀（默认根据输入文件名生成）')
    args = parser.parse_args()

    # 检查数据文件是否存在
    if not os.path.isfile(args.data):
        print(f"[错误] 数据文件 '{args.data}' 不存在！")
        print("请先运行 scan_cutoff.sh 或 scan_rel_cutoff.sh 生成数据文件。")
        sys.exit(1)

    # 读取数据
    params, energies = read_data(args.data)
    
    if len(params) == 0:
        print("[错误] 数据文件中没有有效数据！")
        sys.exit(1)

    print(f"读取到 {len(params)} 个数据点")

    # 计算能量差
    deltas_ha, deltas_ev = compute_energy_deltas(params, energies)

    # 确定标签和输出前缀
    if args.rel_cutoff:
        param_label = 'REL_CUTOFF (Ry)'
        output_prefix = args.output or 'rel_cutoff_convergence'
    else:
        param_label = 'CUTOFF (Ry)'
        output_prefix = args.output or 'cutoff_convergence'

    # 文本输出（始终显示）
    plot_text(params, energies, deltas_ha, deltas_ev, param_label)

    # 绘图
    if HAS_MATPLOTLIB:
        plot_with_matplotlib(params, energies, deltas_ha, deltas_ev,
                             param_label, args.title, output_prefix)
    else:
        print("\n[提示] 未安装 matplotlib，无法生成图片。")
        print("安装方法: pip install matplotlib")
        print("文本结果已输出到终端。")


if __name__ == '__main__':
    main()
