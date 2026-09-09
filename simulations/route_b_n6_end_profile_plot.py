"""Render the measured end-profile continuation and its perturbative predictions."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def render():
    data = json.loads((ROOT/'simulations/results/route_b_n6_end_profile_probe.json').read_text())
    colors = {'one': '#1679a5', 'even': '#db7235', 'odd': '#8356a2'}
    labels = {'one': 'One end', 'even': 'Equal ends (reflection kept)', 'odd': 'Opposite ends'}
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fig, (ax, scaling) = plt.subplots(1, 2, figsize=(12.4, 5.7))
    q0 = complex(data['q0']['real'], data['q0']['imag'])
    for name, profile in data['profiles'].items():
        color = colors[name]
        for branch, style, marker in ((0, '-', 'o'), (1, '--', 's')):
            rows = [r for r in profile['rows'] if r['branch'] == branch]
            qs = [q0]+[complex(r['q']['real'], r['q']['imag']) for r in rows]
            ax.plot(np.real(qs), np.imag(qs), style, color=color, lw=1.9)
            ax.scatter(qs[-1].real, qs[-1].imag, marker=marker, color=color, s=35, zorder=4)
        x = np.array(profile['scaling']['epsilon'])
        y = np.array(profile['scaling']['q_branch_separation'])
        scaling.loglog(x, y, 'o-', color=color, markersize=4, lw=1.7)
        slopes = np.array(profile['leading_slopes']['real'])+1j*np.array(profile['leading_slopes']['imag'])
        if name != 'even':
            prediction = abs(slopes[0]-slopes[1])*x**profile['order']
            scaling.loglog(x, prediction, ':', color='#555555', lw=1.3)
    ax.scatter(q0.real, q0.imag, marker='*', color='#20242b', s=125, zorder=5)
    ax.annotate('Original diabolic seed', (q0.real, q0.imag), xytext=(-185, -3),
                textcoords='offset points', fontsize=9)
    ax.set(xlabel='Re q', ylabel='Im q', title='Two candidate paths per profile, up to ε = 0.05')
    ax.legend(handles=[Line2D([0], [0], color='#555555', ls='-', marker='o', label='Branch 1'),
                       Line2D([0], [0], color='#555555', ls='--', marker='s', label='Branch 2')],
              loc='lower right', fontsize=9)
    scaling.set(xlabel='End-profile amplitude ε', ylabel='Separation |q₁ − q₂|',
                title='Linear versus quadratic opening')
    scaling.text(0.00014, 0.0013, '∝ ε', color='#555555', fontsize=12)
    scaling.text(0.0055, 0.00006, '∝ ε²', color='#555555', fontsize=12)
    scaling.legend(handles=[Line2D([0], [0], color='#555555', ls=':',
                label='Leading perturbation prediction')], loc='upper left', fontsize=9)
    for axis in (ax, scaling):
        axis.grid(alpha=0.17)
        axis.spines[['top', 'right']].set_visible(False)
    fig.suptitle('Reflection symmetry does not preserve this crossing', fontsize=16, y=0.98)
    fig.text(0.5, 0.91, 'N6-E-A2-T-007 · full 90D block · complex q · numerical continuation',
             ha='center', color='#555555')
    fig.legend(handles=[Line2D([0], [0], color=colors[n], lw=2, label=labels[n])
                        for n in colors], loc='lower center', ncol=3, frameon=False)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.81, bottom=0.17, wspace=0.29)
    for suffix in ('png', 'svg'):
        fig.savefig(ROOT/f'visualizations/route_b_n6_end_profiles.{suffix}', dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    render()
