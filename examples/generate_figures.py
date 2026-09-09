import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import numpy as np, json, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from src.config import A, C
from src.simulate import sim_median, sim_stubborn, burst

BLUE, ORANGE, PURPLE = "#2563a8", "#d97706", "#6b4fa0"
plt.rcParams.update({'font.size': 6.6, 'font.family': 'serif', 'axes.linewidth': .5,
    'xtick.major.width': .5, 'ytick.major.width': .5, 'xtick.major.size': 2, 'ytick.major.size': 2,
    'axes.spines.top': False, 'axes.spines.right': False, 'lines.linewidth': 1.0,
    'xtick.major.pad': 1.5, 'ytick.major.pad': 1.5})

def main():
    L1 = np.load('Lopt1.npy')
    L2 = np.load('Lopt2.npy')

    e0 = np.array([1., -1., .8])
    h1 = [e0] + [np.zeros(3)] * 2
    h2 = [e0] + [np.zeros(3)] * 4
    fl = lambda v: np.maximum(v, 1e-16)

    fig, ax = plt.subplots(2, 2, figsize=(3.45, 2.02))
    LG = dict(fontsize=5.6, handlelength=1.5, borderpad=.12, labelspacing=.2, handletextpad=.4,
              frameon=True, facecolor='white', edgecolor='none', framealpha=.85)

    T = 250
    k = np.arange(T)
    d31 = burst(T, 1, 3, 50., np.random.default_rng(3))
    d51 = burst(T, 1, 5, 50., np.random.default_rng(3))

    a = ax[0, 0]
    a.semilogy(k, fl(sim_median(A, C, L1, 1, d31, h1, T)), color=BLUE, label=r'median $q\!=\!1$')
    a.semilogy(k, fl(sim_stubborn(A, C, L1, 3, 0.8, d31, e0, T)), color=ORANGE, ls='-.', label='stubborn')
    a.semilogy(k, fl(sim_stubborn(A, C, L1, 3, 0.8, d51, e0, T)), color=ORANGE, ls=':', label=r'stubborn $P\!=\!5$')
    a.set_title(r'(a) $b=1$, $P=3$', fontsize=6.6, pad=2)
    a.legend(loc='lower left', **LG)

    T2 = 600
    k2 = np.arange(T2)
    d25 = burst(T2, 2, 5, 50., np.random.default_rng(3))

    a = ax[0, 1]
    a.semilogy(k2, fl(sim_median(A, C, L2, 2, d25, h2, T2)), color=BLUE, label=r'median $q\!=\!2$')
    a.semilogy(k2, fl(sim_median(A, C, L1, 1, d25, h1, T2)), color=PURPLE, ls='--', label=r'median $q\!=\!1$')
    a.semilogy(k2, fl(sim_stubborn(A, C, L1, 3, 0.8, d25, e0, T2)), color=ORANGE, ls='-.', label='stubborn')
    a.set_title(r'(b) $b=2$, $P=5$', fontsize=6.6, pad=2)
    a.legend(loc='lower left', **LG)

    for a_ in ax[0]:
        a_.set_xlabel(r'$k$', labelpad=0)
        a_.set_ylim(1e-12, 3e2)
        a_.set_ylabel(r'$|e_k|$', labelpad=0)
        a_.grid(True, lw=.25, color='#dcdcd8')
        a_.set_axisbelow(True)
        a_.set_yticks([1e-12, 1e-8, 1e-4, 1e0])

    with open('fig2c.json') as f:
        D = json.load(f)
    qs = np.array(D['q'])

    a = ax[1, 0]
    for av, c, mk in [(0.5, BLUE, 'o'), (0.7, PURPLE, 's'), (0.9, ORANGE, '^')]:
        y = np.array([np.nan if v is None else v for v in D['a%s' % av]])
        a.plot(qs, y, color=c, marker=mk, ms=2.6, label=r'$a=%.1f$' % av)
        a.axhline(av, color=c, lw=.5, ls='--')
    a.set_xlabel(r'$q$', labelpad=0)
    a.set_ylabel(r'$\alpha^\star$', labelpad=0)
    a.set_title(r'(c) certified rate', fontsize=6.6, pad=2)
    a.set_ylim(0.4, 1.0)
    a.legend(loc='upper left', **LG)

    a = ax[1, 1]
    ceil = np.array(D['ceil'])
    a.plot(qs, 1 + 1 / (2 * qs), color=PURPLE, ls='--', marker='s', ms=2.6, label=r'bound $1+\frac{1}{2q}$')
    a.plot(qs, ceil, color=BLUE, marker='o', ms=2.6, label=r'LMI (22)')
    a.axhline(1.0, color='#999', lw=.5)
    a.set_xlabel(r'$q$', labelpad=0)
    a.set_ylabel(r'$\max\,a$', labelpad=0)
    a.set_title(r'(d) feasibility ceiling', fontsize=6.6, pad=2)
    a.legend(loc='upper right', **LG)

    for a_ in ax[1]:
        a_.grid(True, lw=.25, color='#dcdcd8')
        a_.set_axisbelow(True)
        a_.set_xticks(qs)

    fig.tight_layout(pad=.2, w_pad=.9, h_pad=.7)
    import os
    os.makedirs('Fig', exist_ok=True)
    fig.savefig('Fig/all.pdf', bbox_inches='tight', pad_inches=0.01)
    print('Fig/all.pdf saved successfully.')

if __name__ == "__main__":
    main()
