"""
Figures of the technical notes on LiF (Notes/main.tex), made from the runs in NL_Chi, which are on the cluster only.

Usage, on ismhpc from ~/work/LiF/Notes:

    ~/miniconda3/bin/python make_figures.py [figure ...] [--out DIR]

Without figure names all the figures are made. The jpeg files are written in DIR (default: Figures). To avoid
untracked files in the cluster copy of the repository, write them outside it (e.g. --out ~/tmp_claude/notes_figures)
and copy them to Notes/Figures on the laptop, where they are committed.
"""
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mppi import Parsers as P
from mppi.Utilities import Constants as Const

HERE = os.path.dirname(os.path.abspath(__file__))
NL_CHI = os.path.join(HERE, '..', 'NL_Chi')
RT_DIR = os.path.join(NL_CHI, 'kx16_nb100_NoTr_E100')    # RT pump-probe runs (YamboNL_Analysis.ipynb)
T_PROBES = [13.0, 14.5, 16.0, 17.5, 19.0]                  # probe times of the RT runs (fs)

plt.rcParams.update({'font.size': 9, 'axes.labelsize': 9, 'legend.fontsize': 8, 'savefig.dpi': 200,
                     'savefig.bbox': 'tight'})


def rt_id(bands, pump_int=None, t_probe=None, neg=False):
    # run id of the RT runs (bands as a string, e.g. '3-6'; neg: pump along -x), as in YamboNL_Analysis.ipynb
    idd = 'rt-ip-bands_%s' % bands
    if pump_int is not None:
        idd += '-EP_%s' % ('%.0e' % pump_int).replace('e+0', 'e') + ('-neg' if neg else '')
    if t_probe is not None:
        idd += '-tp_%.2f' % t_probe
    return idd + '-' + {(True, True): 'Pp', (True, False): 'P', (False, True): 'p'}[(pump_int is not None,
                                                                                       t_probe is not None)]


def load_nl(run_id):
    return P.YamboNLDBParser(os.path.join(RT_DIR, run_id, 'ndb.Nonlinear'), verbose=False)


def pump_center(run):
    # time of the maximum of the QSFIELD 5SIGMA pump envelope and its width sigma (fs), as in yambo RT_small_a
    pump = run.Efield2[0]
    wP, sigma = pump['freq_range'][0], pump['damping'] / (2 * np.sqrt(2 * np.log(2)))   # au, damping = FWHM
    t_P = pump['initial_time'] + np.pi / wP * np.rint(wP / np.pi * 5 * sigma)
    return t_P / Const.FsToAu, sigma / Const.FsToAu


def fig_rt_setup(out):
    """
    Pump field of the RT runs (6e7 kW/cm^2, bands 3-6), its envelope and the probe times.
    """
    run = load_nl(rt_id('3-6', 6e7))                                  # pump-only run (pump as Field1)
    t = np.array(run.IO_TIME_points) / Const.FsToAu
    t_P, sigma = pump_center(run)
    E0P = run.Efield2[0]['amplitude']
    f = np.exp(-(t - t_P)**2 / (2 * sigma**2))
    t_p = [load_nl(rt_id('3-6', 6e7, tp)).Efield[0]['initial_time'] / Const.FsToAu for tp in T_PROBES]

    fig, ax = plt.subplots(figsize=(6.5, 2.6))
    ax.plot(t, np.array(run.E_ext[0])[0].real / E0P, 'k', lw=0.8, label=r'pump field $E_{\mathrm{P}}(t)/E_{0\mathrm{P}}$')
    ax.plot(t, f, 'b--', lw=1, label=r'envelope $\pm f(t-t_{\mathrm{P}})$')
    ax.plot(t, -f, 'b--', lw=1)
    for tp in t_p:
        ax.axvline(tp, color='r', lw=1)
    ax.axvline(t_P, color='gray', ls=':', lw=1, label=r'pump maximum $t_{\mathrm{P}}$ = %.2f fs' % t_P)
    ax.plot([], [], 'r', lw=1, label=r'probe times $t_{\mathrm{p}}$')
    ax.set_xlim(4, 28)
    ax.set_ylim(-1.15, 1.15)
    ax.set_xlabel('time (fs)')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False)
    top = ax.twiny()                                                  # delays tau = t_p - t_P on the top axis
    top.set_xlim(ax.get_xlim())
    top.set_xticks(t_p)
    top.set_xticklabels(['%+.1f' % (tp - t_P) for tp in t_p], color='r')
    top.set_xlabel(r'delay $\tau = t_{\mathrm{p}}-t_{\mathrm{P}}$ (fs)', color='r')
    fig.savefig(os.path.join(out, 'rt_setup.jpeg'))
    plt.close(fig)
    print('rt_setup: t_P = %.3f fs, sigma = %.3f fs, tau =' % (t_P, sigma), np.round(np.array(t_p) - t_P, 3))


FIGURES = {'rt_setup': fig_rt_setup}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Figures of the technical notes on LiF')
    parser.add_argument('figures', nargs='*', help='names of the figures (default: all): %s' % ', '.join(FIGURES))
    parser.add_argument('--out', default=os.path.join(HERE, 'Figures'), help='output folder of the jpeg files')
    args = parser.parse_args()
    out = os.path.expanduser(args.out)
    os.makedirs(out, exist_ok=True)
    for name in args.figures or FIGURES:
        FIGURES[name](out)
