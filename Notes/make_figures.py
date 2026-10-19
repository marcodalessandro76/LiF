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
from mppi import Optics as O
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


def rt_eps(run, ref=None, eta=0.6):
    # energy (eV) and eps along x from the delta probe of the run, minus the polarization of ref (pump-only run), with a
    # Lorentzian smoothing eta (eV)
    return O.Linear_Response(np.array(run.IO_TIME_points), run.Polarization[0], run.Efield[0],
                             pol_ref=None if ref is None else ref.Polarization[0], eta=eta)


def fig_rt_even_odd(out, t_probe=16.0, pump_int=6e7, erange=(10, 25)):
    """
    Even and odd parts in the pump field of the RT Delta eps (bands 3-6) at one probe time, from the runs with the pump
    along +x and -x.
    """
    E, eps_eq = rt_eps(load_nl(rt_id('3-6', None, 10.0)))                 # probe-only run
    ref_p, ref_m = load_nl(rt_id('3-6', pump_int)), load_nl(rt_id('3-6', pump_int, neg=True))
    run_p, run_m = load_nl(rt_id('3-6', pump_int, t_probe)), load_nl(rt_id('3-6', pump_int, t_probe, neg=True))
    eps_p, eps_m = rt_eps(run_p, ref_p)[1], rt_eps(run_m, ref_m)[1]
    even, odd = (eps_p + eps_m) / 2 - eps_eq, (eps_p - eps_m) / 2
    tau = run_p.Efield[0]['initial_time'] / Const.FsToAu - pump_center(ref_p)[0]
    s = (E >= erange[0]) & (E <= erange[1])

    fig, axes = plt.subplots(1, 2, figsize=(7, 2.6), sharex=True)
    for ax, part, lab in zip(axes, [np.real, np.imag], ['Re', 'Im']):
        ax.plot(E[s], 1e3 * part(even[s]), 'k', lw=1.2, label='even part (third order)')
        ax.plot(E[s], 1e3 * part(odd[s]), 'r--', lw=1.2, label='odd part (second order)')
        ax.axhline(0, color='gray', lw=0.5)
        ax.set_xlabel('probe energy (eV)')
        ax.set_ylabel(r'%s $\Delta\varepsilon$ ($10^{-3}$)' % lab)
        ax.set_xlim(erange)
        ax.set_xticks(np.arange(erange[0], erange[1] + 0.1, 5))
    fig.suptitle(r'pump $%s\cdot10^{%d}$ kW/cm$^2$, $\tau$ = %+.2f fs' % (('%.0e' % pump_int)[0],
                 int(np.log10(pump_int)), tau), fontsize=9, y=0.97)
    fig.tight_layout()
    fig.legend(*axes[0].get_legend_handles_labels(), loc='upper center', bbox_to_anchor=(0.5, 0.02), ncol=2,
               frameon=False)
    fig.savefig(os.path.join(out, 'rt_even_odd.jpeg'))
    plt.close(fig)
    print('rt_even_odd: tau = %+.2f fs, max|even| %.2e, max|odd| %.2e' % (tau, np.abs(even[s]).max(),
                                                                          np.abs(odd[s]).max()))


GRIDS = {'kx8': 'kx8_nb100_NoTr_E100', 'kx12': 'kx12_nb100_NoTr_E100', 'kx16': 'kx16_nb100_NoTr_E100'}
NAME_SINE = 'pulse-E1_1e3-nlenrange_10.0-25.016-nlensteps_155-bands_3-6-damp_0.3-sin-nltime_100'
NAME_PP = 'Pp-Ep_1e3-EP_1e6-PE_1.55-nlenrange_10.0-25.016-nlensteps_155-bands_3-6-damp_0.3-sin-nltime_100'
MIX_KEYS = [(1, 0), (1, 1), (1, -1), (1, 2), (1, -2)]


def chi_kgrids(out, broadening=0.3):
    """
    chi1 (sine runs) and the mixing keys (P&p runs, pump 1e6 kW/cm^2) of the three k grids, with the a posteriori
    broadening (eV) added to the damping 0.3 eV of the dynamics. Cached in out/chi_kgrids_<broadening>.npz.
    """
    cache = os.path.join(out, 'chi_kgrids_%.2f.npz' % broadening)
    if os.path.exists(cache):
        return dict(np.load(cache))
    res = {}
    for g, folder in GRIDS.items():
        d_sine = P.YamboNLDBParser(os.path.join(NL_CHI, folder, NAME_SINE, 'ndb.Nonlinear'), verbose=False)
        d_pp = P.YamboNLDBParser(os.path.join(NL_CHI, folder, NAME_PP, 'ndb.Nonlinear'), verbose=False)
        fm = O.Xn_frequency_mixing(d_pp, X_order=(1, 2), verbose=False)
        mix = fm.compute_Xn(broadening=broadening)[0]
        res[g + '_chi1'] = O.Xn_single_frequency(d_sine, X_order=3, verbose=False).compute_Xn(broadening=broadening)[0][1]
        for key in MIX_KEYS:
            res[g + '_%d_%d' % key] = mix[key]
        res[g + '_EP'] = d_pp.Efield2[0]['amplitude']
        res['freqs'] = fm.probe_freqs * Const.HaToeV
        print('chi_kgrids: %s done' % g, flush=True)
    np.savez(cache, **res)
    return res


def fig_chi_kgrid(out, erange=(10, 25)):
    """
    Dependence on the k grid of the spurious second order key chi_(1,1) and of the third order term chi_0 at the probe
    frequency (total broadening 0.6 eV).
    """
    r = chi_kgrids(out)
    f = r['freqs']
    s = (f >= erange[0]) & (f <= erange[1])
    w = (f >= 10) & (f <= 20)
    chi0 = {g: (r[g + '_1_0'] - r[g + '_chi1']) / (r[g + '_EP']**2 / 4) for g in GRIDS}
    style = {'kx8': ('g', ':'), 'kx12': ('b', '--'), 'kx16': ('k', '-')}

    fig, axes = plt.subplots(1, 2, figsize=(7, 2.7))
    for g in GRIDS:
        c, ls = style[g]
        axes[0].semilogy(f[s], np.abs(r[g + '_1_1'][s]), color=c, ls=ls, lw=1.2, label=g)
        axes[1].plot(f[s], chi0[g][s].imag, color=c, ls=ls, lw=1.2, label=g)
    axes[0].set_ylabel(r'$|\chi_{(1,1)}|$ (au)')
    axes[1].set_ylabel(r'Im $\chi^{(3)}_0$ (au)')
    axes[1].axhline(0, color='gray', lw=0.5)
    for ax in axes:
        ax.set_xlabel('probe energy (eV)')
        ax.set_xlim(erange)
        ax.set_xticks(np.arange(erange[0], erange[1] + 0.1, 5))
    fig.tight_layout()
    fig.legend(*axes[1].get_legend_handles_labels(), loc='upper center', bbox_to_anchor=(0.5, 0.02), ncol=3,
               frameon=False)
    fig.savefig(os.path.join(out, 'chi_kgrid.jpeg'))
    plt.close(fig)

    print('chi_kgrid: medians of |chi_(1,+-1)| in 10-20 eV (au) and ratios between the grids')
    for key in [(1, 1), (1, -1)]:
        med = {g: np.median(np.abs(r[g + '_%d_%d' % key][w])) for g in GRIDS}
        print('  %-7s' % str(key), ' '.join('%s %.2e' % (g, m) for g, m in med.items()),
              '  kx8/kx12 %.2f  kx12/kx16 %.2f' % (med['kx8'] / med['kx12'], med['kx12'] / med['kx16']))
    print('chi_kgrid: chi_0, |chi_g - chi_kx16|/max|chi_kx16| (mean/max) and ratio of the maxima, below/above 16 eV')
    for g in ['kx8', 'kx12']:
        for lab, sel in [('10-16', (f >= 10) & (f < 16)), ('16-21', (f >= 16) & (f <= 21))]:
            d = np.abs(chi0[g][sel] - chi0['kx16'][sel]) / np.abs(chi0['kx16'][sel]).max()
            print('  %-5s %s eV: %.2f / %.2f, ratio of the maxima %.2f' % (g, lab, d.mean(), d.max(),
                  np.abs(chi0[g][sel]).max() / np.abs(chi0['kx16'][sel]).max()))


FIGURES = {'rt_setup': fig_rt_setup, 'rt_even_odd': fig_rt_even_odd, 'chi_kgrid': fig_chi_kgrid}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Figures of the technical notes on LiF')
    parser.add_argument('figures', nargs='*', help='names of the figures (default: all): %s' % ', '.join(FIGURES))
    parser.add_argument('--out', default=os.path.join(HERE, 'Figures'), help='output folder of the jpeg files')
    args = parser.parse_args()
    out = os.path.expanduser(args.out)
    os.makedirs(out, exist_ok=True)
    for name in args.figures or FIGURES:
        FIGURES[name](out)
