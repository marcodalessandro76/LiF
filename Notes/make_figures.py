"""
Figures of the technical notes on LiF (Notes/main.tex), made from the runs in NL_Chi, which are on the cluster only.

Usage, on ismhpc from ~/work/LiF/Notes:

    ~/miniconda3/bin/python make_figures.py [figure ...] [--out DIR]

Without figure names all the figures are made. The jpeg files are written in DIR (default: Figures). To avoid
untracked files in the cluster copy of the repository, write them outside it (e.g. --out ~/tmp_claude/notes_figures)
and copy them to Notes/Figures on the laptop, where they are committed.
"""
import os
import glob
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scipy.optimize import minimize_scalar
from scipy.interpolate import CubicSpline
from mppi import Parsers as P
from mppi import Optics as O
from mppi.Optics.Utils import time_weighted
from mppi.Utilities import Constants as Const

HERE = os.path.dirname(os.path.abspath(__file__))
NL_CHI = os.environ.get('LIF_NL_CHI', os.path.join(HERE, '..', 'NL_Chi'))   # data folder (env variable for tests)
RT_DIR = os.path.join(NL_CHI, 'kx16_nb100_NoTr_E100')    # RT pump-probe runs (YamboNL_Analysis.ipynb)
T_PROBES = [13.0, 14.5, 16.0, 17.5, 19.0]                  # probe times of the RT runs (fs)

plt.rcParams.update({'font.size': 9, 'axes.labelsize': 9, 'legend.fontsize': 8, 'savefig.dpi': 200,
                     'savefig.bbox': 'tight'})


def rt_id(bands, pump_int=None, t_probe=None, neg=False, level='ip'):
    # run id of the RT runs (bands as a string, e.g. '3-6'; neg: pump along -x; level 'ip' or 'hsex' (HARTREE+SEX)), as
    # in YamboNL_Analysis.ipynb
    idd = 'rt-%s-bands_%s' % (level, bands)
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
        axes[0].plot(f[s], np.abs(r[g + '_1_1'][s]), color=c, ls=ls, lw=1.2, label=g)
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


SCAN39 = 'nlenrange_10.0-25.113-nlensteps_39'                  # 39 frequencies of the test grid (step wP/4)
NAMES_PUMP_KX12 = {1e6: NAME_PP,
                   4e6: 'Pp-Ep_1e3-EP_4e6-PE_1.55-nlenrange_10.0-25.016-nlensteps_155-bands_3-6-damp_0.3-sin-nltime_100',
                   1.6e7: 'Pp-Ep_1e3-EP_1.6e7-PE_1.55-%s-bands_3-6-damp_0.3-sin-nltime_100' % SCAN39,
                   6.4e7: 'Pp-Ep_1e3-EP_6.4e7-PE_1.55-%s-bands_3-6-damp_0.3-sin-nltime_100' % SCAN39,
                   2.56e8: 'Pp-Ep_1e3-EP_2.56e8-PE_1.55-%s-bands_3-6-damp_0.3-sin-nltime_100' % SCAN39}


def chi_pumps_kx12(out):
    """
    Third order terms extracted from the kx12 P&p runs with the pump intensities of NAMES_PUMP_KX12 (no a posteriori
    broadening: the runs above 4e6 kW/cm^2 have only 39 frequencies). The runs not yet done are skipped. Cached in
    out/chi_pumps_kx12.npz: the intensities already in the cache are not computed again.
    """
    cache = os.path.join(out, 'chi_pumps_kx12.npz')
    res = dict(np.load(cache)) if os.path.exists(cache) else {}
    folder = os.path.join(NL_CHI, GRIDS['kx12'])
    todo = {I: name for I, name in NAMES_PUMP_KX12.items() if '%.1e_chi0' % I not in res
            and os.path.exists(os.path.join(folder, name, 'ndb.Nonlinear'))}
    if not todo:
        return res
    d_sine = P.YamboNLDBParser(os.path.join(folder, NAME_SINE, 'ndb.Nonlinear'), verbose=False)
    sf = O.Xn_single_frequency(d_sine, X_order=3, verbose=False)
    chi1, f1 = sf.compute_Xn()[0][1], sf.fields_freqs * Const.HaToeV
    for I, name in todo.items():
        d = P.YamboNLDBParser(os.path.join(folder, name, 'ndb.Nonlinear'), verbose=False)
        fm = O.Xn_frequency_mixing(d, X_order=(1, 2), verbose=False)
        mix = fm.compute_Xn()[0]
        f = fm.probe_freqs * Const.HaToeV
        idx = np.array([np.argmin(np.abs(f1 - x)) for x in f])
        assert np.allclose(f1[idx], f, atol=1e-6)
        EP = d.Efield2[0]['amplitude']
        tag = '%.1e' % I
        res[tag + '_freqs'] = f
        res[tag + '_chi0'] = (mix[(1, 0)] - chi1[idx]) / (EP**2 / 4)
        res[tag + '_chip'] = mix[(1, 2)]
        res[tag + '_chim'] = mix[(1, -2)]
        print('chi_pumps_kx12: pump %.1e done (E0P %.3e au, %d frequencies)' % (I, EP, len(f)), flush=True)
    np.savez(cache, **res)
    return res


def fig_chi_pump(out, erange=(10, 25)):
    """
    Third order term chi_0 extracted with different pump intensities (kx12, damping 0.3 eV) and its change w.r.t. the
    lowest intensity: the fifth order terms, projected on chi_0, change it proportionally to the intensity.
    """
    r = chi_pumps_kx12(out)
    ref, f = r['1.0e+06_chi0'], r['1.0e+06_freqs']
    w = (f >= 10) & (f <= 20)
    norm = np.abs(ref[w]).max()
    style = {'1.0e+06': ('k', '-', None), '4.0e+06': ('b', '--', None), '1.6e+07': ('r', 'none', 'o'),
             '6.4e+07': ('m', 'none', 's'), '2.6e+08': ('g', 'none', '^')}
    label = {'1.0e+06': r'$10^6$', '4.0e+06': r'$4\cdot10^6$', '1.6e+07': r'$1.6\cdot10^7$',
             '6.4e+07': r'$6.4\cdot10^7$', '2.6e+08': r'$2.56\cdot10^8$'}

    fig, axes = plt.subplots(1, 2, figsize=(7, 2.7))
    for tag, (c, ls, m) in style.items():
        if tag + '_chi0' not in r:
            continue
        x, y = r[tag + '_freqs'], r[tag + '_chi0']
        s = (x >= erange[0]) & (x <= erange[1])
        axes[0].plot(x[s], y[s].imag, color=c, ls=ls, marker=m, ms=3, lw=1.2, label=label[tag] + ' kW/cm$^2$')
        if tag != '1.0e+06':
            idx = np.array([np.argmin(np.abs(f - xx)) for xx in x])
            d = 100 * np.abs(y - ref[idx]) / norm
            axes[1].plot(x[s], d[s], color=c, ls='-' if m else ls, marker=m, ms=3, lw=1)
            ww = (x >= 10) & (x <= 20)
            rel = np.abs(y - ref[idx]) / np.abs(ref[idx])
            print('chi_pump: %s, chi_0 change w.r.t. 1e6 in 10-20 eV: mean %.2f%% of the max, median relative %.2f%%'
                  % (tag, d[ww].mean(), 100 * np.median(rel[ww])))
            for key in ['chip', 'chim']:
                rk = np.abs(r[tag + '_' + key] - r['1.0e+06_' + key][idx]) / np.abs(r['1.0e+06_' + key][idx])
                print('          %s: median relative change %.2f%%' % (key, 100 * np.median(rk[ww])))
    axes[0].set_ylabel(r'Im $\chi^{(3)}_0$ (au)')
    axes[0].axhline(0, color='gray', lw=0.5)
    axes[1].set_ylabel(r'$|\chi^{(3)}_0-\chi^{(3)}_0(10^6)|$ (% of the max)')
    for ax in axes:
        ax.set_xlabel('probe energy (eV)')
        ax.set_xlim(erange)
        ax.set_xticks(np.arange(erange[0], erange[1] + 0.1, 5))
    fig.tight_layout()
    fig.legend(*axes[0].get_legend_handles_labels(), loc='upper center', bbox_to_anchor=(0.5, 0.02), ncol=5,
               frameon=False)
    fig.savefig(os.path.join(out, 'chi_pump.jpeg'))
    plt.close(fig)


HB = Const.Planck_reduced_ev_fs      # hbar (eV fs)
WP_NL = 1.55                          # pump frequency of the yambo_nl runs (eV)


def invint_correction(freqs, chi, dt=0.01):
    # removes the effects of the INVINT integrator of the yambo_nl runs (time step dt, fs) from chi on the grid freqs (eV)
    # of the output polarization: chi(w_eff) exp(-i w_eff dt/2), w_eff = (2 hbar/dt) arctan(w dt/2 hbar) (red shift of
    # the resonances and delay dt/2 of the fields)
    weff = 2 * HB / dt * np.arctan(freqs * dt / (2 * HB))
    chi = np.asarray(chi) * np.exp(-1j * freqs / HB * dt / 2)
    return np.array([np.interp(weff, freqs, y.real) + 1j * np.interp(weff, freqs, y.imag)
                     for y in chi.reshape(-1, chi.shape[-1])]).reshape(chi.shape)


def cranknic_correction(freqs, chi, dt=0.01):
    # removes the static energy map of the CRANKNIC integrator of the yambo_nl runs (time step dt, fs) from chi on the
    # grid freqs (eV): chi(w_eff), w_eff = (2 hbar/dt) arctan(w dt/2 hbar) (no delay of the fields, unlike INVINT)
    weff = 2 * HB / dt * np.arctan(freqs * dt / (2 * HB))
    chi = np.asarray(chi)
    return np.array([np.interp(weff, freqs, y.real) + 1j * np.interp(weff, freqs, y.imag)
                     for y in chi.reshape(-1, chi.shape[-1])]).reshape(chi.shape)


HSEX_CHI = 'nlenrange_7.85-20.250-nlensteps_128-bands_3-6-damp_0.6-cn-nltime_40'   # HSEX chi runs (CRANKNIC)
_tr_hsex = []


def chi_hsex():
    # Transient_Response of the HARTREE+SEX chi (kx16, bands 3-6, damping 0.6 eV in the dynamics, CRANKNIC 0.01 fs,
    # 128 probe frequencies 7.85-20.15 eV, pump 1e6 kW/cm^2), read once
    if not _tr_hsex:
        sine, pp = load_nl('hsex-pulse-E1_1e3-' + HSEX_CHI), load_nl('hsex-Pp-Ep_1e3-EP_1e6-PE_1.55-' + HSEX_CHI)
        _tr_hsex.append(O.Transient_Response.from_Xn(O.Xn_frequency_mixing(pp, X_order=(1, 2), verbose=False),
                                                     O.Xn_single_frequency(sine, X_order=3, verbose=False),
                                                     centrosymmetric=True, verbose=False))
    return _tr_hsex[0]


def fig_eq_hsex(out, erange=(8, 20)):
    """
    Equilibrium eps with HARTREE+SEX: RT probe-only run (smoothing 0.6 eV) vs 1 + 4 pi chi1 of the HSEX chi runs (damping
    0.6 eV, CRANKNIC correction), with the IP RT probe-only run for reference.
    """
    tr = chi_hsex()
    f = tr.freqs
    E, eq = rt_eps(load_nl(rt_id('3-6', None, 10.0, level='hsex')))
    E_ip, eq_ip = rt_eps(load_nl(rt_id('3-6', None, 10.0)))
    eps = {lab: 1 + 4 * np.pi * c for lab, c in [('no correction', tr.chi1),
                                                ('CRANKNIC correction', cranknic_correction(f, tr.chi1))]}
    ref = np.interp(f, E, eq.real) + 1j * np.interp(f, E, eq.imag) - 1
    for lab, y in eps.items():
        print('eq_hsex: 1 + 4 pi chi1 (%s) vs eps_eq of the RT probe-only run: in 10-20 eV rel. rms %.4f, ratio of the '
              'maxima %.3f; best shift of the edge (11-13.5 eV) %+.3f eV' % ((lab,) + compare(f, y - 1, ref, (10, 20))[:2]
                                                                            + compare(f, y - 1, ref, (11, 13.5))[2:]))
    s, s_ip = (E >= erange[0]) & (E <= erange[1]), (E_ip >= erange[0]) & (E_ip <= erange[1])
    fig, axes = plt.subplots(1, 2, figsize=(7, 2.6))
    for ax, part, lab in zip(axes, [np.imag, np.real], ['Im', 'Re']):
        ax.plot(E_ip[s_ip], part(eq_ip[s_ip]), color='gray', ls=':', lw=1, label='RT, IP')
        ax.plot(E[s], part(eq[s]), 'k', lw=1.2, label='RT, HSEX')
        ax.plot(f, part(eps['CRANKNIC correction']), 'r--', lw=1, label=r'$1+4\pi\chi^{(1)}$, HSEX')
        ax.set_xlabel('energy (eV)')
        ax.set_ylabel(r'%s $\varepsilon$' % lab)
        ax.set_xlim(*erange)
    axes[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'eq_hsex.jpeg'))
    plt.close(fig)


def dchi_local(tr, taus, E0P, sigma, wP):
    """
    As tr.dchi(taus, E0P, pump='accumulated', sigma=sigma, wP=wP) of the MPPI Transient_Response tr, with the local weight
    f^2(tau+t) (the pump acts at the emission time; not in MPPI, the notes quote its result): each third order term is
    weighted in time on the probe grid of the runs and read on the output grid tr.w.
    """
    f, n = tr.freqs / Const.HaToeV, int(round(2 * tr.wP / tr.step))
    chi0 = (tr.keys[(1, 0)] - tr.chi1) / (tr.E0P**2 / 4)
    out = []
    for tau in taus:
        W = lambda t: np.exp(-((tau * Const.FsToAu + np.maximum(t, 0)) / (sigma * Const.FsToAu))**2)
        ph = np.exp(-2j * wP / HB * tau)
        out.append(E0P**2 / 4 * (time_weighted(f, chi0, W)[n:-n] - time_weighted(f, tr.keys[(1, 2)], W)[:-2 * n] * ph
                                 - time_weighted(f, tr.keys[(1, -2)], W)[2 * n:] / ph))
    return tr.w, np.array(out)


def compare(w, y, ref, erange):
    # rel. rms distance, ratio of the maxima of the moduli and best rigid shift (eV, > 0 to the right) of the complex
    # spectrum y(w) w.r.t. ref(w) in erange
    s = (w >= erange[0]) & (w <= erange[1])
    shifted = lambda a: np.interp(w[s] - a, w, y.real) + 1j * np.interp(w[s] - a, w, y.imag)
    shift = minimize_scalar(lambda a: np.sum(np.abs(shifted(a) - ref[s])**2), bounds=(-1, 1), method='bounded').x
    return (np.linalg.norm(y[s] - ref[s]) / np.linalg.norm(ref[s]), np.abs(y[s]).max() / np.abs(ref[s]).max(), shift)


def rt_even(pump_int, level='ip'):
    """
    Even part in the pump field of the RT Delta eps (bands 3-6) at the probe times T_PROBES, from the runs with the pump
    along +x and -x. Returns the energies (eV), eps_eq of the probe-only run, the list of the even parts, the delays
    (fs), t_P and sigma (fs), E0P (au) and the pump frequency (eV).
    """
    E, eps_eq = rt_eps(load_nl(rt_id('3-6', None, 10.0, level=level)))
    ref_p, ref_m = load_nl(rt_id('3-6', pump_int, level=level)), load_nl(rt_id('3-6', pump_int, neg=True, level=level))
    t_P, sigma = pump_center(ref_p)
    pump = ref_p.Efield2[0]
    even, taus = [], []
    for tp in T_PROBES:
        run_p, run_m = load_nl(rt_id('3-6', pump_int, tp, level=level)), load_nl(rt_id('3-6', pump_int, tp, True, level))
        even.append((rt_eps(run_p, ref_p)[1] + rt_eps(run_m, ref_m)[1]) / 2 - eps_eq)
        taus.append(run_p.Efield[0]['initial_time'] / Const.FsToAu - t_P)
    return E, eps_eq, even, taus, t_P, sigma, pump['amplitude'], pump['freq_range'][0] * Const.HaToeV


def fig_dchi_rt(out, pump_int=6e7, erange=(13.6, 21.3), name='dchi_rt', level='ip'):
    """
    Transient spectrum Delta eps at the five probe times: RT response (even part in the pump field, bands 3-6) vs
    4 pi dchi from the kx16 chi (total broadening 0.6 eV, INVINT correction) with the envelope at the probe time and
    with the accumulated envelope. For a pump different from 6e7 kW/cm^2 (beyond the third order) only the accumulated
    envelope is plotted, and the RT response is also compared (printed numbers) with the third order RT response, the
    even part at 6e7 rescaled by pump_int/6e7. With level='hsex' the RT runs and the chi are at the HARTREE+SEX level
    (chi of chi_hsex, damping 0.6 eV in the dynamics, CRANKNIC correction; no local weight).
    """
    E, eps_eq, even, taus, t_P, sigma, E0P, wP = rt_even(pump_int, level)
    if pump_int != 6e7:
        third = [y * pump_int / 6e7 for y in rt_even(6e7, level)[2]]
        lab_rt = 'RT (even part)'
    else:
        third, lab_rt = None, 'RT, third order (even part)'

    if level == 'hsex':
        tr, correction, model_names = chi_hsex(), cranknic_correction, ['probe time', 'accumulated']
    else:
        # chi of the kx16 grid, total broadening 0.6 eV
        r = chi_kgrids(out)
        tr = O.Transient_Response(r['freqs'], r['kx16_chi1'], {key: r['kx16_%d_%d' % key] for key in MIX_KEYS},
                                  r['kx16_EP'], WP_NL, centrosymmetric=True, verbose=False)
        correction, model_names = invint_correction, ['probe time', 'local', 'accumulated']
    f = tr.freqs
    models = {'probe time': ('r', '--', 'envelope at the probe time'), 'accumulated': ('b', '-', 'accumulated envelope')}
    deps = {}
    for m in model_names:
        # the phases use the RT pump frequency wP, the shifted grids the one of the yambo_nl runs
        w, d = dchi_local(tr, taus, E0P, sigma, wP) if m == 'local' else tr.dchi(taus, E0P, pump=m, sigma=sigma, wP=wP)
        deps[m] = 4 * np.pi * correction(w, d)
    rt = [np.interp(w, E, y.real) + 1j * np.interp(w, E, y.imag) for y in even]
    print('%s: pump %.1e kW/cm^2, E0P %.4e au (%.4f times the yambo_nl pump, sqrt(I ratio) %.4f), wP %.3f eV, '
          't_P %.3f fs, sigma %.3f fs' % (name, pump_int, E0P, E0P / tr.E0P, np.sqrt(pump_int / 1e6), wP, t_P, sigma))
    if level != 'hsex':
        eq = np.interp(f, E, eps_eq.real) + 1j * np.interp(f, E, eps_eq.imag) - 1
        for lab, c in [('no correction', tr.chi1), ('INVINT correction', invint_correction(f, tr.chi1))]:
            print('dchi_rt: 4 pi chi1 (%s) vs eps_eq - 1 of the RT probe-only run: in 10-25 eV rel. rms %.3f, ratio of '
                  'the maxima %.3f; best shift of the edge (12-17.3 eV) %+.3f eV' % ((lab,) + compare(f, 4 * np.pi * c,
                  eq, (10, 25))[:2] + compare(f, 4 * np.pi * c, eq, (12, 17.3))[2:]))
    print('%s: chi vs RT even part in %.1f-%.1f eV: rel. rms / ratio of the maxima / best shift (eV)' % ((name,) + erange))
    for k, tau in enumerate(taus):
        print('  tau %+.2f fs | ' % tau + ' | '.join('%s %.2f / %.2f / %+.2f' % ((m,) + compare(w, deps[m][k], rt[k],
                                                                                              erange)) for m in deps))
    if third is not None:
        rt3 = [np.interp(w, E, y.real) + 1j * np.interp(w, E, y.imag) for y in third]
        print('%s: RT third order (6e7 x %d) vs RT %.0e and accumulated chi vs RT third order: rel. rms / ratio of '
              'the maxima / best shift (eV)' % (name, round(pump_int / 6e7), pump_int))
        for k, tau in enumerate(taus):
            print('  tau %+.2f fs | third order vs RT %.2f / %.2f / %+.2f | chi vs third order %.2f / %.2f / %+.2f'
                  % ((tau,) + compare(w, rt3[k], rt[k], erange) + compare(w, deps['accumulated'][k], rt3[k], erange)))

    if third is not None:
        models = {'accumulated': models['accumulated']}
    s = (E >= w[0]) & (E <= w[-1])
    fig, axes = plt.subplots(len(taus), 2, figsize=(7, 1.75 * len(taus)), sharex=True)
    for k, tau in enumerate(taus):
        for ax, part, lab in zip(axes[k], [np.real, np.imag], ['Re', 'Im']):
            ax.plot(E[s], 1e3 * part(even[k][s]), 'k', lw=1.2, label=lab_rt)
            for m, (c, ls, ml) in models.items():
                ax.plot(w, 1e3 * part(deps[m][k]), color=c, ls=ls, lw=1, label=r'from the $\chi$, ' + ml)
            ax.axhline(0, color='gray', lw=0.5)
            ax.set_ylabel(r'%s $\Delta\varepsilon$ ($10^{-3}$)' % lab)
            ax.text(0.02, 0.92, r'$\tau$ = %+.2f fs' % tau, transform=ax.transAxes, va='top', fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel('probe energy (eV)')
        ax.set_xlim(*((11, 17) if level == 'hsex' else (13, 22)))
    fig.tight_layout()
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc='upper center', bbox_to_anchor=(0.5, 0.005),
               ncol=len(models) + 1, frameon=False)
    fig.savefig(os.path.join(out, name + '.jpeg'))
    plt.close(fig)


DAVIDE_DIR = os.path.join(NL_CHI, 'Davide_data', '20250312_Report_simulations_data_kx16_ball')
DAVIDE_FILE = 'o-abs_0-50eV_sm0.6eV.YPP-eps_along_E'


def ypp_fields(fname):
    # (initial time (fs), width (fs), frequency (eV), intensity (kW/cm^2)) of the fields (probe, pump) in the header of a
    # ypp file
    hdr = [l for l in open(fname) if l.startswith('#')]
    get = lambda key: [float(l.split(':')[1].split()[0]) for l in hdr if '[FIELDs] ' + key in l]
    return list(zip(get('Initial time'), get('Width'), get('Frequency range'), get('Intensity')))


def load_davide(prefix='td-hsex-cv-ber'):
    """
    RT runs of D. Sangalli of the set prefix (ypp spectra, smoothing 0.6 eV), sorted by the delay tau = t_p - t_P:
    energies (eV), eps along x (n_delays, n_energies) and delays (fs). The first one (-10.66 fs) is the reference.
    """
    folders = glob.glob(os.path.join(DAVIDE_DIR, prefix + '*delta*fs*'))
    t_p = np.array([ypp_fields(os.path.join(f, DAVIDE_FILE))[0][0] for f in folders])
    t_start, sigma, wP, _ = ypp_fields(os.path.join(folders[0], DAVIDE_FILE))[1]
    t_P = t_start + np.pi * HB / wP * np.rint(wP / (np.pi * HB) * 5 * sigma)   # maximum of the QSFIELD 5SIGMA pump
    # one file of the HSEX set starts with a stray 'ok' line
    data = [np.loadtxt(os.path.join(folders[i], DAVIDE_FILE), comments=['#', 'ok']) for i in np.argsort(t_p)]
    return data[0][:, 0], np.array([d[:, 4] + 1j * d[:, 1] for d in data]), np.sort(t_p) - t_P


def highpass(y, t, s=0.55):
    # y (delays, energies) minus its Gaussian smoothing in the delay t (width s, fs): removes the part that follows the
    # pump envelope and keeps the oscillations at 2wP (the smoothing keeps exp(-(2 wP s)^2/2) = 3% of them)
    K = np.exp(-(t[:, None] - t[None, :])**2 / (2 * s**2))
    return y - (K @ y) / K.sum(1)[:, None]


_hsex_ac = {}


def hsex_ac_data():
    """
    Data of the figures maps_hsex, ac_hsex and ac_traces_hsex (HSEX, pump 6e9 kW/cm^2), computed once: Im Delta eps of
    the RT runs of D. Sangalli (bands 3-12, pump along +x only) and of the HSEX chi (bands 3-6, accumulated envelope,
    pump of the RT runs) at the same delays; AC part of the chi (chi_+- terms only) with the two envelope models; RT
    after the high-pass in the delay; complex amplitude a(w) of the 2wP component, y = f^2(tau) Re[a e^{2i wP tau}].
    """
    if _hsex_ac:
        return _hsex_ac
    E, eps, taus_all = load_davide()
    taus = taus_all[1:]
    tr = chi_hsex()
    w = tr.w
    pump = load_nl(rt_id('3-6', 6e9, level='hsex'))                   # pump of the RT runs (our pump-only run)
    E0P, sigma, wP = pump.Efield2[0]['amplitude'], pump_center(pump)[1], pump.Efield2[0]['freq_range'][0] * Const.HaToeV
    tr_ac = O.Transient_Response(tr.freqs, tr.chi1, {(1, 0): tr.chi1, (1, 2): tr.keys[(1, 2)], (1, -2): tr.keys[(1, -2)]},
                                 tr.E0P, tr.wP, verbose=False)        # chi_0 = 0: only the chi_+- terms

    def im_deps(t, trr, model='accumulated'):
        _, d = trr.dchi(t, E0P, pump=model, sigma=sigma, wP=wP)
        return (4 * np.pi * cranknic_correction(w, d)).imag

    # maps: change w.r.t. the reference delay, the chi on a 0.02 eV grid (cubic spline), 11-17 eV
    Em = np.arange(w[0], w[-1] + 1e-9, 0.02)
    Em = Em[(Em >= 11) & (Em <= 17)]
    full = CubicSpline(w, im_deps(taus_all, tr), axis=1)(Em)
    d = dict(Em=Em, taus=taus, w=w, sigma=sigma, wP=wP, dim_chi=full[1:] - full[0],
             dim_rt=np.array([np.interp(Em, E, y - eps[0].imag) for y in eps[1:].imag]))
    # oscillations at 2wP on the output grid w
    d['hp_rt'] = highpass(np.array([np.interp(w, E, y - eps[0].imag) for y in eps[1:].imag]), taus)
    d['ac'] = {m: im_deps(taus, tr_ac, m) for m in ['probe time', 'accumulated']}
    d['taus_fine'] = np.arange(taus[0], taus[-1] + 1e-9, 0.02)
    d['ac_fine'] = {m: im_deps(d['taus_fine'], tr_ac, m) for m in ['probe time', 'accumulated']}
    W2 = 2 * wP / HB                                                    # 2 wP (rad/fs)

    def demod(y, sel=None):
        # complex amplitude a(w) of the 2wP component, fit of f^2 Re[a e^{2i wP tau}] on the delays sel (|tau| < 6 fs)
        sel = np.abs(taus) < 6 if sel is None else sel
        g = np.exp(-taus[sel]**2 / sigma**2)
        return 2 * (y[sel] * (g * np.exp(-1j * W2 * taus[sel]))[:, None]).sum(0) / (g**2).sum()

    ac, hp_ac = d['ac']['accumulated'], highpass(d['ac']['accumulated'], taus)
    d['a_rt'], d['a_chi'] = demod(d['hp_rt']), demod(ac)
    a_hp, a_full = demod(hp_ac), demod(highpass(im_deps(taus, tr), taus))
    big = np.abs(d['a_chi']) > 0.5 * np.abs(d['a_chi']).max()

    # printed numbers
    s, m3 = (w >= 11) & (w <= 17), np.abs(taus) < 3
    wins = {'all': np.ones(len(taus), bool), 'tau < -1.5': taus < -1.5, '|tau| < 1.5': np.abs(taus) < 1.5,
            'tau > 1.5': taus > 1.5}
    rel = lambda x, y, m: np.linalg.norm((x - y)[m][:, s]) / np.linalg.norm(y[m][:, s])
    print('hsex_ac: %d delays %.2f..%.2f fs, E0P %.4e au, sigma %.3f fs, wP %.3f eV' % (len(taus), taus[0], taus[-1], E0P,
                                                                                      sigma, wP))
    print('maps, Delta Im eps in 11-17 eV, RT vs chi: rel. rms / ratio of the maxima RT/chi: ' + ' | '.join(
          '%s %.2f / %.2f' % (k, np.linalg.norm(d['dim_rt'][m] - d['dim_chi'][m]) / np.linalg.norm(d['dim_chi'][m]),
                              np.abs(d['dim_rt'][m]).max() / np.abs(d['dim_chi'][m]).max()) for k, m in wins.items()))
    print('filter on the chi map (|tau| < 3 fs): rel. rms w.r.t. the AC part %.3f (high-pass of the AC part), %.3f '
          '(high-pass of the full map); 2wP amplitude kept %.3f, phase changed by <= %.3f rad' % (
          rel(hp_ac, ac, m3), rel(highpass(im_deps(taus, tr), taus), ac, m3), np.median(np.abs(a_full / d['a_chi'])[big]),
          np.abs(np.angle(a_full / d['a_chi'])[big]).max()))
    print('RT high-pass vs chi AC part (|tau| < 3 fs): rel. rms %.2f, ratio of the maxima %.2f' % (
          rel(d['hp_rt'], ac, m3), np.abs(d['hp_rt'][m3][:, s]).max() / np.abs(ac[m3][:, s]).max()))
    print('2wP component:  w (eV) | |a| RT | |a| chi | RT/chi (same filter) | arg a RT - chi (rad)')
    for wi in [12.3, 12.9, 13.5, 13.95, 14.6, 15.0, 15.5, 16.0, 16.5]:
        i = np.argmin(abs(w - wi))
        print('                %6.2f | %.4f | %.4f | %.2f | %+.2f' % (w[i], abs(d['a_rt'][i]), abs(d['a_chi'][i]),
              abs(d['a_rt'][i] / a_hp[i]), np.angle(d['a_rt'][i] / d['a_chi'][i])))
    cent = lambda y: (taus[:, None] * y[:, s]**2).sum() / (y[:, s]**2).sum()
    width = lambda y: np.sqrt((taus[:, None]**2 * y[:, s]**2).sum() / (y[:, s]**2).sum() - cent(y)**2)
    print('envelope of the fringes vs the delay (11-17 eV): centroid / rms width of the AC power; RT high-pass %+.2f / '
          '%.2f fs' % (cent(d['hp_rt']), width(d['hp_rt'])))
    for m, y in d['ac'].items():
        y = highpass(y, taus)
        print('  chi %-11s %+.2f / %.2f fs | RT vs chi (same filter), rel. rms / ratio of the norms: ' % (m, cent(y),
              width(y)) + ' | '.join('%s %.2f / %.2f' % (k, rel(d['hp_rt'], y, mm), np.linalg.norm(d['hp_rt'][mm][:, s])
              / np.linalg.norm(y[mm][:, s])) for k, mm in wins.items() if k != 'all'))
    _hsex_ac.update(d)
    return _hsex_ac


def fig_maps_hsex(out, tau_cut=0.):
    """
    Delta Im eps vs the delay and the probe energy (HSEX, pump 6e9 kW/cm^2): RT runs of D. Sangalli (bands 3-12) and the
    chi (bands 3-6, accumulated envelope), same color scale, and the spectra at the delay nearest to tau_cut.
    """
    d = hsex_ac_data()
    i = np.argmin(np.abs(d['taus'] - tau_cut))
    v = max(np.abs(d['dim_rt']).max(), np.abs(d['dim_chi']).max())
    fig = plt.figure(figsize=(7, 2.9))
    gs = fig.add_gridspec(1, 5, width_ratios=[2, 2, 0.1, 0.45, 1.2], wspace=0.08)
    axes = [fig.add_subplot(gs[0])]
    axes += [fig.add_subplot(gs[k], sharey=axes[0]) for k in [1, 4]]
    for ax, y, t in zip(axes[:2], [d['dim_rt'], d['dim_chi']], ['RT', r'from the $\chi$']):
        pc = ax.pcolormesh(d['taus'], d['Em'], y.T, cmap='bwr', norm=TwoSlopeNorm(0, -v, v), shading='gouraud',
                           rasterized=True)
        ax.axvline(d['taus'][i], color='k', lw=0.8, ls='--')
        ax.set_xlabel(r'delay $\tau$ (fs)')
        ax.set_title(t, fontsize=9)
    fig.colorbar(pc, cax=fig.add_subplot(gs[2])).ax.set_title(r'$\Delta$Im $\varepsilon$', fontsize=8)
    axes[0].set_ylabel('probe energy (eV)')
    axes[2].plot(d['dim_rt'][i], d['Em'], 'k', lw=1.2, label='RT')
    axes[2].plot(d['dim_chi'][i], d['Em'], 'b--', lw=1, label=r'$\chi$')
    axes[2].axvline(0, color='gray', lw=0.5)
    axes[2].set_xlabel(r'$\Delta$Im $\varepsilon$, $\tau$ = %+.2f fs' % d['taus'][i])
    axes[2].legend(frameon=False)
    for ax in axes[1:]:
        ax.tick_params(labelleft=False)
    fig.savefig(os.path.join(out, 'maps_hsex.jpeg'))
    plt.close(fig)


def fig_ac_hsex(out):
    """
    Oscillations at 2wP (HSEX, pump 6e9 kW/cm^2): RT of D. Sangalli after the high-pass in the delay and AC part of the
    chi (chi_+- terms, accumulated envelope), same color scale; amplitude |a| and phase arg a of the 2wP component.
    """
    d = hsex_ac_data()
    w, ac = d['w'], d['ac']['accumulated']
    v = np.abs(ac).max()
    fig = plt.figure(figsize=(7, 3.0))
    gs = fig.add_gridspec(1, 6, width_ratios=[1, 1, 0.07, 0.3, 0.8, 0.8], wspace=0.1)
    axes = [fig.add_subplot(gs[0])]
    axes += [fig.add_subplot(gs[k], sharey=axes[0]) for k in [1, 4, 5]]
    for ax, y, t in zip(axes[:2], [d['hp_rt'], ac], ['RT, high-pass', r'$\chi$, AC part']):
        pc = ax.pcolormesh(d['taus'], w, y.T, cmap='bwr', norm=TwoSlopeNorm(0, -v, v), shading='gouraud',
                           rasterized=True)
        ax.set_xlim(-6, 6)
        ax.set_xlabel(r'delay $\tau$ (fs)')
        ax.set_title(t, fontsize=9)
    fig.colorbar(pc, cax=fig.add_subplot(gs[2])).ax.set_title(r'$\Delta$Im $\varepsilon$', fontsize=8)
    axes[0].set_ylabel('probe energy (eV)')
    axes[0].set_ylim(w[0], w[-1])
    axes[2].plot(np.abs(d['a_rt']), w, 'k', lw=1.2, label='RT')
    axes[2].plot(np.abs(d['a_chi']), w, 'b--', lw=1, label=r'$\chi$')
    axes[2].set_xlabel(r'amplitude $|a|$')
    axes[2].legend(frameon=False)
    axes[3].plot(np.unwrap(np.angle(d['a_rt'])), w, 'k', lw=1.2)
    axes[3].plot(np.unwrap(np.angle(d['a_chi'])), w, 'b--', lw=1)
    axes[3].set_xlabel(r'phase arg $a$ (rad)')
    for ax in axes[1:]:
        ax.tick_params(labelleft=False)
    fig.savefig(os.path.join(out, 'ac_hsex.jpeg'))
    plt.close(fig)


def fig_ac_traces_hsex(out, energies=(12.3, 16.0)):
    """
    Oscillations at 2wP vs the delay at two probe energies (HSEX, pump 6e9 kW/cm^2): RT of D. Sangalli after the
    high-pass and AC part of the chi with the accumulated envelope and with the envelope at the probe time.
    """
    d = hsex_ac_data()
    fig, axes = plt.subplots(1, len(energies), figsize=(7, 2.5), sharey=False)
    for ax, wi in zip(axes, energies):
        i = np.argmin(np.abs(d['w'] - wi))
        ax.plot(d['taus_fine'], d['ac_fine']['accumulated'][:, i], 'b', lw=1, label=r'$\chi$, accumulated envelope')
        ax.plot(d['taus_fine'], d['ac_fine']['probe time'][:, i], 'r--', lw=0.8,
                label=r'$\chi$, envelope at the probe time')
        ax.plot(d['taus'], d['hp_rt'][:, i], 'ko', ms=2.5, label='RT, high-pass')
        ax.axhline(0, color='gray', lw=0.5)
        ax.axvline(0, color='gray', lw=0.5, ls=':')
        ax.set_xlim(-7.5, 7.5)
        ax.set_xlabel(r'delay $\tau$ (fs)')
        ax.set_title('%.2f eV' % d['w'][i], fontsize=9)
    axes[0].set_ylabel(r'$\Delta$Im $\varepsilon$, oscillating part')
    fig.tight_layout()
    fig.legend(*axes[0].get_legend_handles_labels(), loc='upper center', bbox_to_anchor=(0.5, 0.02), ncol=3,
               frameon=False)
    fig.savefig(os.path.join(out, 'ac_traces_hsex.jpeg'))
    plt.close(fig)


FIGURES = {'rt_setup': fig_rt_setup, 'rt_even_odd': fig_rt_even_odd, 'chi_kgrid': fig_chi_kgrid,
           'chi_pump': fig_chi_pump, 'dchi_rt': fig_dchi_rt,
           'dchi_rt_6e9': lambda out: fig_dchi_rt(out, pump_int=6e9, name='dchi_rt_6e9'), 'eq_hsex': fig_eq_hsex,
           'dchi_rt_hsex': lambda out: fig_dchi_rt(out, erange=(11, 17), name='dchi_rt_hsex', level='hsex'),
           'dchi_rt_hsex_6e9': lambda out: fig_dchi_rt(out, pump_int=6e9, erange=(11, 17), name='dchi_rt_hsex_6e9',
                                                       level='hsex'),
           'maps_hsex': fig_maps_hsex, 'ac_hsex': fig_ac_hsex, 'ac_traces_hsex': fig_ac_traces_hsex}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Figures of the technical notes on LiF')
    parser.add_argument('figures', nargs='*', help='names of the figures (default: all): %s' % ', '.join(FIGURES))
    parser.add_argument('--out', default=os.path.join(HERE, 'Figures'), help='output folder of the jpeg files')
    args = parser.parse_args()
    out = os.path.expanduser(args.out)
    os.makedirs(out, exist_ok=True)
    for name in args.figures or FIGURES:
        FIGURES[name](out)
