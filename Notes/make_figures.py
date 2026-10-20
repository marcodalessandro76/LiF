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
from scipy.special import erf
from scipy.optimize import minimize_scalar
from mppi import Parsers as P
from mppi import Optics as O
from mppi.Utilities import Constants as Const

HERE = os.path.dirname(os.path.abspath(__file__))
NL_CHI = os.environ.get('LIF_NL_CHI', os.path.join(HERE, '..', 'NL_Chi'))   # data folder (env variable for tests)
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


def envelope_weight(tau, t, sigma, model):
    # intensity envelope (|E_P|^2 in units of E0P^2/4) that weights the response at the time t (fs) after a probe at the
    # delay tau (fs): 'probe time' f^2(tau), 'local' f^2(tau+t) (pump at the emission time), 'accumulated'
    # W(tau,t) = (1/t) int_0^t f^2(tau+t') dt', with f^2 = exp(-t^2/sigma^2)
    f2 = np.exp(-(tau / sigma)**2)
    if model == 'probe time':
        return f2 * np.ones_like(t)
    if model == 'local':
        return np.exp(-((tau + np.maximum(t, 0)) / sigma)**2)
    tt = np.where(t > 1e-9, t, 1.)
    W = np.sqrt(np.pi) * sigma / (2 * tt) * (erf((tau + tt) / sigma) - erf(tau / sigma))
    return np.where(t > 1e-9, W, f2)


def time_weighted(freqs, chi, W):
    # chi(w) on the uniform grid freqs (eV) -> int dt e^{iwt} W(t) R(t), with R(t) the response function of chi (discrete
    # Fourier transforms on the grid: the period in t is 2 pi hbar/step, 42.7 fs for the step wP/16)
    t = np.fft.fftfreq(len(freqs), (freqs[1] - freqs[0]) / (2 * np.pi * HB))
    E = np.exp(-1j * np.outer(t, freqs) / HB)
    return E.conj().T @ ((E @ chi) * W(t)) / len(freqs)


def dchi_from_chi(freqs, chi0, key_p2, key_m2, taus, E0P, sigma, wP, model):
    """
    Third order dchi(w,tau) (au) for a delta probe at the delay tau from the maximum of the pump E0P f sin(wP(t-t_P)):
    (E0P^2/4) int dt e^{iwt} W(tau,t) [R_0 - R_+ e^{-2i wP tau} - R_- e^{2i wP tau}] (Sec. 1.5 of the notes), with
    chi_+(w) = chi_(1,2)(w-2wP) and chi_-(w) = chi_(1,-2)(w+2wP) read on the probe grid freqs with the pump of the yambo_nl
    runs (2 WP_NL = n steps), the phases with the pump frequency wP (eV) of the RT runs. Returns the output frequencies
    freqs[n:-n] and dchi, shape (len(taus), len(freqs)-2n).
    """
    step = freqs[1] - freqs[0]
    n = int(round(2 * WP_NL / step))
    assert abs(n * step - 2 * WP_NL) < 1e-3 * step, '2 wP is not a multiple of the frequency step'
    out = []
    for tau in taus:
        W = lambda t: envelope_weight(tau, t, sigma, model)
        ph = np.exp(-2j * wP / HB * tau)
        # the time weight of the keys (1,+-2) is applied on their own grid and then shifted: a shift in frequency is a
        # phase in time, which does not change W R
        out.append(E0P**2 / 4 * (time_weighted(freqs, chi0, W)[n:-n] - time_weighted(freqs, key_p2, W)[:-2 * n] * ph
                                 - time_weighted(freqs, key_m2, W)[2 * n:] / ph))
    return freqs[n:-n], np.array(out)


def compare(w, y, ref, erange):
    # rel. rms distance, ratio of the maxima of the moduli and best rigid shift (eV, > 0 to the right) of the complex
    # spectrum y(w) w.r.t. ref(w) in erange
    s = (w >= erange[0]) & (w <= erange[1])
    shifted = lambda a: np.interp(w[s] - a, w, y.real) + 1j * np.interp(w[s] - a, w, y.imag)
    shift = minimize_scalar(lambda a: np.sum(np.abs(shifted(a) - ref[s])**2), bounds=(-1, 1), method='bounded').x
    return (np.linalg.norm(y[s] - ref[s]) / np.linalg.norm(ref[s]), np.abs(y[s]).max() / np.abs(ref[s]).max(), shift)


def rt_even(pump_int):
    """
    Even part in the pump field of the RT Delta eps (bands 3-6) at the probe times T_PROBES, from the runs with the pump
    along +x and -x. Returns the energies (eV), eps_eq of the probe-only run, the list of the even parts, the delays
    (fs), t_P and sigma (fs), E0P (au) and the pump frequency (eV).
    """
    E, eps_eq = rt_eps(load_nl(rt_id('3-6', None, 10.0)))
    ref_p, ref_m = load_nl(rt_id('3-6', pump_int)), load_nl(rt_id('3-6', pump_int, neg=True))
    t_P, sigma = pump_center(ref_p)
    pump = ref_p.Efield2[0]
    even, taus = [], []
    for tp in T_PROBES:
        run_p, run_m = load_nl(rt_id('3-6', pump_int, tp)), load_nl(rt_id('3-6', pump_int, tp, neg=True))
        even.append((rt_eps(run_p, ref_p)[1] + rt_eps(run_m, ref_m)[1]) / 2 - eps_eq)
        taus.append(run_p.Efield[0]['initial_time'] / Const.FsToAu - t_P)
    return E, eps_eq, even, taus, t_P, sigma, pump['amplitude'], pump['freq_range'][0] * Const.HaToeV


def fig_dchi_rt(out, pump_int=6e7, erange=(13.6, 21.3), name='dchi_rt'):
    """
    Transient spectrum Delta eps at the five probe times: RT response (even part in the pump field, bands 3-6) vs
    4 pi dchi from the kx16 chi (total broadening 0.6 eV, INVINT correction) with the envelope at the probe time and
    with the accumulated envelope. For a pump different from 6e7 kW/cm^2 (beyond the third order) only the accumulated
    envelope is plotted, and the RT response is also compared (printed numbers) with the third order RT response, the
    even part at 6e7 rescaled by pump_int/6e7.
    """
    E, eps_eq, even, taus, t_P, sigma, E0P, wP = rt_even(pump_int)
    if pump_int != 6e7:
        third = [y * pump_int / 6e7 for y in rt_even(6e7)[2]]
        lab_rt = 'RT (even part)'
    else:
        third, lab_rt = None, 'RT, third order (even part)'

    # chi of the kx16 grid, total broadening 0.6 eV
    r = chi_kgrids(out)
    f = r['freqs']
    chi0 = (r['kx16_1_0'] - r['kx16_chi1']) / (r['kx16_EP']**2 / 4)
    models = {'probe time': ('r', '--', 'envelope at the probe time'), 'accumulated': ('b', '-', 'accumulated envelope')}
    deps = {}
    for m in ['probe time', 'local', 'accumulated']:
        w, d = dchi_from_chi(f, chi0, r['kx16_1_2'], r['kx16_1_-2'], taus, E0P, sigma, wP, m)
        deps[m] = 4 * np.pi * invint_correction(w, d)
    rt = [np.interp(w, E, y.real) + 1j * np.interp(w, E, y.imag) for y in even]
    print('dchi_rt: pump %.1e kW/cm^2, E0P %.4e au (%.4f times the yambo_nl pump, sqrt(I ratio) %.4f), wP %.3f eV, '
          't_P %.3f fs, sigma %.3f fs' % (pump_int, E0P, E0P / r['kx16_EP'], np.sqrt(pump_int / 1e6), wP, t_P, sigma))
    eq = np.interp(f, E, eps_eq.real) + 1j * np.interp(f, E, eps_eq.imag) - 1
    for lab, c in [('no correction', r['kx16_chi1']), ('INVINT correction', invint_correction(f, r['kx16_chi1']))]:
        print('dchi_rt: 4 pi chi1 (%s) vs eps_eq - 1 of the RT probe-only run: in 10-25 eV rel. rms %.3f, ratio of the '
              'maxima %.3f; best shift of the edge (12-17.3 eV) %+.3f eV' % ((lab,) + compare(f, 4 * np.pi * c, eq,
              (10, 25))[:2] + compare(f, 4 * np.pi * c, eq, (12, 17.3))[2:]))
    print('dchi_rt: chi vs RT even part in %.1f-%.1f eV: rel. rms / ratio of the maxima / best shift (eV)' % erange)
    for k, tau in enumerate(taus):
        print('  tau %+.2f fs | ' % tau + ' | '.join('%s %.2f / %.2f / %+.2f' % ((m,) + compare(w, deps[m][k], rt[k],
                                                                                              erange)) for m in deps))
    if third is not None:
        rt3 = [np.interp(w, E, y.real) + 1j * np.interp(w, E, y.imag) for y in third]
        print('dchi_rt: RT third order (6e7 x %d) vs RT %.0e and accumulated chi vs RT third order: rel. rms / ratio of '
              'the maxima / best shift (eV)' % (round(pump_int / 6e7), pump_int))
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
        ax.set_xlim(13, 22)
    fig.tight_layout()
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc='upper center', bbox_to_anchor=(0.5, 0.005),
               ncol=len(models) + 1, frameon=False)
    fig.savefig(os.path.join(out, name + '.jpeg'))
    plt.close(fig)


FIGURES = {'rt_setup': fig_rt_setup, 'rt_even_odd': fig_rt_even_odd, 'chi_kgrid': fig_chi_kgrid,
           'chi_pump': fig_chi_pump, 'dchi_rt': fig_dchi_rt,
           'dchi_rt_6e9': lambda out: fig_dchi_rt(out, pump_int=6e9, name='dchi_rt_6e9')}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Figures of the technical notes on LiF')
    parser.add_argument('figures', nargs='*', help='names of the figures (default: all): %s' % ', '.join(FIGURES))
    parser.add_argument('--out', default=os.path.join(HERE, 'Figures'), help='output folder of the jpeg files')
    args = parser.parse_args()
    out = os.path.expanduser(args.out)
    os.makedirs(out, exist_ok=True)
    for name in args.figures or FIGURES:
        FIGURES[name](out)
