"""
Figures of the article (PumpProbe-NL-Chi/main.tex), made from the runs in NL_Chi, which are on the cluster only.

The data and the calculations are NOT duplicated here: they come from the figure script of the technical notes
(Notes/make_figures.py, loaded as the module `notes`), so that the article and the notes show the same numbers. This
script only sets the style of the article (column width, fonts) and the layout of the panels. A new calculation goes in
Notes/make_figures.py (as a data function, e.g. `hsex_ac_data`, `ip_maps_data`) and is plotted here.

Usage, on ismhpc from ~/work/LiF/PumpProbe-NL-Chi:

    LIF_NL_CHI=$HOME/work/LiF/NL_Chi ~/miniconda3/bin/python make_figures.py [figure ...] [--out DIR]

Without figure names all the figures are made. As for the notes, write the jpeg files outside the cluster copy of the
repository (e.g. --out ~/tmp_claude/paper_figures) and copy them to PumpProbe-NL-Chi/Figures on the laptop, where they
are committed.
"""
import os
import argparse
import importlib.util
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))


def load_notes():
    # the figure script of the notes as the module `notes` (not with `import make_figures`: this file has the same name)
    spec = importlib.util.spec_from_file_location('notes', os.path.join(HERE, '..', 'Notes', 'make_figures.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


notes = load_notes()

# Style of the article (Phys. Rev. B, revtex4-2 reprint: column width 3.375 in, page width 7.0 in), set after loading
# the notes, which set their own rcParams
COLUMN, PAGE = 3.375, 7.0
plt.rcParams.update({'font.size': 8, 'axes.labelsize': 8, 'legend.fontsize': 7, 'xtick.labelsize': 7,
                     'ytick.labelsize': 7, 'savefig.dpi': 300, 'savefig.bbox': 'tight'})

MODELS = {'probe time': ('r', '--', 'envelope at the probe time'), 'accumulated': ('b', '-', 'accumulated envelope')}


def fig_eq_hsex(out, erange=(8, 20)):
    """
    Fig. eq_hsex: Im eps at equilibrium with HSEX (kx16, bands 3-6): RT probe-only run (smoothing 0.6 eV), 1 + 4 pi chi1
    of the chi runs (damping 0.6 eV, CRANKNIC correction) and the IP RT probe-only run for reference (as Fig. 9 of the
    notes, imaginary part only).
    """
    tr = notes.chi_hsex()
    f = tr.freqs
    E, eq = notes.rt_eps(notes.load_nl(notes.rt_id('3-6', None, 10.0, level='hsex')))
    E_ip, eq_ip = notes.rt_eps(notes.load_nl(notes.rt_id('3-6', None, 10.0)))
    eps = 1 + 4 * np.pi * notes.cranknic_correction(f, tr.chi1)
    s, s_ip = (E >= erange[0]) & (E <= erange[1]), (E_ip >= erange[0]) & (E_ip <= erange[1])
    fig, ax = plt.subplots(figsize=(COLUMN, 2.3))
    ax.plot(E_ip[s_ip], eq_ip[s_ip].imag, color='gray', ls=':', lw=1, label='RT, independent particles')
    ax.plot(E[s], eq[s].imag, 'k', lw=1.2, label='RT, HSEX')
    ax.plot(f, eps.imag, 'r--', lw=1, label=r'$1+4\pi\chi^{(1)}$, HSEX')
    ax.set_xlabel('energy (eV)')
    ax.set_ylabel(r'Im $\varepsilon$')
    ax.set_xlim(*erange)
    ax.set_ylim(0, 1.45 * eq[s].imag.max())                           # room for the legend above the curves
    ax.legend(frameon=False, loc='upper right')
    fig.savefig(os.path.join(out, 'eq_hsex.jpeg'))
    plt.close(fig)


def fig_chi3_hsex(out):
    """
    Fig. chi3_hsex: Im chi_0, chi_+, chi_- with HSEX (CRANKNIC correction) on the output grid of the transient response,
    vs the frequency of the spectrum, with the exciton E_x (maximum of Im chi1) and E_x + 2wP (as Fig. 10 of the notes).
    """
    tr = notes.chi_hsex()
    w, comp = tr.w, tr.components()
    chi1 = notes.cranknic_correction(w, np.interp(w, tr.freqs, tr.chi1.real) + 1j * np.interp(w, tr.freqs, tr.chi1.imag))
    Ex = w[np.argmax(chi1.imag)]
    labels = {'chi_0': r'$\chi^{(3)}_0$', 'chi_+': r'$\chi^{(3)}_+$', 'chi_-': r'$\chi^{(3)}_-$'}
    lines = {'chi_0': [(Ex, '--')], 'chi_+': [(Ex, ':'), (Ex + 2 * tr.wP, '--')], 'chi_-': [(Ex, ':')]}
    fig, axes = plt.subplots(1, 3, figsize=(PAGE, 2.0), sharex=True)
    for ax, key in zip(axes, labels):
        for x, ls in lines[key]:
            ax.axvline(x, color='k', lw=0.6, ls=ls)
        ax.plot(w, notes.cranknic_correction(w, comp[key]).imag / 1e2, 'r', lw=1)
        ax.axhline(0, color='gray', lw=0.5)
        ax.text(0.04, 0.93, labels[key], transform=ax.transAxes, va='top', fontsize=9)
        ax.set_xlabel('energy (eV)')
        ax.set_xlim(w[0], w[-1])
    axes[0].set_ylabel(r'Im $\chi$ ($10^2$ au)')
    fig.tight_layout(w_pad=0.8)
    fig.savefig(os.path.join(out, 'chi3_hsex.jpeg'))
    plt.close(fig)
    print('chi3_hsex: exciton %.2f eV, E_x + 2wP %.2f eV' % (Ex, Ex + 2 * tr.wP))


def fig_dchi_cuts(out, t_probes=(13.0, 16.0, 19.0), erange=(11, 17)):
    """
    Fig. dchi_cuts: Im Delta eps at three delays (before, at and after the pump maximum), HSEX, pump 6e7 kW/cm^2: RT even
    part vs 4 pi dchi from the chi with the envelope at the probe time and with the accumulated envelope (extract of
    Fig. 11 of the notes, imaginary part only).
    """
    E, eps_eq, even, taus, t_P, sigma, E0P, wP = notes.rt_even(6e7, 'hsex', list(t_probes))
    tr = notes.chi_hsex()
    deps = {}
    for m in MODELS:
        w, d = tr.dchi(taus, E0P, pump=m, sigma=sigma, wP=wP)
        deps[m] = 4 * np.pi * notes.cranknic_correction(w, d)
    s = (E >= erange[0]) & (E <= erange[1])
    fig, axes = plt.subplots(len(taus), 1, figsize=(COLUMN, 1.45 * len(taus) + 0.6), sharex=True)
    for ax, k in zip(axes, range(len(taus))):
        ax.plot(E[s], 1e3 * even[k][s].imag, 'k', lw=1.2, label='RT, third order (even part)')
        for m, (c, ls, lab) in MODELS.items():
            ax.plot(w, 1e3 * deps[m][k].imag, color=c, ls=ls, lw=1, label=r'from the $\chi$, ' + lab)
        ax.axhline(0, color='gray', lw=0.5)
        ax.set_ylabel(r'$\Delta$Im $\varepsilon$ ($10^{-3}$)')
        ax.text(0.03, 0.92, r'$\tau$ = %+.2f fs' % taus[k], transform=ax.transAxes, va='top', fontsize=8)
        rt = np.interp(w, E, even[k].real) + 1j * np.interp(w, E, even[k].imag)
        print('dchi_cuts: tau %+.2f fs | ' % taus[k] + ' | '.join('%s %.2f / %.2f / %+.2f' % ((m,) + notes.compare(
              w, deps[m][k], rt, erange)) for m in MODELS))
    axes[-1].set_xlabel('probe energy (eV)')
    axes[-1].set_xlim(*erange)
    fig.tight_layout(h_pad=0.3, rect=(0, 0.11, 1, 1))
    fig.legend(*axes[0].get_legend_handles_labels(), loc='lower center', bbox_to_anchor=(0.55, 0.0), ncol=1,
               frameon=False)
    fig.savefig(os.path.join(out, 'dchi_cuts.jpeg'))
    plt.close(fig)


def fig_maps(out, tau_cut=0.):
    """
    Fig. maps: Delta Im eps vs the delay and the probe energy (HSEX, bands 3-6, pump 6e7 kW/cm^2, 41 delays): RT even part
    and the chi with the accumulated envelope, same color scale, and the spectra at the delay nearest to tau_cut (as
    Fig. 13 of the notes, with the same approximations in the RT and in the chi).
    """
    d = notes.hsex_dense_data()
    taus, Em = d['taus'], d['Em']
    rt, chi = 1e3 * d['dim_rt'], 1e3 * d['dim_chi']['accumulated']
    i = np.argmin(np.abs(taus - tau_cut))
    v = max(np.abs(rt).max(), np.abs(chi).max())
    fig = plt.figure(figsize=(PAGE, 2.7))
    gs = fig.add_gridspec(1, 5, width_ratios=[2, 2, 0.09, 0.42, 1.2], wspace=0.08)
    axes = [fig.add_subplot(gs[0])]
    axes += [fig.add_subplot(gs[k], sharey=axes[0]) for k in [1, 4]]
    for ax, y, t in zip(axes[:2], [rt, chi], ['RT (even part)', r'$\chi$, accumulated envelope']):
        pc = ax.pcolormesh(taus, Em, y.T, cmap='bwr', norm=TwoSlopeNorm(0, -v, v), shading='gouraud', rasterized=True)
        ax.axvline(taus[i], color='k', lw=0.7, ls='--')
        ax.set_xlabel(r'delay $\tau$ (fs)')
        ax.set_title(t, fontsize=8)
    cb = fig.colorbar(pc, cax=fig.add_subplot(gs[2]))
    cb.set_label(r'$\Delta$Im $\varepsilon$ ($10^{-3}$)', fontsize=7)
    axes[0].set_ylabel('probe energy (eV)')
    axes[2].plot(rt[i], Em, 'k', lw=1.2, label='RT')
    axes[2].plot(chi[i], Em, 'b--', lw=1, label=r'$\chi$')
    axes[2].axvline(0, color='gray', lw=0.5)
    axes[2].set_xlabel(r'$\Delta$Im $\varepsilon$ ($10^{-3}$)')
    axes[2].set_title(r'$\tau$ = %+.2f fs' % taus[i], fontsize=8)
    axes[2].legend(frameon=False, loc='lower right')
    for ax in axes[1:]:
        ax.tick_params(labelleft=False)
    axes[0].set_ylim(Em[0], Em[-1])
    fig.savefig(os.path.join(out, 'maps.jpeg'))
    plt.close(fig)


# name: function(out), filled as the figures of the article are defined (Outline.md, "Figures")
FIGURES = {'eq_hsex': fig_eq_hsex, 'chi3_hsex': fig_chi3_hsex, 'dchi_cuts': fig_dchi_cuts, 'maps': fig_maps}

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Figures of the article on LiF')
    parser.add_argument('figures', nargs='*', help='names of the figures (default: all): %s' % ', '.join(FIGURES))
    parser.add_argument('--out', default=os.path.join(HERE, 'Figures'), help='output folder of the jpeg files')
    args = parser.parse_args()
    out = os.path.expanduser(args.out)
    os.makedirs(out, exist_ok=True)
    for name in args.figures or FIGURES:
        FIGURES[name](out)
