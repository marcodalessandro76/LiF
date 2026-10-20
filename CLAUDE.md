# LiF pump-probe project — notes for Claude Code

The user (Marco D'Alessandro) writes in Italian: reply in Italian. Code, docstrings and commit messages stay in English.
This is a research project: the work is on notebooks and analysis, not on the development of MPPI. MPPI (the python
package used for the computations and the analysis) has its own repository and its own notes
(`D:\Projects\Research\MPPI\CLAUDE.md`): changes to MPPI classes go there, in a separate session.

## Scientific goal
Compute the transient (non-equilibrium) absorption of bulk LiF in a pump-probe setup from the equilibrium non-linear
susceptibilities computed with `yambo_nl`, and compare it with the direct real-time (RT) simulations of Davide Sangalli
(one RT run for each pump-probe delay). Pump: IR, 1.55 eV (~800 nm). Probe: XUV/broad (around the F 2p / Li 1s
transitions, the analysis window is ~10-20 eV in the yambo_nl runs).

Theory (the reference PDFs of the work go in `References/`, tracked in git):
- `References/Pp_and_non_linear_functions.pdf` (Attaccalite, Sangalli; draft): chi_eff = P(w)/E_p(w);
  chi^neq,(1)[E_P](w) = (P_Pp(w) - P_P(w))/E_p(w) (pump-probe minus pump-only polarization, divided by the probe), and
  its expansion in equilibrium susceptibilities, Eq. (11):
  chi^neq - chi^eq,(1) = [chi2(w+wP,-wP) E_p(w+wP)/E_p(w) + chi2(w-wP,wP) E_p(w-wP)/E_p(w)] E_P(wP)
                       + [chi3(w,wP,-wP) + chi3(w-2wP,wP,wP) E_p(w-2wP)/E_p(w) + chi3(w+2wP,-wP,-wP) E_p(w+2wP)/E_p(w)] E_P^2(wP) + O(E_P^3)
  For a delta probe the ratios E_p(w')/E_p(w) are 1 (Eq. 12); for a monochromatic probe only chi3(w;w,wP,-wP)
  survives (Eq. 13); for a broad pulsed probe the shape of the probe enters. The note itself says that the frequency
  arguments and the treatment of the +-wP combinations still have to be fixed (red N.B.), and it uses real cos
  amplitudes (E_w = E_-w), not the convention of MPPI.
- `NL_Chi/Attivita nuova.txt`: the original plan (starting SAVE of Sangalli, how the YamboPy
  `o-*.YamboPy-SF_probe_order_n_m` files store chi at n times the probe and m times the pump frequency).
- `References/Report_LiF_Davide.pdf`: Davide's report on the RT transient absorption of LiF.
- `References/RT_Analisi_LiF.txt` (old, Jan 2026): plan of the RT transient absorption runs (Sangalli's folders and
  yambo_nl executable, probe-only / pump-probe / pump-only runs, pulse parameters, delay at pump-probe overlap).

### Translation of Eq. (11) into MPPI quantities (MPPI >= 1.3)
`Xn_frequency_mixing(data, X_order=(1,2)).compute_Xn()[0][(n,m)]` at probe frequency w is the ratio between the
component of P at n*w + m*wP and the product of the field amplitudes E(+-w)^|n| E(+-wP)^|m| (signs of n and m), with
E(w) = (i E0/2) exp(i w t0), E(-w) = E(w)^* (Boyd, Nonlinear Optics, Sec. 1.3). It therefore INCLUDES Boyd's
degeneracy factor D: chi(1,+-1) = 2 chi2, chi(1,+-2) = 3 chi3, third order part of chi(1,0) = 6 chi3(w;w,wP,-wP)|E_P(wP)|^2.
In these quantities the pump induced change of the response at the frequency w is

  dchi(w) = sum_{m=+-1,+-2} chi_(1,m)(w - m wP) * E_p(w - m wP)/E_p(w) * E_P(sign(m) wP)^|m|  +  [chi_(1,0)(w) - chi1(w)]

i.e. the key (1,m) must be read at the probe frequency w - m wP (shifted grid), the pump enters with its complex
amplitude E_P(+wP) = (i E0P/2) exp(i wP t0P) and E_P(-wP) = E_P(wP)^* (the delay enters through t0P), and
chi_(1,0) - chi1 already contains |E_P|^2 (chi1 from a probe-only run). To be checked against `eval_dchi_neq`
(see the open points).

## Repository and data
- Laptop: `D:\RICERCA\DFT AND MANY BODY\SIMULATIONS\LiF` (path with spaces: quote it). Cluster: `~/work/LiF` on ismhpc.
  Same git repo (GitHub `marcodalessandro76/LiF`, branch `main`), synced only through git (commit/push on one
  machine, pull on the other). In git there are only README, this file, `Attivita nuova.txt`, the notebooks, `References/` and the
  technical notes `Notes/` (LaTeX sources and figures):
  the data (~11 GB on the cluster) are NOT in git and live only on the cluster (the yambo runs) or on the laptop
  (Davide's reports and tarballs in `NL_Chi/Davide_Google_drive`, ~320 MB).
  Never `git add` data folders, `.ipynb_checkpoints` or tarballs (see `.gitignore`).
- `NL_Chi/` (current work):
  - `YamboNL_Analysis.ipynb`: the yambo_nl datasets (MPPI YamboInput/YamboCalculator/Dataset, slurm, all appended
    with `skip=True`), no analysis: inversion check of the k sampling (`YamboDftParser.expand_IBZ_kpoints`), kx8 delta,
    delta runs on kx12/kx16, damping 0.3 eV sine + P&p runs on kx8/kx12 (`study_eta`) and kx16, the 39-frequency test
    grid, kx12 pump tests (1.6e7 at 39 frequencies, 4e6 at 155), restart and chunk repair of the corrupted frequencies.
    The sections of the runs deleted on 2026-10-12 (damping 0.1 eV production, chi(1,+-1) diagnostics, kx8 intensity
    scan, kx12 4e6 at 39 frequencies) were removed (git history, e.g. 3fd6f06).
  - `NL-Chi_Analysis.ipynb` (simplified 2026-10-12): linear response from the delta runs (kx8/kx12/kx16, eta 0.3 and
    0.5 eV: kx8 not adequate), chi from the damping 0.3 eV sine and P&p runs on kx12/kx16 (linear check, non-linear
    keys and their k dependence with 0.5 eV, relative sizes, quadratic/cubic ratio R and I_eq), kx12 pump intensity
    test (1e6, 4e6, 1.6e7), conclusions. The old diagnostics of chi(1,+-1) (kx8, damping 0.1 eV: X_order (1,3), raw
    spectrum of P(t), intensity scaling, integrator, first k tests), the old k convergence and the broadening
    validation on the oscillator were removed (summarized in one markdown cell; full version in git commit 3fd6f06).
  - `Dchi-neq_Analysis.ipynb` (2026-10-12): RT IP runs of Davide (reading, spectra, maps), kx16 damping 0.3 eV chi with
    `compute_Xn(broadening=0.3)` (total 0.6 eV = RT smoothing), formula and function `dchi_neq` (third order terms only,
    replaces `eval_dchi_neq`, which read the (1,+-2) keys on the wrong shifted grid):
    dchi(w,tau) = |E_P(tau)|^2 [x3(w) - chi_(1,2)(w-2wP) e^{-2i wP tau} - chi_(1,-2)(w+2wP) e^{2i wP tau}],
    |E_P|^2 = E0P^2/4 exp(-tau^2/sigma^2) (delta probe, QSFIELD PHPI90 pump = -E0 f(t) sin(wP(t-t_P)): minus sign from
    (+-i)^2; quasi-monochromatic pump, envelope at the probe time), w in 13.10-21.82 eV (exact index shift 32 = 2wP/step)
    (INVINT phase and energy shift corrected since 2026-10-13, see below). First
    comparison: good below 17.5 eV (bleach 17 eV ~1.5x larger than RT), shapes differ above.
    2026-10-13: a rigid shift of ~0.3 eV to the right of the chi spectra improves the agreement (rms 0.93 -> 0.44 in
    13.6-17.5 eV); it is static (same at all delays) and already in the equilibrium spectrum (edge of Im eps 0.26 eV
    below the RT earliest delay). QSFIELD PHPI90 checked in the Lumen source (`RT_small_a.F`): E = -E0 f [sin(wP tau) +
    tau cos(wP tau)/(wP sigma^2)], so the formula (sine, minus sign) is right. Same scissor (Davide's folder
    `sh6.0+str0.05`). IP optics runs with `yambo -o c` (section "IP linear response with yambo optics" at the end of
    `YamboNL_Analysis.ipynb`, `study_ip_bands`, kx16, bands 3-6/3-15/1-15, damping 0.6 eV, run in
    `kx16_nb100_NoTr_E100/ip-optics-bands_*-damp_0.6`): optics 1-15 = RT earliest delay (edge 0.03 eV: yambo 5.3 and
    Lumen agree); bands 7-15 move the edge by 0.07 eV only (bands 1-2 nothing) but add absorption above 19 eV and +0.5 to
    Re eps. `invint_correction` = chi(w_eff) exp(-i w_eff dt/2) (energy map + phase of the field delay dt/2, the phase is
    the larger effect): yambo_nl chi1 = optics 3-6 (edge 0.01 eV, eps 2%; delta run 0.8%). The equilibrium shift is
    solved. dEps: rms 0.93 -> 0.68, residual best shift 0.34 -> 0.22 eV, now dependent on the pump envelope (0.15-0.27
    near the maximum, 0.04-0.14 in the tail) with RT/chi3 amplitude 0.6-0.85 at the maximum: likely beyond third order
    (pump induced edge shift comparable to the 0.6 eV broadening; RT pump 1500x our highest perturbative test).
    2026-10-14, our RT runs (`yambo_nl`, Davide's IP setup, section "RT pump-probe runs" of `YamboNL_Analysis.ipynb`,
    `kx16_nb100_NoTr_E100/rt-ip-*`, 26 runs: bands 3-6 with pump 6e9/6e8/6e7 at t_p 13-19 fs, bands 1-15 6e9 at 16 fs
    and 6e7 at 14.5-17.5 fs, pump-only and probe-only references; analysis in the last section of
    `Dchi-neq_Analysis.ipynb`): our 1-15 6e9 run = Davide's (3e-7); probe-only CRANKNIC = optics (0.8%). Bands 1-15 vs
    3-6: same dEps in 13.6-17.5 eV (rms 0.03-0.08), very different above 18 eV. Intensity: 6e9 is beyond third order
    (fifth order 25-45% of the third near the maximum), 6e7 has the spurious linear-in-E_P term (k discretization)
    comparable with the third order: the orders are separated at each frequency with the three intensities
    (a e + b e^2 + c e^4). The third order part b agrees with chi3 near the pump maximum (tau -1.41/0.09 fs: rms
    0.28/0.19, amplitude 1.04/1.16, shift <= 0.1 eV, also above 17.5 eV): the residual 0.2 eV shift is from the orders
    above the third. Away from the maximum chi3/b = 0.6 (tau -2.9) and 1.4-1.7 (+1.6, +3.1): limit of the formula with
    the envelope at the probe time (probe polarization lives ~1.1 fs, sigma 3.2 fs); next: a formula with the time
    dependence of the envelope, chi with more bands for 18-22 eV.
    Reversed pump (`study_rt_neg`, bands 3-6, 6e7 along -x, `pump_sign=-1` in `build_rt_input`/`rt_id`, ids `-EP_6e7-neg`):
    the even part [dEps(+E_P)+dEps(-E_P)]/2 is the clean third order RT (odd part = spurious second order, as large as
    the even one at 6e7). It equals the fit b for tau >= 0 (rms 0.01-0.06), and chi3 reproduces it with no shift
    (<= 0.06 eV near the maximum) over 13-22 eV: rms 0.17/0.23, amplitude 0.95/1.20 at tau -1.41/0.09 fs; chi3/RT 0.69
    at -2.9 fs and 1.44/1.70 at +1.6/+3.1 fs (envelope approximation of the formula).
    Envelope correction (2026-10-14, section "Time dependence of the pump envelope", `dchi_neq_t`): the formula above
    takes chi3 at Omega = 0 (envelope at the probe time); the pump acts during the life of the probe coherence (~1.1 fs).
    Model: transform of the constant-envelope dchi to s = t - t_p, times W(tau,s) = (1/s) int_0^s f^2(tau+s') ds'
    (accumulated envelope: adiabatic Stark/Franz-Keldysh shift accumulated on the coherence), back to w; same W applied
    to the (1,+-2) terms (an extension, not derived). chi3 vs RT third order (even 6e7) at all 5 delays: amplitude
    1.04-1.07 (was 0.69-1.70), rms 0.17-0.26; the local weight f^2(tau+s) overcorrects. Reversed pump at 6e9 (6 more
    runs `-EP_6e9-neg`): odd part 6-13% of the even one (spurious term small); the even part at 6e9 departs from the
    third order up to +1.6 fs (rms 0.61-0.25, amplitude 0.74-0.82, spectrum up to 0.13 eV higher), coincides at +3.1 fs:
    the RT of Davide is beyond the perturbative regime near/before the pump maximum.
    Both notebooks come from the split (2026-10-05) of the old `Transient_Abs_NL-Chi.ipynb` (in the git history).
  - cluster only, `kx8_nb100_NoTr_E100/`, `kx12_nb100_NoTr_E100/`, `kx16_nb100_NoTr_E100/`: copies of Sangalli's
    NoTr_E100 SAVEs (`/work/sangalli/simulations/LiF/yambo/e80_kx6_pbe_sr/kx{8,12,16}_nb100/NoTr_E100/SAVE`: same
    lattice, 100 bands, the same 8 symmetries without inversion and time reversal, 100/294/648 IBZ points), copied
    without the yambo 5.2 setup databases and with the setup redone (`yambo` without arguments) with Lumen 2.1.0.
    Runs kept after the cleanup of 2026-10-12 (~8.5 GB deleted on /work, not in the home quota): delta
    (`lresponse-bands_3-6-delta`, NLstep 0.01 fs) on the three grids; damping 0.3 eV sine and P&p 1e6
    (`...nlenrange_10.0-25.016-nlensteps_155-...-damp_0.3-sin-nltime_100`) on kx8 (1 node), kx12 (P&p = the 2-node
    run) and kx16 (repaired); kx12 P&p 4e6 (155 frequencies) and 1.6e7 (39 frequencies, with its `_corrupted` and
    `_restart_NL_CPU_4.8` fragments); one corrupted run, the kx12 1-node P&p `..._corrupted` (case B of
    `YamboNL_Corrupted_Frequencies.ipynb`); the `_corrupted*` fragments inside the kx16 runs (cases D-F); the job outputs
    of the deleted 2-node runs (`job_..._2nodes.out`, kx12 `job_pulse..._corrupted.out`: MXM segfaults). Deleted: the
    damping 0.1 eV production and 5-frequency diagnostic runs (kx8, kx12), the kx8 intensity scan, the kx12 4e6 run at
    39 frequencies, the 2-node kx8 runs and the kx12 2-node sine `_corrupted` (case A, results kept as text).
  - cluster only, `Davide_data/`: Davide's RT results at kx16 in `20250312_Report_simulations_data_kx16_ball/` (plus
    tarballs): one folder per probe time `delta<t_p>fs` with `o-abs_0-50eV_sm0.6eV.YPP-eps_along_E` (E, Im eps x/y/z,
    Re eps x/y/z; ypp, Lorentzian smoothing 0.6 eV). Two sets: `td-hsex-cv-ber_..._3-12b` (HARTREE+SEX, 61 delays,
    10-60 eV) and `td-ip-ber_..._1-15b` (IP, 62 delays, 0-50 eV); yambo 5.3, CRANKNIC 3 as, 40 fs, no dephasing; probe
    DELTA along x; pump QSFIELD PHPI90 5SIGMA 1.56 eV, 6e9 kW/cm^2, sigma 3.219 fs, switched on at 0.003 fs. Pump
    maximum at t_P = 0.003 + (pi/wP) nint(5 sigma wP/pi) = 15.909 fs (yambo `RT_small_a`), delay tau = t_p - t_P
    (-7.41..+7.34 fs, step 0.25 fs; probe times read from the headers). No probe-only run: the reference is the earliest
    delay (t_p = 5.25 fs, tau = -10.66 fs), as in Davide's plotting script (`20250312_Report_simulations_data_kx16.py`:
    t0 = 15.9 fs, reference = first delay, Fresnel R at 72 deg; his "Neq." curve and right panel are averages over
    |tau - 0.1| <= 2 fs, his left panel shows -X). `Dchi-neq_Analysis.ipynb` uses only the IP set (`load_rt_runs`,
    `rt_quantity`, `plot_rt_spectra`, `plot_rt_map`; the dR/R map at 72 deg s-pol reproduces his figure). N.B. the RT
    pump is at 1.56 eV, the yambo_nl pump at 1.55 eV.
- `RT_Transient_Absorption/` (earlier work, transient absorption from RT simulations): `Transient_Absorption.ipynb`
  (fit of the experimental IR pump, MPPI Gaussian pulses vs ypp_rt fields, yambo_nl RT pump-probe dynamics with a 1.56 eV pump of fwhm 7.58 fs and a
  35 eV probe of fwhm 0.42 fs, transient spectra and reflectivity, delays 7.0 and 9.12 fs) and `Sum_frequency.ipynb`.
  Cluster only: `NoTr_E100` (SAVE + collisions `coll-par-cv_1-15B_X59RL-50B_noSM_H339RL_kpar`, fixsym for the 100
  field), `Field_test`, `nl_input_template`, three `td-hsex-*` runs.
- The old version of MPPI's `Analysis_Optics.ipynb`, based on the LiF runs (now replaced in MPPI by an analytical model),
  is in the MPPI git history: `git -C <MPPI> show b785feb:sphinx_source/tutorials/Analysis_Optics.ipynb`. It can be
  moved here.

## Technical notes (`Notes/`, LaTeX, started 2026-10-14)
Technical notes on this activity: NOT the paper, but a complete and as concise as possible explanation of what is
done (theory, implementation, computations) and of the notation and conventions. Written in English, in LaTeX.
- Files: `Notes/main.tex` (main file: `\documentclass[10pt,a4paper]{paper}`, then `\input{header.tex}`),
  `Notes/header.tex` (packages and new commands: put there the macros of the recurring symbols, e.g. the MPPI keys
  chi_(n,m), Delta eps, E_P, the pump envelope f, so that the notation stays uniform), figures in `Notes/Figures/`
  as jpeg (`\includegraphics{Figures/<name>.jpeg}`). Style of the user's example notes: header with inputenc utf8,
  babel (english here), amsmath/amsfonts/amssymb, bm (`\bm{}` for vectors), geometry (margin 1in), booktabs,
  graphicx, subcaption, tikz, `\parindent` 4pt; `\section`/`\subsection`, numbered equations with `\label{eq:...}` and
  `\eqref`, figures `[h!]` with `\label{fig:...}` and `Fig.~\ref{...}`, subfigures side by side for panels.
- Workflow: the notes are written on the laptop (LaTeX is there, not on the cluster) with the Claude of VS Code and
  compiled locally (the PDF is not tracked, `*.pdf` is in `.gitignore`; the LaTeX auxiliary files too). The data and
  the notebook executions are only on the cluster: the figures are made there (Claude CLI session in tmux `lif`, from
  the notebooks/data with a script kept in `Notes/`), saved as jpeg in `Notes/Figures` and pushed; the laptop session
  asks for them through the user (or checks small things with ssh, see "Cluster ismhpc"). Sync only through git, with
  a file edited on one side at a time.
- User's guidelines (2026-10-15): notes for the user and the colleagues (technical, explicit derivations). Section 1
  is the theory (P vs chi, delta probe, chi^neq, dchi for monochromatic and pulsed pump, accumulated envelope), written
  without MPPI notation (chi_0 = 6 chi3(w;w,wP,-wP), chi_+- = 3 chi3(w;w-+2wP,+-wP,+-wP)): the MPPI keys and the
  implementation come in the later sections. Vectors in bold (`\bm`) with Cartesian indices until the fields are
  taken along x. Only a delta probe (no monochromatic probe). Display equations must fit the line (split with align,
  check the log for Overfull \hbox). Bibliography in `Notes/biblio.bib` (bibtex, unsrt); do NOT cite the
  Attaccalite-Sangalli draft (these notes replace it). After each change Claude builds `Notes/main.pdf` in place
  with latexmk (user's request, 2026-10-19) and checks the log (errors, Overfull \hbox, undefined references); if the
  PDF is locked by the viewer (Acrobat: "I can't write on file main.pdf"), ask the user to close it. After a failed
  run latexmk does not recompile unchanged sources and reports the old errors (the log is not rewritten): use
  `latexmk -g` (force) once the lock is gone.
- Sources: this file (conventions, results, numbers) and the markdown cells of the notebooks (tracked in git, they can
  be read on the laptop): `Dchi-neq_Analysis.ipynb` ("Formula used for dchi^neq", "Frequency grids", "Correction of
  the INVINT integrator", "Our RT runs with the bands 3-6", "Pump intensity: third order RT response", "High pump
  intensity: orders above the third", "Time dependence of the pump envelope"), `NL-Chi_Analysis.ipynb` (chi
  extraction, chi(1,+-1), k convergence, pump intensity test), `YamboNL_Analysis.ipynb` (inputs of all the runs),
  `YamboNL_Corrupted_Frequencies.ipynb`; `References/Pp_and_non_linear_functions.pdf` (Attaccalite-Sangalli, Eq. 11),
  `References/Report_LiF_Davide.pdf`; Boyd, Nonlinear Optics, ch. 1 (conventions of the fields and of chi).
- Outline agreed as a starting point (to be refined with the user):
  1. Introduction: transient absorption of LiF from the equilibrium non-linear chi (yambo_nl) vs direct RT pump-probe
     runs (D. Sangalli's and ours).
  2. Conventions: E(t) = sum_n E(w_n) e^{-i w_n t}, E(-w) = E(w)^*, sine field E0 sin(w(t-t0)) -> E(w) = (i E0/2)
     e^{i w t0}; MPPI keys chi_(n,m) with Boyd's degeneracy factors (chi_(1,+-1) = 2 chi2, chi_(1,+-2) = 3 chi3,
     third order part of chi_(1,0) = 6 chi3(w;w,wP,-wP)|E_P|^2); units (yambo kWLm2 = kW/cm^2, E0 = sqrt(4 pi I/c)
     au, 1e6 kW/cm^2 -> 1.194e-4 au); Delta eps = 4 pi dchi.
  3. Theory: chi^neq expansion (Eq. 11 of the reference), exact third order expression for a delta probe
     (dchi(w,t_p) = 3 int dOmega e^{-i Omega t_p} int dw1 chi3(w; w-Omega, w1, Omega-w1) E_P(w1) E_P(Omega-w1), three
     groups Omega ~ 0, +-2wP), monochromatic limit = the formula with the envelope at the probe time
     (dchi = |E_P(tau)|^2 [x3(w) - chi_(1,2)(w-2wP) e^{-2i wP tau} - chi_(1,-2)(w+2wP) e^{2i wP tau}], minus sign
     from the sine pump), accumulated envelope model (W(tau,s) = (1/s) int_0^s f^2(tau+s') ds', Stark shift
     accumulated on the probe coherence; extension to the (1,+-2) terms not derived).
  4. Computation of chi with yambo_nl: SAVE and setup (kx16, bands 3-6, no inversion/time reversal symmetries,
     Lumen 2.1.0), sine and pump&probe runs (155 frequencies 10-25 eV, step wP/16, damping 0.3 eV, NLtime 100 fs,
     NLstep 0.01 fs, pump 1e6), harmonic analysis (`Xn_frequency_mixing`, X_order (1,2)), a posteriori broadening
     0.3 eV, INVINT correction, shifted frequency grids (2wP = 32 steps).
  5. Numerical issues: spurious chi_(1,+-1) (k discretization of the Berry coupling, decreases with the k grid,
     dropped), k convergence of the cubic terms, pump intensity test (4e6 perturbative, 1.6e7 not), silently corrupted
     frequencies of yambo_nl (short: detection and repair).
  6. RT runs: setup (CRANKNIC 3 as, 40 fs, no damping; DELTA probe; QSFIELD PHPI90 5SIGMA pump 1.56 eV, sigma 3.219 fs,
     t_P = 15.909 fs; delays), Delta eps = eps[P_Pp - P_P] - eps_eq with the 0.6 eV Lorentzian smoothing of the FT
     (a window of the analysis, not a dephasing), check vs D. Sangalli and vs yambo optics, bands 3-6 vs 1-15, pump
     intensities 6e9/6e7, reversed pump and even/odd parts.
  7. Results: chi3 vs RT third order (accumulated envelope: amplitude 1.04-1.07, rms 0.17-0.26), 6e9 beyond the
     perturbative regime (orders above the third, odd part small).
  8. Open points (k convergence of the cubic terms, bands above 6 for 18-22 eV, residual 5% / 0.07 eV).
- State (2026-10-19): Sec. 1 (theory) of `Notes/main.tex` is done and agreed with the user; Sec. 2 (extraction of the
  chi from the RT polarization with sine fields: switch-on and damping, harmonic fit for one and two fields, mixing
  keys chi_(n,m), field amplitudes; a general description of the method of MPPI, no implementation details) is done.
  Sec. 3 "Transient absorption of LiF" (structure given by the user, 2026-10-19): 3.1 RT transient absorption (3.1.1
  setup of our RT runs at 6e7 kW/cm^2, kx16, bands 3-6, with the figure of the pump and of the probe times; 3.1.2 even
  and odd parts in the pump field with the reversed pump, origin of the spurious second order term, the comparison
  with the chi uses the even part), 3.2 Non-linear chi analysis (done: extraction from the keys, Eq. chi_from_keys;
  3.2.1 convergence with the k grid, Fig. chi_kgrid; 3.2.2 convergence with the pump intensity, kx12, five intensities
  1e6-2.56e8 kW/cm^2, Fig. chi_pump and Tab. chi_pump with the maximum changes of chi_0, chi_+, chi_- in 10-20 eV of
  the probe grid: 0.20/0.41/0.18% at 4e6, 0.91/1.84/0.73% at 1.6e7, 3.89/7.45/3.11% at 6.4e7, 14.3/25.2/11.4% at
  2.56e8, growth ~ I_P up to 6.4e7 (fifth order), less at 2.56e8; the RT runs at 6e7 have a few % beyond the third
  order near the pump maximum), 3.3 transient absorption from the chi and comparison with the RT results (done
  2026-10-20: Fig. dchi_rt = RT even part at 6e7 vs 4 pi dchi from the kx16 chi at the 5 delays, envelope at the probe
  time and accumulated; details: window 13.1-21.8 eV, keys with 2wP = 3.10 eV and phases with 1.56 eV, derivative term
  of the pump neglected, INVINT correction (with it 1+4pi chi1 = eps_eq of the RT probe-only run: 4.4% rms in 10-25 eV,
  edge -0.015 eV; without it 15.6% and 0.29 eV); in 13.6-21.3 eV rel. rms / ratio of the maxima: probe time
  0.17-0.77 / 0.69-1.70, accumulated 0.17-0.26 / 1.04-1.07, local weight 0.65-1.38; residual ~5% and a shift of
  0.06-0.07 eV at/after the maximum; 3.3.1 third order regime, 3.3.2 high pump intensity: Fig. dchi_rt_6e9 = RT even
  part at 6e9 vs chi with the accumulated envelope only (user's choice): rel. rms 0.96/0.76/0.61/0.42/0.32, chi larger
  by 12-42%, 0.12-0.16 eV lower; the RT 6e9 vs the third order RT (6e7 even x 100, printed only): rms 0.83 -> 0.09,
  20-25% smaller and 0.14 eV higher up to the maximum, third order again at +3.09 fs). Nothing after Sec. 3 is defined:
  do not refer to sections that do not exist. Figures: `Notes/make_figures.py` (registry FIGURES, run on ismhpc from
  ~/work/LiF/Notes with `--out ~/tmp_claude/notes_figures`, then scp the jpeg to Notes/Figures on the laptop and commit
  there; the cluster copy of the repo never commits). Our RT runs at 6e7: only 5 probe times (13-19 fs, tau -2.91 ...
  +3.09 fs); more delays later (each run ~50 min on 1 node, two runs per delay with the reversed pump). A LaTeX comment at the
  end of main.tex keeps notes for later (comparison of the envelope weights; relation between the keys chi_(1,m) and
  the chi of Sec. 1, shifted grids, chi0 from the difference of two runs).

### Reminder for the article (2026-10-20, NOT for the technical notes)
The technical notes are for the discussion with the colleagues; the user does not want these elements there. To be
included in the document that prepares the article:
1. Introduction (goal, the two approaches, LiF, IR pump and XUV probe).
2. Computational details (table or appendix): DFT (functional, lattice, 100 bands, scissor 6.0 eV + stretch 0.05),
   k grids (IBZ points, the 8 symmetries without inversion/time reversal), codes and versions (Lumen 2.1.0, yambo 5.3,
   MPPI), parameters of the chi runs and of the RT runs, intensity -> amplitude E0 = sqrt(4 pi I/c) (1e6 kW/cm^2 ->
   1.194e-4 au); the notes never name the codes.
3. Conclusions and open points (k convergence of the cubic terms above 16 eV, bands above 6 for 18-22 eV, more delays
   at 6e7, the +-2wP terms of the envelope model not derived, residual 5% / 0.07 eV, higher order chi at 6e9).
4. Relation with the simulations of D. Sangalli: our 1-15 6e9 run = his (3e-7); bands 1-15 vs 3-6 (same below 17.5 eV,
   different above 18 eV).
5. Validation of the extraction: AnharmonicOscillator (1e-6-1e-3), chi1 sine vs delta (1-2%).
6. A posteriori broadening: chi(w + i Delta), only the terms linear in the probe (now a comment at the end of main.tex).
7. Silently corrupted frequencies of yambo_nl (appendix, or only the dedicated notebook).
8. Outlook: Delta R/R (the experimental observable, Davide's 72 deg), excitonic (HSEX) runs.

### Model of Sec. 1 of the notes (reference for checking and implementing it)
Atomic units, E(t) = int dw E(w) e^{-iwt} (E(w) = (2pi)^-1 int dt E(t)e^{iwt}), fields along x, delta probe
E_p(t) = E0p delta(t-t_p) -> E_p(w-Omega)/E_p(w) = e^{-i Omega t_p}. chi^neq = (P_Pp - P_P)/E_p = chi1 + dchi with
dchi^(2) = 2 int dOmega e^{-i Omega t_p} chi2(w;w-Omega,Omega) E_P(Omega) and
dchi^(3) = 3 int dOmega e^{-i Omega t_p} int dw1 chi3(w;w-Omega,w1,Omega-w1) E_P(w1)E_P(Omega-w1).
- Monochromatic sine pump E0P sin(wP(t-t_P)), E_P(+-wP) = +-(iE0P/2)e^{+-iwP t_P}, tau = t_p - t_P:
  dchi = (iE0P/2)[chi2_+ e^{-iwP tau} - chi2_- e^{iwP tau}] + (E0P^2/4)[chi_0 - chi_+ e^{-2iwP tau} - chi_- e^{2iwP tau}],
  chi2_+-(w) = 2chi2(w;w-+wP,+-wP), chi_0(w) = 6chi3(w;w,wP,-wP), chi_+-(w) = 3chi3(w;w-+2wP,+-wP,+-wP). In MPPI keys:
  chi_0 = (chi_(1,0)-chi1)/(E0P_run^2/4), chi_+-(w) = chi_(1,+-2)(w-+2wP), chi2_+-(w) = chi_(1,+-1)(w-+wP).
- Gaussian pulse E0P f(t-t_P) sin(wP(t-t_P)), f(t) = e^{-t^2/2sigma^2}: E_P = E^+ + E^-, E^+-(w) = +-(iE0P/2)
  (sigma/sqrt(2pi)) e^{iw t_P} e^{-(w-+wP)^2 sigma^2/2}. Second order: exact, two terms with E^+-. Third order: chi3 taken
  out of the integral over w1 at w1 = +-wP (pump arguments off resonance), Omega dependence kept:
  dchi^(3) = int dOmega e^{-i Omega t_p}[6chi3(w;w-Omega,wP,Omega-wP)S_0 + 3chi3(w;w-Omega,wP,Omega-wP)S_+2
  + 3chi3(w;w-Omega,-wP,Omega+wP)S_-2], S_0 = int dw1 E^+(w1)E^-(Omega-w1) = (E0P^2/4)(sigma/2sqrt(pi)) e^{i Omega t_P}
  e^{-Omega^2 sigma^2/4}, S_+-2 = -(same with Omega -+ 2wP).
- Envelope at the probe time (chi independent of Omega over hbar/sigma): the monochromatic result with f(tau) (second
  order) and f^2(tau) (third order). Valid if hbar/eta << sigma (pump envelope constant during the probe dynamics,
  eta = width of the transitions = total broadening, 0.6 eV in our comparison); underestimates before the pump
  maximum, overestimates after it.
- Accumulated envelope (third order): independent two-level transitions, pump shift d eps(t) = c_eps (E0P^2/4) f^2(t-t_P)
  (dynamical Stark / Franz-Keldysh), response function R_eps(t) = i|d|^2 e^{-i eps t - eta t}; accumulated phase ->
  dR_eps(t) = -i c_eps (E0P^2/4) R_eps(t) int_0^t f^2(tau+t')dt'; the unknown sum_eps(-i c_eps)R_eps(t) = R_0(t)/t is
  fixed by requiring the monochromatic chi_0 (calibration, not a physical case), with chi_0(w) = int_0^inf dt R_0(t)
  e^{iwt}, R_0(t) = (2pi)^-1 int dw chi_0(w) e^{-iwt}. Result dR = (E0P^2/4) W(tau,t) R_0(t),
  W(tau,t) = (1/t) int_0^t f^2(tau+t')dt' = (sqrt(pi) sigma/2t)[erf((tau+t)/sigma) - erf(tau/sigma)], and
  dchi^(3)(w,tau) = (E0P^2/4) int_0^inf dt e^{iwt} W(tau,t)[R_0 - R_+ e^{-2iwP tau} - R_- e^{2iwP tau}] (same W on the
  +-2wP terms: extension, not derived). Limits: W -> 1 (sigma -> inf, monochromatic), W -> f^2(tau) (hbar/eta << sigma).
  Equivalent in frequency: 6chi3(w;w-Omega,wP,Omega-wP) = (1/Omega) int_{w-Omega}^{w} chi_0(w')dw' (version discarded
  from the notes: the user found the time form clearer; useful as an independent numerical check).

### To do after the notes (not started: the user wants to go on with the notes first)
1. Check that `dchi_neq` and `dchi_neq_t` (+ `envelope_weight`, `envelope_conv`) of `Dchi-neq_Analysis.ipynb` implement
   the model above: chi_0 = x3 = (chi_(1,0)-chi1_sin)/(EP**2/4) with the same broadening on both; shifted grids
   (n = 2wP/step = 32, keys (1,+-2) sliced [:-2n] / [2n:]); signs and phases of the +-2wP terms (minus from the sine
   pump, phases with the RT pump frequency `wP_phase`); f^2 = exp(-tau^2/sigma^2) with sigma from the ypp header (field
   width) and t_P of QSFIELD 5SIGMA; W of `envelope_weight('accumulated')` = Eq. above (for s < 0 it returns f^2(tau));
   `envelope_conv` = inverse DFT on the probe grid, times W, direct DFT (period 2 pi hbar/step = 42.7 fs; for the
   (1,+-2) keys the transform is done on the input grid and then shifted: equivalent because a frequency shift is a
   phase in time, to be verified); the INVINT correction applied to dchi at the output frequency (for the shifted
   keys the probe enters at w-+2wP: check whether the energy map and the dt/2 phase are applied consistently); the
   neglected derivative term of the QSFIELD pump. A numerical check: the frequency form above must give the same result.
   Partial check (2026-10-20): `fig_dchi_rt` of `Notes/make_figures.py` implements the formulas of the notes
   independently (`dchi_from_chi`, `envelope_weight`, `time_weighted`, `invint_correction`) and gives the same numbers
   as `dchi_neq`/`dchi_neq_t` (rms and ratios of the maxima of all the weights to the second decimal).
2. DONE (2026-10-20, MPPI f0f8add): class `Transient_Response` in MPPI Optics, see the handoff below.

### Handoff (2026-10-20): use the MPPI class `Transient_Response` in the LiF notebooks (for the session on the cluster)
Goal: replace the local implementations of dchi (`dchi_neq`, `dchi_neq_t`, `envelope_weight`, `envelope_conv` of
`NL_Chi/Dchi-neq_Analysis.ipynb`; `dchi_from_chi`, `envelope_weight`, `time_weighted` of `Notes/make_figures.py`) with
the class of MPPI, which implements the formulas of Sec. 1 of the notes. MPPI f0f8add is already pulled in
`~/Applications/MPPI` (85 tests pass on the cluster). The tutorial is the last section of MPPI `Analysis_Optics.ipynb`.
- API (MPPI units: energies in eV, times in fs, fields in au; `from mppi import Optics as O`):
  `tr = O.Transient_Response.from_Xn(mix, single, broadening=0.3, centrosymmetric=True)` with `mix` the
  `Xn_frequency_mixing` of the P&p runs (X_order (1,2) or (1,3)) and `single` the `Xn_single_frequency` of the sine runs
  on the SAME probe grid (true for the 155-frequency runs). It computes chi1 and the keys with the same broadening (as
  `chi_kgrids` of make_figures), takes E0P_run and wP (1.55 eV) of the runs, checks that the step divides 2wP (32
  steps) and builds chi_0, chi_+- on the output grid `tr.w` (13.1-21.8 eV). `centrosymmetric=True` drops the keys
  (1,+-1) (the k discretization error). `tr.components()` gives chi_0, chi_+, chi_- on `tr.w`.
  `w, dchi = tr.dchi(taus, E0P, pump='accumulated' | 'probe time' | 'monochromatic', sigma=sigma, wP=wP_RT)`:
  taus = t_p - t_P (fs), E0P the RT pump amplitude (au), sigma the FIELD envelope width (fs, FWHM/(2 sqrt(2 ln 2)) of
  the ypp/yambo header, as in `pump_center`), wP = 1.56 eV of the RT runs for the phases only (the grid shift always
  uses the 1.55 eV of the chi runs, as `wP_phase` did). dchi has shape (len(taus), len(w)), Delta eps = 4 pi dchi.
  Arrays constructor for the cases without a common grid (e.g. the 39-frequency pump tests vs the 155 sine run):
  `O.Transient_Response(freqs_eV, chi1, keys, E0P_run, wP_eV, centrosymmetric=True)` with chi1 at the matching
  frequencies (as in `chi_pumps_kx12`; same broadening on chi1 and keys).
- NOT in MPPI (user's decisions): the INVINT correction (keep `invint_correction` here and apply it to the output as
  now, `4*np.pi*invint_correction(w, dchi)`, dt = 0.01 fs); the derivative term of the QSFIELD pump (neglected); the
  'local' weight f^2(tau+t) (the notes quote its result: keep a local version only if that number is still computed);
  the accumulated envelope for the second order terms (ValueError; irrelevant here, LiF is centrosymmetric). The class
  uses the pump E0P f sin(wP(t-t_P)): the minus sign of QSFIELD PHPI90 does not change the third order.
- Already checked by the MPPI session (2026-10-20, on the cache `~/tmp_claude/notes_figures/chi_kgrids_0.30.npz`,
  kx16, arrays constructor, taus -2.91 ... +3.09 fs, E0P 9.25e-4, sigma 3.219 fs, wP 1.56 eV): the class gives the
  dchi of `dchi_from_chi` on the same output grid (91 frequencies, 13.10-21.82 eV) with max |difference| / max |dchi|
  2e-9 (probe time) and 3.5e-9 (accumulated). So the numbers of Sec. 3.3 of the notes do not change. In the notebook,
  check also `from_Xn` against the current chi (same broadening) and the old `dchi_neq`/`dchi_neq_t`, and report to the
  user before deleting the old functions; keep the RT side (`load_rt_runs`, `rt_eps`, `rt_even`, `compare`,
  `reflectivity`, plots) here.
- Workflow: edit and run the notebooks on the cluster (analysis notebooks may be executed on the login node, see
  "Cluster ismhpc"); commit and push from the cluster only with the user's go; the laptop session does not edit these
  files in the meantime.

### Handoff (2026-10-19): kx12 pump test at higher intensities (for the CLI session on the cluster)
Prepared by the laptop session; the user launches it from the cluster session. Goal: extend Fig. 4 of the notes
(Sec. 3.2.2, chi_0 extracted with different pump intensities, kx12) beyond 1.6e7 kW/cm^2, to see how the fifth (and
higher) order terms projected on chi_0 grow. New section of `YamboNL_Analysis.ipynb` after the kx12 pump test,
"Pump intensity test on the kx12 grid at higher intensities": dataset `study_pump_kx12_high`, two P&p runs on kx12
with the 39 frequencies of the test grid and the pump 6.4e7 and 2.56e8 kW/cm^2 (8 E_P and 16 E_P), `num_tasks=2`
(two nodes on all12h at the same time, 1-2 h each; check that no more than 3 nodes are in use). `fmt_int` now uses 3
significant digits (ids `EP_6.4e7` and `EP_2.56e8`; the ids of the existing runs do not change).
1. `cd ~/work/LiF && git pull`; dry run: `cd NL_Chi && ~/miniconda3/bin/python ~/tmp_claude/launch_dataset.py
   study_pump_kx12_high --dry` (prints the two ids and skip flags, submits nothing).
2. Launch: `nohup ~/miniconda3/bin/python ~/tmp_claude/launch_dataset.py study_pump_kx12_high >
   ~/tmp_claude/pump_kx12_high.log 2>&1 &`; monitor with `squeue -u $USER` and the job outputs.
3. Check the corrupted frequencies of both runs (kx12 runs had silently wrong groups of frequencies, see the section
   "Restart of the corrupted frequencies"): the old criterion |chi_(1,0) - chi1|/|chi1| > 1e-2 does not work at these
   intensities (the third order term is large). Use max|P(t)| of each frequency (`YamboNLDBParser(...).Polarization`):
   it is dominated by the response to the pump and almost the same at all the frequencies; a corrupted group (indices
   with the same remainder modulo the number of `NL_CPU` w-groups, usually from its 3rd frequency on) has it reduced by
   20-50%. Cross-check with chi_0 at the same frequencies of the 1.6e7 run (the change must grow smoothly with I_P).
   Repair with the restart mechanism of that section (`backup_run_files`, `remove_fragments`, `missing_fragments`,
   restart with `skip=False, clean_restart=False, activate_BeeOND=False`), adding the cells to the notebook.
4. Figure: `cd ~/work/LiF/Notes && LIF_NL_CHI=$HOME/work/LiF/NL_Chi ~/miniconda3/bin/python make_figures.py chi_pump
   --out ~/tmp_claude/notes_figures` (it adds the new intensities to the cache `chi_pumps_kx12.npz` and prints the
   median changes of chi_0, chi_+ and chi_- w.r.t. 1e6). Report the numbers to the user; the laptop session copies
   `chi_pump.jpeg` and updates the notes (the right panel may need a log scale with five intensities).
5. If the notebook was modified on the cluster (restart cells), commit and push it from the cluster with the user's
   go; the laptop session does not edit `YamboNL_Analysis.ipynb` in the meantime.

Outcome (2026-10-19, CLI session on the cluster): both runs done (jobs 59718/59719, 1h19). 6.4e7 correct at all the 39
frequencies (max|P| within 0.6% of the neighbours). 2.56e8: group 2 corrupted from its 3rd frequency (indices 10, 14,
..., 38, max|P| 63%), as in the 1.6e7 run; the restart with NL_CPU 4.8 (job 59720) repaired 10 and 14 but the group
was corrupted again from its 3rd frequency (18-38, 75%); the remaining 6 were recomputed by the chunk repair (2 per
job, jobs 59721-3), all correct (within 0.2%). Wrong fragments kept as `_corrupted` / `_corrupted_repair1`, backups
`_backup_orig` / `_backup_repair1`; cells in the section "Repair of the kx12 run with pump 2.56e8" of the notebook
(`chunk_repair` is now defined in a cell of its own, which `launch_dataset.py` always executes). `chi_pump.jpeg` in
`~/tmp_claude/notes_figures` (cache updated). Changes w.r.t. 1e6 in 10-20 eV (chi_0 mean % of the max / median
relative; chi_+ and chi_- median relative): 4e6 0.05/0.19%, 0.25%, 0.55%; 1.6e7 0.25/1.0%, 1.1%, 2.3%; 6.4e7
0.93/3.8%, 3.6%, 6.1%; 2.56e8 3.35/12.5%, 11.0%, 21.4%. From 1.6e7 on a factor 4 in I_P multiplies the change by
3.5-3.8 (~ I_P: fifth order projected on the chi); the change is concentrated in 17-21 eV (peak ~18.5 eV, 14% of the
max at 2.56e8); chi_- changes about twice as much as chi_0 and chi_+. The linear scale of the right panel is still
readable (the 4e6 curve is flat at zero).

## Cluster ismhpc
- Reached from the laptop with `ssh -o BatchMode=yes -o ClearAllForwardings=yes ismhpc '...'` (host in
  `~/.ssh/config`, jump through `narro`; ClearAllForwardings avoids the LocalForward 4444 clash). CentOS 7 /
  glibc 2.17: VS Code Remote-SSH does not work. The connection sometimes times out: just retry later.
- Never run heavy computations on the login node `frontend`: use slurm. At most three nodes in use at the same time (cluster rule,
  2026-10-13): MPPI datasets of 1-node runs with `num_tasks` <= 3, counting the other running jobs. Always use the partition `all12h` (32 cores
  per node), also for the short test/diagnostic runs: it runs on all the nodes and gives no problems (user's choice,
  2026-10-06; the `debug` partition is not used). Runs longer than 12 h (e.g. the kx16 yambo_nl runs on one node,
  ~14 h) use `slownodes` (time limit 2 days, wnode01-02, 04-05, 07-08): `yambo_nl_code(nodes, partition)` in the first
  cell of `YamboNL_Analysis.ipynb`. Home quota 19.5 GB (11 GB used on 2026-10-05: watch the size of new runs,
  BeeOND scratch is used for the runs).
- Python: `~/miniconda3` base, python 3.13; numpy, scipy, matplotlib, netCDF4 and the Jupyter stack are pip-installed
  (update with pip). Notebooks are executed in place with
  `jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.kernel_name=python3 --NotebookClient.record_timing=False <nb>.ipynb`
  The user also works interactively: the LocalForward 4444 of the ssh config is the tunnel used to run a Jupyter server
  on the cluster and open the notebooks in the browser of the laptop.
- MPPI is installed in editable mode from `~/Applications/MPPI` (branch master): `git pull` there to get the latest
  version (laptop copy: `D:\Projects\Research\MPPI`, python `C:/Users/Marco/miniconda3/python.exe`).
- QE and Yambo need different environments (pw.x: Intel MPI + MKL; yambo: OpenMPI), loaded by the module files
  `~/module_script/qe_module` and `~/module_script/yambo_module`, passed to MPPI as `RunRules(pre_processing=...)`
  (included in the slurm scripts, sourced before direct runs, output discarded). Without them pw.x fails with
  `libmkl_gf_lp64.so: cannot open shared object file`. `code.show_environment()` prints the modules and executable.
  `export OMPI_MCA_btl=^openib` (in .bashrc and yambo_module) only hides a harmless OpenFabrics warning.
- Production RunRules (always `partition='all12h', memory='125000'` and BeeOND):
  - QE: `time='11:59:00', ntasks_per_node=16, cpus_per_task=2, omp_num_threads=2,
    pre_processing='/home/dalessandro/module_script/qe_module'`, `QeCalculator(rr, activate_BeeOND=True)`
  - Yambo: `ntasks_per_node=32, cpus_per_task=1, omp_num_threads=1,
    pre_processing='/home/dalessandro/module_script/yambo_module'`,
    `YamboCalculator(rr, executable='yambo_nl', activate_BeeOND=True)` (or `yambo_rt`, `ypp`, `yambo`)
- Claude submits slurm jobs only when the user asks explicitly, or for repairs/tests agreed with the user (user's rule,
  2026-10-09); otherwise the user runs the notebooks that launch them from Jupyter. Jobs of the production notebook are
  launched without Jupyter with `~/tmp_claude/launch_dataset.py` (see the Handoff section). The analysis notebooks
  (e.g. `NL-Chi_Analysis.ipynb`) Claude may execute directly on the login node:
  `nohup jupyter nbconvert ... --execute --inplace ... > ~/tmp_claude/<nb>.log 2>&1 &`.
- Claude can run either on the laptop (desktop app; cluster commands through
  `ssh -o BatchMode=yes -o ClearAllForwardings=yes -o ConnectTimeout=30 ismhpc`, allowed in `.claude/settings.local.json`)
  or directly on the cluster (Claude Code CLI in the conda env `tools`, inside `tmux new -s lif`, started in `~/work/LiF`,
  Remote Control on): there the commands are local (no ssh prefix) and long monitoring survives the laptop.
- Ask the user before launching p2y/yambo/yambo_nl runs (cost, quota). The user prefers clean runs from scratch to
  reusing old outputs, unless stated otherwise.

## QuantumESPRESSO and Yambo: things learned
- QE 7.0 (`pw.x`). Yambo: the Lumen 2.1.0 fork, rebuilt on 2026-10-01 in `~/Applications/Lumen` (core, nl-project,
  rt-project; executables yambo, ypp, yambo_nl, yambo_rt, p2y in `gpl-gcc_10.2/bin`).
- The nscf used for p2y must use `force_symmorphic=True` (default of MPPI PwInput >= 1.3): with non-symmorphic
  symmetries this yambo silently activates no runlevel (inputs with only setup variables, report named `r-..._ypp`).
- p2y always ends with MPI_ABORT after writing the wavefunctions, but the SAVE is fine.
- MPPI YamboInput: keep `reformat=True` (yambo regenerates the input with the `-V` verbosity of args, dropping the
  other variables). Units are set with `set_array_variables(units=...)` (eV, fs, kWLm2, ...).
- ypp bands: use `BANDS_kpts` (not `BANDS_path` labels) and `INTERP_mode='BOLTZ'` (NN gives step-like bands).
- yambo_nl: the DELTA field is written as a single time step of value E0/dt (so int E dt = E0); SIN/SOFTSIN fields start
  at `FieldX_Tstart`; the transient after the switch on decays with the dephasing (NL_damping).
- yambo field intensities (`FieldX_Int`, input unit `kWLm2`, reports `[kWCMm2]`) are in kW/cm^2: E0 = sqrt(4 pi I/c) au
  with I in au (AU2KWCMm2 = 6.4364086e12), 1e6 -> 1.194e-4 au. The "kW/m^2" written in these notes and in the notebooks
  for the yambo_nl runs are kW/cm^2. Davide's RT pump is 6e9 kW/cm^2 (E0 = 9.25e-3 au, 6000 times our 1e6 pump).
- yambo vector alat of a fcc cell is alat/2 (YamboDftParser `rescale=True` lattice differs by 2 from PwParser).
- SAVEs made with older yambo versions: the setup databases (`ndb.gops`, `ndb.kindx`, `ndb.kpts`) written by yambo 5.2
  give NaN overlaps (`ndb.Overlap` DIP_S, `ndb.dipoles` DIP_iR) and a NaN Berry polarization at t=0 with Lumen 2.1.0
  ("Found NaN in carr. Dynamics stopped"): remove them and redo the setup (`yambo` without arguments in the run folder).
- yambo_nl INVINT integrator (default): Cayley step with the Hamiltonian at the beginning of the step. For the static
  part a transition energy E evolves as E_eff = (2hbar/dt) arctan(E dt/2hbar): dt = 0.05 fs red shifts the 10-20 eV
  spectrum by ~2 eV (17.4 -> 15.4 eV), dt = 0.01 fs by ~0.1 eV; the field acts with a delay dt/2 (phase w dt/2, also
  the delta kick). CRANKNIC (midpoint Hamiltonian, 2x cost) has the same static phase error. Use NLstep = 0.01 fs.
- yambo truncates the job string (-J) at 100 characters: keep the run ids of MPPI datasets within 100 characters.
- `UseDipoles` (fixed dipoles instead of the Berry coupling) gives only the linear response correctly.

## MPPI for the analysis (user side)
- Parsers: `YamboNLDBParser(<run>/ndb.Nonlinear)` reads the yambo_nl database and its fragments (attributes
  `IO_TIME_points`, `Polarization` [one (3,nt) array per run], `E_ext`, `Efield`, `Efield_general`, `n_runs`,
  `NL_damping`, `get_time()`); `YamboOutputParser.from_file(o-file)` reads o- files (also `NL_pol_F1` and the
  `YPP-eps_along_E` spectra, columns `col1`, `col2`, ...).
- Optics: `Linear_Response(time, pol, efield, eta=...)` (delta field: eps = 1 + 4 pi P(w)/E(w));
  `Xn_single_frequency(data, X_order=n)` (one SIN field per run); `Xn_frequency_mixing(data, X_order=(n_probe,
  n_pump))` (probe = field 1 varying with the run, pump = field 2 fixed). Methods `check_harmonic_reliability`,
  `eval_Pw`, `eval_Ew`, `compute_Xn(set_units_of_measure=False)`. chi in atomic units (P = chi E), esu with
  `set_units_of_measure=True`.
- Conventions of MPPI >= 1.3 (changed on 2026-10-02, see the MPPI notes and Boyd ch. 1):
  - negative orders use the conjugate field: the keys (1,-1), (1,-2) changed in sign/phase w.r.t. the results computed
    with MPPI < 1.3 (and differ from YamboPy, which does not conjugate);
  - the chi include the degeneracy factor D (see above); the zero-th order of `Xn_single_frequency` is divided by |E|^2;
  - `Linear_Response` now multiplies the FFT by dt (before: factor 2, eps-1 changes by dt/2, ~3% for the 0.05 fs IO step);
  - dephasing threshold 12/damp in both Xn classes (6/damp gave ~3% errors on the third order keys);
  - all the harmonics with a sizeable amplitude must be in the fit (e.g. the third harmonic of the pump, X_order=(1,3)),
    otherwise they spoil the weak terms; if the estimated time window is longer than the simulation the fit starts at
    the dephasing time with a warning.
  - Every result of the LiF notebooks computed before these changes must be recomputed.
- MPPI >= a7aac3f (2026-10-07, review of the Optics classes driven by LiF, see the MPPI CLAUDE.md):
  - `compute_Xn(..., broadening=None, renormalize=True)` of both Xn classes and `O.Utils.lorentzian_broadening(freqs,
    chi, delta)`: a posteriori Lorentzian broadening (eV) of the keys linear in the probe (key 1 of the sine, keys (1,m)
    of the mixing; the others are returned unchanged). Use it instead of the prototype in the notebook.
  - the harmonic fits are stored per frequency (eval_Pw/eval_Ew/compute_Xn no longer repeat them: LiF 201 frequencies
    13 s -> 2 s, identical results); the fit window of the mixing never starts before the dephasing time (no LiF
    frequency was affected); the time sampling warnings are now one summary line (the old ones on LiF were spurious).
  - `generate_frequencies` drops a key whose frequency coincides with another one (e.g. (1,-3) when w = 6 wP).
  - INVINT phase: the fields act with a delay dt/2 (phase w dt/2 = 0.13 rad at 17 eV with NLstep 0.01 fs), the same for
    the delta and the sine/P&p runs, so it cancels in the internal comparisons. It can be removed by adding dt/2 to the
    `initial_time` of the fields (`efield['initial_time']` for `Linear_Response`, and consistently for the Xn classes);
    to be done for the comparison with the RT runs of Davide (dt = 3 as, negligible delay).
- MPPI >= f0f8add (2026-10-20): `O.Transient_Response` (transient response dchi from the chi, monochromatic pump,
  envelope at the probe time, accumulated envelope; see the handoff of 2026-10-20 above), `O.Utils.time_weighted`,
  `AnharmonicOscillator.delta_probe_response` (exact P_Pp - P_P of the oscillator for a delta probe with any pump).
- MPPI tutorials useful as reference: `Analysis_Optics.ipynb` and `Model_AnharmonicOscillator.ipynb` (classical
  anharmonic oscillator with analytical chi: the Optics module reproduces chi1, chi2, chi3 and the mixing keys to
  1e-6-1e-3), `Tutorial_YamboNLDBParser.ipynb`. `mppi.Models.AnharmonicOscillator` can be used to test any new analysis
  formula (e.g. dchi^neq) against exact results before applying it to LiF.

## Status and open points (2026-10-06)
Results so far (details and comments in `NL-Chi_Analysis.ipynb`):
- chi1 from delta, sine and P&p (1,0) agree (kx8: 2.7% delta vs sine after fixing NLstep and the setup). The third
  order term is `(chi(1,0) - chi1)/(EP**2/4)` (EP = `Efield2[0]['amplitude']`, E0 in au), checked also against the
  pump intensity scaling.
- chi(1,+-1) is NOT zero although LiF is centrosymmetric: |chi(1,+-1)||E_P|/|chi1| ~ 1e-3-1e-2, 10 times the third
  order term, resonant with the probe. Excluded: harmonic fit (X_order (1,3), raw spectrum of P(t)), k sampling of the
  ground state (k/-k closure, E and dipoles), nonlinearity/noise (exact E_p*E_P scaling), time integrator (CRANKNIC).
  Remaining candidate: discretization of the Berry coupling in k. The kx12 test (5 frequencies) was inconclusive
  because chi1 itself is not k converged there.
- k convergence of chi1 (delta runs): with damping/broadening 0.1 eV no grid is converged above the gap (kx8 vs kx16
  19% mean, kx12 vs kx16 11%); with 0.3 eV kx12 is converged within 2% (kx8 6%, artifact peak at 17.5 eV); with
  0.6 eV also kx8 (2.6%).
- A posteriori broadening of chi1 and of the keys (1,m) (linear in the probe): Lorentzian convolution on the probe
  grid = chi(w + i Delta); validated on the oscillator (~3%, limited by the window tails) and on LiF (sine broadened by
  0.2 eV vs delta with eta 0.3 eV: 3%); reliable up to Delta ~ 0.3 eV in a 10-20 eV window. Not the same as a larger
  dynamical damping (the pump-only denominators are not broadened), but analogous to Davide's smoothing (he uses
  dephasing 0 in the dynamics and sm 0.6 eV in post-processing).
### GOAL OF THE CURRENT PHASE (2026-10-07)
Understand why the chi2 terms chi(1,+-1), which must vanish in centrosymmetric LiF, are dominant w.r.t. the chi3 terms
((1,0)-chi1 and (1,+-2)). Every analysis of this phase must be oriented to this question; the comparison with Davide
(transient absorption, eval_dchi_neq) is postponed. Results with damping 0.3 eV (final runs, 2026-10-08, 10-20 eV, all
frequencies): quadratic/cubic ratio (median) 22 on kx8 and 9 on kx12 (vs (1,+-2): 40 and 15); |chi(1,+-1)||E_P|/|chi1|
max 2.8e-3 -> 1.3e-3 and median 3e-4/2e-4 -> 1.5e-4/1.2e-4 from kx8 to kx12 (factor ~2, close to (12/8)^2 = 2.25),
while the cubic terms change little (max |chi(1,0)-chi1|/|chi1| 1.1e-4 -> 8.6e-5, median (1,+-2) unchanged): chi2
decreases with the k grid, the cubic terms do not. The chi(1,+-1) of the two grids have different shapes (mean rel.
diff 40-60%). Single field SHG |chi_2||E_p|/chi1 < 3e-6: per unit field and without the degeneracy factor (chi(1,+-1) = 2 chi2) the
probe SHG chi2 is ~5-7 times smaller than the mixing chi2 (medians, different frequency arguments).

Analyses to implement in `NL-Chi_Analysis.ipynb` (new part) when the single-node runs are available:
a. Revision of the new part on the 1-node runs (same names as before); replace the broadening prototype of the old
   part with `compute_Xn(broadening=...)` of MPPI.
b. 1-node vs 2-node runs (`..._2nodes`, `..._corrupted` folders): are the "good" 2-node runs identical to the 1-node
   ones? (relative differences of all the keys, in particular chi(1,+-1)).
c. Signal to noise map: for each frequency and key, residual of the harmonic fit (rms of P(t) - fitted signal in the
   fit window, from the stored fits of `Xn_frequency_mixing`/`eval_Pw`, see `perform_harmonic_analysis`, results dict
   + B0 + residual) vs the amplitude of the components (1,0), (1,+-1), (1,+-2) and of the difference (1,0)-chi1
   (difference of two runs: estimate its noise from the residuals of both). Flag the frequencies where a key is not
   above the noise.
d. Intensity scan (runs `study_scan` of `YamboNL_Analysis.ipynb`, kx8, damping 0.3 eV, 39 frequencies = every 4th of
   the 155 grid, step wP/4: pump 2.5e5 and 4e6 kW/m^2, probe 4e3 kW/m^2; reference = the 1-node kx8 P&p run at the
   same frequencies, indices 0,4,...,152): exponents a (probe) and b (pump) of |P(key)| for each key and frequency
   (expected a=1 and b=0,1,2 for (1,0),(1,+-1),(1,+-2)); stability of the normalized chi between the intensities
   (criterion: <~1% and exponents within ~0.01); third order (1,0)-chi1 from different intensity pairs. Purpose: check
   that the extracted chi are not affected by noise (low field) or higher orders (high field) and whether a larger
   pump can be used to improve the S/N of the cubic terms.
e. Trend of the chi2/chi3 ratio with the k grid (kx8 vs kx12, possibly kx16): the main hypothesis for the chi2 terms is
   the k discretization of the Berry coupling in the dynamics with two fields.
Done (2026-10-08): b (kx8 1 node = 2 nodes bit by bit) and d (section "Intensity scan (kx8)" of the notebook, with
conclusions). Outcome: chi(1,+-1) is exactly bilinear (a = b = 1.000, normalized chi stable to ~1e-3 over x16 pump and
x4 probe intensity): a genuine second order response of the discretized dynamics, ~Dk^2 (kx8 -> kx12). Its
dominance is a statement on the signals at our pump field, not on the chi: the quadratic component grows as E_P and the
cubic one as E_P^2, so with fixed chi the ratio R ~ 1/E_P (median 41.5, 21.2, 10.6 at 2.5e5, 1e6, 4e6 kW/m^2), with no
breakdown of perturbation theory (that would show up as chi changing with the field). The intrinsic quantity is the
crossover intensity: cubic signal larger above I_eq = I_P R^2 = 4.7e8 (kx8), 7.3e7 (kx12) kW/m^2 (median 10-20 eV; RT pump 1e8). Third
order: (1,+-2) stable (b = 2.00); (1,0)-chi1 limited by a floor ~1e-6 chi1 on (1,0) (20% error at pump/4, 1% between
1e6 and 4e6). Next runs: pump >~4e6 kW/m^2 (still perturbative) and kx12/kx16; chi(1,+-1) to be dropped in dchi^neq
(zero by symmetry) but it measures the k discretization error (the cubic keys differ 35-50% between kx8 and kx12).
Still open: c (S/N map from the fit residuals, partly superseded by the scan) and e on kx16.
Broadening 0.5 eV (2026-10-08: 0.3 eV damping + `compute_Xn(broadening=0.2)`, shown side by side with 0.3 eV in the
notebook): chi1 and the cubic keys get closer between kx8 and kx12 (chi1 4.8% -> 2.3% mean; chi(1,1,-1) 48% -> 35%,
ratio of the maxima 1.67 -> 1.24; delta kx12 vs kx16 0.8%), while chi(1,+-1) does not (ratio of the maxima kx8/kx12
3.5-4.2 -> 3.3-3.6, mean difference 60-90%) and R does not decrease (median 22 -> 29 kx8, 8.5 -> 8.9 kx12; I_eq
8.6e8 / 7.8e7 kW/m^2). So chi(1,+-1) is not a slowly converging spectral sum but a k discretization error of the
coupling (Berry finite differences), R ~ Dk^(2-3) / E_P: only a finer k grid and a stronger pump reduce it.

Next:
1. Run the damping 0.3 eV sine and P&p runs (kx8 and kx12, 10-25 eV, 155 frequencies with step wP/16, NLtime 100 fs)
   and revise the whole `NL-Chi_Analysis.ipynb` on them (chi extraction, chi(1,+-1) relative size, third order).
2. `eval_dchi_neq(omega, omega_P, tau, EP, x11, x1m1, x11m1, x12, x1m2)`: check it against the formula above (shifted
   grids, now exact with the step wP/16; complex pump amplitude i*EP/2 instead of the real EP; phases exp(+-i wP t0P)
   and the delay convention; E_p(w')/E_p(w) ratios of the probe used in the RT runs). Test it on the oscillator.
3. Done in MPPI (a7aac3f): use `compute_Xn(broadening=...)` in the notebook instead of the prototype. Before the
   comparison with Davide, remove the INVINT phase w dt/2 consistently from all the chi (see the MPPI section).
5. The damping 0.3 eV runs are being repeated on a single node (2026-10-07): the 2-node runs showed MXM segmentation
   faults (hundreds of 400 MB core dumps, deleted) and the kx12 sine run had 36 wrong frequencies. The 2-node runs
   are kept with the suffixes `_2nodes` / `_corrupted` for comparison. Use 1 node for yambo_nl until the inter-node
   communication is understood.
   Outcome (2026-10-08): the kx8 1-node runs are bitwise identical to the 2-node ones; the kx12 1-node sine is good,
   but the kx12 1-node P&p (wnode07) is corrupted at 37 frequencies (indices 10,14,...,154, |P| 76-89% of the right
   value, no error message), while the 2-node kx12 P&p (wnode02,05) is good at all the frequencies (identical to the
   1-node one elsewhere, (1,0) consistent with the sine and the delta). So the kx12 P&p used is the 2-node run, renamed
   to the original name; the 1-node one is kept as `..._corrupted`. Both corrupted runs ran on wnode07 (the kx8 2-node
   P&p on wnode02,07 is fine): suspected faulty node, intermittent and silent. wnode07 is excluded from all the
   yambo_nl jobs (`exclude_nodes = 'wnode07'` in the RunRules cell of `YamboNL_Analysis.ipynb`, passed as
   `RunRules(..., exclude=exclude_nodes)` to rr, rr_debug and rr_2nodes; option of MPPI >= 1e65f4c) and a new run is
   always checked against an independent one (sine vs (1,0), delta).
   UPDATE (2026-10-08): the kx12 pump test run with 1.6e7 kW/m^2 (39 frequencies, wnode02, NL_CPU 4.8) has the same
   corruption: indices 10,14,...,38 (|P| ~80%), i.e. the same frequency group (indices 2 mod 4) from its third
   frequency, as the corrupted kx12 1-node P&p (indices 10,14,...,154); the 2-node kx12 sine had two groups (1,2 mod 8)
   from their third frequency. So it is NOT a faulty node but an intermittent yambo_nl problem, so far only on kx12,
   always starting at the third frequency of a group. The wnode07 exclusion is harmless but not the fix: check every run
   (|chi(1,0)-chi1|/|chi1| > 1e-2 flags the corrupted frequencies, or the pattern of max|P| vs the neighbours).
   2026-10-09: the kx16 sine (groups 1 and 2 from indices 9 and 10) and P&p (group 1 from 9, group 2 from 50) runs are
   also corrupted; the delta run of the same grid is the arbiter (exact below the gap). Restart tests on the kx12 1.6e7
   run (section "Restart of the corrupted frequencies" of `YamboNL_Analysis.ipynb`): recomputing the 8 corrupted
   frequencies with the same NL_CPU 4.8 (same group) gives correct values, and NL_CPU 1.32 gives fragments identical bit
   by bit: the error is random/intermittent, correct runs are exactly reproducible. NL_CPU 1.32 (32 tasks on k) is ~2.5x
   slower than 4.8. Repair: yambo_nl recomputes only the frequencies whose `ndb.Nonlinear_fragment_<index+1>` is missing
   (`NL_start_and_restart`, "RESTART, remaining frequencies"): rename the corrupted fragments, back up job out/report/LOG,
   rerun the same input with `skip=False, clean_restart=False, activate_BeeOND=False` (with BeeOND MPPI copies the job
   folder content in the scratch root while yambo looks in <scratch>/<jobname>, so the fragments would not be found).
   Helpers `backup_run_files`, `remove_fragments`, `missing_fragments` in the notebook; launcher script
   `~/tmp_claude/launch_dataset.py <dataset> [--prepare "<first line of a '# restart preparation' cell>"] [--dry]`.
   kx16 repair launched 2026-10-09 (jobs 59071/59072, slownodes, NL_CPU 4.8). A report for the yambo developers is to
   be written in a dedicated notebook.
4. Decide whether to move the old MPPI Analysis_Optics (LiF version) here.

## Outcome of the handoff (2026-10-10, CLI session on the cluster)
- kx16 repairs done: P&p 64/64 (job 59072), sine 70/70 in chunks (index 154 flagged at 0.976 only because the health
  reference at the grid edge is the constant value of index 152: recomputed twice bit-identical, restored by hand,
  copies `_good_chunk*`). Full check: kx16 sine vs delta mean 1.4% (as kx12), P&p (1,0) vs sine max 4e-5, no
  corrupted frequency. Steps 1-3 of the handoff below are done; step 4 (commit) is pending (user's go).
- `NL-Chi_Analysis.ipynb` executed with kx16 and the pump tests (conclusions cell rewritten). Results (10-20 eV,
  medians): chi(1,+-1) decrease kx12 -> kx16 by 2.6-3.6 (|chi| ~ dk^p, p = 3.4-4.5; kx8 -> kx12 p = 1.0-2.9), single
  field SHG too: k discretization error, drop them. R = 21.7/8.5/5.2 (0.3 eV) on kx8/kx12/kx16, I_eq kx16 2-2.7e7
  kW/m^2 (cubic ~2x quadratic at the RT pump 1e8). Cubic keys: size stable (p ~ 0) but NOT converged in shape in
  16-21 eV (kx12 vs kx16 29-46% of the max, maxima ratio 1.4-2.1; chi1 converged at 1-2%). Pump 4e6 perturbative,
  1.6e7 not (fifth order, changes x5). Open: k convergence of the cubic terms (kx20/kx24 or 0.6 eV total broadening).
- `NL-Chi_Analysis.ipynb` simplified (2026-10-12, user's request, script `~/tmp_claude/simplify_nlchi.py`): no damping
  0.1 eV in the new part, kx8 only in the delta section (eta 0.3/0.5, shows kx8 is not adequate), non-linear terms
  plotted with 0.5 eV on kx12/kx16 (table with p from the kx12/kx16 medians), the kx8 intensity scan and the "dependence
  on the k grid" section removed, one kx12 pump test section (exponent b 1e6->4e6 on 155 frequencies, changes at 1.6e7).
  Then the "Previous analysis" part was removed and the transient absorption part moved to `Dchi-neq_Analysis.ipynb`.
- Cleanup of the runs (2026-10-12, user's request, `~/tmp_claude/cleanup_runs.sh`): see the list in the NL_Chi data
  section; `YamboNL_Analysis.ipynb` (sections of the deleted runs removed, script `~/tmp_claude/clean_yambonl_nb.py`)
  and `YamboNL_Corrupted_Frequencies.ipynb` (case A dropped from the computation, kept as text) updated and checked.
- `YamboNL_Corrupted_Frequencies.ipynb` executed on the final data (cases A-F; 9 of 10 corrupted groups start at the
  3rd frequency of the group; a corrupted group never recovers).

## Handoff (2026-10-09, from the laptop session to the Claude CLI session on the cluster)
The laptop session ("NL-CHI") hands the work over to a Claude CLI session on the cluster (`tmux` session `lif`,
started in `~/work/LiF`). State at 18:20 of 2026-10-09:

Running (independent of any Claude session; check with `squeue -u $USER` and the logs):
- **kx16 P&p repair**, job 59072 (slownodes, wnode01, NL_CPU 4.8, no BeeOND): recomputes 64 corrupted frequencies
  (group 1 indices 9..153 and group 2 indices 50..154) of
  `kx16_nb100_NoTr_E100/Pp-Ep_1e3-EP_1e6-PE_1.55-nlenrange_10.0-25.016-nlensteps_155-bands_3-6-damp_0.3-sin-nltime_100`.
  22/64 done at 18:20, all correct; ~5 h left. Its launcher process was killed (only the slurm job runs): nothing
  to do at the end but the checks below. Health: `python ~/tmp_claude/check_repair_progress.py Pp` (new/ref of max|P|
  w.r.t. the good groups 0 and 3, `*` if off by more than 2%).
- **kx16 sine chunk repair**: driver on the login node (`launch_dataset.py none --exec "# chunk repair: kx16 sine"`,
  log `~/tmp_claude/chunk_repair_kx16_sine.log`, one line per chunk with the health ratio of each frequency) that
  recomputes the 70 remaining corrupted frequencies of
  `kx16_nb100_NoTr_E100/pulse-E1_1e3-nlenrange_10.0-25.016-nlensteps_155-bands_3-6-damp_0.3-sin-nltime_100` in 18
  chunks of 2 frequencies per group (jobs on all12h, ~30 min each + queue). Chunks 0 and 1 (17,18,21,22,25,26,29,30)
  correct. Frequencies 9, 10, 13, 14 were already repaired by job 59071 (cancelled when its groups got corrupted again
  from their 3rd frequency). Fragments: `_corrupted` (original run), `_corrupted_repair1` (job 59071),
  `_corrupted_chunk<k>_pass<p>` (chunks off tolerance, redone in a later pass). The driver ends with "all the 70
  frequencies done" or "stop after 3 passes". If it has to be stopped: `scancel` its job, kill the two python
  processes (MPPI runs the task in a child process, SIGINT does not work), then `python ~/tmp_claude/clean_placeholders.py`
  (removes the placeholder fragments, identical to their `_corrupted` copy) before starting it again.

Helper scripts (in `~/tmp_claude`, run with `~/miniconda3/bin/python`): `launch_dataset.py <dataset|none> [--prepare
"<first line of a '# restart preparation' cell>"] [--exec "<first line of a '# chunk repair' cell>"] [--dry]` (executes
the first cells and the damping 0.3 eV section of `YamboNL_Analysis.ipynb`, then runs one dataset; use `nohup ... &`);
`check_repair_progress.py [sine|Pp]`; `clean_placeholders.py`; `check_kx16_sine.py` (sine vs P&p vs delta kx16 per
frequency group); `check_restart.py`; `plot_corrupted.py` (figure good vs corrupted chi1 of the kx12 P&p run, used for
the report). Backups of job outputs/reports/LOG of the restarts: `<run>_backup_<tag>` folders next to the runs.

Next steps, in this order:
1. When both repairs end: full check of the kx16 runs (sine first harmonic vs P&p (1,0) at every frequency, < 1e-3
   on good frequencies; both vs the kx16 delta, exact below the gap; `check_kx16_sine.py`); redo any frequency still
   wrong with the same restart mechanism. Then remove the placeholder logic from the plan and note the outcome here.
2. `NL-Chi_Analysis.ipynb`: add kx16 to the damping 0.3/0.5 eV section (chi1, the keys (1,+-1), (1,+-2), chi(1,1,-1)
   on kx8/kx12/kx16): convergence kx12 -> kx16 of the cubic terms (the key question: are the chi3 converged at kx12/kx16?),
   power law in Dk of chi(1,+-1) with three grids, quadratic/cubic ratio and I_eq; add the pump test (4e6 and 1.6e7 on
   kx12, results above: 4e6 perturbative, 1.6e7 not) and the kx12 P&p run with pump 4e6 (155 frequencies, good:
   third order = 1e6 result within 0.2%). Keep the analysis concise, oriented to the goal of the phase.
   Prepared (2026-10-09, tested on a copy with kx8/kx12): `~/tmp_claude/add_kx16_analysis.py` (backs up the notebook,
   adds kx16 to the damping 0.3/0.5 eV section, a "Dependence of the second and third order terms on the k grid" cell
   with the exponent p of |chi| ~ dk^p, the sections "Pump intensity test (kx12)" and "kx12 P&p run with pump 4e6").
   To apply when the kx16 repairs end, then execute the notebook and rewrite the conclusions cell on the numbers.
3. Dedicated notebook for the yambo developers on the corrupted frequencies: list of the cases (grid, nodes, NL_CPU,
   group, first corrupted frequency = 3rd computed by the group in all cases but one (13th), size of the error),
   the figure of `plot_corrupted.py` (P(t) wrong from the first time step, the whole response reduced by 30-50%, as
   if part of the k points were missing in the group), the restart tests 1a/1b (correct runs reproducible bit by bit,
   error intermittent), the chunk repair (2 frequencies per group never corrupted), the MXM segfaults of the 2-node
   runs, the input and the yambo version (Lumen 2.1.0).
   Draft done (2026-10-09): `NL_Chi/YamboNL_Corrupted_Frequencies.ipynb`, generated by
   `~/tmp_claude/devnb/make_dev_notebook.py` (cases A-F, fragment-by-fragment comparison wrong vs correct). Findings:
   P(t) identical bit by bit up to the switch-on of the field, wrong from the first step with the field (step 12);
   projection on the correct P(t) constant within a group (e.g. 0.564 on 37 frequencies of the kx12 P&p); first wrong
   frequency = 3rd of its group in 9 of 10 corrupted groups. To be re-executed (and the summary revised on the final
   numbers) when the kx16 repairs end: before that the kx16 cases D/E compare with placeholders.
4. Commit the notebooks (`YamboNL_Analysis.ipynb`, `NL-Chi_Analysis.ipynb`) at the end of the phase (never data).
