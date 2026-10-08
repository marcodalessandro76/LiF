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
  machine, pull on the other). In git there are only README, this file, `Attivita nuova.txt`, the notebooks and `References/`:
  the data (~13 GB on the cluster) are NOT in git and live only on the cluster (the yambo runs) or on the laptop
  (Davide's reports and tarballs in `NL_Chi/Davide_Google_drive`, ~320 MB).
  Never `git add` data folders, `.ipynb_checkpoints` or tarballs (see `.gitignore`).
- `NL_Chi/` (current work):
  - `YamboNL_Analysis.ipynb`: the yambo_nl datasets (MPPI YamboInput/YamboCalculator/Dataset, slurm, all appended
    with `skip=True`), no analysis: inversion check of the k sampling (`YamboDftParser.expand_IBZ_kpoints`), kx8
    production runs (delta, sine, P&p, damping 0.1 eV), diagnostic runs for the chi(1,+-1) terms (5 probe frequencies
    17.41-17.61 eV: reference, pump/4, probe*4, NLstep/2, CRANKNIC, kx12 reference), delta runs on kx12/kx16, and the
    damping 0.3 eV sine + P&p runs on kx8 and kx12 (single dataset `study_eta`, 4 runs in series, 2 nodes each).
  - `NL-Chi_Analysis.ipynb`: chi extraction (linear in three ways, second order), diagnostics of chi(1,+-1),
    k grid convergence of chi1 from the delta runs, a posteriori Lorentzian broadening on the probe grid
    (`lorentzian_broadening`, validated on the anharmonic oscillator), then the not yet revised part (third order,
    `eval_dchi_neq`, `build_delta_dict`, comparison with Davide). To be revised on the damping 0.3 eV runs.
    Both notebooks come from the split (2026-10-05) of the old `Transient_Abs_NL-Chi.ipynb` (in the git history).
  - cluster only, `kx8_nb100_NoTr_E100/`, `kx12_nb100_NoTr_E100/`, `kx16_nb100_NoTr_E100/`: copies of Sangalli's
    NoTr_E100 SAVEs (`/work/sangalli/simulations/LiF/yambo/e80_kx6_pbe_sr/kx{8,12,16}_nb100/NoTr_E100/SAVE`: same
    lattice, 100 bands, the same 8 symmetries without inversion and time reversal, 100/294/648 IBZ points), copied
    without the yambo 5.2 setup databases and with the setup redone (`yambo` without arguments) with Lumen 2.1.0.
    kx8: production runs (`lresponse-bands_3-6-delta` with NLstep 0.01 fs, `pulse-E1_1e3-nlenrange_10.0-20.0-
    nlensteps_201-bands_3-6-damp_0.1-sin`, `Pp-Ep_1e3-EP_1e6-PE_1.55-nlenrange_10.0-20.0-nlensteps_201-bands_3-6-
    damp_0.1-sin`, the last two with Lumen 2.0.0) and the diagnostic runs `...nlensteps_5...`; kx12: delta,
    diagnostic reference; kx16: delta. The damping 0.3 eV runs (`...nlenrange_10.0-25.016-nlensteps_155-...-damp_0.3-
    sin-nltime_100`) go in the kx8 and kx12 folders.
  - cluster only, `Davide_data/`: Davide's RT results at kx16 (`o-abs_0-50eV_sm0.6eV.YPP-eps_along_E` for each delay,
    in folders named `...delta<value>fs...`, read by `build_delta_dict`), plus tarballs.
- `RT_Transient_Absorption/` (earlier work, transient absorption from RT simulations): `Transient_Absorption.ipynb`
  (fit of the experimental IR pump, MPPI Gaussian pulses vs ypp_rt fields, yambo_nl RT pump-probe dynamics with a 1.56 eV pump of fwhm 7.58 fs and a
  35 eV probe of fwhm 0.42 fs, transient spectra and reflectivity, delays 7.0 and 9.12 fs) and `Sum_frequency.ipynb`.
  Cluster only: `NoTr_E100` (SAVE + collisions `coll-par-cv_1-15B_X59RL-50B_noSM_H339RL_kpar`, fixsym for the 100
  field), `Field_test`, `nl_input_template`, three `td-hsex-*` runs.
- The old version of MPPI's `Analysis_Optics.ipynb`, based on the LiF runs (now replaced in MPPI by an analytical model),
  is in the MPPI git history: `git -C <MPPI> show b785feb:sphinx_source/tutorials/Analysis_Optics.ipynb`. It can be
  moved here.

## Cluster ismhpc
- Reached from the laptop with `ssh -o BatchMode=yes -o ClearAllForwardings=yes ismhpc '...'` (host in
  `~/.ssh/config`, jump through `narro`; ClearAllForwardings avoids the LocalForward 4444 clash). CentOS 7 /
  glibc 2.17: VS Code Remote-SSH does not work. The connection sometimes times out: just retry later.
- Never run heavy computations on the login node `frontend`: use slurm. Always use the partition `all12h` (32 cores
  per node), also for the short test/diagnostic runs: it runs on all the nodes and gives no problems (user's choice,
  2026-10-06; the `debug` partition is not used). Home quota 19.5 GB (11 GB used on 2026-10-05: watch the size of new runs,
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
- Claude does not submit sbatch/slurm jobs (the notebooks that launch them are run by the user from Jupyter). The
  analysis notebooks (e.g. `NL-Chi_Analysis.ipynb`) Claude may execute directly on the login node, as the user does
  through the Jupyter tunnel: `nohup jupyter nbconvert ... --execute --inplace ... > ~/tmp_claude/<nb>.log 2>&1 &`.
  Light ssh commands (git pull, ls, reading outputs) use the prefix
  `ssh -o BatchMode=yes -o ClearAllForwardings=yes -o ConnectTimeout=30 ismhpc` (allowed in `.claude/settings.local.json`).
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
4. Decide whether to move the old MPPI Analysis_Optics (LiF version) here.
