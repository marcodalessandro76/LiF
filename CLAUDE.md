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
  (Davide's reports and tarballs in `NL_Chi/Davide_Google_drive`, ~320 MB, and
  `RT_Transient_Absorption/Report_LiF_Davide.pdf`).
  Never `git add` data folders, `.ipynb_checkpoints` or tarballs (see `.gitignore`).
- `NL_Chi/` (current work):
  - `Transient_Abs_NL-Chi.ipynb`: yambo_nl datasets (built with MPPI YamboInput/YamboCalculator/Dataset, slurm),
    extraction of the chi with `Xn_single_frequency` / `Xn_frequency_mixing`, `eval_dchi_neq` and `build_delta_dict`,
    comparison with Davide's spectra for each delay. 89 cells, ~14 empty; written with MPPI < 1.3.
  - cluster only, `kx8_nb100_NoTr_E100/`: a copy of Sangalli's SAVE (`/work/sangalli/simulations/LiF/yambo/
    e80_kx6_pbe_sr/kx8_nb100/NoTr_E100`, kx8 grid, 100 bands, no time reversal, field along 100) and the yambo_nl runs
    (bands 3-6, damping 0.1-0.2 eV, SIN fields, Tstart 0.1 fs): `lresponse-bands_3-6-delta` (linear response, delta),
    `pulse-*` (probe only, 5-25 eV/100 steps and 10-20 eV/201 steps, intensity 1e3 kW/m^2), `Pp-Ep_1e3-EP_1e6-PE_1.55-*`
    (pump 1.55 eV at 1e6 + probe at 1e3, same two probe grids), `pump_1.55-*` (older tests), slurm `job_*.sh/.out`.
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
- Never run heavy computations on the login node `frontend`: use slurm. Partition `debug` (2h) for short tests,
  `all12h` for production (32 cores per node). Home quota 19.5 GB (11 GB used on 2026-10-05: watch the size of new runs,
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
- MPPI tutorials useful as reference: `Analysis_Optics.ipynb` and `Model_AnharmonicOscillator.ipynb` (classical
  anharmonic oscillator with analytical chi: the Optics module reproduces chi1, chi2, chi3 and the mixing keys to
  1e-6-1e-3), `Tutorial_YamboNLDBParser.ipynb`. `mppi.Models.AnharmonicOscillator` can be used to test any new analysis
  formula (e.g. dchi^neq) against exact results before applying it to LiF.

## Open points (2026-10-05)
1. In `Transient_Abs_NL-Chi.ipynb` the third order term is computed as `x11m1 = (chi_freqmix[(1,0)] - chi_lr[1])/Ew[(0,2)]`:
   wrong, Ew[(0,2)] = E_P(wP)^2; it must be divided by |E_P(wP)|^2 = EP**2/4 (EP = pump amplitude in au).
2. `eval_dchi_neq(omega, omega_P, tau, EP, x11, x1m1, x11m1, x12, x1m2)`: check it against the formula above (shifted
   grids, complex pump amplitude i*EP/2 instead of the real EP, phases exp(+-i wP t0P) and the delay convention,
   E_p(w')/E_p(w) ratios of the probe used in the RT runs; the docstring mentions x10 but the argument is x11m1;
   `round(omega_P/domega)` assumes wP is a multiple of the grid step). Test it on the anharmonic oscillator.
3. LiF (rock salt, Fm-3m) is centrosymmetric: all even order chi vanish, so chi(1,+-1) and the F2 term (linear in E_P)
   of `eval_dchi_neq` must be zero within the numerical accuracy; check |chi(1,+-1) E_P| << |chi(1,+-2) E_P^2| on the
   yambo_nl data (on the oscillator the even keys are ~1e-9 of the linear one). The leading pump induced terms are
   then chi(1,0)-chi1 and chi(1,+-2).
4. Recompute all the chi with MPPI >= 1.3 and check the old warnings ("time sampling starts before the dephasing time":
   some probe frequencies need a long time window because of near-degenerate harmonics).
5. Clean up `Transient_Abs_NL-Chi.ipynb` (empty cells; the old "Field 3 not found" prints disappear with MPPI 1.3) and
   decide whether to move the old MPPI Analysis_Optics (LiF version) here.
