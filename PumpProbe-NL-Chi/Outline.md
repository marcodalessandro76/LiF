# Understanding pump and probe signal from excitonic non linear response functions

Working outline of the article (provisional title). Started on 2026-11-23 from the user's draft (kept verbatim at the
end of this file) and the discussion with Claude. Target: Phys. Rev. B (revtex4-2, can be changed later).

References to the technical notes (`Notes/main.tex`) are given with the numbers of 2026-11-23 and with the LaTeX
labels, which do not change when the notes are edited.

## Main message

The RT pump-probe simulations give the transient response Delta eps(w, tau) exactly (for a given Hamiltonian) but they
do not explain it. Built from the equilibrium non-linear susceptibilities, the same signal is a sum of terms with a
fixed dependence on the pump field, each one tied to the poles of chi3: the features of the signal can be traced back to
properties of the system. Shown on the exciton of LiF (HARTREE+SEX dynamics), IR pump at 1.55 eV, delta (attosecond)
probe.

## Structure

### 1. Introduction
Written by D. Sangalli. It contains the motivation of the approach (decided 2026-11-24; a short placeholder in
`main.tex`): the established approach is the RT simulation with explicit pump and probe, exact but
with no direct physical interpretation; the chi are properties of the system and the field enters in a way fixed by
the formalism, so building the same response from the chi explains the origin of its features.

### 2. Theory: RT vs chi approach to dchi^neq (provisional title)
Sec. 2 holds the theoretical frame of the article only (the motivation is in the Introduction).
- Rule (agreed): Sec. 2 is general (no LiF, no codes); everything specific to LiF and to the computations goes in Sec. 3.
- Subsections:
  - 2.1 Pump-probe response (DRAFTED 2026-11-24): conventions, chi^neq = (P_Pp - P_P)/E_p (Eq. 7 of the notes,
    `eq:chineq`) and how the RT computes it, delta probe, Delta eps.
  - 2.2 Expansion in the equilibrium susceptibilities (DRAFTED 2026-11-24): P^(n), dchi^(2), dchi^(3).
    Removed from 2.2 (user, 2026-11-24), to be introduced where they are needed: (a) the term with n pump fields is
    multiplied by (-1)^n when E_P -> -E_P, so the odd and even terms are separated with the two signs of the pump; with
    inversion symmetry dchi^(2) = 0 and the response is even in E_P (needed in Sec. 3, even part of the RT); (b) if
    chi3 depends only on w, the convolution theorem gives dchi^(3) = 3 chi3(w) E_P(t_p)^2, the response to the pump
    field at the probe time (possible link to the envelope at the probe time in 2.3).
  - 2.3 Pulsed pump (DRAFTED 2026-11-24, restructured on the user's request): general formulas for the Gaussian pump,
    no approximation on the envelope and no mention of the envelope at the probe time; no comment on inversion
    symmetry. Revised 2026-11-26: opening sentence (Eqs. 5-6 specialized to a pulsed pump, ranges of Omega that
    contribute, evaluated in II D); pump and its +-
    frequency terms (time and spectrum); second order (exact, sum over E^+-); third order (the four products E^+E^-,
    E^-E^+ around Omega = 0, E^+E^+ around 2wP, E^-E^- around -2wP, no a,b indices; cross terms equal by
    w1 -> Omega - w1 and permutation symmetry; chi3 at w1 = +-wP; S_0, S_+-2). Second order exact; third order with the only approximation of a pump far from resonance
    (pump arguments of chi3 fixed at +-wP). The
    monochromatic chi2_+-, chi_0, chi_+- are NOT defined here (user, 2026-11-26: useless at this level).
  - 2.4 Approximations for the pump envelope (structure of 2026-11-26, user's choice: the two approximations are
    subsubsections of one subsection, II D 1 and II D 2). Short opening: for a pulse the chi enter also away from the
    centers, the Omega dependence must be treated, two distinct approximations follow.
  - 2.4.1 Envelope at the probe time (DRAFTED 2026-11-24): the simple case, which connects to expressions already in
    the literature (to be cited: TODO). Opens with the chi taken out at the centers, where they are those of a
    monochromatic pump: definitions of chi2_+-, chi_0, chi_+- (`eq:chi_mono`, moved here 2026-11-26; reused in 2.4.2).
    Second and third order results by the convolution theorem, monochromatic limit, validity (eta >> 1/sigma; in time
    1/eta << sigma; under/overestimate before/after the maximum).
  - 2.4.2 Accumulated envelope (to do): the relevant case; salient points and result in the text, detailed derivation
    in Appendix A.
- Reduced version of Secs. 1.1-1.3 of the notes: conventions, delta probe, chi^neq and its expansion in the
  susceptibilities (Eqs. 1-11).
- Sec. 1.4 of the notes (monochromatic pump) is not a separate part: from the general third order expression (Eq. 9,
  `eq:dchi3`) directly to the pulsed pump. Its result stays in the text: the structure
  |E_P|^2 [chi_0 - chi_+ e^{-2i wP tau} - chi_- e^{2i wP tau}] and the definitions of chi_0, chi_+- (Eqs. 16-17,
  `eq:chi3_mono`, `eq:dchi_mono`) appear as the limit of the envelope at the probe time (Eq. 26, `eq:dchi_probetime`).
  They are the basis of the DC/AC analysis of Sec. 4. Second order terms: one sentence (zero in a centrosymmetric
  system). The meaning of the DC (chi_0) and AC (chi_+-) parts is NOT discussed here but in Sec. 3, with the results for
  LiF (user, 2026-11-24).
- Envelope at the probe time: short, in the text (standard approximation; its failure is the reason for the accumulated
  envelope). Derivation in Appendix A. To be revisited while writing.
- Accumulated envelope (Sec. 1.5 of the notes, up to Eq. 32, `eq:dchi_accumulated`): probably the main theoretical
  development of the work. Salient points in the text, details of the model in Appendix A.

### 3. Pump and probe signal at the exciton peak in LiF
HSEX level only (the IP level is excluded: at most one sentence or the Supplemental Material).
- Equilibrium: the exciton, 1 + 4 pi chi1 vs the RT probe-only run (as Fig. 9 of the notes, `fig:eq_hsex`).
- Meaning of the DC part (chi_0: follows the envelope, cycle average) and of the AC part (chi_+-: oscillations at 2wP,
  seen only with a probe shorter than the pump cycle), with the results for LiF (moved here from Sec. 2).
- Computational details of the two approaches, without heavy details: DFT (functional, lattice, 100 bands, scissor
  6.0 eV + stretch 0.05), kx16 (the 8 symmetries without inversion/time reversal), bands 3-6, HSEX kernel (collisions),
  codes and versions (named here, unlike the notes); RT: CRANKNIC 3 as, 40 fs, delta probe, QSFIELD pump; chi: sine
  and pump-probe runs, harmonic analysis (Sec. 2 of the notes, short with references, details in Appendix B),
  CRANKNIC correction (energy map).
- Spurious second order: the keys chi_(1,+-1) are not zero and the RT response has an odd part in E_P. Same origin
  (k discretization of the Berry coupling), one paragraph: removed by dropping the keys and by taking the even part of
  the RT (reversed pump).
- k convergence: PENDING (user's decision). The tests (k grid, chi2, pump intensity) are at the IP level only; at the
  HSEX level only kx16 and bands 3-6 (and the bands above 6 matter: 36% on Im eps at equilibrium, 3-6 vs 3-12). For the
  comparison of the two approaches it is enough that grid and bands are the same in both.
- Maximum pump intensity compatible with the chi up to third order (pump test of the chi; RT 6e9 vs 6e7 x 100).
- Reference figure: 2D map Delta eps(w, tau), RT vs chi, as Fig. 13 of the notes (`fig:maps_hsex`) but with the same
  approximations in both: bands 3-6, HSEX, kx16, pump 6e7 kW/cm^2 (third order regime), even part of the RT, accumulated
  envelope for the chi. 41 delays: t_p = 11.00-21.00 fs with step 0.25 fs, tau = -4.91..+5.09 fs (the 0.25 fs grid of
  the RT runs of D. Sangalli). New RT runs: handoff of 2026-11-23 in `CLAUDE.md`.
- Line cuts at a few delays, envelope at the probe time vs accumulated (as Fig. 11 of the notes, `fig:dchi_rt_hsex`):
  it justifies the model of Sec. 2.

### 4. Analysis and discussion of the non-equilibrium signal
Once the chi signal is known, its terms are analysed to understand their physical origin.
- Literature: "model + delta pump" (to be clarified by the user).
- Proposals of Claude (to be discussed):
  - DC part first (98% of the map; it survives with a probe longer than the pump cycle), then the AC part (it needs a
    probe shorter than the pump cycle).
  - DC part: change of the line shape of the exciton near its equilibrium energy. From the fit of the notes it is mainly
    a broadening (symmetric shape +, -, +) with a small blue shift: avoid calling it a "derivative feature".
    - RT: phenomenological interpretation as a change of the exciton parameters dE, dGamma, dS. In the notes the fit is
      done on chi_0 (Sec. 3.5.3, `sec:lif_delay_bands`); here do it on the DC part of the RT signal (what one can do
      with the RT, or an experiment, alone).
    - chi: microscopic reason, the poles of chi3. Why does it happen? Idea to be verified (new work, not in the notes):
      chi3(w; w, wP, -wP) has terms with two denominators at the same frequency w, a double pole at the exciton that
      gives the shape dchi1/dw; the coefficient is a sum over the states reached with one pump photon (E_x +- hbar wP):
      its real part gives the Stark shift, its imaginary part the broadening when those states are real.
  - AC part: oscillations at 2wP (Sec. 3.5.2 of the notes, `sec:lif_delay_2wp`): poles of chi_+-, replica of the
    exciton at E_x + 2 hbar wP, interference picture.
  - Message of the section: each feature of the map corresponds to a pole of chi3.

### 5. Conclusions
With the outlook: Delta R/R (the experimental observable), intensities beyond the third order, ...

### Appendices
- A. Envelope of the pump: derivation of the envelope at the probe time and of the accumulated envelope.
- B. Non-linear susceptibilities from the real-time polarization (Sec. 2 of the notes).
- (C. Correction of the integrator, if it does not fit in Sec. 3.)

## Figures (provisional)
| # | Content | Source |
|---|---------|--------|
| 1 | Equilibrium spectrum (HSEX): RT probe-only vs 1 + 4 pi chi1 | notes Fig. 9, `eq_hsex` |
| 2 | Map Delta eps(w, tau), RT even part vs chi, 6e7, 41 delays | new, needs the new RT runs |
| 3 | Line cuts at a few delays, envelope at the probe time vs accumulated | notes Fig. 11, `dchi_rt_hsex` |
| 4 | Im chi_0, chi_+, chi_- (poles) | notes Fig. 10, `chi3_hsex` |
| 5 | DC part: exciton line fit | notes Figs. 16-17, `dc_hsex`, `line_fit_hsex` |
| 6 | AC part: 2wP amplitude and phase, traces vs the delay | notes Figs. 14-15, `ac_hsex`, `ac_traces_hsex` |
| ? | Setup: pump envelope and probe times | notes Fig. 1, `rt_setup` |

Figures 4-6 will be redone on the new runs where possible (pump 6e7, even part, bands 3-6). Figure code:
`make_figures.py` of this folder, with the data and the calculations taken from `Notes/make_figures.py`.

## Open points and decisions pending
- k (and bands) convergence at the HSEX level: pending (user).
- "literature: model + delta pump" in Sec. 4: to be clarified.
- Envelope at the probe time: in the text for now, to be revisited while writing.
- Authors: M. D'Alessandro, C. Attaccalite, M. Grüning, D. Sangalli (order of 2026-11-26); affiliations still placeholders in
  `main.tex`.
- Material listed in the "Reminder for the article" of `CLAUDE.md` (computational details, codes and versions, relation
  with the simulations of D. Sangalli, validation of the extraction, a posteriori broadening, corrupted frequencies of
  yambo_nl, outlook): to be placed.
- RT runs of D. Sangalli (6e9, bands 3-12): at most for the comment on the pump intensity.

## Original draft (user, 2026-11-23)

```text
TITOLO : UNDERSTANDING PUMP AND PROBE SIGNAL FROM EXCITONIC NON LINEAR RESPONSE FUNCTIONS

1)Introduction

......la scrivera Davide successivamente….

2) Theroy : RT vs Chi approach for the dneq_chi (titolo provvisorio e indicativo)

Qui dobbiamo dire che l'approccio consolidato per calcolare il dchi_neq è quello delle simulazioni RT con il pump and probe
esplicito, EQ 7 delle note nella quale si dice che le P_Pp e P_p di omega sono le FT delle polarizzazioni funzioni del tempo
calcolate esplicitamente. Questo approccio però ha il limite che non consente una intepretazione fisica diretta del risultato ottenuto.
Da questo punto di vista riuscire a costruire la stessa risposta con le non-linear chi puo consentire di comprendere e spiegare
l'origine delle features del segnale calcolato perché le chi rappresentano una caratteristica del sistema e il campo esterno interviene
in un modo ben definito dal formalismo (che poi spieghiamo dopo)


In questa sezione vanno inseriti anche gli aspetti tecnici per la costruzione del dchi_neq con le chi di non equilibrio. 
Quanti e quali dettagli è una cosa che devo ancora capire bene, visto che dobbiamo scrivere un articolo (probabilemte un PRB)
e non delle note didattiche. Bisogna sicuramente mettere una versione ridotta dei paragrafi da 1.1 fino a 1.3 per fissare notazioni/convenzioni
e definire le formule fondamentali. 
La sezione 1.4 invece non la metterei nel flusso di questa sezione del lavoro, puo andare al limite in appendice.
Va invece discussa la sezione 1.5 che rappresenta probabilmewnte lo sviluppo teorico piu importantedel lavoro e bisogna arrivare alla eq 32
che rappresenta il punto finale per il calcolo delle dchi_neq. Nel testo metterei i punti salienti e rimandeiamo all'appendici per i dettagli del
modello di accumulated envelope. Per quanto riguarda il caso semplice di envelope at the probe time, metterei in appendice anche quello, ma non 
sono sicuro, vediamo in corso di stesura

3) Pump and probe signal at the ex peak in LiF. 

La figura di riferimento di questa analisi è la  : Fig 13 (2d map RT vs NLCHI - delay -5,5 step 0.25 Even contribution per campo ordine 6*10^7). 
Che andra ricostruita usando le stesse approssimazioni (numero di bande, valore del campo, selezione del contributo even in 
entrambi gli schemi di calcolo). 

Dopo aver copmletato l'analisi teorica in 2 mostriamo il confronto computazionale dei due approcci (solo per il caso HSEX).
In questa sezione ci sono quindi i dettagli computazioanli per ottenere i risultati nei 2 approcci. Qui si puo parlare (senza eccessivi dettagli che
rendono il discorso pesante) delle convergenze rispetto a k, del problema delle chi2 che non sono subito zero e di contro del fatto che il segnale
odd calcolato in RT non è zero e che quindi va rimosso. Anche un commento sull'intensità massima del campo comptabilie con le chi fino al terzo ordine


4)Analysis and discussion of the NEQ signal

Once that the NLCHI signal is known we can analize the various terms to understang their physical origin

	literature : model +  delta pump

Signal anaylsis :

	AC PART : oscillations
	DC PART : derivative feature close to the equilibrium excitonic energy. 
			RT -> interpretation in term of the shift of the exciton paramters
			CHI -> properties of chi3 poles. Why this happens ? 
```
