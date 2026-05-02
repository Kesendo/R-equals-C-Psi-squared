# On the Q Axis and the PTF Lesson

**Status:** Reflection. Synthesises the F86 universal-shape finding (2026-05-02) with the existing chiral classification of the Liouvillian (Π class AIII per [PT_SYMMETRY_ANALYSIS](../experiments/PT_SYMMETRY_ANALYSIS.md)), the Inside-Observability theorem ([PRIMORDIAL_QUBIT §9](../hypotheses/PRIMORDIAL_QUBIT.md)), and the two-time reading ([ON_TWO_TIMES](ON_TWO_TIMES.md)). The methodology is calibrated against the PTF retraction arc (EQ-014, 2026-04-19) which had structurally the same convergence-of-many-threads pattern.

**Date:** 2026-05-02
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Depends on:** [F86 in ANALYTICAL_FORMULAS.md](../docs/ANALYTICAL_FORMULAS.md), [PROOF_F86_QPEAK](../docs/proofs/PROOF_F86_QPEAK.md), [PT_SYMMETRY_ANALYSIS](../experiments/PT_SYMMETRY_ANALYSIS.md), [FRAGILE_BRIDGE](../hypotheses/FRAGILE_BRIDGE.md), [PERSPECTIVAL_TIME_FIELD](../hypotheses/PERSPECTIVAL_TIME_FIELD.md), [ON_TWO_TIMES](ON_TWO_TIMES.md), [PRIMORDIAL_QUBIT](../hypotheses/PRIMORDIAL_QUBIT.md), EQ-014 and EQ-022 in [EMERGING_QUESTIONS](../review/EMERGING_QUESTIONS.md), `memory/project_retraction_lesson_ptf_to_f86`.

---

The Q axis kept opening doors today. Each one led to a structural anchor that already existed in the repo, and the convergence of all of them on the EP at the slowest rate-channel pair felt, for a few hours, like enough evidence to claim a closed form for Q_peak. It was not. The closed-form claims got retracted within the same session, one after the other, against extended-N data. What survived is the *shape* around the peak, not the peak's location, and it survived because it is a symmetry, not a number.

This is the same arc as PTF, ten days earlier. The structural lesson of both arcs is the same. Worth writing down before the next door opens.

---

## 1. Two retractions, one survivor each

### PTF (April 18 to April 27, 2026)

The Perspectival Time Field reading (`hypotheses/PERSPECTIVAL_TIME_FIELD.md`) proposed that under a J-coupling defect in an N=7 XY chain, each site has its own time-rescaling rate α_i, and that the seven α_i across sites satisfy a state-independent closure law

    Σ_i ln(α_i)  =  0     (within ±0.05 across all tested initial states)

The empirical pattern was striking. The painters-around-a-mountain image carried the interpretation. Under the surface, the closure looked like a theorem.

EQ-014 closed it. Direct first-order extraction `Σ f_i = lim_{δJ → 0} Σ ln(α_i)/δJ` from exact RK4 evolution gave +0.05 to +1.29 across initial states — state-dependent, nonzero. The ±0.05 closure at δJ = 0.1 was a fortuitous combination of small first-order coefficients (some states) and partial second-order cancellation (others), not an exact identity. The theorem route closed.

What survived was the chiral mirror law

    Σ f_i(ψ_k)  =  Σ f_i(ψ_{N+1−k})     (machine-precision exact across N=5, 7, 8)

This is symmetry-enforced: K_1 = diag((−1)^i) is the chiral sublattice symmetry of the open XY chain, and the J-defect perturbation V_L is K_1-odd. Eigenvector mixing peaks at the K_1-fixed point ψ_4 at N=7, where Σ f_i = +2.14 (more than twice the outermost pair). Symmetry forces the equality. The closure law was the wrong-level claim; the chiral mirror is the right-level one.

### F86 (this session, 2026-05-02)

After Q_SCALE established the c-specific N-invariant Q_peak constants (c=3 → 1.6, c≥4 → 1.8, per relative-J convention), the morning's per-bond fine-grid scan suggested two clean closed forms:

- Q_peak(Endpoint, N) = csc(π/(N+1))  (chain-edge sine-mode anchor via F2b)
- Q_peak(Interior, c=3) → csc(π/5)   (pentagonal asymptote, golden-ratio family)

Both matched at N=7 within sub-percent. The csc(π/5) value at c=3 N=7 was striking enough to warrant a memory entry, a framework primitive, and a section in F86. The pentagonal anchor felt structurally meaningful — at N=4 the OBC dispersion is exactly ±φ, ±1/φ, and the n_XY=2 sector Liouvillian Im(λ) values decompose into integer combinations of {φ, 1/φ, 1, √5} (verified bit-exact in [`eq018_golden_ratio_check.py`](../simulations/eq018_golden_ratio_check.py)).

Then the N=8 data arrived. The c=3 trajectory was 1.566 → 1.689 → 1.743 → **1.750** — it crossed csc(π/5) = 1.7013 between N=6 and N=7 and continued past it. The Endpoint formula was even simpler: at N=5 the observed value 2.40 differs by 20 % from csc(π/6) = 2.0, at N=6 by 9 %, at N=7 by −3 %, at N=8 by −13 %. The "1.4 % match at N=7" reported earlier was itself a dQ=0.05 grid-snap artefact; with parabolic interpolation the N=7 Endpoint value is 2.53, not 2.65, and csc(π/8) = 2.6131 is *off* in the other direction.

Both closed forms retracted. The framework primitives `q_peak_endpoint(N)` and `Q_PEAK_INTERIOR_C3_ANCHOR` that had been promoted earlier in the same session were removed.

What survived was the universal shape under relative-Q normalisation:

    K(Q) / |K|_max  =  f(Q/Q_EP)

The pairwise residual under `x = (Q − Q_peak)/Q_peak, y = K(Q)/|K|_max` is 21× tighter than under absolute-Q shift. Across all six tested cases (c=3 N=5..8, c=4 N=7,8) the y values at common x cluster within 1-3 %. The half-width on the left side has

    HWHM_left / Q_peak  ≈  0.756 ± 0.005

universal across the tested range. The position is chain-specific; the shape is not.

---

## 2. The structural hierarchy

The convergent threads, organised vertically by abstraction level:

**Top: existing classification.** Π is class AIII chiral (Altland-Zirnbauer framework), established at length in [`experiments/PT_SYMMETRY_ANALYSIS.md`](../experiments/PT_SYMMETRY_ANALYSIS.md). The classification is order-4, linear (not anti-linear), and gives Π·L·Π⁻¹ = −L − 2Σγ·I for all N (F1). The fragile-bridge Hopf bifurcation in [`hypotheses/FRAGILE_BRIDGE.md`](../hypotheses/FRAGILE_BRIDGE.md) is the same chiral classification at a different scale: Petermann factor K = 403 signals an exceptional point in the *complex* γ plane.

**Middle: 2-level local instance (today's reduction).** For two adjacent rate channels at HD = 2k−1 and HD = 2k+1 in the (n, n+1) coherence block, the effective Liouvillian takes the form

    L_eff − (trace/2)·I  =  [ −Δ/2     +iJ·g_eff ]
                            [ +iJ·g_eff   +Δ/2  ]      with Δ = 4γ₀

The same-sign-imaginary off-diagonal pattern admits an EP at finite J·g_eff = 2γ₀, with degenerate eigenvalue Re(λ) = −4γ₀·k. This is "PT-phenomenology-like" (EP at finite coupling) but algebraically belongs in the chiral class above. The opposite-sign pattern (+iJg, −iJg) gives discriminant 4γ₀² + J²g² with no EP — verified numerically. The local-2-level-EP at Q_EP = 2/g_eff is the rate-channel instance of the chiral classification; the global Hopf bifurcation is the complex-γ-plane instance. Whether the two are connected algebraically is open.

**Lower: universal shape (today's finding).** The 2-level eigenvector rotation `tan θ = J·g_eff / 2γ₀ = Q/Q_EP` makes every probe-overlap observable a function of Q/Q_EP alone. Hence K(Q) / |K|_max = f(Q/Q_EP) for some universal function f. Q_peak ≈ Q_EP for the slowest channel pair; g_eff varies with chain (c, N, bond position) but the shape in `(Q − Q_peak)/Q_peak` coordinates does not. HWHM_left/Q_peak ≈ 0.756 is the numerical witness.

**Interpretive: the inside view.** Inside-Observability (`hypotheses/PRIMORDIAL_QUBIT.md` §9) says Q is the only inside-measurable quantity; γ₀ alone is invisible. The shape universality reads cleanly through this: Q_peak chain-specific = the chain knowing itself through its own peak position; universal shape = the framework-level fact about the chain that the chain-as-its-own-instrument cannot identify in absolute units. t_peak = 1/(4γ₀) is the felt-time horizon (`reflections/ON_TWO_TIMES.md`) — γ₀-time flows invisibly, but the EP-time-scale that lives inside the standing-wave envelope is set by it.

The hierarchy reads top-down (algebra → local form → universal shape → reading), or bottom-up (reading → universal shape → local form → algebra). The same structure recurs across levels. That's part of why it kept opening doors.

---

## 3. What the PT-like form is and is not

The "PT-symmetric-like" framing of the 2-level effective is half-right and needs the qualifier the existing classification provides. Π is *linear*, classical Bender-Boettcher PT requires *anti-linear* (P linear, T anti-linear, PT anti-linear). This was established at length in `experiments/PT_SYMMETRY_ANALYSIS.md` (lines 31, 105, 111, 392): Π is class AIII chiral, not PT.

The 2-level off-diagonal pattern same-sign-imaginary (+iJg_eff, +iJg_eff) is non-Hermitian and admits an EP at finite coupling. This *phenomenology* matches PT-symmetric Hamiltonians (real eigenvalues below EP, complex pair above), but the algebraic mechanism is different: the same-sign-imaginary structure is what the L matrix actually does in the rate-channel-effective basis, and it sits inside the chiral classification of the full L. The phenomenology is PT-like; the algebra is chiral.

This matters because the chiral classification gives us infrastructure the PT label does not: F1 spectral mirror, Π² parity sectoring, Hopf bifurcation in FRAGILE_BRIDGE, complex-γ-plane EP via Petermann K=403. The 2-level EP at Q_EP = 2/g_eff is a local instance of all of this. The phenomenological resemblance to Bender PT is suggestive but should not guide further work; the structural classification should.

---

## 4. The named structural law: EP-rotation universality

The F86 analog of PTF's chiral mirror law. PTF's surviving Tier-1 law:

    Σ f_i(ψ_k)  =  Σ f_i(ψ_{N+1−k})     (K_1 sublattice symmetry of the open XY chain)

F86's surviving Tier-1 *candidate*:

    K(Q) / |K|_max  =  f(Q/Q_EP)         (2-level EP rotation symmetry)

Both are symmetry-enforced. PTF's symmetry is chiral on sites (K_1 = diag((−1)^i)); F86's is the 2-level EP rotation `tan θ = Q/Q_EP` on rate channels. PTF's law is exact at machine precision; F86's law is a Tier-1 candidate witnessed by 21× residual tightening across six (c, N) cases and HWHM_left/Q_peak ≈ 0.756 stability to 1.8 %.

The closed-form value of HWHM_left/Q_peak (0.756 looks numerically close to 3/4 = 0.75 and 2/√7 = 0.7559) follows from the 2-level eigenstructure but has not been derived analytically. Per the lesson below: do not promote the numerical match to a closed form before deriving f(x).

---

## 5. The methodological closing

The lesson is recorded in `memory/project_retraction_lesson_ptf_to_f86`. Decision rule for promoting any closed-form claim to Tier 1:

1. The claim must hold at more than one (c, N) anchor point. A single-point match is a trajectory crossing.
2. Extend the sample at least one step in N (or in the relevant parameter) before promoting.
3. Distinguish empirical pattern at precision P from exact structural identity. A law that works to 5 % at one parameter point is not a theorem.
4. Look for a symmetry that explains the pattern. If you can name the symmetry (K_1 chiral mirror, EP rotation, ...) the survivor is structural; if you cannot, the closed-form claim is fragile.
5. Numerical-fit-looking closed forms need analytical derivation from the underlying eigenstructure before promotion. The numerical witness is necessary but not sufficient.

PTF spent eleven days (April 18 to April 27) between the closure-law over-claim and the chiral-mirror-law survivor. F86 compressed the same arc into a few hours within one session. The repetition is the lesson. When many threads converge — class AIII chiral, 2-level EP, universal shape, Inside-Observability, two-time reading, F1 palindrome, F2b OBC dispersion, golden-ratio fixed point at N=4 — the convergence itself becomes pressure to claim more than the data supports. Naming the surviving symmetry is the way to absorb the convergence without over-claiming.

The Q axis is the inside-observable. The shape around its peak is what the chain cannot tell about itself. The peak's location is what the chain *can* tell about itself. Two facts, one axis, both inside-readable in the right coordinate. The framework is consistent with the reading; the work is to keep the right level.
