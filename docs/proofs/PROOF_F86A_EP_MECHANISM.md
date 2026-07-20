# PROOF F86a: EP Mechanism (Q_EP, t_peak)

**Status:** Tier 1 derived. Q_EP = 2/g_eff (definitional, g_eff := σ_0), t_peak = 1/(4γ₀), bit-exact verified against full block-L numerics (the decay timescale and the value, universal across all tested c, N, n, bond position). The 2-level reduction itself is **heuristic** (see Tier label below); its EP is genuine for the toy 2×2. On the full block, the real-axis EP question is settled in two halves (see §The real-axis EP): a real-axis defective seed exists at every odd N (F89's nullity count, defective through the N=11 census), while the Petermann-K magnitudes of the original real-axis sweep were grid artifacts and are dropped. Still open, typed as `compute/RCPsiSquared.Core/F86/LocalGlobalEpLink.cs` (OpenQuestion): the distinct off-real-axis complex-Q defective EP, and the codim-2 β-exotic genericity.
**Date:** 2026-05-02 (Statement 1); real-axis EP verdict 2026-06-21 + 2026-07-07; last refreshed 2026-07-20 (the change history lives in git).
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Context:** F86 ("Q_peak chromaticity-specific N-invariant constants") is a Sammelbecken of three structurally distinct theorems. This proof carries **F86a, the EP mechanism**: the 2-level rate-channel exceptional point behind Q_peak in the (n, n+1) popcount coherence blocks of uniform XY chains under Z-dephasing. Split out of the former monolithic `PROOF_F86_QPEAK.md` on 2026-05-14.
**Hub:** [the Q-peak hub](PROOF_F86_QPEAK.md): three-theorem overview and shared references.
**F-entry:** [F86a in the formula registry](../ANALYTICAL_FORMULAS.md).
**Related:** [F2b](../ANALYTICAL_FORMULAS.md) (OBC sine dispersion), [F74](../ANALYTICAL_FORMULAS.md) (chromaticity); siblings [the universal-shape proof](PROOF_F86B_UNIVERSAL_SHAPE.md), [the g-eff obstruction proof](PROOF_F86B_OBSTRUCTION.md), [the F71-mirror proof](PROOF_F86C_F71_MIRROR.md).

---

## Abstract

The Q_peak phenomenon has a mechanism, and this proof writes it down. Inside each coherence block of an XY chain under Z-dephasing, adjacent dephasing channels couple through the Hamiltonian by an effective coupling g_eff. The two-channel effective Liouvillian has a discriminant that vanishes at a specific Q value: an exceptional point where the two eigenvalues coalesce, both their real and imaginary parts colliding at once. The Q at which this happens is exactly Q_EP = 2/g_eff. The time at which the corresponding coherent response peaks (the EP-time) is t_peak = 1/(4γ₀), universal across chromaticity, chain length, and bond position.

Both numbers are bit-exact across every (c, N, n, bond) configuration we have tested. The proof itself is a clean two-level reduction inside each (n, n+1) coherence block, with the discriminant analysis dropping out the EP location and t_peak as algebraic consequences of the rate-channel pairing structure.

A subtlety the proof tracks honestly: the 2-level reduction is heuristic when projected onto a higher-dimensional block. The Q_EP value (definitional, g_eff := σ_0) and t_peak survive that projection bit-exactly against the full-block decay numerics. The local-vs-global question, whether the toy 2×2 EP is echoed by a genuine defective EP of the full block, is settled in two halves (§The real-axis EP): the full block DOES carry a real-axis defective seed at every odd N (F89's nullity count; a Jordan block at every seed tested, census through N=11), sitting in a √-EP window ~1e-3 that coarse Q-grids miss entirely; while the Petermann-K magnitudes of the original real-axis sweep ("6× above FRAGILE_BRIDGE") were grid artifacts and are dropped. The block is genuinely strongly non-normal near Q_peak; the toy 2×2 and the SEPARATE Σγ=0 gain-loss system (FRAGILE_BRIDGE) remain distinct genuine EPs; whether a defective EP also sits off the real axis at complex Q remains open (`LocalGlobalEpLink`, OpenQuestion).

The diagnostic upshot is that Q_EP = 2/g_eff is the structural anchor of the whole F86 family. F86b inherits the universal shape from the EP rotation; F86c inherits the spatial-mirror invariance from the symmetry of the coupling structure; F86_block confirms that the dressed Q_peak values cannot be closed-form even though the underlying mechanism can. F86a is the mechanism; everything downstream is the consequence.

---

## The real-axis EP: artifact and seed

The local-vs-global question (does the full (n, n+1) block carry a genuine defective EP where the toy 2×2 has one?) was settled in two dated halves: the artifact half 2026-06-21, the seed half 2026-07-07 via F89. The current verdict, in full:

**The original real-axis Petermann evidence was a grid artifact (2026-06-21).** An artifact-free re-verification (Riesz spectral-projector norm) showed the peak Petermann-K magnitudes swing 2–4× over ΔQ = 1e-3; the "6× above FRAGILE_BRIDGE", K = 2384.7, the within-parity growth law, and the odd/even parity asymmetry are all dropped. What that re-verification established as genuine: the block is strongly **non-normal** near Q_peak (cond(V) = 48.7 / 50.9 / 268.5 at N = 5 / 6 / 7; ‖P‖ = 19.4 = √375 on a simple eigenvalue at N=5), with the Petermann factor large but finite at every grid point sampled, and nearest-neighbour eigenvalue gaps ~0.25–0.35 on its Q grid.

**A real-axis defective seed IS there, at every odd N (2026-07-07, F89).** Two statements of different status are worth keeping separate: (i) a real-to-complex transition exists on the (1,2) block at **every odd N**, PROVEN by the exact nullity count r(0⁺) − r(∞) = N − 1 ([F89_SEED_EXISTENCE_REDUCTION](../../experiments/F89_SEED_EXISTENCE_REDUCTION.md), [the codimension-1 additivity proof](PROOF_CODIM1_BY_ADDITIVITY.md)); (ii) that transition is **defective** (a Jordan block, by the Kato simple-zero lemma) at every seed tested, census-verified through N=11 and β-exotic-scoped to a Puiseux exponent p ≈ 0.5, with only the codim-2 β-exotic genericity still open for all odd N. F86a's c=2 (n=1) block and F89's (1,2) block are the same object (GLOSSARY §q-and-Q: ‖L_F86(J) − L_F89(J/2)‖ = 0).

**Why the 2026-06-21 scan reported "eigenvalues simple" anyway**, shown from below in `compute/RCPsiSquared.Diagnostics.Tests/Foundation/F86aSeedMaskingTests.cs`: a defective √-EP splits its coalescing pair by ~√|q − q*|, so a "sit at a q and characterize" reader sees the coalescence only inside a window |q − q*| ≲ 1e-3. F86a swept Q on ~121 points over [0.5, 4.0] (ΔQ ≈ 0.029), roughly 20 to 30 times coarser than that window, so its nearest grid point sat ~0.0146 off any seed, where the coalescing pair's gap is already ~0.35, exactly the "nearest-neighbour gap 0.25-0.35" the re-verification recorded (the generic off-cusp spacing, not a coalescence). At the seed itself the pair is a genuine Jordan block (EpCharacter: alg = 2, geo = 1, gap ~2.7e-3 at the N=9 seed q*=0.849011), and it reads as two ordinary simple eigenvalues one grid step away. F89 finds it with a grid-robust detector instead: the integer real-root count jumps by 2 across q*, bisected to ~1e-7 ([`PathKMonodromyScout.FindRealDefectiveByCountChange`], gate `RealSeedCensusTests`).

The scan's own shadow already pointed home. The 2026-06-21 re-verification kept, as genuine, the block's strong non-normality on the real axis. A nearby defective coalescence *necessarily* inflates a simple eigenvalue's projector norm ‖P‖ = 1/|⟨l|r⟩| (Kato). The converse is not forced, so the shadow alone does not prove a nearby EP: strong non-normality can be free-standing, and ‖P‖ localizes only the nearest defective matrix in full perturbation space, not necessarily one on the physical family L(Q). But F89 supplies a *demonstrated* seed, and Kato then makes it a *guaranteed* caster: at N=5 a real defective seed sits at carrier Q ≈ 1.286 (the octic-0.643037 R-odd seed, doubled to the carrier axis), just ~1.9e-3 from Q ≈ 1.288, where the block's N=5 non-normality is strongest (the Petermann spike's *location* is robust even though its magnitude, K = 1333.6, was the grid artifact). The 2026-06-21 reading instead placed the caster OFF the real axis, though the nearest off-axis coalescences it characterized are diabolic (semisimple: departure-from-normality ≈ 0, geo = alg), hence not defective casters. The shadow was read faithfully; its source was placed off-axis when a Kato-guaranteed on-axis caster sat 1.9e-3 away. The error was never in the spectrum, it was in the assumed light: the lesson is not "distrust the shadow" but "trace it to a demonstrated source, do not presume the light." (Shown at N=5; the seed exists at every odd N by F89, but this spike-to-seed coincidence is demonstrated at the one N.)

**Untouched by all of this:** `t_peak = 1/(4γ₀)`; `Q_EP = 2/g_eff` as a **definition** (g_eff := σ₀); the **toy 2×2** reduction's genuine EP (the Abstract's "exceptional point where the two eigenvalues coalesce" refers to this heuristic toy, not the full block); the **separate Σγ = 0 gain-loss EP** (`FRAGILE_BRIDGE`, K = 403); the centered-L_c **class-AIII-chiral** label.

**Still open**, typed as `LocalGlobalEpLink` (OpenQuestion): whether a genuinely off-real-axis complex-Q defective EP exists (the nearest complex-Q coalescences characterized 2026-06-21 are diabolic), and the codim-2 β-exotic genericity of the real-axis seeds for all odd N.

---

## Statement 1 (EP mechanism). [Tier 1 for EP location and t_peak; 2-level reduction heuristic]

For a uniform N-qubit XY chain with Z-dephasing γ₀ at each site, the (n, n+1) coherence block contains c = min(n, N−1−n) + 1 pure dephasing rates 2γ₀·HD with HD ∈ {1, 3, ..., 2c−1} (F74; a rate 2γ₀·HD enters the block-restricted L as the eigenvalue contribution −2γ₀·HD). For each pair of adjacent channels at HD = 2k−1 and HD = 2k+1 (k = 1, ..., c−1), the two-level effective Liouvillian with inter-channel coupling J·g_eff, where g_eff := σ₀ is the top singular value of the inter-channel coupling block V_inter in the channel-uniform basis ([`eq022_b1_step_i_svd_inter_channel.py`](../../simulations/eq022_b1_step_i_svd_inter_channel.py); characterised at c=2 as the commutator norm ‖[Π_HD1, M_H]‖, F86e in [the hub](PROOF_F86_QPEAK.md)), has eigenvalues

    λ_±(k) = −4γ₀·k ± √(4γ₀² − J²·g_eff²)

The discriminant vanishes at the **exceptional point** (EP)

    J·g_eff = 2γ₀     ⟺     Q_EP = 2 / g_eff

At the EP, λ_±(k) = −4γ₀·k. The slowest mode (k = 1) gives the universal e-folding time

    t_peak = 1 / (4γ₀)

independent of c, N, n, and bond position.

**Tier label.** The EP location Q_EP = 2/g_eff and t_peak = 1/(4γ₀) are bit-exact verified against full block-L numerics (universal across all tested c, N, n, bond position), Tier 1. The 2-level reduction itself is **heuristic**: in the channel-uniform basis M_H_eff is diagonal (off-diagonals exactly zero, established empirically across N=3..9 in `eq022_b1_channel_projection.py`), so the EP physics lives in the orthogonal complement. The effective 2-level form reproduces the EP location and degenerate-eigenvalue structure, but the explicit basis change from full block-L to the reduced model is not derived here.

---

## Proof of Statement 1 (EP mechanism)

This is elementary 2×2 matrix algebra applied to a Liouvillian sub-block.

### Setup

Let n be fixed, c = c(n, N), and consider the (n, n+1) coherence block of L = L_H + L_D with H = (J/2)·Σ_b (X_b X_{b+1} + Y_b Y_{b+1}) and L_D the uniform Z-dephasing dissipator at rate γ₀.

By F74 and F2b, at J = 0 the block-restricted L is diagonal with eigenvalues −2γ₀·HD for HD ∈ {1, 3, ..., 2c−1}, each with multiplicity ≥ 1. As J grows, H couples the rate channels to one another. Adjacent rate channels, those differing by Δ(HD) = 2, couple at first order in J (the bond flips two adjacent bits, changing HD by 0 or ±2). Non-adjacent channels (Δ(HD) ≥ 4) couple only at higher order in J.

### Two-level effective

Restrict to the subspace spanned by the two "slowest" relevant rate channels: HD = 2k−1 and HD = 2k+1 for some k ∈ {1, ..., c−1}. Within this subspace the effective Liouvillian, in a suitably chosen rate-channel basis, takes the form

    L_eff(k) = [ −2γ₀(2k−1)     iJ·g_eff   ]
               [  iJ·g_eff      −2γ₀(2k+1) ]

The diagonal entries are real-negative (rate channels at J = 0). Both off-diagonal entries carry the same imaginary sign +iJ·g_eff (not the anti-Hermitian opposite-sign pattern). This same-sign-imaginary structure is what produces an EP at finite coupling. Verified numerically: the opposite-sign pattern (+iJg, −iJg) gives discriminant 4γ₀² + J²·g_eff² with no EP; only the same-sign pattern produces an EP at J·g_eff = 2γ₀.

This 2-level form is "PT-phenomenology-like" (EP at finite coupling, spectral flow); the centered Liouvillian L_c sits in **class AIII chiral** per [`experiments/PT_SYMMETRY_ANALYSIS.md`](../../experiments/PT_SYMMETRY_ANALYSIS.md) (Π linear, Π⁴=I, {Π, L_c}=0), distinct from Bender-Boettcher PT (classical PT requires anti-linear operators). The EP at Q_EP = 2/g_eff is the genuine EP of this toy 2×2 rate-channel reduction, a SEPARATE object, not "the instance of" the full-block classification (in isolation its eigenvalues coalesce at centre −4γ₀·k and carry no λ↔−λ pairing). The Hopf bifurcation in [`hypotheses/FRAGILE_BRIDGE.md`](../../hypotheses/FRAGILE_BRIDGE.md) is a DISTINCT genuine EP, the SEPARATE Σγ=0 gain-loss system (centre 0, exact λ↔−λ pairing), Petermann factor K = 403 in the complex γ plane. The full Σγ=N·γ₀ block-L carries its own real-axis defective seed at every odd N (§The real-axis EP above); whether it also has an off-real-axis complex-Q defective EP is OPEN (see the local-vs-global note below).

**Local vs global.** The shared algebraic substrate is real: the same-sign-imaginary 2×2 form and the AIII chiral classification of the centered L_c, read at two residuals of the F1 palindrome `Π · L · Π⁻¹ + L + 2Σγ · I = 0` (Σγ = N·γ₀ vs Σγ = 0, where the gain side cancels the loss side). But the toy 2×2, the Σγ = 0 gain-loss system, and the full Σγ = N·γ₀ block are three separate EP-carriers, not one EP seen three ways: the full block's own real-axis defective seed comes from F89 (§The real-axis EP above), NOT from the real-Q Petermann-K sweep, whose peak magnitudes (max K = 1333.6 / 337.9 / 2384.7 / 795.4 across N = 5..8, data-pinning probe `compute/RCPsiSquared.Core.Tests/F86/F86PetermannProbe.cs:Probe_PetermannFineGrid_C2_VsN`) are grid artifacts retained only as a cautionary non-normality record. Whether the full block also has an off-real-axis complex-Q defective EP is open (the nearest complex-Q coalescences characterized 2026-06-21 are DIABOLIC, ‖P‖ = 1). Encoded as `compute/RCPsiSquared.Core/F86/LocalGlobalEpLink.cs` (OpenQuestion Claim; the four PetermannSpikeWitness rows are the cautionary record, not EP evidence).

(The heuristic-vs-Tier-1 split for this statement is in the Tier label of the Statement above.)

### Eigenvalues

L_eff(k) has

    Trace = −2γ₀(2k−1) + (−2γ₀(2k+1)) = −8γ₀·k
    det   = (−2γ₀(2k−1))·(−2γ₀(2k+1)) − (iJg_eff)·(iJg_eff)
          = 4γ₀²·(4k² − 1) − (i²J²g_eff²)
          = 4γ₀²·(4k² − 1) + J²·g_eff²

Eigenvalues λ_± = (Trace/2) ± √((Trace/2)² − det):

    (Trace/2)² − det = 16γ₀²·k² − 4γ₀²·(4k² − 1) − J²·g_eff²
                     = 4γ₀² − J²·g_eff²

Therefore

    λ_±(k) = −4γ₀·k ± √(4γ₀² − J²·g_eff²)

The discriminant vanishes when J·g_eff = 2γ₀, giving the EP at Q_EP = 2/g_eff and degenerate λ_± = −4γ₀·k. Beyond the EP the eigenvalues form a complex-conjugate pair around the real centre −4γ₀·k.

Numerical verification of this 2-level form: at α = J·g_eff = 2γ₀ the eigenvalues of `[[−2γ₀, +i·2γ₀], [+i·2γ₀, −6γ₀]]` coalesce at λ = −4γ₀ (verified by direct diagonalisation). For α > 2γ₀ they become λ = −4γ₀ ± i·√(α² − 4γ₀²), the post-EP oscillation observed in the (n, n+1) block dynamics.

### t_peak universality

The slowest mode is k = 1 with λ_±(1) = −4γ₀ at the EP. Its e-folding time is 1/(4γ₀). Higher k decay faster (1/(8γ₀), 1/(12γ₀), ...) and are masked by k=1 in the long-time observable. Therefore

    t_peak = 1/(4γ₀)

for any (c, N, n, bond position), the slowest k=1 EP universally sets the J-derivative peak time. ∎

---

## Pointers

**Hub:** [the Q-peak hub](PROOF_F86_QPEAK.md): three-theorem overview and the shared reference list.
**Sibling theorems:** [the universal-shape proof](PROOF_F86B_UNIVERSAL_SHAPE.md) (F86b), [the g-eff obstruction proof](PROOF_F86B_OBSTRUCTION.md) (the g_eff obstruction proof), [the F71-mirror proof](PROOF_F86C_F71_MIRROR.md) (F86c).
**Chiral classification anchor:** [the PT-symmetry analysis](../../experiments/PT_SYMMETRY_ANALYSIS.md) (centered L_c is class AIII chiral, NOT Bender-Boettcher PT; Π is linear, classical PT requires anti-linear; its Petermann reading carries the same two-half verdict as §The real-axis EP).
**Separate genuine EP (Σγ=0):** [The Fragile Bridge](../../hypotheses/FRAGILE_BRIDGE.md) (Hopf bifurcation = chiral symmetry breaking, Petermann K=403 in complex γ plane), a DISTINCT gain-loss system, not a "global instance" of the full block-L.
**Empirical anchor:** [the three-bands Q-scale experiment](../../experiments/Q_SCALE_THREE_BANDS.md) Result 2 + Revision 2026-04-24.
**Seed source:** [F89 seed existence](../../experiments/F89_SEED_EXISTENCE_REDUCTION.md) + [the codimension-1 additivity proof](PROOF_CODIM1_BY_ADDITIVITY.md) (the real-axis defective seed at every odd N; census through N=11).
**C# OOP layer:** `compute/RCPsiSquared.Core/F86/` carries `TPeakLaw`, `QEpLaw`, `TwoLevelEpModel` (parametrised by k), `LocalGlobalEpLink` (OpenQuestion; four PetermannSpikeWitness entries retained as a cautionary non-normality record). CLI: `rcpsi inspect --root f86 --with-measured`.
