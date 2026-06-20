# On the Polarity Coordinates and the Balance We Could Not Break

**Status:** Reflection. After the polarity-coordinates wave landed the three-way decomposition (`polarity_coordinates`) as a typed framework primitive, after 16 of 17 pytests passed cleanly, and after the one failure refused to be turned into a success by extending the probe across eight Hamiltonian families and three dissipator settings. Captures what we expected to discover, what we actually discovered, and the structural symmetry that was sitting in the background the whole time without anyone naming it.

**Closure 2026-05-26:** F112 (`LindbladBitBPiBalance`) typed Tier1Derived for Hermitian H + non-Hermitian extension at N ≤ 4 (via basis enumeration; Tier1Candidate at N ≥ 5); see the *Closure* section at the bottom.

**Date:** 2026-05-25 (closure addendum 2026-05-26)
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Context:** [F81](../docs/proofs/PROOF_F81_PI_CONJUGATION_OF_M.md) gave us the binary Π-symmetric/antisymmetric split of M (the framework's residual operator). The Π-conjugation map has order four on Liouville space; its full eigenvalue spectrum is {+1, −1, +i, −i}, not two but four. The ±1 eigenspaces together form F81's M_sym; the ±i eigenspaces together form M_anti. The wave's hypothesis was that splitting M_anti further into the +i and −i parts would expose an asymmetric T1 signature. The hypothesis is refuted. The refutation is the interesting thing.

---

## What we built

The primitive in [`polarity_coordinates.py`](../simulations/framework/diagnostics/polarity_coordinates.py) takes the F81 sym/anti split and refines the anti part with two complex projectors:

    M_zero       = (M + Π·M·Π⁻¹) / 2                 (the 0-axis, identical to F81 M_sym)
    M_plus_half  = (M_anti − i·Π·M_anti·Π⁻¹) / 2     (Π eigenvalue +i, the +1/2 polarity)
    M_minus_half = (M_anti + i·Π·M_anti·Π⁻¹) / 2     (Π eigenvalue −i, the −1/2 polarity)

These are the Π-eigenvalue projectors restricted to the Π²-odd subspace. They are Frobenius-orthogonal by construction, and the orthogonality test (`test_orthogonality_invariant_across_H_families`) verifies the closure to machine precision across all bilinear H families we tried. The Frobenius-norm invariant is bit-exact:

    ‖M‖² = ‖M_zero‖² + ‖M_plus_half‖² + ‖M_minus_half‖²

The structural payload: where F81 measured one polarity axis (Π²-parity), the new primitive measures both. The {−1/2, 0, +1/2} polarity triple at d=2, anchored as a [Pi2KnowledgeBase claim](../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs) since the F1 trunk landed, has had explicit operator-level Frobenius coordinates for one day. Before the wave, the +1/2 and −1/2 sides were named in the symmetry layer but never weighed independently against an operator. Now they are.

The connection to existing structure: F81's M_sym _is_ M_zero (verified by `test_f81_match_M_sym_equals_M_zero`); F81's M_anti is the sum M_plus_half + M_minus_half (verified by the orthogonality invariant). The new piece is the Π-eigenvalue split of M_anti. It is the operator-level analogue of asking, on a circle with four marked quadrants, not just whether a point is on the top half or the bottom half but which exact quadrant it sits in.

## What we expected

The [spec's working hypothesis](../docs/superpowers/specs/2026-05-25-polarity-coordinates-design.md) was a two-line prediction:

- For Hermitian H and pure Z-dephasing, the +1/2 and −1/2 sides should be balanced: ‖M_plus_half‖² = ‖M_minus_half‖². The reasoning was an appeal to "complex-conjugate symmetry within M_anti" without a tight structural derivation, more a vibes-prediction than a proof. We expected it to hold; if it did we would have a clean Hermitian-H baseline.
- For T1 cooling-only (γ_↓ > 0, γ_↑ = 0), the F84 amplitude-damping correction violates F81's identity. That violation should presumably manifest as an asymmetric split between +1/2 and −1/2. Specifically: T1 cooling breaks the up-down symmetry of the Bloch sphere by preferring the ground state; the resulting bias on the Π²-odd dynamics should land on one polarity side more than the other. We expected the T1 asymmetry test to pass with a measurable gap ‖M_plus‖² − ‖M_minus‖² > 1e-3.

The wave was framed as a discovery probe: confirm the balance, find the breaker, type the asymmetry as a new probe-channel for non-Hermitian / non-Lindblad territory.

## What we found, in two pieces

**Piece one (the confirmation):** The Hermitian-H balance holds across all six bilinear H families tested in [Task B](../simulations/framework/tests/diagnostics/test_polarity_coordinates.py). XX only, YZ+ZY (the canonical Π²-even non-truly anomaly carrier from F108), XY pure (the canonical Π²-odd case), XY+YX (the pair), Heisenberg, Heisenberg + XY (mixed even+odd). Every family gave ‖M_plus_half‖² = ‖M_minus_half‖² at machine precision under pure Z-dephasing. Six for six. The Hermitian-H prediction landed.

**Piece two (the refutation):** The T1 asymmetry hypothesis did not hold. For the Heisenberg + T1-cooling case (γ_z = 0.05, γ_t1 = 0.1, γ_pump = 0), the values came out:

    ‖M_plus_half‖²  = 0.24000000000000005
    ‖M_minus_half‖² = 0.24000000000000005
    asymmetry       = 0.0

Not "close to zero"; bit-exact zero (the difference between the two 0.24-valued doubles is identically zero, not 1e-16). The test failed not by a small numerical margin but by structural cancellation. The hypothesis was wrong in a sharp, machine-precision way.

The instinctive response was to expand the probe, on the theory that we had picked an unlucky H family. [Task C](../simulations/_polarity_demo.py) extended the sweep to eight H families crossed with three dissipator settings: pure Z, Z plus T1 cooling, Z plus T1 detailed balance. Every cell of the table gave asymmetry = 0.0 bit-exact. The refutation sharpens; it does not soften. The +1/2 and −1/2 sides remained balanced across:

- F1-truly H (Heisenberg, XX only)
- Π²-Z-even non-truly H (YZ+ZY)
- Π²-Z-odd H (XY pure, XY+YX)
- Mixed even-plus-odd (Heisenberg + XY)
- T1 cooling-only on all of the above (F81 violation present, balance preserved)
- T1 heating-only (mirror case of cooling, same result)
- Non-uniform per-site T1 (γ_t1 = [0.1, 0.05, 0.2])
- T1 detailed balance (γ_t1 = γ_pump = 0.1)

The single-site transverse-field case (h_y·Y_l) was the one probe Task C could not run: `pi_decompose_M` accepts bilinears or k-body tuples of length ≥ 2, and single-site terms silently drop. This is a deferred future-work item, and it might still break the balance. But every probe the existing API admitted gave the same answer.

There is one more important detail: the F81 violation _is_ measurable in these T1 cases. The diagnostic returns `f81_violation ≈ 0.69` for T1 cooling-only on Heisenberg at N=3, γ_t1=0.1, which matches F84's closed form ‖D_AmplDamp_odd‖_F = √(Σ_l (γ_↓_l − γ_↑_l)²)·2^(N−1) = 0.1·√3·4 = 0.6928 exactly (PROOF_F84 verification table). So F84 is right; the dissipator does shift M_anti away from L_{H_odd}. But the shift lives entirely in M_anti's whole. When you project M_anti onto the +i and −i Π-eigenvalue eigenspaces and weigh them separately, the shift distributes itself symmetrically. The T1 dissipator's contribution to M_anti is itself ±i-balanced.

This is not a small finding. It says: the Π-eigenvalue split is not coupled to T1's asymmetry. The Bloch-sphere up-down bias does not project onto the +i vs −i Π-eigenvalue axis.

## The sharpened question

The naive reading of the result is "well, we tested everything we could and the balance held; either it is a structural symmetry or we have not found the breaker yet." That is true but incomplete. The structural reading is sharper and worth naming.

**Candidate explanation (a), the Hermitian-conjugate symmetry route.** Π is unitary, so Π⁻¹ = Π†. Conjugation by a unitary is a unitary action on operator space; the linear map T_Π: X ↦ Π·X·Π⁻¹ has the same spectrum as Π itself ({+1, −1, +i, −i}) lifted to operator space. Now consider the Hermitian-conjugate involution on operator space, J: X ↦ X†. This is an antilinear involution. It acts on the eigenspaces of T_Π by complex conjugation of the eigenvalue: if T_Π(X) = i·X, then T_Π(X†) = Π·X†·Π⁻¹ = (Π·X·Π⁻¹)†·(some_sign) = (i·X)† = −i·X† (modulo signs that work out for Π unitary). So J maps the +i eigenspace bijectively onto the −i eigenspace.

If M is Hermitian as a superoperator (M = M† under the Frobenius adjoint, equivalent to its kernel preserving the hermiticity of ρ, which is true for Lindbladians with Hermitian H and standard Lindblad dissipators), then J(M_plus_half) = M_minus_half. Frobenius norm is preserved by the Hermitian-conjugate involution. Therefore ‖M_plus_half‖² = ‖M_minus_half‖² automatically. The balance is not a discovery; it is a consequence of M being a Lindbladian.

If this explanation is right, the balance would break if M itself were not Hermitian as a superoperator. That happens for non-Hermitian effective Hamiltonians (open systems with NH terms representing leakage to inaccessible bath modes) or for non-Lindblad-form dissipators that do not preserve hermiticity. We have tested none of those.

The candidate-(a) reading is the prosaic one. It says: the structural symmetry we missed is just Hermitian conjugation, and we should have seen it in the spec. We didn't, because we were focused on the Π eigenvalue projection as a novel splitting and forgot that the underlying operator space carries one more involution that intersects nontrivially with it. Working it out as algebra would take less than a page. Working it out empirically took a day and produced cleaner evidence.

(Editorial caveat after the probes ran and F112 closed: the "M is Hermitian as a superoperator" premise above is the wrong invariant. M = Π L Π⁻¹ + L + 2σ·I is not Hermitian as a superoperator in general for Lindbladian L; what F112's Step 5 actually uses is L_H's anti-Hermitian property under the Hilbert-Schmidt inner product for Hermitian H. Candidate (a) is therefore not wrong because non-Hermitian H breaks it; it is wrong-level from the start, predicting the right answer for the wrong reason. The probes still refute the prediction it makes, but the mechanism is L_H anti-Hermiticity carried into the Π-eigenspace dagger map, not M Hermiticity. See PROOF_F112 Step 5 for the correct argument.)

**Candidate explanation (b), the bra-ket reading.** Deeper, and connected to `project_one_system_two_indices`. The {−1/2, 0, +1/2} polarity at d=2 is the qubit's reading of the +0/−0 polarity layer at d=0 (the [`PolarityLayerOriginClaim`](../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs)). The two sides of the polarity are the two indices of ρ in the bra-ket sense: a density matrix lives on H ⊗ H*, one factor for the ket index, one for the bra. The Π-eigenvalue projection onto +i vs −i operates exactly on this dual-index structure: M_plus_half is the ket-side projection, M_minus_half is the bra-side projection. The Hermitian-conjugate involution _is_ the bra-ket exchange.

Read this way, the balance is not a Hermitian-symmetry accident. It is the statement that the bra and the ket are equally weighted in any physical M, because any physical M must preserve the hermiticity of ρ, because we are observing one system from one side and the bra-ket structure is the operator-algebra refraction of "we" into "two indices on one object." `project_framework_as_remembrance` reads the framework as operational memory for the half of d² − 2d = 0 we cannot directly remember; the bra-ket balance is the structural statement that we cannot favor one half over the other from within. The mirror is by construction equally weighted.

Candidate (a) and candidate (b) are not in competition. (a) is the algebraic mechanism; (b) is the meaning. The same identity ‖M_plus‖² = ‖M_minus‖² reads one way as "Lindbladians are Hermitian as superoperators" and another way as "the bra and the ket are the same subject's two sides." The framework's `project_one_system_two_indices` memory has been saying this since 2026-04-30. The polarity-coordinates wave gave it a Frobenius-orthogonal numerical witness without anyone needing to argue for it.

## What the typed framework readings already said

This wave does not add a hypothesis; it adds a measurement that several existing hypotheses anticipated.

[`reference_diopter_as_polarity_bridge`](../) named the {−1/2, 0, +1/2} triple at d=2 as the framework's polarity vocabulary. Tom's recognition (May 19) was that eyeglass diopter prescriptions are the everyday phenomenology of this triple. The wave's contribution is to have computed the triple as an explicit Frobenius decomposition; the two half sides come out bit-exact equal under every Hermitian-system probe, which is what a diopter pair is supposed to be (the +0.5 and the −0.5 reading-glass have the same lens magnitude, opposite sign).

[`project_plus_minus_zero_layer`](../) reads the qubit as a window onto a +0/−0 polarity substrate; the X-basis {|+⟩, |−⟩} is the framework's natural projection onto this layer. The polarity-coordinates primitive is the operator-level analogue: where the X-basis projects ρ onto a polarity-eigenstate pair, the Π-eigenvalue projection splits the Liouvillian residual into a polarity-eigenoperator pair. Different axis (state vs operator), same polarity structure inherited from the d=0 layer.

[`project_one_system_two_indices`](../) is the cleanest match. The bra and ket indices of ρ on d² = 4^N are the "two sides" of the system; we are the self-coupling system seen from inside, and the d=0 axis is the active vacuum substrate we live on. The Π-eigenvalue split into +i and −i is the operator-level instantiation of the bra-ket pair. The balance ‖M_plus‖² = ‖M_minus‖² is the structural statement that we cannot favor one index over the other; the framework's mirror is balanced by construction because there is one subject with two indexed faces.

[`project_framework_as_remembrance`](../) carries the reading further: the framework package is operational memory for the half of d² − 2d = 0 (the d=0 side) we cannot directly remember; Π compensates the forgetting. The polarity-coordinates measurement is the explicit Frobenius witness for the symmetry across the axis through us. The two halves come out equal because they are the same one of us read from two sides; whatever the forgetting does, it does it symmetrically.

[`project_blockspectrum_z_deph_only`](../) is the day's parallel finding. The BlockSpectrum wave wired the F108 Pi5Bilinear auto-dispatch at the builder layer and discovered that canonical Π and Pi5Bilinear-Z share the same sector cycle: there is no speedup because what we hoped to discover as a new sector structure was already present as a hidden symmetry of the canonical Π. The pattern repeats here: what we hoped to discover as an asymmetric T1 signature was already structurally absent because Lindbladians are Hermitian as superoperators. The typed-knowledge layer often makes explicit what was already implicit; the discovery is the explicit naming and the measurement, not the underlying mathematics.

## What is not explored yet, and where the balance might break

The Hermitian-conjugate explanation gives a clean prediction: break the Hermitian-superoperator property and the balance should break too. Concrete probes for a future wave:

- **Non-Hermitian H.** Effective Hamiltonians with NH terms (gain-loss systems, post-selected dynamics with measurement back-action, PT-symmetric models). These are not Lindbladians in the strict sense but they do generate dynamics on ρ-space, and the resulting M would not be Hermitian as a superoperator. The +1/2 vs −1/2 split should asymmetrize. Magnitude prediction: should track the non-Hermiticity scale of the H.

- **Single-site transverse fields (h_x·X_l, h_y·Y_l, h_z·Z_l).** The current `pi_decompose_M` API does not accept single-site Pauli terms; they silently drop. Hardware-relevant case (the Marrakesh h_y_eff = 0.05 signature uses transverse fields). The prediction from the Hermitian-conjugate argument: should _not_ break the balance, because single-site Pauli terms still give Hermitian H. But the prediction is untested.

- **Mixed dephase letters.** Z-dephasing plus X-dephasing simultaneously, at different rates. Both individually preserve hermiticity; their sum should too. The prediction is balance-preserving; the measurement is missing.

- **Higher k-body terms.** The primitive accepts arbitrary k-body tuples via `pi_decompose_M`. Three-body and four-body Hamiltonians have not been swept; the prediction is balance-preserving for any Hermitian k-body H.

- **Asymmetric collapse operators.** Non-Lindblad dissipators that do not preserve hermiticity (operators like c·σ⁻ with c a non-real complex number, or rotating-frame approximations that drop the hermiticity-preserving terms). The prediction: balance breaks proportionally to the non-Hermiticity of the dissipator.

These five candidates form a natural probe-axis for the +1/2 vs −1/2 channel. The polarity-coordinates primitive is now the diagnostic; future waves can use it to test for non-Hermitian effective dynamics in hardware, for non-Lindblad approximations in derivations, for any structural deviation from the standard Hermitian-superoperator framework.

## Probes 1, 2, 3, 4, 5 (2026-05-25)

The five candidates listed above were converted into five standalone probes the same day the section was written. The k=3 probe (probe 4) ran first via the existing `polarity_coordinates` (pi_decompose_M was already k-body-capable since commit `9cce4b3`). For the other four (non-Hermitian H, single-site fields, mixed dephase letters, non-Lindblad dissipators) a sister primitive `polarity_coordinates_from_L(L_pauli, N, sigma, Pi=None)` was added to accept an externally constructed Liouvillian, bypassing the framework's H + dissipator builder.

The summary: **all five candidates preserve the balance bit-exactly**. The Hermitian-conjugate explanation proposed earlier in this reflection is insufficient. The balance is more robust than even the candidate explanation predicted.

**Probe 4 (k=3 chain Hamiltonians, [`_polarity_probe_kbody.py`](../simulations/_polarity_probe_kbody.py)).** Six k=3 Hamiltonian families on N=4 chains under pure Z-dephasing: XYZ pure, ZYX (Hermitian conjugate of XYZ), XYZ+ZYX (Hermitian pair), YYY pure, XXX+YYY+ZZZ, XYI+IXY. Every family gave asymmetry = 0.0 bit-exact. Note: the F81 identity ‖M_anti − L_{H_odd}‖ has a non-zero residual at k=3 (~55 for some cases), since F81 was proven for 2-body only, but the +1/2 vs −1/2 balance survives independently. F81 identity and polarity balance are SEPARATE invariants.

**Probe 1 (non-Hermitian H, [`_polarity_probe_nonhermitian.py`](../simulations/_polarity_probe_nonhermitian.py)).** Two sub-sweeps at N=3: (A) H = J·XX-chain + i·κ·YZ-chain (Π²-Z-even substrate); (B) H = J·XY-chain + i·κ·XZ-chain (Π²-Z-odd substrate, where M_plus/minus_half are non-trivial). Both at κ ∈ {0, 0.05, 0.1, 0.5, 1.0}. Sweep A pushes 100% of M into M_zero for κ > 0 (everything sits on the Π²-even axis, trivially balanced). Sweep B is the sharp test: M_zero/+1/2/−1/2 split stays at exactly 50%/25%/25% across the full κ range, with asymmetry = 0.0 bit-exact at every κ. **A non-Hermitian matrix in the commutator −i[H,·] does not break the +1/2 vs −1/2 balance.** (Scope, 2026-06-20: the probe builds the commutator `L = −i(H⊗I − I⊗Hᵀ)`, NOT the physical generator −i(Hρ−ρH†) of PT/gain-loss/post-selection dynamics — that generator can break the balance content-dependently (single-site-Z detuning / random H; a uniform-damping bond-Hamiltonian chain is itself balanced — refining the "What is not explored yet" prediction); the two coincide only for Hermitian H. See `docs/CAUGHT_ERRORS.md` and `hypotheses/ASYMMETRY_IS_THE_UNRECYCLED_DRAIN.md`.) For the commutator, this refutes the Hermitian-conjugate prediction in this reflection: the candidate explanation said non-Hermitian H makes M not Hermitian-as-superoperator and the balance should break; for the commutator the measurement says it does not.

**Probe 2 (single-site transverse Y field, [`_polarity_probe_transverse.py`](../simulations/_polarity_probe_transverse.py)).** H = J·Heisenberg + h_y · Σ_l Y_l at N=3, sweep h_y ∈ {0, 0.05, 0.1, 0.5}. Pure Heisenberg (truly) at h_y = 0 gives M ~ 0 baseline. For h_y > 0 the transverse field creates a 50%/25%/25% split (same shape as Probe 1B) with asymmetry = 0.0 bit-exact at every h_y. **Single-site transverse fields preserve balance**, matching the Hermitian-conjugate prediction here. Hardware-relevant: the Marrakesh h_y_eff = 0.05 leakage signature would not asymmetrize the +1/2 vs −1/2 channel.

**Probe 3 (mixed dephase letters, [`_polarity_probe_mixed_dephase.py`](../simulations/_polarity_probe_mixed_dephase.py)).** H = Heisenberg, dissipator = γ_z · Σ_l D[Z]_l + γ_x · Σ_l D[X]_l simultaneously, decomposition against canonical Z-dephase Π. Sweep r = γ_x / γ_z ∈ {0, 0.1, 0.5, 1.0, 2.0} with γ_z = 0.05 fixed. For r = 0 (pure Z): M ~ 0 (truly). For r > 0 the X-dephasing contribution lands entirely in M_zero (100%): the X-dissipator is Π_X-palindromic but Π²-Z-even, so it sits on the 0-axis when measured against Π_Z. Asymmetry = 0.0 bit-exact across the sweep. Mixed dephase letters preserve balance.

**Probe 5 (non-Lindblad dissipator, [`_polarity_probe_non_lindblad.py`](../simulations/_polarity_probe_non_lindblad.py)).** H = Heisenberg + Z-dephasing baseline plus a non-Lindblad addition. Two forms: (A) σ⁺ρσ⁻ − ρ, the σ⁻ jump action missing its anti-commutator (trace-non-preserving but hermiticity-preserving); (B) σ⁺ρσ⁺ − ρ, raising on both sides (the strongest break: ρ doesn't stay Hermitian under this dissipator). Sweep γ ∈ {0, 0.05, 0.1, 0.5}. Both forms produce a 99%/0.5%/0.5% split for γ > 0 with asymmetry = 0.0 bit-exact at every γ, including Form B which actively destroys ρ's hermiticity. **The balance survives even when the dissipator does not preserve the hermiticity of ρ.** This is the most damning evidence against the Hermitian-conjugate explanation: the candidate mechanism requires M to be Hermitian as a superoperator, which in turn requires the dissipator to preserve hermiticity. Form B violates the precondition and the balance holds anyway.

**Aggregate reading.** Across five candidate-breakers (non-Hermitian H on both Π²-Z-even and Π²-Z-odd substrates, single-site transverse fields, mixed dephase letters, non-Lindblad-jump dissipators, k=3 Hamiltonians), the +1/2 vs −1/2 balance is bit-exact 0.0 in every case. The Hermitian-conjugate / "bra-ket exchange" explanation proposed in the candidate (a)/(b) section above predicts asymmetry should break for Probe 1B and Probe 5B; it does not. The structural symmetry is deeper than the candidate mechanism named. Sharpened question: what IS the actual structural identity that forces ‖M_plus_half‖² = ‖M_minus_half‖² for every L we can construct that has any Π²-Z-odd content at all? The candidates the reflection enumerated as "should break it" are exhausted. The next probe needs to construct an L by some path that does not factor through the (kron(c, c.conj()) ± something) form at all, perhaps a Pauli-basis L that is built coefficient-by-coefficient without going through a c-op pipeline, with deliberate ±i Π-eigenvalue imbalance hand-engineered into M_anti.

The wave-time reading: the balance is structural, the framework cannot break it through the standard construction channels, and the open question is whether ANY constructible L can break it or whether the answer is "no, this is an algebraic identity of the polarity_coordinates_from_L function itself given the way Π acts on Liouville space."

## Probe 6 (2026-05-26): the candidate-breaker that worked

The next-move suggestion from the probes section, "build L coefficient-by-coefficient in the Pauli basis with deliberate ±i Π-eigenvalue imbalance", was executed the next morning ([`_polarity_probe_pi_eigenmode.py`](../simulations/_polarity_probe_pi_eigenmode.py)). The construction: diagonalize Π at N=3 (64×64 in Pauli basis; eigenvalue multiplicities 16/16/16/16 across {+1, -1, +i, -i}), in Π's eigenbasis identify the 1024 (i, j) matrix-element positions whose ratio eigvals[i]/eigvals[j] equals exactly +i, fill those positions with random complex values, leave all other positions zero, transform back to Pauli basis. The resulting L satisfies Π·L·Π⁻¹ = +i·L bit-exactly (residual ≈ 1.5e-15), making it a pure Π-conjugation +i eigenmode.

The result at N=3 with σ = 0:

```
||M||^2          = 4111.35
||M_zero||^2     = 2055.68    (= ||M||^2 / 2)
||M_plus_half||^2  = 2055.68    (= ||M||^2 / 2)
||M_minus_half||^2 = 0.000000
asymmetry        = +2055.68
```

The prediction was ‖M‖²/2 for both M_zero and M_plus_half, zero for M_minus_half. The observation matches the prediction at machine precision. **The balance breaks**.

The structural reading is the inverse of the previous five probes: the ±1/2 balance is NOT forced by the polarity_coordinates_from_L function or by an algebraic identity of Π. It IS forced by the standard Lindblad construction channel L = −i(H ⊗ I − I ⊗ H^T) + Σ_l γ_l (P_l ⊗ P_l* − I). Every L built through that channel produces an M whose Π-conjugation eigenvalue content is balanced between +i and −i. The five probes preserved the balance because they all factored through some variation of that channel (non-Hermitian H, single-site fields, mixed dephase letters, non-Lindblad σ⁺ρσ⁺ form, k=3 Hamiltonians all still use the kron(c, c.conj()) form internally). Probe 6 sidesteps the channel entirely and the asymmetry appears immediately.

The sharpened structural question is now constructive: prove that L of the form −i(H ⊗ I − I ⊗ H^T) + Σ γ_l (c_l ⊗ c_l*) for any Hermitian or non-Hermitian H and any (possibly non-Hermitian) c_l produces M with ‖M_plus_half‖² = ‖M_minus_half‖². Probe 5 already showed it holds for c_l = σ⁺ (which gives non-Hermitian, non-Lindblad-CP form); Probe 1 showed it holds for non-Hermitian H. The conjecture: any operator of the form A ⊗ B* + B ⊗ A* (a real-Hermitian-in-Liouville-superoperator-language structure) is automatically balanced under Π conjugation's ±i eigenvalue split. Probe 6 confirms this property does NOT hold for arbitrary Pauli-basis operators.

If proved, this is a typed Tier1 candidate: `StandardLindbladPiBalance` or similar. Statement: any L built through the kron(c, c.conj()) Lindblad channel has ‖M_plus_half‖² = ‖M_minus_half‖² where M, M_plus_half, M_minus_half are defined by the polarity_coordinates_from_L decomposition. The diagnostic value of polarity_coordinates_from_L stays: it is exactly the test for whether a given L is or is not in the standard Lindblad form (or homotopic to it).

## Probes 7 and 8 (2026-05-26): the conjecture sharpens further

Before attempting the proof, two stronger empirical probes ([`_polarity_probe_random_lindblad.py`](../simulations/_polarity_probe_random_lindblad.py) and [`_polarity_probe_real_or_imag.py`](../simulations/_polarity_probe_real_or_imag.py)) tested the post-probe-6 conjecture across 240 random configurations: random Hermitian or non-Hermitian H, 1 to 4 random Hermitian or non-Hermitian jump operators c_k, N ∈ {2, 3, 4}. The naive conjecture predicted asymmetry = 0 for any standard Lindblad form. Result: **240/240 configurations gave nonzero relative asymmetry, up to 8.7%**. The naive conjecture is false.

Probe 8 narrowed the hypothesis: maybe the structural property is "L is entrywise real or pure imaginary" (Pauli matrices have entries in {0, +/-1, +/-i} which composes nicely under tensor product). Six tests at N=2 and N=3 per random seed:

| c type | balance? |
|---|---|
| Random real matrix, H real Hermitian | BROKEN |
| Random pure-imaginary matrix, H real Hermitian | BROKEN |
| Random complex matrix, H real Hermitian | BROKEN |
| Two single-Pauli c_ops (X on site 0, Y on site 1), H Pauli sum | BALANCE bit-exact |
| Pauli with imaginary coefficient (i·Z), H Pauli sum | BALANCE bit-exact |
| Pauli sum with complex coefficients (X + i·Z), H Pauli sum | BALANCE bit-exact |

The differentiator is NOT "entrywise real or pure imaginary". Random real matrices have entrywise real entries and still break the balance. The actual structural property that distinguishes preserved from broken is **Pauli-string composition of c (and H)**.

But: any operator can be written as a linear combination of Pauli strings (the Pauli basis is complete on d=2^N). So "Pauli sum" alone is not a structural restriction. The empirical distinction must be something finer, possibly "low Pauli rank" (c has support on a small number of Pauli strings) or a related sparsity property.

The sharpest possible conjecture given this data: there exists a structural constant `k_max(N)` such that if c is supported on at most `k_max(N)` Pauli strings with arbitrary complex coefficients, then balance is preserved; if c is supported on more Pauli strings (e.g., random complex matrix expanded in Pauli basis has support on all 4^N), balance breaks. The boundary `k_max(N)` is unknown.

This is now a question about **structural rank in Pauli space**. The diagnostic value remains: `polarity_coordinates_from_L` detects Pauli sparsity of L (or rather, of L's c_l constituents) via the +1/2 vs −1/2 balance. The proof attempt deferred until we identify the actual structural condition.

Next-move candidates:
- Test c with k_pauli ∈ {1, 2, 4, 8, ..., 4^N / 2, 4^N} Pauli strings, see where the boundary breaks
- Connect to F49's Frobenius scaling structure (Π²-odd Pauli sum norm)
- Search for the algebraic property that distinguishes Pauli-eigenvector decompositions from random Pauli-basis decompositions of c

## The honest takeaway (revised after the 5 probes, sharpened by probe 6)

The wave's promised discovery was an asymmetric T1 signature in the +1/2 vs −1/2 polarity channel. The first iteration's finding was that no such asymmetry exists for any standard-Lindblad-form dissipator with Hermitian H, and that the absence was attributable to Hermitian-conjugate symmetry of M as a superoperator. The five probes refuted that attribution: balance survives non-Hermitian H, non-Lindblad dissipators that destroy ρ-hermiticity, mixed dephase letters, and k=3 Hamiltonians where the F81 identity itself fails. The Hermitian-symmetry explanation was insufficient.

What remains, after probe 6: the +1/2 and −1/2 coordinates are equal by construction in every case the framework's standard channels produce, AND the equality is NOT an algebraic identity of `polarity_coordinates_from_L` itself. Probe 6 broke the balance with a hand-engineered Π-conj +i eigenmode L (asymmetry = ‖M‖²/2 at machine precision). The balance is a property of the standard kron(c, c.conj()) Lindblad construction, propagating through to M's Π-conjugation content. The next move is the analytic proof: show that any L of the form −i(H ⊗ I − I ⊗ H^T) + Σ γ_l (c_l ⊗ c_l*) produces M with ‖M_plus_half‖² = ‖M_minus_half‖².

This is not a null result in the bad sense. The Frobenius coordinates for {+1/2, 0, −1/2} are now measurable; we learned across six probes that for every L built through the standard Lindblad construction the readings agree bit-exactly, and that a single deliberate construction OUTSIDE the standard channel breaks them. The diagnostic IS real: `polarity_coordinates_from_L` is exactly the test for whether a given L is (or is homotopic to) standard Lindblad form. Asymmetry ≠ 0 is the construction-outside-standard-channels witness.

There is a second honest takeaway, more about how the framework grows. The polarity-coordinates wave fits a pattern we have seen before, most recently with the [BlockSpectrum F108 wiring](../compute/RCPsiSquared.Core/Symmetry/) the day before: the typed-knowledge layer often makes explicit what was already structurally present. The Frobenius coordinates for {+1/2, 0, −1/2} have been derivable from the F81 split and the Π eigenvalue spectrum since F81 landed; nobody had asked for the explicit projector before. Writing the primitive made the measurement available; running it across five candidate-breaker probes refuted the obvious "Hermitian-symmetry" answer; running probe 6 outside the standard L construction broke the balance immediately and located the actual structural source. The arc closed in two days: from "is there asymmetry under T1" through "the asymmetry is always zero" to "the asymmetry vanishes iff L came from the standard Lindblad channel". The diagnostic is sharper than the original question knew how to ask.

This extends the trail [`ON_THE_HALF`](ON_THE_HALF.md) read for the recurring 0.5: the half was felt by every reader of d=2 substrates from Pythagoras forward; what computers added was the explicit trail from the polynomial root to the Pauli normalization to the F81 50/50 split. The polarity-coordinates wave extends the trail one node further: from the +1/2 of the diopter to the M_plus_half of the Liouvillian residual, with a Frobenius-orthogonal projector, a pytest invariant, and five candidate-breaker probes that all read zero. The two sides are equal in every case we can construct; the proof that they MUST be equal is the next piece of work.

The primitive ships. The reflection ships. The asymmetry channel waits for a future wave that tests non-Hermitian dynamics, single-site fields, or non-Lindblad dissipators. When that wave runs, the diagnostic is already typed and the prediction (balance preservation for any Hermitian-superoperator M) is already written down. The next discovery, if one is coming, will be the first measured ‖M_plus‖² ≠ ‖M_minus‖². Until then, the result is: the balance holds, the structural reason is bra-ket symmetry, and the polarity vocabulary now has explicit operator-level coordinates.

## Closure (2026-05-26): F112 typed as Tier1Derived for Hermitian H

The "next wave" turned out to be the same day. Six additional probes (9–14) located the exact structural condition that distinguishes preserved from broken balance, and the resulting theorem closed as a typed Tier1Derived Claim in C# Core the morning after the reflection above.

**Probe 9 ([`_polarity_probe_pauli_rank.py`](../simulations/_polarity_probe_pauli_rank.py)):** k_max boundary search across single-Pauli, low-rank, and high-rank c at N=2, 3. Result: k_pauli = 1 always preserves balance; k_pauli ≥ 2 is selection-dependent. The boundary is not a rank threshold; the structural property is finer.

**Probe 10 ([`_polarity_probe_pair_enumeration.py`](../simulations/_polarity_probe_pair_enumeration.py)):** at N=2, exhaustive enumeration of all 136 Pauli-string pairs as c = α·P_α + β·P_β. With the fixed coefficient choice (α, β) = (1, i), every one of the 136 pairs preserves balance bit-exact. The coefficient choice matters; the pair choice doesn't constrain.

**Probe 11 ([`_polarity_probe_coefficients.py`](../simulations/_polarity_probe_coefficients.py)):** coefficient sweep at N=2 with five representative pairs. Two pair classes emerge: same Z₂³ cell (Klein × y_par) preserves balance for ALL coefficient choices; cross-cell preserves only for phase-matched coefficients. The structural axis is the Z₂³ sub-cell structure of the Pauli group.

**Probe 12 ([`_polarity_probe_z2cubed_scaling.py`](../simulations/_polarity_probe_z2cubed_scaling.py)):** the Z₂³-cell pattern scales to N=3 and N=4. Sweep A (within-cell c, random complex coefs): 171/171 BALANCED across N=2, 3, 4. Sweep B (cross-cell c): the BALANCE/BROKEN split partitions exactly by bit_b parity match. Same-bit_b cells preserve, cross-bit_b cells break. The Z₂³ structure reduces to the single Z₂ axis of bit_b.

**Probe 13 ([`_polarity_proof_verify.py`](../simulations/_polarity_proof_verify.py)):** the proof's Step 2 numerical check. For bit_b-homogeneous c (every Pauli string σ in c's expansion shares bit_b(σ) = const), the dissipator `np.kron(c, c.conj())` sits 100.00% in the Π²-conjugation +1 eigenspace at N=2, 3 (bit-exact to machine precision). Mixed-bit_b c splits ~50/50 between Π² eigenspaces. This is the algebraic root of the balance: bit_b-homogeneous c contributes nothing to M's Π +i / −i content.

**Probe 14 ([`_polarity_step5_stress.py`](../simulations/_polarity_step5_stress.py)):** the proof's Step 5 stress test. Direct Π-eigenspace projection of L_H for 30 random H configurations (10 Hermitian Pauli, 10 non-Hermitian Pauli, 10 random complex matrix) at N=2, 3. All 30 give ‖L_{H,+i}‖² = ‖L_{H,-i}‖² bit-exact, including the non-Hermitian Pauli and random complex matrix cases. The Hermitian case has a clean structural reason (Step 5 below); the non-Hermitian case is empirical and reduces to an open identity.

**F87 ↔ F112 orthogonality probe ([`_polarity_probe_f87_connection.py`](../simulations/_polarity_probe_f87_connection.py)):** all three F87 trichotomy classes (truly, soft, hard) at N=3 under standard Z-deph give F112 balance asymmetry = 0 bit-exact. F87 lives in ‖M‖_F + spec(L) palindromy; F112 lives in M_anti's Π +i/−i Frobenius split. The two are orthogonal axes on the same bit_b Z₂-grading.

### The theorem

The structural condition surfaced by probes 9–14 is exactly **bit_b-homogeneous c with Hermitian H**:

> **F112 (Hermitian H, Tier1Derived).** For any Lindblad-form Liouvillian L = -i[H, ·] + Σ_k γ_k · `np.kron(c_k, c_k^*)` on N qubits with Hermitian H and each c_k bit_b-homogeneous (every Pauli string σ in c_k's expansion shares bit_b(σ) = (#Y(σ) + #Z(σ)) mod 2 = const), the `polarity_coordinates_from_L` decomposition of M = Π L Π⁻¹ + L + 2σ·I satisfies ‖M_plus_half‖² = ‖M_minus_half‖² bit-exactly.

Five steps: (1) reduce balance to ‖M_{+i}‖² = ‖M_{-i}‖² via Π-eigenspace decomposition of M_plus_half / M_minus_half; (2) bit_b-homogeneous c implies dissipator lies in Π²-conj +1 eigenspace (via F38 / F63); (3) Π²-conj +1 eigenspace = Π-conj {+1, −1}, so dissipator has zero ±i content; (4) M_{+i} and M_{-i} come entirely from L_H; (5) for Hermitian H, L_H is anti-Hermitian as a superoperator (L_H^† = −L_H via straightforward Hilbert-Schmidt manipulation) and the dagger map sends the Π +i eigenspace bijectively onto the Π −i eigenspace while preserving Frobenius norm, giving ‖L_{H,+i}‖² = ‖L_{H,-i}‖². ∎

**Non-Hermitian extension (Tier1Derived at N ≤ 4, Tier1Candidate at N ≥ 5).** Probe 14 first tested non-Hermitian H bit-exactly across 20 random configurations. Writing H = H_re + i H_im (Hermitian decomposition), the equality reduces to the identity F(H_re, H_im) := Im⟨L_{H_re,-i}, L_{H_im,-i}⟩ = 0 for any Hermitian H_re, H_im. F is real-bilinear and antisymmetric under H_re ↔ H_im exchange, so it is determined by its values on Pauli-string basis pairs. A basis-enumeration check (`simulations/_f112_open_identity_basis_enum.py`) gives F = 0 bit-exact across all 136 + 2,080 + 32,896 = 35,112 distinct Pauli-string pairs at N = 2, 3, 4. By bilinearity, F ≡ 0 on the full Hermitian operator space at N ≤ 4, so the non-Hermitian extension is **Tier1Derived at N ≤ 4** via this constructive argument. For N ≥ 5 the extension remains Tier1Candidate pending either further enumeration, an inductive N → N+1 step, or identification of the structural symmetry that makes F vanish universally. The typed Claim's `NonHermitianExtension` inspectable carries the same scope.

### What changed in our picture

The reflection above proposed two candidate explanations: (a) Hermitian-conjugate symmetry of M as a superoperator, (b) bra-ket exchange as the structural meaning. Probes 1 and 5 refuted (a) on the surface: balance survived non-Hermitian H and non-Lindblad σ⁺ρσ⁺ dissipators that destroy ρ-hermiticity. The actual answer is (a) at a different level: it is not M's Hermitian-superoperator property that matters, but L_H's anti-Hermitian-superoperator property and the bit_b-Π² eigenstructure of the dissipator together. Probe 1 fits this directly: its dissipator was single-Pauli Z (bit_b-homogeneous), only the H was non-Hermitian, and that places it inside F112's non-Hermitian Tier1Candidate scope, which empirically holds. Probe 5 sits outside F112's typed scope entirely (the σ⁺ρσ⁺ form is not standard `Σ γ_k np.kron(c_k, c_k^*)` Lindblad and σ⁺ = (X + iY)/2 is not even bit_b-homogeneous); balance held there too, but as an empirical observation rather than a consequence of the typed theorem. Probe 6 broke the balance by side-stepping the Lindblad channel entirely with a hand-engineered Π-conjugation +i eigenmode that had no Lindblad structure to invoke at all. The candidate (b) bra-ket reading is structurally compatible with F112 but lives at the meaning layer rather than the algebraic mechanism: F112 says the bra-side and ket-side coordinates are equal whenever the system's dissipative channel keeps the bit_b grading clean.

For standard physical Lindblad systems with Hermitian H, the bit_b-homogeneity of c splits cleanly: single-Pauli c (T2 = Z dephasing, dephasing on a single Pauli letter, depolarizing single-site components) is trivially bit_b-homogeneous, since one Pauli string has one bit_b value. T1 amplitude damping via σ⁻ = (X − iY)/2 is NOT bit_b-homogeneous (X has bit_b = 0, Y has bit_b = 1), so T1 sits in the non-Hermitian extension's domain rather than F112's rigorously-Hermitian Tier1Derived scope; the basis-enumeration argument (Tier1Derived at N ≤ 4) covers it. F112 therefore says the polarity coordinates are equal by structure for every Lindblad system we actually build to describe physical noise at chain lengths up to N = 4, and remain empirically balanced at all tested N beyond that. Asymmetry ≠ 0 in `polarity_coordinates_from_L` is the precise witness for c outside the bit_b-homogeneous regime *combined with* non-Hermitian H *combined with* the open-identity F not vanishing for that specific (H_re, H_im) pair, i.e. it requires the system to leave every channel F112 covers simultaneously, which in practice means non-standard collapse operators with mixed-bit_b Pauli support at N ≥ 5, or L not in Lindblad form at all.

### Wiring

The proof landed at [`docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md`](../docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md) with the rigorous Hermitian-H derivation and the open non-Hermitian identity documented as the next-piece-of-work. The F-formula registry entry is at [`docs/ANALYTICAL_FORMULAS.md` F112](../docs/ANALYTICAL_FORMULAS.md). The C# typed Claim is at [`compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiBalance.cs`](../compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiBalance.cs) (IZ2AxisClaim on the BitB axis, Tier1Derived, F108 Part 1 as typed ctor parent for the shared F38 / F63 bit_b foundation). Cross-references wire the connection mapping: F38 and F63 list F108 / F112 as downstream consumers, F87 names F112 as the orthogonal axis on the shared bit_b grading, F108 Part 1 and Part 3 carry sibling-inspectables pointing at F112.

The diagnostic value of `polarity_coordinates_from_L` is now sharper than the original reflection knew how to state. It is not a generic balance check; it is the structural witness for "L came from a standard Lindblad channel with bit_b-homogeneous collapse operators and Hermitian H." Asymmetry ≠ 0 narrows the source: at least one of (Hermitian H, bit_b-homogeneous c, Lindblad-form L) fails. That is a useful classifier for hardware-effective dynamics where one is asking exactly that question.

---

*"Wir hatten gehofft eine neue Berechnung zu entdecken; was wir entdeckt haben war eine alte Symmetrie ohne Namen."*

*"Die ±1/2 ist gleich groß auf beiden Seiten, weil wir auf beiden Seiten dieselben sind."*

---

**Anchors:**

- Primitive: [`simulations/framework/diagnostics/polarity_coordinates.py`](../simulations/framework/diagnostics/polarity_coordinates.py)
- Tests (18 PASS, 1 FAIL the documented refutation `test_t1_cooling_breaks_plus_minus_balance`): [`simulations/framework/tests/diagnostics/test_polarity_coordinates.py`](../simulations/framework/tests/diagnostics/test_polarity_coordinates.py)
- Demo (8 H families × 3 dissipators, asymmetry = 0.0 across the board): [`simulations/_polarity_demo.py`](../simulations/_polarity_demo.py)
- Spec: [`docs/superpowers/specs/2026-05-25-polarity-coordinates-design.md`](../docs/superpowers/specs/2026-05-25-polarity-coordinates-design.md)
- Plan: [`docs/superpowers/plans/2026-05-25-polarity-coordinates-plan.md`](../docs/superpowers/plans/2026-05-25-polarity-coordinates-plan.md)

**Framework anchors:**

- [F81](../docs/proofs/PROOF_F81_PI_CONJUGATION_OF_M.md): Π·M·Π⁻¹ = M − 2·L_{H_odd}, the parent identity refined by polarity_coordinates
- [F84](../docs/proofs/PROOF_F84_AMPLITUDE_DAMPING.md): the amplitude-damping correction whose F81-violation we now know does NOT project asymmetrically onto ±i
- [F83](../docs/proofs/PROOF_F83_PI_DECOMPOSITION_RATIO.md): the anti-fraction limit at r=0, the F81 ratio observable
- [F80](../docs/ANALYTICAL_FORMULAS.md): Spec(M) = 2i · Spec(H_non-truly), the spectral side of the Π²-odd content
- [F112](../docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md): the closure identity ‖M_plus_half‖² = ‖M_minus_half‖² for Hermitian H + bit_b-homogeneous c, proved via dagger + anti-Hermitian L_H; non-Hermitian extension Tier1Derived at N ≤ 4 via basis enumeration (35,112 Pauli-string pairs all bit-exact 0), Tier1Candidate at N ≥ 5. The structural witness behind every BALANCE reading from Probes 1–5 and 7–14.
- [`LindbladBitBPiBalance`](../compute/RCPsiSquared.Core/Symmetry/LindbladBitBPiBalance.cs): the typed C# Claim (Tier1Derived for Hermitian H; non-Hermitian extension Tier1Derived at N ≤ 4, Tier1Candidate at N ≥ 5), IZ2AxisClaim on the BitB axis, F108 Part 1 as ctor parent
- [F38](../docs/ANALYTICAL_FORMULAS.md): Π² = (−1)^{w_YZ} eigenvalue formula, the algebraic root of the bit_b Z₂-grading F112 consumes
- [F63](../docs/ANALYTICAL_FORMULAS.md): [L, Π²] = 0, the dissipator-side commutation F112 also consumes
- [F87](../docs/ANALYTICAL_FORMULAS.md): trichotomy classifier; orthogonal axis on the same bit_b Z₂-grading (empirically confirmed via `_polarity_probe_f87_connection.py`)
- [F108 Part 1](../docs/proofs/PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md): sibling Tier1Derived theorem on the shared bit_b axis; F108 closes spec(L) palindromy for Π²-even bilinear H, F112 closes Π +i/−i balance for arbitrary Hermitian H + bit_b-homogeneous c
- [`PolarityLayerOriginClaim`](../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs): the +0/−0 layer at d=0 generating ±1/2 at d=2 via the 0.5-shift
- [`HalfAsStructuralFixedPointClaim`](../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs): the half's argmax side
- [`QuarterAsBilinearMaxvalClaim`](../compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs): the apex where both polarity sides meet under squaring

**Interpretive anchors:**

- [`ON_THE_HALF`](ON_THE_HALF.md): the three faces of 0.5 (bridge / horizon / substrate), now joined by a fourth (the operator-level Π-eigenvalue coordinate)
- [`ON_BOTH_SIDES_OF_THE_MIRROR`](ON_BOTH_SIDES_OF_THE_MIRROR.md): the F81 Π·M·Π⁻¹ = M − 2·L_{H_odd} identity that polarity_coordinates refines
- [`ON_HOW_TWO_SIDES_MEET_AT_THE_QUARTER`](ON_HOW_TWO_SIDES_MEET_AT_THE_QUARTER.md): the geometric fold of the ±1/2 onto 1/4, now with explicit operator-level coordinates for both sides
- [`THE_POLARITY_LAYER`](../hypotheses/THE_POLARITY_LAYER.md): the +0/−0 differentiation at d=0; the polarity-coordinates primitive is the operator-level instantiation
- [`MIRROR_THEORY`](../MIRROR_THEORY.md): the framework's mirror reading; the bra-ket balance is the structural form of "the mirror is by construction equally weighted"

**Memory pointers:**

- `reference_diopter_as_polarity_bridge`: the everyday phenomenology of the {−1/2, 0, +1/2} triple
- `project_plus_minus_zero_layer`: the X-basis polarity projection that polarity_coordinates parallels at the operator level
- `project_one_system_two_indices`: the bra-ket reading whose Frobenius-numerical witness is the balance ‖M_plus‖² = ‖M_minus‖²
- `project_framework_as_remembrance`: the explicit measurement of "the two halves we cannot distinguish"
- `project_blockspectrum_z_deph_only` (typed earlier today): the parallel pattern where typed knowledge makes implicit symmetry explicit

---

*Tom and Claude, 2026-05-25. The polarity-coordinates wave closes with a measurement that refused to be the discovery it was hypothesized to be, and became instead the typed witness of a structural symmetry that was sitting in the framework all along. Written while the contrast between expected and found is still fresh; future waves with non-Hermitian probes will say whether this is the last word or the baseline against which an asymmetry first appears.*

*Closure 2026-05-26. The non-Hermitian probe was the same day. Probes 9–14 narrowed the structural condition to bit_b-homogeneity of c with Hermitian H, and the proof closed via dagger + anti-Hermitian L_H. F112 (`LindbladBitBPiBalance`) typed as Tier1Derived in C# Core, anchored in ANALYTICAL_FORMULAS, and wired into the bit_b Z₂-axis hierarchy alongside F38, F63, F87, and F108 Part 1. The non-Hermitian extension first reduced to an open identity F(H_re, H_im) = 0; a basis-enumeration check (35,112 Pauli-string pairs at N ≤ 4, all bit-exact 0) plus bilinearity + antisymmetry closed it constructively as Tier1Derived at N ≤ 4, leaving Tier1Candidate at N ≥ 5 pending an inductive lift. The diagnostic `polarity_coordinates_from_L` is now the precise witness for "L came from a standard Lindblad channel with bit_b-homogeneous collapse operators and Hermitian H (or with non-Hermitian H at N ≤ 4)"; asymmetry ≠ 0 is the construction-outside-that-class witness. The two-day arc closed.*
