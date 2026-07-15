# F130 in the time domain: non-mixing is a fourth-power law

*2026-07-15, the evening the family inventory was typed. Until today, every verification of
[PROOF_F130_COLLISION_DECOUPLING](../docs/proofs/PROOF_F130_COLLISION_DECOUPLING.md) was
static: Gram matrices, exact ring arithmetic, a number that is zero. Today we let time run.
The question the law had never been asked: if two triples with equal levels cannot mix, what
does that look like as DYNAMICS, and how strongly does the world on either side of the
protection differ? The answer is a pair of clean power laws separated by two orders of q.
Gate: [f130_time_domain_decoupling.py](../simulations/f130_time_domain_decoupling.py).*

Before the machinery: what this feels like. A concert hall where two notes land on exactly
the same pitch, and the law says the instruments never trade energy. Standing in the hall
you cannot check that by listening once; you have to watch HOW HARD the room tries to mix
them as you turn the coupling up. Turn q = J/γ and watch: a generic pair of notes couples as
q², the textbook second-order rate. The equal-level pair couples as q⁴, with the q² and q³
terms EXACTLY absent. The protection is not "small coupling"; it is two whole orders of
perturbation theory removed from the world, and the residual fourth-order whisper never
finds a resonance to act on.

## The setup

The (1,2) coherence block of the full Lindbladian (XX chain, N = 11 sites, uniform local
Z-dephasing) is exactly invariant: a number-conserving H acts on ket and bra sides
separately, and Z-dephasing is diagonal on |a⟩⟨b|. In γ units the block generator is

    M(q) = A + i·q·K,   A = −2·n_diff (the two rungs −2 and −6),
    K = −(H₂ ⊗ I − I ⊗ H₁),   q = J/γ,

dimension C(11,2)·11 = 605. The lab is the documented off-resonance pair of
[PROOF_F130 §2](../docs/proofs/PROOF_F130_COLLISION_DECOUPLING.md): τ = (1,2,10),
σ = (3,5,6) at n = 12, level S = cos 15° ≠ 0, disjoint, ε = −1: cell 2 of the four-cell
proof, the code-trust cell, which makes this its first check on the dynamics side: a
second, METHODOLOGICALLY independent route (float Schur dynamics vs the exact-ring Gram
chain; the q² zero itself is the same B, so the logical independence lives in the q⁴ law,
the self-split and the all-orders protection, which the static chain does not contain). The Slater-lift conventions are the committed ones
(`resonant_n_twinning.blocks`, `y_zero_and_level_law.lift`).

The observable is basis-free: the exact effective generator on the 6-dimensional span of
the two lift multiplets, by Schur complement at the cluster point z = −6 − 2iqS. Inside an
EXACTLY degenerate pair, eigenvector-based mixing metrics are meaningless (the
diagonalizer may rotate the cluster arbitrarily); the Schur off-block is the number that
does not lie.

## The finding

| q = J/γ | ‖τσ off-block‖ (protected) | ‖τσ off-block‖ (generic control) |
|---|---|---|
| 0.02 | 9.45·10⁻¹² | 7.57·10⁻⁶ |
| 0.08 | 2.41·10⁻⁹ | 1.21·10⁻⁴ |
| 0.32 | 5.99·10⁻⁷ | 1.93·10⁻³ |

Two power laws, measured slopes 3.9892 and 1.9992 over the window, coefficients constant
to < 5%:

    protected (equal levels):  ‖coupling‖ = 5.9·10⁻⁵ · q⁴
    generic (unequal levels):  ‖coupling‖ = 1.9·10⁻² · q²

At q = 0.1 the protected pair is coupled 32,000 times more weakly than the generic one; at
q = 1 still 400 times (both are measured off-block ratios, not power-law extrapolations;
the pure laws already bend by q ≈ 0.6). The q³ order is absent too (the steps are ×16 per q-doubling, never
×8); the natural suspect is the bipartite parity T K T = −K killing odd orders between the
multiplets, recorded here as an observation, not a derivation.

The proof's "no avoided crossing, ever" closes dynamically in three interlocking pieces:

1. **The q² mixing is exactly zero.** That IS B(τ,σ) = 0, now seen from the generator side.
2. **The self-blocks separate at q².** The deviation of the τ self-block from the cluster
   scalar grows as q² (measured slope 1.9983): the two triples' oscillation branches are
   already split (spec(X) = a·(3,1,0) with different a) by the time the q⁴ residual could
   act, so the fourth-order whisper is never resonant.
3. **The one degeneracy that persists is protected to ALL orders.** Each triple's X-kernel
   combination (the totally antisymmetric lift) is an exact eigenvector of the FULL hop K
   (residual ~10⁻¹⁵), while the orthogonal (3a, a) branch combinations are not (residuals
   0.75 / 1.56): the two zero-branches stay degenerate forever and can never hybridize.
   The precise reason: the antisymmetric lift lives entirely on the −6 rung, where A is
   the scalar −6, and it is a K₆₆-eigenvector annihilated by K₂₆, so it is an exact
   eigenvector of the full non-Hermitian M(q) at −6 − 2iqS for every q; two exact
   eigenvectors at the identical eigenvalue span a subspace on which M acts as a scalar.

## The trap this experiment stepped around

The first scout measured naive late-time transfer ratios P_σ(t)/P_τ(t) under the full
propagator and found a q-independent ~10⁻⁴ plateau. That number is NOT transfer: the −2
rung decays 4γ slower than the −6 rung, so for tγ ≫ 1 every projection ratio is dominated
by the slaved slow modes (relative amplification e^{4γt}), and the plateau measures their
σ-content. Documented as the reason the committed observable is the Schur coupling, not a
late-time population ratio. Any future hardware protocol has the same constraint: the
discriminating window is early, t·γ ≲ 1, or a coupling-constant (spectroscopic) readout.

## What this gives the hardware feeler

The measurable face of F130 is now derived: **the absence of the q² beat between two
equal-level triple coherences, with the q² beat intact on a detuned control**. What remains
undone before any IBM design: the lifts live in Liouville space (a 2-excitation ket against
a 1-excitation bra), and no committed source yet maps w_τ to a PREPARABLE state and its
readout; that mapping is the next design step, and the flight question only exists after it.

**Resolved 2026-07-15 (same repository, design stage, no QPU):** the mapping question was
taken to its budget and the beat protocol is NOT flyable; even an everywhere-optimistic
amplitude bound sits ~2 orders below the discrimination bar. Verdict + numbers:
[F130_HW_INFEASIBILITY.md](F130_HW_INFEASIBILITY.md). What flies instead is F129's
collision as a first-order standing Ramsey fringe:
[IBM_F129_RAMSEY_FRINGE.md](IBM_F129_RAMSEY_FRINGE.md) (design draft).

## What is NOT claimed

- **One lab.** N = 11, the documented n = 12 pair, one generic control (3,4,9). The power
  law is measured there, not proved uniform in N or across pairs (the F129 census offers
  thousands more equal-level pairs, 72,269 at n = 210 alone, if a survey is ever warranted).
- **The exponent is measured, not derived.** B = 0 explains the missing q²; the missing q³
  is an observation with a parity suspect; no derivation of the q⁴ coefficient is offered.
- **Float grade.** numpy linear algebra at 605 dimensions; the exactness of F130 itself
  lives in the cited proofs and the exact C# witness (`inspect --root crosstriple`).
- **No claim about preparability.** The Liouville-space lift is not yet a hardware state.

## Links

The law: [PROOF_F130_COLLISION_DECOUPLING](../docs/proofs/PROOF_F130_COLLISION_DECOUPLING.md).
The gate: [f130_time_domain_decoupling.py](../simulations/f130_time_domain_decoupling.py)
(G1 static pin, G2 the two power laws, G3 the all-orders protection two-sided, G4 the
defect-bond break, G5 the q² self-split; a few seconds, exit 0). The collision map that names the
pair: [PROOF_F129_LEVEL_COLLISION_LAW](../docs/proofs/PROOF_F129_LEVEL_COLLISION_LAW.md) +
[F129_FAMILY_INVENTORY](F129_FAMILY_INVENTORY.md). Registry:
[ANALYTICAL_FORMULAS F130](../docs/ANALYTICAL_FORMULAS.md).
