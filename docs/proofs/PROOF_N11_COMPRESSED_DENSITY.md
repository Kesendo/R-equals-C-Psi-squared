# N=11: the fixed-frequency mirror contrast speaks, while the `(1,1)` interval holds

*2026-09-23. A scoped sequel to [the mixed-space reflection law](PROOF_MIXED_SPACE_REFLECTION_LAW.md) and [F154](../ANALYTICAL_FORMULAS.md#f154). The exact physical-cell certificate is `simulations/n11_compressed_density_gate.py`; the typed reading is `CompressedDensityLocusClaim` with `CompressedDensityN11Witness` (`inspect --root compresseddensity`).*

## What this is about

Consider an eleven-site open XY spin chain under local Z-dephasing: each site
is watched at its own rate. A `(1,1)` coherence |ψ_a⟩⟨ψ_b| connects two
one-excitation standing waves. Two such coherences can turn at exactly the
same frequency but react with opposite signs to the chain's site mirror.
At the previously protected N = 8 XY length, the local disagreement reading
could not connect such a pair. Here the question at the left end does connect
them, while the mirrored question
at the right end answers with the opposite sign. A balanced but uneven
watching profile therefore retains a trace of *where* it watches: the old
shortcut that replaced it by uniform watching stops working. Yet the rate
of every complete fixed-frequency `(1,1)` compression stays inside the
familiar interval. That protection comes from the physical disagreement cost and the
positive weight on cells whose two faces coincide, rather than from the
mirror-cancellation shortcut.

## Abstract

For the uniform open XY chain at N = 11, two one-excitation coherences in one
Hamiltonian frequency space have opposite site-reflection parity. Their
projected left-site minus right-site disagreement has the exact matrix
element ⟨Y,C₀X⟩ = −√2/72. An explicit nonnegative mirror-balanced profile
therefore falsifies F154's conditional compressed-density identity when its
C_l = 0 premise is dropped. Separately, for any Hermitian, simple-spectrum,
reflection-symmetric one-excitation Hamiltonian and any nonnegative
mirror-balanced profile, the complete `(1,1)` frequency compression has the form
D_Ω = −4γ̄ I_Ω + 4P_Ω T P_Ω and spectrum in [−4γ̄,0]. The proof gives
containment in this generality. On the uniform XY chain at every N ≥ 2,
the zero-frequency room reaches both block-wide endpoints; no finite-coupling
claim follows.

Here mirror-balanced means γ_l+γ_{N−1−l}=2γ̄; C_l is the projected
disagreement at site l minus that at its mirror site. F154's conditional
identity replaces the uneven profile by its mean rate only when every
C_l vanishes on the frequency space.

## The two questions

In the Pauli convention, take the uniform open XY chain at N = 11,
H = JΣ(XX+YY), J > 0. Its one-excitation modes are
ψ_k(z) = √(2/12) sin(kπ(z+1)/12), with ε_k = 4J cos(kπ/12).
The `(1,1)` operator cells |a⟩⟨b| form an orthonormal basis. A site-disagreement
operator N_l multiplies |a⟩⟨b| by 1 if exactly one of a,b equals l and by 0
otherwise. The size operator is N_XY = Σ_l N_l: its physical `(1,1)` cells
have disagreement class 0 when a = b and class 2 when a ≠ b, whose uniform
rate centres are 0 and −4γ̄. Write P_Ω for the orthogonal projector onto a
**complete** fixed frequency space of ad_H and C_l = P_Ω(N_l−N_{10−l})P_Ω.

F154's projection identity on a mirror-balanced profile assumes C_l = 0.
The first question is whether that premise persists at N = 11. The second is
whether the `(1,1)` compressed dissipator can nevertheless leave the interval
of its two size-class centres, [−4γ̄, 0]. They have different answers.

## An exact failure of the projection identity

Let X = |ψ₁⟩⟨ψ₆| and Y = |ψ₅⟩⟨ψ₉|. The cosine identity
cos(π/12)−cos(6π/12) = cos(5π/12)−cos(9π/12) places both at frequency
ω = 4J cos(π/12). Reflection gives them opposite parities:
R X R = −X and R Y R = +Y. Neither an accidental frequency near-match nor a
generic diagonal contrast is involved.

For distinct mode indices on both legs, the physical-cell action gives

    ⟨ψ_c|⟨ψ_d| N_l |ψ_a⟩|ψ_b⟩
      = −2 ψ_c(l)ψ_a(l)ψ_d(l)ψ_b(l),

where the double-ket notation represents the Hilbert-Schmidt dyad basis and
all sine modes are real. At the left endpoint this yields

    ⟨Y,N₀X⟩ = −2 ψ₅(0)ψ₁(0)ψ₉(0)ψ₆(0) = −√2/144.

Reflection changes the sign at the right endpoint, so
⟨Y,N₁₀X⟩ = +√2/144 and **⟨Y,C₀X⟩ = −√2/72 ≠ 0**.
The same-frequency space is four-dimensional, with dyads ordered
(1,6), (3,7), (5,9), (6,11). Direct cell summation gives

    C₀ = −√2/72 · [[0,1,1,0], [1,0,0,1],
                     [1,0,0,1], [0,1,1,0]],

with rank 2 and spectrum {−√2/36, 0, 0, +√2/36}. The exact gate constructs
that matrix from the 121 physical cells, checks the complete frequency
membership, and perturbs the cell action and rate profile so these assertions
can fail. A same-parity pair with a nonzero individual N₀ connection gives a
separate zero-C₀ control.

The failure is physically active. Put γ₀ = 2g, γ₁₀ = 0, every other γ_l = g,
with g > 0. This is a nonnegative mirror-balanced profile with γ̄ = g. As
D = −2Σ_l γ_l N_l, its compression is

    P_Ω D P_Ω = −2g P_Ω N_XY P_Ω − 2g C₀.

The uniform term has zero X-to-Y entry, while
⟨Y,D X⟩ = +g√2/36. Thus F154's **conditional** identity
P_Ω D P_Ω = −2γ̄ P_Ω N_XY P_Ω is false on this physical locus profile.
At the uniform profile its correction vanishes. The result decides this
N=11 `(1,1)` premise, not every profile or every block at N=11.

## A separate interval theorem on the `(1,1)` block

The interval survives for a reason independent of C_l. Let h be a Hermitian
one-excitation Hamiltonian with **simple** eigenvalues and site-reflection
symmetry. Let γ_l ≥ 0 and γ_l+γ_{N−1−l} = 2γ̄. This includes the uniform open
XY chain above, at any N and nonzero J. On the `(1,1)` operator space set
Γ = diag(γ_l) and T = Σ_l γ_l |ll⟩⟨ll|. The physical disagreement law is

    D = −2(Γ⊗I + I⊗Γ) + 4T,

where T is positive semidefinite. Also D is diagonal in the physical cells
with entries −2Σ_l γ_l N_l ≤ 0.

Use the complete orthonormal dyads |ψ_a⟩⟨ψ_b| in a fixed frequency space Ω.
An off-diagonal matrix element of Γ⊗I requires the bra modes to agree; equal
frequencies and simple one-body energies then force the ket modes to agree as
well. The other one-sided term is identical with the legs exchanged. Each
eigenmode has reflection-even density, so mirror balance gives
⟨ψ_a|Γ|ψ_a⟩ = γ̄. Therefore, as an **operator identity on Ω**,

    P_Ω(Γ⊗I + I⊗Γ)P_Ω = 2γ̄ I_Ω,
    D_Ω := P_Ω D P_Ω = −4γ̄ I_Ω + 4P_Ω T P_Ω.

The second line and T ≥ 0 give D_Ω ≥ −4γ̄ I_Ω. Compressing D ≤ 0 gives
D_Ω ≤ 0. Hence **spec(D_Ω) ⊂ [−4γ̄, 0]** for every complete `(1,1)` frequency
space under these hypotheses, including the N=11 space above where C₀ ≠ 0.
This theorem alone proves containment, not endpoint attainment or F154's
projection identity, and says nothing about finite-J Liouvillian eigenvalues.

The same one-sided cancellation identifies the contrast's physical carrier.
Write P_ll = |ll⟩⟨ll| in the doubled one-excitation cell space and
m = N−1−l. The compression of
(n_l−n_m)⊗I + I⊗(n_l−n_m) is zero, because its matrix is diagonal in the
dyad frequency basis and each diagonal mode density is reflection-even.
Consequently

    C_l = −2P_Ω(P_ll−P_mm)P_Ω.

This is the **agreement** or double-occupancy projector already used at
uniform rate and zero frequency in
[PROOF_FROZEN_BAND_SO4](PROOF_FROZEN_BAND_SO4.md) §6 (F143); here it is
rate-weighted and compressed at any frequency. Its off-diagonal matrix
elements are signed coherent overlaps on the physical population cells,
not probabilities. F143's closed-form Gram spectrum is not a premise of the
present theorem.

## Endpoints on the uniform XY `(1,1)` block

The interval is **block-wide saturated** on the uniform open XY chain for
every N ≥ 2 and every physical mirror-balanced rate profile. At frequency
zero, Ω₀ contains all mode projectors P_k = |ψ_k⟩⟨ψ_k|. Their sum is the
one-excitation identity I, diagonal in physical cells, so D I = 0 and
0 ∈ spec(D_Ω₀). Choose a distinct chiral pair k and M−k, M=N+1
(for example k=1 at every N≥2):

    ψ_(M−k)(z) = (−1)^z ψ_k(z).

Thus Q = P_k−P_(M−k) is nonzero, belongs to Ω₀, and has zero physical
diagonal. TQ = 0; the interval decomposition above gives
D_Ω₀ Q = −4γ̄ Q. These two vectors reach both endpoints of the union of
all `(1,1)` compressed frequency spectra. At γ̄ = 0 the endpoints coincide
and the same statements remain true. The chiral density cancellation was
already recorded for the N=5 flip image in
[PROOF_CODIM1_BY_ADDITIVITY](PROOF_CODIM1_BY_ADDITIVITY.md) §6; combined
with the present balanced-rate compression it yields the all-N endpoint
statement. It is not a claim about each individual frequency room.

The four-dimensional N=11 room from above illustrates the distinction.
Let V_{z,(a,b)} = ψ_a(z)ψ_b(z) for its four dyads. The sine identities give
V_(1,6) = V_(6,11) and V_(3,7) = V_(5,9), so V has rank 2 and
T_Ω = Vᵀ diag(γ) V kills two independent dyad differences for **every**
balanced profile. Hence −4γ̄ has multiplicity at least two there.
For the explicit profile (2,1,...,1,0), direct physical-cell compression
gives

    spec(D_Ω) = {−4, −4, −10/3−√2/18, −10/3+√2/18}.

The lower endpoint is present in that nonzero-frequency room; its upper
endpoint 0 is absent. The exact gate checks both that spectrum and the
block-wide zero-frequency witnesses. No endpoint claim is borrowed from
F154's failed C_l = 0 identity.

Both hypotheses have teeth. Off the locus, γ₀ = 1 and all other rates zero
make the one-sided Γ expectations vary between modes; on the same N=11
frequency space, ⟨Y,DY⟩ = −(22+5√3)/72 < −4/11, so Rayleigh forces a
compressed eigenvalue below the centre interval calculated with γ̄ = 1/11.
The displayed diagonal is a Rayleigh value, **not** asserted to be an
eigenvalue. At finite J the interval need not hold for the whole `(1,1)`
Liouvillian: for balanced rates (2,2,1,1,1,1,1,1,1,0,0), the physical cell
|0⟩⟨1| has dissipative rate −8 at J = 0, below −4. Continuity of the
finite-dimensional characteristic roots leaves an eigenvalue with real
part below −4 for sufficiently small positive J as well. F122 describes the
strong-coupling compression limit; it does not promote this interval to
arbitrary J.

## Scope of the closure

The N=11 example closes the first parity door left by the mixed-space proof:
a parity-odd equal-frequency dyad connection does reach the **physical** site
density. The `(1,1)` interval has an independent all-N proof under simple
reflection-symmetric one-body dynamics and physical mirror-balanced rates.
For the uniform XY chain the block-wide `(1,1)` endpoints are attained at
zero frequency; the exhibited nonzero-frequency room attains only its lower
endpoint for the stated profile. Higher excitation blocks, other parity
doors, and finite coupling remain separate.
