using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The Klein-V₄ subgroup of unitary involutions on the 4^N Pauli basis that
/// realizes the dephase-letter Klein V₄ {I, Z↔Y, Z↔X, Y↔X} on the operator-space
/// F1 palindrome family {Π_Z, Π_X, Π_Y}.
///
/// <para><b>Welle 12 closure (2026-05-27).</b> Surfaced during the Welle 10d audit when
/// the F112 sparse-rep refactor uncovered the basis-convention twist D · L_natural · D
/// in the standard codebase pipeline; the diagonal D turned out to BE the operator-
/// space lift of the Z↔Y dephase-letter swap, not just a basis change. Welle 12 Tasks 1
/// and 2 closed the universal-N structural proofs.</para>
///
/// <para><b>Three involutions on operator space</b> (all order-2 unitary, all built as
/// per-site Kronecker products of 4×4 single-site matrices on basis (I, X, Z, Y),
/// matching the <see cref="PauliLetter"/> packing convention <c>i = a + 2·b</c>):</para>
/// <list type="bullet">
///   <item><b>D</b> = ⊗_l d_l with d_l = diag(1, 1, 1, −1) on (I, X, Z, Y). Diagonal
///         entry at flat index k equals (−1)^n_Y(k) where n_Y(k) counts Y letters in
///         the k-th Pauli string. Realizes Z↔Y: D · Π_Z · D = Π_Y. See
///         <c>docs/proofs/PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md</c>.</item>
///   <item><b>H</b> = ⊗_l h with h the 4×4 basis-permutation matrix that fixes the I
///         and Y basis vectors and swaps the X and Z basis indices per site. H is a
///         pure permutation (signs all +1), self-inverse. Realizes Y↔X: H · Π_Y · H
///         = Π_X (and H · Π_X · H = Π_Y). NOTE: h does NOT swap the physical Pauli
///         operators σ_X ↔ σ_Z; it permutes their basis indices in the Pauli-string
///         enumeration (see "Convention note" in
///         <c>PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md</c>).</item>
///   <item><b>Q_zx</b> = H · D = ⊗_l (h · d_l). Order-2 involution. Realizes Z↔X:
///         Q_zx · Π_Z · Q_zx = Π_X.</item>
/// </list>
///
/// <para><b>Klein-V₄ closure</b>: D · Q_zx · H = I (identity); D, Q_zx, H commute
/// pairwise. The three operators together with I form a Klein four-group ≃ Z₂ × Z₂
/// in U(4^N) for any N, faithfully representing the V₄ on dephase letters.</para>
///
/// <para><b>Universal N via tensor-product factorization</b>: each identity reduces
/// via the Kronecker mixed-product property to a 4×4 per-site check. See
/// <c>docs/proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md</c> for the full
/// structural proof.</para>
///
/// <para><b>Significance</b>: F1-family identities transfer between dephase
/// letters by two distinct routes, each enabled by a different subset of the
/// Klein-V₄ generators.</para>
///
/// <para><b>Route 1 (per-axis structural re-run)</b>: identities whose proof
/// reduces to (a) the F38 Π_d² eigenvalue formula on Pauli strings and (b)
/// Pauli-basis matrix-support structure transfer between all three dephase
/// letters by direct re-derivation, substituting axis_d for the relevant Z₂
/// axis (bit_b for d ∈ {Y, Z}, bit_a for d = X). The Welle 11 F112
/// non-Hermitian extension is the prototype; the Welle 13 closure of F112-X
/// and F112-Y (docs/proofs/PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md) gives
/// the full cross-dephase F112 family this way. F108 Parts 1, 2, 3 are
/// Klein-V₄ equivalent per Welle 14
/// (docs/proofs/PROOF_F108_KLEIN_V4_EQUIVALENCE.md): Part 3 = Part 1 via D
/// (operator-space, bit_b-axis-preserving) and Part 2 = Part 1 via Hadamard
/// transport (Hilbert-space, bit_a ↔ bit_b axis swap). On Π_5bilinear the
/// operator-space Klein-V₄ assignment is PERMUTED relative to canonical Π_d:
/// the {I, D, H} subgroup acts (D: Z↔Y, and H: Z↔X via
/// H · Π_5b(Z) · H = Π_5b(X) bit-exact, closed 2026-05-29 in
/// simulations/f108_bita_d_search.py), while Q_zx (= H·D) is the one that
/// leaves the canonical Π_5b set. The canonical-Π pairings fail (Q_zx on Z↔X,
/// H on Y↔X), which is why the original Welle-14 sweep read it as the {I, D}
/// subgroup alone. L-transport (the Lindbladian, not the operator) still needs
/// the Hilbert-space Hadamard, carried only by Q_zx (Route 2 mechanism).</para>
///
/// <para><b>Route 2 (Hadamard transport via Q_zx)</b>: Q_zx is the
/// operator-space lift of the Hilbert-space unitary U_H^⊗N (Hadamard ⊗N),
/// so Q_zx-conjugation maps Lindblad-form L_Z to Lindblad-form L_X with
/// rotated (H, c_k). This gives F112-Z → F112-X transport explicitly via
/// a Lindblad-form-preserving similarity. The D and Q_yx (= H) involutions
/// are operator-space-only: they do NOT lift to Hilbert-space unitaries
/// (D-conjugation would require V such that V·Y·V⁻¹ = −Y, V·X·V⁻¹ = X,
/// V·Z·V⁻¹ = Z, which is impossible by Pauli algebra). So D and H
/// intertwine the Π_d operators but do NOT transport L between dephase
/// representations. F112-Y cannot be obtained "for free" from F112-Z via
/// D-conjugation; it requires Route 1 (direct bit_b axis re-run with d = Y
/// phase).</para>
///
/// <para>Bottom line: Klein-V₄ equivariance is genuine for Route 1
/// (universal across all three letters) and partial for Route 2 (only the
/// Hadamard {I, Q_zx} subgroup transports L). F1-family identities
/// depending only on Π_d / Π_d² eigenvalue structure and Pauli-support
/// disjointness transfer freely between all three dephase letters;
/// identities depending on the Lindblad-form L itself transfer only
/// between (Z, X) via Q_zx.</para>
///
/// <para>Tier1Derived universal N per the two PROOF docs. No <see cref="IZ2AxisClaim"/>
/// implementation: Klein-V₄ does not sit on a single Z₂ axis cleanly. It is the
/// cross-axis primitive that intertwines the bit_a and bit_b axes of the Klein
/// V₄ on the Pauli group itself.</para></summary>
public sealed class Pi2KleinV4DephaseSwapGroup : Claim
{
    /// <summary>The Klein-V₄ closure theorem in one line.</summary>
    public string Theorem =>
        "On the 4^N Pauli basis for any N >= 1, the three unitary involutions " +
        "D, H, Q_zx = H · D form a faithful Klein V₄ subgroup of U(4^N) that " +
        "realizes the dephase-letter Klein V₄ {I, Z↔Y, Z↔X, Y↔X} on the F1 " +
        "palindrome family: D · Π_Z · D = Π_Y, H · Π_X · H = Π_Y, " +
        "Q_zx · Π_Z · Q_zx = Π_X.";

    /// <summary>Build the D operator (diagonal involution, realizes Z↔Y).
    /// D = ⊗_l diag(1, 1, 1, −1) on basis order (I, X, Z, Y). Each diagonal entry
    /// at flat index k equals (−1)^n_Y(k) where n_Y(k) counts Y letters in the
    /// k-th Pauli string (Y is letter index 3 under the <c>a + 2·b</c> packing).
    /// </summary>
    /// <param name="N">number of sites; must be >= 1 and small enough that 4^N
    /// fits in <see cref="int"/>.</param>
    /// <returns>a 4^N × 4^N real diagonal unitary involution as a sparse
    /// <see cref="ComplexMatrix"/>; D² = I.</returns>
    public static ComplexMatrix BuildD(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be >= 1");
        long dim = 1L << (2 * N);
        if (dim > int.MaxValue)
            throw new ArgumentOutOfRangeException(nameof(N), N, $"4^N = {dim} overflows int32");
        var M = Matrix<Complex>.Build.Sparse((int)dim, (int)dim);
        for (long k = 0; k < dim; k++)
        {
            var letters = PauliIndex.FromFlat(k, N);
            int nY = 0;
            foreach (var l in letters)
                if (l == PauliLetter.Y) nY++;
            M[(int)k, (int)k] = (nY & 1) == 0 ? Complex.One : -Complex.One;
        }
        return M;
    }

    /// <summary>Build the H operator (basis-permutation involution, realizes Y↔X).
    /// H = ⊗_l h where h is the 4×4 permutation that fixes the I and Y basis
    /// vectors and swaps the X and Z basis indices (letter index 1 ↔ 2) per site.
    /// H is a pure permutation matrix (all entries 0 or +1), self-inverse.
    ///
    /// <para>On a multi-site Pauli string σ = σ_1 ⊗ ... ⊗ σ_N, H acts by per-site
    /// X↔Z letter swap (I and Y unchanged), no phase factor.</para></summary>
    /// <param name="N">number of sites; must be >= 1 and small enough that 4^N
    /// fits in <see cref="int"/>.</param>
    /// <returns>a 4^N × 4^N real permutation matrix as a sparse
    /// <see cref="ComplexMatrix"/>; H² = I.</returns>
    public static ComplexMatrix BuildH(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be >= 1");
        long dim = 1L << (2 * N);
        if (dim > int.MaxValue)
            throw new ArgumentOutOfRangeException(nameof(N), N, $"4^N = {dim} overflows int32");
        var M = Matrix<Complex>.Build.Sparse((int)dim, (int)dim);
        for (long k = 0; k < dim; k++)
        {
            var letters = PauliIndex.FromFlat(k, N);
            var swapped = new PauliLetter[N];
            for (int i = 0; i < N; i++)
            {
                swapped[i] = letters[i] switch
                {
                    PauliLetter.X => PauliLetter.Z,
                    PauliLetter.Z => PauliLetter.X,
                    _ => letters[i],
                };
            }
            long newK = PauliIndex.ToFlat(swapped);
            M[(int)newK, (int)k] = Complex.One;
        }
        return M;
    }

    /// <summary>Build the Q_zx operator = H · D (realizes Z↔X).
    /// Q_zx = ⊗_l (h · d_l) per site. Order-2 involution because h and d_l commute
    /// per-site, so (h · d_l)² = h² · d_l² = I.</summary>
    /// <param name="N">number of sites; must be >= 1 and small enough that 4^N
    /// fits in <see cref="int"/>.</param>
    /// <returns>a 4^N × 4^N real signed-permutation matrix as a
    /// <see cref="ComplexMatrix"/>; Q_zx² = I.</returns>
    public static ComplexMatrix BuildQzx(int N) => BuildH(N) * BuildD(N);

    public Pi2KleinV4DephaseSwapGroup()
        : base("Pi2 Klein V₄ dephase-swap group: D, H, Q_zx form Klein V₄ subgroup of U(4^N), realizing the dephase-letter Klein V₄ on the F1 palindrome family {Π_Z, Π_X, Π_Y}; Tier1Derived universal N",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md + " +
               "docs/proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md + " +
               "simulations/_d_pi_z_swap_verify.py + " +
               "simulations/_klein_dephase_swap_explore.py + " +
               "reflections/D_PI_Z_EQUALS_PI_Y.md + " +
               "compute/RCPsiSquared.Core/Symmetry/PiOperator.cs")
    {
    }

    public override string DisplayName =>
        "Pi2 Klein V₄ dephase-swap group (D · Π_Z · D = Π_Y; H · Π_Y · H = Π_X; Q_zx · Π_Z · Q_zx = Π_X)";

    public override string Summary =>
        $"Klein V₄ subgroup {{I, D, H, Q_zx}} of U(4^N) realizing dephase-letter swaps " +
        $"{{Z↔Y, Y↔X, Z↔X}} on the F1 palindrome family. D, H, Q_zx all order-2 " +
        $"involutions; D · Q_zx · H = I; pairwise commuting. Per-site factorization " +
        $"gives universal N. ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem", summary: Theorem);
            yield return new InspectableNode("D involution (Z↔Y swap)",
                summary: "D = ⊗ d_l with d_l = diag(1, 1, 1, −1) on (I, X, Z, Y); D² = I; " +
                         "D · Π_Z · D = Π_Y. Pure diagonal involution. See " +
                         "PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md (Welle 12 Task 1).");
            yield return new InspectableNode("H involution (Y↔X swap)",
                summary: "H = ⊗ h with h the 4×4 X↔Z basis-index permutation (I, Y fixed); " +
                         "H² = I; H · Π_Y · H = Π_X (and H · Π_X · H = Π_Y). Pure permutation " +
                         "matrix; all signs +1. NOTE: h permutes basis indices, NOT the " +
                         "physical Pauli operators σ_X ↔ σ_Z.");
            yield return new InspectableNode("Q_zx involution (Z↔X swap)",
                summary: "Q_zx = H · D = ⊗ (h · d_l). Q_zx² = I (because h and d_l commute " +
                         "per-site); Q_zx · Π_Z · Q_zx = Π_X.");
            yield return new InspectableNode("Klein-V₄ closure",
                summary: "D · Q_zx · H = I; D, Q_zx, H commute pairwise. The three involutions " +
                         "plus identity form Klein V₄ ≃ Z₂ × Z₂ in U(4^N) for any N. Any two " +
                         "determine the third.");
            yield return new InspectableNode("Universal-N reduction",
                summary: "Each identity reduces to a per-site 4×4 check via the Kronecker " +
                         "mixed-product property: (⊗ a_l)(⊗ b_l)(⊗ c_l) = ⊗(a_l · b_l · c_l). " +
                         "Verified bit-exact at N = 1, 2, 3, 4 via numpy + sympy symbolic.");
            yield return new InspectableNode("Significance",
                summary: "F1-family identities transfer between dephase letters by two routes: " +
                         "Route 1 (per-axis re-derivation): identities depending only on F38 Π_d² " +
                         "eigenvalue structure + Pauli-support disjointness transfer universally " +
                         "(all three d), substituting axis_d (bit_b for Z/Y, bit_a for X). The " +
                         "F112 Welle 13 closure (PROOF_F112_CROSS_DEPHASE_VIA_KLEIN_V4.md) is the " +
                         "first formalization. Route 2 (Hadamard transport via Q_zx): identities " +
                         "depending on L itself transfer only between (Z, X) via Q_zx, which is " +
                         "the operator-space lift of U_H^⊗N. D and Q_yx (= H) are operator-space-" +
                         "only and do NOT transport L between dephase representations: D is the " +
                         "transpose ρ ↦ ρ^T, so D·L_Z·D⁻¹ stays a valid Z-dephasing Lindbladian " +
                         "(= L_Z with H → −H^T) — it keeps Lindblad form but does not rotate the " +
                         "dephasing axis (no Hilbert-space lift for D). F108 Welle 14 closure " +
                         "(PROOF_F108_KLEIN_V4_EQUIVALENCE.md): F108 Part 3 follows from Part 1 via " +
                         "OPERATOR-SPACE D-conjugation (D · Π_5b(Z) · D = Π_5b(Y) bit-exact; bilinear " +
                         "set fixed); F108 Part 2 follows from Part 1 via HILBERT-SPACE Hadamard " +
                         "transport (U_op = U_H^⊗N ⊗ (U_H^⊗N)^*; bilinear-set bijection). " +
                         "PERMUTED operator-space Klein-V₄ for Π_5b: the canonical-Π pairings fail " +
                         "(Q_zx · Π_5b(Z) · Q_zx ≠ ±Π_5b(X), H · Π_5b(Y) · H ≠ ±Π_5b(X), residual 2.0 " +
                         "at N=1,2,3), but H serves the Z↔X pairing instead: H · Π_5b(Z) · H = Π_5b(X) " +
                         "bit-exact (closed 2026-05-29). So {I, D, H} acts on Π_5b at the operator level " +
                         "(D: Z↔Y, H: Z↔X); only Q_zx (= H·D) leaves the canonical set. F112-Y cannot be " +
                         "obtained 'for free' from F112-Z via D-conjugation; the correct path is Route 1.");
            yield return new InspectableNode("No IZ2AxisClaim implementation",
                summary: "Klein-V₄ does not sit on a single Z₂ axis cleanly. It is the cross-axis " +
                         "primitive that intertwines the bit_a and bit_b axes of the Klein V₄ on " +
                         "the Pauli group itself: X = (1, 0) and Z = (0, 1) are the single-axis " +
                         "generators, and swapping their basis indices via h is precisely the V₄ " +
                         "involution (a, b) ↔ (b, a) that maps the bit_a-axis to the bit_b-axis.");
            yield return new InspectableNode("Action precision: four oriented operators, not three letters (2026-06-20)",
                summary: "The conjugation action is the REGULAR V₄-action on the four oriented " +
                         "palindrome operators {Π_X, Π_X⁻¹, Π_Y = Π_Z⁻¹, Π_Z}, not a permutation " +
                         "of the three letters {X, Y, Z} (three transpositions would generate S₃, " +
                         "and the bare-letter map is ill-defined for Q_zx and H). Π_Y = Π_Z⁻¹ " +
                         "exactly at all N (per-site π_Y·π_Z = I), so the three letters carry only " +
                         "two independent operators and the identity-letter slot is realized by " +
                         "Π_X⁻¹ (I ↦ Π_X⁻¹, X ↦ Π_X, Y ↦ Π_Y, Z ↦ Π_Z). Each generator is a " +
                         "fixed-point-free double transposition (e.g. D = (Π_X Π_X⁻¹)(Π_Y Π_Z)); " +
                         "ambient group B₂ ≅ D₄. The four oriented operators are a copy of the " +
                         "Pauli Klein group {I, X, Y, Z} itself — which is why 'Klein-V₄' is the " +
                         "right name. (The three dephasing-DIAGONAL operators Q_v carry the " +
                         "complementary S₃; the V₄-vs-S₃ difference is the Π-vs-Q distinction.)");
            yield return new InspectableNode("Downstream sibling: F114 commutator D-conjugation parity",
                summary: "CommutatorDConjugationSign (F114, Tier1Derived, ctor child of this Claim, " +
                         "Welle 15 closure 2026-05-27) characterizes D's action on the H-commutator " +
                         "superoperator L_σ = −i[σ, ·] with closed form ε(σ) = (−1)^{n_Y(σ) + 1}. " +
                         "Whereas this Claim makes D the Π swap-operator across {Z, Y} dephase letters " +
                         "(D · Π_Z · D = Π_Y), F114 makes D the L_H sign-flip-operator with per-term " +
                         "n_Y bookkeeping. Together the two Claims characterize D's action on the two " +
                         "main dephase-letter-sensitive structures (Π and L_H).");
        }
    }
}
