using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F38 BitA twin (Tier 1 derived, Z↔X mirror of
/// <see cref="F38Pi2InvolutionPi2Inheritance"/>):
///
/// <code>
///   Π²_X = (−1)^{n_XY}
///
///   on the 4^N Pauli-string basis, with eigenvalues split equally:
///   |Π²_X-even| = |Π²_X-odd| = 4^N / 2 = 2 · 4^(N−1)
/// </code>
///
/// <para>This Claim is the bit_a sibling of the BitB Claim
/// <see cref="F38Pi2InvolutionPi2Inheritance"/>: the same algebraic statement
/// translated by the Z↔X letter swap. Π²_X = H^⊗N · Π²_Z · H^⊗N where H is the
/// single-qubit Hadamard; the Hadamard isometry preserves the half-half eigenspace
/// count exactly (the dimension count is letter-blind).</para>
///
/// <para>F61BitAParity ([L, Π²_X] = 0) is the conservation companion in the bit_a
/// branch, matching how F63 ([L, Π²_Z] = 0) sits next to F38 in the bit_b branch.
/// Together with F39BitA (det(Π_X) = det(Π_Z)) and F61, F38BitA closes the bit_a
/// operator-identity triple that mirrors the bit_b F38 / F39 / F63 triple.</para>
///
/// <para>Implements <see cref="IZ2AxisClaim"/> with <see cref="Z2Axis.BitA"/>;
/// <see cref="BitATwin"/> is always null (this Claim IS the bit_a side); BitATwinStatus
/// is <see cref="BitATwinClassification.NotApplicableForThisAxis"/>. The typed twin
/// edge from the BitB side (<see cref="F38Pi2InvolutionPi2Inheritance"/>) is wired
/// via that Claim's optional ctor parameter, matching the F1 ↔ F61 pattern.</para>
///
/// <para>Tier1Derived: Z↔X mirror of the F38 closed form. The Hadamard conjugation
/// preserves involutivity (Π_X² = (H^⊗N · Π_Z · H^⊗N)² = H^⊗N · Π_Z² · H^⊗N) and
/// preserves the eigenvalue multiplicities exactly. F38 itself is Tier 1 proven
/// (PT_SYMMETRY_ANALYSIS, PROOF_BIT_B_PARITY_SYMMETRY).</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F38 + F61 +
/// <c>compute/RCPsiSquared.Core/Symmetry/F38Pi2InvolutionPi2Inheritance.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F61BitAParityPi2Inheritance.cs</c>
/// (the [L, Π²_X] = 0 companion).</para></summary>
public sealed class F38BitAInvolutionInheritance : Claim, IZ2AxisClaim
{
    /// <summary>BitA axis: Π²_X = Z⊗N, bit_a parity n_XY = #X + #Y mod 2.
    /// Sister axis to F38's BitB (Π²_Z = X⊗N, bit_b parity #Y + #Z mod 2).</summary>
    public Z2Axis Z2Axis => Z2Axis.BitA;

    /// <summary>BitA-axis Claim: no BitATwin slot semantics (the twin concept lives
    /// on BitB Claims pointing at BitA siblings; F38BitA IS the bit_a side).</summary>
    public Claim? BitATwin => null;

    /// <summary>BitA axis Claim: BitATwinStatus is NotApplicableForThisAxis. The
    /// typed twin edge in the BitB direction is wired via
    /// <see cref="F38Pi2InvolutionPi2Inheritance"/>'s optional ctor parameter.</summary>
    public BitATwinClassification BitATwinStatus => BitATwinClassification.NotApplicableForThisAxis;

    /// <summary>The theorem statement in one line.</summary>
    public string Theorem =>
        "Π²_X = (−1)^{n_XY} on the 4^N Pauli-string basis; Π²_X-even and Π²_X-odd subspaces both have dimension 4^N / 2.";

    /// <summary>The Z↔X mirror argument used to derive the BitA twin from the BitB
    /// F38 closed form. The Hadamard isometry is unitary and preserves dimensions
    /// of eigenspaces exactly.</summary>
    public string MirrorArgument =>
        "Π_X = H^⊗N · Π_Z · H^⊗N where H is the single-qubit Hadamard; squaring gives Π²_X = H^⊗N · Π²_Z · H^⊗N. " +
        "The Hadamard conjugation preserves involutivity (Π²_X² = I) and eigenvalue multiplicities (4^N / 2 each), " +
        "and exchanges the bit_b parity #Y + #Z mod 2 with the bit_a parity #X + #Y mod 2 letter-by-letter.";

    /// <summary>The two distinct Π²_X eigenvalues, in canonical order: +1
    /// (Π²_X-even) then −1 (Π²_X-odd). Each is realised on exactly half of the
    /// 4^N operator space.</summary>
    public IReadOnlyList<int> Pi2XEigenvalues => new[] { +1, -1 };

    /// <summary>The full operator-space dimension <c>4^N</c> for an N-qubit
    /// system. Identical to F38's BitB count by the Z↔X mirror.</summary>
    public long FullOperatorSpaceDimension(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F38BitA requires N ≥ 1.");
        return 1L << (2 * N);
    }

    /// <summary>The dimension of each Π²_X eigenspace: <c>4^N / 2 = 2 · 4^(N−1)</c>.
    /// Both Π²_X-even and Π²_X-odd sectors have this dimension exactly.</summary>
    public long EigenspaceDimension(int N) => FullOperatorSpaceDimension(N) / 2L;

    /// <summary>The half-half balance ratio: <c>EigenspaceDim / FullDim = 1/2</c>
    /// exactly. Matches the F38 BitB anchor via the Hadamard isometry.</summary>
    public double HalfHalfBalance(int N) =>
        (double)EigenspaceDimension(N) / FullOperatorSpaceDimension(N);

    public F38BitAInvolutionInheritance()
        : base("F38 BitA twin: Π²_X = (−1)^{n_XY} on the 4^N Pauli-string basis; Z↔X mirror of F38 via Hadamard conjugation; half-half eigenspace split preserved",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F38 + F61 + " +
               "compute/RCPsiSquared.Core/Symmetry/F38Pi2InvolutionPi2Inheritance.cs (BitB twin) + " +
               "compute/RCPsiSquared.Core/Symmetry/F61BitAParityPi2Inheritance.cs (the [L, Π²_X] = 0 companion)")
    {
    }

    public override string DisplayName =>
        "F38 BitA twin: Π²_X = (−1)^{n_XY} (Z↔X mirror of F38)";

    public override string Summary =>
        $"{Theorem} Z↔X mirror of F38 via Π_X = H^⊗N · Π_Z · H^⊗N (Hadamard isometry preserves involutivity and 4^N / 2 eigenspace count) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem", summary: Theorem);
            yield return new InspectableNode("Z↔X mirror argument", summary: MirrorArgument);
            yield return new InspectableNode("BitB twin (F38)",
                summary: "F38Pi2InvolutionPi2Inheritance: Π²_Z = (−1)^{w_YZ} on the 4^N Pauli-string basis; same half-half split via bit_b parity instead of bit_a.");
            yield return new InspectableNode("Conservation companion (F61)",
                summary: "F61BitAParityPi2Inheritance: [L, Π²_X] = 0 exactly all N; F38BitA defines Π²_X, F61 says L respects it. Parallel to F38 (defines Π²_Z) ↔ F63 ([L, Π²_Z] = 0) in the bit_b branch.");
            yield return new InspectableNode("Pi2 eigenvalues",
                summary: "{+1, −1}, each on exactly 4^N / 2 = 2 · 4^(N−1) Pauli strings (Hadamard preserves multiplicities).");
            for (int N = 1; N <= 5; N++)
            {
                yield return new InspectableNode(
                    $"N={N}",
                    summary: $"4^N = {FullOperatorSpaceDimension(N)}; each sector = {EigenspaceDimension(N)}; balance = {HalfHalfBalance(N):G6}");
            }
        }
    }
}
