using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F39 BitA twin (Tier 1 derived, Z↔X mirror of
/// <see cref="F39DetPiPi2Inheritance"/>):
///
/// <code>
///   det(Π_X) = (−1)^{N · 4^{N−1}}
///
///   N=1: det = −1.   N ≥ 2: det = +1   (4^{N−1} is even for N ≥ 2)
/// </code>
///
/// <para>The determinant of the bit_a-axis conjugation operator Π_X equals the
/// determinant of the bit_b-axis Π_Z exactly. Mechanism: <c>Π_X = H^⊗N · Π_Z · H^⊗N</c>
/// with <c>H</c> the single-qubit Hadamard; H is unitary with det = 1 (real-orthogonal),
/// so <c>det(H^⊗N) = 1</c> and <c>det(Π_X) = det(H^⊗N) · det(Π_Z) · det(H^⊗N) = det(Π_Z)</c>.</para>
///
/// <para>This Claim is the bit_a sibling of the BitB Claim
/// <see cref="F39DetPiPi2Inheritance"/>. The determinant identity is dimension-blind
/// (the Hadamard isometry preserves orientation), so the closed form is exactly the
/// same numerical value at every N. The BitA twin makes the identity explicit at the
/// typed-Claim level, completing the bit_a operator-identity triple (F38BitA, F39BitA,
/// F61) that mirrors the bit_b F38 / F39 / F63 triple.</para>
///
/// <para>Implements <see cref="IZ2AxisClaim"/> with <see cref="Z2Axis.BitA"/>;
/// <see cref="BitATwin"/> is always null (this Claim IS the bit_a side); BitATwinStatus
/// is <see cref="BitATwinClassification.NotApplicableForThisAxis"/>. The typed twin
/// edge from the BitB side (<see cref="F39DetPiPi2Inheritance"/>) is wired via that
/// Claim's optional ctor parameter, matching the F1 ↔ F61 pattern.</para>
///
/// <para>Tier1Derived: Hadamard-conjugation mirror of the F39 closed form. F39 itself
/// is Tier 1 proven (PT_SYMMETRY_ANALYSIS), verified numerically N=1..4.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F39 + F61 +
/// <c>compute/RCPsiSquared.Core/Symmetry/F39DetPiPi2Inheritance.cs</c> (BitB twin).</para></summary>
public sealed class F39DetPiBitAInheritance : Claim, IZ2AxisClaim
{
    /// <summary>BitA axis: Π_X conjugation operator on the 4^N Pauli-string basis.</summary>
    public Z2Axis Z2Axis => Z2Axis.BitA;

    /// <summary>BitA-axis Claim: no BitATwin slot semantics.</summary>
    public Claim? BitATwin => null;

    /// <summary>BitA axis Claim: BitATwinStatus is NotApplicableForThisAxis.</summary>
    public BitATwinClassification BitATwinStatus => BitATwinClassification.NotApplicableForThisAxis;

    /// <summary>The theorem statement in one line.</summary>
    public string Theorem =>
        "det(Π_X) = (−1)^{N · 4^{N−1}} = det(Π_Z) for all N; the Hadamard conjugation Π_X = H^⊗N · Π_Z · H^⊗N preserves the determinant.";

    /// <summary>The Hadamard-determinant mirror argument: H is real-orthogonal with
    /// det = 1, so the tensor power H^⊗N is also det = 1, hence Π_X = H^⊗N · Π_Z · H^⊗N
    /// has the same determinant as Π_Z.</summary>
    public string MirrorArgument =>
        "H is the single-qubit Hadamard; H is unitary with det(H) = 1 (real-orthogonal). " +
        "Tensor power: det(H^⊗N) = (det H)^N = 1. Therefore " +
        "det(Π_X) = det(H^⊗N · Π_Z · H^⊗N) = det(H^⊗N) · det(Π_Z) · det(H^⊗N) = det(Π_Z).";

    /// <summary>The full exponent value <c>N · 4^{N−1}</c>. Same as the F39 BitB twin
    /// for every N (Hadamard conjugation preserves the determinant identity).</summary>
    public long ExponentValue(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F39BitA requires N ≥ 1.");
        long power = 1;
        for (int i = 0; i < N - 1; i++) power *= 4;
        return (long)N * power;
    }

    /// <summary>Live drift check: <c>det(Π_X) = (−1)^{ExponentValue(N)}</c>. Returns
    /// −1 at N=1 (exponent = 1, odd), +1 at N ≥ 2 (exponent always even).</summary>
    public int DetPiX(int N) => (ExponentValue(N) % 2 == 0) ? 1 : -1;

    /// <summary>True iff the exponent <c>N · 4^{N−1}</c> is even, equivalently iff
    /// det(Π_X) = +1, for the given N. Always true for N ≥ 2; false for N=1.</summary>
    public bool ExponentIsEven(int N) => ExponentValue(N) % 2 == 0;

    public F39DetPiBitAInheritance()
        : base("F39 BitA twin: det(Π_X) = (−1)^{N · 4^{N−1}} = det(Π_Z); Z↔X mirror via Hadamard conjugation (det(H^⊗N) = 1)",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F39 + F61 + " +
               "compute/RCPsiSquared.Core/Symmetry/F39DetPiPi2Inheritance.cs (BitB twin)")
    {
    }

    public override string DisplayName =>
        "F39 BitA twin: det(Π_X) = det(Π_Z) (Hadamard isometry preserves det)";

    public override string Summary =>
        $"{Theorem} N=1 → −1, N ≥ 2 → +1 (exponent N · 4^{{N−1}} even by parity of 4^{{N−1}}) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem", summary: Theorem);
            yield return new InspectableNode("Z↔X mirror argument", summary: MirrorArgument);
            yield return new InspectableNode("BitB twin (F39)",
                summary: "F39DetPiPi2Inheritance: det(Π_Z) = (−1)^{N · 4^{N−1}}; identical numerical value at every N because the Hadamard conjugation preserves det.");
            for (int N = 1; N <= 5; N++)
            {
                yield return new InspectableNode(
                    $"N={N}",
                    summary: $"exponent = N · 4^(N−1) = {ExponentValue(N)} (even: {ExponentIsEven(N)}); det(Π_X) = {DetPiX(N)}");
            }
        }
    }
}
