using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F5 closed form for depolarizing noise: <c>error = γ · 2·(N−2)/3</c>
/// (Tier 1 proven). Linear in γ and N, Hamiltonian-independent. Each constant
/// in this closed form sits on the Pi2-Foundation:
///
/// <list type="bullet">
///   <item><b>"2" multiplier</b>: <see cref="Pi2DyadicLadderClaim.Term"/>(0) =
///         <c>a_0</c> = d (qubit dimension; number-anchor side of d=2).</item>
///   <item><b>"3" denominator</b>: <c>d² − 1 = a_{−1} − 1 = 4 − 1 = 3</c>. The
///         number of non-identity Pauli operators per qubit (X, Y, Z) — i.e. all
///         operators in the per-site Pauli basis except the identity. Direct
///         consequence of <see cref="Pi2OperatorSpaceMirrorClaim"/> at N=1: total
///         operator-space d² = 4 minus the identity = 3.</item>
///   <item><b>"(N−2)" linear factor</b>: same chain-derivation overhead shift as
///         <see cref="F49Pi2Inheritance"/> uses ("(N−2) qubits" reading); the
///         "−2 shift" is structural in the F1-residual scaling derivation, NOT a
///         free parameter.</item>
/// </list>
///
/// <para>The composite ratio <c>2/3 = d/(d²−1) = a_0/(a_{−1}−1)</c> is the
/// fraction of off-diagonal Pauli mass in the depolarizing channel: 2 of the 3
/// non-identity Paulis (X, Y) are anti-commute-with-Z (off-diagonal in the
/// dephasing eigenbasis at d=2), 1 (Z) commutes (diagonal). The d=2 split into
/// 2:2 (immune:decaying, <see cref="QubitNecessityPi2Inheritance"/>) becomes
/// I-removed split into 1:2 (immune-Z : decaying-X-Y) for the depolarizing
/// channel at unit budget. Hence 2/3.</para>
///
/// <para>Tier outcome: <b>Tier1Derived</b>. F5 is Tier 1 proven in
/// <c>experiments/DEPOLARIZING_PALINDROME.md</c> + ANALYTICAL_FORMULAS F5; the
/// Pi2-Foundation anchoring of its three constants is algebraic-trivial
/// composition once F5 is granted.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F5 +
/// <c>experiments/DEPOLARIZING_PALINDROME.md</c> +
/// <c>compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs</c> (the parent F1
/// identity which BREAKS for depolarizing, with this measured residual) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs</c>.</para></summary>
public sealed class F5DepolarizingErrorPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly Pi2OperatorSpaceMirrorClaim _mirror;

    /// <summary>The "2" multiplier in F5's <c>2·(N−2)/3</c>. Exactly equal to
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(0) = <c>a_0</c> = d.</summary>
    public double DCoefficient => _ladder.Term(0);

    /// <summary>The "3" denominator: <c>d² − 1 = 4 − 1 = 3</c>. The number of
    /// non-identity Pauli operators per qubit. Pi2-derived from
    /// <see cref="Pi2OperatorSpaceMirrorClaim"/>: total operator-space at N=1 is
    /// 4 (= d²); subtracting the identity gives 3 non-identity Paulis (X, Y, Z).</summary>
    public double DSquaredMinusOne =>
        (_mirror.PairAt(1)?.OperatorSpace ?? 4.0) - 1.0;

    /// <summary>The composite ratio <c>2/3 = d/(d²−1) = DCoefficient / DSquaredMinusOne</c>.
    /// The fraction of off-diagonal Pauli mass in the depolarizing channel.</summary>
    public double TwoOverThree => DCoefficient / DSquaredMinusOne;

    /// <summary>The chain-derivation-overhead linear factor <c>(N − 2)</c> in F5's
    /// closed form. Same shift as <see cref="F49Pi2Inheritance"/>'s "(N−2) qubits"
    /// reading: structural F1-residual derivation, NOT a free parameter.</summary>
    public int NShiftFactor(int N) => N - 2;

    /// <summary>Live re-composition of F5's closed-form coefficient at γ=1:
    /// <c>2·(N−2)/3 = DCoefficient · NShiftFactor(N) / DSquaredMinusOne</c>.
    /// Drift between the algebraic Pi2-derived value and the F5 closed form
    /// surfaces here.</summary>
    public double LiveCoefficient(int N)
    {
        if (N < 2)
            throw new ArgumentOutOfRangeException(nameof(N), N, "F5 closed form requires N ≥ 2.");
        return DCoefficient * NShiftFactor(N) / DSquaredMinusOne;
    }

    public F5DepolarizingErrorPi2Inheritance(Pi2DyadicLadderClaim ladder, Pi2OperatorSpaceMirrorClaim mirror)
        : base("F5 depolarizing error coefficients (2 · (N−2) / 3) inherit from Pi2-Foundation",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F5 + " +
               "experiments/DEPOLARIZING_PALINDROME.md + " +
               "compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _mirror = mirror ?? throw new ArgumentNullException(nameof(mirror));
    }

    public override string DisplayName =>
        "F5 depolarizing error coefficients as Pi2-Foundation inheritance";

    public override string Summary =>
        $"error = γ·2·(N−2)/3: 2 = a_0 = d; 3 = a_{{-1}} − 1 = d² − 1; (N−2) = derivation-overhead shift " +
        $"(same as F49); 2/3 = d/(d²−1) = {TwoOverThree:F4} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F5 closed form",
                summary: "error = γ · 2·(N−2)/3 (Tier 1 proven; linear in γ and N; Hamiltonian-independent)");
            yield return InspectableNode.RealScalar("DCoefficient (= a_0 = d)", DCoefficient);
            yield return InspectableNode.RealScalar("DSquaredMinusOne (= a_{-1} − 1 = d² − 1 = 3)", DSquaredMinusOne);
            yield return InspectableNode.RealScalar("TwoOverThree (= d / (d² − 1))", TwoOverThree);
            yield return new InspectableNode("interpretation",
                summary: "2 = number of off-diagonal Paulis (X, Y); 3 = total non-identity Paulis (X, Y, Z); 2/3 is the off-diagonal Pauli mass fraction in the depolarizing channel");
            yield return new InspectableNode("(N−2) shift",
                summary: "same chain-derivation-overhead as F49Pi2Inheritance; structural, not a free parameter");
            for (int N = 2; N <= 6; N++)
            {
                yield return new InspectableNode(
                    $"chain N={N}",
                    summary: $"(N−2) = {NShiftFactor(N)}; coefficient 2·(N−2)/3 = {LiveCoefficient(N):F4}");
            }
        }
    }
}
