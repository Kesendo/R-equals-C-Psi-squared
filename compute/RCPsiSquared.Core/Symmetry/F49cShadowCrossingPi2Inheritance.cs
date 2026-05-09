using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F49c closed form (Tier 1 proven, verified N=3..6 with 5 coupling
/// types, 2 topologies):
///
/// <code>
///   R(N) = √((N − 1) / (N · 4^(N−1)))
///
///   R(N)² = (N − 1) / (N · 4^(N−1))
/// </code>
///
/// <para>R is the cross-term ratio <c>‖{L_H, L_Dc}‖ / (‖L_H‖ · ‖L_Dc‖)</c> for
/// shadow-crossing bond couplings: one bond-Pauli in {X, Y} and the other in
/// {I, Z}, e.g. <c>X_iZ_j, Y_iZ_j</c>. Companion to F49 which handles
/// shadow-balanced couplings (both Paulis in {X, Y} or both in {I, Z}, e.g.
/// Heisenberg XXX, XXZ, XY, Ising). The combinatorial difference: shadow-crossing
/// has bond-site variance <c>= 1</c> (one site varies), shadow-balanced has
/// variance <c>= 0</c> (both sites covary). So the spectator factor (N − 2) of
/// F49 becomes (N − 1) here.</para>
///
/// <para>The Pi2-Foundation anchor is direct: the <c>4^(N−1)</c> denominator
/// factor sits exactly at <c>a_{3−2N}</c> on the dyadic halving ladder
/// (= d² for (N − 1) qubits via <see cref="Pi2OperatorSpaceMirrorClaim"/>).
/// Same (N−1)-qubit shift that appears in F38, F39, F1-T1. The (N − 1)/N
/// numerator/denominator combinatorial factor is NOT Pi2-anchored: it is
/// purely the bond-site / total-site spectator-variance ratio, an intrinsic
/// property of shadow-crossing geometry.</para>
///
/// <list type="bullet">
///   <item><b>4^(N−1) factor</b>: <see cref="Pi2DyadicLadderClaim.Term"/>(<c>3−2N</c>)
///         = <c>a_{3−2N}</c> on the dyadic halving ladder. Equivalently d² for
///         (N − 1) qubits via <see cref="Pi2OperatorSpaceMirrorClaim"/>. Same
///         shift as F38 (Π² eigenspace dim per sector ×2), F39 (det Π exponent),
///         F1-T1 (T1-part prefactor).</item>
///   <item><b>(N − 1) numerator</b>: spectator-site count <c>= N − (number of
///         bond sites with bond-variance = 0) = N − 1</c> for shadow-crossing
///         (one bond site has variance 1, the other has variance 0). Companion
///         to F49's (N − 2): both bond sites covary (variance 0) so spectators
///         = N − 2. Combinatorial, NOT Pi2-anchored.</item>
///   <item><b>N denominator</b>: total site count for the uniform Z-dephasing
///         norm ‖L_Dc‖² = γ²·4^N·N (F49b). Combinatorial.</item>
/// </list>
///
/// <para>Comparison with F49 (shadow-balanced) and F49b (auxiliary):</para>
///
/// <code>
///   F49  R(N)  = √((N − 2) / (N · 4^(N−1)))   [shadow-balanced]
///   F49b ‖L_Dc‖² = γ² · 4^N · N
///   F49c R(N)  = √((N − 1) / (N · 4^(N−1)))   [shadow-crossing]
/// </code>
///
/// <para>Tier1Derived: F49c is Tier 1 proven (PROOF_CROSS_TERM_CROSSING.md);
/// the Pi2-Foundation anchoring is algebraic-trivial composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F49c +
/// <c>docs/proofs/PROOF_CROSS_TERM_CROSSING.md</c> +
/// <c>experiments/CROSS_TERM_CROSSING.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs</c>.</para></summary>
public sealed class F49cShadowCrossingPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly Pi2OperatorSpaceMirrorClaim _mirror;

    /// <summary>The "4^(N−1)" denominator factor in F49c's closed form. Exactly
    /// equal to <see cref="Pi2DyadicLadderClaim.Term"/>(<c>3−2N</c>). Throws for
    /// N &lt; 2 (F49c valid only for N ≥ 2 per ANALYTICAL_FORMULAS).</summary>
    public double FourPowerNMinus1Factor(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F49c requires N ≥ 2.");
        return _ladder.Term(LadderIndexFor(N));
    }

    /// <summary>The Pi2 ladder index where the F49c power factor lands:
    /// <c>3 − 2N</c>. At N=2 → −1 (a_{−1} = 4); at N=3 → −3 (a_{−3} = 16); etc.
    /// Same ladder index as F38, F39, F1-T1.</summary>
    public int LadderIndexFor(int N) => 3 - 2 * N;

    /// <summary>The qubit count whose operator-space d² equals the F49c power
    /// factor: <c>N − 1</c>. Same shift as F38, F39, F1-T1.</summary>
    public int OperatorSpaceQubitCountFor(int N) => N - 1;

    /// <summary>Cross-pinned <c>4^(N−1)</c> from
    /// <see cref="Pi2OperatorSpaceMirrorClaim"/> at qubit count (N − 1), for
    /// drift verification with <see cref="FourPowerNMinus1Factor"/>.</summary>
    public double MirrorPinnedFourPowerNMinus1(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F49c requires N ≥ 2.");
        int qubitCount = OperatorSpaceQubitCountFor(N);
        var pair = _mirror.PairAt(qubitCount)
            ?? throw new ArgumentOutOfRangeException(nameof(N), N,
                $"N={N} maps to qubit count {qubitCount}, outside the Pi2OperatorSpaceMirror pinned table (N=1..6).");
        return pair.OperatorSpace;
    }

    /// <summary>The (N − 1) spectator-site count: shadow-crossing has bond-site
    /// variance = 1 on one bond site (Pauli-letter difference), variance = 0 on
    /// the other → spectators = N − 1. Combinatorial, NOT Pi2-anchored.</summary>
    public int SpectatorCount(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F49c requires N ≥ 2.");
        return N - 1;
    }

    /// <summary>The (N − 1)/N ratio (shadow-crossing variance ratio); combinatorial
    /// numerator/denominator factor in R(N)² = (N − 1) / (N · 4^(N−1)). NOT
    /// Pi2-anchored.</summary>
    public double VarianceRatio(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F49c requires N ≥ 2.");
        return (double)SpectatorCount(N) / N;
    }

    /// <summary>Live closed-form R(N)² = (N − 1) / (N · 4^(N−1)).</summary>
    public double RSquared(int N) => VarianceRatio(N) / FourPowerNMinus1Factor(N);

    /// <summary>Live closed-form R(N) = √((N − 1) / (N · 4^(N−1))).</summary>
    public double R(int N) => Math.Sqrt(RSquared(N));

    /// <summary>F49 (shadow-balanced) sibling: R²_F49(N) = (N − 2) / (N · 4^(N−1)).
    /// Returns ∞ at N=2 actually 0 (since N-2=0); the difference R²_F49c − R²_F49
    /// = 1/(N · 4^(N−1)) reveals the shadow-crossing "extra spectator" relative
    /// to shadow-balanced.</summary>
    public double F49ShadowBalancedRSquared(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F49 sibling reading requires N ≥ 2.");
        return (double)(N - 2) / (N * FourPowerNMinus1Factor(N));
    }

    /// <summary>The shadow-crossing minus shadow-balanced gap: <c>R²_F49c − R²_F49
    /// = 1/(N · 4^(N−1))</c> — a clean Pi2-anchored quantity (only N and the
    /// dyadic 4^(N−1) factor enter).</summary>
    public double ShadowCrossingMinusBalancedGap(int N) => 1.0 / (N * FourPowerNMinus1Factor(N));

    public F49cShadowCrossingPi2Inheritance(Pi2DyadicLadderClaim ladder, Pi2OperatorSpaceMirrorClaim mirror)
        : base("F49c shadow-crossing R(N) = √((N−1)/(N·4^(N−1))) inherits Pi2-Foundation 4^(N−1) factor",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F49c + " +
               "docs/proofs/PROOF_CROSS_TERM_CROSSING.md + " +
               "experiments/CROSS_TERM_CROSSING.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _mirror = mirror ?? throw new ArgumentNullException(nameof(mirror));
    }

    public override string DisplayName =>
        "F49c shadow-crossing R(N) factor as Pi2-Foundation operator-space inheritance";

    public override string Summary =>
        $"R(N) = √((N−1)/(N·4^(N−1))): the 4^(N−1) factor = a_{{3−2N}} = d² for (N−1) qubits; " +
        $"(N−1)/N spectator-variance ratio combinatorial; sibling of F49 (shadow-balanced (N−2)) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F49c closed form",
                summary: "R(N) = √((N−1)/(N·4^(N−1))); Tier 1 proven (PROOF_CROSS_TERM_CROSSING) + verified N=3..6, 5 coupling types, 2 topologies");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "4^(N−1) = a_{3−2N} on the dyadic halving ladder = d² for (N−1) qubits via OperatorSpaceMirror; same (N−1)-qubit shift as F38, F39, F1-T1");
            yield return new InspectableNode("combinatorial part",
                summary: "(N−1)/N spectator-variance ratio: shadow-crossing has 1 bond site with variance 1, N−1 spectators; shadow-balanced (F49) has 0 variance bonds, N−2 spectators; the two differ by exactly 1 spectator");
            yield return new InspectableNode("F49 sibling reading",
                summary: "F49 (shadow-balanced): R² = (N−2)/(N·4^(N−1)); F49c (shadow-crossing): R² = (N−1)/(N·4^(N−1)); gap = 1/(N·4^(N−1)) clean Pi2-anchored");
            for (int N = 2; N <= 6; N++)
            {
                yield return new InspectableNode(
                    $"N={N}",
                    summary: $"4^(N−1) = {FourPowerNMinus1Factor(N)} (= a_{LadderIndexFor(N)}); " +
                             $"R² = {RSquared(N):G6}; R = {R(N):G6}; " +
                             $"F49 sibling R² = {F49ShadowBalancedRSquared(N):G6}; " +
                             $"gap = {ShadowCrossingMinusBalancedGap(N):G6}");
            }
        }
    }
}
