using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F69 closed form (Tier 1, sextic minimal polynomial, verified N=3):
///
/// <code>
///   |ψ(α)⟩ = α |GHZ_3⟩ + √(1-α²) |W_3⟩,  α ∈ [0, 1]
///
///   C(α)      = Tr(ρ_AB²)   = -5α⁴/18 + 2α²/9 + 5/9
///   L1_off(α) = √6 α √(1-α²) + (2/3)(1-α²)
///   CΨ(α)     = C(α) · L1_off(α) / 3
///
///   Stationarity dCΨ/dα = 0, in x = α², gives the integer-coefficient sextic:
///   P(x) = 2900 x⁶ - 8060 x⁵ + 4211 x⁴ + 3832 x³ - 2428 x² - 512 x + 300 = 0
///
///   Irreducible over ℚ (sympy Poly.is_irreducible = True).
/// </code>
///
/// <para>F69 is the structural finding that the GHZ+W superposition at N=3
/// can lift pair-CΨ(0) ABOVE the 1/4 fold, even though both pure components
/// (F60 GHZ and F62 W) sit strictly below it. The optimum α²_opt is a degree-6
/// algebraic number with NO radical form in nested square roots; the sextic
/// is irreducible over ℚ, so no closed-form expression in radicals exists.</para>
///
/// <para><b>Optimum (computed from sextic root at 25-digit precision):</b></para>
/// <list type="bullet">
///   <item>α²_opt ≈ 0.375420720711069</item>
///   <item>α_opt ≈ 0.612715856422101</item>
///   <item>min pair-CΨ(0) ≈ 0.320411541127025</item>
///   <item>ratio to 1/4 fold ≈ 1.281646× (28% above the bilinear-apex maxval)</item>
///   <item>3-tangle τ_ABC ≈ 0.799453 (near-GHZ limit)</item>
///   <item>pair concurrence C(A,B) ≈ 0.0210 (essentially zero)</item>
/// </list>
///
/// <para>Pi2-Foundation anchors:</para>
/// <list type="bullet">
///   <item><b>FoldPosition = 1/4 = a_3</b>: same dyadic-ladder anchor that
///         <see cref="F60GhzBornBelowFoldPi2Inheritance.FoldPosition"/> and
///         <see cref="F62WStateBornBelowFoldPi2Inheritance.FoldPosition"/> reference.
///         F60 GHZ_3 and F62 W_3 both sit below this fold; F69 says the mix lifts above.</item>
///   <item><b>F60 sibling (GHZ baseline)</b>: GHZ_3 has CΨ(0) = 1/(2³−1) = 1/7 ≈ 0.143
///         (F60.CPsiAtZeroForGhz(3)); strictly below fold.</item>
///   <item><b>F62 sibling (W baseline)</b>: W_3 has CΨ(0) = 10/81 ≈ 0.1235
///         (F62.CPsiAtZeroForWState(3)); strictly below fold.</item>
///   <item><b>Mix lifts above (F69 finding)</b>: the two-parameter slice α·GHZ_3 +
///         √(1-α²)·W_3 admits an interior maximum that exceeds 1/4. The mechanism
///         is interference between GHZ and W via the L1-off √6·α·√(1-α²) cross
///         term, which has no equivalent at N ≥ 4.</item>
/// </list>
///
/// <para><b>Scope (N=3 only):</b> The same GHZ+W slice peaks at 0.167 (N=4),
/// 0.146 (N=5), 0.134 (N=6), 0.125 (N=7), 0.118 (N=8); all below 1/4. The
/// GHZ-purity contribution scales as 1/(2^N − 1) and collapses too fast for
/// the W-contribution to lift the peak back above the fold.</para>
///
/// <para><b>Full Dicke subspace:</b> Pair-CΨ has no non-product local maxima
/// on the permutation-symmetric Dicke sphere at any tested N ∈ {3..8}. The
/// only non-product stationary points are Dicke basis elements (max ≈ 0.123
/// at N=3, ≤ 1/12 for larger N) and the GHZ+W family optimum itself; all
/// saddles on the full sphere with escape Δpair-CΨ ≈ 0.68 on 1% c_2
/// perturbation at N=3.</para>
///
/// <para>Tier1Derived: closed form (sextic + numerical optimum at 25-digit
/// precision) verified algebraically. Pi2-Foundation anchoring is the typed
/// statement that F69's mixed state lifts above F60/F62's shared 1/4 fold.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F69 (line 1467) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>
/// (1/4 = a_3 fold reference) +
/// <c>compute/RCPsiSquared.Core/Symmetry/F60GhzBornBelowFoldPi2Inheritance.cs</c>
/// (GHZ sibling, pure-state below fold) +
/// <c>compute/RCPsiSquared.Core/Symmetry/F62WStateBornBelowFoldPi2Inheritance.cs</c>
/// (W-state sibling, pure-state below fold).</para></summary>
public sealed class F69GhzWSexticAboveFoldPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly F60GhzBornBelowFoldPi2Inheritance _f60;
    private readonly F62WStateBornBelowFoldPi2Inheritance _f62;

    /// <summary>The 1/4 fold position. Live from Pi2DyadicLadder a_3.
    /// Same anchor that F60 and F62 use as their FoldPosition.</summary>
    public double FoldPosition => _ladder.Term(3);

    /// <summary>Sextic coefficients for P(x) = 2900x⁶ − 8060x⁵ + 4211x⁴ + 3832x³ − 2428x² − 512x + 300,
    /// the integer-coefficient minimal polynomial for x = α²_opt. Irreducible over ℚ;
    /// optimum is a degree-6 algebraic number with no radical form.</summary>
    public IReadOnlyList<int> SexticCoefficients { get; } = new[] { 2900, -8060, 4211, 3832, -2428, -512, 300 };

    /// <summary>Numerical optimum α²_opt = 0.375420720711069 (25-digit precision; root of
    /// SexticCoefficients above the bilinear-apex fold).</summary>
    public const double AlphaSquaredOptimum = 0.375420720711069;

    /// <summary>Numerical optimum α_opt = 0.612715856422101 (= √α²_opt).</summary>
    public const double AlphaOptimum = 0.612715856422101;

    /// <summary>Numerical pair-CΨ(0) at the GHZ+W optimum: 0.320411541127025.
    /// Maximum value attained by any non-product permutation-symmetric state on N=3 qubits.</summary>
    public const double OptimumPairCpsi = 0.320411541127025;

    /// <summary>Lift ratio above the 1/4 fold: ≈ 1.281646. F69's 28% excess is
    /// the structural quantity; underlying value is OptimumPairCpsi / FoldPosition.</summary>
    public double LiftRatioAboveFold => OptimumPairCpsi / FoldPosition;

    /// <summary>True iff F69's optimum exceeds the 1/4 fold (the structural claim of F69).</summary>
    public bool OptimumExceedsFold => OptimumPairCpsi > FoldPosition;

    /// <summary>True iff the F60 GHZ baseline at N is strictly below the fold (sibling check).
    /// F69 requires both F60 and F62 below to make the "mix lifts above" finding non-trivial.</summary>
    public bool GhzBaselineBelowFold(int N) => _f60.CPsiAtZeroForGhz(N) < FoldPosition;

    /// <summary>True iff the F62 W baseline at N is strictly below the fold (sibling check).</summary>
    public bool WBaselineBelowFold(int N) => _f62.CPsiAtZeroForWState(N) < FoldPosition;

    /// <summary>The GHZ+W superposition mechanism is N=3-specific (verified by enumeration
    /// at N=4..8). The same family peaks at 0.167 (N=4), 0.146 (N=5), 0.134 (N=6),
    /// 0.125 (N=7), 0.118 (N=8); all below 1/4. F69's lift exists only at N=3.</summary>
    public bool LiftExistsAtN(int N) => N == 3;

    /// <summary>3-tangle τ_ABC at the optimum (N=3). The high value 0.799 indicates
    /// the optimum lies near the GHZ limit (τ_ABC = 1 for pure GHZ).</summary>
    public const double ThreeTangleAtOptimum = 0.799453;

    /// <summary>Pair concurrence C(A, B) at the optimum. Essentially zero (0.021),
    /// indicating no two-qubit entanglement; the entanglement is purely tripartite.</summary>
    public const double PairConcurrenceAtOptimum = 0.0210;

    public F69GhzWSexticAboveFoldPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        F60GhzBornBelowFoldPi2Inheritance f60,
        F62WStateBornBelowFoldPi2Inheritance f62)
        : base("F69 GHZ+W mix lifts pair-CΨ(0) above 1/4 fold at N=3 (sextic minimal polynomial, irreducible over ℚ); siblings F60 + F62 both born below fold",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F69 + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/F60GhzBornBelowFoldPi2Inheritance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/F62WStateBornBelowFoldPi2Inheritance.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _f60 = f60 ?? throw new ArgumentNullException(nameof(f60));
        _f62 = f62 ?? throw new ArgumentNullException(nameof(f62));
    }

    public override string DisplayName =>
        "F69 GHZ+W sextic above-fold (N=3) as Pi2-Foundation a_3 + F60 + F62 inheritance";

    public override string Summary =>
        $"|ψ(α)⟩ = α|GHZ_3⟩ + √(1-α²)|W_3⟩ optimum sits at pair-CΨ ≈ 0.3204 > 1/4 (lift ratio {LiftRatioAboveFold:F4}); sextic is irreducible over ℚ; N=3 only ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F69 sextic",
                summary: "P(x) = 2900x⁶ − 8060x⁵ + 4211x⁴ + 3832x³ − 2428x² − 512x + 300; irreducible over ℚ (sympy.Poly.is_irreducible = True); α²_opt is a degree-6 algebraic number with no radical form in nested square roots");
            yield return InspectableNode.RealScalar("FoldPosition (= a_3 = 1/4)", FoldPosition);
            yield return InspectableNode.RealScalar("AlphaSquaredOptimum", AlphaSquaredOptimum);
            yield return InspectableNode.RealScalar("OptimumPairCpsi", OptimumPairCpsi);
            yield return InspectableNode.RealScalar("LiftRatioAboveFold", LiftRatioAboveFold);
            yield return new InspectableNode("F60/F62 siblings below fold",
                summary: $"F60 GHZ_3: CΨ(0) = 1/(2³-1) = {_f60.CPsiAtZeroForGhz(3):F6} < 1/4; F62 W_3: CΨ(0) = 10/81 = {_f62.CPsiAtZeroForWState(3):F6} < 1/4. Both below; mix lifts above (mechanism: interference via L1-off √6·α·√(1-α²) cross term).");
            yield return new InspectableNode("3-tangle near-GHZ at optimum",
                summary: $"τ_ABC = {ThreeTangleAtOptimum} (close to 1 = pure GHZ); pair concurrence C(A,B) = {PairConcurrenceAtOptimum} ≈ 0 (no bipartite entanglement, all tripartite)");
            yield return new InspectableNode("N=3 unique",
                summary: "GHZ+W slice peaks at N=4..8: 0.167, 0.146, 0.134, 0.125, 0.118; all below 1/4 fold. The lift exists only at N=3 because GHZ-purity scales as 1/(2^N − 1) and collapses for N ≥ 4.");
            yield return new InspectableNode("full Dicke saddles",
                summary: "Pair-CΨ has no non-product local maxima on the permutation-symmetric Dicke sphere at N ∈ {3..8}; only Dicke basis elements (max ≈ 0.123 at N=3) and the GHZ+W family optimum itself; all saddles with escape Δpair-CΨ ≈ 0.68 on 1% c_2 perturbation at N=3");
        }
    }
}
