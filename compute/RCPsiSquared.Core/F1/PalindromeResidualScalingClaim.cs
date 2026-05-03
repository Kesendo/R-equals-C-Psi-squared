using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;

namespace RCPsiSquared.Core.F1;

/// <summary>Closed-form Frobenius-norm scaling of the F1-palindrome residual M for
/// non-truly Hamiltonian classes (Tier 1 derived; verified bit-exact across chain, ring,
/// star, K_N at N = 4, 5).
///
/// <code>
///     ‖M(N, G)‖²_F  =  c_H · F(N, G)
///
///     main class          F(N, G) = B(G) · 4^(N−2)            (B = bond count)
///     single-body class   F(N, G) = (D2(G) / 2) · 4^(N−2)     (D2 = Σ_i deg_G(i)²)
/// </code>
///
/// <para>For chains, B = N − 1 and D2 = 4N − 6, giving F_main(N) = (N−1)·4^(N−2) and
/// F_single(N) = (2N−3)·4^(N−2). The adjacent-N ratio ‖M(N+1)‖²/‖M(N)‖² is a function of
/// N alone, independent of γ, J, topology, and equals 4·N/(N−1) (main) or
/// 4·(2N−1)/(2N−3) (single-body).</para>
///
/// <para>Source: <c>experiments/OPERATOR_RIGIDITY_ACROSS_CUSP.md</c>. The scaling
/// computation lives in <see cref="PalindromeResidualScaling"/>.</para>
/// </summary>
public sealed class PalindromeResidualScalingClaim : Claim
{
    public int N { get; }
    public HamiltonianClass HClass { get; }

    /// <summary>Optional bond count (graph-aware). When null, the chain default B = N − 1
    /// is used.</summary>
    public int? BondCount { get; }

    /// <summary>Optional degree-squared sum (graph-aware). When null, the chain default
    /// D2 = 4N − 6 is used.</summary>
    public int? DegreeSquaredSum { get; }

    public double Factor => (BondCount, DegreeSquaredSum) switch
    {
        (null, null) => PalindromeResidualScaling.FactorChain(N, HClass),
        ({ } b, { } d2) => PalindromeResidualScaling.FactorFromGraph(N, b, d2, HClass),
        _ => throw new InvalidOperationException(
            "BondCount and DegreeSquaredSum must both be provided for graph-aware scaling, or both omitted for the chain default."),
    };

    public double AdjacentRatio => PalindromeResidualScaling.AdjacentRatio(N, HClass);

    public PalindromeResidualScalingClaim(int N, HamiltonianClass hClass,
        int? bondCount = null, int? degreeSquaredSum = null)
        : base($"F1 residual ‖M‖² scaling, {hClass} class",
               Tier.Tier1Derived,
               "RCPsiSquared.Core.Lindblad.PalindromeResidualScaling + experiments/OPERATOR_RIGIDITY_ACROSS_CUSP.md")
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), $"N must be ≥ 2; got {N}");
        this.N = N;
        HClass = hClass;
        BondCount = bondCount;
        DegreeSquaredSum = degreeSquaredSum;
    }

    public override string DisplayName =>
        $"‖M‖² scaling ({HClass}, N={N}): F = {Factor:G6}, adjacent ratio {AdjacentRatio:G6}";

    public override string Summary => HClass switch
    {
        HamiltonianClass.Main =>
            $"main class: ‖M‖² = c_H · B · 4^(N−2); chain B = N−1 → F({N}) = {Factor:G6}; ratio ‖M({N + 1})‖²/‖M({N})‖² = {AdjacentRatio:G6}",
        HamiltonianClass.SingleBody =>
            $"single-body class: ‖M‖² = c_H · (D2/2) · 4^(N−2); chain D2 = 4N−6 → F({N}) = {Factor:G6}; ratio ‖M({N + 1})‖²/‖M({N})‖² = {AdjacentRatio:G6}",
        _ => $"unknown Hamiltonian class {HClass}",
    };

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("N", N);
            yield return new InspectableNode("Hamiltonian class", summary: HClass.ToString());
            yield return InspectableNode.RealScalar("F(N, G)", Factor, "G6");
            yield return InspectableNode.RealScalar("adjacent ratio ‖M(N+1)‖²/‖M(N)‖²", AdjacentRatio, "G6");
            if (BondCount is { } b) yield return InspectableNode.RealScalar("bond count B", b);
            if (DegreeSquaredSum is { } d2) yield return InspectableNode.RealScalar("degree² sum D2", d2);
            yield return new InspectableNode("verification",
                summary: "bit-exact at N = 4, 5 across chain, ring, star, K_N");
        }
    }
}
