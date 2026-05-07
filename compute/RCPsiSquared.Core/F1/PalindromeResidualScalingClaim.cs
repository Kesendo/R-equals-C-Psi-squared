using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;

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
public sealed class PalindromeResidualScalingClaim : Claim, IDriftCheckable
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

    public DriftReport Verify()
    {
        // Operational verification: anchor the Hamiltonian-dependent constant c_H from
        // a small system (N=2, 16x16 Pauli-basis) and use the closed form to predict
        // ‖M‖² at the registered N. Then build the actual Liouvillian at N, compute
        // ‖M(N)‖² independently, and assert it matches the prediction. Same H structure
        // (XX+YZ chain, J=1) at both N=2 and N to keep c_H consistent.
        //
        // Limitations: SingleBody verification needs a single-body H builder we do not
        // currently expose; graph-aware mode (BondCount + DegreeSquaredSum supplied)
        // would need a custom topology builder. Both fall back to a labelled OK report
        // until those builders land. N is capped at 5 because PalindromeResidual.Build
        // internally constructs a 4^N x 4^N transform; at N=5 that is 1024x1024 (~16MB
        // dense complex), at N=6 it is 4096x4096 (~256MB) which stalls under default xunit
        // memory budgets.
        if (BondCount is not null || DegreeSquaredSum is not null)
            return new DriftReport(
                ClaimType: typeof(PalindromeResidualScalingClaim),
                IsDrift: false,
                Description: $"graph-aware mode skips operational drift check at (N={N}, {HClass}).",
                Magnitude: null);

        if (HClass != HamiltonianClass.Main)
            return new DriftReport(
                ClaimType: typeof(PalindromeResidualScalingClaim),
                IsDrift: false,
                Description: $"SingleBody operational verification not yet implemented at N={N}; needs a single-body H builder.",
                Magnitude: null);

        if (N < 2 || N > 5)
            return new DriftReport(
                ClaimType: typeof(PalindromeResidualScalingClaim),
                IsDrift: false,
                Description: $"Operational verification limited to N=2..5 (4^N Pauli-basis transform); N={N} skipped.",
                Magnitude: null);

        double mNormSquaredAt(int n)
        {
            var bonds = new Bond[n - 1];
            for (int i = 0; i < n - 1; i++) bonds[i] = new Bond(i, i + 1, 1.0);
            var terms = new (PauliLetter, PauliLetter, Complex)[]
            {
                (PauliLetter.X, PauliLetter.X, Complex.One),
                (PauliLetter.Y, PauliLetter.Z, Complex.One),
            };
            var H = PauliHamiltonian.Bilinear(n, bonds, terms).ToMatrix();
            const double gammaZ = 0.05;
            var gammaList = Enumerable.Repeat(gammaZ, n).ToArray();
            var L = PauliDephasingDissipator.BuildZ(H, gammaList);
            var M = PalindromeResidual.Build(L, n, n * gammaZ, PauliLetter.Z);
            double frob = M.FrobeniusNorm();
            return frob * frob;
        }

        // Anchor c_H at N=2 (one bond), where ‖M(2)‖² = c_H · F(2) and F(2) = 1 for Main
        // chain (B=1, 4^0=1). So c_H = ‖M(2)‖² directly.
        double cHObserved = mNormSquaredAt(2);
        if (cHObserved < 1e-12)
            return new DriftReport(
                ClaimType: typeof(PalindromeResidualScalingClaim),
                IsDrift: false,
                Description: $"c_H anchor at N=2 is below 1e-12 (Hamiltonian behaved truly under F1); operational verification not applicable.",
                Magnitude: null);

        if (N == 2)
            // c_H comes from N=2 itself; verification is reflexive at the anchor point.
            return new DriftReport(
                ClaimType: typeof(PalindromeResidualScalingClaim),
                IsDrift: false,
                Description: $"F73 operational anchor at N=2 (H=XX+YZ chain): ‖M(2)‖² = {cHObserved:G10} = c_H by definition.",
                Magnitude: 0.0);

        double observed = mNormSquaredAt(N);
        double predicted = cHObserved * Factor;
        double diff = Math.Abs(observed - predicted);
        // Tolerance: allow for accumulated floating-point error across the 4^N basis
        // transform. Empirically the error is 12+ orders of magnitude below predicted
        // for N=3..5; 1e-6 relative is generous.
        double tolerance = 1e-6 * predicted;
        bool drift = diff > tolerance;

        return new DriftReport(
            ClaimType: typeof(PalindromeResidualScalingClaim),
            IsDrift: drift,
            Description: drift
                ? $"F73 operational drift at (N={N}, {HClass}, H=XX+YZ chain): observed ‖M({N})‖² = {observed:G10}, predicted c_H · F = {cHObserved:G6} · {Factor:G6} = {predicted:G10}, |delta| = {diff:G3}."
                : $"F73 operational OK at (N={N}, {HClass}, H=XX+YZ chain): observed {observed:G10} matches predicted c_H · F = {predicted:G10} to {diff:G3} (tolerance {tolerance:G3}).",
            Magnitude: diff);
    }
}
