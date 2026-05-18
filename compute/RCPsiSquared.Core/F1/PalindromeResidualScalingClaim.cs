using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.F1;

/// <summary>Closed-form Frobenius-norm scaling of the F1-palindrome residual M for
/// non-truly Hamiltonian classes (Tier 1 derived; verified bit-exact across arbitrary
/// connected and disconnected graphs at N = 4..7, including random Erdős-Rényi and
/// weighted edges).
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
/// <para><b>γ-independence (uniform AND non-uniform γ, closed 2026-05-18).</b> F(N, G)
/// holds verbatim for arbitrary per-site γ_l: the dissipator-block residual
/// M_D = Π·L_D·Π⁻¹ + L_D + 2Σγ·I vanishes per Pauli string, so ‖M‖²_F = ‖M_H‖²_F and
/// M_H carries no γ dependence by construction. Closes the earlier F1 OpenQuestion that
/// conjectured a Σγ_l² replacement of (Σγ)² in F(N, G); no formula change required.
/// Anchor: <c>docs/proofs/PROOF_F1_NONUNIFORM_GAMMA.md</c>;
/// verification: <c>simulations/_f1_nonuniform_gamma_verify.py</c>.</para>
///
/// <para><b>General-topology universality (closed 2026-05-18).</b> The (B, D2)
/// parameterisation extends bit-exactly to disconnected graphs (B and D2 sum across
/// components), weighted edges (B → Σ_b J²_b), random connected Erdős-Rényi graphs at
/// N=5, 6, and the F1 palindromic-pairing identity holds at N=7 across chain, ring,
/// star, and K_4 + disjoint-3-chain via the
/// <see cref="BlockSpectrum.LiouvillianBlockSpectrum"/> dogfood path. The analytic
/// content was already in <c>docs/proofs/PROOF_CROSS_TERM_FORMULA.md</c> Lemma 3 +
/// Corollary (bond-disjointness independent of connectivity); the verification record
/// lives in the Tier-2 <see cref="F1GeneralTopologyVerifiedClaim"/>. Synthesis proof:
/// <c>docs/proofs/PROOF_F1_GENERAL_TOPOLOGY.md</c>; verification: Python script
/// <c>simulations/_f1_general_topology_verify.py</c> + C# test class
/// <c>compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN7BlockSpectrumTests.cs</c>.</para>
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

    private DriftReport? _cachedReport;

    public DriftReport Verify() => _cachedReport ??= ComputeVerify();

    private DriftReport ComputeVerify()
    {
        // Operational verification: anchor the Hamiltonian-dependent constant c_H from
        // a small system (N=2, 16x16 Pauli-basis) and use the closed form to predict
        // ‖M‖² at the registered N. Then build the actual Liouvillian at N, compute
        // ‖M(N)‖² independently, and assert it matches the prediction. Same H structure
        // (XX+YZ chain, J=1) at both N=2 and N to keep c_H consistent.
        //
        // Result is memoised on the claim instance: the Hamiltonian, gamma, and N are
        // immutable, so subsequent Verify() calls return the cached DriftReport in O(1)
        // instead of paying the ~14s Liouvillian-build cost again.
        //
        // Limitations: SingleBody verification needs a single-body H builder we do not
        // currently expose; graph-aware mode (BondCount + DegreeSquaredSum supplied)
        // would need a custom topology builder. Both fall back to a labelled OK report
        // until those builders land. N is capped at 5 because PalindromeResidual.Build
        // internally constructs a 4^N x 4^N transform; at N=5 that is 1024x1024 (~16MB
        // dense complex), at N=6 it is 4096x4096 (~256MB) which stalls under default xunit
        // memory budgets.
        //
        // For graph-aware operational verification at N=5, see
        // compute/RCPsiSquared.Core.Tests/F1/F1GeneralTopologyN7BlockSpectrumTests.cs
        // (4 tests: chain, ring, star, triangle + disjoint bond), plus the Python
        // simulations/_f1_general_topology_verify.py for N=5, 6 named/random/disconnected/
        // weighted coverage. The F1 palindromic-pairing identity at N=7 is exercised in
        // the same test class via LiouvillianBlockSpectrum.ComputeSpectrumPerBlock.
        if (BondCount is not null || DegreeSquaredSum is not null)
            return new DriftReport(
                ClaimType: typeof(PalindromeResidualScalingClaim),
                IsDrift: false,
                Description: $"graph-aware mode skips operational drift check at (N={N}, {HClass}); see F1GeneralTopologyN7BlockSpectrumTests for graph-aware verification.",
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
                Description: $"PalindromeResidualScaling operational anchor at N=2 (H=XX+YZ chain): ‖M(2)‖² = {cHObserved:G10} = c_H by definition.",
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
                ? $"PalindromeResidualScaling operational drift at (N={N}, {HClass}, H=XX+YZ chain): observed ‖M({N})‖² = {observed:G10}, predicted c_H · F = {cHObserved:G6} · {Factor:G6} = {predicted:G10}, |delta| = {diff:G3}."
                : $"PalindromeResidualScaling operational OK at (N={N}, {HClass}, H=XX+YZ chain): observed {observed:G10} matches predicted c_H · F = {predicted:G10} to {diff:G3} (tolerance {tolerance:G3}).",
            Magnitude: diff);
    }
}
