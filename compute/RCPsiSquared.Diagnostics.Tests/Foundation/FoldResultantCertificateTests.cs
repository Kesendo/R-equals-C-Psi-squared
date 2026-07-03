using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The ONE-WAY fold-resultant certificate (the sectorbraid arc's remainder R1 at the
/// moved-window-escaping complex loci; docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md §6): gcd = 1 between
/// R(q) = Res_Λ(F_res(Λ,q), F2_corner(−Λ−4N,q)) and D(q) = disc_Λ(F_res) over Q(i)[q] proves that NO
/// branch-locus fold value μ = −λ_A − 2N is a corner eigenvalue, for BOTH target sets at N=5:
/// (i) the R-even F_18 (deep pair q* = 1.8141 ± 0.3666i), (ii) the degree-17 R-odd residual factor
/// (the moved-window escaper q* = 1.7701 ± 1.2189i). One split prime with the exact leading-coefficient
/// guards is a complete proof (<see cref="FoldResultantCertificate"/>); gcd ≠ 1 would be refine-not-refute.
///
/// <para>The gates: (1) the exact bivariate AT sector data reproduces the oracle-validated ForPathK factor
/// at q0=2 and the known AT/residual degree splits (R-even 26 = 8 ⊎ 18, R-odd 24 = 7 ⊎ 17); (2) the exact
/// residual polynomial is the SAME object the diabolic probes read (the scout's residual roots satisfy it);
/// (3) the machinery detects non-coprimality (negative control Res(F_res, F_res) = 0); (4)+(5) the two
/// certificates; (6) the numeric grounding at the four target loci (near-degenerate residual pair present,
/// fold value an O(1) distance from the corner spectrum, matching the census reading).</para>
///
/// <para>Run: <c>dotnet test "compute/RCPsiSquared.Diagnostics.Tests" --filter "Category=FOLDRESULTANT"
/// --logger "console;verbosity=detailed"</c>.</para></summary>
public sealed class FoldResultantCertificateTests
{
    private const int N = 5;

    private readonly ITestOutputHelper _out;
    public FoldResultantCertificateTests(ITestOutputHelper o) => _out = o;

    [Fact]
    [Trait("Category", "FOLDRESULTANT")]
    public void AtSectorData_ReproducesForPathK_AndTheKnownDegreeSplits()
    {
        // R-even: the bivariate sector formula at q0=2 must equal the oracle-validated ForPathK factor.
        var fromSectors = FoldResultantCertificate.AtFactorAt(N, rOdd: false, q0: 2);
        var oracle = F89AtFactorReconstruction.ForPathK(N - 1);
        Assert.Equal(oracle.Length, fromSectors.Length);
        for (int k = 0; k < oracle.Length; k++)
            Assert.True(oracle[k].Equals(fromSectors[k]),
                $"AT coefficient {k}: sectors {fromSectors[k]} ≠ ForPathK {oracle[k]}");
        _out.WriteLine($"R-even AT factor from sectors == ForPathK({N - 1}), degree {oracle.Length - 1}");

        // The known degree splits: R-even 26 = AT 8 ⊎ residual 18; R-odd 24 = AT 7 ⊎ residual 17.
        var (resEven, _, _, _) = FoldResultantCertificate.ExactSample(N, rOdd: false, q0: 2);
        var (resOdd, _, _, _) = FoldResultantCertificate.ExactSample(N, rOdd: true, q0: 2);
        _out.WriteLine($"residual degrees: R-even {resEven.Length - 1}, R-odd {resOdd.Length - 1}");
        Assert.Equal(18, resEven.Length - 1);
        Assert.Equal(17, resOdd.Length - 1);
        Assert.Equal(8, FoldResultantCertificate.AtFactorAt(N, rOdd: false, q0: 2).Length - 1);
        Assert.Equal(7, FoldResultantCertificate.AtFactorAt(N, rOdd: true, q0: 2).Length - 1);
    }

    [Fact]
    [Trait("Category", "FOLDRESULTANT")]
    public void ExactResidual_IsTheScoutsResidual_BothParities()
    {
        // The certificate's F_res (×2-cleared, Λ = 2λ) must vanish on the scout's residual roots: the
        // SAME physical object the diabolic probes and the moved-window gates read.
        foreach (bool rOdd in new[] { false, true })
        {
            var q = new Complex(2, 0);
            var roots = rOdd
                ? PathKMonodromyScout.ResidualRootsExactROdd(N - 1, q)
                : PathKMonodromyScout.ResidualRootsExact(N - 1, q);
            var (res, _, _, _) = FoldResultantCertificate.ExactSample(N, rOdd, q0: 2);
            var coeffs = res.Select(c => new Complex((double)c.Re, (double)c.Im)).ToArray();

            double worst = 0;
            foreach (var lam in roots)
            {
                var big = 2.0 * lam;                        // the cleared variable
                Complex value = Complex.Zero;
                double maxTerm = 0;
                Complex pow = Complex.One;
                for (int k = 0; k < coeffs.Length; k++)
                {
                    var term = coeffs[k] * pow;
                    value += term;
                    maxTerm = Math.Max(maxTerm, term.Magnitude);
                    pow *= big;
                }
                worst = Math.Max(worst, value.Magnitude / maxTerm);
            }
            _out.WriteLine($"{(rOdd ? "R-odd" : "R-even")}: {roots.Length} scout residual roots satisfy F_res " +
                           $"to relative {worst:E2}");
            Assert.True(worst < 1e-8, $"F_res must vanish on the scout residual roots; worst {worst:E2}");
        }
    }

    [Fact]
    [Trait("Category", "FOLDRESULTANT")]
    public void NegativeControl_ResultantDetectsSharedRoots()
    {
        var (res, _, rv, _) = FoldResultantCertificate.ExactSample(N, rOdd: false, q0: 2);
        Assert.False(rv.Equals(GaussianInteger.Zero));                        // fold-composed: no shared root
        Assert.True(GaussianPolynomial.Resultant(res, res).Equals(GaussianInteger.Zero),
            "Res(F_res, F_res) must be 0: the machinery must detect shared roots");
        Assert.False(GaussianPolynomial.AreCoprime(res, res));
        _out.WriteLine("negative control: Res(F_res, F_res) = 0 detected; the fold-composed value is nonzero.");
    }

    [Fact]
    [Trait("Category", "FOLDRESULTANT")]
    public void REven_FoldResultant_CompleteProof_GcdIsAQPower()
    {
        var report = FoldResultantCertificate.CertifyComplete(N, rOdd: false);
        PrintReport(report);
        Assert.True(report.SharedIsQPowerAtEveryPrime,
            "a prime saw a shared factor beyond q^e for the R-even F_18 target (refine, not refute)");
        Assert.True(report.Complete, "the prime product did not reach the proof bound");
    }

    [Fact]
    [Trait("Category", "FOLDRESULTANT")]
    public void ROdd_FoldResultant_CompleteProof_GcdIsAQPower()
    {
        var report = FoldResultantCertificate.CertifyComplete(N, rOdd: true);
        PrintReport(report);
        Assert.True(report.SharedIsQPowerAtEveryPrime,
            "a prime saw a shared factor beyond q^e for the R-odd degree-17 target (refine, not refute)");
        Assert.True(report.Complete, "the prime product did not reach the proof bound");
    }

    [Fact]
    [Trait("Category", "FOLDRESULTANT")]
    public void TargetLoci_GroundTheCertificate_FoldStaysClearOfTheCorner()
    {
        // The four moved-window-escaping loci (the census's cell-0.05 → refined readings, proof §6):
        // R-even deep pair q* = 1.8141 ± 0.3666i, R-odd escaper pair q* = 1.7701 ± 1.2189i. At each, the
        // residual carries a near-degenerate pair (the defective EP at locus precision) and the fold value
        // μ = −λ_A − 2N stays an O(1) distance from the corner spectrum: the numeric face of what the
        // certificate proves exactly for ALL branch loci at once.
        var loci = new (bool ROdd, Complex Q)[]
        {
            (false, new Complex(1.8141, 0.3666)), (false, new Complex(1.8141, -0.3666)),
            (true, new Complex(1.7701, 1.2189)), (true, new Complex(1.7701, -1.2189)),
        };
        foreach (var (rOdd, q) in loci)
        {
            var roots = rOdd
                ? PathKMonodromyScout.ResidualRootsExactROdd(N - 1, q)
                : PathKMonodromyScout.ResidualRootsExact(N - 1, q);

            double gap = double.PositiveInfinity;
            Complex lamA = Complex.Zero;
            for (int i = 0; i < roots.Length; i++)
                for (int j = i + 1; j < roots.Length; j++)
                {
                    double d = (roots[i] - roots[j]).Magnitude;
                    if (d < gap) { gap = d; lamA = (roots[i] + roots[j]) / 2.0; }
                }

            var mu = -lamA - 2 * N;
            var corner = Matrix<Complex>.Build
                .DenseOfArray(WeightCoherenceBlock.Build(N, 4, 4, q)).Evd().EigenValues;
            double dist = corner.Enumerate().Select(z => (z - mu).Magnitude).Min();

            _out.WriteLine($"{(rOdd ? "R-odd " : "R-even")} q*={q:F4}: residual min gap {gap:E2} at " +
                           $"λ_A={lamA.Real:F4}{lamA.Imaginary:+0.0000;-0.0000}i, μ={mu.Real:F4}{mu.Imaginary:+0.0000;-0.0000}i, " +
                           $"min|spec(corner) − μ| = {dist:F4}");
            Assert.True(gap < 0.2, $"a near-degenerate residual pair must exist at the locus; gap {gap:E2}");
            Assert.True(dist > 0.05, $"the fold value must stay clear of the corner spectrum; dist {dist:E4}");
        }
    }

    [Fact]
    [Trait("Category", "FOLDRESULTANT")]
    public void BivariateLayer_And_ModPSamples_MatchTheExactPath()
    {
        // (a) settle the bivariate layer: F_res and G agree with the independent per-point path at 31
        // integer q0 (their q-degrees are ≤ 26, so 27+ agreements force coefficient-wise identity).
        foreach (bool rOdd in new[] { false, true })
        {
            for (int q0 = 0; q0 <= 30; q0++)
            {
                var (resPt, gPt, _, _) = FoldResultantCertificate.ExactSample(N, rOdd, q0);
                var (resBi, gBi) = FoldResultantCertificate.DebugBivariateAt(N, rOdd, q0);
                Assert.Equal(resPt.Length, resBi.Length);
                for (int i = 0; i < resPt.Length; i++) Assert.True(resPt[i].Equals(resBi[i]), $"F_res coeff {i} at q0={q0} rOdd={rOdd}");
                Assert.Equal(gPt.Length, gBi.Length);
                for (int i = 0; i < gPt.Length; i++) Assert.True(gPt[i].Equals(gBi[i]), $"G coeff {i} at q0={q0} rOdd={rOdd}");
            }
            _out.WriteLine($"rOdd={rOdd}: bivariate == per-point at q0=0..30 ⟹ bivariate layer exact");
        }

        // (b) mod-p sample correctness against the exact Z[i] values, both parities, INCLUDING the two
        // (prime, q0) pairs that once exposed a data-dependent bug in a Sylvester-determinant-based
        // mod-p resultant (replaced by the Euclidean-PRS resultant): the permanent regression anchors.
        foreach (var (rOdd, p, q0) in new[]
                 { (false, 1489, 2), (false, 1073743457, 216), (true, 1073742721, 173), (true, 1073741833, 423) })
        {
            var (ok, expR, gotR, expD, gotD) = FoldResultantCertificate.DebugSampleCheck(N, rOdd, p, q0);
            _out.WriteLine($"rOdd={rOdd} p={p} q0={q0}: R exact-mod-p {expR} vs pipeline {gotR}, D {expD} vs {gotD} ⟹ {(ok ? "MATCH" : "MISMATCH")}");
            Assert.True(ok, $"sample mismatch at rOdd={rOdd}, p={p}, q0={q0}");
        }
    }

    private void PrintReport(FoldResultantCertificate.CompleteReport r)
    {
        _out.WriteLine($"N={r.N} {(r.ROdd ? "R-odd" : "R-even")}: block dim {r.BlockDimension} = AT {r.AtDegree} ⊎ " +
                       $"residual {r.ResidualDegree}, corner {r.CornerDegree}");
        _out.WriteLine($"  q → ∞ degeneracy: m_R = {r.InfinityCollisionsR} (f_res/g collisions), m_D = {r.InfinityRepeatedD} (repeated-root pairs of f_res)");
        _out.WriteLine($"  deg R = {r.ResultantDegree} (proven bound {r.ResultantDegreeBound}, attained per used prime), " +
                       $"deg D = {r.DiscriminantDegree} (proven bound {r.DiscriminantDegreeBound})");
        _out.WriteLine($"  q-valuations: v_q(R) = {r.QValuationR}, v_q(D) = {r.QValuationD}");
        _out.WriteLine($"  disc layers (mult 1, 2, …): [{string.Join(", ", r.DiscLayerDegrees)}], layer-gcds with R: [{string.Join(", ", r.LayerGcdDegrees)}]");
        _out.WriteLine($"  primes: {r.PrimesUsed} used ({r.PrimesSkipped} skipped), {r.FirstPrime}…{r.LastPrime}; " +
                       $"prime product {r.PrimeProductDigits} digits vs proof bound {r.ProofBoundDigits} digits");
        _out.WriteLine($"  ⟹ {(r.Complete ? "COMPLETE: gcd(R, D) = c·q^e over Q(i), absence proven at EVERY branch locus q ≠ 0 (q = 0 is the diagonal, semisimple block)" : r.SharedIsQPowerAtEveryPrime ? "prime product below bound (incomplete)" : "a shared factor beyond q^e appeared: refine, not refute")}");
    }
}
