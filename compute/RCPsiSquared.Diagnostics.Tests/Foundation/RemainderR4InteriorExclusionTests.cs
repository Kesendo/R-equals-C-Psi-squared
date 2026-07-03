using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Remainder R4 (the exclusion half of diamond membership, PROOF_CODIM1_BY_ADDITIVITY §7): the
/// non-member interior cores do not carry the braid. The four interior non-members (1,1),(4,4),(1,4),(4,1)
/// reduce by the Klein full flip (spec(4,4)=spec(1,1), spec(4,1)=spec(1,4)) and the F89d bra/ket fold
/// (spec(1,4)(q̄) = −conj spec(1,1)(q) − 2N) to the single core (1,1) and its two shared values
/// {λ_A, μ = −λ_A − 2N}. (1,1)×μ IS R1's corner certificate (the corner (p_c+1,p_c+1) = (4,4) = (1,1) by
/// Klein); (1,1)×λ_A is the ONE new certificate here, the identity-composition sibling of the fold-resultant
/// certificate: gcd(Res_Λ(F_res, F2_{(1,1)}(Λ)), disc_Λ(F_res)) over Q(i) is a pure q-power, so no (1,2)
/// residual root (the braid λ_A) is a (1,1) eigenvalue at any branch locus q ≠ 0. With R1's corner + these
/// symmetries this closes the complete N=5 interior-core shared-λ exclusion, promoting the N=5 verdict.
///
/// <para>Run: <c>dotnet test "compute/RCPsiSquared.Diagnostics.Tests" --filter "Category=R4INTERIOR"
/// --logger "console;verbosity=detailed"</c>.</para></summary>
public sealed class RemainderR4InteriorExclusionTests
{
    private readonly ITestOutputHelper _out;
    public RemainderR4InteriorExclusionTests(ITestOutputHelper o) => _out = o;

    // ── The ONE new certificate: (1,1) does not carry the braid eigenvalue λ_A (identity composition). ──
    // The analytic degree bound rBound = resDeg·targetDeg − mR is TIGHT for the corner-fold (deg_q R = 422)
    // but LOOSE here (true deg_q R = 412 < 422; R-odd 384 < 394): the identity leading forms share a common
    // factor beyond the isolated-root collisions mR counts. The engine therefore guards on the CERTIFIED
    // empirical degree trueDegR = max_p DegP(R mod p): sampling more distinct split primes than can divide
    // lc_q(R) (∏ N(π_p) ≤ N(lc) ≤ ‖R‖², the Hadamard bound ⟹ ≤ 2·log₂‖R‖/30 of the p ≥ 2³⁰ primes)
    // certifies the max IS deg_q R exactly; used primes attain it, so π_p ∤ lc_q(R) and the Gauss lift is
    // unchanged. The pins below (412/384, sampled > bound) hold the certificate trail visible.
    [Fact]
    [Trait("Category", "R4INTERIOR")]
    public void Interior11_ExcludesTheBraidEigenvalue_REven_CompleteCertificate()
    {
        var r = FoldResultantCertificate.CertifyBlockExclusion(
            5, rOdd: false, tWKet: 1, tWBra: 1, GaussianInteger.One, GaussianInteger.Zero);
        _out.WriteLine($"(1,1)×λ_A R-even: complete={r.Complete}, resDeg={r.ResidualDegree}, targetDeg={r.CornerDegree}, "
            + $"rBound={r.ResultantDegreeBound}, trueDegR={r.ResultantDegree}, primesUsed={r.PrimesUsed}, "
            + $"sampled={r.PrimesSampled} (lc-divisor bound {r.LcDivisorBound}), proofDigits={r.ProofBoundDigits}");
        Assert.Equal(412, r.ResultantDegree);                       // the certified true degree, 10 below the bound 422
        Assert.True(r.PrimesSampled > r.LcDivisorBound, "the degree certificate needs more sampled primes than can divide lc_q(R)");
        Assert.True(r.Complete, "the (1,1) core must exclude every (1,2) R-even residual root (the braid λ_A)");
    }

    [Fact]
    [Trait("Category", "R4INTERIOR")]
    public void Interior11_ExcludesTheBraidEigenvalue_ROdd_CompleteCertificate()
    {
        var r = FoldResultantCertificate.CertifyBlockExclusion(
            5, rOdd: true, tWKet: 1, tWBra: 1, GaussianInteger.One, GaussianInteger.Zero);
        _out.WriteLine($"(1,1)×λ_A R-odd: complete={r.Complete}, resDeg={r.ResidualDegree}, trueDegR={r.ResultantDegree}, "
            + $"primesUsed={r.PrimesUsed}, sampled={r.PrimesSampled} (lc-divisor bound {r.LcDivisorBound})");
        Assert.Equal(384, r.ResultantDegree);                       // the certified true degree, 10 below the bound 394
        Assert.True(r.PrimesSampled > r.LcDivisorBound, "the degree certificate needs more sampled primes than can divide lc_q(R)");
        Assert.True(r.Complete, "the (1,1) core must exclude every (1,2) R-odd residual root");
    }

    // ── The symmetry reduction: interior four → (1,1) via Klein + F89d (so ONE certificate covers all four). ──
    [Fact]
    [Trait("Category", "R4INTERIOR")]
    public void InteriorFour_ReduceTo_Core11_ByKleinAndF89d()
    {
        const int N = 5, twoN = 10;
        foreach (var q in new[] { new Complex(0.620878, 0), new Complex(1.8141, 0.3666), new Complex(1.7701, 1.2189) })
        {
            var qb = Complex.Conjugate(q);
            var s11 = Spec(N, 1, 1, q);
            Assert.True(MultisetDist(Spec(N, 4, 4, q), s11) < 1e-9, $"Klein (4,4)=(1,1) at q={q}");
            Assert.True(MultisetDist(Spec(N, 4, 1, q), Spec(N, 1, 4, q)) < 1e-9, $"Klein (4,1)=(1,4) at q={q}");
            Assert.True(MultisetDist(Spec(N, 1, 4, qb), s11.Select(z => -Complex.Conjugate(z) - twoN).ToArray()) < 1e-9,
                $"F89d (1,4)(q̄) = −conj spec(1,1)(q) − 2N at q={q}");
        }
    }

    // ── Grounding: λ_A ∉ spec(1,1) with a clear gap at every N=5 locus (the 2nd-locus 0.208 the proof flagged
    //    as census-only is now the certificate's job); the fold-partner μ is R1's corner, already certified. ──
    [Fact]
    [Trait("Category", "R4INTERIOR")]
    public void BraidEigenvalue_NotIn_Core11Spectrum_AtEveryLocus()
    {
        const int N = 5;
        foreach (var (q, guess) in new[] {
            (new Complex(0.620878, 0), new Complex(-4.6189, 0)),
            (new Complex(1.077615, 0), new Complex(-3.7917, 0)),
            (new Complex(0.9938, 0.1825), new Complex(-4.7121, -0.8242)),
            (new Complex(1.8141, 0.3666), new Complex(-6.1586, -9.0359)),
            (new Complex(1.7701, 1.2189), new Complex(-3.7562, -0.5065)) })
        {
            var s12 = Spec(N, 1, 2, q);
            var s11 = Spec(N, 1, 1, q);
            var lamA = s12.OrderBy(z => (z - guess).Magnitude).First();
            double d = s11.Min(z => (lamA - z).Magnitude);
            _out.WriteLine($"q={q}: dist(λ_A={lamA:F3}, spec(1,1)) = {d:F4}");
            Assert.True(d > 0.05, $"λ_A must not be a (1,1) eigenvalue at q={q}; dist {d:F4}");
        }
    }

    // ── The propagation gate: the complete N=5 interior-four exclusion, assembled from its parts. ──
    // At every branch locus q* ≠ 0 (both parities): cert A ((1,1)×λ_A, identity composition; the two facts
    // above) proves NO residual root is a (1,1) eigenvalue; cert B (R1's corner-fold, CertifyComplete, gate
    // FOLDRESULTANT) proves NO fold image −λ−2N is a (4,4) eigenvalue; the Klein full flip spec(4,4) =
    // spec(1,1) (same q, q-independent unitary) makes A+B cover BOTH shared values {λ_A, μ = −λ_A−2N} on
    // BOTH diagonal cores. For the off-diagonal pair (1,4)/(4,1): composing F89d (spec(1,4)(q̄) =
    // −conj spec(1,1)(q) − 2N) with the pencil reality conj L(q) = L(−q̄) (A real, C entrywise imaginary)
    // cancels the two conjugations HOLOMORPHICALLY, spec(1,4)(q) = −spec(1,1)(−q) − 2N, and the bipartite
    // gauge (the chain's alternating-sign diagonal flips q, so spec(1,1)(−q) = spec(1,1)(q)) turns it into
    //   spec(1,4)(q) = −spec(1,1)(q) − 2N at EVERY q
    // (the (1,4) sibling of the §7 diamond fold spec(3,3) = −spec(2,3) − 2N). Then λ_A ∈ spec(1,4)(q*) ⟺
    // μ ∈ spec(1,1)(q*), excluded by B + Klein, and μ ∈ spec(1,4)(q*) ⟺ λ_A ∈ spec(1,1)(q*), excluded by A;
    // (4,1) follows by Klein. This fact pins the two symmetry ingredients (q-evenness, the composed
    // holomorphic fold) at the three locus classes and re-runs all four certificates: it alone IS the
    // closure of remainder R4's interior-core exclusion at N=5 (PROOF_CODIM1_BY_ADDITIVITY §7 scoping).
    [Fact]
    [Trait("Category", "R4INTERIOR")]
    public void InteriorFourExclusion_Propagates_FromTheTwoCertificates_KleinAndTheHolomorphicFold()
    {
        const int N = 5, twoN = 10;

        // (1) the two symmetry ingredients, at a real locus, the R-even deep locus, the R-odd escaper locus
        foreach (var q in new[] { new Complex(0.620878, 0), new Complex(1.8141, 0.3666), new Complex(1.7701, 1.2189) })
        {
            var s11 = Spec(N, 1, 1, q);
            Assert.True(MultisetDist(Spec(N, 1, 1, -q), s11) < 1e-9, $"q-evenness of spec(1,1) at q={q}");
            Assert.True(MultisetDist(Spec(N, 1, 4, q), s11.Select(z => -z - twoN).ToArray()) < 1e-9,
                $"holomorphic fold spec(1,4)(q) = −spec(1,1)(q) − 2N at q={q}");
        }
        _out.WriteLine("q-evenness + the composed holomorphic (1,4)-fold hold at all three locus classes");

        // (2) the four certificates: cert A (identity, both parities) + cert B (corner-fold, both parities)
        foreach (bool rOdd in new[] { false, true })
        {
            var a = FoldResultantCertificate.CertifyBlockExclusion(
                N, rOdd, tWKet: 1, tWBra: 1, GaussianInteger.One, GaussianInteger.Zero);
            var b = FoldResultantCertificate.CertifyComplete(N, rOdd);
            _out.WriteLine($"{(rOdd ? "R-odd " : "R-even")}: cert A (1,1)×λ_A complete={a.Complete} " +
                           $"(trueDegR={a.ResultantDegree}), cert B corner-fold complete={b.Complete} " +
                           $"(trueDegR={b.ResultantDegree})");
            Assert.True(a.Complete, $"cert A must complete ({(rOdd ? "R-odd" : "R-even")})");
            Assert.True(b.Complete, $"cert B must complete ({(rOdd ? "R-odd" : "R-even")})");
        }
    }

    private static Complex[] Spec(int n, int wk, int wb, Complex q) =>
        Matrix<Complex>.Build.DenseOfArray(WeightCoherenceBlock.Build(n, wk, wb, q)).Evd().EigenValues.ToArray();

    private static double MultisetDist(Complex[] a, Complex[] b)
    {
        if (a.Length != b.Length) return double.NaN;
        var used = new bool[b.Length];
        double worst = 0;
        foreach (var z in a)
        {
            int best = -1; double bd = double.PositiveInfinity;
            for (int j = 0; j < b.Length; j++)
                if (!used[j]) { double dd = (z - b[j]).Magnitude; if (dd < bd) { bd = dd; best = j; } }
            used[best] = true; worst = Math.Max(worst, bd);
        }
        return worst;
    }
}
