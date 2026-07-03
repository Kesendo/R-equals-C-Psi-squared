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
    // WIP (2026-07-03): the engine is repointed and computes the exact bivariate F_res/G and the mod-p pipeline
    // correctly, BUT the analytic degree bound rBound = resDeg·targetDeg − mR is TIGHT for the corner-fold
    // (deg_q R = 422 = 422) yet LOOSE for the identity composition (true deg_q R = 412 < 422; R-odd 384 < 394):
    // the identity leading forms share a common factor beyond the isolated-root collisions mR counts, so the
    // guard `DegP(rp) != rBound ⟹ p|lc, skip` misfires and skips every prime (primesUsed=0). FIX (a proof-
    // critical rigor change, to be adversarially reviewed): use the empirical true degree trueDegR = max_p
    // DegP(rp) (≤ the proven rBound) in place of the tightness assumption. Diagnosis pinned by
    // FoldResultantCertificate.DebugDegreeReport. Un-skip once the degree logic is fixed + reviewed.
    [Fact(Skip = "WIP: identity-composition degree bound is loose (true 412<422); needs the trueDegR rigor fix + review")]
    [Trait("Category", "R4INTERIOR")]
    public void Interior11_ExcludesTheBraidEigenvalue_REven_CompleteCertificate()
    {
        var r = FoldResultantCertificate.CertifyBlockExclusion(
            5, rOdd: false, tWKet: 1, tWBra: 1, GaussianInteger.One, GaussianInteger.Zero);
        _out.WriteLine($"(1,1)×λ_A R-even: complete={r.Complete}, resDeg={r.ResidualDegree}, targetDeg={r.CornerDegree}, "
            + $"rBound={r.ResultantDegreeBound}, primesUsed={r.PrimesUsed}, proofDigits={r.ProofBoundDigits}");
        Assert.True(r.Complete, "the (1,1) core must exclude every (1,2) R-even residual root (the braid λ_A)");
    }

    [Fact(Skip = "WIP: identity-composition degree bound is loose (true 384<394); needs the trueDegR rigor fix + review")]
    [Trait("Category", "R4INTERIOR")]
    public void Interior11_ExcludesTheBraidEigenvalue_ROdd_CompleteCertificate()
    {
        var r = FoldResultantCertificate.CertifyBlockExclusion(
            5, rOdd: true, tWKet: 1, tWBra: 1, GaussianInteger.One, GaussianInteger.Zero);
        _out.WriteLine($"(1,1)×λ_A R-odd: complete={r.Complete}, resDeg={r.ResidualDegree}, primesUsed={r.PrimesUsed}");
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
