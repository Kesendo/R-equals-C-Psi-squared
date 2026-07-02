using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The EXACT holomorphic-fold identity (premise 1 of the complex-q resultant route, remainder R1):
/// the diagonal core (p_c,p_c) = (3,3) carries the HOLOMORPHIC fold μ = −λ − 2N of every (1,2) eigenvalue λ,
/// exactly and at every q, NOT the antiunitary −conj(λ) − 2N. The moved-window gate (ComplexQRateWindowTests)
/// only read this to locus precision (≤ 10⁻²) at a cell-0.05 scan; the resultant premise F_corner(−λ − 2N)
/// needs it as an exact identity.
///
/// <para>Finding (the two gates below): the core (3,3) carries −λ − 2N for every (1,2) RESIDUAL (F_18) root,
/// exact to ~10⁻¹³ in exact-arithmetic-at-integer-q (Gaussian-integer blocks), NOT the antiunitary −conj(λ) − 2N.
/// The fold is RESIDUAL-SPECIFIC: the full-spectrum identity spec(3,3) = −spec(band) − 2N is FALSE (second test,
/// ≈ 44), so only the residual roots fold, consistent with the diamond being a residual-braid phenomenon. It thus
/// carries the same status as remainder R3's gap byte-identity: OBSERVED EXACTLY, not yet derived from the
/// intertwiner. This upgrades the moved-window's 10⁻² locus-precision reading to the EXACT holomorphic identity
/// the resultant's F_corner(−λ − 2N) composition needs; the remaining analytic step is a from-first-principles
/// derivation via the diamond fold maps (the naive "two conjugations cancel to a full-spectrum identity" chain is
/// refuted by the second test, so the derivation must track the residual sub-structure specifically).</para>
///
/// <para>Run: <c>dotnet test "compute/RCPsiSquared.Diagnostics.Tests" --filter "Category=FOLDIDENTITY"
/// --logger "console;verbosity=detailed"</c>.</para></summary>
public sealed class HolomorphicFoldIdentityTests
{
    private readonly ITestOutputHelper _out;
    public HolomorphicFoldIdentityTests(ITestOutputHelper o) => _out = o;

    // Premise 1, the resultant-relevant form: the core (3,3) carries −λ − 2N for every (1,2) RESIDUAL root
    // (the F_18 roots the fold-resultant composes against), exact at integer q.
    [Fact]
    [Trait("Category", "FOLDIDENTITY")]
    public void N5_Core33_CarriesHolomorphicFold_OfEveryResidualRoot_ExactAtIntegerQ()
    {
        const int N = 5, twoN = 2 * N;
        foreach (int q0 in new[] { 2, 3, 5, 7 })
        {
            var q = new Complex(q0, 0);
            var residual = PathKMonodromyScout.ResidualRootsExact(N - 1, q);         // the F_18 (1,2) residual roots
            var core33 = Matrix<Complex>.Build.DenseOfArray(WeightCoherenceBlock.Build(N, 3, 3, q)).Evd().EigenValues;

            double worst = 0;
            foreach (var lam in residual)
            {
                var mu = -lam - twoN;                                               // the HOLOMORPHIC fold
                double d = core33.Enumerate().Select(z => (z - mu).Magnitude).Min();
                worst = Math.Max(worst, d);
            }
            _out.WriteLine($"q0={q0}: {residual.Length} residual roots; worst min|spec(3,3) − (−λ−2N)| = {worst:E2}");
            Assert.True(worst < 1e-8,
                $"core (3,3) must carry −λ−2N for every (1,2) residual root at q0={q0}; worst {worst:E2}");
        }
    }

    // The fold is RESIDUAL-SPECIFIC, not a full-spectrum identity: spec(3,3) ≠ −spec(band) − 2N for either band
    // orientation (both dim 100), so it is the F_18 residual roots that fold into the core, not the whole (1,2)
    // spectrum. Report-only: this records that the exact fold identity (the passing test above) belongs to the
    // residual/diamond family (the same "observed, not yet derived from the intertwiner" status as remainder R3's
    // gap byte-identity), NOT to a full block-spectral similarity.
    [Fact]
    [Trait("Category", "FOLDIDENTITY")]
    public void N5_Fold_IsResidualSpecific_NotFullSpectrum()
    {
        const int N = 5, twoN = 2 * N;
        var q = new Complex(3, 0);
        var core33 = Matrix<Complex>.Build.DenseOfArray(WeightCoherenceBlock.Build(N, 3, 3, q)).Evd().EigenValues
            .Enumerate().OrderBy(z => z.Real).ThenBy(z => z.Imaginary).ToArray();
        double best = double.PositiveInfinity;
        foreach (var (wk, wb) in new[] { (2, 3), (3, 2) })
        {
            var fold = Matrix<Complex>.Build.DenseOfArray(WeightCoherenceBlock.Build(N, wk, wb, q)).Evd().EigenValues
                .Enumerate().Select(z => -z - twoN).OrderBy(z => z.Real).ThenBy(z => z.Imaginary).ToArray();
            double worst = core33.Length == fold.Length
                ? core33.Zip(fold, (a, b) => (a - b).Magnitude).Max() : double.NaN;
            best = Math.Min(best, worst);
            _out.WriteLine($"full-spectrum band (wKet={wk},wBra={wb}) dim {fold.Length}: ‖sort spec(3,3) − sort(−spec−2N)‖_∞ = {worst:E2}");
        }
        _out.WriteLine($"⟹ the fold is RESIDUAL-SPECIFIC: no full-spectrum band identity (best {best:E2} ≫ 0), " +
                       "only the F_18 residual roots fold into the core (the passing residual test above), R3-family status.");
        Assert.True(best > 1.0, $"the fold must be residual-specific, not a full-spectrum identity; got best {best:E2}");
    }
}
