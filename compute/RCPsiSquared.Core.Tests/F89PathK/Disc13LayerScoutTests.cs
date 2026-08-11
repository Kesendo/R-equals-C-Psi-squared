using System;
using RCPsiSquared.Core.F89PathK;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>The (1,3)@N=6 disc-layer scout (sideways_spin_ladder arc, the census's exact follow-up):
/// the weight-general sector pencil + derived-rate AT clearing
/// (<see cref="FoldResultantCertificate.WeightSectorPencil"/> /
/// <see cref="F89AtFactorReconstruction.ClearedAtSectorsFromPencil"/>) validated against the
/// committed (1,2) scout pins, then run on the census block. The separating object is the
/// multiplicity-2 layer (doubled defective = double disc root, doubled diabolic = quadruple, away
/// from Λ = −2N), with the Λ = −2N stratum read separately as
/// deg gcd(F_res(−2N, q), ∂Λ F_res(−2N, q)) mod p, the full-evenness read, and the composition
/// identity chaining the two reads live. Every pin measured 2026-08-11 and cross-checked at two
/// split primes; independent per-point reproductions by two review rounds.</summary>
public class Disc13LayerScoutTests
{
    private readonly ITestOutputHelper _out;
    public Disc13LayerScoutTests(ITestOutputHelper output) => _out = output;

    /// <summary>The generalized path must reproduce the committed (1,2)@N=4 scout pins
    /// (EvenNDiscLayerScoutTests.N4_Control_KnownOcticFactorization: 52, 24, [20, 4]) from the
    /// WeightCoherenceBlock route, both parities, and read the N=4 disc REAL mod p (the exact
    /// form of the committed N4_ReImControl). Fast (~1 s).</summary>
    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    [Trait("Category", "EVEN_DISC13_SCOUT")]
    public void N4_Control_GeneralPathMatchesCommitted12Pins(bool rOdd)
    {
        var (p, degD, vD, layers, atDeg, resDeg, discReal, m2Gcd) =
            FoldResultantCertificate.WeightDiscLayersAtNthPrime(4, 1, 2, rOdd, 0);
        _out.WriteLine($"N4 (1,2) rOdd={rOdd}: p={p} degD={degD} vD={vD} layers=[{string.Join(", ", layers)}] at={atDeg} res={resDeg} real={discReal} m2Gcd={m2Gcd}");
        Assert.Equal(52, degD);
        Assert.Equal(24, vD);
        Assert.Equal(new[] { 20, 4 }, layers);
        Assert.True(discReal, "the N=4 (1,2) disc is real (the self-fold), mod p included");
        Assert.Equal(-1, m2Gcd);   // f squarefree, no evenness: the u-analysis is void here
    }

    /// <summary>The generalized path must reproduce the committed (1,2)@N=6 scout pins
    /// (EvenNDiscLayerScoutTests via the certificate: 926, 536, [124, 133]), both parities, and
    /// read the (1,2)@N=6 disc genuinely COMPLEX mod p (the two embeddings disagree).</summary>
    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    [Trait("Category", "SLOW_EVEN_DISC13")]
    public void N6_Control_GeneralPathMatchesCommitted12Pins(bool rOdd)
    {
        var (p, degD, vD, layers, atDeg, resDeg, discReal, _) =
            FoldResultantCertificate.WeightDiscLayersAtNthPrime(6, 1, 2, rOdd, 0);
        _out.WriteLine($"N6 (1,2) rOdd={rOdd}: p={p} degD={degD} vD={vD} layers=[{string.Join(", ", layers)}] at={atDeg} res={resDeg} real={discReal}");
        Assert.Equal(926, degD);
        Assert.Equal(536, vD);
        Assert.Equal(new[] { 124, 133 }, layers);
        Assert.False(discReal, "the N=6 (1,2) sector disc is genuinely complex");
    }

    /// <summary>The (1,3)@N=6 layer profile, both parities, cross-checked at the first two split
    /// primes (1073741833, 1073741857; measured identical 2026-08-11, ~2.3 s per read). The pins:
    /// NO odd-multiplicity layer exists (layers[0] = layers[2] = 0) — the census's Klein even-order
    /// theorem exact mod p, the reason a sign-change detector is structurally blind here; the
    /// multiplicity-2 layer (108 / 132 distinct roots) is the stratum carrying the census loci
    /// (8 + 5 real ones found there, a floor closed at 9 + 6 = fifteen by the exact Sturm count,
    /// Disc13SturmTests; most mult-2 roots sit at complex q), the
    /// multiplicity-4 layer (120 / 168) the doubled diabolic crossings; and the sector disc is
    /// REAL mod p (both embeddings agree), the exact form of the census's ≤ 5·10⁻¹² float read.</summary>
    [Theory]
    [InlineData(false, 0)]
    [InlineData(true, 0)]
    [InlineData(false, 1)]
    [InlineData(true, 1)]
    [Trait("Category", "SLOW_EVEN_DISC13")]
    public void N6_Scout_Disc13Layers(bool rOdd, int nth)
    {
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var (p, degD, vD, layers, atDeg, resDeg, discReal, m2Gcd) =
            FoldResultantCertificate.WeightDiscLayersAtNthPrime(6, 1, 3, rOdd, nth);
        _out.WriteLine($"N6 (1,3) rOdd={rOdd} nth={nth}: p={p} degD={degD} vD={vD} layers=[{string.Join(", ", layers)}] at={atDeg} res={resDeg} real={discReal} m2Gcd={m2Gcd}  [{sw.Elapsed}]");
        Assert.Equal(rOdd ? 2124 : 1572, degD);
        Assert.Equal(rOdd ? 1188 : 876, vD);
        Assert.Equal(rOdd ? new[] { 0, 132, 0, 168 } : new[] { 0, 108, 0, 120 }, layers);
        Assert.Equal(rOdd ? 12 : 18, atDeg);
        Assert.Equal(rOdd ? 48 : 42, resDeg);
        Assert.True(discReal, "the (1,3) sector disc is real: gauge + fold both act within the sector");
        // the aggregation identity: deg − v = Σ (k+1)·layers[k]
        int agg = 0;
        for (int k = 0; k < layers.Length; k++) agg += (k + 1) * layers[k];
        Assert.Equal(degD - vD, agg);
        // u | L₂ (deg gcd(u, mult-2 layer) = deg u): gcd(u, disc_M G) = 1, the hypothesis that
        // (i) excludes the 1+1 case (two coincident defective EPs on μ = 0 read as an order-2
        // f-root too; they would push a u-root into the mult-4 layer) and (ii) fixes the mult-2
        // split exactly (18/24 of its roots are the μ = 0 diabolics). A review round found both
        // conclusions silently resting on this.
        Assert.Equal(rOdd ? 24 : 18, m2Gcd);
    }

    /// <summary>The Λ = −12 stratum of (1,3)@N=6, both parities, and the (1,2)@N=4 control at its
    /// own line Λ = −8. Pins measured 2026-08-11: at (1,3) the Λ-derivative of F_res at −12 is
    /// IDENTICALLY ZERO mod p (GDeg = −1) — F_res is a polynomial in (Λ+12)², the within-sector
    /// fold at the ALGEBRA level (gauge-real coefficients + fold evenness), so gcd(f, ∂Λ) = f
    /// whole; and f = F_res(−12, q) is a PERFECT SQUARE (f-layers [0, 18] / [0, 24]): μ and −μ
    /// strands touch the μ = 0 line only together. The N=4 (1,2) control discriminates on every
    /// read: its disc is real (the antiunitary fold) but its sector is NOT a real family (the
    /// gauge anticommutes with R at (p+q)(N−1) = 9 odd), so GDeg = 7 (no evenness), gcd = 0 and
    /// f squarefree ([8]) — single strands cross the line alone, the s-partner living in the
    /// other sector.</summary>
    [Theory]
    [InlineData(6, 1, 3, false, -12, 0, 36, -1, 36, new[] { 0, 18 }, true)]
    [InlineData(6, 1, 3, true, -12, 0, 48, -1, 48, new[] { 0, 24 }, true)]
    [InlineData(6, 1, 3, false, -12, 1, 36, -1, 36, new[] { 0, 18 }, true)]
    [InlineData(6, 1, 3, true, -12, 1, 48, -1, 48, new[] { 0, 24 }, true)]
    [InlineData(4, 1, 2, false, -8, 0, 8, 7, 0, new[] { 8 }, false)]
    [InlineData(4, 1, 2, true, -8, 0, 8, 7, 0, new[] { 8 }, false)]
    [Trait("Category", "SLOW_EVEN_DISC13")]
    public void LambdaStratum_FoldEvenness_And_SquareStratum(
        int n, int wKet, int wBra, bool rOdd, int lambdaCleared, int nth,
        int expFDeg, int expGDeg, int expGcdDeg, int[] expFLayers, bool expEven)
    {
        var (p, fDeg, gDeg, gcdDeg, gcdVal, fLayers, evenShifted) =
            FoldResultantCertificate.WeightLambdaStratumAtNthPrime(n, wKet, wBra, rOdd, lambdaCleared, nth);
        _out.WriteLine($"N{n} ({wKet},{wBra}) rOdd={rOdd} Lam0={lambdaCleared} nth={nth}: p={p} fDeg={fDeg} gDeg={gDeg} gcdDeg={gcdDeg} gcdQVal={gcdVal} fLayers=[{string.Join(", ", fLayers)}] even={evenShifted}");
        Assert.Equal(expFDeg, fDeg);
        Assert.Equal(expGDeg, gDeg);
        Assert.Equal(expGcdDeg, gcdDeg);
        Assert.Equal(expFLayers, fLayers);
        Assert.Equal(0, gcdDeg < 0 ? 0 : gcdVal);
        // EVERY odd Taylor coefficient of F_res at Λ₀ vanishing EXACTLY over ℤ[i] (not just a₁ =
        // the gDeg read, and not merely mod p; two review rounds tightened this in turn) ⟺
        // F_res ∈ ℤ[i][(Λ−Λ₀)², q]: the fold as algebra, a theorem of the exact residual, present
        // at (1,3)@N=6 (gauge-real + within-sector fold), absent at the N=4 control. The mod-p
        // reads beside it (gcd of degree 0 means a nonzero CONSTANT gcd; f-layers) get their
        // second prime from the nth = 1 rows.
        Assert.Equal(expEven, evenShifted);
    }

    /// <summary>The composition identity disc_Λ(F) = ±4^m·f·disc_M(G)² for F = G((Λ+12)²), chaining
    /// the LIVE layer profile and the LIVE stratum read against each other (a first draft asserted
    /// it over hardcoded literals, which a review round noted cannot catch a drifted
    /// re-measurement): deg disc G = (deg D − deg f)/2, v(disc G) = v(D)/2 (v(f) = 0, so v(D) is
    /// FORCED even), and the stripped-layer identity 1·L₂ + 2·L₄ = deg u + (deg disc G − v disc G)
    /// with f = c·u². R-even: 768, 438, 108 + 2·120 = 348 = 18 + 330; R-odd: 1038, 594,
    /// 132 + 2·168 = 468 = 24 + 444.</summary>
    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    [Trait("Category", "SLOW_EVEN_DISC13")]
    public void CompositionIdentity_ChainsLayersAndStratum(bool rOdd)
    {
        var (_, degD, vD, layers, _, _, _, _) =
            FoldResultantCertificate.WeightDiscLayersAtNthPrime(6, 1, 3, rOdd, 0);
        var (_, fDeg, _, _, _, fLayers, evenShifted) =
            FoldResultantCertificate.WeightLambdaStratumAtNthPrime(6, 1, 3, rOdd, -12, 0);
        Assert.True(evenShifted);
        Assert.Equal(4, layers.Length);
        Assert.Equal(2, fLayers.Length);
        Assert.Equal(0, fLayers[0]);
        int uDeg = fLayers[1];
        Assert.Equal(fDeg, 2 * uDeg);
        Assert.Equal(0, vD % 2);
        int degDiscG = (degD - fDeg) / 2;
        int vDiscG = vD / 2;
        Assert.Equal(uDeg + (degDiscG - vDiscG), layers[1] + 2 * layers[3]);
    }
}
