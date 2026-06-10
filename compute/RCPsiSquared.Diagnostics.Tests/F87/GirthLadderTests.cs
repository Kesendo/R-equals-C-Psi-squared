using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F87;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>Pins the girth-ladder predictor (<see cref="GirthLadder.Forecast"/>) to the exact
/// witnesses of the Python anchors (<c>simulations/f87_windowed_monomial_converse.py</c> canonical
/// pairs + the two k=4 pairs of <c>simulations/f87_girth_dichotomy.py</c> Block 3): the deg-1
/// outright branch with its closed-form coefficient (IXXZ+XIXZ 573440, IIZ+IZI 9216 == the P31
/// identity 6·4^N·Σc_l²), the higher-rung branch on all three ℓ-kinds (K3, flux, multi-Z lift,
/// the γ⁵ witness IIXY+ZXZY), and the bipartite control. Everything here is dense-H-only
/// (2^N ≤ 32), no Liouvillian, so the whole class runs in milliseconds.</summary>
public class GirthLadderTests
{
    private static PauliTerm T(string label) => new(PauliLabel.Parse(label), Complex.One);
    private static List<PauliTerm> Pair(params string[] labels) => labels.Select(T).ToList();

    private const double RelTol = 1e-6;

    private static void AssertRelEqual(double expected, double actual) =>
        Assert.True(Math.Abs(actual - expected) <= RelTol * Math.Abs(expected),
            $"expected {expected} within rel {RelTol}, got {actual}");

    [Fact]
    public void Forecast_IXXZplusXIXZ_N5_IsDeg1Outright_MStar7_Coefficient573440()
    {
        // The first k=4 deg-1 pure-cycle representative (girth dichotomy Block 3): ℓ = 3 cycle,
        // t₃ = [0, 0, 0, 64, 0], so m* = 7 exact and the coefficient is the closed form
        // P_{7,1} = 7·C(6,3)·Σt₃² = 7·20·4096 = 573440 (matches the exact CRT p₇).
        var f = GirthLadder.Forecast(Pair("IXXZ", "XIXZ"), n: 5);

        Assert.Equal(3, f.Ell);
        Assert.Equal("3-cycle", f.EllKind);
        Assert.Equal(GirthLadderBranch.Deg1Outright, f.Branch);
        Assert.Equal(7, f.MStar);
        Assert.False(f.MStarIsLowerBound);
        AssertRelEqual(573440.0, f.Coefficient);

        // consistency with 7·C(6,3)·Σt₃² from the reported per-site moments
        double sumSq = f.GirthMoments.Sum(t => t * t);
        AssertRelEqual(4096.0, sumSq);
        AssertRelEqual(7 * 20 * sumSq, f.Coefficient);
    }

    [Fact]
    public void Forecast_K3_XXZplusXZX_N4_IsHigherRung_Bound9()
    {
        // K3 triangle: ℓ = 3 cycle with t₃ ≡ 0 (deg-3 branch); the exact engine fires at
        // m* = 9 = 2·3+3, which the lower bound reaches exactly here.
        var f = GirthLadder.Forecast(Pair("XXZ", "XZX"), n: 4);

        Assert.Equal(3, f.Ell);
        Assert.Equal("3-cycle", f.EllKind);
        Assert.Equal(GirthLadderBranch.HigherRung, f.Branch);
        Assert.Equal(9, f.MStar);
        Assert.True(f.MStarIsLowerBound);
        Assert.Equal(0.0, f.Coefficient);
        Assert.All(f.GirthMoments, t => Assert.True(Math.Abs(t) < 1e-9));
    }

    [Fact]
    public void Forecast_Flux_IXYplusXIY_N4_ComplexH_IsHigherRung_Bound9()
    {
        // The flux pair (odd #Y) has COMPLEX H with Gaussian-integer entries; the signed 3-walk
        // cancels (t₃ ≡ 0) although the unsigned 3-cycle exists, so ℓ = 3 with the higher rung
        // (exact engine: m* = 9, γ³).
        var f = GirthLadder.Forecast(Pair("IXY", "XIY"), n: 4);

        Assert.Equal(3, f.Ell);
        Assert.Equal("3-cycle", f.EllKind);
        Assert.Equal(GirthLadderBranch.HigherRung, f.Branch);
        Assert.Equal(9, f.MStar);
        Assert.True(f.MStarIsLowerBound);
        Assert.All(f.GirthMoments, t => Assert.True(Math.Abs(t) < 1e-9));
    }

    [Fact]
    public void Forecast_MultiZ_XXZplusZZZ_N4_IsDiagonalLift_HigherRung_Bound5()
    {
        // ZZZ lifts the diagonal, so ℓ = 1 without any hopping cycle; but there is no
        // single-site-Z component, so t₁ ≡ 0 and the deg-1 rung is dead (exact engine:
        // m* = 5 = 2·1+3, γ³, p₅ = 61440·γ³).
        var f = GirthLadder.Forecast(Pair("XXZ", "ZZZ"), n: 4);

        Assert.Equal(1, f.Ell);
        Assert.Equal("diagonal-lift", f.EllKind);
        Assert.Equal(GirthLadderBranch.HigherRung, f.Branch);
        Assert.Equal(5, f.MStar);
        Assert.True(f.MStarIsLowerBound);
        Assert.Equal(0.0, f.Coefficient);
        Assert.All(f.GirthMoments, t => Assert.True(Math.Abs(t) < 1e-9));
    }

    [Fact]
    public void Forecast_IIXYplusZXZY_N5_IsHigherRung_Bound9_ActualRungIsEleven()
    {
        // The γ⁵ rung witness (girth dichotomy Block 3): ℓ = 3, t₃ ≡ 0, AND the γ³ class also
        // dies, so the exact engine fires only at m* = 11 = 2·3+5 (p₁₁ = 86507520·γ⁵). The
        // forecast's 2ℓ+3 = 9 is therefore a strict LOWER bound here, which is exactly what
        // MStarIsLowerBound carries.
        var f = GirthLadder.Forecast(Pair("IIXY", "ZXZY"), n: 5);

        Assert.Equal(3, f.Ell);
        Assert.Equal("3-cycle", f.EllKind);
        Assert.Equal(GirthLadderBranch.HigherRung, f.Branch);
        Assert.Equal(9, f.MStar);
        Assert.True(f.MStarIsLowerBound);
        Assert.Equal(0.0, f.Coefficient);
    }

    [Fact]
    public void Forecast_DiagLift_IIZplusIZI_N4_IsDeg1Outright_MStar3_Coefficient9216_MatchesP31ClosedForm()
    {
        // The single-site-Z lift (the anchor's DIAG_LIFT): ℓ = 1 diagonal lift with
        // t₁ = [0, 16, 32, 16] = 2^N·c_l (c = [0, 1, 2, 1]), m* = 3 exact, and the two closed
        // forms agree: P_{3,1} = 3·C(2,1)·Σt₁² = 9216 = 6·4^N·Σc_l².
        const int n = 4;
        var f = GirthLadder.Forecast(Pair("IIZ", "IZI"), n);

        Assert.Equal(1, f.Ell);
        Assert.Equal("diagonal-lift", f.EllKind);
        Assert.Equal(GirthLadderBranch.Deg1Outright, f.Branch);
        Assert.Equal(3, f.MStar);
        Assert.False(f.MStarIsLowerBound);
        AssertRelEqual(9216.0, f.Coefficient);

        // per-site girth moments are pinned: t₁ = [0, 16, 32, 16]
        Assert.Equal(new[] { 0.0, 16.0, 32.0, 16.0 }, f.GirthMoments);

        // cross-check the P31 identity: c_l = t₁^(l) / 2^N, P_{3,1} = 6·4^N·Σc_l²
        double d = 1 << n;
        double sumCSq = f.GirthMoments.Sum(t => (t / d) * (t / d));
        AssertRelEqual(6 * Math.Pow(4, n) * sumCSq, f.Coefficient);
    }

    [Fact]
    public void Forecast_Bipartite_XXZplusZXX_N4_IsBipartite()
    {
        // The soft control: no diagonal lift and the hopping graph is bipartite (no odd cycle),
        // so every odd power sum vanishes exactly and the pair is soft at every γ.
        var f = GirthLadder.Forecast(Pair("XXZ", "ZXX"), n: 4);

        Assert.Equal(GirthLadderBranch.Bipartite, f.Branch);
        Assert.Equal(0, f.Ell);
        Assert.Equal("bipartite", f.EllKind);
        Assert.Equal(0, f.MStar);
        Assert.False(f.MStarIsLowerBound);
        Assert.Equal(0.0, f.Coefficient);
        Assert.Empty(f.GirthMoments);
    }
}
