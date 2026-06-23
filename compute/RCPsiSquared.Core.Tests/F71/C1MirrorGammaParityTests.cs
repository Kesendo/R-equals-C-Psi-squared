using RCPsiSquared.Core.F71;
using RCPsiSquared.Core.Knowledge;
using Xunit;

namespace RCPsiSquared.Core.Tests.F71;

/// <summary>Tests for F101 <see cref="C1MirrorGammaParity"/>: the F71 c₁ bond-mirror
/// deviation D(b) = c₁(b) − c₁(N−2−b) is exactly odd in the F71-anti-palindromic
/// γ-component (graceful breakdown of F71 under non-uniform per-site dephasing).
///
/// <para>These tests cover the typed claim's Tier and its γ_sym / γ_anti decomposition
/// helpers (site-mirror l ↔ N−1−l). The numerical verification lives in the witness
/// script <c>simulations/f71_nonuniform_gamma_verification.py</c> (oddness + palindromic
/// survival, N=3,4,5), recorded in
/// <c>docs/proofs/PROOF_F101_C1_MIRROR_GAMMA_PARITY.md</c>.</para></summary>
public sealed class C1MirrorGammaParityTests
{
    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, new C1MirrorGammaParity().Tier);
    }

    [Theory]
    [InlineData(new double[] { 0.05, 0.05, 0.05 })]          // uniform
    [InlineData(new double[] { 0.07, 0.05, 0.07 })]          // non-uniform palindromic
    [InlineData(new double[] { 0.03, 0.06, 0.06, 0.03 })]    // palindromic, N=4
    public void IsPalindromic_TrueForPalindromicGamma(double[] siteGamma)
    {
        Assert.True(C1MirrorGammaParity.IsPalindromic(siteGamma));
        Assert.True(C1MirrorGammaParity.PalindromicDeviation(siteGamma) < 1e-12);
    }

    [Theory]
    [InlineData(new double[] { 0.05, 0.05, 0.09 })]          // γ_0 != γ_2
    [InlineData(new double[] { 0.03, 0.05, 0.07 })]          // linear ramp
    public void IsPalindromic_FalseForNonPalindromicGamma(double[] siteGamma)
    {
        Assert.False(C1MirrorGammaParity.IsPalindromic(siteGamma));
        Assert.True(C1MirrorGammaParity.PalindromicDeviation(siteGamma) > 1e-3);
    }

    [Fact]
    public void Decomposition_GammaSymPlusGammaAnti_RecoversGamma()
    {
        var siteGamma = new double[] { 0.03, 0.05, 0.07 };
        var sym = C1MirrorGammaParity.PalindromicComponent(siteGamma);
        var anti = C1MirrorGammaParity.AntiPalindromicComponent(siteGamma);
        for (int l = 0; l < siteGamma.Length; l++)
            Assert.Equal(siteGamma[l], sym[l] + anti[l], precision: 12);
    }

    [Fact]
    public void PalindromicComponent_IsItselfPalindromic()
    {
        var sym = C1MirrorGammaParity.PalindromicComponent(new double[] { 0.03, 0.05, 0.07 });
        for (int l = 0; l < sym.Length; l++)
            Assert.Equal(sym[l], sym[sym.Length - 1 - l], precision: 12);
    }

    [Fact]
    public void AntiPalindromicComponent_IsItselfAntiPalindromic()
    {
        var anti = C1MirrorGammaParity.AntiPalindromicComponent(new double[] { 0.03, 0.05, 0.07 });
        for (int l = 0; l < anti.Length; l++)
            Assert.Equal(anti[l], -anti[anti.Length - 1 - l], precision: 12);
    }

    [Fact]
    public void AntiPalindromicComponent_OfLinearRamp_MatchesClosedForm()
    {
        // γ = [0.03, 0.05, 0.07], F71-reverse = [0.07, 0.05, 0.03]
        // γ_anti = (γ − F71(γ)) / 2 = [-0.02, 0.0, 0.02]
        var anti = C1MirrorGammaParity.AntiPalindromicComponent(new double[] { 0.03, 0.05, 0.07 });
        Assert.Equal(-0.02, anti[0], precision: 12);
        Assert.Equal(0.0, anti[1], precision: 12);
        Assert.Equal(0.02, anti[2], precision: 12);
    }
}
