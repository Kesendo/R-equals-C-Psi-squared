using RCPsiSquared.Core.F71;
using RCPsiSquared.Core.Knowledge;
using Xunit;

namespace RCPsiSquared.Core.Tests.F71;

/// <summary>Tests for F100 <see cref="C1QPeakMirrorJParity"/>: the F71 c₁/Q_peak
/// bond-mirror deviation D(b) = c₁(b) − c₁(N−2−b) is exactly odd in the
/// F71-anti-palindromic J-component (graceful breakdown of F71 under non-uniform J).
///
/// <para>These tests cover the typed claim's Tier and its J_sym / J_anti decomposition
/// helpers. The numerical verification lives in the witness scripts: c₁ in
/// <c>simulations/f71_nonuniform_j_verification.py</c> (oddness + palindromic survival,
/// N=3,4,5), Q_peak in <c>simulations/f100_qpeak_nonuniform_j_verification.py</c>
/// (N=4,5,6 plus a c=3 spot check). Both are recorded in
/// <c>docs/proofs/PROOF_F100_C1_QPEAK_MIRROR_J_PARITY.md</c>.</para></summary>
public sealed class C1QPeakMirrorJParityTests
{
    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, new C1QPeakMirrorJParity().Tier);
    }

    [Theory]
    [InlineData(new double[] { 1.0, 1.0, 1.0 })]        // uniform
    [InlineData(new double[] { 1.2, 0.8, 1.2 })]        // non-uniform palindromic
    [InlineData(new double[] { 0.3, 0.5, 0.5, 0.3 })]   // palindromic, N=5
    public void IsPalindromic_TrueForPalindromicJ(double[] bondJ)
    {
        Assert.True(C1QPeakMirrorJParity.IsPalindromic(bondJ));
        Assert.True(C1QPeakMirrorJParity.PalindromicDeviation(bondJ) < 1e-12);
    }

    [Theory]
    [InlineData(new double[] { 1.0, 1.0, 1.4 })]        // J_0 != J_2
    [InlineData(new double[] { 0.7, 1.0, 1.3 })]        // linear ramp
    public void IsPalindromic_FalseForNonPalindromicJ(double[] bondJ)
    {
        Assert.False(C1QPeakMirrorJParity.IsPalindromic(bondJ));
        Assert.True(C1QPeakMirrorJParity.PalindromicDeviation(bondJ) > 1e-3);
    }

    [Fact]
    public void Decomposition_JSymPlusJAnti_RecoversJ()
    {
        var bondJ = new double[] { 0.7, 1.0, 1.3 };
        var sym = C1QPeakMirrorJParity.PalindromicComponent(bondJ);
        var anti = C1QPeakMirrorJParity.AntiPalindromicComponent(bondJ);
        for (int b = 0; b < bondJ.Length; b++)
            Assert.Equal(bondJ[b], sym[b] + anti[b], precision: 12);
    }

    [Fact]
    public void PalindromicComponent_IsItselfPalindromic()
    {
        var sym = C1QPeakMirrorJParity.PalindromicComponent(new double[] { 0.7, 1.0, 1.3 });
        for (int b = 0; b < sym.Length; b++)
            Assert.Equal(sym[b], sym[sym.Length - 1 - b], precision: 12);
    }

    [Fact]
    public void AntiPalindromicComponent_IsItselfAntiPalindromic()
    {
        var anti = C1QPeakMirrorJParity.AntiPalindromicComponent(new double[] { 0.7, 1.0, 1.3 });
        for (int b = 0; b < anti.Length; b++)
            Assert.Equal(anti[b], -anti[anti.Length - 1 - b], precision: 12);
    }

    [Fact]
    public void AntiPalindromicComponent_OfLinearRamp_MatchesClosedForm()
    {
        // J = [0.7, 1.0, 1.3], F71-reverse = [1.3, 1.0, 0.7]
        // J_anti = (J - F71(J)) / 2 = [-0.3, 0.0, 0.3]
        var anti = C1QPeakMirrorJParity.AntiPalindromicComponent(new double[] { 0.7, 1.0, 1.3 });
        Assert.Equal(-0.3, anti[0], precision: 12);
        Assert.Equal(0.0, anti[1], precision: 12);
        Assert.Equal(0.3, anti[2], precision: 12);
    }
}
