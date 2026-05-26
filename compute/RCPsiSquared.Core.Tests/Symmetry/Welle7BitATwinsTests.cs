using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Welle 7 (2026-05-26): structural tests for the 4 new BitA twin Claims +
/// 4 BitATwinStatus reclassifications.
///
/// <para>Track A net effect on the PolarityCubeMap inventory: BitA Claim count
/// 2 → 6 (added F38BitA, F39BitA, F63BitAReference, ZGlobalEigenstateMirrorBitA);
/// BitBSpecific status count 2 → 6 (added F82, F84, F91, F93); NeedsDerivation
/// count drops by 4.</para></summary>
public class Welle7BitATwinsTests
{
    // ------------------------------------------------------------------
    // A1. F38BitAInvolutionInheritance
    // ------------------------------------------------------------------

    [Fact]
    public void F38BitA_Z2Axis_IsBitA()
    {
        var claim = new F38BitAInvolutionInheritance();
        Assert.Equal(Z2Axis.BitA, claim.Z2Axis);
    }

    [Fact]
    public void F38BitA_Tier_IsTier1Derived()
    {
        var claim = new F38BitAInvolutionInheritance();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void F38BitA_BitATwinStatus_IsNotApplicableForThisAxis()
    {
        IZ2AxisClaim claim = new F38BitAInvolutionInheritance();
        Assert.Equal(BitATwinClassification.NotApplicableForThisAxis, claim.BitATwinStatus);
    }

    [Fact]
    public void F38BitA_BitATwin_IsNull()
    {
        var claim = new F38BitAInvolutionInheritance();
        Assert.Null(claim.BitATwin);
    }

    [Fact]
    public void F38BitA_Theorem_IsNonEmpty()
    {
        var claim = new F38BitAInvolutionInheritance();
        Assert.False(string.IsNullOrWhiteSpace(claim.Theorem));
    }

    [Theory]
    [InlineData(1, 4L, 2L)]
    [InlineData(2, 16L, 8L)]
    [InlineData(3, 64L, 32L)]
    [InlineData(5, 1024L, 512L)]
    public void F38BitA_DimensionsMatchF38(int N, long fullDim, long eigenspaceDim)
    {
        var claim = new F38BitAInvolutionInheritance();
        Assert.Equal(fullDim, claim.FullOperatorSpaceDimension(N));
        Assert.Equal(eigenspaceDim, claim.EigenspaceDimension(N));
        Assert.Equal(0.5, claim.HalfHalfBalance(N), precision: 14);
    }

    // ------------------------------------------------------------------
    // A2. F39DetPiBitAInheritance
    // ------------------------------------------------------------------

    [Fact]
    public void F39BitA_Z2Axis_IsBitA()
    {
        var claim = new F39DetPiBitAInheritance();
        Assert.Equal(Z2Axis.BitA, claim.Z2Axis);
    }

    [Fact]
    public void F39BitA_Tier_IsTier1Derived()
    {
        var claim = new F39DetPiBitAInheritance();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void F39BitA_BitATwinStatus_IsNotApplicableForThisAxis()
    {
        IZ2AxisClaim claim = new F39DetPiBitAInheritance();
        Assert.Equal(BitATwinClassification.NotApplicableForThisAxis, claim.BitATwinStatus);
    }

    [Fact]
    public void F39BitA_Theorem_IsNonEmpty()
    {
        var claim = new F39DetPiBitAInheritance();
        Assert.False(string.IsNullOrWhiteSpace(claim.Theorem));
    }

    [Theory]
    [InlineData(1, -1)]
    [InlineData(2, 1)]
    [InlineData(3, 1)]
    public void F39BitA_DetPiX_MatchesF39(int N, int expected)
    {
        var claim = new F39DetPiBitAInheritance();
        Assert.Equal(expected, claim.DetPiX(N));
    }

    // ------------------------------------------------------------------
    // A3. F63BitAReference
    // ------------------------------------------------------------------

    [Fact]
    public void F63BitARef_Z2Axis_IsBitA()
    {
        var claim = new F63BitAReference();
        Assert.Equal(Z2Axis.BitA, claim.Z2Axis);
    }

    [Fact]
    public void F63BitARef_Tier_IsTier1Derived()
    {
        var claim = new F63BitAReference();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void F63BitARef_BitATwinStatus_IsNotApplicableForThisAxis()
    {
        IZ2AxisClaim claim = new F63BitAReference();
        Assert.Equal(BitATwinClassification.NotApplicableForThisAxis, claim.BitATwinStatus);
    }

    [Fact]
    public void F63BitARef_Theorem_IsNonEmpty()
    {
        var claim = new F63BitAReference();
        Assert.False(string.IsNullOrWhiteSpace(claim.Theorem));
    }

    // ------------------------------------------------------------------
    // A4. ZGlobalEigenstateMirrorBitAInheritance
    // ------------------------------------------------------------------

    [Fact]
    public void ZGlobalMirrorBitA_Z2Axis_IsBitA()
    {
        var claim = new ZGlobalEigenstateMirrorBitAInheritance(new HalfAsStructuralFixedPointClaim());
        Assert.Equal(Z2Axis.BitA, claim.Z2Axis);
    }

    [Fact]
    public void ZGlobalMirrorBitA_Tier_IsTier1Derived()
    {
        var claim = new ZGlobalEigenstateMirrorBitAInheritance(new HalfAsStructuralFixedPointClaim());
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void ZGlobalMirrorBitA_BitATwinStatus_IsNotApplicableForThisAxis()
    {
        IZ2AxisClaim claim = new ZGlobalEigenstateMirrorBitAInheritance(new HalfAsStructuralFixedPointClaim());
        Assert.Equal(BitATwinClassification.NotApplicableForThisAxis, claim.BitATwinStatus);
    }

    [Fact]
    public void ZGlobalMirrorBitA_AlphaAtMirror_IsZero()
    {
        Assert.Equal(0.0, ZGlobalEigenstateMirrorBitAInheritance.AlphaAtMirror);
    }

    [Fact]
    public void ZGlobalMirrorBitA_GammaAtMirror_IsOne()
    {
        Assert.Equal(1.0, ZGlobalEigenstateMirrorBitAInheritance.GammaAtMirror);
    }

    [Theory]
    [InlineData(1.0, 0.0)]
    [InlineData(-1.0, 0.0)]
    public void ZGlobalMirrorBitA_AlphaFromGamma_MatchesUniversalShape(double gammaX, double expectedAlpha)
    {
        Assert.Equal(expectedAlpha, ZGlobalEigenstateMirrorBitAInheritance.AlphaFromGammaAtMirror(gammaX), precision: 14);
    }

    [Fact]
    public void ZGlobalMirrorBitA_NullHalf_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new ZGlobalEigenstateMirrorBitAInheritance(null!));
    }

    // ------------------------------------------------------------------
    // BitB siblings: BitATwin populated when wired with the new ctor parameter
    // ------------------------------------------------------------------

    private static F38Pi2InvolutionPi2Inheritance BuildF38(F38BitAInvolutionInheritance? bitA = null) =>
        new F38Pi2InvolutionPi2Inheritance(
            new Pi2DyadicLadderClaim(),
            new Pi2OperatorSpaceMirrorClaim(),
            new Pi2I4MemoryLoopClaim(),
            new HalfAsStructuralFixedPointClaim(),
            bitA);

    [Fact]
    public void F38_WithBitATwin_BitATwin_IsF38BitA()
    {
        var bitA = new F38BitAInvolutionInheritance();
        var f38 = BuildF38(bitA);
        Assert.Same(bitA, f38.BitATwin);
        Assert.Equal(BitATwinClassification.Filled, f38.BitATwinStatus);
    }

    [Fact]
    public void F38_WithoutBitATwin_StaysTrivialNotYetTyped()
    {
        var f38 = BuildF38();
        Assert.Null(f38.BitATwin);
        Assert.Equal(BitATwinClassification.TrivialNotYetTyped, f38.BitATwinStatus);
    }

    [Fact]
    public void F39_WithBitATwin_BitATwin_IsF39BitA()
    {
        var bitA = new F39DetPiBitAInheritance();
        var f39 = new F39DetPiPi2Inheritance(
            new Pi2DyadicLadderClaim(),
            new Pi2OperatorSpaceMirrorClaim(),
            bitA);
        Assert.Same(bitA, f39.BitATwin);
        Assert.Equal(BitATwinClassification.Filled, f39.BitATwinStatus);
    }

    [Fact]
    public void F39_WithoutBitATwin_StaysTrivialNotYetTyped()
    {
        var f39 = new F39DetPiPi2Inheritance(
            new Pi2DyadicLadderClaim(),
            new Pi2OperatorSpaceMirrorClaim());
        Assert.Null(f39.BitATwin);
        Assert.Equal(BitATwinClassification.TrivialNotYetTyped, f39.BitATwinStatus);
    }

    [Fact]
    public void F63_WithBitATwin_BitATwin_IsF63BitAReference()
    {
        var bitA = new F63BitAReference();
        var f38 = BuildF38();
        var f63 = new F63LCommutesPi2Pi2Inheritance(f38, new Pi2DyadicLadderClaim(), bitA);
        Assert.Same(bitA, f63.BitATwin);
        Assert.Equal(BitATwinClassification.Filled, f63.BitATwinStatus);
    }

    [Fact]
    public void F63_WithoutBitATwin_StaysTrivialNotYetTyped()
    {
        var f38 = BuildF38();
        var f63 = new F63LCommutesPi2Pi2Inheritance(f38, new Pi2DyadicLadderClaim());
        Assert.Null(f63.BitATwin);
        Assert.Equal(BitATwinClassification.TrivialNotYetTyped, f63.BitATwinStatus);
    }

    [Fact]
    public void XGlobalEigenstateMirror_WithBitATwin_BitATwin_IsZGlobalMirrorBitA()
    {
        var half = new HalfAsStructuralFixedPointClaim();
        var bitA = new ZGlobalEigenstateMirrorBitAInheritance(half);
        var x = new XGlobalEigenstateMirrorPi2Inheritance(half, bitA);
        Assert.Same(bitA, x.BitATwin);
        Assert.Equal(BitATwinClassification.Filled, x.BitATwinStatus);
    }

    [Fact]
    public void XGlobalEigenstateMirror_WithoutBitATwin_StaysTrivialNotYetTyped()
    {
        var half = new HalfAsStructuralFixedPointClaim();
        var x = new XGlobalEigenstateMirrorPi2Inheritance(half);
        Assert.Null(x.BitATwin);
        Assert.Equal(BitATwinClassification.TrivialNotYetTyped, x.BitATwinStatus);
    }

    // ------------------------------------------------------------------
    // A5-A8: 4 reclassifications now report BitBSpecific
    // ------------------------------------------------------------------

    [Fact]
    public void F82_BitATwinStatus_IsBitBSpecific()
    {
        var f82 = BuildF82();
        Assert.Equal(BitATwinClassification.BitBSpecific, ((IZ2AxisClaim)f82).BitATwinStatus);
    }

    [Fact]
    public void F84_BitATwinStatus_IsBitBSpecific()
    {
        var f82 = BuildF82();
        var f84 = new F84ThermalAmplitudeDampingPi2Inheritance(new Pi2DyadicLadderClaim(), f82);
        Assert.Equal(BitATwinClassification.BitBSpecific, ((IZ2AxisClaim)f84).BitATwinStatus);
    }

    private static F82T1AmplitudeDampingPi2Inheritance BuildF82()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f1Identity = new RCPsiSquared.Core.F1.F1PalindromeIdentity();
        var f81 = new F81Pi2Inheritance(f1Identity, ladder, mirror, memoryLoop);
        return new F82T1AmplitudeDampingPi2Inheritance(ladder, f81);
    }

    [Fact]
    public void F91_BitATwinStatus_IsBitBSpecific()
    {
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f91 = new F91Pi2Inheritance(memoryLoop);
        Assert.Equal(BitATwinClassification.BitBSpecific, ((IZ2AxisClaim)f91).BitATwinStatus);
    }

    [Fact]
    public void F93_BitATwinStatus_IsBitBSpecific()
    {
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f93 = new F93Pi2Inheritance(memoryLoop);
        Assert.Equal(BitATwinClassification.BitBSpecific, ((IZ2AxisClaim)f93).BitATwinStatus);
    }
}
