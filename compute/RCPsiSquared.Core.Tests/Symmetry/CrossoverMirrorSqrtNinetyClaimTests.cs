using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class CrossoverMirrorSqrtNinetyClaimTests
{
    private static CrossoverMirrorSqrtNinetyClaim BuildClaim() =>
        new CrossoverMirrorSqrtNinetyClaim(new NinetyDegreeMirrorMemoryClaim());

    [Fact]
    public void Tier_IsTier1Candidate()
    {
        // Bit-exact numerical witness; analytical derivation of S_light²=σ_x↔σ_y is the
        // promotion step, so the claim is a Tier-1 candidate, not yet derived.
        Assert.Equal(Tier.Tier1Candidate, BuildClaim().Tier);
    }

    [Fact]
    public void Turn_S_IsBlockDiagonal_PureRotationNoDarkLightSwap()
    {
        Assert.True(BuildClaim().TurnIsBlockDiagonal());
    }

    [Fact]
    public void LightPlaneSquare_IsTheNinetyDegreeRotation_BitExact()
    {
        // The heart: S_light² = [[0,−1],[1,0]] = σ_x↔σ_y 90° = the NinetyDegreeMirror.
        // So S_light = √(NinetyDegreeMirror); the crossover mirror is Π turned by half the
        // 90° angle-anchor.
        var claim = BuildClaim();
        Assert.True(claim.LightPlaneSquareIsNinetyDegree());
        Assert.True(claim.LightPlaneSquareResidual() < 1e-9,
            $"S_light² should equal the σ_x↔σ_y 90° rotation; residual {claim.LightPlaneSquareResidual():E2}");
    }

    [Fact]
    public void CrossoverMirror_IsOrderFour_MSquaredIsMinusIdentity()
    {
        Assert.True(BuildClaim().CrossoverMirrorIsOrderFour());
    }

    [Fact]
    public void Constructor_RejectsNullParent()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new CrossoverMirrorSqrtNinetyClaim(null!));
    }

    [Fact]
    public void TypedParent_NinetyDegree_IsExposed()
    {
        Assert.NotNull(BuildClaim().NinetyDegree);
    }

    [Fact]
    public void Anchor_References_LocalityResultParentAndReflection()
    {
        var f = BuildClaim();
        Assert.Contains("PI_OPERATOR_ENTANGLEMENT.md", f.Anchor);
        Assert.Contains("crossover_mirror_sqrt_ninety.py", f.Anchor);
        Assert.Contains("ON_THE_SQUARE_ROOT_OF_THE_MIRROR.md", f.Anchor);
        Assert.Contains("NinetyDegreeMirrorMemory", f.Anchor);
    }
}
