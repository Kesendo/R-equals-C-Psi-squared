using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class ThreeLadderHingeClaimTests
{
    private static ThreeLadderHingeClaim MakeClaim()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var absorption = new AbsorptionTheoremClaim(ladder);

        // MomentTowerPumpChannelClaim's parent chain (mirrors MomentTowerPumpChannelClaimTests.MakeParents).
        var part2 = new F108Part2Pi2XEvenAlwaysPalindromic();
        var f113 = new LindbladBitBPiBreakMagnitude(
            new LindbladBitBPiBalance(
                new F108Part1Pi2EvenAlwaysPalindromic(part2),
                new LindbladBitAPiBalance(part2)));
        var f81 = new F81Pi2Inheritance(
            new F1PalindromeIdentity(),
            ladder,
            new Pi2OperatorSpaceMirrorClaim(),
            new Pi2I4MemoryLoopClaim());
        var f84 = new F84ThermalAmplitudeDampingPi2Inheritance(
            ladder,
            new F82T1AmplitudeDampingPi2Inheritance(ladder, f81));
        var moment = new MomentTowerPumpChannelClaim(f113, f84);

        return new ThreeLadderHingeClaim(absorption, moment);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, MakeClaim().Tier);
    }

    [Fact]
    public void Battery_AllCasesPass()
    {
        var claim = MakeClaim();
        Assert.NotEmpty(claim.Cases);
        foreach (var c in claim.Cases)
            Assert.True(c.Passes, $"battery case '{c.Name}' failed: expected {c.Expected}, got {c.Actual}");
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }

    [Fact]
    public void TypedParents_AreWired()
    {
        var claim = MakeClaim();
        Assert.NotNull(claim.Absorption);
        Assert.NotNull(claim.Moment);
    }

    [Fact]
    public void Anchor_ReferencesProofAndWitness()
    {
        var claim = MakeClaim();
        Assert.Contains("PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md", claim.Anchor);
        Assert.Contains("LadderHingeWitness", claim.Anchor);
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        var absorption = new AbsorptionTheoremClaim(new Pi2DyadicLadderClaim());
        Assert.Throws<ArgumentNullException>(() => new ThreeLadderHingeClaim(null!, null!));
        Assert.Throws<ArgumentNullException>(() => new ThreeLadderHingeClaim(absorption, null!));
    }
}
