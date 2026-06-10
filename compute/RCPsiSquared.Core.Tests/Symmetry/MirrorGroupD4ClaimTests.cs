using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class MirrorGroupD4ClaimTests
{
    private static MirrorGroupD4Claim MakeClaim()
    {
        var kleinV4 = new Pi2KleinV4DephaseSwapGroup();
        return new MirrorGroupD4Claim(
            new KleinEightCellClaim(new KleinFourCellClaim()),
            new CommutatorDConjugationSign(kleinV4),
            kleinV4);
    }

    // ------------------------------------------------------------------
    // Claim metadata
    // ------------------------------------------------------------------

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var claim = MakeClaim();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void Anchor_ReferencesProofAndVerifier()
    {
        var claim = MakeClaim();
        Assert.Contains("PROOF_PI_FACTORS_AS_R_TIMES_D.md", claim.Anchor);
        Assert.Contains("mirror_inventory_d4.py", claim.Anchor);
        Assert.Contains("PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md", claim.Anchor);
        Assert.Contains("F114", claim.Anchor);
    }

    [Fact]
    public void TypedParents_AreWired()
    {
        var claim = MakeClaim();
        Assert.NotNull(claim.Cube);
        Assert.NotNull(claim.F114);
        Assert.NotNull(claim.KleinV4);
        Assert.Same(claim.KleinV4, claim.F114.KleinV4);
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        var klein4 = new KleinFourCellClaim();
        var cube = new KleinEightCellClaim(klein4);
        var kleinV4 = new Pi2KleinV4DephaseSwapGroup();
        var f114 = new CommutatorDConjugationSign(kleinV4);

        Assert.Throws<ArgumentNullException>(() => new MirrorGroupD4Claim(null!, f114, kleinV4));
        Assert.Throws<ArgumentNullException>(() => new MirrorGroupD4Claim(cube, null!, kleinV4));
        Assert.Throws<ArgumentNullException>(() => new MirrorGroupD4Claim(cube, f114, null!));
    }

    // ------------------------------------------------------------------
    // Self-check battery (N = 2 dense superoperators, built in the ctor)
    // ------------------------------------------------------------------

    [Fact]
    public void Battery_AllCasesPass()
    {
        var claim = MakeClaim();
        Assert.Equal(10, claim.Cases.Count);
        foreach (var c in claim.Cases)
            Assert.True(c.Passes, $"battery case '{c.Name}' failed: expected {c.Expected}, got {c.Actual}");
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }

    [Fact]
    public void Battery_GroupClosure_IsEight()
    {
        var claim = MakeClaim();
        var closure = claim.Cases.Single(c => c.Name.Contains("|⟨R, D⟩| = 8"));
        Assert.Equal("8", closure.Expected);
        Assert.Equal("8", closure.Actual);
    }

    [Fact]
    public void Battery_Factorization_IsExact()
    {
        var claim = MakeClaim();
        var fact = claim.Cases.Single(c => c.Name.Contains("Π_Z = R·D"));
        Assert.True(fact.Passes, $"factorization dev out of tolerance: {fact.Actual}");
    }

    [Fact]
    public void Battery_SpineV4_IsKleinSubgroup()
    {
        var claim = MakeClaim();
        var spine = claim.Cases.Single(c => c.Name.Contains("spine V₄"));
        Assert.Equal("4/4", spine.Actual);
    }

    [Fact]
    public void Battery_CubeCharacters_AllSixteenStrings()
    {
        var claim = MakeClaim();
        foreach (var axis in new[] { "bit_a = char(Ad_{Z^⊗N})", "bit_b = char(Ad_{X^⊗N})", "y_par = char(θ)" })
        {
            var c = claim.Cases.Single(x => x.Name.StartsWith(axis, StringComparison.Ordinal));
            Assert.Equal("16/16", c.Actual);
        }
    }

    [Fact]
    public void Summary_CarriesTheFactorizationAndTheAntiautomorphismAxis()
    {
        var claim = MakeClaim();
        Assert.Contains("Π_Z = R·D", claim.Summary);
        Assert.Contains("antiautomorphism axis", claim.Summary);
        Assert.Contains($"{claim.Cases.Count}/{claim.Cases.Count} battery PASS", claim.Summary);
    }
}
