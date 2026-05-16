using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Verifies the F89 ↔ F88 bridge anchors: XX, YY bond terms in F88-Pp,
/// ρ_cc initial state classified by F88PopcountPairLens(N, 1, 2).</summary>
public class F89F88KleinPpAnchorTests
{
    [Fact]
    public void F89XxEntry_KleinCell_IsPp()
    {
        Assert.Equal("Pp", F89F88KleinPpAnchor.F89XxEntry.Cell);
        Assert.Equal(+1, F89F88KleinPpAnchor.F89XxEntry.Pi2ZEigenvalue);
        Assert.Equal(+1, F89F88KleinPpAnchor.F89XxEntry.Pi2XEigenvalue);
    }

    [Fact]
    public void F89YyEntry_KleinCell_IsPp()
    {
        Assert.Equal("Pp", F89F88KleinPpAnchor.F89YyEntry.Cell);
        Assert.Equal(+1, F89F88KleinPpAnchor.F89YyEntry.Pi2ZEigenvalue);
        Assert.Equal(+1, F89F88KleinPpAnchor.F89YyEntry.Pi2XEigenvalue);
    }

    [Fact]
    public void Constants_MatchKleinPpEigenvalues()
    {
        Assert.Equal("Pp", F89F88KleinPpAnchor.F89BondKleinCell);
        Assert.Equal(+1, F89F88KleinPpAnchor.F89BondPi2ZEigenvalue);
        Assert.Equal(+1, F89F88KleinPpAnchor.F89BondPi2XEigenvalue);
        Assert.Equal(1, F89F88KleinPpAnchor.F89InitialStateNp);
        Assert.Equal(2, F89F88KleinPpAnchor.F89InitialStateNq);
    }

    [Fact]
    public void InitialStateLens_AtN3_IsPopcountMirror()
    {
        var lens = F89F88KleinPpAnchor.InitialStateLens(N: 3);
        Assert.Equal(ConfigurationKind.PopcountMirror, lens.Kind);
        Assert.Equal(0.0, lens.Alpha, precision: 12);
    }

    [Fact]
    public void InitialStateLens_AtN4_IsKIntermediate()
    {
        var lens = F89F88KleinPpAnchor.InitialStateLens(N: 4);
        Assert.Equal(ConfigurationKind.KIntermediate, lens.Kind);
        // α = C(4,2) / (2·(C(4,1) + C(4,2))) = 6 / (2·(4+6)) = 6/20 = 0.3
        Assert.Equal(0.3, lens.Alpha, precision: 12);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(8)]
    public void InitialStateLens_AtNAtLeast5_IsGeneric(int N)
    {
        var lens = F89F88KleinPpAnchor.InitialStateLens(N);
        Assert.Equal(ConfigurationKind.Generic, lens.Kind);
        Assert.Equal(0.5, lens.Alpha, precision: 12);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var anchor = new F89F88KleinPpAnchor(
            new KleinFourCellClaim(),
            new F89TopologyOrbitClosure(new Pi2DyadicLadderClaim()));
        Assert.Equal(Tier.Tier1Derived, anchor.Tier);
    }

    [Fact]
    public void Constructor_NullKleinFourCell_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => new F89F88KleinPpAnchor(
            null!,
            new F89TopologyOrbitClosure(new Pi2DyadicLadderClaim())));
    }

    [Fact]
    public void Constructor_NullF89_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => new F89F88KleinPpAnchor(
            new KleinFourCellClaim(),
            null!));
    }

    [Fact]
    public void OtherBondTerms_NonPp_AreCorrectlyDifferentCells()
    {
        // Sanity: F87-Truly XY+YX (both Mp) should NOT land in Pp.
        var xy = new KleinBilinearEntry(PauliLetter.X, PauliLetter.Y);
        var yz = new KleinBilinearEntry(PauliLetter.Y, PauliLetter.Z);
        Assert.NotEqual("Pp", xy.Cell);
        Assert.NotEqual("Pp", yz.Cell);
    }
}
