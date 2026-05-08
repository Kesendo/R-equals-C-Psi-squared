using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class Pi2OperatorSpaceMirrorClaimTests
{
    private readonly ITestOutputHelper _out;

    public Pi2OperatorSpaceMirrorClaimTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        Assert.Equal(Tier.Tier1Derived, mirror.Tier);
    }

    [Fact]
    public void Pairs_HasSixEntriesForNOneToSix()
    {
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        Assert.Equal(6, mirror.Pairs.Count);
        Assert.Equal(new[] { 1, 2, 3, 4, 5, 6 }, mirror.Pairs.Select(p => p.N).ToArray());
    }

    [Theory]
    [InlineData(1, -1, 4.0,    3, 0.25)]
    [InlineData(2, -3, 16.0,   5, 0.0625)]
    [InlineData(3, -5, 64.0,   7, 1.0 / 64.0)]
    [InlineData(4, -7, 256.0,  9, 1.0 / 256.0)]
    [InlineData(5, -9, 1024.0, 11, 1.0 / 1024.0)]
    [InlineData(6, -11, 4096.0, 13, 1.0 / 4096.0)]
    public void Pair_HasCorrectIndicesAndValues(int N, int upperIdx, double opSpace,
        int lowerIdx, double mirrorMass)
    {
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        var pair = mirror.PairAt(N);
        Assert.NotNull(pair);
        Assert.Equal(upperIdx, pair!.UpperIndex);
        Assert.Equal(opSpace, pair.OperatorSpace, precision: 12);
        Assert.Equal(lowerIdx, pair.LowerIndex);
        Assert.Equal(mirrorMass, pair.MirrorMass, precision: 12);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void MirrorProduct_AlwaysEqualsOne(int N)
    {
        // The Tier1Derived inversion identity 4^N · 1/4^N = 1, surfaced as a typed
        // Schicht-1-fact per qubit count.
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        var pair = mirror.PairAt(N)!;
        Assert.Equal(1.0, pair.OperatorSpace * pair.MirrorMass, precision: 12);
        Assert.Equal(1.0, mirror.LiveMirrorProduct(N), precision: 12);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void Pair_AgreesWithLadderTermsAtBothIndices(int N)
    {
        // Cross-verification: the pinned values must match Pi2DyadicLadder.Term at the
        // same indices. Drift between the pair and the ladder is caught here.
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        var ladder = new Pi2DyadicLadderClaim();
        var pair = mirror.PairAt(N)!;

        Assert.Equal(ladder.Term(pair.UpperIndex), pair.OperatorSpace, precision: 12);
        Assert.Equal(ladder.Term(pair.LowerIndex), pair.MirrorMass, precision: 12);
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(6)]
    public void Pair_LadderIndicesAreMirrorPartners(int N)
    {
        // The two ladder indices in each pair must be mirror-partners under the ladder's
        // inversion symmetry: 2 - UpperIndex = LowerIndex (and vice versa).
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        var ladder = new Pi2DyadicLadderClaim();
        var pair = mirror.PairAt(N)!;

        Assert.Equal(pair.LowerIndex, ladder.MirrorPartnerIndex(pair.UpperIndex));
        Assert.Equal(pair.UpperIndex, ladder.MirrorPartnerIndex(pair.LowerIndex));
    }

    [Fact]
    public void OperatorSpace_EqualsFourPowerN()
    {
        // Literal reading: each upper-side value is 4^N. Pauli basis cardinality.
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        foreach (var p in mirror.Pairs)
            Assert.Equal(Math.Pow(4.0, p.N), p.OperatorSpace, precision: 10);
    }

    [Fact]
    public void MirrorMass_EqualsOneOverFourPowerN()
    {
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        foreach (var p in mirror.Pairs)
            Assert.Equal(1.0 / Math.Pow(4.0, p.N), p.MirrorMass, precision: 14);
    }

    [Fact]
    public void Anchor_References_MirrorTheory_AndZeroIsTheMirror()
    {
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        Assert.Contains("MIRROR_THEORY.md", mirror.Anchor);
        Assert.Contains("ZERO_IS_THE_MIRROR.md", mirror.Anchor);
    }

    [Fact]
    public void Reconnaissance_EmitsMirrorTable()
    {
        // Documents the per-N stride-2 mirror substructure of the dyadic ladder as a
        // queryable Schicht-1 table.
        var mirror = new Pi2OperatorSpaceMirrorClaim();
        _out.WriteLine("  N | upper a_{-(2N-1)} | OperatorSpace 4^N | lower a_{2N+1} | MirrorMass 1/4^N | product");
        _out.WriteLine("  --|-------------------|-------------------|-----------------|-------------------|--------");
        foreach (var p in mirror.Pairs)
        {
            _out.WriteLine($"  {p.N} | a_{p.UpperIndex,-15} | {p.OperatorSpace,17:F0} | a_{p.LowerIndex,-13} | {p.MirrorMass,17:G8} | {p.OperatorSpace * p.MirrorMass:F1}");
        }
    }
}
