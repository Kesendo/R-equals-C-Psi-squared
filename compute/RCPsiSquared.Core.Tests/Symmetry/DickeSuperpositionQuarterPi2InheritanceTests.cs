using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class DickeSuperpositionQuarterPi2InheritanceTests
{
    private readonly ITestOutputHelper _out;

    public DickeSuperpositionQuarterPi2InheritanceTests(ITestOutputHelper output) => _out = output;

    private static DickeSuperpositionQuarterPi2Inheritance Build() =>
        new DickeSuperpositionQuarterPi2Inheritance(
            new Pi2DyadicLadderClaim(),
            new QuarterAsBilinearMaxvalClaim(),
            new HalfAsStructuralFixedPointClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        // Theorem 1 + Theorem 2 are algebraic identities (Cauchy-Schwarz + AM-GM +
        // Dicke amplitude uniformity). Pure Tier1Derived F86 closed form.
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Fact]
    public void QuarterCeiling_IsExactlyOneQuarter_FromLadderTermThree()
    {
        // The 1/4 ceiling is exactly a_3 on the Pi2 dyadic ladder
        // (= QuarterAsBilinearMaxval = (1/d)² for d=2).
        var d = Build();
        Assert.Equal(0.25, d.QuarterCeiling, precision: 14);
    }

    [Fact]
    public void SectorBalance_IsExactlyOneHalf_FromLadderTermTwo()
    {
        // The AM-GM saturation balance p_n = p_{n+1} = 1/2 is exactly a_2 on the Pi2
        // dyadic ladder (= HalfAsStructuralFixedPoint = 1/d for d=2).
        var d = Build();
        Assert.Equal(0.5, d.SectorBalance, precision: 14);
    }

    [Fact]
    public void SectorBalanceSquared_EqualsQuarterCeiling()
    {
        // The AM-GM identity: (1/2)² = 1/4. The 1/4 ceiling is the square of the
        // sector balance — direct manifest of "1/4 = (1/2)² unifies framework
        // quarter-boundaries" (memory project_quarter_as_polarity_squared).
        var d = Build();
        Assert.Equal(d.QuarterCeiling, d.SectorBalanceSquared, precision: 14);
    }

    [Theory]
    [InlineData(3, 0, 3)]      // C(3,0)·C(3,1) = 1·3 = 3
    [InlineData(3, 1, 9)]      // C(3,1)·C(3,2) = 3·3 = 9
    [InlineData(4, 1, 24)]     // C(4,1)·C(4,2) = 4·6 = 24
    [InlineData(5, 2, 100)]    // C(5,2)·C(5,3) = 10·10 = 100
    [InlineData(6, 2, 300)]    // C(6,2)·C(6,3) = 15·20 = 300
    public void BlockDimension_EqualsBinomialProduct(int N, int n, long expected)
    {
        Assert.Equal(expected, Build().BlockDimension(N, n));
    }

    [Theory]
    [InlineData(3, 0)]
    [InlineData(3, 1)]
    [InlineData(4, 1)]
    [InlineData(5, 2)]
    [InlineData(6, 2)]
    [InlineData(8, 3)]
    public void LiveBlockCpsiAtZero_EqualsExactlyOneQuarter(int N, int n)
    {
        // Theorem 1 mechanism: M_block · |amplitude|² = M_block · 1/(4·M_block) = 1/4.
        // The M_block factor cancels exactly. Live drift check.
        var d = Build();
        Assert.Equal(0.25, d.LiveBlockCpsiAtZero(N, n), precision: 12);
    }

    [Theory]
    [InlineData(3, 0)]   // C_block = 1/4 at N=3, n=0
    [InlineData(8, 3)]   // C_block = 1/4 at N=8, n=3
    public void LiveBlockCpsiAtZero_EqualsQuarterCeiling(int N, int n)
    {
        // Cross-verification: live block-CΨ at t=0 equals the Pi2 ladder's a_3 anchor.
        var d = Build();
        Assert.Equal(d.QuarterCeiling, d.LiveBlockCpsiAtZero(N, n), precision: 12);
    }

    [Fact]
    public void BlockDimension_NLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => Build().BlockDimension(0, 0));
    }

    [Fact]
    public void BlockDimension_NPlusOneOutOfRange_Throws()
    {
        // n+1 > N is invalid (the (n+1) sector does not exist).
        Assert.Throws<ArgumentOutOfRangeException>(() => Build().BlockDimension(5, 5));
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var half = new HalfAsStructuralFixedPointClaim();
        Assert.Throws<ArgumentNullException>(() =>
            new DickeSuperpositionQuarterPi2Inheritance(null!, quarter, half));
        Assert.Throws<ArgumentNullException>(() =>
            new DickeSuperpositionQuarterPi2Inheritance(ladder, null!, half));
        Assert.Throws<ArgumentNullException>(() =>
            new DickeSuperpositionQuarterPi2Inheritance(ladder, quarter, null!));
    }

    [Fact]
    public void TypedParents_AreExposed()
    {
        var d = Build();
        Assert.NotNull(d.Quarter);
        Assert.NotNull(d.Half);
    }

    [Fact]
    public void Anchor_References_BlockCpsiQuarterProof_AndPi2Foundation()
    {
        var d = Build();
        Assert.Contains("PROOF_BLOCK_CPSI_QUARTER.md", d.Anchor);
        Assert.Contains("Pi2KnowledgeBaseClaims.cs", d.Anchor);
        Assert.Contains("Pi2DyadicLadderClaim.cs", d.Anchor);
    }

    [Fact]
    public void Reconnaissance_EmitsTheoremTable()
    {
        var d = Build();
        _out.WriteLine("");
        _out.WriteLine("    Theorem 1: C_block(0) = 1/4 on (|D_n⟩+|D_{n+1}⟩)/√2 (Tier1Derived)");
        _out.WriteLine("    Theorem 2: 1/4 is universal block-purity ceiling (Cauchy-Schwarz + AM-GM)");
        _out.WriteLine("");
        _out.WriteLine($"    QuarterCeiling     = a_3 = {d.QuarterCeiling}");
        _out.WriteLine($"    SectorBalance      = a_2 = {d.SectorBalance}");
        _out.WriteLine($"    SectorBalance²     = {d.SectorBalanceSquared} (= QuarterCeiling)");
        _out.WriteLine("");
        _out.WriteLine("     N | n | M_block | |amp|²       | C_block(0)");
        _out.WriteLine("    ---|---|---------|--------------|------------");
        foreach (var (N, n) in new[] { (3, 0), (4, 1), (5, 2), (6, 2), (8, 3) })
        {
            long mb = d.BlockDimension(N, n);
            double amp = d.PerEntryAmplitude(N, n);
            double live = d.LiveBlockCpsiAtZero(N, n);
            _out.WriteLine($"     {N} | {n} | {mb,7} | {amp * amp,12:G6} | {live:F4}");
        }
    }
}
