using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F49bCenteredDissipatorPi2InheritanceTests
{
    private readonly ITestOutputHelper _out;

    public F49bCenteredDissipatorPi2InheritanceTests(ITestOutputHelper output) => _out = output;

    private static F49bCenteredDissipatorPi2Inheritance Build() =>
        new F49bCenteredDissipatorPi2Inheritance(new Pi2DyadicLadderClaim(), new Pi2OperatorSpaceMirrorClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Theory]
    [InlineData(1, 4.0)]      // N=1: 4^1 = 4 = a_{-1} (= d² for 1 qubit)
    [InlineData(2, 16.0)]     // N=2: 4^2 = 16 = a_{-3}
    [InlineData(3, 64.0)]     // N=3: 4^3 = 64 = a_{-5}
    [InlineData(4, 256.0)]    // N=4
    [InlineData(5, 1024.0)]   // N=5
    [InlineData(6, 4096.0)]   // N=6
    public void FourPowerNFactor_EqualsFourToTheN(int N, double expected)
    {
        var f = Build();
        Assert.Equal(expected, f.FourPowerNFactor(N), precision: 12);
    }

    [Theory]
    [InlineData(1, -1)]    // N=1: ladder index 1 - 2·1 = -1
    [InlineData(2, -3)]
    [InlineData(3, -5)]
    [InlineData(6, -11)]
    public void LadderIndexFor_EqualsOneMinusTwoN(int N, int expected)
    {
        Assert.Equal(expected, Build().LadderIndexFor(N));
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void MirrorPinnedFourPowerN_AgreesWithFourPowerNFactor(int N)
    {
        // F49b's 4^N IS the operator-space dimension d² for N qubits — direct anchor,
        // no derivation overhead. Cross-check between the ladder reading and the
        // pinned mirror table.
        var f = Build();
        Assert.Equal(f.FourPowerNFactor(N), f.MirrorPinnedFourPowerN(N), precision: 12);
    }

    [Theory]
    [InlineData(1, 4.0)]       // N=1: 4·1 = 4
    [InlineData(2, 32.0)]      // N=2: 16·2 = 32
    [InlineData(3, 192.0)]     // N=3: 64·3 = 192
    [InlineData(4, 1024.0)]    // N=4: 256·4 = 1024
    [InlineData(5, 5120.0)]    // N=5: 1024·5 = 5120
    public void LiveScalingAtUnitGamma_EqualsFourPowerNTimesN(int N, double expected)
    {
        var f = Build();
        Assert.Equal(expected, f.LiveScalingAtUnitGamma(N), precision: 10);
    }

    [Fact]
    public void FourPowerNFactor_NLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => Build().FourPowerNFactor(0));
    }

    [Fact]
    public void MirrorPinnedFourPowerN_OutsidePinnedRange_Throws()
    {
        // Pi2OperatorSpaceMirror table covers N=1..6; N=7 is outside.
        Assert.Throws<ArgumentOutOfRangeException>(() => Build().MirrorPinnedFourPowerN(7));
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F49bCenteredDissipatorPi2Inheritance(null!, new Pi2OperatorSpaceMirrorClaim()));
        Assert.Throws<ArgumentNullException>(() =>
            new F49bCenteredDissipatorPi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }

    [Fact]
    public void Anchor_References_F49b_AndCrossTermProof_AndPi2Foundation()
    {
        var f = Build();
        Assert.Contains("ANALYTICAL_FORMULAS.md", f.Anchor);
        Assert.Contains("PROOF_CROSS_TERM_FORMULA.md", f.Anchor);
        Assert.Contains("Pi2DyadicLadderClaim.cs", f.Anchor);
        Assert.Contains("Pi2OperatorSpaceMirrorClaim.cs", f.Anchor);
    }

    [Fact]
    public void Reconnaissance_EmitsScalingTable()
    {
        var f = Build();
        _out.WriteLine("");
        _out.WriteLine("    F49b closed form: ‖L_Dc‖² = γ² · 4^N · N (Tier 1 proven)");
        _out.WriteLine("    Pi2-Anker: 4^N = a_{1−2N} = d² für N qubits (direkt, no shift)");
        _out.WriteLine("");
        _out.WriteLine("     N | 4^N    | ladder    | scaling 4^N·N | OperatorSpaceMirror d²(N)");
        _out.WriteLine("    ---|--------|-----------|---------------|---------------------------");
        for (int N = 1; N <= 6; N++)
        {
            _out.WriteLine($"     {N} | {f.FourPowerNFactor(N),6:F0} | a_{f.LadderIndexFor(N),-8} | {f.LiveScalingAtUnitGamma(N),13:F0} | {f.MirrorPinnedFourPowerN(N),6:F0}");
        }
        _out.WriteLine("");
        _out.WriteLine("    F-formula 4^k power-shift comparison:");
        _out.WriteLine("      F49  uses 4^(N−2) = a_{5−2N}    (residual scaling)");
        _out.WriteLine("      F39  uses 4^(N−1) = a_{3−2N}    (det Π exponent)");
        _out.WriteLine("      F1-T1 uses 4^(N−1) = a_{3−2N}   (T1-part prefactor)");
        _out.WriteLine("      F49b uses 4^N     = a_{1−2N}    (full operator-space, no shift)");
    }
}
