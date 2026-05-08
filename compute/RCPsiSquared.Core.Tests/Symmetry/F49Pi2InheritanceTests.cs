using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F49Pi2InheritanceTests
{
    private readonly ITestOutputHelper _out;

    public F49Pi2InheritanceTests(ITestOutputHelper output) => _out = output;

    private static F49Pi2Inheritance Build() =>
        new F49Pi2Inheritance(new Pi2DyadicLadderClaim(), new Pi2OperatorSpaceMirrorClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Theory]
    [InlineData(2, 1.0)]    // N=2 → 4^0 = 1 = a_1 (the trivial identity scale)
    [InlineData(3, 4.0)]    // N=3 → 4^1 = 4 = a_{-1} = d² for 1 qubit
    [InlineData(4, 16.0)]   // N=4 → 4^2 = 16 = a_{-3} = d² for 2 qubits
    [InlineData(5, 64.0)]   // N=5 → 4^3 = 64 = a_{-5} = d² for 3 qubits
    [InlineData(6, 256.0)]  // N=6 → 4^4 = 256 = a_{-7} = d² for 4 qubits
    [InlineData(7, 1024.0)] // N=7 → 4^5 = 1024 = a_{-9} = d² for 5 qubits
    [InlineData(8, 4096.0)] // N=8 → 4^6 = 4096 = a_{-11} = d² for 6 qubits
    public void PowerFactor_EqualsFourToTheNMinusTwo(int chainN, double expected)
    {
        var f = Build();
        Assert.Equal(expected, f.PowerFactor(chainN), precision: 12);
    }

    [Theory]
    [InlineData(2, 1)]    // N=2 → ladder index 1 (self-mirror pivot, a_1=1)
    [InlineData(3, -1)]   // N=3 → ladder index -1 (a_{-1}=4)
    [InlineData(4, -3)]   // N=4 → ladder index -3 (a_{-3}=16)
    [InlineData(5, -5)]
    [InlineData(8, -11)]
    public void LadderIndexFor_EqualsFiveMinusTwoN(int chainN, int expected)
    {
        var f = Build();
        Assert.Equal(expected, f.LadderIndexFor(chainN));
    }

    [Theory]
    [InlineData(3, 1)]    // N_chain=3 → qubit count 1
    [InlineData(4, 2)]    // N_chain=4 → qubit count 2
    [InlineData(5, 3)]
    [InlineData(8, 6)]
    public void OperatorSpaceQubitCountFor_EqualsChainNMinusTwo(int chainN, int expected)
    {
        var f = Build();
        Assert.Equal(expected, f.OperatorSpaceQubitCountFor(chainN));
    }

    [Theory]
    [InlineData(3)]   // chainN=3 → qubit count 1, in mirror table (N=1..6)
    [InlineData(4)]
    [InlineData(8)]   // chainN=8 → qubit count 6, in mirror table
    public void IsInPinnedMirrorTable_TrueForChainNThreeToEight(int chainN)
    {
        var f = Build();
        Assert.True(f.IsInPinnedMirrorTable(chainN));
    }

    [Fact]
    public void IsInPinnedMirrorTable_FalseForChainNNine()
    {
        // chainN=9 → qubit count 7, outside the pinned table (currently N=1..6)
        var f = Build();
        Assert.False(f.IsInPinnedMirrorTable(9));
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void MirrorPinnedPowerFactor_AgreesWithPowerFactor(int chainN)
    {
        // Cross-verification: the operator-space dimension pinned in the mirror table
        // equals the F49 power factor 4^(N−2) computed via the ladder. Drift between
        // the two surfaces here.
        var f = Build();
        Assert.Equal(f.PowerFactor(chainN), f.MirrorPinnedPowerFactor(chainN), precision: 12);
    }

    [Fact]
    public void PowerFactor_ChainNLessThanTwo_Throws()
    {
        var f = Build();
        Assert.Throws<ArgumentOutOfRangeException>(() => f.PowerFactor(1));
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F49Pi2Inheritance(null!, new Pi2OperatorSpaceMirrorClaim()));
        Assert.Throws<ArgumentNullException>(() =>
            new F49Pi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }

    [Fact]
    public void Anchor_References_OperatorRigidity_AndPalindromeResidualScaling_AndPi2Foundation()
    {
        var f = Build();
        Assert.Contains("OPERATOR_RIGIDITY_ACROSS_CUSP.md", f.Anchor);
        Assert.Contains("PalindromeResidualScalingClaim.cs", f.Anchor);
        Assert.Contains("Pi2DyadicLadderClaim.cs", f.Anchor);
        Assert.Contains("Pi2OperatorSpaceMirrorClaim.cs", f.Anchor);
    }

    [Fact]
    public void Reconnaissance_EmitsInheritanceTable()
    {
        var f = Build();
        _out.WriteLine("");
        _out.WriteLine("    F49: ‖M(N, G)‖²_F = c_H · B(G) · 4^(N−2)");
        _out.WriteLine("    inheritance: 4^(N−2) = d² for (N−2) qubits");
        _out.WriteLine("");
        _out.WriteLine("    chain N | 4^(N−2) | ladder index | qubit count | in pinned table");
        _out.WriteLine("    --------|---------|--------------|-------------|----------------");
        for (int chainN = 2; chainN <= 9; chainN++)
        {
            string powerStr = chainN >= 2 ? $"{f.PowerFactor(chainN),6:F0}" : "(N/A)";
            _out.WriteLine($"    {chainN,7} | {powerStr} | a_{f.LadderIndexFor(chainN),-9} | {f.OperatorSpaceQubitCountFor(chainN),-11} | {f.IsInPinnedMirrorTable(chainN)}");
        }
    }
}
