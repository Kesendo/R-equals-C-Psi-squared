using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F39DetPiPi2InheritanceTests
{
    private readonly ITestOutputHelper _out;

    public F39DetPiPi2InheritanceTests(ITestOutputHelper output) => _out = output;

    private static F39DetPiPi2Inheritance Build() =>
        new F39DetPiPi2Inheritance(new Pi2DyadicLadderClaim(), new Pi2OperatorSpaceMirrorClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        // F39 is Tier 1 proven (PT_SYMMETRY_ANALYSIS, verified N=1..4); the Pi2
        // anchoring is algebraic-trivial composition.
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Theory]
    [InlineData(1, 1.0)]      // N=1: 4^0 = 1 = a_1 (trivial identity scale)
    [InlineData(2, 4.0)]      // N=2: 4^1 = 4 = a_{-1} (= d² for 1 qubit)
    [InlineData(3, 16.0)]     // N=3: 4^2 = 16 = a_{-3} (= d² for 2 qubits)
    [InlineData(4, 64.0)]     // N=4: 4^3 = 64 = a_{-5} (= d² for 3 qubits)
    [InlineData(5, 256.0)]    // N=5: 4^4 = 256 = a_{-7} (= d² for 4 qubits)
    public void PowerNMinus1Factor_EqualsFourToTheNMinus1(int N, double expected)
    {
        var f = Build();
        Assert.Equal(expected, f.PowerNMinus1Factor(N), precision: 12);
    }

    [Theory]
    [InlineData(1, 1)]      // N=1: ladder index 3 - 2 = 1 (a_1 = 1)
    [InlineData(2, -1)]
    [InlineData(3, -3)]
    [InlineData(5, -7)]
    public void LadderIndexFor_EqualsThreeMinusTwoN(int N, int expected)
    {
        Assert.Equal(expected, Build().LadderIndexFor(N));
    }

    [Theory]
    [InlineData(1, 0)]      // N=1: qubit count 0 (degenerate; a_1 = 1 = trivial scale)
    [InlineData(2, 1)]      // N=2: qubit count 1
    [InlineData(3, 2)]      // N=3: qubit count 2
    [InlineData(5, 4)]
    public void OperatorSpaceQubitCountFor_EqualsNMinus1(int N, int expected)
    {
        Assert.Equal(expected, Build().OperatorSpaceQubitCountFor(N));
    }

    [Theory]
    [InlineData(1, 1L)]        // N=1: 1 · 4^0 = 1
    [InlineData(2, 8L)]        // N=2: 2 · 4 = 8
    [InlineData(3, 48L)]       // N=3: 3 · 16 = 48
    [InlineData(4, 256L)]      // N=4: 4 · 64 = 256
    [InlineData(5, 1280L)]     // N=5: 5 · 256 = 1280
    public void ExponentValue_EqualsNTimesFourPowerNMinus1(int N, long expected)
    {
        Assert.Equal(expected, Build().ExponentValue(N));
    }

    [Theory]
    [InlineData(1, false)]    // exponent 1 is odd
    [InlineData(2, true)]     // 8 is even
    [InlineData(3, true)]     // 48 is even
    [InlineData(4, true)]
    [InlineData(5, true)]
    public void ExponentIsEven_TrueIffNGreaterEqualTwo(int N, bool expected)
    {
        Assert.Equal(expected, Build().ExponentIsEven(N));
    }

    [Theory]
    [InlineData(1, -1)]    // det(Π) = -1 at N=1
    [InlineData(2, 1)]     // det(Π) = +1 for N ≥ 2
    [InlineData(3, 1)]
    [InlineData(4, 1)]
    [InlineData(5, 1)]
    public void DetPi_FollowsExponentParity(int N, int expected)
    {
        // F39 closed form: det(Π) = (-1)^(N·4^(N-1)). The Pi2 ladder anchor 4^(N-1)
        // being even for N ≥ 2 forces det = +1 there; at N=1 (4^0 = 1, exponent = 1)
        // det = -1.
        Assert.Equal(expected, Build().DetPi(N));
    }

    [Fact]
    public void PowerNMinus1Factor_NLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => Build().PowerNMinus1Factor(0));
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F39DetPiPi2Inheritance(null!, new Pi2OperatorSpaceMirrorClaim()));
        Assert.Throws<ArgumentNullException>(() =>
            new F39DetPiPi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }

    [Fact]
    public void Anchor_References_F39_AndPtSymmetryAnalysis_AndPi2Foundation()
    {
        var f = Build();
        Assert.Contains("ANALYTICAL_FORMULAS.md", f.Anchor);
        Assert.Contains("PT_SYMMETRY_ANALYSIS.md", f.Anchor);
        Assert.Contains("Pi2DyadicLadderClaim.cs", f.Anchor);
        Assert.Contains("Pi2OperatorSpaceMirrorClaim.cs", f.Anchor);
    }

    [Fact]
    public void Reconnaissance_EmitsDeterminantTable()
    {
        var f = Build();
        _out.WriteLine("");
        _out.WriteLine("    F39 closed form: det(Π) = (−1)^(N · 4^(N−1))");
        _out.WriteLine("    Pi2-Anker: 4^(N−1) = a_{3−2N} = d² für (N−1) qubits");
        _out.WriteLine("");
        _out.WriteLine("     N | 4^(N−1) | ladder | qubit-count | exponent N·4^(N−1) | even? | det(Π)");
        _out.WriteLine("    ---|---------|--------|-------------|--------------------|-------|--------");
        for (int N = 1; N <= 5; N++)
        {
            _out.WriteLine($"     {N} | {f.PowerNMinus1Factor(N),7:F0} | a_{f.LadderIndexFor(N),-5} | {f.OperatorSpaceQubitCountFor(N),-11} | {f.ExponentValue(N),18} | {f.ExponentIsEven(N),5} | {f.DetPi(N),+3}");
        }
    }
}
