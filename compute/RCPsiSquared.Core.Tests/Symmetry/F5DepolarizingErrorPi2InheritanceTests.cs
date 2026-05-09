using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F5DepolarizingErrorPi2InheritanceTests
{
    private readonly ITestOutputHelper _out;

    public F5DepolarizingErrorPi2InheritanceTests(ITestOutputHelper output) => _out = output;

    private static F5DepolarizingErrorPi2Inheritance Build() =>
        new F5DepolarizingErrorPi2Inheritance(new Pi2DyadicLadderClaim(), new Pi2OperatorSpaceMirrorClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        // F5 is Tier 1 proven; the Pi2-Foundation anchoring is algebraic-trivial
        // composition. Tier1Derived inheritance.
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Fact]
    public void DCoefficient_IsExactlyTwo_FromLadderTermZero()
    {
        // The "2" in F5's 2·(N−2)/3 is a_0 = d (qubit dimension).
        var f = Build();
        Assert.Equal(2.0, f.DCoefficient, precision: 14);
    }

    [Fact]
    public void DSquaredMinusOne_IsExactlyThree()
    {
        // The "3" denominator: d² − 1 = 4 − 1 = 3 (number of non-identity Paulis).
        // Pi2-derived from operator-space-mirror N=1 anchor.
        var f = Build();
        Assert.Equal(3.0, f.DSquaredMinusOne, precision: 14);
    }

    [Fact]
    public void TwoOverThree_IsExactly2DividedBy3()
    {
        // The composite Pi2-derived ratio: d / (d² − 1) = 2/3.
        var f = Build();
        Assert.Equal(2.0 / 3.0, f.TwoOverThree, precision: 14);
    }

    [Theory]
    [InlineData(2, 0)]    // N=2: (N−2) = 0; F5 error vanishes (degenerate chain)
    [InlineData(3, 1)]
    [InlineData(4, 2)]
    [InlineData(5, 3)]
    [InlineData(8, 6)]
    public void NShiftFactor_EqualsNMinusTwo(int N, int expected)
    {
        Assert.Equal(expected, Build().NShiftFactor(N));
    }

    [Theory]
    [InlineData(2, 0.0)]              // N=2: 2·0/3 = 0
    [InlineData(3, 2.0 / 3.0)]        // N=3: 2·1/3 = 2/3
    [InlineData(4, 4.0 / 3.0)]        // N=4: 2·2/3 = 4/3
    [InlineData(5, 6.0 / 3.0)]        // N=5: 2·3/3 = 2
    [InlineData(6, 8.0 / 3.0)]        // N=6: 2·4/3 = 8/3
    public void LiveCoefficient_EqualsClosedForm(int N, double expected)
    {
        // Cross-verification: live composition of the Pi2-anchored constants matches
        // F5's closed form 2·(N−2)/3 bit-exact.
        var f = Build();
        Assert.Equal(expected, f.LiveCoefficient(N), precision: 12);
    }

    [Fact]
    public void LiveCoefficient_NLessThanTwo_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => Build().LiveCoefficient(1));
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F5DepolarizingErrorPi2Inheritance(null!, new Pi2OperatorSpaceMirrorClaim()));
        Assert.Throws<ArgumentNullException>(() =>
            new F5DepolarizingErrorPi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }

    [Fact]
    public void Anchor_References_F5_AndDepolarizingPalindromeExperiment_AndPi2Foundation()
    {
        var f = Build();
        Assert.Contains("ANALYTICAL_FORMULAS.md", f.Anchor);
        Assert.Contains("DEPOLARIZING_PALINDROME.md", f.Anchor);
        Assert.Contains("F1PalindromeIdentity.cs", f.Anchor);
        Assert.Contains("Pi2DyadicLadderClaim.cs", f.Anchor);
        Assert.Contains("Pi2OperatorSpaceMirrorClaim.cs", f.Anchor);
    }

    [Fact]
    public void Reconnaissance_EmitsDecomposition()
    {
        var f = Build();
        _out.WriteLine("");
        _out.WriteLine("    F5 closed form: error = γ · 2·(N−2)/3 (Tier 1 proven)");
        _out.WriteLine("");
        _out.WriteLine($"    \"2\" multiplier  = a_0 = d           = {f.DCoefficient}");
        _out.WriteLine($"    \"3\" denominator = a_(-1) - 1 = d²-1 = {f.DSquaredMinusOne}");
        _out.WriteLine($"    \"2/3\" ratio    = d/(d²-1)          = {f.TwoOverThree:F6}");
        _out.WriteLine("");
        _out.WriteLine("    interpretation: 2 = off-diag Paulis (X,Y); 3 = total non-identity Paulis (X,Y,Z)");
        _out.WriteLine("                    2/3 = off-diag fraction in depolarizing channel");
        _out.WriteLine("");
        _out.WriteLine("     N | (N−2) | coefficient");
        _out.WriteLine("    ---|-------|------------");
        for (int N = 2; N <= 6; N++)
            _out.WriteLine($"     {N} |   {f.NShiftFactor(N)}   | {f.LiveCoefficient(N):F4}");
    }
}
