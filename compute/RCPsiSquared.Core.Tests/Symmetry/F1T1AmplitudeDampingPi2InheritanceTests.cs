using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F1T1AmplitudeDampingPi2InheritanceTests
{
    private readonly ITestOutputHelper _out;

    public F1T1AmplitudeDampingPi2InheritanceTests(ITestOutputHelper output) => _out = output;

    private static F1T1AmplitudeDampingPi2Inheritance Build() =>
        new F1T1AmplitudeDampingPi2Inheritance(new Pi2DyadicLadderClaim(), new Pi2OperatorSpaceMirrorClaim());

    [Fact]
    public void Tier_IsTier1Candidate()
    {
        // The closed form is empirically verified N=3..6 but the analytic derivation is
        // open (F1 open-question Item 2). The Pi2-Foundation anchoring is conditional on
        // the formula being correct; Tier1Candidate matches the parent formula's tier.
        Assert.Equal(Tier.Tier1Candidate, Build().Tier);
    }

    [Fact]
    public void FourMultiplier_IsExactlyFour_FromLadderTermMinusOne()
    {
        // The "4" multiplier in 4·(Σγ_T1)² is exactly a_{-1} = 4 = d² for N=1 qubit.
        var f = Build();
        Assert.Equal(4.0, f.FourMultiplier, precision: 14);
    }

    [Fact]
    public void ThreeMultiplier_IsExactlyThree()
    {
        // The "3" in 3·Σγ_T1² — small integer from T1 dissipator algebra. Documented
        // but not claimed as Pi2-Foundation anchor.
        Assert.Equal(3, F1T1AmplitudeDampingPi2Inheritance.ThreeMultiplier);
    }

    [Theory]
    [InlineData(2, 16.0)]    // N=2: 2^(2+2) = 2^4 = 16 = a_{-3}
    [InlineData(3, 32.0)]    // N=3: 2^5 = 32 = a_{-4}
    [InlineData(4, 64.0)]    // N=4: 2^6 = 64 = a_{-5}
    [InlineData(5, 128.0)]   // N=5: 2^7 = 128 = a_{-6}
    [InlineData(6, 256.0)]   // N=6: 2^8 = 256 = a_{-7}
    public void HPartPrefactor_EqualsTwoToTheNPlusTwo(int N, double expected)
    {
        var f = Build();
        Assert.Equal(expected, f.HPartPrefactor(N), precision: 12);
    }

    [Theory]
    [InlineData(2, -3)]    // N=2: ladder index -(2+1) = -3
    [InlineData(3, -4)]
    [InlineData(4, -5)]
    [InlineData(6, -7)]
    public void HPartLadderIndex_EqualsMinusNPlusOne(int N, int expected)
    {
        Assert.Equal(expected, Build().HPartLadderIndex(N));
    }

    [Theory]
    [InlineData(2, 4.0)]      // N=2: 4^1 = 4 = a_{-1} = d² for 1 qubit
    [InlineData(3, 16.0)]     // N=3: 4^2 = 16 = a_{-3} = d² for 2 qubits
    [InlineData(4, 64.0)]     // N=4: 4^3 = 64 = a_{-5} = d² for 3 qubits
    [InlineData(5, 256.0)]    // N=5: 4^4 = 256 = a_{-7} = d² for 4 qubits
    [InlineData(6, 1024.0)]   // N=6: 4^5 = 1024 = a_{-9} = d² for 5 qubits
    public void T1PartPrefactor_EqualsFourToTheNMinusOne(int N, double expected)
    {
        var f = Build();
        Assert.Equal(expected, f.T1PartPrefactor(N), precision: 12);
    }

    [Theory]
    [InlineData(2, -1)]    // N=2: ladder index 3-2*2 = -1
    [InlineData(3, -3)]
    [InlineData(4, -5)]
    [InlineData(6, -9)]
    public void T1PartLadderIndex_EqualsThreeMinusTwoN(int N, int expected)
    {
        Assert.Equal(expected, Build().T1PartLadderIndex(N));
    }

    [Theory]
    [InlineData(2, 1)]    // N_chain=2 → qubit count 1
    [InlineData(3, 2)]
    [InlineData(4, 3)]
    [InlineData(6, 5)]
    public void T1PartOperatorSpaceQubitCount_EqualsNMinusOne(int N, int expected)
    {
        Assert.Equal(expected, Build().T1PartOperatorSpaceQubitCount(N));
    }

    [Fact]
    public void HPartPrefactor_NLessThanTwo_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => Build().HPartPrefactor(1));
    }

    [Fact]
    public void T1PartPrefactor_NLessThanTwo_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => Build().T1PartPrefactor(1));
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F1T1AmplitudeDampingPi2Inheritance(null!, new Pi2OperatorSpaceMirrorClaim()));
        Assert.Throws<ArgumentNullException>(() =>
            new F1T1AmplitudeDampingPi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }

    [Fact]
    public void Anchor_References_F82_AndF1OpenQuestions_AndPi2Foundation()
    {
        var f = Build();
        Assert.Contains("ANALYTICAL_FORMULAS.md", f.Anchor);
        Assert.Contains("F1OpenQuestions.cs", f.Anchor);
        Assert.Contains("F1PalindromeIdentity.cs", f.Anchor);
        Assert.Contains("Pi2DyadicLadderClaim.cs", f.Anchor);
        Assert.Contains("Pi2OperatorSpaceMirrorClaim.cs", f.Anchor);
    }

    [Fact]
    public void Reconnaissance_EmitsFactorTable()
    {
        var f = Build();
        _out.WriteLine("");
        _out.WriteLine("    F1 T1-amplitude-damping closed form:");
        _out.WriteLine("    ‖M‖² = 2^(N+2)·n_YZ·‖H‖²_F  +  4^(N−1)·[3·Σγ² + 4·(Σγ)²]");
        _out.WriteLine("");
        _out.WriteLine("     N | 2^(N+2) | ladder | 4^(N−1) | ladder | qubit-count");
        _out.WriteLine("    ---|---------|--------|---------|--------|------------");
        for (int N = 2; N <= 6; N++)
        {
            _out.WriteLine($"     {N} | {f.HPartPrefactor(N),7:F0} | a_{f.HPartLadderIndex(N),-5} | " +
                          $"{f.T1PartPrefactor(N),7:F0} | a_{f.T1PartLadderIndex(N),-5} | {f.T1PartOperatorSpaceQubitCount(N)}");
        }
        _out.WriteLine("");
        _out.WriteLine($"    constants: ThreeMultiplier = {F1T1AmplitudeDampingPi2Inheritance.ThreeMultiplier} (T1 algebra)");
        _out.WriteLine($"               FourMultiplier  = {f.FourMultiplier} = a_{{-1}}");
    }
}
