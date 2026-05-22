using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F80PiCommutatorAnticommutatorIdentityTests
{
    private static F80PiCommutatorAnticommutatorIdentity Build() =>
        new F80PiCommutatorAnticommutatorIdentity(new RCPsiSquared.Core.F1.F1PalindromeIdentity());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Theory]
    [InlineData('X', 1)]
    [InlineData('Y', -1)]
    [InlineData('Z', 1)]
    public void Epsilon_MatchesPerSiteSign(char letter, int expected)
    {
        // F80 Step 5 per-site identities: μ(X·a)=+c_X(a)·X·μ(a), μ(Y·a)=−c_Y(a)·Y·μ(a),
        // μ(Z·a)=+c_Z(a)·Z·μ(a). The leading sign ε_P is +1 for X and Z, −1 for Y.
        Assert.Equal(expected, Build().Epsilon(letter));
    }

    [Fact]
    public void Epsilon_InvalidLetter_Throws()
    {
        Assert.Throws<ArgumentException>(() => Build().Epsilon('I'));
    }

    [Theory]
    [InlineData('X', 'Y', 1)]
    [InlineData('Y', 'X', 1)]
    [InlineData('X', 'Z', -1)]
    [InlineData('Z', 'X', -1)]
    public void SignFor_MatchesProvenStepFiveTable(char p, char q, int expected)
    {
        // Π·[H,·]·Π⁻¹ = s·{H,·}, s = −ε_P·ε_Q. (X,Y),(Y,X): s=+1 (M = −2i·H⊗I_bra);
        // (X,Z),(Z,X): s=−1 (M = +2i·I_ket⊗Hᵀ). Verified bit-exact N=3,4,5.
        Assert.Equal(expected, Build().SignFor(p, q));
    }

    [Theory]
    [InlineData('X', 'Y')]
    [InlineData('X', 'Z')]
    [InlineData('Y', 'X')]
    [InlineData('Z', 'X')]
    public void SignFor_EqualsMinusEpsilonProduct(char p, char q)
    {
        // The proof's closed form for the sign: s = −ε_P·ε_Q.
        var f = Build();
        Assert.Equal(-f.Epsilon(p) * f.Epsilon(q), f.SignFor(p, q));
    }

    [Theory]
    [InlineData('X', 'X')]   // no Y/Z partner
    [InlineData('Y', 'Z')]   // no X
    [InlineData('Z', 'Y')]   // no X
    public void SignFor_NonPi2OddBondPair_Throws(char p, char q)
    {
        // F80 Step 5 covers exactly the four Π²-odd 2-body pairs (one X, one Y/Z).
        Assert.Throws<ArgumentException>(() => Build().SignFor(p, q));
    }

    [Fact]
    public void Constructor_RejectsNullParent()
    {
        Assert.Throws<ArgumentNullException>(
            () => new F80PiCommutatorAnticommutatorIdentity(null!));
    }

    [Fact]
    public void Anchor_ReferencesF80ProofAndF1Parent()
    {
        var f = Build();
        Assert.Contains("PROOF_F80_BLOCH_SIGNWALK.md", f.Anchor);
        Assert.Contains("F1PalindromeIdentity.cs", f.Anchor);
    }
}
