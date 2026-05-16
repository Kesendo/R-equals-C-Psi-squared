using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Diagnostics.F87;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>Bridge claim verification: F89's bond Hamiltonian H_b = J·(XX+YY) is
/// F87-Truly across every chain N, satisfying the F89 AT-lock precondition.</summary>
public class F89F87TrulyInheritanceTests
{
    private static ChainSystem MakeChain(int N) => new(N, J: 1.0, GammaZero: 0.05);

    // F87 trichotomy is N-invariant for Pauli-pair bond Hamiltonians (the per-bond
    // Π² eigenvalue is bit-pattern arithmetic, not L-spectral); we sample N=3..5 to
    // keep the full-L Evd cheap (each test builds a 4^N × 4^N complex Liouvillian).
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void ClassifyF89Hamiltonian_AcrossN_IsTrulyEveryTime(int N)
    {
        Assert.Equal(TrichotomyClass.Truly,
            F89F87TrulyInheritance.ClassifyF89Hamiltonian(MakeChain(N)));
    }

    [Theory]
    [InlineData(3)]
    [InlineData(5)]
    public void IsAtLockPreconditionSatisfied_AtAnyN_True(int N)
    {
        Assert.True(F89F87TrulyInheritance.IsAtLockPreconditionSatisfied(MakeChain(N)));
    }

    [Fact]
    public void F89BondTerms_ExpectedShape_XXPlusYY()
    {
        var terms = F89F87TrulyInheritance.F89BondTerms;
        Assert.Equal(2, terms.Count);
        Assert.Equal((Core.Pauli.PauliLetter.X, Core.Pauli.PauliLetter.X), (terms[0].LetterA, terms[0].LetterB));
        Assert.Equal((Core.Pauli.PauliLetter.Y, Core.Pauli.PauliLetter.Y), (terms[1].LetterA, terms[1].LetterB));
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var trichotomy = new F87TrichotomyClassification();
        var bridge = new F89F87TrulyInheritance(trichotomy);
        Assert.Equal(Tier.Tier1Derived, bridge.Tier);
    }

    [Fact]
    public void Constructor_NullTrichotomy_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => new F89F87TrulyInheritance(null!));
    }
}
