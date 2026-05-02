using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F77;

namespace RCPsiSquared.Diagnostics.Tests.F77;

public class PauliPairTrichotomyTests
{
    private static ChainSystem MakeChain(int N) => new(N, J: 1.0, GammaZero: 0.05);

    [Fact]
    public void EmptyTerms_AreTruly()
    {
        Assert.Equal(TrichotomyClass.Truly,
            PauliPairTrichotomy.Classify(MakeChain(3), Array.Empty<PauliPairBondTerm>()));
    }

    [Fact]
    public void XXplusYY_IsTruly_UnderZDephasing()
    {
        // Canonical XY chain (XX + YY) is "truly" under Z-dephasing — the F1 anchor.
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.X),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Y),
        };
        Assert.Equal(TrichotomyClass.Truly,
            PauliPairTrichotomy.Classify(MakeChain(3), terms));
    }

    [Fact]
    public void Heisenberg_XX_YY_ZZ_IsTruly()
    {
        // Heisenberg = XX+YY+ZZ. ZZ is bit_b sum even, also "truly".
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.X),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Y),
            new PauliPairBondTerm(PauliLetter.Z, PauliLetter.Z),
        };
        Assert.Equal(TrichotomyClass.Truly,
            PauliPairTrichotomy.Classify(MakeChain(3), terms));
    }

    [Fact]
    public void YZplusZY_IsSoft_UnderZDephasing()
    {
        // YZ+ZY: bit_b sum is (1+1)+(1+1) = 4 → Π²-even, but bond-flipped form. Soft per F77.
        // This is the EQ-030 Marrakesh-confirmed soft Hamiltonian (drop=28).
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Z),
            new PauliPairBondTerm(PauliLetter.Z, PauliLetter.Y),
        };
        Assert.Equal(TrichotomyClass.Soft,
            PauliPairTrichotomy.Classify(MakeChain(3), terms));
    }

    [Fact]
    public void XXplusXY_IsHard()
    {
        // XX+XY: mixes Π²-even with Π²-odd. Hard per F77 trichotomy.
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.X),
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.Y),
        };
        Assert.Equal(TrichotomyClass.Hard,
            PauliPairTrichotomy.Classify(MakeChain(3), terms));
    }

    [Fact]
    public void XYplusYX_IsSoft_BondFlip()
    {
        // XY+YX: bond-flipped Z-free pair. Soft per Marrakesh skeleton-trace test.
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.Y),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.X),
        };
        Assert.Equal(TrichotomyClass.Soft,
            PauliPairTrichotomy.Classify(MakeChain(3), terms));
    }
}
