using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F87;

namespace RCPsiSquared.Diagnostics.Tests.F87;

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
        // Canonical XY chain (XX + YY) is "truly" under Z-dephasing, the F1 anchor.
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
        // YZ+ZY: bit_b sum is (1+1)+(1+1) = 4 → Π²-even, but bond-flipped form. Soft per F87.
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
        // XX+XY: mixes Π²-even with Π²-odd. Hard per F87 trichotomy.
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

    [Fact]
    public void Classify_VerdictIsGammaInvariant_HardAndSoftWitnesses()
    {
        // Windowed all-γ theorem (WindowedConverseAllGammaClaim, closed 2026-06-10): a Hard
        // verdict for a windowed diagonal-cell pair is γ-universal, hard at one γ is hard at
        // every γ > 0 (Pascal-Gram positivity F117, no residual). Truly/Soft are γ-independent
        // by construction. Sweep both witnesses across the γ axis; the verdict must not move.
        //
        // Hard witness: XXZ+XZX, a windowed diagonal-cell Klein-(0,1) k=3 pair (X/Y masks
        // 011 vs 101 have different (1+x)-valuations, hard per PalindromeHardSweepTests).
        // No k=2 diagonal-cell pair is hard (both diagonal-cell Mixed k=2 strings share one
        // mask), so the hard witness goes through the k-body overload at N=4.
        var hardTemplates = new[]
        {
            new PauliTerm(new[] { PauliLetter.X, PauliLetter.X, PauliLetter.Z }, Complex.One),
            new PauliTerm(new[] { PauliLetter.X, PauliLetter.Z, PauliLetter.X }, Complex.One),
        };
        // Soft witness: YZ+ZY, the EQ-030 Marrakesh-confirmed bond-flip pair from above.
        var softTerms = new[]
        {
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Z),
            new PauliPairBondTerm(PauliLetter.Z, PauliLetter.Y),
        };
        foreach (double gamma in new[] { 0.05, 0.3, 1.0 })
        {
            Assert.Equal(TrichotomyClass.Hard,
                PauliPairTrichotomy.Classify(new ChainSystem(N: 4, J: 1.0, GammaZero: gamma), hardTemplates));
            Assert.Equal(TrichotomyClass.Soft,
                PauliPairTrichotomy.Classify(new ChainSystem(N: 3, J: 1.0, GammaZero: gamma), softTerms));
        }
    }
}
