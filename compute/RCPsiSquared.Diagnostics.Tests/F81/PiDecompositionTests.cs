using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F49;
using RCPsiSquared.Diagnostics.F81;

namespace RCPsiSquared.Diagnostics.Tests.F81;

public class PiDecompositionTests
{
    [Fact]
    public void TrulyHamiltonian_HasZero_M_And_F81_HoldsTrivially()
    {
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.X),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Y),
        };
        var d = PiDecomposition.Decompose(chain, terms);
        Assert.True(d.MNormSquared < 1e-18);
        Assert.True(d.F81Violation < 1e-9);
    }

    [Fact]
    public void PurePi2OddHamiltonian_F81_HoldsBitExactly()
    {
        // F81 Tier 1: for pure Π²-odd H (no T1), M_anti = L_{H_odd} bit-exactly.
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.Y),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.X),
        };
        var d = PiDecomposition.Decompose(chain, terms);
        Assert.True(d.F81Violation < 1e-9, $"F81 violation = {d.F81Violation:E3}");
    }

    [Fact]
    public void PurePi2EvenNonTruly_M_anti_IsZero_LHOdd_IsZero()
    {
        // YZ+ZY is Π²-even non-truly. H_odd = 0 (no Π²-odd terms).
        // F81 says M_anti = L_{H_odd} = 0 → M_anti = 0.
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Z),
            new PauliPairBondTerm(PauliLetter.Z, PauliLetter.Y),
        };
        var d = PiDecomposition.Decompose(chain, terms);
        Assert.True(d.LHOddNormSquared < 1e-18);
        Assert.True(d.MAntiNormSquared < 1e-9, $"M_anti norm² = {d.MAntiNormSquared:E3}");
        Assert.True(d.F81Violation < 1e-9);
    }

    [Fact]
    public void Decomposition_MSymPlusMAnti_EqualsM()
    {
        // M = M_sym + M_anti by construction.
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.Y),
        };
        var d = PiDecomposition.Decompose(chain, terms);
        var sum = d.MSym + d.MAnti;
        var diff = sum - d.M;
        Assert.True(diff.FrobeniusNorm() < 1e-12);
    }

    [Fact]
    public void Decomposition_NormSquared_Decomposes_Orthogonally()
    {
        // Frobenius orthogonality: ‖M‖² = ‖M_sym‖² + ‖M_anti‖².
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.Y),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Z),
        };
        var d = PiDecomposition.Decompose(chain, terms);
        Assert.Equal(d.MNormSquared, d.MSymNormSquared + d.MAntiNormSquared, 8);
    }

    [Fact]
    public void T1_BreaksF81_ViolationIsPositive()
    {
        // T1 amplitude damping breaks F81; the violation is the quantitative diagnostic.
        var chain = new ChainSystem(N: 2, J: 1.0, GammaZero: 0.05);
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.X),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Y),
        };
        double[] gammaT1 = { 0.02, 0.02 };
        var d = PiDecomposition.Decompose(chain, terms, gammaT1);
        Assert.True(d.F81Violation > 1e-3, $"T1 should produce F81 violation > 0; got {d.F81Violation}");
    }
}
