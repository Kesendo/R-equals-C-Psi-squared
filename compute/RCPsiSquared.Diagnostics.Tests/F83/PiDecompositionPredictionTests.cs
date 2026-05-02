using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F49;
using RCPsiSquared.Diagnostics.F83;

namespace RCPsiSquared.Diagnostics.Tests.F83;

public class PiDecompositionPredictionTests
{
    [Fact]
    public void TrulyHamiltonian_Predicts_MAllZero_AntiFraction_Zero()
    {
        // XX+YY truly → all terms drop → ‖M‖² = 0, anti-fraction defaults to 0.
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.X),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Y),
        };
        var d = PiDecompositionPrediction.Predict(chain, terms);
        Assert.Equal(0.0, d.MSquared, 12);
        Assert.Equal(0.0, d.MAntiSquared, 12);
        Assert.Equal(0.0, d.AntiFraction, 12);
    }

    [Fact]
    public void PurePi2OddHamiltonian_HasAntiFractionOneHalf()
    {
        // F81 50/50 finding: pure Π²-odd → anti-fraction = 1/2.
        // XY+YX is pure Π²-odd (each has bit_b sum = 0+1 = 1).
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.Y),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.X),
        };
        var d = PiDecompositionPrediction.Predict(chain, terms);
        Assert.Equal(0.5, d.AntiFraction, 10);
    }

    [Fact]
    public void PurePi2EvenNonTrulyHamiltonian_HasAntiFractionZero()
    {
        // F81 100/0 finding: pure Π²-even non-truly → anti-fraction = 0 (M_anti = 0).
        // YZ+ZY: bit_b sum (1+1)+(1+1)=4 → both terms Π²-even, non-truly.
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Z),
            new PauliPairBondTerm(PauliLetter.Z, PauliLetter.Y),
        };
        var d = PiDecompositionPrediction.Predict(chain, terms);
        Assert.Equal(0.0, d.AntiFraction, 10);
    }

    [Fact]
    public void EqualMix_Pi2OddAndEven_GivesAntiFractionOneSixth()
    {
        // r = ‖H_even_nontruly‖²/‖H_odd‖² = 1 → anti-fraction = 1/(2+4·1) = 1/6.
        // Use XY+YX (Π²-odd) + YZ+ZY (Π²-even non-truly), both with same coefficient.
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.Y),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.X),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Z),
            new PauliPairBondTerm(PauliLetter.Z, PauliLetter.Y),
        };
        var d = PiDecompositionPrediction.Predict(chain, terms);
        Assert.Equal(1.0, d.R, 10);
        Assert.Equal(1.0 / 6.0, d.AntiFraction, 10);
    }

    [Fact]
    public void Decomposition_Sums_Anti_PlusSym_EqualsTotal()
    {
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Z),
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.Y),
        };
        var d = PiDecompositionPrediction.Predict(chain, terms);
        Assert.Equal(d.MSquared, d.MAntiSquared + d.MSymSquared, 10);
    }
}
