using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F49;

namespace RCPsiSquared.Diagnostics.Tests.F49;

public class FrobeniusScalingTests
{
    [Fact]
    public void TrulyHamiltonian_Predicts_ZeroResidualNorm()
    {
        // XX+YY (truly) → all terms drop in F85 sum → predicted ‖M‖² = 0.
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.X, PauliLetter.X),
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Y),
        };
        Assert.Equal(0.0, FrobeniusScaling.PredictNormSquaredFromTerms(chain, terms), 12);
    }

    [Fact]
    public void Pi2OddSoftHamiltonian_Predicts_NonzeroNorm_MatchesActualResidual()
    {
        // YZ+ZY soft Hamiltonian: bit_b sum (1+1)+(1+1)=4 → Π²-even? Wait that's both sites.
        // Per term: YZ has bit_b=(1,1) sum=2 even, ZY has bit_b=(1,1) sum=2 even. Π²-even.
        // But not truly (each has #Y=1 odd). So Π²-even non-truly → c=2.
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var terms = new[]
        {
            new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Z),
            new PauliPairBondTerm(PauliLetter.Z, PauliLetter.Y),
        };
        double predicted = FrobeniusScaling.PredictNormSquaredFromTerms(chain, terms);
        Assert.True(predicted > 0, "soft non-truly Hamiltonian should give nonzero ‖M‖²");

        // Cross-check: build L explicitly and measure ‖M‖_F², compare.
        var bilinearSpec = terms.Select(t =>
            (t.LetterA, t.LetterB, (System.Numerics.Complex)chain.J)).ToList();
        var H = PauliHamiltonian.Bilinear(chain.N, chain.Bonds, bilinearSpec).ToMatrix();
        var L = PauliDephasingDissipator.BuildZ(H, Enumerable.Repeat(chain.GammaZero, chain.N).ToArray());
        var M = PalindromeResidual.Build(L, chain.N, chain.SigmaGamma, PauliLetter.Z);
        double actual = Math.Pow(M.FrobeniusNorm(), 2);
        // Match within 5% (small numerical tolerance for Frobenius computation through the basis transform)
        Assert.InRange(predicted, 0.95 * actual, 1.05 * actual);
    }

    [Fact]
    public void T1Contribution_IsHamiltonianIndependent()
    {
        // T1 part: 4^(N-1) · [3·Σγ² + 4·(Σγ)²], no Hamiltonian dependence.
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        double[] gammaT1 = { 0.02, 0.03, 0.05 };
        double withTerms = FrobeniusScaling.PredictNormSquaredFromTerms(chain,
            new[] { new PauliPairBondTerm(PauliLetter.X, PauliLetter.X),
                    new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Y) },
            gammaT1);
        double withoutTerms = FrobeniusScaling.PredictNormSquaredFromTerms(chain,
            Array.Empty<PauliPairBondTerm>(), gammaT1);
        // Both should equal the T1-only contribution since H is "truly" → z_part = 0.
        double sumG = gammaT1.Sum();
        double sumG2 = gammaT1.Sum(g => g * g);
        double expected = Math.Pow(4, chain.N - 1) * (3 * sumG2 + 4 * sumG * sumG);
        Assert.Equal(expected, withTerms, 10);
        Assert.Equal(expected, withoutTerms, 10);
    }

    [Theory]
    [InlineData(2, HamiltonianClass.Main, 1.0)]
    [InlineData(3, HamiltonianClass.Main, 8.0)]
    [InlineData(4, HamiltonianClass.Main, 48.0)]
    public void PredictNormSquared_WithCh_MatchesScalingFormula(int N, HamiltonianClass cls, double expectedFactor)
    {
        var chain = new ChainSystem(N, J: 1.0, GammaZero: 0.05);
        Assert.Equal(expectedFactor, FrobeniusScaling.PredictNormSquared(chain, cH: 1.0, cls), 10);
    }
}
