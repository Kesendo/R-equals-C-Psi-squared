using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Tests.ChainSystems;

public class ChainSystemBuildersTests
{
    [Fact]
    public void BuildHamiltonian_XY_IsHermitian()
    {
        var chain = new ChainSystem(N: 3, J: 1.0, GammaZero: 0.05);
        var H = chain.BuildHamiltonian();
        var diff = H - H.ConjugateTranspose();
        Assert.True(diff.FrobeniusNorm() < 1e-12);
    }

    [Fact]
    public void BuildHamiltonian_Heisenberg_DiffersFromXY()
    {
        var xy = new ChainSystem(3, 1.0, 0.05).BuildHamiltonian();
        var heis = new ChainSystem(3, 1.0, 0.05, HType: HamiltonianType.Heisenberg).BuildHamiltonian();
        Assert.True((xy - heis).FrobeniusNorm() > 1e-3);
    }

    [Fact]
    public void BuildLiouvillian_PassesF1PalindromeIdentity()
    {
        // The whole Foundation in one test: ChainSystem.BuildLiouvillian() satisfies
        // Π·L·Π⁻¹ + L + 2σ·I = 0 bit-exactly.
        var chain = new ChainSystem(N: 2, J: 1.0, GammaZero: 0.05);
        var L = chain.BuildLiouvillian();
        var residual = PalindromeResidual.Build(L, chain.N, chain.SigmaGamma, PauliLetter.Z);
        Assert.True(residual.FrobeniusNorm() < 1e-9);
    }

    [Fact]
    public void SigmaGamma_IsNGammaZero()
    {
        var chain = new ChainSystem(N: 5, J: 1.0, GammaZero: 0.07);
        Assert.Equal(5 * 0.07, chain.SigmaGamma, 12);
    }
}
