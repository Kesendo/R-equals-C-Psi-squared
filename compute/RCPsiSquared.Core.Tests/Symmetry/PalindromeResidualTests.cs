using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class PalindromeResidualTests
{
    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    public void F1_Palindrome_IsExactlyZero_ForZDephasedXYChain(int N)
    {
        // The canonical F1 anchor: Π · L · Π⁻¹ + L + 2σ·I = 0 for Z-dephased XY chain.
        // This is the structural identity behind the entire chiral-classification stack.
        double gamma = 0.05;
        var H = PauliHamiltonian.XYChain(N, J: 1.0).ToMatrix();
        var gammaList = Enumerable.Repeat(gamma, N).ToArray();
        var L = PauliDephasingDissipator.BuildZ(H, gammaList);
        double sigmaGamma = N * gamma;

        var residual = PalindromeResidual.Build(L, N, sigmaGamma, PauliLetter.Z);
        double frob = residual.FrobeniusNorm();
        Assert.True(frob < 1e-9, $"F1 palindrome residual not zero at N={N}: ‖M‖_F = {frob}");
    }

    [Fact]
    public void F1_Palindrome_IsExactlyZero_ForZDephasedHeisenbergChain()
    {
        // Heisenberg = XY + ZZ. ZZ is bit_b = (1+1)=2 even, hence Π²-even, also "truly" → palindrome holds.
        int N = 3;
        double gamma = 0.05;
        var H = PauliHamiltonian.HeisenbergChain(N, J: 1.0).ToMatrix();
        var gammaList = Enumerable.Repeat(gamma, N).ToArray();
        var L = PauliDephasingDissipator.BuildZ(H, gammaList);
        double sigmaGamma = N * gamma;

        var residual = PalindromeResidual.Build(L, N, sigmaGamma, PauliLetter.Z);
        Assert.True(residual.FrobeniusNorm() < 1e-9);
    }

    [Fact]
    public void F1_Palindrome_BreaksFor_T1Dissipator()
    {
        // T1 amplitude damping σ⁻ = (X − iY)/2 carries Y-component (bit_b=1) and breaks F1 trivially.
        int N = 2;
        double gammaZ = 0.05;
        double gammaT1 = 0.02;
        var H = PauliHamiltonian.XYChain(N, J: 1.0).ToMatrix();
        var L = T1Dissipator.Build(H,
            Enumerable.Repeat(gammaZ, N).ToArray(),
            Enumerable.Repeat(gammaT1, N).ToArray());
        // The σ_offset for the T1 case includes only the Z-dephasing part (T1's contribution to σ
        // shows up as a separate residual term). Test against pure Z-dephasing σ: residual ≠ 0.
        var residual = PalindromeResidual.Build(L, N, N * gammaZ, PauliLetter.Z);
        Assert.True(residual.FrobeniusNorm() > 1e-3, "T1 dissipator should give nonzero F1 residual");
    }
}
