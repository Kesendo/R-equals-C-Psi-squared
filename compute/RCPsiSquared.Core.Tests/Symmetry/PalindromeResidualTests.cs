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
    public void F1_Palindrome_IsExactlyZero_ForXDephasedHeisenbergChain()
    {
        // Heisenberg = XX+YY+ZZ is invariant under cyclic Pauli rotation X→Y→Z, so it's
        // "truly" under X-, Y-, and Z-dephasing. F1 must hold bit-exactly with the
        // X-axis-specific Π built by PiOperator.BuildFull(N, PauliLetter.X).
        int N = 3;
        double gamma = 0.05;
        var H = PauliHamiltonian.HeisenbergChain(N, J: 1.0).ToMatrix();
        var gammaList = Enumerable.Repeat(gamma, N).ToArray();
        var L = PauliDephasingDissipator.Build(H, gammaList, PauliLetter.X);
        double sigmaGamma = N * gamma;

        var residual = PalindromeResidual.Build(L, N, sigmaGamma, PauliLetter.X);
        double frob = residual.FrobeniusNorm();
        Assert.True(frob < 1e-9, $"F1 X-dephased Heisenberg residual not zero at N={N}: ‖M‖_F = {frob}");
    }

    [Fact]
    public void F1_Palindrome_IsExactlyZero_ForYDephasedHeisenbergChain()
    {
        // Same SU(2)-cyclic argument as the X-dephasing case: Heisenberg is truly under Y too.
        int N = 3;
        double gamma = 0.05;
        var H = PauliHamiltonian.HeisenbergChain(N, J: 1.0).ToMatrix();
        var gammaList = Enumerable.Repeat(gamma, N).ToArray();
        var L = PauliDephasingDissipator.Build(H, gammaList, PauliLetter.Y);
        double sigmaGamma = N * gamma;

        var residual = PalindromeResidual.Build(L, N, sigmaGamma, PauliLetter.Y);
        double frob = residual.FrobeniusNorm();
        Assert.True(frob < 1e-9, $"F1 Y-dephased Heisenberg residual not zero at N={N}: ‖M‖_F = {frob}");
    }

    [Fact]
    public void F1_Palindrome_IsExactlyZero_ForNonUniformZDephasing()
    {
        // F1 is per-site additive in the Klein parities, so non-uniform γ_l still gives M=0
        // for a truly H. The σ in Π·L·Π⁻¹ + L + 2σ·I = 0 must be the actual Σγ_l, not the
        // chain.SigmaGamma uniform shortcut. Documents the caller pattern when γ is non-uniform.
        int N = 3;
        double[] gammaList = { 0.04, 0.07, 0.06 };  // non-uniform
        var H = PauliHamiltonian.XYChain(N, J: 1.0).ToMatrix();
        var L = PauliDephasingDissipator.BuildZ(H, gammaList);
        double sigmaGamma = gammaList.Sum();  // 0.17, NOT chain.SigmaGamma (which is N·γ₀)

        var residual = PalindromeResidual.Build(L, N, sigmaGamma, PauliLetter.Z);
        double frob = residual.FrobeniusNorm();
        Assert.True(frob < 1e-9, $"F1 non-uniform γ residual not zero at N={N}: ‖M‖_F = {frob}");
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
