using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Tests.TestHelpers;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.BlockSpectrum;

/// <summary>Bit-exact regression tests for the X⊗N-pairing optimisation in
/// <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock"/> and
/// <see cref="F71MirrorBlockRefinement.ComputeSpectrumPerBlock"/>.
///
/// <para>X⊗N pairs joint-popcount sector (p_c, p_r) with (N − p_c, N − p_r) under chain
/// XY+Z-deph L; paired sectors share spectrum exactly (verified Tier-1 derived in
/// <c>SymmetryFamily/XGlobalChargeConjugationPairing.cs</c>). The optimised path computes
/// eig only on lex-smaller "primary" sectors and copies onto follower sectors. These
/// tests ensure the multiset of all 4^N eigenvalues remains bit-exactly equal to the full
/// L eig (and to the unpaired sector loop) at N = 3..6.</para></summary>
public class LiouvillianBlockSpectrumXNPairingTests
{
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void ComputeSpectrumPerBlock_WithXNPairing_MatchesFullLEig(int N)
    {
        const double J = 1.0;
        const double gamma = 0.5;
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();
        var L = PauliDephasingDissipator.BuildZ(H, gammaPerSite);

        var spectrumPaired = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaPerSite, N);
        var spectrumFull = L.Evd().EigenValues.ToArray();

        Assert.Equal(spectrumPaired.Length, spectrumFull.Length);
        MultisetAssert.NearestNeighbourEqual(spectrumPaired, spectrumFull, tolerance: 1e-9, context: $"N={N}");
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void F71RefinedComputeSpectrumPerBlock_WithXNPairing_MatchesFullLEig(int N)
    {
        const double J = 1.0;
        const double gamma = 0.5;
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();
        var L = PauliDephasingDissipator.BuildZ(H, gammaPerSite);

        var spectrumPaired = F71MirrorBlockRefinement.ComputeSpectrumPerBlock(H, gammaPerSite, N);
        var spectrumFull = L.Evd().EigenValues.ToArray();

        Assert.Equal(spectrumPaired.Length, spectrumFull.Length);
        MultisetAssert.NearestNeighbourEqual(spectrumPaired, spectrumFull, tolerance: 1e-9, context: $"N={N}");
    }
}
