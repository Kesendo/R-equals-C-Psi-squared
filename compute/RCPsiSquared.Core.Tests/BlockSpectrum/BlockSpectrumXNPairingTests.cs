using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Tests.TestHelpers;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.BlockSpectrum;

/// <summary>Regression tests for the X⊗N copy
/// (<see cref="LiouvillianBlockSpectrum.SectorPairing.XNCopy"/>) in
/// <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock"/> and
/// <see cref="F71MirrorBlockRefinement.ComputeSpectrumPerBlock"/>.
///
/// <para>X⊗N pairs joint-popcount sector (p_c, p_r) with (N − p_c, N − p_r) under chain
/// XY+Z-deph L; paired sectors share spectrum exactly (verified Tier-1 derived in
/// <c>SymmetryFamily/XGlobalChargeConjugationPairing.cs</c>). The X⊗N copy computes eig
/// on one sector of each pair and copies it onto the other. These tests check the
/// multiset of all 4^N eigenvalues against the full L eig at N = 3..6 (nearest neighbour,
/// within 1e-9).</para></summary>
public class BlockSpectrumXNPairingTests
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

        var spectrumPaired = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(
            H, gammaPerSite, N, LiouvillianBlockSpectrum.SectorPairing.XNCopy);
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

        var spectrumPaired = F71MirrorBlockRefinement.ComputeSpectrumPerBlock(
            H, gammaPerSite, N, LiouvillianBlockSpectrum.SectorPairing.XNCopy);
        var spectrumFull = L.Evd().EigenValues.ToArray();

        Assert.Equal(spectrumPaired.Length, spectrumFull.Length);
        MultisetAssert.NearestNeighbourEqual(spectrumPaired, spectrumFull, tolerance: 1e-9, context: $"N={N}");
    }
}
