using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.SymmetryFamily;
using RCPsiSquared.Core.Tests.TestHelpers;
using Xunit;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.SymmetryFamily;

public sealed class F93DetuningAntiPalindromicSpectralInvarianceTests
{
    [Fact]
    public void Tier_IsTier1Derived()
    {
        var sectors = new JointPopcountSectors();
        var f71 = new F71MirrorBlockRefinement(sectors);
        var inventory = new SymmetryFamilyInventory();
        var f93 = new F93DetuningAntiPalindromicSpectralInvariance(sectors, f71, inventory);
        Assert.Equal(Tier.Tier1Derived, f93.Tier);
    }

    [Theory]
    [InlineData(new double[] { 0.0, 0.0, 0.0 }, 0.0)]
    [InlineData(new double[] { 0.1, 0.2, 0.3, 0.4, 0.5 }, 0.0)]
    [InlineData(new double[] { 0.5, 0.5, 0.5, 0.5 }, 0.0)]
    public void AntiPalindromicDeviation_OfAntiPalindromicProfiles_IsZero(double[] hPerSite, double expected)
    {
        var dev = F93DetuningAntiPalindromicSpectralInvariance.AntiPalindromicDeviation(hPerSite);
        Assert.True(System.Math.Abs(dev - expected) < 1e-7,
            $"Expected dev={expected}, got {dev}");
    }

    [Fact]
    public void AntiPalindromicDeviation_OfPermutedProfile_IsNonzero()
    {
        // Profile [0.7, 0.2, 0.5, 0.3]: pairs (0.7+0.3, 0.2+0.5) = (1.0, 0.7), NOT constant.
        var hPerSite = new double[] { 0.7, 0.2, 0.5, 0.3 };
        var dev = F93DetuningAntiPalindromicSpectralInvariance.AntiPalindromicDeviation(hPerSite);
        Assert.True(dev > 0.05, $"Expected dev > 0.05 for non-anti-pal profile, got {dev}");
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    public void Spectrum_InvariantUnderAntiPalindromicH_ChainXY_UniformGamma(int N)
    {
        const double J = 1.0;
        const double gamma = 0.5;
        const double havg = 0.4;
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();

        // Reference: uniform h
        var Hunif = BuildXYChainPlusZDetuning(N, J, Enumerable.Repeat(havg, N).ToArray());
        var spectrumUnif = F71MirrorBlockRefinement.ComputeSpectrumPerBlock(Hunif, gammaPerSite, N);

        // Anti-palindromic h: h_l + h_{N-1-l} = 2·havg ∀l
        var hAnti = BuildAntiPalindromicH(N, havg);
        var Hanti = BuildXYChainPlusZDetuning(N, J, hAnti);
        var spectrumAnti = F71MirrorBlockRefinement.ComputeSpectrumPerBlock(Hanti, gammaPerSite, N);

        Assert.Equal(spectrumUnif.Length, spectrumAnti.Length);
        MultisetAssert.NearestNeighbourEqual(spectrumUnif, spectrumAnti, tolerance: 1e-9);
    }

    private static double[] BuildAntiPalindromicH(int N, double havg)
    {
        var h = new double[N];
        double mid = (N - 1) / 2.0;
        for (int l = 0; l < N; l++)
            h[l] = havg + 0.2 * (l - mid);
        return h;
    }

    private static ComplexMatrix BuildXYChainPlusZDetuning(int N, double J, double[] hPerSite)
    {
        if (hPerSite.Length != N)
            throw new System.ArgumentException($"hPerSite length {hPerSite.Length} != N = {N}");
        var allTerms = new List<PauliTerm>();
        Complex c = J / 2.0;
        for (int b = 0; b < N - 1; b++)
        {
            allTerms.Add(PauliTerm.TwoSite(N, b, PauliLetter.X, b + 1, PauliLetter.X, c));
            allTerms.Add(PauliTerm.TwoSite(N, b, PauliLetter.Y, b + 1, PauliLetter.Y, c));
        }
        for (int l = 0; l < N; l++)
        {
            if (System.Math.Abs(hPerSite[l]) > 1e-15)
                allTerms.Add(PauliTerm.SingleSite(N, l, PauliLetter.Z, hPerSite[l]));
        }
        return new PauliHamiltonian(N, allTerms).ToMatrix();
    }
}
