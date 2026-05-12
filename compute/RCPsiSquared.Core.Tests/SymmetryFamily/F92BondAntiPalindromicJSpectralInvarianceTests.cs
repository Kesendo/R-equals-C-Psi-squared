using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.SymmetryFamily;
using Xunit;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.SymmetryFamily;

public sealed class F92BondAntiPalindromicJSpectralInvarianceTests
{
    [Fact]
    public void Tier_IsTier1Derived()
    {
        var sectors = new JointPopcountSectors();
        var f71 = new F71MirrorBlockRefinement(sectors);
        var f92 = new F92BondAntiPalindromicJSpectralInvariance(sectors, f71);
        Assert.Equal(Tier.Tier1Derived, f92.Tier);
    }

    [Theory]
    [InlineData(new double[] { 0.5, 0.5, 0.5 }, 0.0)]
    [InlineData(new double[] { 0.3, 0.4, 0.5, 0.6, 0.7 }, 0.0)]
    [InlineData(new double[] { 0.3, 0.5, 0.4, 0.6 }, 0.0)]
    public void AntiPalindromicDeviation_OfAntiPalindromicProfiles_IsZero(double[] bondJ, double expected)
    {
        var dev = F92BondAntiPalindromicJSpectralInvariance.AntiPalindromicDeviation(bondJ);
        Assert.True(System.Math.Abs(dev - expected) < 1e-7,
            $"Expected dev={expected}, got {dev}");
    }

    [Fact]
    public void AntiPalindromicDeviation_OfPermutedProfile_IsNonzero()
    {
        // Profile [0.7, 0.2, 0.5, 0.3, 0.6, 0.4]: same multiset as monotonic but pairs are
        // (0.7+0.4), (0.2+0.6), (0.5+0.3) = (1.1, 0.8, 0.8), NOT constant. Avg = 0.45.
        var bondJ = new double[] { 0.7, 0.2, 0.5, 0.3, 0.6, 0.4 };
        var dev = F92BondAntiPalindromicJSpectralInvariance.AntiPalindromicDeviation(bondJ);
        Assert.True(dev > 0.05, $"Expected dev > 0.05 for non-anti-pal profile, got {dev}");
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    public void Spectrum_InvariantUnderAntiPalindromicJ_ChainXY_UniformGamma(int N)
    {
        // F92 (parallel to F91): the F71-REFINED DIAGONAL-BLOCK spectrum is invariant under
        // any J satisfying J_b + J_{N-2-b} = 2·J_avg. The full L spectrum
        // (LiouvillianBlockSpectrum.ComputeSpectrumPerBlock) generally DIFFERS; anti-palindromy
        // breaks F71 as an L-symmetry, populating the (F71-even ↔ F71-odd) cross blocks.
        const double gamma = 0.5;
        const double Javg = 1.0;
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();

        // Reference: uniform J = Javg
        var Hunif = BuildXYChainWithBondJ(N, Enumerable.Repeat(Javg, N - 1).ToArray());
        var spectrumUnif = F71MirrorBlockRefinement.ComputeSpectrumPerBlock(Hunif, gammaPerSite, N);

        // Anti-palindromic J: J_b + J_{N-2-b} = 2·Javg
        var Janti = BuildAntiPalindromicJ(N, Javg);
        var Hanti = BuildXYChainWithBondJ(N, Janti);
        var spectrumAnti = F71MirrorBlockRefinement.ComputeSpectrumPerBlock(Hanti, gammaPerSite, N);

        Assert.Equal(spectrumUnif.Length, spectrumAnti.Length);
        AssertMultisetEqual(spectrumUnif, spectrumAnti, tolerance: 1e-9);
    }

    private static double[] BuildAntiPalindromicJ(int N, double Javg)
    {
        // For N qubits there are N-1 bonds, b ∈ {0..N-2}. Anti-palindromic: J_b + J_{N-2-b} = 2·Javg.
        // Linear gradient: J_b = Javg + 0.3·(b - (N-2)/2).
        int numBonds = N - 1;
        var J = new double[numBonds];
        double mid = (numBonds - 1) / 2.0;
        for (int b = 0; b < numBonds; b++) J[b] = Javg + 0.3 * (b - mid);
        return J;
    }

    private static ComplexMatrix BuildXYChainWithBondJ(int N, double[] bondJ)
    {
        if (bondJ.Length != N - 1)
            throw new System.ArgumentException($"bondJ length {bondJ.Length} != N-1 = {N - 1}");
        var allTerms = new List<PauliTerm>();
        for (int b = 0; b < N - 1; b++)
        {
            Complex c = bondJ[b] / 2.0;  // F86 convention (J/2) per existing PauliHamiltonian.XYChain
            allTerms.Add(PauliTerm.TwoSite(N, b, PauliLetter.X, b + 1, PauliLetter.X, c));
            allTerms.Add(PauliTerm.TwoSite(N, b, PauliLetter.Y, b + 1, PauliLetter.Y, c));
        }
        return new PauliHamiltonian(N, allTerms).ToMatrix();
    }

    private static void AssertMultisetEqual(IList<Complex> a, IList<Complex> b, double tolerance)
    {
        Assert.Equal(a.Count, b.Count);
        var unmatched = b.ToList();
        foreach (var x in a)
        {
            int bestIdx = -1; double bestDist = double.MaxValue;
            for (int i = 0; i < unmatched.Count; i++)
            {
                double d = (x - unmatched[i]).Magnitude;
                if (d < bestDist) { bestDist = d; bestIdx = i; }
            }
            Assert.True(bestDist < tolerance,
                $"No multiset match within {tolerance}: {x} (best {unmatched[bestIdx]} at {bestDist:E3})");
            unmatched.RemoveAt(bestIdx);
        }
    }
}
