using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.SymmetryFamily;
using Xunit;

namespace RCPsiSquared.Core.Tests.SymmetryFamily;

public sealed class XGlobalChargeConjugationPairingTests
{
    [Fact]
    public void Tier_IsTier1Derived()
    {
        var x = new XGlobalChargeConjugationPairing(new SymmetryFamilyInventory());
        Assert.Equal(Tier.Tier1Derived, x.Tier);
    }

    [Theory]
    [InlineData(3, 0, 0, 3, 3)]
    [InlineData(3, 1, 0, 2, 3)]
    [InlineData(4, 1, 2, 3, 2)]
    [InlineData(5, 2, 3, 3, 2)]
    public void PairSector_MapsCorrectly(int N, int pCol, int pRow, int expectedPairCol, int expectedPairRow)
    {
        var (pcOut, prOut) = XGlobalChargeConjugationPairing.PairSector(N, pCol, pRow);
        Assert.Equal(expectedPairCol, pcOut);
        Assert.Equal(expectedPairRow, prOut);
    }

    [Theory]
    [InlineData(4, 2, 2, true)]
    [InlineData(3, 1, 2, false)]
    [InlineData(4, 1, 3, false)]
    public void IsSelfPaired_DetectsCorrectly(int N, int pCol, int pRow, bool expected)
    {
        Assert.Equal(expected, XGlobalChargeConjugationPairing.IsSelfPaired(N, pCol, pRow));
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void PairedSectors_ShareSpectrum_ChainXY_UniformGamma(int N)
    {
        const double J = 1.0;
        const double gamma = 0.5;
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();
        var spectrum = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaPerSite, N);

        var decomp = JointPopcountSectorBuilder.Build(N);
        var bySector = new Dictionary<(int, int), List<Complex>>();
        int offset = 0;
        foreach (var sector in decomp.SectorRanges)
        {
            var list = new List<Complex>(sector.Size);
            for (int i = 0; i < sector.Size; i++) list.Add(spectrum[offset + i]);
            bySector[(sector.PCol, sector.PRow)] = list;
            offset += sector.Size;
        }

        foreach (var ((pc, pr), eigs) in bySector)
        {
            var (pairPc, pairPr) = XGlobalChargeConjugationPairing.PairSector(N, pc, pr);
            if ((pc, pr) == (pairPc, pairPr)) continue;  // self-paired, trivial
            var pairEigs = bySector[(pairPc, pairPr)];
            Assert.Equal(eigs.Count, pairEigs.Count);
            AssertMultisetEqual(eigs, pairEigs, tolerance: 1e-9);
        }
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
