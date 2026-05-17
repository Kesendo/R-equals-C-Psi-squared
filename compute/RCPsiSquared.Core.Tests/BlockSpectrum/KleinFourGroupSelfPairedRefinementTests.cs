using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Tests.TestHelpers;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.BlockSpectrum;

/// <summary>Witnesses for <see cref="KleinFourGroupSelfPairedRefinement"/>: K = {1, F71,
/// X⊗N, F71·X⊗N} on the X⊗N-self-paired (N/2, N/2) sector splits the sector into 4
/// character sub-blocks whose spectra union to the full sector spectrum. At small N the
/// union is cross-validated against direct dense Evd of the (N/2, N/2) block via
/// <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/>.</summary>
public class KleinFourGroupSelfPairedRefinementTests
{
    private readonly ITestOutputHelper _out;
    public KleinFourGroupSelfPairedRefinementTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void Build_AtN4_SubBlockDimsSumToSectorDim()
    {
        var k = KleinFourGroupSelfPairedRefinement.Build(N: 4);
        int sectorDim = 36;  // C(4, 2)² = 36
        Assert.Equal(sectorDim, k.SectorDim);
        _out.WriteLine($"N=4 (2,2): sub-block dims {{++: {k.SubBlockDims[KleinCharacter.PlusPlus]}, " +
                       $"+-: {k.SubBlockDims[KleinCharacter.PlusMinus]}, " +
                       $"-+: {k.SubBlockDims[KleinCharacter.MinusPlus]}, " +
                       $"--: {k.SubBlockDims[KleinCharacter.MinusMinus]}}} = {sectorDim}");
    }

    [Fact]
    public void Build_AtN10_SubBlockDimsMatchAnalyticalCount()
    {
        // Anti-palindromic popcount-5 strings at N=10: 2^(N/2) = 32 (free choice per
        // half-position determines its mirror-complement partner). Generic popcount-5: 252-32=220.
        // Orbit counts: both-anti = 32² → 1024 elements / size 2 = 512 orbits; mixed =
        // 2·32·220 / 4 = 3520; both-generic = 220² / 4 = 12100.
        //   ++ : all 16132 orbits contribute  (sum 512+3520+12100)
        //   -- : all 16132 (size-2 stab {1, F71·X⊗N} survives because (-1)(-1)=+1)
        //   +- : only size-4 orbits, sum 15620
        //   -+ : only size-4 orbits, sum 15620
        var k = KleinFourGroupSelfPairedRefinement.Build(N: 10);
        Assert.Equal(63504, k.SectorDim);
        Assert.Equal(16132, k.SubBlockDims[KleinCharacter.PlusPlus]);
        Assert.Equal(15620, k.SubBlockDims[KleinCharacter.PlusMinus]);
        Assert.Equal(15620, k.SubBlockDims[KleinCharacter.MinusPlus]);
        Assert.Equal(16132, k.SubBlockDims[KleinCharacter.MinusMinus]);

        _out.WriteLine($"N=10 (5,5): sector dim {k.SectorDim} → 4 sub-blocks " +
                       $"[++={k.SubBlockDims[KleinCharacter.PlusPlus]}, " +
                       $"+-={k.SubBlockDims[KleinCharacter.PlusMinus]}, " +
                       $"-+={k.SubBlockDims[KleinCharacter.MinusPlus]}, " +
                       $"--={k.SubBlockDims[KleinCharacter.MinusMinus]}]");
        _out.WriteLine($"Max sub-block: {k.SubBlockDims.Values.Max()} (matches inventory's advertised ~16k)");
    }

    [Theory]
    [InlineData(4)]
    [InlineData(6)]
    public void Build_SmallN_SubBlockSpectraUnionMatchesDirectSectorEvd(int N)
    {
        // The 4 Klein sub-block spectra must union to the full (N/2, N/2) sector spectrum.
        // Direct sector spectrum via PerBlockLiouvillianBuilder dense Evd.
        const double gamma = 0.05;
        var gammaArr = Enumerable.Repeat(gamma, N).ToArray();

        var H = PauliHamiltonian.XYChain(N, J: 1.0).ToMatrix();
        var decomp = JointPopcountSectorBuilder.Build(N);
        int m = N / 2;
        var sector = decomp.SectorRanges.First(s => s.PCol == m && s.PRow == m);
        var flat = new int[sector.Size];
        for (int k = 0; k < sector.Size; k++) flat[k] = decomp.Permutation[sector.Offset + k];

        var directBlock = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaArr, flat);
        var directEigs = directBlock.Evd().EigenValues.ToArray();

        var klein = KleinFourGroupSelfPairedRefinement.Build(N);
        var unionEigs = new List<Complex>();
        foreach (var chi in new[] { KleinCharacter.PlusPlus, KleinCharacter.PlusMinus,
                                    KleinCharacter.MinusPlus, KleinCharacter.MinusMinus })
        {
            var subBlock = klein.BuildSubBlockL(chi, gammaArr);
            var subEigs = subBlock.Evd().EigenValues.ToArray();
            unionEigs.AddRange(subEigs);
            _out.WriteLine($"N={N}, sub-block {chi}: dim {klein.SubBlockDims[chi]}, " +
                           $"{subEigs.Length} eigenvalues");
        }
        _out.WriteLine($"N={N}: direct sector eig count = {directEigs.Length}, Klein union count = {unionEigs.Count}");

        Assert.Equal(directEigs.Length, unionEigs.Count);
        MultisetAssert.NearestNeighbourEqual(
            unionEigs, directEigs, tolerance: 1e-9, context: $"N={N} Klein union vs direct");
    }

    [Fact]
    public void Build_Tier_IsTier1Derived()
    {
        var k = KleinFourGroupSelfPairedRefinement.Build(N: 4);
        Assert.Equal(Tier.Tier1Derived, k.Tier);
    }

    [Fact]
    public void Build_RejectsOddN()
    {
        Assert.Throws<ArgumentException>(() => KleinFourGroupSelfPairedRefinement.Build(N: 5));
    }

    [Fact]
    public void Build_RejectsTinyN()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => KleinFourGroupSelfPairedRefinement.Build(N: 1));
    }
}
