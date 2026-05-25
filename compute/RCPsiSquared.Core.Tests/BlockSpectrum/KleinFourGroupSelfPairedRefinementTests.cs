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
    [InlineData(8)]
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

    [Theory]
    [InlineData(4, KleinCharacter.PlusPlus)]
    [InlineData(4, KleinCharacter.MinusMinus)]
    [InlineData(6, KleinCharacter.PlusPlus)]
    public void BuildSubBlockL_UniformBondJEqualsOne_MatchesScalarOverload(int N, KleinCharacter chi)
    {
        // Regression: passing bondJ = [1.0, 1.0, ..., 1.0] (length N − 1) through the new
        // per-bond overload must reproduce the scalar J = 1 overload bit-exact (the scalar
        // overload itself just forwards to the per-bond path with a uniform list).
        const double gamma = 0.05;
        var gammaArr = Enumerable.Repeat(gamma, N).ToArray();
        var refinement = KleinFourGroupSelfPairedRefinement.Build(N);

        var bondJUniform = Enumerable.Repeat(1.0, N - 1).ToArray();
        var perBondL = refinement.BuildSubBlockL(chi, gammaArr, bondJUniform);
        var scalarL = refinement.BuildSubBlockL(chi, gammaArr);

        double diff = (perBondL - scalarL).FrobeniusNorm();
        _out.WriteLine($"N={N}, χ={chi}: ‖per-bond[1,..,1] − scalar J=1‖_F = {diff:G3}");
        Assert.True(diff < 1e-12,
            $"Uniform bondJ=[1,..,1] should match scalar J=1 bit-exact; got Frobenius diff {diff:G3}");
    }

    [Fact]
    public void BuildSubBlockL_NonUniformBondJ_DiffersFromUniform()
    {
        // Capability: bondJ = [1.0, 2.0, 1.0] at N=4 should produce a sub-block spectrum
        // that differs structurally from uniform bondJ = [1.0, 1.0, 1.0]. At N=4 the (2, 2)
        // sector has dim 36 split as ++/--/+-/-+ = [10, 10, 8, 8]; we pick the ++ sub-block.
        const int N = 4;
        const double gamma = 0.05;
        var gammaArr = Enumerable.Repeat(gamma, N).ToArray();
        var refinement = KleinFourGroupSelfPairedRefinement.Build(N);

        var bondJUniform = new[] { 1.0, 1.0, 1.0 };
        var bondJNonUniform = new[] { 1.0, 2.0, 1.0 };

        var uniformL = refinement.BuildSubBlockL(KleinCharacter.PlusPlus, gammaArr, bondJUniform);
        var nonUniformL = refinement.BuildSubBlockL(KleinCharacter.PlusPlus, gammaArr, bondJNonUniform);

        var uniformEigs = uniformL.Evd().EigenValues.ToArray();
        var nonUniformEigs = nonUniformL.Evd().EigenValues.ToArray();
        Assert.Equal(uniformEigs.Length, nonUniformEigs.Length);

        // Find at least one uniform eigenvalue whose nearest non-uniform partner is > 1e-3 away
        // (i.e., the multiset differs by more than numerical noise).
        var taken = new bool[nonUniformEigs.Length];
        double maxNearestDiff = 0.0;
        for (int i = 0; i < uniformEigs.Length; i++)
        {
            int bestIdx = -1;
            double bestDist = double.MaxValue;
            for (int j = 0; j < nonUniformEigs.Length; j++)
            {
                if (taken[j]) continue;
                double d = (uniformEigs[i] - nonUniformEigs[j]).Magnitude;
                if (d < bestDist) { bestDist = d; bestIdx = j; }
            }
            if (bestIdx >= 0) taken[bestIdx] = true;
            if (bestDist > maxNearestDiff) maxNearestDiff = bestDist;
        }
        _out.WriteLine($"N=4 χ=++: max nearest-neighbour spectrum diff (uniform vs [1,2,1]) = {maxNearestDiff:G3}");
        Assert.True(maxNearestDiff > 1e-3,
            $"Non-uniform bondJ=[1,2,1] should produce structurally different spectrum from uniform; " +
            $"got max nearest-neighbour diff {maxNearestDiff:G3} (expected > 1e-3)");
    }

    [Fact]
    public void BuildSubBlockL_BondJWrongLength_Throws()
    {
        const int N = 4;
        var gammaArr = Enumerable.Repeat(0.05, N).ToArray();
        var refinement = KleinFourGroupSelfPairedRefinement.Build(N);
        Assert.Throws<ArgumentException>(() =>
            refinement.BuildSubBlockL(KleinCharacter.PlusPlus, gammaArr, new double[2]));
    }

    [Fact]
    public void Build_AtN10_SubBlockBuildTimeIsTractable()
    {
        // Reconnaissance: confirm sub-block construction at N=10 (5, 5) is fast even
        // though dense Evd at that dim would be ~3 h per sub-block. We build just one
        // (++ sub-block, dim 16132) and verify:
        //   - matrix construction completes within a few minutes
        //   - the matrix is sparse-but-dense-storage (per-row nnz ≪ dim)
        //   - the matrix is non-zero (sanity)
        // Full sector spectrum requires 4 sub-block Evds at ~3 h each (~12 h total)
        // — see the SLOW_ test below for that overnight workload.
        const double gamma = 0.05;
        var gammaArr = Enumerable.Repeat(gamma, 10).ToArray();

        var klein = KleinFourGroupSelfPairedRefinement.Build(N: 10);
        Assert.Equal(16132, klein.SubBlockDims[KleinCharacter.PlusPlus]);

        var sw = System.Diagnostics.Stopwatch.StartNew();
        var subBlock = klein.BuildSubBlockL(KleinCharacter.PlusPlus, gammaArr);
        sw.Stop();
        Assert.Equal(16132, subBlock.RowCount);
        Assert.Equal(16132, subBlock.ColumnCount);

        int nnzSample = 0;
        double frobSqSample = 0;
        // Sample 64 random rows for nnz count + Frobenius² (avoiding MathNet FrobeniusNorm
        // which goes through MKL marshaling and trips the int32 array-size limit at this
        // sub-block dim: 16132² × 16 bytes = 4.2 GB > 2 GB managed-array cap).
        var rng = new Random(1);
        int rowsSampled = 64;
        for (int s = 0; s < rowsSampled; s++)
        {
            int row = rng.Next(16132);
            for (int col = 0; col < 16132; col++)
            {
                var v = subBlock[row, col];
                double mag2 = v.Real * v.Real + v.Imaginary * v.Imaginary;
                if (mag2 > 1e-24) nnzSample++;
                frobSqSample += mag2;
            }
        }
        double meanNnzPerRow = (double)nnzSample / rowsSampled;
        double frobNormSampled = Math.Sqrt(frobSqSample * 16132.0 / rowsSampled);

        _out.WriteLine($"N=10 (5,5) Klein ++ sub-block: dim 16132, build time {sw.Elapsed.TotalSeconds:F1} s");
        _out.WriteLine($"  mean nnz per row (sampled 64 rows): {meanNnzPerRow:F1}");
        _out.WriteLine($"  Frobenius norm (estimated from sample): {frobNormSampled:F1}");

        Assert.True(frobNormSampled > 1.0, "++ sub-block should be substantially non-zero");
        Assert.True(meanNnzPerRow < 50,
            $"sub-block should be sparse (≤ ~4N + a few per row); got mean {meanNnzPerRow}");
    }

    [Fact(Skip = "Blocked: dense complex Evd at sub-block dim 16132 trips the .NET/MKL int32 array-size cap (16132² × 16 bytes = 4.2 GB > 2 GB managed marshaling limit; same issue compute/RCPsiSquared.Compute solved at N=8 via NativeMemory + ILP64 LAPACK in MklDirect.cs). Implementing the same NativeMemory + ILP64 path for the Klein sub-block Evd is the natural Phase 3c. Once that exists, the workload itself is reachable (~3 h per sub-block dense Evd, ~12 h total, ~4 GB peak per sub-block — overnight on commodity hardware).")]
    public void SLOW_Build_AtN10_FullSectorSpectrumViaFourSubBlockEvd_F1PalindromeHolds()
    {
        // The Phase 3b target: at N=10 (5, 5) compute the full sector spectrum (63 504
        // eigenvalues) by independently eigendecomposing the 4 Klein character sub-blocks
        // (dims [16132, 15620, 15620, 16132]), then verify F1 palindrome λ → −λ − 2·Σγ on
        // the union. This is the workload the direct dense (5, 5) Evd cannot tackle
        // (~65 GB matrix); the Klein splitting brings it to 4 sub-blocks at ~4 GB each.
        // Per-sub-block dense complex Evd at dim ~16k extrapolates from N=8 timing as
        // (16132/1225)³ ≈ 2280× → ~3 h per sub-block, ~12 h total. For faster top-K
        // recovery (not full spectrum) combine the Klein splitting with the Phase 2 sparse
        // Krylov path (see JwSlaterPairShiftInvertArnoldi for the per-sub-block analogue).
        const double gamma = 0.05;
        const int N = 10;
        var gammaArr = Enumerable.Repeat(gamma, N).ToArray();
        double sumGamma = N * gamma;

        var klein = KleinFourGroupSelfPairedRefinement.Build(N);
        _out.WriteLine($"N=10 (5,5): Klein decomposition → 4 sub-blocks " +
                       $"[++: {klein.SubBlockDims[KleinCharacter.PlusPlus]}, " +
                       $"+-: {klein.SubBlockDims[KleinCharacter.PlusMinus]}, " +
                       $"-+: {klein.SubBlockDims[KleinCharacter.MinusPlus]}, " +
                       $"--: {klein.SubBlockDims[KleinCharacter.MinusMinus]}]");

        var unionEigs = new List<Complex>(63504);
        var sw = System.Diagnostics.Stopwatch.StartNew();
        foreach (var chi in new[] { KleinCharacter.PlusPlus, KleinCharacter.PlusMinus,
                                    KleinCharacter.MinusPlus, KleinCharacter.MinusMinus })
        {
            var swSub = System.Diagnostics.Stopwatch.StartNew();
            var subBlock = klein.BuildSubBlockL(chi, gammaArr);
            var subBuildTime = swSub.Elapsed;
            swSub.Restart();
            var subEigs = subBlock.Evd().EigenValues.ToArray();
            var subEvdTime = swSub.Elapsed;
            unionEigs.AddRange(subEigs);
            _out.WriteLine($"  sub-block {chi} dim {klein.SubBlockDims[chi]}: " +
                           $"build {subBuildTime.TotalSeconds:F1} s, Evd {subEvdTime.TotalSeconds:F1} s, " +
                           $"{subEigs.Length} eigenvalues");
        }
        sw.Stop();
        _out.WriteLine($"Total wall: {sw.Elapsed.TotalMinutes:F1} min, collected {unionEigs.Count} eigenvalues");

        Assert.Equal(63504, unionEigs.Count);

        // F1 palindrome verification: every λ in the collected set has its F1 mirror
        // −λ − 2·Σγ also in the set. Greedy nearest-neighbour matching for robustness
        // against the high degeneracies in the XY-chain spectrum.
        var mirrors = unionEigs.Select(l => new Complex(-l.Real - 2 * sumGamma, -l.Imaginary)).ToArray();
        MultisetAssert.NearestNeighbourEqual(
            mirrors, unionEigs, tolerance: 1e-7, context: "N=10 (5,5) F1 palindrome");
        _out.WriteLine($"F1 palindrome holds on full N=10 (5,5) spectrum within 1e-7");
    }
}
