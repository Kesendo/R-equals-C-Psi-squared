using System.Numerics;
using System.Runtime.InteropServices;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Tests.TestHelpers;

namespace RCPsiSquared.Core.Tests.BlockSpectrum;

/// <summary>Parity witness for the LP64-bridge port (2026-05-18). At small N where both
/// eigensolver paths (MathNet's managed wrapper + the new MklDirect + NativeMemory + ILP64-aware
/// path) are exercisable, the per-block spectrum must agree across the two paths (the
/// spectra within 1e-9, the per-block elements exactly).
///
/// <para>The point of this test class is not to verify the F1 palindromic-pairing identity
/// itself (that's covered by <c>LiouvillianBlockSpectrumTests</c>, <c>F1GeneralTopologyN7…</c>,
/// <c>F1GeneralTopologyN8…</c>, and, through the bridge, <c>F1GeneralTopologyN9…</c>).
/// It's the cross-route witness for the bridge itself: forcing every block
/// through MklDirect on N=3, 5 must produce the same eigenvalues (as a multiset) and the same
/// per-cell native buffer (vs MathNet's column-major storage) the production-default Auto path
/// produces for those same blocks under the MathNet branch.</para>
///
/// <para>Anchors: <see cref="PerBlockLiouvillianBuilder.BuildBlockZIntoNativeMemory"/> (the
/// new native-buffer build path), <see cref="LiouvillianBlockSpectrum.EigenPath"/> (the
/// dispatch enum), <see cref="LiouvillianBlockSpectrum.Lp64ComplexCeiling"/> (the threshold).
/// Reference for the bridge motivation: 2026-05-18 memory entry
/// <c>reference_lp64_2gb_ceiling_bypass.md</c>.</para></summary>
public class PerBlockLiouvillianBuilderNativeMemoryParityTests
{
    /// <summary>BuildBlockZ (MathNet) vs BuildBlockZIntoNativeMemory must produce bit-exact
    /// matrix elements for every (i, j) of every joint-popcount sector at N=3, since both
    /// follow the same accumulation pattern with identical operand order (no associativity
    /// reorder).
    ///
    /// <para>The test walks every sector of <see cref="JointPopcountSectorBuilder.Build"/>(3),
    /// builds each block via both paths, and compares every cell exactly: the loops add
    /// identical contributions in the same order, so a cell that differs at all is a
    /// finding about the builders.</para></summary>
    [Fact]
    public unsafe void BuildBlockZIntoNativeMemory_MatchesBuildBlockZ_PerCell_AtN3()
    {
        const int N = 3;
        // A coupling and rates that round, so a change in the order of accumulation between
        // the two builders could show in some cell.
        const double J = 0.7;
        var H = PauliHamiltonian.XYChain(N, J: J).ToMatrix();
        var gammaPerSite = new[] { 0.3, 0.37, 0.41 };
        var decomp = JointPopcountSectorBuilder.Build(N);

        int totalCellsCompared = 0;
        double maxAbsDelta = 0.0;

        foreach (var sector in decomp.SectorRanges)
        {
            int size = sector.Size;
            if (size == 0) continue;

            var flatIndices = new int[size];
            for (int k = 0; k < size; k++)
                flatIndices[k] = decomp.Permutation[sector.Offset + k];

            var managed = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaPerSite, flatIndices);

            IntPtr nativePtr = PerBlockLiouvillianBuilder.BuildBlockZIntoNativeMemory(
                H, gammaPerSite, flatIndices);
            try
            {
                var nativeData = (Complex*)nativePtr;
                for (int i = 0; i < size; i++)
                {
                    for (int j = 0; j < size; j++)
                    {
                        // Column-major: native[j*size + i] == managed[i, j].
                        var nativeCell = nativeData[(long)j * size + i];
                        var managedCell = managed[i, j];
                        double delta = (nativeCell - managedCell).Magnitude;
                        if (delta > maxAbsDelta) maxAbsDelta = delta;
                        Assert.True(nativeCell == managedCell,
                            $"Block (p_c={sector.PCol}, p_r={sector.PRow}, size={size}) cell ({i}, {j}): " +
                            $"managed={managedCell}, native={nativeCell}, |Δ|={delta:E3}");
                        totalCellsCompared++;
                    }
                }
            }
            finally
            {
                NativeMemory.Free((void*)nativePtr);
            }
        }

        Assert.True(totalCellsCompared > 0, "expected at least one populated sector at N=3");
    }

    /// <summary>Auto-dispatch must equal forced MathNet at small N where both routes pick the
    /// MathNet branch (block sizes well below <see cref="LiouvillianBlockSpectrum.Lp64ComplexCeiling"/>).
    /// This is the "no regression" witness: production code paths that don't trip the ceiling
    /// keep their pre-bridge behaviour exactly.</summary>
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void ComputeSpectrumPerBlock_Auto_MatchesMathNet_AtSmallN(int N)
    {
        var H = PauliHamiltonian.XYChain(N, J: 1.0).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(0.5, N).ToArray();

        var eigsAuto = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(
            H, gammaPerSite, N, LiouvillianBlockSpectrum.EigenPath.Auto);
        var eigsMathNet = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(
            H, gammaPerSite, N, LiouvillianBlockSpectrum.EigenPath.MathNet);

        Assert.Equal(eigsMathNet.Length, eigsAuto.Length);
        // At Auto vs MathNet small-N, both branches hit the same MathNet code path at the
        // same outer parallelism, so the spectra are equal element by element, on a LAPACK
        // build that returns the same bits for the same input.
        for (int i = 0; i < eigsMathNet.Length; i++)
            Assert.Equal(eigsMathNet[i], eigsAuto[i]);
    }

    /// <summary>The cross-route witness. Force every block through MklDirect +
    /// NativeMemory at N=3, 4, 5 where the MathNet path is the production default, then
    /// compare the resulting spectra. The two routes use independent LAPACK invocations
    /// (MathNet's z_eigen wrapper vs MklDirect's zgeev_ direct P/Invoke); the multiset of
    /// eigenvalues must agree to MKL's per-block Evd noise floor.
    ///
    /// <para>1e-9 tolerance matches the existing <c>LiouvillianBlockSpectrumTests</c> envelope
    /// for full-L vs per-block comparisons. The actual N=3, 5 observed residual is typically
    /// 1e-13..1e-14 because the per-block QR/Schur iteration converges to the same eigenvalue
    /// set regardless of which LAPACK entry point invoked it (same algorithm, same input,
    /// modulo internal workspace-size choices).</para></summary>
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void ComputeSpectrumPerBlock_MklDirectNative_MatchesMathNet_AtSmallN(int N)
    {
        var H = PauliHamiltonian.XYChain(N, J: 1.0).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(0.5, N).ToArray();

        var eigsMathNet = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(
            H, gammaPerSite, N, LiouvillianBlockSpectrum.EigenPath.MathNet);
        var eigsMklDirect = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(
            H, gammaPerSite, N, LiouvillianBlockSpectrum.EigenPath.MklDirectNative);

        Assert.Equal(eigsMathNet.Length, eigsMklDirect.Length);
        MultisetAssert.NearestNeighbourEqual(
            eigsMklDirect, eigsMathNet, tolerance: 1e-9, context: $"MklDirect vs MathNet N={N}");
    }

    /// <summary>Independent of the cross-route check: the forced MklDirect path must agree
    /// with the dense full-L baseline at N=3, 4, 5 to the same tolerance the parent
    /// <see cref="LiouvillianBlockSpectrum"/> claim uses for its witness. This
    /// closes the loop: MklDirect-per-block ≈ MathNet-per-block ≈ full dense L Evd, all
    /// within 1e-9.</summary>
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void ComputeSpectrumPerBlock_MklDirectNative_MatchesFullL_AtSmallN(int N)
    {
        var H = PauliHamiltonian.XYChain(N, J: 1.0).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(0.5, N).ToArray();
        var L = PauliDephasingDissipator.BuildZ(H, gammaPerSite);

        var eigsFull = L.Evd().EigenValues.ToArray();
        var eigsMklDirect = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(
            H, gammaPerSite, N, LiouvillianBlockSpectrum.EigenPath.MklDirectNative);

        Assert.Equal(eigsFull.Length, eigsMklDirect.Length);
        MultisetAssert.NearestNeighbourEqual(
            eigsMklDirect, eigsFull, tolerance: 1e-9, context: $"MklDirect vs full-L N={N}");
    }

    /// <summary>Threshold sanity check: the bridge constant must be reachable from a public
    /// surface (the N=9 test references it via class-local copy
    /// <c>F1GeneralTopologyN9BlockSpectrumChainTests.Lp64EvdSquareMatrixCeiling</c>), so a
    /// regression that drifts the central constant gets caught at build time inside this
    /// assembly.</summary>
    [Fact]
    public void Lp64ComplexCeiling_Is11585()
    {
        Assert.Equal(11_585, LiouvillianBlockSpectrum.Lp64ComplexCeiling);
    }
}
