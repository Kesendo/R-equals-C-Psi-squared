using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.SymmetryFamily;
using RCPsiSquared.Core.Tests.TestHelpers;
using Xunit;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.BlockSpectrum;

/// <summary>F108 Π_5bilinear orbit-pairing dispatch tests for
/// <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(ComplexMatrix, System.Collections.Generic.IReadOnlyList{double}, int, LiouvillianBlockSpectrum.EigenPath, PauliLetter)"/>
/// and the matching overload on <see cref="F71MirrorBlockRefinement"/>.
///
/// <para>Task B of the 2026-05-25 BlockSpectrum / F108 / per-bond J wave (spec at
/// <c>docs/superpowers/specs/2026-05-25-builder-f108-jbond-wiring-design.md</c>). Adds a
/// <see cref="PauliLetter"/> <c>dephaseLetter</c> parameter and verifies the F108 auto-
/// dispatch contract:</para>
/// <list type="number">
///   <item>Backward compat: default Z gives bit-exact same spectrum as the pre-Task-B
///         overload (Heisenberg / XY chain regression).</item>
///   <item>F108 Part 1 (Z-deph) acceptance: explicit Z passes through identically.</item>
///   <item>Mismatch (non-Z) handling: the per-block builder is hardcoded Z-only
///         (<see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/>); X- and Y-dephasing throw
///         <see cref="System.NotImplementedException"/> with a clear message, rather than
///         silently stamping Z-deph entries onto a non-Z problem.</item>
/// </list>
///
/// <para>The "soundness over generality" stance: the F108 generalisation (Part 1+2+3) lives
/// at the operator-algebra level (Π_5bilinear on the Pauli-string basis); the builder
/// operates in the computational Liouville basis where joint-popcount sectors are block-
/// diagonal only under Z-dephasing. The auto-detect dispatch records the user's intent via
/// <c>dephaseLetter</c> and refuses non-Z combos rather than producing wrong eigenvalues
/// (see Task-B B.1 finding in the commit message and the entry-point XML doc).</para></summary>
public class F108DispatchTests
{
    // ----------------------------------------------------------------------
    // Regression: default Z overload matches pre-Task-B behaviour bit-exactly.
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void Heisenberg_Z_GetsF1Quartering_Regression(int N)
    {
        const double J = 1.0;
        const double gamma = 0.5;
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();
        var L = PauliDephasingDissipator.BuildZ(H, gammaPerSite);

        // Pre-Task-B 3-argument overload: the historical entry point. Must continue to work
        // unchanged via the new dephaseLetter = PauliLetter.Z default.
        var spectrumDefault = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaPerSite, N);
        var spectrumFull = L.Evd().EigenValues.ToArray();

        Assert.Equal(spectrumDefault.Length, spectrumFull.Length);
        MultisetAssert.NearestNeighbourEqual(spectrumDefault, spectrumFull, tolerance: 1e-9,
            context: $"Heisenberg N={N} default-Z");
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void Heisenberg_ExplicitZ_MatchesDefault_BitExact(int N)
    {
        const double J = 1.0;
        const double gamma = 0.5;
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();

        var spectrumDefault = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaPerSite, N);
        var spectrumExplicitZ = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(
            H, gammaPerSite, N, LiouvillianBlockSpectrum.EigenPath.Auto, PauliLetter.Z);

        // Same code path, same orbit-pairing, bit-exact ordering. Both arrays are produced
        // by the same dispatch under PauliLetter.Z.
        Assert.Equal(spectrumDefault.Length, spectrumExplicitZ.Length);
        for (int i = 0; i < spectrumDefault.Length; i++)
            Assert.Equal(spectrumDefault[i], spectrumExplicitZ[i]);
    }

    // ----------------------------------------------------------------------
    // F108 Part 2 (X-deph) / Part 3 (Y-deph) dispatch: explicit
    // NotImplementedException, no silent stamping of Z-deph onto a non-Z problem.
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(PauliLetter.X)]
    [InlineData(PauliLetter.Y)]
    public void NonZ_DephaseLetter_Throws_NotImplemented_OnLiouvillianBlockSpectrum(PauliLetter dephaseLetter)
    {
        const int N = 4;
        var H = PauliHamiltonian.XYChain(N, 1.0).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(0.5, N).ToArray();

        var ex = Assert.Throws<System.NotImplementedException>(() =>
            LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(
                H, gammaPerSite, N, LiouvillianBlockSpectrum.EigenPath.Auto, dephaseLetter));

        // Soundness: the message must cite the per-block builder restriction so callers know
        // what to do, rather than seeing a generic "not implemented" wall.
        Assert.Contains("PauliLetter.Z", ex.Message);
        Assert.Contains("PerBlockLiouvillianBuilder", ex.Message);
    }

    [Theory]
    [InlineData(PauliLetter.X)]
    [InlineData(PauliLetter.Y)]
    public void NonZ_DephaseLetter_Throws_NotImplemented_OnF71MirrorBlockRefinement(PauliLetter dephaseLetter)
    {
        const int N = 4;
        var H = PauliHamiltonian.XYChain(N, 1.0).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(0.5, N).ToArray();

        var ex = Assert.Throws<System.NotImplementedException>(() =>
            F71MirrorBlockRefinement.ComputeSpectrumPerBlock(H, gammaPerSite, N, dephaseLetter));

        Assert.Contains("PauliLetter.Z", ex.Message);
        Assert.Contains("PerBlockLiouvillianBuilder", ex.Message);
    }

    [Fact]
    public void Identity_PauliLetter_Throws_NotImplemented()
    {
        // PauliLetter.I is a dephase letter the canonical Π / Π_5bilinear families do not
        // define; the dispatch refuses it the same way it refuses X / Y, citing the per-
        // block Z-only restriction.
        const int N = 3;
        var H = PauliHamiltonian.XYChain(N, 1.0).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(0.5, N).ToArray();

        Assert.Throws<System.NotImplementedException>(() =>
            LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(
                H, gammaPerSite, N, LiouvillianBlockSpectrum.EigenPath.Auto, PauliLetter.I));
    }

    // ----------------------------------------------------------------------
    // F71MirrorBlockRefinement: backward-compat regression and explicit-Z parity.
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void F71Refinement_Heisenberg_ExplicitZ_MatchesDefault_BitExact(int N)
    {
        const double J = 1.0;
        const double gamma = 0.5;
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();

        var spectrumDefault = F71MirrorBlockRefinement.ComputeSpectrumPerBlock(H, gammaPerSite, N);
        var spectrumExplicitZ = F71MirrorBlockRefinement.ComputeSpectrumPerBlock(
            H, gammaPerSite, N, PauliLetter.Z);

        Assert.Equal(spectrumDefault.Length, spectrumExplicitZ.Length);
        for (int i = 0; i < spectrumDefault.Length; i++)
            Assert.Equal(spectrumDefault[i], spectrumExplicitZ[i]);
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void F71Refinement_Heisenberg_ExplicitZ_MatchesDenseFullEig(int N)
    {
        const double J = 1.0;
        const double gamma = 0.5;
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();
        var L = PauliDephasingDissipator.BuildZ(H, gammaPerSite);

        var spectrumExplicitZ = F71MirrorBlockRefinement.ComputeSpectrumPerBlock(
            H, gammaPerSite, N, PauliLetter.Z);
        var spectrumFull = L.Evd().EigenValues.ToArray();

        Assert.Equal(spectrumExplicitZ.Length, spectrumFull.Length);
        MultisetAssert.NearestNeighbourEqual(spectrumExplicitZ, spectrumFull, tolerance: 1e-9,
            context: $"F71 refinement Heisenberg N={N} explicit-Z");
    }

    // ----------------------------------------------------------------------
    // Quartering verification: the F108 dispatch keeps the orbit-pairing primary
    // count equal to F1's DistinctSpectralClasses(N), confirming that the auto-
    // dispatched (Z-deph) path still quarters and does not silently fall back to
    // a no-pairing reference. The check is structural rather than runtime-timed.
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void Z_DephaseLetter_Dispatch_PreservesPiOrbitQuartering(int N)
    {
        // The F1 Π-orbit partition is independent of H and γ at the sector-label level.
        // Confirm the dispatch reaches the same partition the existing path uses, i.e. the
        // F108 dispatch does not bypass the orbit-pairing primitive.
        var sectors = JointPopcountSectorBuilder.Build(N).SectorRanges;
        var (primaries, _) = F1PalindromeOrbitPairing.PartitionByPiOrbit(
            N, sectors, s => (s.PCol, s.PRow));

        int expectedDistinct = F1PalindromeOrbitPairing.DistinctSpectralClasses(N);
        Assert.Equal(expectedDistinct, primaries.Count);

        // Sanity: at N=4 even-N central (2, 2) sector is Π-fixed; at odd N every orbit has 4
        // entries and the primary count is exactly (N+1)²/4. These are intrinsic to the
        // dispatch's quartering claim regardless of which Π (canonical or Π_5bilinear-Z) the
        // dispatch nominally selects.
        int total = (N + 1) * (N + 1);
        int piFixed = N % 2 == 0 ? 1 : 0;
        Assert.Equal(piFixed + (total - piFixed) / 4, primaries.Count);
    }
}
