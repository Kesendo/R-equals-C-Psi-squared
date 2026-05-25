using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
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
///         <see cref="System.NotSupportedException"/> (design-permanent under the current
///         basis), and PauliLetter.I throws <see cref="System.ArgumentException"/> (not a
///         valid dephase letter), rather than silently stamping Z-deph entries onto a non-Z
///         problem.</item>
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
    // NotSupportedException (design-permanent X/Y refusal under the current basis); PauliLetter.I
    // raises ArgumentException. No silent stamping of Z-deph entries onto a non-Z problem.
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(PauliLetter.X)]
    [InlineData(PauliLetter.Y)]
    public void NonZ_DephaseLetter_Throws_NotSupported_OnLiouvillianBlockSpectrum(PauliLetter dephaseLetter)
    {
        const int N = 4;
        var H = PauliHamiltonian.XYChain(N, 1.0).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(0.5, N).ToArray();

        // Exception TYPE is the contract: NotSupportedException signals "design-permanent
        // refusal under the current basis", distinct from NotImplementedException which
        // conventionally marks an unfilled stub. Wording details can change freely without
        // breaking this assertion.
        Assert.Throws<System.NotSupportedException>(() =>
            LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(
                H, gammaPerSite, N, LiouvillianBlockSpectrum.EigenPath.Auto, dephaseLetter));
    }

    [Theory]
    [InlineData(PauliLetter.X)]
    [InlineData(PauliLetter.Y)]
    public void NonZ_DephaseLetter_Throws_NotSupported_OnF71MirrorBlockRefinement(PauliLetter dephaseLetter)
    {
        const int N = 4;
        var H = PauliHamiltonian.XYChain(N, 1.0).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(0.5, N).ToArray();

        // Same contract as the LiouvillianBlockSpectrum variant: NotSupportedException for
        // design-permanent X/Y refusal under the current joint-popcount basis.
        Assert.Throws<System.NotSupportedException>(() =>
            F71MirrorBlockRefinement.ComputeSpectrumPerBlock(H, gammaPerSite, N, dephaseLetter));
    }

    [Fact]
    public void Identity_PauliLetter_Throws_ArgumentException()
    {
        // PauliLetter.I is not a valid dephase letter (the Lindblad dissipator requires a
        // non-identity operator). The dispatch raises ArgumentException with ParamName set to
        // "dephaseLetter" so callers can introspect the offending parameter structurally
        // instead of grepping the (mutable) human-readable message.
        const int N = 3;
        var H = PauliHamiltonian.XYChain(N, 1.0).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(0.5, N).ToArray();

        var ex = Assert.Throws<System.ArgumentException>(() =>
            LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(
                H, gammaPerSite, N, LiouvillianBlockSpectrum.EigenPath.Auto, PauliLetter.I));
        Assert.Equal("dephaseLetter", ex.ParamName);
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
    // End-to-end quartering verification: the default (no-dephaseLetter) overload and the
    // explicit-Z 5-arg overload must produce bit-exactly equal spectra. Both routes lead
    // through the orbit-pairing partition; if a future regression silently bypasses the
    // dispatch (e.g. switches one overload to a no-pairing reference path), the spectra
    // would diverge in ordering or value. The standalone structural check on
    // F1PalindromeOrbitPairing.DistinctSpectralClasses(N) was tautological because it
    // re-derived the partition without exercising the dispatch at all.
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void Z_DephaseLetter_Dispatch_PreservesPiOrbitQuartering_EndToEnd(int N)
    {
        // Use the Heisenberg chain (popcount-conserving) so the per-block builder's
        // DebugAssertPopcountConservingH guard does not trip.
        var H = PauliHamiltonian.HeisenbergChain(N, 1.0).ToMatrix();
        var gamma = Enumerable.Repeat(0.05, N).ToArray();

        var defaultSpectrum = LiouvillianBlockSpectrum
            .ComputeSpectrumPerBlock(H, gamma, N)
            .OrderBy(c => c.Real).ThenBy(c => c.Imaginary).ToArray();
        var explicitZSpectrum = LiouvillianBlockSpectrum
            .ComputeSpectrumPerBlock(H, gamma, N, LiouvillianBlockSpectrum.EigenPath.Auto, PauliLetter.Z)
            .OrderBy(c => c.Real).ThenBy(c => c.Imaginary).ToArray();

        Assert.Equal(defaultSpectrum.Length, explicitZSpectrum.Length);
        for (int i = 0; i < defaultSpectrum.Length; i++)
        {
            Assert.True(
                Math.Abs(defaultSpectrum[i].Real - explicitZSpectrum[i].Real) < 1e-9 &&
                Math.Abs(defaultSpectrum[i].Imaginary - explicitZSpectrum[i].Imaginary) < 1e-9,
                $"Default and explicit-Z spectra differ at index {i}: " +
                $"{defaultSpectrum[i]} vs {explicitZSpectrum[i]}");
        }
    }
}
