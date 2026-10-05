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

/// <summary>Dephase-letter tests for
/// <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(ComplexMatrix, System.Collections.Generic.IReadOnlyList{double}, int, LiouvillianBlockSpectrum.EigenPath, PauliLetter)"/>
/// and the matching overload on <see cref="F71MirrorBlockRefinement"/>.
///
/// <para>Task B of the 2026-05-25 BlockSpectrum / F108 / per-bond J wave (spec at
/// <c>docs/superpowers/specs/2026-05-25-builder-f108-jbond-wiring-design.md</c>). Adds a
/// <see cref="PauliLetter"/> <c>dephaseLetter</c> parameter and verifies its contract:</para>
/// <list type="number">
///   <item>Backward compat: the default Z overload matches the full eigensolver within
///         1e-9 (XY chain regression).</item>
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
/// diagonal only under Z-dephasing. The <c>dephaseLetter</c> parameter records the caller's
/// intent, and the builder refuses non-Z combos rather than producing wrong eigenvalues
/// (see Task-B B.1 finding in the commit message and the entry-point XML doc).</para></summary>
public class F108DispatchTests
{
    // ----------------------------------------------------------------------
    // Regression: the default Z overload against the full eigensolver.
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void XyChain_DefaultZ_MatchesFullEigensolver_Regression(int N)
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
            context: $"XY chain N={N} default-Z");
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void XyChain_ExplicitZ_MatchesDefault_BitExact(int N)
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
    public void F71Refinement_XyChain_ExplicitZ_MatchesDefault_BitExact(int N)
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
    public void F71Refinement_XyChain_ExplicitZ_MatchesDenseFullEig(int N)
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
            context: $"F71 refinement XY chain N={N} explicit-Z");
    }

    // ----------------------------------------------------------------------
    // End to end: the default (no-dephaseLetter) overload and the explicit-Z 5-arg overload
    // run one code path and agree bit for bit on a LAPACK build that returns the same bits for
    // the same input, as this one does, so a switch of either to another pairing would show
    // here as a difference; which pairing each overload applies is
    // checked bit for bit in BlockSpectrumPairingGuardTests
    // (DefaultOverloads_UseThePiOrbitPairing for these two, RequestedPairing_ReachesBothEngines
    // for the SectorPairing overloads).
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void Heisenberg_DefaultAndExplicitZ_Agree_EndToEnd(int N)
    {
        // Use the Heisenberg chain (popcount-conserving) so the per-block builder's exact
        // popcount check passes.
        var H = PauliHamiltonian.HeisenbergChain(N, 1.0).ToMatrix();
        var gamma = Enumerable.Repeat(0.05, N).ToArray();

        var defaultSpectrum = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gamma, N);
        var explicitZSpectrum = LiouvillianBlockSpectrum
            .ComputeSpectrumPerBlock(H, gamma, N, LiouvillianBlockSpectrum.EigenPath.Auto, PauliLetter.Z);

        Assert.Equal(defaultSpectrum.Length, explicitZSpectrum.Length);
        for (int i = 0; i < defaultSpectrum.Length; i++)
            Assert.Equal(defaultSpectrum[i], explicitZSpectrum[i]);
    }
}
