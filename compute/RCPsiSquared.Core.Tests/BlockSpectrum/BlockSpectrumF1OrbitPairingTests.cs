using System.Collections.Generic;
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

/// <summary>Bit-exact integration tests for the F1 Π-orbit pairing optimisation in
/// <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock"/> and
/// <see cref="F71MirrorBlockRefinement.ComputeSpectrumPerBlock"/>.
///
/// <para>The F1 palindrome conjugation Π is order-4 and on joint-popcount labels acts as
/// the whole-sector cycle (p_c, p_r) ↦ (N − p_r, p_c), grouping the (N+1)² sectors into
/// orbits of 4. One eigendecomposition per orbit feeds three followers: the Π²-image
/// (X⊗N partner) by a verbatim copy, the Π/Π³-images by the F1 reflection λ ↦ −2Σγ − λ.
/// These tests pin that the union of all 4^N eigenvalues stays bit-exactly equal to the
/// dense full-L eig (N=3..6) and to the un-halved per-block loop (N=5,6,7) where dense L
/// is infeasible.</para></summary>
public class BlockSpectrumF1OrbitPairingTests
{
    // ----------------------------------------------------------------------
    // Path A: LiouvillianBlockSpectrum.ComputeSpectrumPerBlock vs dense L.Evd
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void ComputeSpectrumPerBlock_WithF1OrbitPairing_MatchesFullLEig(int N)
    {
        const double J = 1.0;
        const double gamma = 0.5;
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();
        var L = PauliDephasingDissipator.BuildZ(H, gammaPerSite);

        var spectrumOrbit = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaPerSite, N);
        var spectrumFull = L.Evd().EigenValues.ToArray();

        Assert.Equal(spectrumOrbit.Length, spectrumFull.Length);
        MultisetAssert.NearestNeighbourEqual(spectrumOrbit, spectrumFull, tolerance: 1e-9, context: $"N={N}");
    }

    // ----------------------------------------------------------------------
    // Path B: F71-refined ComputeSpectrumPerBlock vs dense L.Evd
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void F71RefinedComputeSpectrumPerBlock_WithF1OrbitPairing_MatchesFullLEig(int N)
    {
        const double J = 1.0;
        const double gamma = 0.5;
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();
        var L = PauliDephasingDissipator.BuildZ(H, gammaPerSite);

        var spectrumOrbit = F71MirrorBlockRefinement.ComputeSpectrumPerBlock(H, gammaPerSite, N);
        var spectrumFull = L.Evd().EigenValues.ToArray();

        Assert.Equal(spectrumOrbit.Length, spectrumFull.Length);
        MultisetAssert.NearestNeighbourEqual(spectrumOrbit, spectrumFull, tolerance: 1e-9, context: $"N={N}");
    }

    // Path B non-uniform γ guard: the F1 reflection of both the F71-even and F71-odd
    // sub-arrays must use the genuine Σγ. The γ-list here is non-uniform (so Σγ ≠ N·γ,
    // catching a wrong N·γ reflection constant) but PALINDROMIC across the chain mirror
    // (γ_l = γ_{N−1−l}). The palindromic constraint is required: F71 spatial-mirror
    // refinement is exact iff γ_l = γ_{N−1−l} (Tier-1 InhomogeneousGammaF71BreakingWitness);
    // a non-palindromic γ breaks the F71 refinement itself, independently of the pairing.
    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void F71RefinedComputeSpectrumPerBlock_WithF1OrbitPairing_PalindromicNonUniformGamma_MatchesFullLEig(int N)
    {
        const double J = 1.0;
        // γ_l = γ_{N−1−l}: a palindromic ramp toward the chain centre, non-uniform.
        var gammaPerSite = Enumerable.Range(0, N)
            .Select(l => 0.1 + 0.17 * Math.Min(l, N - 1 - l))
            .ToArray();
        Assert.Equal(0.0, InhomogeneousGammaF71BreakingWitness.F71AsymmetryNorm(gammaPerSite), 12);
        Assert.True(gammaPerSite.Distinct().Count() > 1, "γ must be non-uniform for the Σγ guard.");

        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var L = PauliDephasingDissipator.BuildZ(H, gammaPerSite);

        var spectrumOrbit = F71MirrorBlockRefinement.ComputeSpectrumPerBlock(H, gammaPerSite, N);
        var spectrumFull = L.Evd().EigenValues.ToArray();

        Assert.Equal(spectrumOrbit.Length, spectrumFull.Length);
        MultisetAssert.NearestNeighbourEqual(spectrumOrbit, spectrumFull, tolerance: 1e-9, context: $"N={N}");
    }

    // Path B decisive test: F1-orbit-halved F71-refined path vs the un-halved F71-refined
    // per-block reference. At N=7 dense L.Evd is infeasible.
    [Theory]
    [InlineData(5, 1e-9)]
    [InlineData(6, 1e-9)]
    [InlineData(7, 1e-7)]
    public void F71RefinedComputeSpectrumPerBlock_WithF1OrbitPairing_MatchesUnHalvedReference(
        int N, double tolerance)
    {
        const double J = 1.0;
        const double gamma = 0.5;
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();

        var spectrumOrbit = F71MirrorBlockRefinement.ComputeSpectrumPerBlock(H, gammaPerSite, N);
        var spectrumReference = ComputeF71RefinedSpectrumNoPairing(H, gammaPerSite, N);

        Assert.Equal(spectrumOrbit.Length, spectrumReference.Length);
        MultisetAssert.NearestNeighbourEqual(
            spectrumOrbit, spectrumReference, tolerance, context: $"N={N}");
    }

    // ----------------------------------------------------------------------
    // Decisive test: F1-orbit-halved path vs an un-halved per-block reference.
    // At N=7 dense L.Evd is infeasible (4^7 = 16384, ~4 GB), so the reference is
    // the Phase-2 loop over ALL sectors with no follower derivation.
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(5, 1e-9)]
    [InlineData(6, 1e-9)]
    [InlineData(7, 1e-7)]
    public void ComputeSpectrumPerBlock_WithF1OrbitPairing_MatchesUnHalvedPerBlockReference(
        int N, double tolerance)
    {
        const double J = 1.0;
        const double gamma = 0.5;
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();

        var spectrumOrbit = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaPerSite, N);
        var spectrumReference = ComputeSpectrumPerBlockNoPairing(H, gammaPerSite, N);

        Assert.Equal(spectrumOrbit.Length, spectrumReference.Length);
        MultisetAssert.NearestNeighbourEqual(
            spectrumOrbit, spectrumReference, tolerance, context: $"N={N}");
    }

    // Non-uniform γ guard: the F1 reflection constant must be the genuine Σγ, not N·γ.
    // A wrong constant survives a uniform-γ test (Σγ = N·γ there) but fails here.
    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    public void ComputeSpectrumPerBlock_WithF1OrbitPairing_NonUniformGamma_MatchesFullLEig(int N)
    {
        const double J = 1.0;
        // Distinct per-site rates so Σγ ≠ N·γ for any single γ.
        var gammaPerSite = Enumerable.Range(0, N).Select(l => 0.1 + 0.13 * l).ToArray();
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var L = PauliDephasingDissipator.BuildZ(H, gammaPerSite);

        var spectrumOrbit = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaPerSite, N);
        var spectrumFull = L.Evd().EigenValues.ToArray();

        Assert.Equal(spectrumOrbit.Length, spectrumFull.Length);
        MultisetAssert.NearestNeighbourEqual(spectrumOrbit, spectrumFull, tolerance: 1e-9, context: $"N={N}");
    }

    // ----------------------------------------------------------------------
    // N=7 F1 self-consistency: the full spectrum equals its own λ ↦ −2Σγ − λ image.
    // ----------------------------------------------------------------------

    [Fact]
    public void ComputeSpectrumPerBlock_AtN7_IsClosedUnderF1Reflection()
    {
        const int N = 7;
        const double J = 1.0;
        const double gamma = 0.5;
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();
        double sumGamma = gammaPerSite.Sum();

        var spectrum = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaPerSite, N);
        var reflected = spectrum
            .Select(z => new Complex(-2.0 * sumGamma - z.Real, -z.Imaginary))
            .ToArray();

        // F1: the spectrum is a multiset palindrome about −Σγ. The reflected multiset must
        // equal the original multiset.
        MultisetAssert.NearestNeighbourEqual(spectrum, reflected, tolerance: 1e-7, context: "N=7 F1");
    }

    // ----------------------------------------------------------------------
    // Follower-derivation test: each follower sector's eigenvalue multiset equals
    // the F1-reflected (or verbatim-copied) primary, depending on the follower kind.
    // ----------------------------------------------------------------------

    [Fact]
    public void ComputeSpectrumPerBlock_AtN5_FollowerSectors_DeriveFromTheirOrbitPrimary()
    {
        const int N = 5;
        const double J = 1.0;
        const double gamma = 0.5;
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();
        double sumGamma = gammaPerSite.Sum();

        var spectrum = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaPerSite, N);

        // Slice the flat spectrum back into per-sector eigenvalue lists.
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

        // For every follower sector, verify it equals the F1-reflected (or copied) primary.
        var (_, followerToPrimary) = F1PalindromeOrbitPairing.PartitionByPiOrbit(
            N, decomp.SectorRanges, s => (s.PCol, s.PRow));
        int checkedReflect = 0, checkedCopy = 0;
        foreach (var (followerIdx, follower) in followerToPrimary)
        {
            var followerLabel = (decomp.SectorRanges[followerIdx].PCol, decomp.SectorRanges[followerIdx].PRow);
            var primaryLabel = (decomp.SectorRanges[follower.PrimaryIndex].PCol, decomp.SectorRanges[follower.PrimaryIndex].PRow);
            var followerEigs = bySector[followerLabel];
            var primaryEigs = bySector[primaryLabel];
            Assert.Equal(primaryEigs.Count, followerEigs.Count);

            if (follower.Kind == F1PalindromeOrbitPairing.F1FollowerKind.F1Reflect)
            {
                var expected = primaryEigs
                    .Select(z => new Complex(-2.0 * sumGamma - z.Real, -z.Imaginary))
                    .ToList();
                MultisetAssert.NearestNeighbourEqual(
                    followerEigs, expected, tolerance: 1e-9,
                    context: $"F1Reflect follower {followerLabel} from primary {primaryLabel}");
                checkedReflect++;
            }
            else
            {
                MultisetAssert.NearestNeighbourEqual(
                    followerEigs, primaryEigs, tolerance: 1e-9,
                    context: $"XnCopy follower {followerLabel} from primary {primaryLabel}");
                checkedCopy++;
            }
        }
        // N=5 (odd): every orbit has size 4 ⇒ 3 followers per orbit, 2 F1Reflect + 1 XnCopy.
        Assert.True(checkedReflect > 0, "no F1Reflect followers were exercised");
        Assert.True(checkedCopy > 0, "no XnCopy followers were exercised");
        Assert.Equal(2 * checkedCopy, checkedReflect);
    }

    // ----------------------------------------------------------------------
    // F71-commutation guard (Path B prerequisite). The F1Reflect followers in
    // F71MirrorBlockRefinement reflect BOTH the F71-even and F71-odd sub-arrays;
    // this is only valid if Π commutes with the F71 mirror pair P_F71 ⊗ P_F71, so
    // that the Π sector-permutation maps even→even and odd→odd. Verify ‖ΠF − FΠ‖_F = 0
    // on the computational Liouville space, the basis F71MirrorBlockRefinement uses.
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void Pi_CommutesWith_F71MirrorPair(int N)
    {
        int d = 1 << N;
        int liouvilleDim = d * d;

        // Π in the computational Liouville basis. PiOperator.BuildFull gives Π in the
        // 4^N Pauli-string basis; the change of basis to the computational vec(ρ) basis
        // is W, whose column k is vec(σ_k) (row-major flat = row·d + col). The Pauli
        // strings are Hilbert-Schmidt orthogonal with ⟨σ_j,σ_k⟩ = 2^N·δ_jk, so W is
        // unitary up to scale: W†W = 2^N·I, hence W⁻¹ = W†/2^N and
        // Π_Liouville = W · Π_Pauli · W⁻¹.
        var piPauli = global::RCPsiSquared.Core.Symmetry.PiOperator.BuildFull(N);
        var W = Matrix<Complex>.Build.Dense(liouvilleDim, liouvilleDim);
        for (long k = 0; k < liouvilleDim; k++)
        {
            var sigma = PauliString.Build(PauliIndex.FromFlat(k, N));
            for (int row = 0; row < d; row++)
                for (int col = 0; col < d; col++)
                    W[row * d + col, (int)k] = sigma[row, col];
        }
        var Winv = W.ConjugateTranspose() / new Complex(d, 0);
        var piLiouville = W * (piPauli * Winv);

        // P_F71 ⊗ P_F71 in the computational Liouville basis: the bit-reversal site
        // mirror on each Hilbert side, flat ↦ mirror(row)·d + mirror(col). Permutation
        // matrix F with F[mirror(flat), flat] = 1.
        var mirrorBits = F71MirrorIndexHelper.BuildHilbertMirrorLookup(N);
        var F = Matrix<Complex>.Build.Sparse(liouvilleDim, liouvilleDim);
        for (int flat = 0; flat < liouvilleDim; flat++)
        {
            int mflat = mirrorBits[flat / d] * d + mirrorBits[flat % d];
            F[mflat, flat] = Complex.One;
        }

        double commutatorNorm = (piLiouville * F - F * piLiouville).FrobeniusNorm();
        Assert.True(commutatorNorm < 1e-10,
            $"N={N}: ‖Π·F − F·Π‖_F = {commutatorNorm:E3} (Π must commute with P_F71⊗P_F71 " +
            "for the F1Reflect followers to act even→even, odd→odd in Path B).");
    }

    // ----------------------------------------------------------------------
    // Test-local un-halved per-block reference. This is the Phase-2 loop of
    // LiouvillianBlockSpectrum.ComputeSpectrumPerBlock run over EVERY joint-popcount
    // sector, with no F1-orbit / X⊗N follower derivation. It is the independent
    // reference for the N=5,6,7 decisive test (dense L.Evd is infeasible at N=7).
    // ----------------------------------------------------------------------

    private static Complex[] ComputeSpectrumPerBlockNoPairing(
        ComplexMatrix H, IReadOnlyList<double> gammaPerSite, int N)
    {
        int liouvilleDim = 1 << (2 * N);
        var decomp = JointPopcountSectorBuilder.Build(N);
        var perm = decomp.Permutation;
        var spectrum = new Complex[liouvilleDim];

        int sectorCount = decomp.SectorRanges.Count;
        var writeOffsets = new int[sectorCount];
        int cum = 0;
        for (int i = 0; i < sectorCount; i++)
        {
            writeOffsets[i] = cum;
            cum += decomp.SectorRanges[i].Size;
        }

        for (int sIdx = 0; sIdx < sectorCount; sIdx++)
        {
            var sector = decomp.SectorRanges[sIdx];
            int size = sector.Size;
            if (size == 0) continue;
            var flatIndices = new int[size];
            for (int k = 0; k < size; k++)
                flatIndices[k] = perm[sector.Offset + k];
            var block = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaPerSite, flatIndices);
            var blockEigs = block.Evd().EigenValues;
            int write = writeOffsets[sIdx];
            for (int i = 0; i < size; i++) spectrum[write + i] = blockEigs[i];
        }
        return spectrum;
    }

    // Test-local un-halved F71-refined per-block reference: the Phase-2 loop of
    // F71MirrorBlockRefinement.ComputeSpectrumPerBlock run over EVERY joint-popcount sector
    // with no F1-orbit / X⊗N follower derivation. Independent reference for the Path-B
    // decisive N=5,6,7 test (dense L.Evd is infeasible at N=7).
    private static Complex[] ComputeF71RefinedSpectrumNoPairing(
        ComplexMatrix H, IReadOnlyList<double> gammaPerSite, int N)
    {
        int d = 1 << N;
        int liouvilleDim = d * d;
        var mirrorBits = F71MirrorIndexHelper.BuildHilbertMirrorLookup(N);
        int Mirror(int flat) => mirrorBits[flat / d] * d + mirrorBits[flat % d];

        var baseDecomp = JointPopcountSectorBuilder.Build(N);
        var spectrum = new Complex[liouvilleDim];

        int sectorCount = baseDecomp.SectorRanges.Count;
        var writeOffsets = new int[sectorCount];
        int cum = 0;
        for (int i = 0; i < sectorCount; i++)
        {
            writeOffsets[i] = cum;
            cum += baseDecomp.SectorRanges[i].Size;
        }

        for (int sIdx = 0; sIdx < sectorCount; sIdx++)
        {
            var sector = baseDecomp.SectorRanges[sIdx];
            int size = sector.Size;
            if (size == 0) continue;

            var sectorFlat = new int[size];
            for (int k = 0; k < size; k++) sectorFlat[k] = baseDecomp.Permutation[sector.Offset + k];
            var (fixedPoints, pairs) = F71MirrorIndexHelper.FindOrbitsInSector(sectorFlat, Mirror);

            int nFix = fixedPoints.Count;
            int nPairs = pairs.Count;
            int unionSize = nFix + 2 * nPairs;
            var unionFlat = new int[unionSize];
            for (int i = 0; i < nFix; i++) unionFlat[i] = fixedPoints[i];
            for (int k = 0; k < nPairs; k++) unionFlat[nFix + k] = pairs[k].S;
            for (int k = 0; k < nPairs; k++) unionFlat[nFix + nPairs + k] = pairs[k].Ps;

            var unionBlock = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaPerSite, unionFlat);
            var rotated = F71MirrorBlockRefinement.RotateUnionBlockF71InPlace(unionBlock, nFix, nPairs);

            int evenSize = nFix + nPairs;
            int oddSize = nPairs;
            int write = writeOffsets[sIdx];
            if (evenSize > 0)
            {
                var evenEigs = rotated.SubMatrix(0, evenSize, 0, evenSize).Evd().EigenValues;
                for (int i = 0; i < evenSize; i++) spectrum[write + i] = evenEigs[i];
                write += evenSize;
            }
            if (oddSize > 0)
            {
                var oddEigs = rotated.SubMatrix(evenSize, oddSize, evenSize, oddSize).Evd().EigenValues;
                for (int i = 0; i < oddSize; i++) spectrum[write + i] = oddEigs[i];
            }
        }
        return spectrum;
    }
}
