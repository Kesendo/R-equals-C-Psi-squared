using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.BlockSpectrum;

/// <summary>Tests for <see cref="InhomogeneousGammaF71BreakingWitness"/>: GR-analogy
/// γ-asymmetry test. Predicts that joint-popcount block-diagonality is γ-blind, F71
/// refinement is exact iff γ palindromic, and F1 palindrome center stays at -Σγ regardless.
/// </summary>
public class InhomogeneousGammaF71BreakingWitnessTests
{
    // ----------------------------------------------------------------------
    // Claim metadata
    // ----------------------------------------------------------------------

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var claim = new InhomogeneousGammaF71BreakingWitness(
            new JointPopcountSectors(),
            new F71MirrorBlockRefinement(new JointPopcountSectors()));
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void DisplayName_IsNonEmpty()
    {
        var claim = new InhomogeneousGammaF71BreakingWitness(
            new JointPopcountSectors(),
            new F71MirrorBlockRefinement(new JointPopcountSectors()));
        Assert.False(string.IsNullOrWhiteSpace(claim.DisplayName));
    }

    [Fact]
    public void Constructor_NullSectors_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new InhomogeneousGammaF71BreakingWitness(
                null!, new F71MirrorBlockRefinement(new JointPopcountSectors())));
    }

    [Fact]
    public void Constructor_NullF71_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new InhomogeneousGammaF71BreakingWitness(new JointPopcountSectors(), null!));
    }

    // ----------------------------------------------------------------------
    // F71AsymmetryNorm: closed-form correctness on known γ-distributions
    // ----------------------------------------------------------------------

    [Fact]
    public void F71AsymmetryNorm_NullList_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            InhomogeneousGammaF71BreakingWitness.F71AsymmetryNorm(null!));
    }

    [Fact]
    public void F71AsymmetryNorm_Uniform_IsZero()
    {
        // [0.5, 0.5, 0.5, 0.5, 0.5] — γ_l = γ_{N-1-l} trivially.
        var gammas = new[] { 0.5, 0.5, 0.5, 0.5, 0.5 };
        double norm = InhomogeneousGammaF71BreakingWitness.F71AsymmetryNorm(gammas);
        Assert.True(norm < 1e-12, $"Uniform γ: expected asymmetry norm ≈ 0, got {norm:E3}.");
    }

    [Fact]
    public void F71AsymmetryNorm_PalindromicNonUniform_IsZero()
    {
        // [0.3, 0.4, 0.5, 0.4, 0.3] — γ_l = γ_{N-1-l} symbol-by-symbol.
        var gammas = new[] { 0.3, 0.4, 0.5, 0.4, 0.3 };
        double norm = InhomogeneousGammaF71BreakingWitness.F71AsymmetryNorm(gammas);
        Assert.True(norm < 1e-12, $"Palindromic γ: expected asymmetry norm ≈ 0, got {norm:E3}.");
    }

    [Fact]
    public void F71AsymmetryNorm_LinearMonotonic_N5_MatchesClosedForm()
    {
        // [0.3, 0.4, 0.5, 0.6, 0.7] at N=5.
        //   l=0: (0.3 − 0.7)² = 0.16
        //   l=1: (0.4 − 0.6)² = 0.04
        //   l=2: (0.5 − 0.5)² = 0
        //   l=3: (0.6 − 0.4)² = 0.04
        //   l=4: (0.7 − 0.3)² = 0.16
        // Σ = 0.40, sqrt = sqrt(0.4) ≈ 0.6324555.
        var gammas = new[] { 0.3, 0.4, 0.5, 0.6, 0.7 };
        double norm = InhomogeneousGammaF71BreakingWitness.F71AsymmetryNorm(gammas);
        double expected = Math.Sqrt(0.4);
        Assert.Equal(expected, norm, precision: 12);
    }

    // ----------------------------------------------------------------------
    // Joint-popcount γ-blindness: holds for ANY γ-list
    // ----------------------------------------------------------------------

    // The three sizes the claim names; each asymmetric list is NOT γ_l = γ_{N-1-l}.
    public static TheoryData<int, double[]> AsymmetricGammas => new()
    {
        { 4, new[] { 0.1, 0.5, 0.8, 0.3 } },
        { 4, new[] { 0.2, 0.9, 0.5, 0.4 } },
        { 5, new[] { 0.3, 0.4, 0.5, 0.6, 0.7 } },
        { 5, new[] { 0.9, 0.2, 0.5, 0.6, 0.3 } },
        { 6, new[] { 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 } },
        { 6, new[] { 0.8, 0.3, 0.6, 0.2, 0.5, 0.4 } },
    };

    public static TheoryData<int, double[]> PalindromicGammas => new()
    {
        { 4, new[] { 0.5, 0.5, 0.5, 0.5 } },
        { 4, new[] { 0.3, 0.4, 0.4, 0.3 } },
        { 5, new[] { 0.5, 0.5, 0.5, 0.5, 0.5 } },
        { 5, new[] { 0.3, 0.4, 0.5, 0.4, 0.3 } },
        { 6, new[] { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 } },
        { 6, new[] { 0.3, 0.4, 0.5, 0.5, 0.4, 0.3 } },
    };

    [Theory]
    [MemberData(nameof(AsymmetricGammas))]
    public void JointPopcountBlockDiagonal_AsymmetricGamma_IsExact(int N, double[] gammaPerSite)
    {
        // Joint-popcount sectors must still be exactly block-diagonal (U(1)×U(1) holds
        // for any γ_l ≥ 0, since per-site popcount conservation doesn't depend on γ).

        double asym = InhomogeneousGammaF71BreakingWitness.F71AsymmetryNorm(gammaPerSite);
        Assert.True(asym > 0.1, $"Test γ-list must be F71-asymmetric; got asymmetry norm = {asym:E3}.");

        var H = PauliHamiltonian.XYChain(N, J: 1.0).ToMatrix();
        var L = PauliDephasingDissipator.BuildZ(H, gammaPerSite);
        var decomp = JointPopcountSectorBuilder.Build(N);

        double offBlockFroSq = 0.0;
        var sectors = decomp.SectorRanges;
        for (int si = 0; si < sectors.Count; si++)
        {
            for (int sj = 0; sj < sectors.Count; sj++)
            {
                if (si == sj) continue;
                var rs = sectors[si];
                var cs = sectors[sj];
                for (int rOff = 0; rOff < rs.Size; rOff++)
                {
                    int rowFlat = decomp.Permutation[rs.Offset + rOff];
                    for (int cOff = 0; cOff < cs.Size; cOff++)
                    {
                        int colFlat = decomp.Permutation[cs.Offset + cOff];
                        var z = L[rowFlat, colFlat];
                        offBlockFroSq += z.Real * z.Real + z.Imaginary * z.Imaginary;
                    }
                }
            }
        }
        // Exact route: the dissipator is diagonal in the cells and H conserves both popcounts,
        // so no off-block entry is ever written; anything but 0.0 is a construction finding.
        Assert.True(offBlockFroSq == 0.0,
            $"N={N}: joint-popcount off-block Frobenius² = {offBlockFroSq:R} under asymmetric γ; expected exactly 0 (γ-blind).");
    }

    // ----------------------------------------------------------------------
    // F71 refinement: exact iff γ palindromic
    // ----------------------------------------------------------------------

    [Theory]
    [MemberData(nameof(PalindromicGammas))]
    public void F71OffBlockFrobenius_PalindromicGamma_IsZero(int N, double[] gammaPerSite)
    {
        double offBlockFro = ComputeF71OffBlockFrobenius(N, gammaPerSite);
        Assert.True(offBlockFro < 1e-10,
            $"N={N} palindromic γ: F71 off-block Frobenius = {offBlockFro:E3}; expected ~0.");
    }

    [Theory]
    [MemberData(nameof(AsymmetricGammas))]
    public void F71OffBlockFrobenius_AsymmetricGamma_IsSignificant(int N, double[] gammaPerSite)
    {
        // Off-block Frobenius (refined basis) against the on-block one, for ratio scaling.
        double offBlockFro = ComputeF71OffBlockFrobenius(N, gammaPerSite);
        double onBlockFro = ComputeF71OnBlockFrobenius(N, gammaPerSite);

        Assert.True(offBlockFro > 0.1,
            $"N={N} asymmetric γ: F71 off-block Frobenius = {offBlockFro:E3}; expected > 0.1.");
        Assert.True(offBlockFro / onBlockFro > 0.05,
            $"N={N} asymmetric γ: off/on Frobenius ratio = {offBlockFro / onBlockFro:E3}; expected > 0.05.");
    }

    [Theory]
    [MemberData(nameof(AsymmetricGammas))]
    public void F71OffBlockFrobenius_IsProportionalToAsymmetryNorm(int N, double[] gammaPerSite)
    {
        // The derivation is in the claim's doc comment (PROOF_F91 Eq. 9 summed over all 4^N
        // cells); two non-parallel δ per N make the gate see isotropy, not one direction.
        // Error model, a worst-case bound rather than a fit: each cross entry carries a few
        // ulp from the ±1/√2 rotation, and a naive sum of n = 4^N positive terms adds at most
        // n·eps relative, so |ratio/2^(N−1) − 1| ≤ c·4^N·eps with c covering the per-term and
        // the final square-root error. c = 4 gives ~3.6e-12 at N = 6; a wrong constant or a
        // non-isotropic form misses by far more.
        double ratio = ComputeF71OffBlockFrobenius(N, gammaPerSite)
                     / InhomogeneousGammaF71BreakingWitness.F71AsymmetryNorm(gammaPerSite);
        double expected = 1 << (N - 1);
        double bound = 4.0 * (1L << (2 * N)) * Math.Pow(2, -52);
        Assert.True(Math.Abs(ratio / expected - 1.0) < bound,
            $"N={N}: off-block / asymmetry norm = {ratio:R}; expected 2^(N−1) = {expected} within {bound:E1}.");
    }

    // ----------------------------------------------------------------------
    // F1 palindrome center stays at -Σγ regardless of γ asymmetry
    // ----------------------------------------------------------------------

    [Theory]
    [MemberData(nameof(AsymmetricGammas))]
    public void F1PalindromeCenter_AsymmetricGamma_StaysAtMinusSumGamma(int N, double[] gammaPerSite)
    {
        // Palindrome center = −Σγ → mirror map λ → −2Σγ − λ.
        double sumGamma = gammaPerSite.Sum();

        // Compute spectrum via the per-block path (matches full L by prior tests).
        var H = PauliHamiltonian.XYChain(N, J: 1.0).ToMatrix();
        var spectrum = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaPerSite, N);

        // Greedy nearest-neighbour matching: each eigenvalue paired with its mirror under
        // λ → −2·Σγ − λ.
        AssertSpectrumPalindromicAroundCenter(spectrum, sumGamma, tol: 1e-7, N: N);
    }

    // ----------------------------------------------------------------------
    // Helpers
    // ----------------------------------------------------------------------

    /// <summary>F71 off-block Frobenius in the refined basis: sum across joint-popcount
    /// sectors of the Frobenius² of the (F71-even × F71-odd) cross blocks. Mirror of
    /// the CLI's ComputeF71OffBlockNorm helper.</summary>
    private static double ComputeF71OffBlockFrobenius(int N, double[] gammaPerSite)
    {
        var H = PauliHamiltonian.XYChain(N, J: 1.0).ToMatrix();
        int d = 1 << N;
        var mirrorBits = F71MirrorIndexHelper.BuildHilbertMirrorLookup(N);
        int Mirror(int flat) => mirrorBits[flat / d] * d + mirrorBits[flat % d];
        double invSqrt2 = 1.0 / Math.Sqrt(2.0);

        var baseDecomp = JointPopcountSectorBuilder.Build(N);
        double offBlockFroSq = 0.0;

        foreach (var sector in baseDecomp.SectorRanges)
        {
            int size = sector.Size;
            if (size == 0) continue;

            var sectorFlat = new int[size];
            for (int k = 0; k < size; k++) sectorFlat[k] = baseDecomp.Permutation[sector.Offset + k];
            var (fixedPoints, pairs) = F71MirrorIndexHelper.FindOrbitsInSector(sectorFlat, Mirror);
            int nFix = fixedPoints.Count;
            int nPairs = pairs.Count;
            int unionSize = nFix + 2 * nPairs;
            if (unionSize == 0) continue;

            var unionFlat = new int[unionSize];
            for (int i = 0; i < nFix; i++) unionFlat[i] = fixedPoints[i];
            for (int k = 0; k < nPairs; k++) unionFlat[nFix + k] = pairs[k].S;
            for (int k = 0; k < nPairs; k++) unionFlat[nFix + nPairs + k] = pairs[k].Ps;

            var unionBlock = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaPerSite, unionFlat);

            var R = Matrix<Complex>.Build.Dense(unionSize, unionSize);
            for (int i = 0; i < nFix; i++) R[i, i] = Complex.One;
            for (int k = 0; k < nPairs; k++)
            {
                int evenCol = nFix + k;
                R[nFix + k, evenCol] = invSqrt2;
                R[nFix + nPairs + k, evenCol] = invSqrt2;
            }
            int oddOffset = nFix + nPairs;
            for (int k = 0; k < nPairs; k++)
            {
                int oddCol = oddOffset + k;
                R[nFix + k, oddCol] = invSqrt2;
                R[nFix + nPairs + k, oddCol] = -invSqrt2;
            }
            var rotated = R.Transpose() * unionBlock * R;

            int evenSize = nFix + nPairs;
            int oddSize = nPairs;
            for (int i = 0; i < evenSize; i++)
                for (int j = 0; j < oddSize; j++)
                {
                    var z = rotated[i, evenSize + j];
                    offBlockFroSq += z.Real * z.Real + z.Imaginary * z.Imaginary;
                }
            for (int i = 0; i < oddSize; i++)
                for (int j = 0; j < evenSize; j++)
                {
                    var z = rotated[evenSize + i, j];
                    offBlockFroSq += z.Real * z.Real + z.Imaginary * z.Imaginary;
                }
        }
        return Math.Sqrt(offBlockFroSq);
    }

    /// <summary>F71 on-block Frobenius: sum of Frobenius² of the (even × even) and
    /// (odd × odd) diagonal blocks in the refined basis. Used to normalize the off-block
    /// scale.</summary>
    private static double ComputeF71OnBlockFrobenius(int N, double[] gammaPerSite)
    {
        var H = PauliHamiltonian.XYChain(N, J: 1.0).ToMatrix();
        int d = 1 << N;
        var mirrorBits = F71MirrorIndexHelper.BuildHilbertMirrorLookup(N);
        int Mirror(int flat) => mirrorBits[flat / d] * d + mirrorBits[flat % d];
        double invSqrt2 = 1.0 / Math.Sqrt(2.0);

        var baseDecomp = JointPopcountSectorBuilder.Build(N);
        double onBlockFroSq = 0.0;

        foreach (var sector in baseDecomp.SectorRanges)
        {
            int size = sector.Size;
            if (size == 0) continue;

            var sectorFlat = new int[size];
            for (int k = 0; k < size; k++) sectorFlat[k] = baseDecomp.Permutation[sector.Offset + k];
            var (fixedPoints, pairs) = F71MirrorIndexHelper.FindOrbitsInSector(sectorFlat, Mirror);
            int nFix = fixedPoints.Count;
            int nPairs = pairs.Count;
            int unionSize = nFix + 2 * nPairs;
            if (unionSize == 0) continue;

            var unionFlat = new int[unionSize];
            for (int i = 0; i < nFix; i++) unionFlat[i] = fixedPoints[i];
            for (int k = 0; k < nPairs; k++) unionFlat[nFix + k] = pairs[k].S;
            for (int k = 0; k < nPairs; k++) unionFlat[nFix + nPairs + k] = pairs[k].Ps;

            var unionBlock = PerBlockLiouvillianBuilder.BuildBlockZ(H, gammaPerSite, unionFlat);

            var R = Matrix<Complex>.Build.Dense(unionSize, unionSize);
            for (int i = 0; i < nFix; i++) R[i, i] = Complex.One;
            for (int k = 0; k < nPairs; k++)
            {
                int evenCol = nFix + k;
                R[nFix + k, evenCol] = invSqrt2;
                R[nFix + nPairs + k, evenCol] = invSqrt2;
            }
            int oddOffset = nFix + nPairs;
            for (int k = 0; k < nPairs; k++)
            {
                int oddCol = oddOffset + k;
                R[nFix + k, oddCol] = invSqrt2;
                R[nFix + nPairs + k, oddCol] = -invSqrt2;
            }
            var rotated = R.Transpose() * unionBlock * R;

            int evenSize = nFix + nPairs;
            int oddSize = nPairs;
            for (int i = 0; i < evenSize; i++)
                for (int j = 0; j < evenSize; j++)
                {
                    var z = rotated[i, j];
                    onBlockFroSq += z.Real * z.Real + z.Imaginary * z.Imaginary;
                }
            for (int i = 0; i < oddSize; i++)
                for (int j = 0; j < oddSize; j++)
                {
                    var z = rotated[evenSize + i, evenSize + j];
                    onBlockFroSq += z.Real * z.Real + z.Imaginary * z.Imaginary;
                }
        }
        return Math.Sqrt(onBlockFroSq);
    }

    /// <summary>F1 palindrome verification: greedy nearest-neighbour multiset matching of
    /// the spectrum against its mirror image under λ → −2·Σγ − λ. Identical to the CLI's
    /// MultisetPalindromeOk + MultisetMatches helpers.</summary>
    private static void AssertSpectrumPalindromicAroundCenter(
        IReadOnlyList<Complex> spectrum, double sumGamma, double tol, int N)
    {
        int n = spectrum.Count;
        var mirrored = new Complex[n];
        for (int i = 0; i < n; i++)
            mirrored[i] = new Complex(-2 * sumGamma - spectrum[i].Real, spectrum[i].Imaginary);

        var taken = new bool[n];
        for (int i = 0; i < n; i++)
        {
            double bestDist = double.PositiveInfinity;
            int bestJ = -1;
            for (int j = 0; j < n; j++)
            {
                if (taken[j]) continue;
                double dist = (spectrum[i] - mirrored[j]).Magnitude;
                if (dist < bestDist) { bestDist = dist; bestJ = j; }
            }
            Assert.True(bestJ >= 0 && bestDist < tol,
                $"N={N}: palindrome match for eigenvalue {i} ({spectrum[i]}) failed; nearest mirror dist = {bestDist:E3} (tol={tol:E1}).");
            taken[bestJ] = true;
        }
    }
}
