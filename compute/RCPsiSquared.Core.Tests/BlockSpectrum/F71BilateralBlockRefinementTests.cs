using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.BlockSpectrum;

/// <summary>Tests for <see cref="F71BilateralBlockRefinement"/>: Z₂ × Z₂ irrep frame for
/// the two independent F71 mirror actions on Liouville space (F71_col = I⊗Mirror,
/// F71_row = Mirror⊗I). Strictly finer in basis structure than the diagonal-Z₂
/// <see cref="F71MirrorBlockRefinement"/>.
///
/// <para><b>Spectral status (empirical, this test file is the witness):</b> the standard
/// chain XY + local Z-dephasing Liouvillian L is NOT block-diagonalised by the bilateral
/// basis Q even at uniform γ — the local Z_l ⊗ Z_l dissipator transforms to Z_{N-1-l} ⊗ Z_l
/// under Mirror ⊗ I, which is a different operator. The off-sub-block Frobenius is O(γ N)
/// at uniform γ. The diagonal-Z₂ subset (= F71MirrorBlockRefinement) IS exact at
/// palindromic γ, recovered as a sanity check.</para>
///
/// <para><b>What is tested:</b> (1) Tier metadata + null-arg guards; (2) Q is real
/// orthogonal with the prescribed entry alphabet; (3) sector counts and ordering;
/// (4) algebraic max-block projection at N=8 lands at expected ≈ 1225 vs joint-popcount
/// 4900 and diagonal-F71 ≈ 2450; (5) off-sub-block Frobenius IS NONZERO for chain XY+local-Z-deph
/// at uniform / palindromic γ (documenting that the bilateral split fails for this L);
/// (6) off-sub-block Frobenius IS ZERO for L = pure −i[H, ·] (γ=0), where the bilateral
/// symmetry is exact (positive control); (7) the diagonal-F71 subset of the bilateral
/// off-blocks IS zero at palindromic γ (independent recovery of the
/// <see cref="F71MirrorBlockRefinement"/> result through the bilateral basis).</para></summary>
public class F71BilateralBlockRefinementTests
{
    // ----------------------------------------------------------------------
    // Claim metadata
    // ----------------------------------------------------------------------

    [Fact]
    public void Tier_IsTier2Empirical()
    {
        var sectors = new JointPopcountSectors();
        var diag = new F71MirrorBlockRefinement(sectors);
        var claim = new F71BilateralBlockRefinement(sectors, diag);
        Assert.Equal(Tier.Tier2Empirical, claim.Tier);
    }

    [Fact]
    public void DisplayName_IsNonEmpty()
    {
        var sectors = new JointPopcountSectors();
        var claim = new F71BilateralBlockRefinement(sectors, new F71MirrorBlockRefinement(sectors));
        Assert.False(string.IsNullOrWhiteSpace(claim.DisplayName));
    }

    [Fact]
    public void Summary_IsNonEmpty()
    {
        var sectors = new JointPopcountSectors();
        var claim = new F71BilateralBlockRefinement(sectors, new F71MirrorBlockRefinement(sectors));
        Assert.False(string.IsNullOrWhiteSpace(claim.Summary));
    }

    [Fact]
    public void Constructor_NullSectors_Throws()
    {
        var diag = new F71MirrorBlockRefinement(new JointPopcountSectors());
        Assert.Throws<ArgumentNullException>(() => new F71BilateralBlockRefinement(null!, diag));
    }

    [Fact]
    public void Constructor_NullDiagonal_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F71BilateralBlockRefinement(new JointPopcountSectors(), null!));
    }

    [Fact]
    public void BuildRefinement_NullDecomp_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => F71BilateralBlockRefinement.RefineBilateral(null!));
    }

    // ----------------------------------------------------------------------
    // Orthogonality + structural invariants of Q
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void BasisChange_IsRealOrthogonal(int N)
    {
        var refined = F71BilateralBlockRefinement.BuildRefinement(N);
        // Q is real-orthogonal: Q^T · Q = I (entries 0, ±1, ±1/√2, ±1/2).
        var qtq = refined.BasisChange.ConjugateTranspose() * refined.BasisChange;
        int dim = refined.D * refined.D;
        var ident = Matrix<Complex>.Build.DenseIdentity(dim);
        double err = (qtq - ident).FrobeniusNorm();
        Assert.True(err < 1e-10,
            $"N={N}: Q^T·Q deviation from identity = {err:E3} (expected real-orthogonal Q).");
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void BasisChange_EntriesInExpectedAlphabet(int N)
    {
        // Q entries should be drawn from {0, ±1, ±1/√2, ±1/2}.
        var refined = F71BilateralBlockRefinement.BuildRefinement(N);
        double invSqrt2 = 1.0 / Math.Sqrt(2.0);
        var allowed = new[] { 0.0, 1.0, -1.0, invSqrt2, -invSqrt2, 0.5, -0.5 };
        int dim = refined.D * refined.D;
        for (int r = 0; r < dim; r++)
        {
            for (int c = 0; c < dim; c++)
            {
                var z = refined.BasisChange[r, c];
                Assert.True(Math.Abs(z.Imaginary) < 1e-12, $"N={N}: Q[{r},{c}] has nonzero imag = {z.Imaginary:E3}");
                bool found = allowed.Any(a => Math.Abs(z.Real - a) < 1e-12);
                Assert.True(found, $"N={N}: Q[{r},{c}] = {z.Real:E3} not in {{0, ±1, ±1/√2, ±1/2}}");
            }
        }
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void SubBlockSizes_SumTo4ToTheN(int N)
    {
        var refined = F71BilateralBlockRefinement.BuildRefinement(N);
        long total = refined.SectorRanges.Sum(s => (long)s.Size);
        Assert.Equal(1L << (2 * N), total);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void SubBlockOffsets_ContiguousNonEmpty(int N)
    {
        var refined = F71BilateralBlockRefinement.BuildRefinement(N);
        int expectedOffset = 0;
        foreach (var sub in refined.SectorRanges)
        {
            if (sub.Size == 0) continue;
            Assert.Equal(expectedOffset, sub.Offset);
            expectedOffset += sub.Size;
        }
        Assert.Equal(refined.D * refined.D, expectedOffset);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void SectorCount_Is4xJointPopcount(int N)
    {
        // Every joint-popcount sector emits 4 entries (one per parity sector); some may be
        // empty (Size=0). Total = 4 · (N+1)².
        var refined = F71BilateralBlockRefinement.BuildRefinement(N);
        Assert.Equal(4 * JointPopcountSectors.SectorCount(N), refined.SectorRanges.Count);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void SectorEntries_AreOrderedByParity(int N)
    {
        // Within each joint-popcount sector the 4 entries appear in canonical order:
        // (++)/(+−)/(−+)/(−−).
        var refined = F71BilateralBlockRefinement.BuildRefinement(N);
        Assert.Equal(0, refined.SectorRanges.Count % 4);
        for (int g = 0; g < refined.SectorRanges.Count; g += 4)
        {
            Assert.Equal((+1, +1), (refined.SectorRanges[g].ColParity, refined.SectorRanges[g].RowParity));
            Assert.Equal((+1, -1), (refined.SectorRanges[g + 1].ColParity, refined.SectorRanges[g + 1].RowParity));
            Assert.Equal((-1, +1), (refined.SectorRanges[g + 2].ColParity, refined.SectorRanges[g + 2].RowParity));
            Assert.Equal((-1, -1), (refined.SectorRanges[g + 3].ColParity, refined.SectorRanges[g + 3].RowParity));
        }
    }

    // ----------------------------------------------------------------------
    // Sector-count comparison: bilateral has strictly more non-empty sub-sectors than
    // the diagonal F71 refinement
    // ----------------------------------------------------------------------

    [Fact]
    public void SectorCount_Bilateral_StrictlyExceedsDiagonalF71_N4()
    {
        const int N = 4;
        var baseDecomp = JointPopcountSectorBuilder.Build(N);
        var diag = F71MirrorBlockRefinement.RefineWithF71(baseDecomp);
        var bilateral = F71BilateralBlockRefinement.RefineBilateral(baseDecomp);

        int diagNonEmpty = diag.SectorRanges.Count(s => s.Size > 0);
        int bilateralNonEmpty = bilateral.SectorRanges.Count(s => s.Size > 0);

        // Bilateral splits each diag sub-block further when orbit-4 entries exist; must be
        // strictly greater than diag non-empty count.
        Assert.True(bilateralNonEmpty > diagNonEmpty,
            $"N={N}: bilateral non-empty sub-sectors = {bilateralNonEmpty} not strictly > diagonal F71 non-empty = {diagNonEmpty}.");

        // Joint-popcount has 25 sectors at N=4. Diagonal F71 emits 50 entries, bilateral 100.
        Assert.Equal(25, JointPopcountSectors.SectorCount(N));
        Assert.Equal(50, diag.SectorRanges.Count);
        Assert.Equal(100, bilateral.SectorRanges.Count);

        // Bilateral non-empty count is bounded by 100 (some sectors yield no orbit-4 vectors).
        Assert.True(bilateralNonEmpty <= 100);
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void SectorCount_Bilateral_ConsistentWithJointPopcount(int N)
    {
        var bilateral = F71BilateralBlockRefinement.BuildRefinement(N);
        Assert.Equal(4 * JointPopcountSectors.SectorCount(N), bilateral.SectorRanges.Count);
        // Sub-block sizes must sum to 4^N and be non-negative.
        long sum = 0;
        foreach (var s in bilateral.SectorRanges)
        {
            Assert.True(s.Size >= 0);
            sum += s.Size;
        }
        Assert.Equal(1L << (2 * N), sum);
    }

    // ----------------------------------------------------------------------
    // Max-block factor projection at N=8 (algebraic only, no L built)
    // ----------------------------------------------------------------------

    [Fact]
    public void N8_MaxBilateralSubBlockSize_ApproxQuartersMaxSectorSize()
    {
        // Algebraic projection: enumerate orbits in the (4, 4) joint-popcount sector
        // (largest at N=8; size 4900). Without building L, count orbit sizes 1/2/4 and
        // project the bilateral sub-block sizes from the orbit counts.
        const int N = 8;
        const int d = 1 << N;
        const int mid = N / 2;

        var palindromic = new bool[d];
        for (int x = 0; x < d; x++)
        {
            int m = 0;
            for (int i = 0; i < N; i++)
                if (((x >> i) & 1) != 0) m |= 1 << (N - 1 - i);
            palindromic[x] = (m == x);
        }

        // Enumerate all (col, row) ∈ (4, 4) sector pairs and bucket by orbit type.
        int orbit1 = 0;        // both palindromic (size 1, contributes to ++ only)
        int orbit2ColFix = 0;  // col palindromic, row not (size 2; contributes to ++ and +−)
        int orbit2RowFix = 0;  // row palindromic, col not (size 2; contributes to ++ and −+)
        int orbit4 = 0;        // neither palindromic (size 4; contributes to all four)
        int total = 0;

        for (int row = 0; row < d; row++)
        {
            if (System.Numerics.BitOperations.PopCount((uint)row) != mid) continue;
            for (int col = 0; col < d; col++)
            {
                if (System.Numerics.BitOperations.PopCount((uint)col) != mid) continue;
                total++;
                bool rowFixed = palindromic[row];
                bool colFixed = palindromic[col];
                if (rowFixed && colFixed) orbit1++;
                else if (colFixed && !rowFixed) orbit2ColFix++;
                else if (rowFixed && !colFixed) orbit2RowFix++;
                else orbit4++;
            }
        }
        Assert.Equal(4900, total);
        Assert.Equal(4900, orbit1 + orbit2ColFix + orbit2RowFix + orbit4);
        // Orbit sizes 2 and 4 must be even/quarter closed.
        Assert.Equal(0, orbit2ColFix % 2);
        Assert.Equal(0, orbit2RowFix % 2);
        Assert.Equal(0, orbit4 % 4);

        // Orbit counts:
        int n1 = orbit1;
        int n2c = orbit2ColFix / 2;
        int n2r = orbit2RowFix / 2;
        int n4 = orbit4 / 4;

        // Bilateral sub-block sizes:
        //   (+,+) = n1 + n2c + n2r + n4
        //   (+,−) = n2c + n4
        //   (−,+) = n2r + n4
        //   (−,−) = n4
        int sizePP = n1 + n2c + n2r + n4;
        int sizePM = n2c + n4;
        int sizeMP = n2r + n4;
        int sizeMM = n4;
        Assert.Equal(4900, sizePP + sizePM + sizeMP + sizeMM);

        int maxBilateral = Math.Max(Math.Max(sizePP, sizePM), Math.Max(sizeMP, sizeMM));

        // Diagonal F71 sub-block sizes for comparison:
        // diag-even = orbit1 + (one half of size-2 orbits) + (two halves of size-4 orbits)
        //           = n1 + n2c + n2r + 2·n4
        // diag-odd  = (other half of size-2 orbits) + (other two halves of size-4 orbits)
        //           = n2c + n2r + 2·n4
        int diagEven = n1 + n2c + n2r + 2 * n4;
        int diagOdd = n2c + n2r + 2 * n4;
        Assert.Equal(4900, diagEven + diagOdd);
        int maxDiag = Math.Max(diagEven, diagOdd);

        Assert.True(maxBilateral <= maxDiag,
            $"N=8 (4,4) sector: bilateral max sub-block = {maxBilateral} should be ≤ diagonal F71 max sub-block = {maxDiag}.");

        // Spec target: bilateral max ≈ 1225 (= 4900 / 4) within orbit-boundary slack.
        // Slack is dominated by n1+n2c+n2r (entries that don't get fully quartered).
        int targetMax = 4900 / 4;  // 1225
        int slack = n1 + n2c + n2r;
        Assert.True(maxBilateral <= targetMax + slack,
            $"N=8 (4,4) sector: bilateral max sub-block = {maxBilateral} exceeds target ≈ {targetMax} by more than slack = {slack} (n1={n1}, n2c={n2c}, n2r={n2r}, n4={n4}).");
        Assert.True(maxBilateral >= targetMax - slack,
            $"N=8 (4,4) sector: bilateral max sub-block = {maxBilateral} below target ≈ {targetMax} by more than slack = {slack}.");

        // The slack must be ≪ 4900 (boundary effect, not bulk).
        Assert.True(slack * 4 < 4900,
            $"N=8 (4,4) sector: orbit-boundary slack {slack} should be ≪ 4900/4; n1+n2c+n2r={slack}, n4={n4}.");
    }

    // ----------------------------------------------------------------------
    // Off-sub-block Frobenius: documenting empirical behaviour for chain XY+local-Z-deph
    //
    // Per the class docstring, the bilateral Z₂ × Z₂ split does NOT block-diagonalise L
    // for the standard chain XY + local Z-dephasing dissipator even at uniform γ. The
    // dissipator term Σ_l γ_l Z_l ⊗ Z_l, under conjugation by Mirror ⊗ I, transforms to
    // Σ_l γ_l Z_{N-1-l} ⊗ Z_l ≠ original. The off-sub-block Frobenius scales as O(γ N).
    // The diagonal Z₂ subset (P_F71 ⊗ P_F71) IS exact at palindromic γ — that is the
    // F71MirrorBlockRefinement result, recoverable through the bilateral basis as a
    // sanity check.
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void OffSubBlockFrobenius_IsNonZero_UniformGamma_ChainXYLocalZDeph(int N)
    {
        // Document empirical breaking: bilateral split fails at uniform γ for chain XY +
        // local Z-deph because the local Z_l ⊗ Z_l dissipator is not invariant under
        // F71_row alone (Mirror ⊗ I).
        const double J = 1.0;
        const double gamma = 0.5;
        var L = BuildXYZDephasingL(N, J, Enumerable.Repeat(gamma, N).ToArray());

        double offBlockFro = F71BilateralBlockRefinement.OffSubBlockFrobenius(L, N);
        Assert.True(offBlockFro > 1.0,
            $"N={N}: bilateral off-sub-block Frobenius at uniform γ = {offBlockFro:E3} (expected > 1.0 — local Z_l ⊗ Z_l dissipator is not bilateral-symmetric; per class docstring the bilateral split fails at uniform γ).");
    }

    [Fact]
    public void OffSubBlockFrobenius_IsNonZero_PalindromicGamma_N4_ChainXYLocalZDeph()
    {
        // Palindromic γ (γ_l = γ_{N-1-l}) suffices for the diagonal F71 refinement to be
        // exact, but NOT for the bilateral Z₂ × Z₂ split (which would require dissipator-side
        // bilateral invariance, not just γ-symmetry).
        const int N = 4;
        var gammaPalindromic = new[] { 0.3, 0.5, 0.5, 0.3 };
        // Sanity: γ-asymmetry norm = 0 (palindromic).
        double asym = InhomogeneousGammaF71BreakingWitness.F71AsymmetryNorm(gammaPalindromic);
        Assert.True(asym < 1e-12,
            $"Sanity: γ-asymmetry norm = {asym:E3} should be ≈ 0 for palindromic γ.");

        var L = BuildXYZDephasingL(N, J: 1.0, gammaPalindromic);
        double offBlockFro = F71BilateralBlockRefinement.OffSubBlockFrobenius(L, N);
        Assert.True(offBlockFro > 1.0,
            $"N={N}: bilateral off-sub-block Frobenius at F71-palindromic γ {string.Join(",", gammaPalindromic)} = {offBlockFro:E3} (expected > 1.0 — γ-palindromy alone does not make the local Z_l ⊗ Z_l dissipator bilateral-symmetric).");
    }

    [Fact]
    public void OffSubBlockFrobenius_IsNonZero_MonotonicGamma_N4_ChainXYLocalZDeph()
    {
        // F71-asymmetric γ: γ_l ≠ γ_{N-1-l}. Both diagonal F71 AND bilateral split leak.
        const int N = 4;
        var gammaMonotonic = new[] { 0.2, 0.4, 0.5, 0.7 };
        double asym = InhomogeneousGammaF71BreakingWitness.F71AsymmetryNorm(gammaMonotonic);
        Assert.True(asym > 0.4, $"Sanity: γ-asymmetry norm should be substantial for {string.Join(",", gammaMonotonic)}; got {asym:E3}.");

        var L = BuildXYZDephasingL(N, J: 1.0, gammaMonotonic);
        double offBlockFro = F71BilateralBlockRefinement.OffSubBlockFrobenius(L, N);
        Assert.True(offBlockFro > 0.1,
            $"N={N}: bilateral off-sub-block Frobenius at monotonic non-palindromic γ {string.Join(",", gammaMonotonic)} = {offBlockFro:E3} (expected > 0.1 — bilateral split also leaks under γ-asymmetry, in addition to dissipator chirality).");
    }

    // ----------------------------------------------------------------------
    // Positive control: at γ = 0 (pure unitary L = -i[H,·]) the bilateral split IS exact.
    // The Hamiltonian L is bilateral-symmetric: each Hilbert side commutes with the
    // chain XY Hamiltonian's F71 mirror.
    // ----------------------------------------------------------------------

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void OffSubBlockFrobenius_IsZero_ZeroGamma_PureChainXYHamiltonian(int N)
    {
        // L_H = -i [H, ·] = -i (H ⊗ I − I ⊗ H^T). Both H ⊗ I and I ⊗ H^T are bilateral-
        // symmetric for chain XY: H is F71-symmetric, so both Hilbert-side independent
        // F71 mirrors commute with L_H.
        var gammaZero = Enumerable.Repeat(0.0, N).ToArray();
        var L = BuildXYZDephasingL(N, J: 1.0, gammaZero);
        double offBlockFro = F71BilateralBlockRefinement.OffSubBlockFrobenius(L, N);
        Assert.True(offBlockFro < 1e-10,
            $"N={N}: bilateral off-sub-block Frobenius at γ=0 (pure -i[H, ·]) = {offBlockFro:E3} (expected ~0 — chain XY Hamiltonian is F71-symmetric on each Hilbert side independently).");
    }

    // ----------------------------------------------------------------------
    // Sanity: the diagonal-F71 subset of the bilateral basis still gives exact
    // block-diagonality at palindromic γ (recovering F71MirrorBlockRefinement through
    // the bilateral basis). Concretely: the (++) ⊕ (−−) and (+−) ⊕ (−+) clusters
    // correspond to F71-even and F71-odd of the diagonal Z₂; sums of cross-cluster
    // entries vanish at palindromic γ even though within-cluster entries do not.
    // ----------------------------------------------------------------------

    [Fact]
    public void DiagonalF71Cluster_OffBlockFrobenius_IsZero_PalindromicGamma_N4()
    {
        const int N = 4;
        var gammaPalindromic = new[] { 0.3, 0.5, 0.5, 0.3 };
        var L = BuildXYZDephasingL(N, J: 1.0, gammaPalindromic);

        var refined = F71BilateralBlockRefinement.BuildRefinement(N);
        var Lp = refined.BasisChange.ConjugateTranspose() * L * refined.BasisChange;

        // Cluster sub-blocks into diagonal-F71 even (++ ⊕ −−) and odd (+− ⊕ −+) groups.
        // Map each bilateral sub-block to a diagonal-F71 cluster index, then accumulate
        // off-cluster Frobenius (entries between EVEN and ODD clusters).
        // Within a joint-popcount sector, sub-blocks appear in (++)/(+−)/(−+)/(−−) order.
        var clusterOfSubblock = new int[refined.SectorRanges.Count];
        for (int i = 0; i < refined.SectorRanges.Count; i++)
        {
            var s = refined.SectorRanges[i];
            int diagParity = s.ColParity * s.RowParity;   // +1 = diagonal F71-even cluster, −1 = odd
            clusterOfSubblock[i] = diagParity > 0 ? 0 : 1;
        }

        double offClusterFroSq = 0.0;
        var subs = refined.SectorRanges;
        for (int si = 0; si < subs.Count; si++)
        {
            for (int sj = 0; sj < subs.Count; sj++)
            {
                if (si == sj) continue;
                if (subs[si].Size == 0 || subs[sj].Size == 0) continue;
                if (clusterOfSubblock[si] == clusterOfSubblock[sj]) continue;
                for (int rOff = 0; rOff < subs[si].Size; rOff++)
                {
                    int rowFlat = subs[si].Offset + rOff;
                    for (int cOff = 0; cOff < subs[sj].Size; cOff++)
                    {
                        int colFlat = subs[sj].Offset + cOff;
                        var z = Lp[rowFlat, colFlat];
                        offClusterFroSq += z.Real * z.Real + z.Imaginary * z.Imaginary;
                    }
                }
            }
        }
        double offClusterFro = Math.Sqrt(offClusterFroSq);
        Assert.True(offClusterFro < 1e-10,
            $"N={N}: at palindromic γ, off-CLUSTER Frobenius (between F71-even cluster {{++,−−}} and F71-odd cluster {{+−,−+}}) = {offClusterFro:E3} (expected ~0 — recovers F71MirrorBlockRefinement diagonal Z₂ through the bilateral basis).");
    }

    // ----------------------------------------------------------------------
    // OffSubBlockFrobenius helper: argument validation
    // ----------------------------------------------------------------------

    [Fact]
    public void OffSubBlockFrobenius_NullL_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => F71BilateralBlockRefinement.OffSubBlockFrobenius(null!, 3));
    }

    [Fact]
    public void OffSubBlockFrobenius_WrongDimension_Throws()
    {
        var wrong = Matrix<Complex>.Build.DenseIdentity(16);
        Assert.Throws<ArgumentException>(() => F71BilateralBlockRefinement.OffSubBlockFrobenius(wrong, 3));
    }

    // ----------------------------------------------------------------------
    // Helpers
    // ----------------------------------------------------------------------

    private static ComplexMatrix BuildXYZDephasingL(int N, double J, double[] gammaPerSite)
    {
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        return PauliDephasingDissipator.BuildZ(H, gammaPerSite);
    }
}
