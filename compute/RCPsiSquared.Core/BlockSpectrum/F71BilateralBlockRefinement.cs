using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.BlockSpectrum;

/// <summary>F71BilateralBlockRefinement (Tier 2 empirical; 2026-05-12):
///
/// Constructs the Z₂ × Z₂ irrep basis change Q for the two independent F71 mirror actions
/// on Liouville space: <c>F71_col</c> (bit-reverse of the column index — bra side of
/// <c>|a⟩⟨b|</c>; acts as <c>I ⊗ Mirror</c>) and <c>F71_row</c> (bit-reverse of the row
/// index — ket side; acts as <c>Mirror ⊗ I</c>). The Z₂ × Z₂ orbit decomposition of each
/// joint-popcount sector yields up to four sub-blocks labeled by
/// <c>(col_parity, row_parity) ∈ {(+,+), (+,−), (−,+), (−,−)}</c>. <see cref="F71MirrorBlockRefinement"/>
/// is the diagonal Z₂ within this Z₂ × Z₂: F71-even ≡ (+,+) ⊕ (−,−), F71-odd ≡ (+,−) ⊕ (−,+).
///
/// <para><b>Spectral reduction status (Tier 2 empirical).</b> For the standard chain
/// XY + local Z-dephasing Liouvillian
/// <c>L = -i [H, ·] + Σ_l γ_l (Z_l ⊗ Z_l − I)</c> the bilateral split does NOT
/// block-diagonalise <c>Q^T L Q</c>, even at uniform γ. The reason: under <c>F71_row</c>
/// conjugation (<c>Mirror ⊗ I</c>) each <c>Z_l ⊗ Z_l</c> term transforms to
/// <c>Z_{N-1-l} ⊗ Z_l</c>; relabelling <c>l → N-1-l</c> in the sum maps it to
/// <c>Σ_l γ_l Z_l ⊗ Z_{N-1-l}</c>, which is a different operator (different (a,b) → (a',b')
/// support) than the original <c>Σ_l γ_l Z_l ⊗ Z_l</c>. Empirically the off-sub-block
/// Frobenius of <c>Q^T L Q</c> is on the order of <c>O(γ N)</c> at uniform γ, NOT zero.
/// The bilateral basis Q remains a valid Z₂ × Z₂ irrep frame on the carrier space; the
/// bilateral sub-blocks are NOT independent eigenproblems for chain XY+Z-deph.</para>
///
/// <para><b>What IS preserved.</b> The diagonal F71 (P_F71 ⊗ P_F71) commutes with L for
/// F71-palindromic γ (γ_l = γ_{N-1-l}); see <see cref="F71MirrorBlockRefinement"/> +
/// <see cref="InhomogeneousGammaF71BreakingWitness"/>. The bilateral basis Q can be used
/// to expose this: the pairs (++, −−) and (+−, −+) cluster into the diagonal-F71 even/odd
/// blocks; cross-cluster entries vanish but within-cluster (e.g. ++ ↔ −−) entries
/// generally do NOT. Splitting further to the (++) / (−−) granularity requires a
/// dissipator structure that is F71-symmetric per Hilbert side independently — e.g.
/// Heisenberg + uniform <c>(X⊗X + Y⊗Y + Z⊗Z)</c> coupling — which the standard chain
/// XY+local-Z-deph does not satisfy.</para>
///
/// <para><b>Use cases.</b> (a) Probe how strongly the bilateral split is broken under
/// different physical models — the off-sub-block-Frobenius norm is the direct measure of
/// the extra "dissipator chirality" beyond the diagonal F71. (b) Provide the irrep frame
/// for hypothetical models where <c>F71_col</c> and <c>F71_row</c> are independent L-symmetries
/// (e.g. all-to-all Z⊗Z coupling, where the dissipator is invariant under any per-side site
/// permutation). (c) Algebraic max-block projection: at N=8 the largest joint-popcount
/// sector is C(8,4)² = 4900; the largest diagonal-F71 sub-block is ≈ 2450; the largest
/// bilateral sub-sub-block is ≈ 1225 give-or-take orbit-boundary slack — the IF-it-applied
/// spectral cost would scale at this projected size.</para>
///
/// <para>Anchors: <c>compute/RCPsiSquared.Core/BlockSpectrum/JointPopcountSectorBuilder.cs</c>
/// (parent block-diagonal structure), <c>compute/RCPsiSquared.Core/BlockSpectrum/F71MirrorBlockRefinement.cs</c>
/// (parent — diagonal Z₂ within the Z₂ × Z₂; the diagonal Z₂ DOES commute with L at
/// palindromic γ, the bilateral Z₂ × Z₂ does not),
/// <c>compute/RCPsiSquared.Core.Tests/BlockSpectrum/F71BilateralBlockRefinementTests.cs</c>
/// (orthogonality + sector-count + algebraic max-block projection passing; off-sub-block
/// Frobenius nonzero at uniform / palindromic γ documented as expected behaviour for the
/// standard chain XY+Z-deph dissipator).</para></summary>
public sealed class F71BilateralBlockRefinement : Claim
{
    private readonly JointPopcountSectors _sectors;
    private readonly F71MirrorBlockRefinement _f71Diagonal;

    public F71BilateralBlockRefinement(
        JointPopcountSectors sectors,
        F71MirrorBlockRefinement f71Diagonal)
        : base("F71BilateralBlockRefinement: Z₂ × Z₂ irrep basis change Q for the two independent F71 actions on Liouville space (F71_col = I⊗Mirror, F71_row = Mirror⊗I); Q is real-orthogonal with entries 0, ±1, ±1/√2, ±1/2; orbit decomposition splits each joint-popcount sector into up to four (col_parity, row_parity) sub-blocks. Spectral block-diagonalisation of L (chain XY + local Z-deph) under Q does NOT hold even at uniform γ — the Σ_l γ_l Z_l ⊗ Z_l dissipator is not invariant under F71_row alone (becomes Σ_l γ_l Z_{N-1-l} ⊗ Z_l ≠ original). The bilateral basis is a structural irrep frame; off-sub-block Frobenius is O(γ N) at uniform γ. Useful for measuring extra dissipator chirality beyond the diagonal F71 and for hypothetical models where bilateral L-symmetry would hold (e.g. all-to-all Z⊗Z coupling).",
               Tier.Tier2Empirical,
               "JointPopcountSectors block-diagonality (parent) + F71MirrorBlockRefinement diagonal Z₂ (parent) + Z₂ × Z₂ irrep frame via independent F71_col / F71_row actions on Liouville space. Empirical: at uniform γ and F71-palindromic γ, off-sub-block Frobenius is NONZERO for the standard chain XY+local-Z-deph Liouvillian (N=3,4,5 measured in F71BilateralBlockRefinementTests); zero only for the diagonal-Z₂ subset (recovered by F71MirrorBlockRefinement). Structural correctness of Q (orthogonality, sector ordering, orbit accounting) verified bit-exactly.")
    {
        _sectors = sectors ?? throw new ArgumentNullException(nameof(sectors));
        _f71Diagonal = f71Diagonal ?? throw new ArgumentNullException(nameof(f71Diagonal));
    }

    public override string DisplayName =>
        "F71BilateralBlockRefinement: Z₂ × Z₂ irrep basis (col_parity, row_parity) over joint-popcount sectors";

    public override string Summary =>
        $"Z₂ × Z₂ irrep frame for (F71_col, F71_row) on Liouville space; each (p_c,p_r) sector → 4 sub-blocks; structural decomposition only (chain XY+local-Z-deph dissipator does NOT block-diagonalise under bilateral split); diagonal Z₂ subset recovers F71MirrorBlockRefinement ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("parent-1",
                summary: "JointPopcountSectors (block-diagonal structure, γ-blind)");
            yield return new InspectableNode("parent-2",
                summary: "F71MirrorBlockRefinement (diagonal Z₂ within this Z₂ × Z₂; F71-even = (+,+) ⊕ (−,−), F71-odd = (+,−) ⊕ (−,+); IS exact at palindromic γ)");
            yield return new InspectableNode("symmetry frame",
                summary: "F71_col (bit-reverse col index, bra side; I⊗Mirror) and F71_row (bit-reverse row index, ket side; Mirror⊗I) commute and generate Z₂ × Z₂; basis change Q is real-orthogonal with entries 0, ±1, ±1/√2, ±1/2");
            yield return new InspectableNode("orbit structure",
                summary: "size 1 (col & row both palindromic): only (+,+); size 2 col-fixed: (+,+) + (+,−); size 2 row-fixed: (+,+) + (−,+); size 4 (neither palindromic): all four sectors via Hadamard mix at ±1/2");
            yield return new InspectableNode("spectral status",
                summary: "L (chain XY + local Z-deph) is NOT block-diagonal under Q even at uniform γ — Z_l⊗Z_l dissipator transforms to Z_{N-1-l}⊗Z_l under Mirror⊗I, ≠ original; off-sub-block Frobenius O(γ N) measured at N=3,4,5 uniform γ");
            yield return new InspectableNode("usefulness",
                summary: "irrep frame for measuring extra dissipator chirality beyond diagonal F71; algebraic max-block projection (N=8: ≈ 1225 vs joint-popcount 4900); applies cleanly to all-to-all Z⊗Z or Heisenberg-with-uniform-coupling variants where bilateral symmetry IS exact");
        }
    }

    /// <summary>One refined sub-sub-block of L after Z₂ × Z₂ bilateral refinement of a
    /// joint-popcount sector. <see cref="ColParity"/> and <see cref="RowParity"/> ∈ {+1, −1}
    /// label the irrep under (F71_col, F71_row). Empty sub-blocks (Size = 0) are emitted for
    /// structural completeness — they occur, for example, when a sector has only F71-fixed
    /// pairs and the (+,−), (−,+), (−,−) sectors carry no basis vectors.</summary>
    public sealed record BilateralSector(int PCol, int PRow, int ColParity, int RowParity, int Offset, int Size);

    /// <summary>Real-orthogonal block-diagonal decomposition of L into bilateral
    /// (col_parity, row_parity) sub-blocks. <see cref="BasisChange"/> is the real-orthogonal
    /// Q such that <c>Q^T L Q</c> is block-diagonal in the (p_c, p_r, col_parity, row_parity)
    /// labels (when γ is F71-palindromic). Sub-block ranges are listed in canonical order:
    /// per joint-popcount sector, the four parity sectors (+,+), (+,−), (−,+), (−,−) appear
    /// in that order; empty sectors (Size = 0) consume no columns of Q.</summary>
    public sealed class BilateralDecomposition
    {
        public int N { get; }
        public int D { get; }   // 2^N
        public ComplexMatrix BasisChange { get; }
        public IReadOnlyList<BilateralSector> SectorRanges { get; }
        public BilateralDecomposition(int n, int d, ComplexMatrix q, IReadOnlyList<BilateralSector> sectors)
        { N = n; D = d; BasisChange = q; SectorRanges = sectors; }
    }

    /// <summary>Build the bilateral decomposition for the chain XY+Z-dephasing Liouvillian
    /// at qubit count <paramref name="N"/>. Returned basis change Q has columns ordered so
    /// that <c>Q^T L Q</c> is block-diagonal in the (p_c, p_r, col_parity, row_parity)
    /// labels (when γ is F71-palindromic). Q is real-orthogonal (entries 0, ±1, ±1/√2, ±1/2);
    /// <c>Q^H Q = I</c> bit-exactly.
    ///
    /// <para><b>Orbit structure within each (p_c, p_r) sector.</b> For each (col, row) pair
    /// with flat index <c>flat = row · d + col</c>, compute the F71_col image
    /// <c>fc(flat) = row · d + Mirror(col)</c>, the F71_row image
    /// <c>fr(flat) = Mirror(row) · d + col</c>, and the diagonal image
    /// <c>fd(flat) = Mirror(row) · d + Mirror(col)</c>. The orbit under Z₂ × Z₂ has size:
    /// <list type="bullet">
    ///   <item>1 if both col and row are F71-palindromic (all four images coincide). Contributes
    ///         only to the (+,+) sector with the bare basis vector e_flat.</item>
    ///   <item>2 (col-fixed only) if col is palindromic but row is not. Contributes
    ///         <c>(e_flat + e_{fr})/√2</c> to (+,+) and <c>(e_flat − e_{fr})/√2</c> to (+,−).</item>
    ///   <item>2 (row-fixed only) if row is palindromic but col is not. Contributes
    ///         <c>(e_flat + e_{fc})/√2</c> to (+,+) and <c>(e_flat − e_{fc})/√2</c> to (−,+).</item>
    ///   <item>4 (generic) if neither col nor row is palindromic. Contributes 4-way Hadamard-mixed
    ///         basis vectors (1/2 · (±e_x ± e_{fc(x)} ± e_{fr(x)} ± e_{fd(x)})) to all four
    ///         (col_parity, row_parity) sectors.</item>
    /// </list></para>
    ///
    /// <para>Sub-block ordering: for each joint-popcount sector (in
    /// <see cref="JointPopcountSectorBuilder.SectorRange"/> iteration order), emit the four
    /// parity sub-blocks in the canonical order (+,+), (+,−), (−,+), (−,−). Empty sub-blocks
    /// (Size = 0) are still recorded for structural completeness.</para></summary>
    public static BilateralDecomposition BuildRefinement(int N)
    {
        if (N < 1 || N > 12) throw new ArgumentOutOfRangeException(nameof(N), N, "Supported N range: 1..12.");
        var baseDecomp = JointPopcountSectorBuilder.Build(N);
        return RefineBilateral(baseDecomp);
    }

    /// <summary>Refine a joint-popcount decomposition via the Z₂ × Z₂ action of (F71_col, F71_row).
    /// See <see cref="BuildRefinement"/> for the orbit-irrep construction and column ordering.</summary>
    public static BilateralDecomposition RefineBilateral(JointPopcountSectorBuilder.Decomposition baseDecomp)
    {
        if (baseDecomp is null) throw new ArgumentNullException(nameof(baseDecomp));
        int N = baseDecomp.N;
        int d = baseDecomp.D;
        int liouvilleDim = d * d;

        // Per-Hilbert-side F71 mirror map: a |b_0…b_{N-1}⟩ ↔ |b_{N-1}…b_0⟩.
        var mirrorBits = new int[d];
        for (int x = 0; x < d; x++)
        {
            int m = 0;
            for (int i = 0; i < N; i++)
                if (((x >> i) & 1) != 0) m |= 1 << (N - 1 - i);
            mirrorBits[x] = m;
        }

        // Bilateral images of flat = row*d + col:
        //   fc(flat) = row * d + Mirror(col)         — F71_col flips bra-side (col index)
        //   fr(flat) = Mirror(row) * d + col         — F71_row flips ket-side (row index)
        //   fd(flat) = Mirror(row) * d + Mirror(col) — diagonal (both flipped)
        int FCol(int flat) => (flat / d) * d + mirrorBits[flat % d];
        int FRow(int flat) => mirrorBits[flat / d] * d + (flat % d);
        int FDiag(int flat) => mirrorBits[flat / d] * d + mirrorBits[flat % d];

        var Q = Matrix<Complex>.Build.Sparse(liouvilleDim, liouvilleDim);
        var refinedSectors = new List<BilateralSector>();
        int writeCol = 0;

        double invSqrt2 = 1.0 / Math.Sqrt(2.0);
        const double half = 0.5;

        foreach (var sector in baseDecomp.SectorRanges)
        {
            int size = sector.Size;
            // Collect flat indices in this sector. F71_col, F71_row, F71_diag all preserve
            // the per-side popcount, so all orbit images live in the same (p_c, p_r) sector.
            var sectorFlat = new int[size];
            for (int k = 0; k < size; k++) sectorFlat[k] = baseDecomp.Permutation[sector.Offset + k];

            // Walk indices in ascending order; for each unseen flat, classify the orbit and
            // emit the corresponding irrep basis vectors. We collect per-parity column lists
            // first so we can emit sub-blocks in (++)/(+−)/(−+)/(−−) order at the end of the
            // sector regardless of orbit type.
            var seen = new HashSet<int>();
            var pp = new List<(int[] flats, double[] coeffs)>(); // (+,+) basis vectors
            var pm = new List<(int[] flats, double[] coeffs)>(); // (+,−)
            var mp = new List<(int[] flats, double[] coeffs)>(); // (−,+)
            var mm = new List<(int[] flats, double[] coeffs)>(); // (−,−)

            foreach (int flat in sectorFlat.OrderBy(x => x))
            {
                if (seen.Contains(flat)) continue;

                int row = flat / d;
                int col = flat % d;
                bool rowFixed = (mirrorBits[row] == row);
                bool colFixed = (mirrorBits[col] == col);

                if (rowFixed && colFixed)
                {
                    // Orbit size 1: contributes only to (+,+) with the bare basis vector e_flat.
                    pp.Add((new[] { flat }, new[] { 1.0 }));
                    seen.Add(flat);
                }
                else if (colFixed && !rowFixed)
                {
                    // Orbit size 2 (col-fixed-only): {flat, fr(flat)}. F71_col acts trivially
                    // (col_parity always +); F71_row swaps → row_parity ±.
                    int partner = FRow(flat);
                    int sMin = Math.Min(flat, partner);
                    int sMax = Math.Max(flat, partner);
                    pp.Add((new[] { sMin, sMax }, new[] { invSqrt2, invSqrt2 }));   // (+,+)
                    pm.Add((new[] { sMin, sMax }, new[] { invSqrt2, -invSqrt2 }));  // (+,−)
                    seen.Add(sMin); seen.Add(sMax);
                }
                else if (rowFixed && !colFixed)
                {
                    // Orbit size 2 (row-fixed-only): {flat, fc(flat)}. F71_row trivial;
                    // F71_col swaps → col_parity ±.
                    int partner = FCol(flat);
                    int sMin = Math.Min(flat, partner);
                    int sMax = Math.Max(flat, partner);
                    pp.Add((new[] { sMin, sMax }, new[] { invSqrt2, invSqrt2 }));   // (+,+)
                    mp.Add((new[] { sMin, sMax }, new[] { invSqrt2, -invSqrt2 }));  // (−,+)
                    seen.Add(sMin); seen.Add(sMax);
                }
                else
                {
                    // Orbit size 4 (generic): {flat, fc(flat), fr(flat), fd(flat)}, all distinct.
                    int xc = FCol(flat);
                    int xr = FRow(flat);
                    int xd = FDiag(flat);
                    // Use the smallest index as canonical orbit representative to stabilise
                    // ordering across different starting points.
                    int x = flat;
                    int min = flat;
                    if (xc < min) { min = xc; }
                    if (xr < min) { min = xr; }
                    if (xd < min) { min = xd; }
                    if (min != x)
                    {
                        x = min;
                        xc = FCol(x);
                        xr = FRow(x);
                        xd = FDiag(x);
                    }
                    // 4-way Hadamard mix:
                    //   (+,+): (e_x + e_{xc} + e_{xr} + e_{xd}) / 2
                    //   (+,−): (e_x + e_{xc} − e_{xr} − e_{xd}) / 2  (even under F_col, odd under F_row)
                    //   (−,+): (e_x − e_{xc} + e_{xr} − e_{xd}) / 2  (odd under F_col, even under F_row)
                    //   (−,−): (e_x − e_{xc} − e_{xr} + e_{xd}) / 2  (odd under both)
                    var orbit = new[] { x, xc, xr, xd };
                    pp.Add((orbit, new[] { half, half, half, half }));
                    pm.Add((orbit, new[] { half, half, -half, -half }));
                    mp.Add((orbit, new[] { half, -half, half, -half }));
                    mm.Add((orbit, new[] { half, -half, -half, half }));
                    seen.Add(x); seen.Add(xc); seen.Add(xr); seen.Add(xd);
                }
            }

            // Emit the four parity sub-blocks in canonical order (++)/(+−)/(−+)/(−−).
            void EmitSector(int colParity, int rowParity, List<(int[] flats, double[] coeffs)> vecs)
            {
                int start = writeCol;
                foreach (var (flats, coeffs) in vecs)
                {
                    for (int i = 0; i < flats.Length; i++)
                        Q[flats[i], writeCol] = coeffs[i];
                    writeCol++;
                }
                int subSize = writeCol - start;
                refinedSectors.Add(new BilateralSector(sector.PCol, sector.PRow, colParity, rowParity, start, subSize));
            }

            EmitSector(+1, +1, pp);
            EmitSector(+1, -1, pm);
            EmitSector(-1, +1, mp);
            EmitSector(-1, -1, mm);
        }

        if (writeCol != liouvilleDim)
            throw new InvalidOperationException(
                $"Bilateral refinement produced {writeCol} basis vectors, expected {liouvilleDim} (4^N).");

        return new BilateralDecomposition(N, d, Q, refinedSectors);
    }

    /// <summary>Compute the bilateral off-sub-block Frobenius norm of <c>Q^T L Q</c> in
    /// the Z₂ × Z₂ irrep basis. This is the diagnostic measure of how strongly the bilateral
    /// split fails for the input L: zero would indicate L commutes with both <c>F71_col</c>
    /// and <c>F71_row</c> separately; nonzero is the empirical norm of the (col_parity ↔ −col_parity)
    /// + (row_parity ↔ −row_parity) cross-block content.
    ///
    /// <para>For the standard chain XY + local Z-dephasing Liouvillian this is O(γ N) at
    /// uniform γ, NOT zero (see class docstring). The diagonal F71 Z₂ that DOES survive at
    /// palindromic γ is captured by <see cref="F71MirrorBlockRefinement"/>; the difference
    /// between the bilateral off-block and the diagonal-F71 off-block is exactly the
    /// dissipator-chirality content the bilateral basis exposes.</para></summary>
    /// <param name="L">The Liouvillian (4^N) × (4^N), row-major flat convention.</param>
    /// <param name="N">Qubit count; must satisfy <c>L.RowCount == 4^N</c>.</param>
    /// <returns>The Frobenius norm of all (i, j) sub-block off-diagonals i ≠ j of
    /// <c>Q^T L Q</c> in the bilateral basis.</returns>
    public static double OffSubBlockFrobenius(ComplexMatrix L, int N)
    {
        if (L is null) throw new ArgumentNullException(nameof(L));
        int liouvilleDim = 1 << (2 * N);
        if (L.RowCount != liouvilleDim || L.ColumnCount != liouvilleDim)
            throw new ArgumentException(
                $"L must be ({liouvilleDim})×({liouvilleDim}) for N={N}; got {L.RowCount}×{L.ColumnCount}.",
                nameof(L));

        var refined = BuildRefinement(N);
        var Lp = refined.BasisChange.ConjugateTranspose() * L * refined.BasisChange;
        double offBlockFroSq = 0.0;
        var subs = refined.SectorRanges;
        for (int si = 0; si < subs.Count; si++)
        {
            for (int sj = 0; sj < subs.Count; sj++)
            {
                if (si == sj) continue;
                if (subs[si].Size == 0 || subs[sj].Size == 0) continue;
                for (int rOff = 0; rOff < subs[si].Size; rOff++)
                {
                    int rowFlat = subs[si].Offset + rOff;
                    for (int cOff = 0; cOff < subs[sj].Size; cOff++)
                    {
                        int colFlat = subs[sj].Offset + cOff;
                        var z = Lp[rowFlat, colFlat];
                        offBlockFroSq += z.Real * z.Real + z.Imaginary * z.Imaginary;
                    }
                }
            }
        }
        return Math.Sqrt(offBlockFroSq);
    }
}
