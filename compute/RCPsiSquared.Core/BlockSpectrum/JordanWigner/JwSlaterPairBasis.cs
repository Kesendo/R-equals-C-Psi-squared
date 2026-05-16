using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.JordanWigner;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.BlockSpectrum.JordanWigner;

/// <summary>JW Slater-pair basis for an arbitrary joint-popcount sector (p_c, p_r) of the
/// XY+Z-dephasing Liouvillian on N qubits. Generalises <see cref="JwBlockBasis"/> (which is
/// hardcoded to F86's c=2 (n_c=1, n_r=2) coherence block) to any (p_c, p_r) sector.
///
/// <para>The JW Slater-pair basis indexes vectors by ordered pairs (K, L) where K ⊂ {1..N}
/// is a p_c-element mode subset on the column side and L ⊂ {1..N} is a p_r-element mode
/// subset on the row side. The basis-transformation matrix element from a computational
/// basis vector |row_state⟩⟨col_state| (with popcount(row_state) = p_r, popcount(col_state)
/// = p_c) to the JW basis vector |L⟩⟨K| is the product of two Slater determinants:</para>
///
/// <code>
///   U[(row, col) flat, (K, L) α] = det[ψ_{L_i}(rowSites_j)]_{i,j=1..p_r}
///                                · det[ψ_{K_i}(colSites_j)]_{i,j=1..p_c}
/// </code>
///
/// <para>where rowSites is the sorted list of excited-bit positions in row_state and ψ_k(j)
/// is the OBC sine mode from <see cref="XyJordanWignerModes"/>.</para>
///
/// <para><b>Tier outcome: <see cref="Tier.Tier1Derived"/>.</b> Standard textbook XY-JW +
/// Slater determinant algebra. Three runtime witnesses verify the construction at
/// <see cref="Tolerance"/> = 1e-10:</para>
/// <list type="bullet">
///   <item><see cref="OrthonormalityResidual"/> = ‖U·U^† − I‖_F.</item>
///   <item><see cref="MhTotalDiagonalityResidual"/> = ‖(U^† · L_H · U)_off‖_F where L_H is
///         the Hamiltonian-only Liouvillian block (γ=0) for this sector. Free-fermion XY
///         is diagonal in Slater-pair basis with eigenvalues −i·(Σε(L) − Σε(K)).</item>
///   <item><see cref="MhTotalEigenvalueMatchResidual"/> = max_α |diag(L_H_JW)[α] −
///         (−i·(Σ_{l∈L} ε_l − Σ_{k∈K} ε_k))|.</item>
/// </list>
///
/// <para>The Z-dephasing dissipator D is NOT diagonal in JW basis (it is quartic in η
/// operators); its sparsity structure in JW basis (couples Slater pairs differing by single
/// 2-fermion swaps) is the basis for the open sparse-eig / Prosen-third-quantization
/// extensions.</para>
///
/// <para>Anchor: generalises <see cref="JwBlockBasis"/> (c=2 only) per the N=10 push of
/// <c>docs/SYMMETRY_FAMILY_INVENTORY.md</c> "requires a new symmetry not yet found".</para>
/// </summary>
public sealed class JwSlaterPairBasis : Claim
{
    public const double Tolerance = 1e-10;

    public int N { get; }
    public int PCol { get; }
    public int PRow { get; }
    public double J { get; }
    public XyJordanWignerModes Modes { get; }
    public IReadOnlyList<int[]> ColumnModeSubsets { get; }
    public IReadOnlyList<int[]> RowModeSubsets { get; }
    public ComplexMatrix U { get; }
    public ComplexMatrix Uinv { get; }
    public IReadOnlyList<int> FlatIndices { get; }
    public double OrthonormalityResidual { get; }
    public double MhTotalDiagonalityResidual { get; }
    public double MhTotalEigenvalueMatchResidual { get; }

    private JwSlaterPairBasis(
        int n, int pCol, int pRow, double j,
        XyJordanWignerModes modes,
        IReadOnlyList<int[]> colSubsets, IReadOnlyList<int[]> rowSubsets,
        ComplexMatrix u, ComplexMatrix uinv,
        IReadOnlyList<int> flatIndices,
        double orthoResidual, double diagResidual, double eigMatchResidual)
        : base($"JW Slater-pair basis for joint-popcount sector (p_c={pCol}, p_r={pRow}) at N={n}: " +
               $"orthonormal U via Slater-determinant construction; Hamiltonian-only L_H diagonal " +
               $"with eigenvalues −i·(Σε(L) − Σε(K)); witnesses {{orth={orthoResidual:G3}, " +
               $"H-diag={diagResidual:G3}, eig-match={eigMatchResidual:G3}}} all < {Tolerance:G3}.",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/F86/JordanWigner/XyJordanWignerModes.cs (sine modes ψ_k, dispersion ε_k) + " +
               "compute/RCPsiSquared.Core/BlockSpectrum/PerBlockLiouvillianBuilder.cs (per-sector L block) + " +
               "compute/RCPsiSquared.Core/F86/JordanWigner/JwBlockBasis.cs (c=2 prior art); " +
               "textbook XY-JW + Slater determinant.")
    {
        N = n; PCol = pCol; PRow = pRow; J = j;
        Modes = modes;
        ColumnModeSubsets = colSubsets;
        RowModeSubsets = rowSubsets;
        U = u; Uinv = uinv;
        FlatIndices = flatIndices;
        OrthonormalityResidual = orthoResidual;
        MhTotalDiagonalityResidual = diagResidual;
        MhTotalEigenvalueMatchResidual = eigMatchResidual;
    }

    public static JwSlaterPairBasis Build(int N, int pCol, int pRow, double J = 1.0)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be ≥ 1.");
        if (pCol < 0 || pCol > N) throw new ArgumentOutOfRangeException(nameof(pCol), pCol, $"pCol must be in [0, {N}].");
        if (pRow < 0 || pRow > N) throw new ArgumentOutOfRangeException(nameof(pRow), pRow, $"pRow must be in [0, {N}].");

        var modes = XyJordanWignerModes.Build(N, J);

        // 1. Enumerate computational basis states with the right popcount.
        var colStates = BlockBasis.PopcountStates(N, pCol);
        var rowStates = BlockBasis.PopcountStates(N, pRow);
        int mc = colStates.Count;
        int mr = rowStates.Count;
        int sectorDim = mc * mr;

        // 2. Enumerate mode subsets (lex-sorted size-p subsets of {1..N}).
        var colSubsets = EnumerateSortedSubsets(N, pCol);
        var rowSubsets = EnumerateSortedSubsets(N, pRow);
        if (colSubsets.Count * rowSubsets.Count != sectorDim)
            throw new InvalidOperationException(
                $"Subset-count mismatch: |colSubsets|·|rowSubsets| = {colSubsets.Count}·{rowSubsets.Count} " +
                $"vs sectorDim {sectorDim}.");

        // 3. Build sector flat-index list (PerBlockLiouvillianBuilder convention: flat = row·d + col).
        int d = 1 << N;
        var flatIndices = new int[sectorDim];
        int rcIdx = 0;
        for (int ri = 0; ri < mr; ri++)
            for (int ci = 0; ci < mc; ci++)
                flatIndices[rcIdx++] = (int)(rowStates[ri] * d + colStates[ci]);

        // 4. Build U[(row, col) flat, (K, L) α] = det[ψ_L(rowSites)] · det[ψ_K(colSites)].
        //    Sector flat-index order: (ri, ci) → ri*mc + ci  (row-major over col, fastest).
        //    JW basis order: (Kidx, Lidx) → Kidx*|L| + Lidx? or (Lidx, Kidx) → Lidx*|K| + Kidx?
        //    Pick: alpha = Lidx*|K| + Kidx (row L outermost matches sector flat outer row).
        var uRaw = new Complex[sectorDim, sectorDim];
        var colSites = new int[pCol];
        var rowSites = new int[pRow];
        for (int ci = 0; ci < mc; ci++)
        {
            ExtractBitsAscending(colStates[ci], N, colSites);
            for (int ri = 0; ri < mr; ri++)
            {
                ExtractBitsAscending(rowStates[ri], N, rowSites);
                int flatRow = ri * mc + ci;  // matches flatIndices[] ordering
                for (int Lidx = 0; Lidx < rowSubsets.Count; Lidx++)
                {
                    double slaterRow = SlaterDeterminant(rowSubsets[Lidx], rowSites, modes);
                    for (int Kidx = 0; Kidx < colSubsets.Count; Kidx++)
                    {
                        double slaterCol = SlaterDeterminant(colSubsets[Kidx], colSites, modes);
                        int alpha = Lidx * colSubsets.Count + Kidx;
                        uRaw[flatRow, alpha] = new Complex(slaterRow * slaterCol, 0.0);
                    }
                }
            }
        }
        var U = Matrix<Complex>.Build.DenseOfArray(uRaw);
        var Uinv = U.ConjugateTranspose();

        // 5. Orthonormality witness: ‖U · U^† − I‖_F.
        var identity = Matrix<Complex>.Build.DenseIdentity(sectorDim);
        double orthoResidual = (U * Uinv - identity).FrobeniusNorm();

        // 6. Hamiltonian-only block via PerBlockLiouvillianBuilder with γ=0.
        var H = PauliHamiltonian.XYChain(N, J).ToMatrix();
        var zeroGamma = new double[N];
        var lhBlock = PerBlockLiouvillianBuilder.BuildBlockZ(H, zeroGamma, flatIndices);

        // 7. Transform L_H to JW basis and verify diagonality + eigenvalue match.
        var lhJw = Uinv * lhBlock * U;
        double diagSqSum = 0.0;
        double maxEigDelta = 0.0;
        for (int alpha = 0; alpha < sectorDim; alpha++)
        {
            int Lidx = alpha / colSubsets.Count;
            int Kidx = alpha % colSubsets.Count;
            double sumEpsL = 0.0;
            foreach (var l in rowSubsets[Lidx]) sumEpsL += modes.Dispersion[l - 1];
            double sumEpsK = 0.0;
            foreach (var k in colSubsets[Kidx]) sumEpsK += modes.Dispersion[k - 1];
            // L_H = -i [H, ρ] gives diag entry -i·(sumEpsL - sumEpsK) on the |L⟩⟨K| basis vector
            // (consistent with the row·d + col convention applied with row-side L, col-side K).
            Complex expected = new Complex(0.0, -(sumEpsL - sumEpsK));
            Complex actual = lhJw[alpha, alpha];
            double delta = (actual - expected).Magnitude;
            if (delta > maxEigDelta) maxEigDelta = delta;
            for (int beta = 0; beta < sectorDim; beta++)
            {
                if (beta == alpha) continue;
                Complex off = lhJw[alpha, beta];
                diagSqSum += off.Real * off.Real + off.Imaginary * off.Imaginary;
            }
        }
        double diagResidual = Math.Sqrt(diagSqSum);

        if (orthoResidual > Tolerance)
            throw new InvalidOperationException(
                $"JwSlaterPairBasis (N={N}, pCol={pCol}, pRow={pRow}): orthonormality residual " +
                $"{orthoResidual:G3} exceeds tolerance {Tolerance:G3}.");
        if (diagResidual > Tolerance)
            throw new InvalidOperationException(
                $"JwSlaterPairBasis (N={N}, pCol={pCol}, pRow={pRow}): L_H off-diagonality residual " +
                $"{diagResidual:G3} exceeds tolerance {Tolerance:G3}.");
        if (maxEigDelta > Tolerance)
            throw new InvalidOperationException(
                $"JwSlaterPairBasis (N={N}, pCol={pCol}, pRow={pRow}): L_H eigenvalue match residual " +
                $"{maxEigDelta:G3} exceeds tolerance {Tolerance:G3}.");

        return new JwSlaterPairBasis(N, pCol, pRow, J, modes, colSubsets, rowSubsets,
            U, Uinv, flatIndices, orthoResidual, diagResidual, maxEigDelta);
    }

    /// <summary>Lex-sorted enumeration of size-k subsets of {1..N}, each subset returned
    /// as a sorted int[] of length k. Output count = C(N, k).</summary>
    public static IReadOnlyList<int[]> EnumerateSortedSubsets(int N, int k)
    {
        var result = new List<int[]>();
        if (k == 0) { result.Add(Array.Empty<int>()); return result; }
        var current = new int[k];
        Recurse(1, 0);
        return result;

        void Recurse(int start, int depth)
        {
            if (depth == k) { result.Add((int[])current.Clone()); return; }
            int remaining = k - depth;
            for (int v = start; v <= N - remaining + 1; v++)
            {
                current[depth] = v;
                Recurse(v + 1, depth + 1);
            }
        }
    }

    /// <summary>Extract ascending bit positions of a state. State uses big-endian convention
    /// (site 0 = MSB) matching <see cref="BlockBasis"/>; bit positions returned are 0-indexed
    /// site numbers in ascending order. Output length = popcount(state) must equal sites.Length.</summary>
    private static void ExtractBitsAscending(long state, int N, int[] sites)
    {
        int written = 0;
        for (int site = 0; site < N; site++)
        {
            // Big-endian: site 0 = MSB = bit (N-1)
            int bit = N - 1 - site;
            if (((state >> bit) & 1L) != 0)
                sites[written++] = site;
        }
        if (written != sites.Length)
            throw new InvalidOperationException(
                $"Bit-extract count {written} != expected {sites.Length} for state {state} at N={N}.");
    }

    /// <summary>Slater determinant det[ψ_{modes_i}(sites_j)]_{i,j=1..n}. For n ≤ 5 (typical
    /// at N ≤ 10) uses Leibniz expansion; falls back to LU for larger n.</summary>
    private static double SlaterDeterminant(int[] modes, int[] sites, XyJordanWignerModes psi)
    {
        int n = modes.Length;
        if (n != sites.Length)
            throw new InvalidOperationException($"Slater dim mismatch: {n} modes vs {sites.Length} sites.");
        if (n == 0) return 1.0;
        if (n == 1) return psi.SineMode(modes[0], sites[0]);
        if (n == 2)
        {
            double a = psi.SineMode(modes[0], sites[0]);
            double b = psi.SineMode(modes[0], sites[1]);
            double c = psi.SineMode(modes[1], sites[0]);
            double dd = psi.SineMode(modes[1], sites[1]);
            return a * dd - b * c;
        }

        // General Leibniz: O(n!·n) — fine for n ≤ ~7 at N ≤ 14.
        var perm = new int[n];
        for (int i = 0; i < n; i++) perm[i] = i;
        double det = 0.0;
        bool first = true;
        double sign = 1.0;
        do
        {
            double term = 1.0;
            for (int i = 0; i < n; i++) term *= psi.SineMode(modes[i], sites[perm[i]]);
            det += (first ? 1.0 : sign) * term;
            first = false;
        } while (NextPermutation(perm, out sign));
        return det;
    }

    /// <summary>In-place next-permutation in lex order; outputs the sign of the resulting
    /// permutation relative to identity. Returns false when wrapping back to identity.</summary>
    private static bool NextPermutation(int[] a, out double sign)
    {
        int n = a.Length;
        int i = n - 2;
        while (i >= 0 && a[i] >= a[i + 1]) i--;
        if (i < 0)
        {
            sign = 1.0;
            return false;
        }
        int j = n - 1;
        while (a[j] <= a[i]) j--;
        (a[i], a[j]) = (a[j], a[i]);
        Array.Reverse(a, i + 1, n - i - 1);
        // Sign: recompute from scratch (small n).
        sign = PermutationSign(a);
        return true;
    }

    private static double PermutationSign(int[] a)
    {
        int inversions = 0;
        for (int i = 0; i < a.Length; i++)
            for (int j = i + 1; j < a.Length; j++)
                if (a[i] > a[j]) inversions++;
        return (inversions & 1) == 0 ? 1.0 : -1.0;
    }

    public override string DisplayName =>
        $"JW Slater-pair basis (p_c={PCol}, p_r={PRow}, N={N}); dim {U.RowCount}";

    public override string Summary =>
        $"orth={OrthonormalityResidual:G3}, H-diag={MhTotalDiagonalityResidual:G3}, " +
        $"eig-match={MhTotalEigenvalueMatchResidual:G3} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("N", N);
            yield return InspectableNode.RealScalar("p_c", PCol);
            yield return InspectableNode.RealScalar("p_r", PRow);
            yield return InspectableNode.RealScalar("J", J);
            yield return InspectableNode.RealScalar("sector dimension", U.RowCount);
            yield return InspectableNode.RealScalar("orthonormality residual", OrthonormalityResidual, "G4");
            yield return InspectableNode.RealScalar("L_H off-diagonality residual", MhTotalDiagonalityResidual, "G4");
            yield return InspectableNode.RealScalar("L_H eigenvalue match residual", MhTotalEigenvalueMatchResidual, "G4");
        }
    }
}
