using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.BlockSpectrum;

/// <summary>Materialises the RAW joint-popcount (p, q̃) block of the XY (Δ=0, free-fermion)
/// chain Liouvillian as a dense complex matrix, plus the full L used to validate it.
///
/// <para><b>Convention.</b> <c>L(q)/γ = −i·q·[Ĥ,·] + D[·]</c>, with
/// <c>Ĥ = Σ_bonds (X_i X_{i+1} + Y_i Y_{i+1})</c> the UNIT-J Pauli XY chain (Δ = 0, NO ZZ
/// term), open chain, <c>q = J/γ</c> (may be complex), <c>γ = 1</c> per site, and <c>D</c> the
/// single-site Z-dephasing dissipator. The q is injected into the commutator by feeding
/// <c>q·Ĥ</c> as the Hamiltonian to the block builder, since <c>−i·q·[Ĥ,·] = −i·[q·Ĥ,·]</c>.
/// The unit-J Ĥ is obtained from <see cref="PauliHamiltonian.XYChain(int, double)"/> with
/// <c>J = 2</c> (that builder carries a J/2 prefactor, so J = 2 gives coefficient 1 on each
/// XX and YY term); it has no ZZ term, so Δ = 0 is automatic.</para>
///
/// <para><b>RAW footing (load-bearing).</b> <see cref="Build"/> returns the genuine principal
/// submatrix of the full L on the sector's flat indices, with NO spatial-symmetry reduction
/// and NO ×2 clearing. It is therefore a DIFFERENT object from the S₂-symmetric-reduced
/// blocks in the F89 path-k monodromy scout; do not cross-validate against those. The only
/// footing-independent validation is the full-L decomposition: L is block-diagonal in the
/// joint-popcount basis, so the multiset union of every (p, q̃) block's spectrum equals the
/// full-L spectrum (guarded by <c>SectorBlockTests.AllSectorBlocks_UnionToFullLiouvillianSpectrum_At_N4</c>).</para>
///
/// <para><b>Reuse.</b> Pure composition of the trusted primitives: sector flat indices from
/// <see cref="JointPopcountSectorBuilder"/>, the computational-basis L block from
/// <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/> (whose D diagonal is
/// <c>Σ_l (−2γ_l)·[bit_l(row) ≠ bit_l(col)]</c> and whose Hamiltonian term already carries the
/// <c>−i</c> commutator sign), and the XY Hamiltonian from <see cref="PauliHamiltonian.XYChain(int, double)"/>.
/// Nothing about the Liouvillian element formula is re-derived here.</para>
///
/// <para><b>Memory.</b> A single (p, q̃) block costs <c>[C(N,p)·C(N,q̃)]² × 16</c> bytes; the full
/// L costs <c>(4^N)² × 16</c> bytes (~256 MB at N=6, ~4 GB at N=7). This materialiser targets
/// small N (N ≤ 6 for the full L); larger N is served by the sector-only per-block paths in
/// <see cref="LiouvillianSectorSweep"/> / <see cref="LiouvillianBlockSpectrum"/>.</para></summary>
public static class SectorBlock
{
    /// <summary>The full L(q)/γ over all 4^N Liouville-space indices in natural row-major order
    /// (<c>flat = row·2^N + col</c>), built via the same trusted
    /// <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/> as the sector blocks. Its
    /// eigenvalues are the full spectrum against which the sector-block union is validated.</summary>
    /// <param name="N">Number of sites (N ≥ 1). Kept small: this materialises a 4^N × 4^N dense matrix.</param>
    /// <param name="q">Complex coupling q = J/γ injected as H = q·Ĥ.</param>
    public static ComplexMatrix BuildFullLiouvillian(int N, Complex q)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be ≥ 1.");
        if (2 * N > 30)
            throw new ArgumentOutOfRangeException(nameof(N), N,
                "BuildFullLiouvillian materialises a 4^N × 4^N dense matrix; N is capped where the flat " +
                "index count 4^N would overflow int32 (and the memory is infeasible far earlier).");

        var H = BuildQHamiltonian(N, q);
        var gamma = UnitGamma(N);
        int fullDim = 1 << (2 * N);   // 4^N
        var flatIndices = new int[fullDim];
        for (int f = 0; f < fullDim; f++) flatIndices[f] = f;
        return PerBlockLiouvillianBuilder.BuildBlockZ(H, gamma, flatIndices);
    }

    /// <summary>The RAW (p, q̃) block: <c>B[i,j] = L(q)[f_i, f_j]</c> where the f's are the
    /// Liouville-space flat indices of the joint-popcount sector (PCol = p, PRow = q̃). Built by
    /// feeding those flat indices and H = q·Ĥ, γ_per_site = 1 to
    /// <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/>.</summary>
    /// <param name="N">Number of sites (N ≥ 1).</param>
    /// <param name="p">Column popcount (0 ≤ p ≤ N).</param>
    /// <param name="qTilde">Row popcount (0 ≤ q̃ ≤ N).</param>
    /// <param name="q">Complex coupling q = J/γ injected as H = q·Ĥ.</param>
    /// <returns>Dense block of size C(N,p)·C(N,q̃).</returns>
    public static Complex[,] Build(int N, int p, int qTilde, Complex q)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be ≥ 1.");
        if (p < 0 || p > N) throw new ArgumentOutOfRangeException(nameof(p), p, $"p must be in [0, {N}].");
        if (qTilde < 0 || qTilde > N) throw new ArgumentOutOfRangeException(nameof(qTilde), qTilde, $"q̃ must be in [0, {N}].");

        var H = BuildQHamiltonian(N, q);
        var gamma = UnitGamma(N);
        var flatIndices = SectorFlatIndices(N, p, qTilde);
        return PerBlockLiouvillianBuilder.BuildBlockZ(H, gamma, flatIndices).ToArray();
    }

    /// <summary>The Liouville-space flat indices of the (PCol = p, PRow = q̃) joint-popcount
    /// sector, read out of <see cref="JointPopcountSectorBuilder"/>'s permutation.</summary>
    public static int[] SectorFlatIndices(int N, int p, int qTilde)
    {
        var decomp = JointPopcountSectorBuilder.Build(N);
        JointPopcountSectorBuilder.SectorRange? range = null;
        foreach (var r in decomp.SectorRanges)
            if (r.PCol == p && r.PRow == qTilde) { range = r; break; }
        if (range is null)
            throw new ArgumentException($"No joint-popcount sector (p={p}, q̃={qTilde}) at N={N}.");

        var flat = new int[range.Size];
        for (int k = 0; k < range.Size; k++) flat[k] = decomp.Permutation[range.Offset + k];
        return flat;
    }

    /// <summary>Unit-J Pauli XY chain scaled by q: <c>q·Ĥ</c>, Ĥ = Σ_bonds (XX + YY), Δ = 0,
    /// open chain. XYChain(N, 2.0) yields coefficient 1 on each XX/YY (its J/2 prefactor), and
    /// scalar-multiplying by the complex q injects q into the commutator.</summary>
    private static ComplexMatrix BuildQHamiltonian(int N, Complex q) =>
        PauliHamiltonian.XYChain(N, 2.0).ToMatrix().Multiply(q);

    private static double[] UnitGamma(int N)
    {
        var gamma = new double[N];
        for (int l = 0; l < N; l++) gamma[l] = 1.0;
        return gamma;
    }
}
