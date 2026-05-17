using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.States;

/// <summary>Partial trace of a multi-qubit density matrix: keep specified qubits,
/// trace out the rest.
///
/// <para>Convention: site 0 = most-significant bit (big-endian), matching
/// <see cref="Pauli.PauliString.SiteOp"/> and the rest of <c>Core</c>. The result
/// is a 2^|keep| × 2^|keep| reduced density matrix.</para>
///
/// <para>Algorithm: precompute (kept_index, traced_bits) per full basis index
/// for an O(d²) sum. Faster than the generic tensor-reshape approach and
/// avoids the numpy-style axis bookkeeping that often introduces sign bugs.</para>
///
/// <para>This is the typed-knowledge-layer companion of the matrix-free
/// partial trace in <c>RCPsiSquared.Propagate.DensityMatrixTools.PartialTrace</c>
/// (which retains the same algorithm and convention for use in the
/// propagation engine). Both should give bit-exact identical output on the
/// same input.</para>
/// </summary>
public static class PartialTrace
{
    /// <summary>Reduce <paramref name="rho"/> (a 2^<paramref name="nQubits"/> × 2^<paramref name="nQubits"/>
    /// density matrix) by tracing out every qubit not in <paramref name="keep"/>.
    /// Returns a 2^|keep| × 2^|keep| dense matrix.</summary>
    /// <param name="rho">The full N-qubit density matrix (must be 2^N × 2^N).</param>
    /// <param name="nQubits">The number of qubits N.</param>
    /// <param name="keep">Sites to keep (0-indexed; in [0, N − 1]; order does
    /// not matter for the trace itself but determines the result's basis
    /// ordering: keep = [0, 2] gives the reduced state on qubits 0, 2 with
    /// site 0 as the most-significant bit).</param>
    public static ComplexMatrix Of(ComplexMatrix rho, int nQubits, IReadOnlyList<int> keep)
    {
        if (rho is null) throw new ArgumentNullException(nameof(rho));
        if (keep is null) throw new ArgumentNullException(nameof(keep));
        int dim = 1 << nQubits;
        if (rho.RowCount != dim || rho.ColumnCount != dim)
            throw new ArgumentException(
                $"rho must be 2^{nQubits} × 2^{nQubits} = {dim} × {dim}; got {rho.RowCount} × {rho.ColumnCount}.",
                nameof(rho));
        foreach (int k in keep)
            if (k < 0 || k >= nQubits)
                throw new ArgumentOutOfRangeException(nameof(keep),
                    $"keep site {k} must be in [0, {nQubits - 1}].");
        if (keep.Distinct().Count() != keep.Count)
            throw new ArgumentException("keep sites must be distinct.", nameof(keep));

        int nKeep = keep.Count;
        int dKeep = 1 << nKeep;
        var traced = Enumerable.Range(0, nQubits).Except(keep).ToArray();
        int nTraced = traced.Length;

        var keptIdx = new int[dim];
        var tracedBits = new int[dim];
        for (int i = 0; i < dim; i++)
        {
            int ki = 0;
            for (int m = 0; m < nKeep; m++)
                ki |= ((i >> (nQubits - 1 - keep[m])) & 1) << (nKeep - 1 - m);
            keptIdx[i] = ki;

            int tb = 0;
            for (int m = 0; m < nTraced; m++)
                tb |= ((i >> (nQubits - 1 - traced[m])) & 1) << (nTraced - 1 - m);
            tracedBits[i] = tb;
        }

        var result = DenseMatrix.Create(dKeep, dKeep, Complex.Zero);
        var dense = rho as DenseMatrix ?? DenseMatrix.OfMatrix(rho);
        var rhoVals = dense.Values;

        for (int i = 0; i < dim; i++)
        {
            int ki = keptIdx[i];
            int tbi = tracedBits[i];
            for (int j = 0; j < dim; j++)
            {
                if (tracedBits[j] != tbi) continue;
                result[ki, keptIdx[j]] += rhoVals[j * dim + i];  // column-major
            }
        }
        return result;
    }
}
