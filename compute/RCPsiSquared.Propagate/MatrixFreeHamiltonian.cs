using System.Numerics;

namespace RCPsiSquared.Propagate;

/// <summary>
/// Matrix-free Hamiltonian for large N (N >= 14).
/// Computes [H, rho] without ever forming the d x d Hamiltonian matrix.
/// Each Heisenberg bond J*(XX + YY + ZZ) is applied via bit manipulation.
///
/// Memory: O(N_bonds) instead of O(d^2) for the Hamiltonian.
/// </summary>
public sealed class MatrixFreeHamiltonian
{
    private readonly int[] _maskA;
    private readonly int[] _maskB;
    private readonly int[] _bitA;
    private readonly int[] _bitB;
    private readonly double[] _J;
    private readonly int _nBonds;
    private readonly int _d;

    public MatrixFreeHamiltonian(int nQubits, Bond[] bonds)
    {
        _d = 1 << nQubits;
        _nBonds = bonds.Length;
        _maskA = new int[_nBonds];
        _maskB = new int[_nBonds];
        _bitA = new int[_nBonds];
        _bitB = new int[_nBonds];
        _J = new double[_nBonds];

        for (int b = 0; b < _nBonds; b++)
        {
            _bitA[b] = nQubits - 1 - bonds[b].QubitA;
            _bitB[b] = nQubits - 1 - bonds[b].QubitB;
            _maskA[b] = 1 << _bitA[b];
            _maskB[b] = 1 << _bitB[b];
            _J[b] = bonds[b].J;
        }
    }

    /// <summary>
    /// Compute [H, rho] into output array. Both arrays are column-major Complex[d*d].
    /// rho[i,j] = rho[j*d + i].
    ///
    /// For each Heisenberg bond (a,b) with coupling J:
    ///   [XX+YY+ZZ, rho][i,j] = J * ((sz_i - sz_j) * rho[i,j]
    ///                               + (1 - sz_i) * rho[i^mask, j]
    ///                               - (1 - sz_j) * rho[i, j^mask])
    /// where sz_x = 1 - 2*((bit_a(x) XOR bit_b(x)) AND 1), mask = maskA | maskB.
    /// </summary>
    public void ApplyCommutator(Complex[] rho, Complex[] output)
    {
        int d = _d;
        int nBonds = _nBonds;
        var maskA = _maskA;
        var maskB = _maskB;
        var bitA = _bitA;
        var bitB = _bitB;
        var J = _J;

        // Parallel over columns for cache-friendly column-major access
        Parallel.For(0, d, j =>
        {
            int jd = j * d;
            for (int i = 0; i < d; i++)
            {
                Complex comm = Complex.Zero;
                for (int b = 0; b < nBonds; b++)
                {
                    int mask = maskA[b] | maskB[b];
                    int iFlip = i ^ mask;
                    int jFlip = j ^ mask;

                    // sz = +1 if bits a,b have same parity, -1 if different
                    int sz_i = 1 - 2 * (((i >> bitA[b]) ^ (i >> bitB[b])) & 1);
                    int sz_j = 1 - 2 * (((j >> bitA[b]) ^ (j >> bitB[b])) & 1);

                    // ZZ part: (sz_i - sz_j) * rho[i,j]
                    int dz = sz_i - sz_j;

                    // XX+YY combined: factor is (1 - sz) which is 0 or 2
                    // Skip random memory read when factor is 0
                    if (dz != 0 || sz_i != 1 || sz_j != 1)
                    {
                        Complex val = dz * rho[jd + i];
                        if (sz_i != 1) // (1 - sz_i) = 2
                            val += 2.0 * rho[jd + iFlip];
                        if (sz_j != 1) // (1 - sz_j) = 2
                            val -= 2.0 * rho[jFlip * d + i];
                        comm += J[b] * val;
                    }
                }
                output[jd + i] = comm;
            }
        });
    }
}
