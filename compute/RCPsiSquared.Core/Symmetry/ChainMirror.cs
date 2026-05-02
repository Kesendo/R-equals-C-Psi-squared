using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F71 chain-mirror spatial symmetry: R |b₀ b₁ … b_{N-1}⟩ = |b_{N-1} … b₁ b₀⟩.
///
/// R is real, symmetric, and involutory (R² = I). It commutes with the uniform XY/Heisenberg
/// chain Hamiltonian and Z-dephasing dissipator, hence with the full Liouvillian L.
/// Proven in docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md.
///
/// Empirically (EQ-024 2026-04-28): F71-eigenstate receivers are capacity-optimal at odd N,
/// capacity-suboptimal at even N. See experiments/J_BLIND_RECEIVER_CLASSES.md.
/// </summary>
public static class ChainMirror
{
    /// <summary>R as a 2^N × 2^N bit-reversal permutation matrix.</summary>
    public static ComplexMatrix Build(int N)
    {
        int d = 1 << N;
        var R = Matrix<Complex>.Build.Sparse(d, d);
        for (int i = 0; i < d; i++)
        {
            int rev = 0;
            int x = i;
            for (int _ = 0; _ < N; _++)
            {
                rev = (rev << 1) | (x & 1);
                x >>= 1;
            }
            R[rev, i] = Complex.One;
        }
        return R;
    }

    /// <summary>Symmetric projector P_sym = (I + R) / 2 onto the F71 +1 eigenspace.</summary>
    public static ComplexMatrix SymmetricProjector(int N)
    {
        int d = 1 << N;
        return (Matrix<Complex>.Build.SparseIdentity(d) + Build(N)) / 2.0;
    }

    /// <summary>Antisymmetric projector P_anti = (I − R) / 2 onto the F71 −1 eigenspace.</summary>
    public static ComplexMatrix AntisymmetricProjector(int N)
    {
        int d = 1 << N;
        return (Matrix<Complex>.Build.SparseIdentity(d) - Build(N)) / 2.0;
    }

    /// <summary>Classify a state vector ψ as F71-eigenstate.
    /// Returns +1, −1, or null (F71-mixed/breaking).</summary>
    public static int? EigenstateClass(ComplexVector psi, double tolerance = 1e-6)
    {
        int N = (int)Math.Round(Math.Log2(psi.Count));
        if ((1 << N) != psi.Count) throw new ArgumentException($"psi length {psi.Count} is not a power of 2");
        var R = Build(N);
        var Rpsi = R * psi;
        double normSq = psi.ConjugateDotProduct(psi).Real;
        if (normSq < tolerance) throw new ArgumentException("psi has zero norm");
        var overlap = psi.ConjugateDotProduct(Rpsi) / normSq;
        if (Math.Abs(overlap.Imaginary) > tolerance) return null;
        if (Math.Abs(overlap.Real - 1.0) < tolerance) return +1;
        if (Math.Abs(overlap.Real + 1.0) < tolerance) return -1;
        return null;
    }

    /// <summary>Counts of sym/asym bond-input subspace dimensions, without building the basis vectors.
    /// Cheap for callers that only need <c>Symmetric</c>/<c>Antisymmetric</c>.Length.</summary>
    public static (int Symmetric, int Antisymmetric) BondMirrorCounts(int N)
    {
        int nBonds = N - 1;
        // For chain-bond mirror b ↔ N-2-b: each (b, N-2-b) pair contributes one sym + one asym;
        // a fixed point (b == N-2-b, exists iff nBonds odd) contributes one sym only.
        int pairs = nBonds / 2;
        int hasSelfMirror = nBonds % 2;
        return (pairs + hasSelfMirror, pairs);
    }

    /// <summary>Bond-mirror basis: orthonormal sym/asym vectors on the (N-1)-dim bond-input space.
    /// R̄: J_b → J_{N-2-b} acts on R^{N-1}. Returns (sym_basis, asym_basis) row-stacked.
    ///
    /// Dimensional structure (verified at N=3..8):
    ///   odd N: balanced k+k split
    ///   even N: unbalanced (k+1)+k split (one self-mirror bond)
    /// </summary>
    public static (double[][] Symmetric, double[][] Antisymmetric) BondMirrorBasis(int N)
    {
        int nBonds = N - 1;
        var sym = new List<double[]>();
        var asym = new List<double[]>();
        var used = new bool[nBonds];
        double inv = 1.0 / Math.Sqrt(2.0);
        for (int b = 0; b < nBonds; b++)
        {
            if (used[b]) continue;
            int bMirror = nBonds - 1 - b;
            if (b == bMirror)
            {
                var v = new double[nBonds];
                v[b] = 1.0;
                sym.Add(v);
                used[b] = true;
            }
            else
            {
                var vS = new double[nBonds]; vS[b] = inv; vS[bMirror] = inv;
                var vA = new double[nBonds]; vA[b] = inv; vA[bMirror] = -inv;
                sym.Add(vS);
                asym.Add(vA);
                used[b] = true;
                used[bMirror] = true;
            }
        }
        return (sym.ToArray(), asym.ToArray());
    }
}
