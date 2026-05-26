using System.Collections.Generic;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Polarity;

/// <summary>F112 non-Hermitian extension verifier: enumerates the Pauli-string basis
/// at length N and checks the open identity Im⟨L_{σ_α,-i}, L_{σ_β,-i}⟩ = 0 for every
/// pair. If all pairs vanish bit-exact (within tolerance), F112 non-Hermitian extension
/// is Tier1Derived at this N by bilinearity + basis spanning.
///
/// <para>C# port of <c>simulations/_f112_open_identity_basis_enum.py</c>. Existing Python
/// anchor: N=2, 3, 4 bit-exact (35,112 distinct ordered pairs total). This primitive
/// extends to N=5 (524,800 pairs, ~16 GB working memory).</para>
///
/// <para>Per pair cost: one Frobenius inner product over a 4^N × 4^N complex matrix.
/// Per-α pre-cache cost: 4^N L matrices each of size (4^N)² complex = 16·4^(2N) bytes.
/// At N=5: 1024 matrices × 16 MB each = 16 GB total cache.</para></summary>
public static class F112NonHermitianBasisEnumeration
{
    /// <summary>Build L_H = -i[H, ·] in the Pauli basis. Mirrors the Python
    /// <c>build_L_H_pauli</c>: L_vec = -i (H ⊗ I − I ⊗ H^T) is the d² × d² vec-basis
    /// commutator superoperator; rotated into the 4^N-dim Pauli basis via
    /// T^† L_vec T / 2^N where T = <see cref="PauliBasis.VecToPauliBasisTransform"/>.
    /// Returns a 4^N × 4^N complex matrix.</summary>
    public static ComplexMatrix BuildLHInPauliBasis(ComplexMatrix H, int N)
    {
        if (H is null) throw new ArgumentNullException(nameof(H));
        int d = 1 << N;
        if (H.RowCount != d || H.ColumnCount != d)
            throw new ArgumentException($"H must be {d}x{d} (2^N x 2^N); got {H.RowCount}x{H.ColumnCount}", nameof(H));

        var Id = Matrix<Complex>.Build.DenseIdentity(d);
        // L_vec = -i (H ⊗ I − I ⊗ H^T)
        var lVec = -Complex.ImaginaryOne * (H.KroneckerProduct(Id) - Id.KroneckerProduct(H.Transpose()));
        var T = PauliBasis.VecToPauliBasisTransform(N);
        return T.ConjugateTranspose() * lVec * T / (1 << N); // 2^N
    }

    /// <summary>Project a 4^N × 4^N operator-space matrix M onto the eigenspace of the
    /// Π-conjugation map X ↦ Π X Π⁻¹ at the requested eigenvalue λ ∈ {+1, -1, +i, -i}.
    /// Π is order-4 on operator space (Π⁴ = I) so the standard idempotent projector is
    /// P_λ(M) = (1/4) Σ_{k=0..3} λ^{-k} Π^k M Π^{-k}.</summary>
    public static ComplexMatrix ProjectOntoPiEigenspace(ComplexMatrix M, ComplexMatrix pi, Complex targetEigenvalue)
    {
        if (M is null) throw new ArgumentNullException(nameof(M));
        if (pi is null) throw new ArgumentNullException(nameof(pi));

        var piInv = pi.ConjugateTranspose();
        var result = Matrix<Complex>.Build.Dense(M.RowCount, M.ColumnCount);
        var curPi = Matrix<Complex>.Build.DenseIdentity(pi.RowCount);
        var curPiInv = Matrix<Complex>.Build.DenseIdentity(pi.RowCount);
        Complex lambdaToK = Complex.One;
        for (int k = 0; k < 4; k++)
        {
            var coef = (1.0 / lambdaToK) / 4.0;
            result = result + coef * (curPi * M * curPiInv);
            curPi = curPi * pi;
            curPiInv = curPiInv * piInv;
            lambdaToK = lambdaToK * targetEigenvalue;
        }
        return result;
    }

    /// <summary>Frobenius inner product ⟨A, B⟩ = Σ A[i,j]* · B[i,j] = Tr(A† B).
    ///
    /// <para>Uses the underlying dense column-major <c>Complex[]</c> storage via
    /// <see cref="DenseMatrix.Values"/> for a tight inner loop. Frobenius inner product
    /// is invariant under any consistent linearization, so the flat array can be walked
    /// in storage order without re-permuting (i, j) indices. This skips MathNet's
    /// per-element method call + bounds check; at N=5 (1024×1024 matrices × 524k pairs)
    /// the difference is hours vs minutes single-threaded.</para>
    ///
    /// <para>Throws <see cref="InvalidCastException"/> if the inputs are not dense
    /// column-major matrices. Deliberately no fallback: a silent slow-path would
    /// regress the N=5 forecast.</para></summary>
    public static Complex FrobeniusInner(ComplexMatrix A, ComplexMatrix B)
    {
        if (A is null) throw new ArgumentNullException(nameof(A));
        if (B is null) throw new ArgumentNullException(nameof(B));
        if (A.RowCount != B.RowCount || A.ColumnCount != B.ColumnCount)
            throw new ArgumentException($"shape mismatch: A is {A.RowCount}x{A.ColumnCount}, B is {B.RowCount}x{B.ColumnCount}");

        var aData = ((DenseMatrix)A).Values;
        var bData = ((DenseMatrix)B).Values;

        Complex sum = Complex.Zero;
        for (int i = 0; i < aData.Length; i++)
            sum += Complex.Conjugate(aData[i]) * bData[i];
        return sum;
    }

    /// <summary>Result record from a single F112 non-Hermitian enumeration run.
    /// Mirrors the Python script's stdout summary.</summary>
    public sealed record EnumerationResult(
        int N,
        int TotalPairs,
        int NonzeroCount,
        double MaxImaginary,
        double MeanAbsImaginary,
        TimeSpan Elapsed,
        IReadOnlyList<(string Alpha, string Beta, double Imag)> NonzeroExamples);

    /// <summary>Enumerate all upper-triangular unordered pairs (σ_α, σ_β) of Pauli
    /// strings at length N, build L_α,-i and L_β,-i for each, compute the Frobenius
    /// inner product ⟨L_α,-i, L_β,-i⟩, and report the maximum |Im|, count of pairs
    /// above tolerance, and up to 10 non-zero examples.
    ///
    /// <para>If all pairs give Im &lt; tolerance, the F112 non-Hermitian extension is
    /// proven at this N via bilinearity + Pauli-basis spanning.</para>
    ///
    /// <para>Memory: pre-caches 4^N L_α,-i matrices, each 4^N × 4^N complex = 16·4^(2N)
    /// bytes. At N=5: 16 GB total cache. At N=6: 1 TB (infeasible without sparse rep).</para>
    /// </summary>
    public static EnumerationResult Enumerate(int N, double tolerance = 1e-10)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be ≥ 1");
        long count = 1L << (2 * N);  // 4^N
        if (count > int.MaxValue)
            throw new ArgumentOutOfRangeException(nameof(N), N, "N too large: 4^N overflows int32");
        int stringCount = (int)count;

        var stopwatch = System.Diagnostics.Stopwatch.StartNew();
        var pi = PiOperator.BuildFull(N, PauliLetter.Z);

        // Pre-compute L_α,-i for every Pauli string α
        var cache = new ComplexMatrix[stringCount];
        for (int k = 0; k < stringCount; k++)
        {
            var letters = PauliIndex.FromFlat(k, N);
            var sigma = PauliString.Build(letters);
            var L = BuildLHInPauliBasis(sigma, N);
            cache[k] = ProjectOntoPiEigenspace(L, pi, -Complex.ImaginaryOne);
        }

        // Iterate upper-triangular pairs, check |Im⟨L_α,-i, L_β,-i⟩|
        int totalPairs = 0;
        int nonzeroCount = 0;
        double maxIm = 0.0;
        double sumAbsIm = 0.0;
        var examples = new List<(string Alpha, string Beta, double Imag)>();
        for (int a = 0; a < stringCount; a++)
        {
            for (int b = a; b < stringCount; b++)
            {
                var inner = FrobeniusInner(cache[a], cache[b]);
                double absIm = Math.Abs(inner.Imaginary);
                if (absIm > maxIm) maxIm = absIm;
                sumAbsIm += absIm;
                if (absIm > tolerance)
                {
                    nonzeroCount++;
                    if (examples.Count < 10)
                    {
                        string alphaName = PauliLabel.Format(PauliIndex.FromFlat(a, N));
                        string betaName = PauliLabel.Format(PauliIndex.FromFlat(b, N));
                        examples.Add((alphaName, betaName, inner.Imaginary));
                    }
                }
                totalPairs++;
            }
        }

        stopwatch.Stop();
        double meanAbsIm = sumAbsIm / totalPairs;
        return new EnumerationResult(N, totalPairs, nonzeroCount, maxIm, meanAbsIm, stopwatch.Elapsed, examples);
    }
}
