using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Pauli;

/// <summary>4^N Pauli-basis decomposition of density matrices.
///
/// The Pauli strings σ_α (α ∈ [0, 4^N)) form an orthogonal basis of d×d matrices
/// (d = 2^N) under the Hilbert-Schmidt inner product: Tr(σ_α σ_β) = d · δ_αβ.
/// Any d×d operator A decomposes as
///   A = Σ_α vec[α] · σ_α,    vec[α] = Tr(σ_α · A) / d.
///
/// For a 1-qubit density matrix this reduces to the standard Bloch-vector form
/// ρ = (I + ⟨X⟩ X + ⟨Y⟩ Y + ⟨Z⟩ Z) / 2 with vec[I] = 1/2, vec[X] = ⟨X⟩/2, etc.
///
/// <para>Also exposes <see cref="VecToPauliBasisTransform"/>, the d² × 4^N matrix M with
/// columns vec_F(σ_α) (column-major flattening). Cached by N. M† M = 2^N · I. Used by
/// Symmetry.PalindromeResidual, Observables.PiProtectedObservables, and the F81 / DZero
/// Π-decomposition diagnostics — the canonical Pauli↔vec round-trip in matrix form.</para>
/// </summary>
public static class PauliBasis
{
    private static readonly System.Collections.Concurrent.ConcurrentDictionary<int, ComplexMatrix> _transformCache = new();

    /// <summary>The d² × 4^N matrix M with M[:, α] = vec_F(σ_α) column-stacked. Cached by N.
    /// Callers should treat the returned matrix as read-only.
    ///
    /// <para><b>Convention note (Welle 10d audit, 2026-05-27)</b>: T columns are vec_F
    /// (column-major) per the BuildVecToPauliBasisTransform loop. The standard codebase
    /// pattern combines T with a vec_R-style commutator superoperator L_vec = -i(H ⊗ I −
    /// I ⊗ H^T) (e.g. Liouvillian.cs, LindbladianBuilder.cs, BondPerturbation.cs,
    /// PauliDephasingDissipator.cs, PalindromeResidual.cs, PiDecomposition.cs,
    /// F112NonHermitianBasisEnumeration.cs). Reading T columns as vec_R implicitly maps
    /// σ_k ↦ σ_k^T, so the resulting Pauli-basis matrix is D · L_natural · D where
    /// D = diag((−1)^n_Y(k)) is a real diagonal unitary involution. A bonus structural
    /// identity surfaced: D · Π_Z · D = Π_Y, so the twist coincides with the Z↔Y
    /// dephasing-letter swap on the operator space.</para>
    ///
    /// <para><b>Refinement (2026-07-29)</b>: the 2026-05-27 note above explained the
    /// invisibility by D-conjugation invariance alone. That is the right account for a
    /// bare sandwich, but NOT for anything that then applies Π. D does not commute with
    /// Π (‖DΠ − ΠD‖ = 11.31 at N=3), so <see cref="Symmetry.PalindromeResidual"/> is not
    /// D-conjugation-invariant and its returned matrix is not D · M_true · D. What
    /// actually holds is stronger and carries a CONJUGATION: since D·Π_Z·D = Π_Y =
    /// conj(Π_Z) and L_natural is REAL in the Pauli basis (a Hermiticity-preserving map
    /// in a Hermitian basis), the returned residual is exactly
    /// <c>M_returned = D · conj(M_true) · D</c>. Frobenius norms, zero patterns and
    /// eigenvalue magnitudes survive; individual entry signs flip where
    /// n_Y(row) + n_Y(col) is odd, and every imaginary part is negated. Verified
    /// entry-wise to 1e-16 at N=2,3 on Heisenberg, transverse-field and T1 Liouvillians;
    /// the first two give a real residual where the conjugation is invisible, the T1 case
    /// gives a genuinely complex one where it is not.</para>
    ///
    /// <para>Consequences a new caller must respect. The Π-symmetric/antisymmetric split
    /// (M ± Π M Π†)/2 DOES commute with M ↦ D·conj(M)·D, so PiDecomposition's norms are
    /// safe. The finer ±i refinement (M_anti ∓ i·Π M_anti Π†)/2 does NOT: conj sends
    /// −i ↦ +i while D·conj(Π) = Π·D sends it back, which SWAPS the two projectors.
    /// So which component is called +1/2 is fixed by the stacking, and
    /// <c>Diagnostics.Polarity.PolarityCoordinates</c> reads its MPlusHalf / MMinusHalf
    /// (hence the sign of Asymmetry) in the column-stack one.</para>
    ///
    /// <para><b>That convention is load-bearing, so do NOT "correct" it.</b> F113
    /// (<see cref="Symmetry.LindbladBitBPiBreakMagnitude"/>) is a signed, nonzero,
    /// hardware-anchored result living on exactly this quantity: for Hermitian H with
    /// single-site Z-drives (ω_l/2)·Z_l plus σ⁻ amplitude damping, asymmetry =
    /// (4^N/2)·Σ_l ω_l·(γ_pump,l − γ_T1,l), so cooling-dominant systems give a NEGATIVE
    /// asymmetry. At ω = 0.13, γ_T1 = 0.001, N = 2 that is −2.08e-3, the Kingston fitted
    /// value. Recomputed from below on 2026-07-29: the column-stack reading reproduces
    /// −2.08e-3 (N=2) and −1.248e-2 (N=3) exactly, and a row-stack transform returns the
    /// same magnitudes with the OPPOSITE sign. Switching this module to row-stack would
    /// therefore flip a hardware-matched sign and silently break F113's agreement with
    /// its own closed form. Note the counterexample domain is Hermitian H; the asymmetry
    /// vanishes for bond Hamiltonians (Heisenberg, XXZ, with or without T1), which is why
    /// a survey over bond families alone reports it as identically zero and misses
    /// this.</para>
    ///
    /// <para>Two callers landed after the 2026-05-27 sweep and were re-checked on
    /// 2026-07-29: <c>Diagnostics.Foundation.SlowManifoldPauliContent</c> reads
    /// this note and deliberately builds its own row-major projection instead (correct by
    /// construction), and DirectSumDecompositionWitness reports only Frobenius and
    /// parity-block norms (invariant; D is diagonal and never moves mass between blocks).
    /// A caller that reads entry SIGNS or compares against a third-party natural L needs
    /// the (−1)^(n_Y(row) + n_Y(col)) correction AND, past a Π, the conjugation (see
    /// <c>F112NonHermitianBasisEnumeration.BuildSparseLSigma</c> for the sign half; it is
    /// a sparse-vs-dense parity fix within the twisted convention, not a correction back
    /// to L_natural).</para></summary>
    public static ComplexMatrix VecToPauliBasisTransform(int N) =>
        _transformCache.GetOrAdd(N, BuildVecToPauliBasisTransform);

    private static ComplexMatrix BuildVecToPauliBasisTransform(int N)
    {
        int d = 1 << N;
        long d2 = 1L << (2 * N);
        var M = Matrix<Complex>.Build.Dense(d * d, (int)d2);
        for (long k = 0; k < d2; k++)
        {
            var letters = PauliIndex.FromFlat(k, N);
            var sigma = PauliString.Build(letters);
            for (int j = 0; j < d; j++)
                for (int i = 0; i < d; i++)
                    M[j * d + i, (int)k] = sigma[i, j];
        }
        return M;
    }

    /// <summary>Decompose ρ into Pauli-basis coefficients: vec[α] = Tr(σ_α · ρ) / 2^N.</summary>
    public static ComplexVector ToPauliVector(ComplexMatrix rho, int N)
    {
        long d2 = 1L << (2 * N);
        double inv = 1.0 / (1 << N);
        var vec = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense((int)d2);
        for (long k = 0; k < d2; k++)
        {
            var letters = PauliIndex.FromFlat(k, N);
            var sigma = PauliString.Build(letters);
            vec[(int)k] = (sigma * rho).Trace() * inv;
        }
        return vec;
    }

    /// <summary>Reconstruct A = Σ_α vec[α] · σ_α from its Pauli-basis coefficients.</summary>
    public static ComplexMatrix FromPauliVector(ComplexVector pauliVec, int N)
    {
        long d2 = 1L << (2 * N);
        if (pauliVec.Count != d2)
            throw new ArgumentException($"expected length 4^N = {d2}; got {pauliVec.Count}");
        int d = 1 << N;
        var rho = Matrix<Complex>.Build.Dense(d, d);
        for (long k = 0; k < d2; k++)
        {
            if (pauliVec[(int)k] == Complex.Zero) continue;
            var letters = PauliIndex.FromFlat(k, N);
            rho = rho + pauliVec[(int)k] * PauliString.Build(letters);
        }
        return rho;
    }
}
