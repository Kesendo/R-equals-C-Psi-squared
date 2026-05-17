using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum.JordanWigner;

namespace RCPsiSquared.Core.BlockSpectrum;

/// <summary>Generic shift-invert Arnoldi on a CSR-stored complex matrix. Lifted from
/// <see cref="JordanWigner.JwSlaterPairShiftInvertArnoldi"/> to a source-agnostic algorithm:
/// inputs are the four raw CSR arrays + Krylov / inner-solve parameters, output is the
/// eigenvalues nearest σ. Reused by the JW Slater-pair path AND by the Klein-4-group
/// self-paired sub-block path (Phase 3c); a single point of algorithmic truth.
///
/// <para>Algorithm: shift-invert spectral transformation with Jacobi-preconditioned BiCGStab
/// inner solve and Modified Gram-Schmidt Arnoldi outer loop. See
/// <see cref="JordanWigner.JwSlaterPairShiftInvertArnoldi"/>'s class docstring for the
/// algorithmic details and the BiCGStab convergence story; this class is the pure-algorithm
/// extraction with no Claim-typing or source-specific metadata.</para></summary>
public static class SparseShiftInvertArnoldi
{
    public const double OuterBreakdownThreshold = 1e-14;
    public const double JacobiFloor = 1e-12;

    /// <summary>Run shift-invert Arnoldi on the CSR matrix (rowPtr, colIdx, values) at dim
    /// <paramref name="dim"/>. Returns the <paramref name="numEig"/> eigenvalues nearest
    /// <paramref name="sigma"/>, sorted by ascending |λ − σ|; plus diagnostic outputs
    /// (deflation flag, terminal Hessenberg subdiagonal magnitude, BiCGStab inner-iter
    /// counts per outer step).
    ///
    /// <para><paramref name="preconditioner"/> selects between Jacobi (diagonal) and Identity
    /// (no preconditioning). Jacobi is the default; Identity is sometimes required when the
    /// CSR matrix is not diagonally dominant (e.g. off-diagonal hopping dominates over
    /// dissipator-only diagonal entries, like in the Klein-projected computational basis at
    /// N=10 (5,5) where many rows have exact-zero diagonal). On such matrices, Jacobi inverse
    /// amplifies near-zero shifted diagonals to ~1/σ ≈ 1000, destabilising BiCGStab.</para></summary>
    public static ShiftInvertResult Run(int dim, int[] rowPtr, int[] colIdx, Complex[] values,
        Complex sigma, int numEig, int numIter, int randomSeed,
        double innerTolerance, int innerMaxIter,
        PreconditionerKind preconditioner = PreconditionerKind.Jacobi)
    {
        if (rowPtr is null) throw new ArgumentNullException(nameof(rowPtr));
        if (colIdx is null) throw new ArgumentNullException(nameof(colIdx));
        if (values is null) throw new ArgumentNullException(nameof(values));
        if (dim < 1) throw new ArgumentOutOfRangeException(nameof(dim), dim, "dim must be ≥ 1.");
        if (rowPtr.Length != dim + 1)
            throw new ArgumentException($"rowPtr length {rowPtr.Length} != dim+1 ({dim + 1})", nameof(rowPtr));
        if (numEig < 1) throw new ArgumentOutOfRangeException(nameof(numEig));
        if (numIter < 1) throw new ArgumentOutOfRangeException(nameof(numIter));
        if (numEig > numIter)
            throw new ArgumentException(
                $"numEig {numEig} must be ≤ numIter {numIter}; Arnoldi cannot return more Ritz values than Krylov dimension.",
                nameof(numEig));
        if (numIter >= dim)
            throw new ArgumentException($"numIter {numIter} must be < dim {dim}.", nameof(numIter));
        if (innerTolerance <= 0) throw new ArgumentOutOfRangeException(nameof(innerTolerance));
        if (innerMaxIter < 1) throw new ArgumentOutOfRangeException(nameof(innerMaxIter));

        var preconditionerDiag = preconditioner == PreconditionerKind.Jacobi
            ? BuildJacobiInverse(dim, rowPtr, colIdx, values, sigma)
            : BuildIdentityInverse(dim);

        var V = new Complex[numIter + 1][];
        V[0] = KrylovOps.RandomNormalized(dim, randomSeed);
        var H = new Complex[numIter + 1, numIter];
        var w = new Complex[dim];
        var innerIters = new List<int>(numIter);

        bool deflated = false;
        double terminalH = 0.0;
        int actualIter = numIter;

        for (int j = 0; j < numIter; j++)
        {
            int innerSteps = SolveShiftedSystem(dim, rowPtr, colIdx, values, sigma, preconditionerDiag,
                V[j], w, innerTolerance, innerMaxIter);
            innerIters.Add(innerSteps);

            for (int i = 0; i <= j; i++)
            {
                Complex hij = KrylovOps.ConjugateDot(V[i], w);
                H[i, j] = hij;
                KrylovOps.AxpyInPlace(w, V[i], -hij);
            }

            double wNorm = Math.Sqrt(KrylovOps.NormSquared(w));
            H[j + 1, j] = new Complex(wNorm, 0.0);

            if (wNorm < OuterBreakdownThreshold)
            {
                deflated = true;
                terminalH = wNorm;
                actualIter = j + 1;
                break;
            }

            if (j + 1 < V.Length)
            {
                V[j + 1] = new Complex[dim];
                double invW = 1.0 / wNorm;
                for (int i = 0; i < dim; i++) V[j + 1][i] = w[i] * invW;
            }
            terminalH = wNorm;
        }

        var Hm = Matrix<Complex>.Build.Dense(actualIter, actualIter, (i, k) => H[i, k]);
        var ritzMu = Hm.Evd().EigenValues.ToArray();

        var recovered = new List<Complex>(ritzMu.Length);
        foreach (var mu in ritzMu)
        {
            if (mu.Magnitude < OuterBreakdownThreshold) continue;
            recovered.Add(sigma + Complex.Reciprocal(mu));
        }
        var sorted = recovered.OrderBy(l => (l - sigma).Magnitude).Take(numEig).ToArray();

        return new ShiftInvertResult(sorted, deflated, terminalH, innerIters.ToArray());
    }

    /// <summary>Apply y = (L − σI)·x using the CSR matvec plus a diagonal shift.</summary>
    private static void ApplyShiftedMatvec(int dim, int[] rowPtr, int[] colIdx, Complex[] values,
        Complex sigma, Complex[] x, Complex[] y)
    {
        Parallel.For(0, dim, alpha =>
        {
            Complex sum = -sigma * x[alpha];
            int start = rowPtr[alpha];
            int end = rowPtr[alpha + 1];
            for (int e = start; e < end; e++)
                sum += values[e] * x[colIdx[e]];
            y[alpha] = sum;
        });
    }

    /// <summary>Right-preconditioned BiCGStab inner solve for (L − σI)·x = b with Jacobi
    /// (diagonal) preconditioning M = diag(L − σI). Initial guess x ≡ 0.</summary>
    private static int SolveShiftedSystem(int n, int[] rowPtr, int[] colIdx, Complex[] values,
        Complex sigma, Complex[] jacobiInv, Complex[] b, Complex[] x, double tol, int maxIter)
    {
        var r = new Complex[n];
        var rt = new Complex[n];
        var p = new Complex[n];
        var v = new Complex[n];
        var s = new Complex[n];
        var t = new Complex[n];
        var z = new Complex[n];
        var y = new Complex[n];

        Array.Clear(x, 0, n);
        Array.Copy(b, r, n);
        Array.Copy(b, rt, n);

        double bNorm = Math.Sqrt(KrylovOps.NormSquared(b));
        if (bNorm < 1e-300) return 0;

        Complex rho = Complex.One, alpha = Complex.One, omega = Complex.One;
        Complex rhoPrev;

        for (int k = 1; k <= maxIter; k++)
        {
            rhoPrev = rho;
            rho = KrylovOps.ConjugateDot(rt, r);
            if (rho.Magnitude < 1e-300) return k;

            if (k == 1)
            {
                Array.Copy(r, p, n);
            }
            else
            {
                if (omega.Magnitude < 1e-300) return k;
                Complex beta = (rho / rhoPrev) * (alpha / omega);
                for (int i = 0; i < n; i++) p[i] = r[i] + beta * (p[i] - omega * v[i]);
            }

            for (int i = 0; i < n; i++) z[i] = jacobiInv[i] * p[i];
            ApplyShiftedMatvec(n, rowPtr, colIdx, values, sigma, z, v);
            Complex rtv = KrylovOps.ConjugateDot(rt, v);
            if (rtv.Magnitude < 1e-300) return k;
            alpha = rho / rtv;

            for (int i = 0; i < n; i++) s[i] = r[i] - alpha * v[i];
            double sNorm = Math.Sqrt(KrylovOps.NormSquared(s));
            if (sNorm / bNorm < tol)
            {
                for (int i = 0; i < n; i++) x[i] += alpha * z[i];
                return k;
            }

            for (int i = 0; i < n; i++) y[i] = jacobiInv[i] * s[i];
            ApplyShiftedMatvec(n, rowPtr, colIdx, values, sigma, y, t);
            Complex tt = KrylovOps.ConjugateDot(t, t);
            if (tt.Magnitude < 1e-300) return k;
            omega = KrylovOps.ConjugateDot(t, s) / tt;

            for (int i = 0; i < n; i++)
            {
                x[i] += alpha * z[i] + omega * y[i];
                r[i] = s[i] - omega * t[i];
            }

            double rNorm = Math.Sqrt(KrylovOps.NormSquared(r));
            if (rNorm / bNorm < tol) return k;
            if (omega.Magnitude < 1e-300) return k;
        }
        return maxIter;
    }

    /// <summary>Identity preconditioner: M = I, M^{−1} = I. Returns all-ones for use in
    /// the same per-element preconditioning step. For matrices that are NOT diagonally
    /// dominant (e.g. Klein-projected computational-basis sub-blocks where the diagonal is
    /// only the dissipator −2γ·hamming, much smaller than off-diagonal hopping ±iJ),
    /// Jacobi preconditioning amplifies near-zero shifted diagonals to ~1/σ and destabilises
    /// BiCGStab; using Identity instead preserves the natural conditioning of (L − σI).</summary>
    private static Complex[] BuildIdentityInverse(int dim)
    {
        var inv = new Complex[dim];
        for (int i = 0; i < dim; i++) inv[i] = Complex.One;
        return inv;
    }

    /// <summary>Pre-compute jacobiInv[α] = 1 / (L[α, α] − σ) with a floor fallback to
    /// avoid division by near-zero on near-singular shifts.</summary>
    private static Complex[] BuildJacobiInverse(int dim, int[] rowPtr, int[] colIdx,
        Complex[] values, Complex sigma)
    {
        var jacobiInv = new Complex[dim];
        for (int alpha = 0; alpha < dim; alpha++)
        {
            Complex d = Complex.Zero;
            int start = rowPtr[alpha];
            int end = rowPtr[alpha + 1];
            for (int e = start; e < end; e++)
            {
                if (colIdx[e] == alpha) { d = values[e]; break; }
            }
            Complex shifted = d - sigma;
            jacobiInv[alpha] = shifted.Magnitude < JacobiFloor ? Complex.One : Complex.Reciprocal(shifted);
        }
        return jacobiInv;
    }
}

/// <summary>Preconditioner selector for <see cref="SparseShiftInvertArnoldi.Run"/>.
/// Jacobi (diagonal) is the default and is appropriate when the CSR matrix has dominant
/// diagonal entries (e.g. JW Slater-pair Liouvillian where −i(ε_L − ε_K) lives on the
/// diagonal). Identity (no preconditioning) is required when the CSR matrix is structurally
/// not diagonally dominant — notably the Klein-projected computational-basis sub-block,
/// where the diagonal is only the dissipator term and many rows have exact-zero diagonal,
/// so Jacobi inverse amplifies BiCGStab residuals catastrophically.</summary>
public enum PreconditionerKind
{
    Jacobi,
    Identity,
}

/// <summary>Output of <see cref="SparseShiftInvertArnoldi.Run"/>: the recovered eigenvalues
/// nearest σ (sorted by ascending |λ − σ|), the deflation flag, the terminal Hessenberg
/// subdiagonal magnitude, and the per-outer-iteration BiCGStab iteration counts.</summary>
public sealed record ShiftInvertResult(
    Complex[] Eigenvalues,
    bool DeflatedEarly,
    double TerminalSubdiagonalMagnitude,
    int[] InnerIterationsPerOuter)
{
    public int MaxInnerIterations => InnerIterationsPerOuter.Length == 0 ? 0 : InnerIterationsPerOuter.Max();
    public double MeanInnerIterations => InnerIterationsPerOuter.Length == 0 ? 0.0 :
        InnerIterationsPerOuter.Average();
}
