using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.BlockSpectrum.JordanWigner;

/// <summary>Shift-invert Arnoldi on the CSR sparse L_JW from
/// <see cref="JwSlaterPairSparseLBuilder"/>. Extracts the eigenvalues nearest a user-
/// supplied complex shift σ by running plain Arnoldi on the operator
/// <c>(L − σI)^(−1)</c>, where each shift-invert step is solved by an in-house BiCGStab
/// on the same CSR matvec — no MathNet sparse conversion needed.
///
/// <para>For σ ≈ 0 this delivers the physically interesting slow modes (Re λ ≈ 0) —
/// the decoherence-time-scale eigenstates that bridge from the steady state into the
/// decaying spectrum. Plain Arnoldi (<see cref="JwSlaterPairArnoldiEig"/>) gives only the
/// largest-magnitude end, which for a Liouvillian is dominated by oscillation rather than
/// slow decay.</para>
///
/// <para>Algorithm: standard shift-invert spectral transformation.
/// <list type="bullet">
///   <item><b>Inner</b> (per Arnoldi step): solve <c>(L − σI)·w = v_j</c> for w using
///         BiCGStab — Bi-Conjugate Gradient Stabilised, robust for moderately non-normal
///         complex systems; convergence to <see cref="InnerTolerance"/> relative residual
///         within <see cref="InnerMaxIterations"/> matvecs.</item>
///   <item><b>Outer</b>: Modified Gram-Schmidt Arnoldi on the shift-invert operator, same
///         pattern as <see cref="JwSlaterPairArnoldiEig"/>. Build Hessenberg H, extract
///         Ritz μ values, recover λ via <c>λ = σ + 1/μ</c>. Largest |μ| corresponds to λ
///         closest to σ in the Euclidean sense.</item>
/// </list></para>
///
/// <para><b>Tier outcome: <see cref="Tier.Tier1Derived"/>.</b> Textbook shift-invert
/// transformation; both BiCGStab and Arnoldi are standard algorithms. Validated against
/// MathNet dense Evd at N = 4 / 5 — recovered eigenvalues match the dense-Evd entries
/// nearest σ within the test-stated tolerance.</para>
///
/// <para>Anchor: <see cref="JwSlaterPairSparseLBuilder"/> (CSR source) +
/// <see cref="JwSlaterPairArnoldiEig"/> (sister plain-Arnoldi primitive); textbook
/// shift-invert Arnoldi (Saad Ch. 6) + BiCGStab (van der Vorst 1992 / Saad Ch. 7).</para>
/// </summary>
public sealed class JwSlaterPairShiftInvertArnoldi : Claim
{
    public const double OuterBreakdownThreshold = 1e-14;

    public JwSlaterPairSparseLBuilder Source { get; }
    public Complex Sigma { get; }
    public int NumEigenvaluesRequested { get; }
    public int NumIterations { get; }
    public double InnerTolerance { get; }
    public int InnerMaxIterations { get; }
    public int RandomSeed { get; }

    /// <summary>Eigenvalues nearest <see cref="Sigma"/>, sorted by ascending |λ − σ|, trimmed
    /// to <see cref="NumEigenvaluesRequested"/> entries. Recovered via <c>λ = σ + 1/μ</c>
    /// from the largest-magnitude Ritz values μ of the shift-invert operator.</summary>
    public Complex[] Eigenvalues { get; }

    public bool DeflatedEarly { get; }
    public double TerminalSubdiagonalMagnitude { get; }

    /// <summary>BiCGStab iteration counts per outer Arnoldi step, length
    /// <see cref="NumIterations"/> (or shorter if deflated). Useful to confirm the inner
    /// solver converges within budget; if any entry == <see cref="InnerMaxIterations"/>,
    /// the corresponding step exhausted its budget without reaching tolerance.</summary>
    public int[] InnerIterationsPerOuter { get; }
    public int MaxInnerIterations => InnerIterationsPerOuter.Length == 0 ? 0 : InnerIterationsPerOuter.Max();
    public double MeanInnerIterations => InnerIterationsPerOuter.Length == 0 ? 0.0 :
        InnerIterationsPerOuter.Average();

    public static JwSlaterPairShiftInvertArnoldi Build(JwSlaterPairSparseLBuilder source,
        Complex sigma, int numEig, int numIter, int randomSeed,
        double innerTolerance, int innerMaxIter)
    {
        if (source is null) throw new ArgumentNullException(nameof(source));
        if (numEig < 1) throw new ArgumentOutOfRangeException(nameof(numEig), numEig, "numEig must be ≥ 1.");
        if (numIter < 1) throw new ArgumentOutOfRangeException(nameof(numIter), numIter, "numIter must be ≥ 1.");
        if (numEig > numIter)
            throw new ArgumentException(
                $"numEig {numEig} must be ≤ numIter {numIter}; Arnoldi cannot return more Ritz values than Krylov dimension.",
                nameof(numEig));

        int dim = source.SectorDim;
        if (numIter >= dim)
            throw new ArgumentException(
                $"numIter {numIter} must be < sectorDim {dim}.", nameof(numIter));
        if (innerTolerance <= 0) throw new ArgumentOutOfRangeException(nameof(innerTolerance));
        if (innerMaxIter < 1) throw new ArgumentOutOfRangeException(nameof(innerMaxIter));

        var rng = new Random(randomSeed);
        var v0 = new Complex[dim];
        double v0NormSq = 0.0;
        for (int i = 0; i < dim; i++)
        {
            v0[i] = new Complex(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5);
            v0NormSq += v0[i].Real * v0[i].Real + v0[i].Imaginary * v0[i].Imaginary;
        }
        double inv0 = 1.0 / Math.Sqrt(v0NormSq);
        for (int i = 0; i < dim; i++) v0[i] *= inv0;

        var V = new Complex[numIter + 1][];
        V[0] = v0;
        var H = new Complex[numIter + 1, numIter];
        var w = new Complex[dim];
        var innerIters = new List<int>(numIter);

        bool deflated = false;
        double terminalH = 0.0;
        int actualIter = numIter;

        for (int j = 0; j < numIter; j++)
        {
            // w = (L − σI)^(−1) V[j], computed by BiCGStab.
            int innerSteps = SolveShiftedSystem(source, sigma, V[j], w, innerTolerance, innerMaxIter);
            innerIters.Add(innerSteps);

            // Modified Gram-Schmidt against V[0..j].
            for (int i = 0; i <= j; i++)
            {
                Complex hij = ConjugateDot(V[i], w);
                H[i, j] = hij;
                AxpyInPlace(w, V[i], -hij);
            }

            double wNormSq = 0.0;
            for (int i = 0; i < dim; i++) wNormSq += w[i].Real * w[i].Real + w[i].Imaginary * w[i].Imaginary;
            double wNorm = Math.Sqrt(wNormSq);
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

        // Hessenberg eigenvalues = Ritz values μ of (L − σI)^(−1).
        var Hm = Matrix<Complex>.Build.Dense(actualIter, actualIter, (i, k) => H[i, k]);
        var ritzMu = Hm.Evd().EigenValues.ToArray();

        // Recover λ = σ + 1/μ; sort by ascending distance from σ.
        var recovered = new List<Complex>(ritzMu.Length);
        foreach (var mu in ritzMu)
        {
            if (mu.Magnitude < OuterBreakdownThreshold) continue;  // skip near-infinite λ
            recovered.Add(sigma + Complex.Reciprocal(mu));
        }
        var sorted = recovered.OrderBy(l => (l - sigma).Magnitude).Take(numEig).ToArray();

        return new JwSlaterPairShiftInvertArnoldi(source, sigma, numEig, numIter, randomSeed,
            innerTolerance, innerMaxIter, sorted, deflated, terminalH, innerIters.ToArray());
    }

    /// <summary>Apply <c>y = (L − σI) x</c> using the CSR matvec plus a diagonal shift.
    /// Parallelised across rows; each thread writes only its own <c>y[alpha]</c>.</summary>
    private static void ApplyShiftedMatvec(JwSlaterPairSparseLBuilder src, Complex sigma, Complex[] x, Complex[] y)
    {
        int dim = src.SectorDim;
        var rowPtr = src.RowPtr;
        var colIdx = src.ColIdx;
        var values = src.Values;
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

    /// <summary>BiCGStab inner solve for <c>(L − σI)·x = b</c>. Initial guess x ≡ 0.
    /// Returns the iteration count at which <c>‖r‖ / ‖b‖ &lt; tol</c> was reached, or
    /// <paramref name="maxIter"/> if convergence failed within budget.</summary>
    private static int SolveShiftedSystem(JwSlaterPairSparseLBuilder src, Complex sigma,
        Complex[] b, Complex[] x, double tol, int maxIter)
    {
        int n = src.SectorDim;
        var r = new Complex[n];
        var rt = new Complex[n];
        var p = new Complex[n];
        var v = new Complex[n];
        var s = new Complex[n];
        var t = new Complex[n];

        // x_0 = 0  →  r_0 = b − A·x_0 = b
        Array.Clear(x, 0, n);
        Array.Copy(b, r, n);
        Array.Copy(r, rt, n);  // r̃_0 = r_0 (van der Vorst's default choice).

        double bNorm = Math.Sqrt(b.Sum(z => z.Real * z.Real + z.Imaginary * z.Imaginary));
        if (bNorm < 1e-300) return 0;

        Complex rho = Complex.One, alpha = Complex.One, omega = Complex.One;
        Complex rhoPrev;

        for (int k = 1; k <= maxIter; k++)
        {
            rhoPrev = rho;
            rho = ConjugateDot(rt, r);
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

            ApplyShiftedMatvec(src, sigma, p, v);
            Complex rtv = ConjugateDot(rt, v);
            if (rtv.Magnitude < 1e-300) return k;
            alpha = rho / rtv;

            for (int i = 0; i < n; i++) s[i] = r[i] - alpha * v[i];
            double sNorm = Math.Sqrt(s.Sum(z => z.Real * z.Real + z.Imaginary * z.Imaginary));
            if (sNorm / bNorm < tol)
            {
                for (int i = 0; i < n; i++) x[i] += alpha * p[i];
                return k;
            }

            ApplyShiftedMatvec(src, sigma, s, t);
            Complex tt = ConjugateDot(t, t);
            if (tt.Magnitude < 1e-300) return k;
            omega = ConjugateDot(t, s) / tt;

            for (int i = 0; i < n; i++)
            {
                x[i] += alpha * p[i] + omega * s[i];
                r[i] = s[i] - omega * t[i];
            }

            double rNorm = Math.Sqrt(r.Sum(z => z.Real * z.Real + z.Imaginary * z.Imaginary));
            if (rNorm / bNorm < tol) return k;
            if (omega.Magnitude < 1e-300) return k;
        }
        return maxIter;
    }

    private static Complex ConjugateDot(Complex[] a, Complex[] b)
    {
        Complex sum = Complex.Zero;
        int n = a.Length;
        for (int i = 0; i < n; i++) sum += Complex.Conjugate(a[i]) * b[i];
        return sum;
    }

    private static void AxpyInPlace(Complex[] y, Complex[] x, Complex alpha)
    {
        int n = y.Length;
        for (int i = 0; i < n; i++) y[i] += alpha * x[i];
    }

    private JwSlaterPairShiftInvertArnoldi(JwSlaterPairSparseLBuilder source,
        Complex sigma, int numEig, int numIter, int randomSeed,
        double innerTol, int innerMaxIter,
        Complex[] eigenvalues, bool deflated, double terminalH, int[] innerIters)
        : base($"Shift-invert Arnoldi at σ=({sigma.Real:G3}, {sigma.Imaginary:G3}) — top-{numEig} eigenvalues " +
               $"nearest σ for sparse L_JW (p_c={source.PCol}, p_r={source.PRow}, N={source.N}, " +
               $"dim={source.SectorDim}); numIter={numIter}, terminal-|H[j+1,j]|={terminalH:G3}" +
               (deflated ? " (deflated early)" : "") + $", inner BiCGStab mean iter={(innerIters.Length == 0 ? 0.0 : innerIters.Average()):F1}.",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/BlockSpectrum/JordanWigner/JwSlaterPairSparseLBuilder.cs (CSR L_JW source) + " +
               "compute/RCPsiSquared.Core/BlockSpectrum/JordanWigner/JwSlaterPairArnoldiEig.cs (sister plain-Arnoldi primitive); " +
               "textbook shift-invert Arnoldi (Saad Ch. 6) and BiCGStab (van der Vorst 1992 / Saad Ch. 7).")
    {
        Source = source;
        Sigma = sigma;
        NumEigenvaluesRequested = numEig;
        NumIterations = numIter;
        RandomSeed = randomSeed;
        InnerTolerance = innerTol;
        InnerMaxIterations = innerMaxIter;
        Eigenvalues = eigenvalues;
        DeflatedEarly = deflated;
        TerminalSubdiagonalMagnitude = terminalH;
        InnerIterationsPerOuter = innerIters;
    }

    public override string DisplayName =>
        $"Shift-invert Arnoldi @ σ=({Sigma.Real:G3}, {Sigma.Imaginary:G3}), top-{NumEigenvaluesRequested} " +
        $"of L_JW (p_c={Source.PCol}, p_r={Source.PRow}, N={Source.N})";

    public override string Summary =>
        $"|λ−σ|_min={(Eigenvalues.Length == 0 ? 0.0 : (Eigenvalues[0] - Sigma).Magnitude):G3}, " +
        $"mean inner BiCGStab iter={MeanInnerIterations:F1}, max={MaxInnerIterations} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return Source;
            yield return new InspectableNode("σ", summary: $"({Sigma.Real:G4}, {Sigma.Imaginary:G4})");
            yield return InspectableNode.RealScalar("numEig requested", NumEigenvaluesRequested);
            yield return InspectableNode.RealScalar("numIter", NumIterations);
            yield return InspectableNode.RealScalar("inner tolerance", InnerTolerance, "G3");
            yield return InspectableNode.RealScalar("inner max iter", InnerMaxIterations);
            yield return InspectableNode.RealScalar("deflated early", DeflatedEarly ? 1 : 0);
            yield return InspectableNode.RealScalar("terminal |H[j+1,j]|", TerminalSubdiagonalMagnitude, "G3");
            yield return InspectableNode.RealScalar("inner-iter mean", MeanInnerIterations, "F1");
            yield return InspectableNode.RealScalar("inner-iter max", MaxInnerIterations);
            for (int i = 0; i < Eigenvalues.Length; i++)
            {
                var e = Eigenvalues[i];
                yield return new InspectableNode($"λ_{i}", summary: $"({e.Real:F6}, {e.Imaginary:F6}), |λ−σ|={(e - Sigma).Magnitude:F6}");
            }
        }
    }
}
