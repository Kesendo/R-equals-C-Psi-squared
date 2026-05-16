using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.BlockSpectrum.JordanWigner;

/// <summary>Top-k largest-magnitude eigenvalues of the sparse L_JW from
/// <see cref="JwSlaterPairSparseLBuilder"/>, computed by managed Arnoldi iteration with
/// Modified Gram-Schmidt orthogonalization. No native dependencies; runs at the N=10
/// max-block where dense LAPACK <c>zgeev</c> is infeasible (~64 GB matrix storage).
///
/// <para>The Arnoldi method builds an orthonormal Krylov basis
/// <c>V = [v_0, A·v_0, A²·v_0, …, A^(m−1)·v_0]</c> (in factored form via Hessenberg
/// reduction) and extracts Ritz values from the small <c>m×m</c> Hessenberg matrix
/// <c>H = V^†·A·V</c>. Largest-magnitude Ritz values converge to the corresponding
/// eigenvalues of <c>A</c> as <c>m</c> grows; tolerance is governed by the residual
/// <c>‖A·v_j − H[:j+1, j]·V[:, :j+1]‖</c>, which is reflected in the subdiagonal entry
/// <c>H[j+1, j]</c>.</para>
///
/// <para><b>Slow-mode trick under uniform γ:</b> the chain XY + Z-dephasing Liouvillian has
/// F1 palindromic spectrum around <c>−Σγ</c>: every eigenvalue λ_fast has a mirror
/// <c>λ_slow = −2·Σγ − λ_fast</c>. Largest-magnitude Arnoldi extracts the fast-decay
/// (most-negative Re λ) end of the spectrum; the slow modes follow by reflection, giving
/// the physically interesting end "for free" without shift-invert.</para>
///
/// <para><b>Tier outcome: <see cref="Tier.Tier1Derived"/>.</b> Textbook non-symmetric
/// Arnoldi; iterations bounded by user-supplied <see cref="NumIterations"/>; deflation on
/// near-zero Hessenberg subdiagonal (lucky breakdown). Validated against MathNet dense Evd
/// at N=4–6 to relative tolerance <c>1e-6</c>; F1 palindrome witness at uniform γ confirms
/// each Arnoldi eigenvalue's mirror is present in the full spectrum to <c>1e-8</c>.</para>
///
/// <para>Anchor: <see cref="JwSlaterPairSparseLBuilder"/> (sparse L_JW source) + F1
/// palindrome (uniform γ slow-mode reflection); textbook Saad "Numerical Methods for Large
/// Eigenvalue Problems" Ch. 6 (non-symmetric Arnoldi with Hessenberg projection).</para>
/// </summary>
public sealed class JwSlaterPairArnoldiEig : Claim
{
    /// <summary>Subdiagonal threshold for lucky-breakdown deflation: when
    /// <c>|H[j+1, j]| &lt; BreakdownThreshold</c>, the Krylov subspace is invariant under A
    /// at depth <c>j+1</c> and the Hessenberg eigenvalues are exact.</summary>
    public const double BreakdownThreshold = 1e-14;

    public JwSlaterPairSparseLBuilder Source { get; }
    public int NumEigenvaluesRequested { get; }
    public int NumIterations { get; }
    public int RandomSeed { get; }

    /// <summary>Ritz values from the Hessenberg projection, sorted by descending magnitude
    /// and trimmed to the first <see cref="NumEigenvaluesRequested"/> entries.</summary>
    public Complex[] Eigenvalues { get; }

    /// <summary>True if the Krylov iteration deflated on a near-zero subdiagonal
    /// (<see cref="BreakdownThreshold"/>) — indicates an invariant subspace was found
    /// exactly; Ritz values from a deflated run are exact for the subspace.</summary>
    public bool DeflatedEarly { get; }

    /// <summary>Subdiagonal entry magnitude at termination: <c>|H[end+1, end]|</c>. Below
    /// <see cref="BreakdownThreshold"/> means deflation; above means the iteration was
    /// truncated at <see cref="NumIterations"/> and convergence is governed by this value
    /// times the Ritz-vector coupling.</summary>
    public double TerminalSubdiagonalMagnitude { get; }

    public static JwSlaterPairArnoldiEig Build(JwSlaterPairSparseLBuilder source,
        int numEig, int numIter, int randomSeed = 0)
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

        // Krylov basis storage: V[k] is the k-th basis vector (Complex[dim]).
        var V = new Complex[numIter + 1][];
        V[0] = KrylovOps.RandomNormalized(dim, randomSeed);
        var H = new Complex[numIter + 1, numIter];
        var w = new Complex[dim];

        bool deflated = false;
        double terminalH = 0.0;
        int actualIter = numIter;

        for (int j = 0; j < numIter; j++)
        {
            SparseMatVec(source, V[j], w);

            // Modified Gram-Schmidt against V[0..j].
            for (int i = 0; i <= j; i++)
            {
                Complex hij = KrylovOps.ConjugateDot(V[i], w);
                H[i, j] = hij;
                KrylovOps.AxpyInPlace(w, V[i], -hij);
            }

            double wNorm = Math.Sqrt(KrylovOps.NormSquared(w));
            H[j + 1, j] = new Complex(wNorm, 0.0);

            if (wNorm < BreakdownThreshold)
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

        // Extract eigenvalues from the upper Hessenberg H[:actualIter, :actualIter] (dense).
        var Hm = Matrix<Complex>.Build.Dense(actualIter, actualIter, (i, k) => H[i, k]);
        var evd = Hm.Evd();
        var ritzValues = evd.EigenValues.ToArray();

        var sorted = ritzValues.OrderByDescending(z => z.Magnitude).Take(numEig).ToArray();

        return new JwSlaterPairArnoldiEig(source, numEig, numIter, randomSeed, sorted, deflated, terminalH);
    }

    private static void SparseMatVec(JwSlaterPairSparseLBuilder src, Complex[] x, Complex[] y)
    {
        int dim = src.SectorDim;
        var rowPtr = src.RowPtr;
        var colIdx = src.ColIdx;
        var values = src.Values;
        // Rows are independent; safe to parallelise (each thread writes only its own y[alpha]).
        Parallel.For(0, dim, alpha =>
        {
            Complex sum = Complex.Zero;
            int start = rowPtr[alpha];
            int end = rowPtr[alpha + 1];
            for (int e = start; e < end; e++)
                sum += values[e] * x[colIdx[e]];
            y[alpha] = sum;
        });
    }

    private JwSlaterPairArnoldiEig(JwSlaterPairSparseLBuilder source,
        int numEig, int numIter, int randomSeed,
        Complex[] eigenvalues, bool deflated, double terminalH)
        : base($"Arnoldi top-{numEig} eigenvalues of sparse L_JW (p_c={source.PCol}, p_r={source.PRow}, " +
               $"N={source.N}, dim={source.SectorDim}); numIter={numIter}, terminal-|H[j+1,j]|={terminalH:G3}" +
               (deflated ? " (deflated early — exact)" : "") + ".",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/BlockSpectrum/JordanWigner/JwSlaterPairSparseLBuilder.cs (sparse L_JW source) + " +
               "F1 palindromic-spectrum theorem (slow modes via mirror reflection under uniform γ); " +
               "textbook non-symmetric Arnoldi with Modified Gram-Schmidt — Saad, 'Numerical Methods for Large Eigenvalue Problems', Ch. 6.")
    {
        Source = source;
        NumEigenvaluesRequested = numEig;
        NumIterations = numIter;
        RandomSeed = randomSeed;
        Eigenvalues = eigenvalues;
        DeflatedEarly = deflated;
        TerminalSubdiagonalMagnitude = terminalH;
    }

    public override string DisplayName =>
        $"Arnoldi top-{NumEigenvaluesRequested} of L_JW (p_c={Source.PCol}, p_r={Source.PRow}, N={Source.N})";

    public override string Summary =>
        $"|λ|_max={Eigenvalues[0].Magnitude:F3}, terminal-|H[j+1,j]|={TerminalSubdiagonalMagnitude:G3} " +
        (DeflatedEarly ? "(deflated) " : "") + $"({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return Source;
            yield return InspectableNode.RealScalar("numEig requested", NumEigenvaluesRequested);
            yield return InspectableNode.RealScalar("numIter", NumIterations);
            yield return InspectableNode.RealScalar("randomSeed", RandomSeed);
            yield return InspectableNode.RealScalar("deflated early", DeflatedEarly ? 1 : 0);
            yield return InspectableNode.RealScalar("terminal |H[j+1,j]|", TerminalSubdiagonalMagnitude, "G3");
            for (int i = 0; i < Eigenvalues.Length; i++)
            {
                var e = Eigenvalues[i];
                yield return new InspectableNode($"λ_{i}", summary: $"({e.Real:F4}, {e.Imaginary:F4}) |λ|={e.Magnitude:F4}");
            }
        }
    }
}
