using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Lindblad;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The Object Manager's read at every θ of a <see cref="DimensionAxis"/>: the marks
/// (the Liouvillian eigenvalue multiset, the fixed contract) against the in-between (how the
/// slow eigen-subspace rotates, the content).
///
/// <para>For the <see cref="DimensionAxis.Crossover"/> axis the two readings split cleanly:
/// <list type="bullet">
///   <item><b>Marks.</b> L(θ) is a similarity transform of L(0), so the sorted eigenvalues are
///         θ-invariant up to the dense-Evd numerical floor.
///         <see cref="DimensionSweepResult.MaxEigenvalueDriftAcrossTheta"/> measures exactly
///         that residual; "small" is the claim, not "exactly zero".</item>
///   <item><b>In-between.</b> The orthonormal projector onto the slow subspace rotates with
///         R_z(θ). <see cref="DimensionSweepResult.SubspaceRotation"/> is the Frobenius distance
///         between consecutive projectors, which is &gt; 0 because R_z(θ) does not commute with
///         the X-containing bond.</item>
/// </list></para>
///
/// <para>The projector is built from the orthonormalized slow eigenvectors precisely so that
/// the metric is gauge-robust: it does not depend on per-eigenvector phase, scale, or the
/// arbitrary basis chosen inside a degenerate eigenspace. Raw eigenvectors are never compared
/// directly.</para>
///
/// <para>The per-step <see cref="DimensionSweepResult.SubspaceRotation"/> is a coarse eyepiece:
/// the chordal (Frobenius) projector distance saturates near its ceiling √(2k) once the subspace
/// has turned a finite amount, so it cannot resolve how far the in-between has rotated.
/// <see cref="DimensionSweepResult.CumulativeRotation"/> is the resolving eyepiece: the largest
/// principal angle between the slow subspace at θ₀ and at θ[p]. Principal angles grow ~linearly
/// for a small rotation, so they keep resolving where the chordal distance has already flattened.
/// They are read from the singular values of Q(θ₀)ᴴ·Q(θ[p]) via
/// <see cref="DimensionSweepResult.SlowBasis"/>, the per-θ orthonormal Q.</para>
/// </summary>
public sealed record DimensionSweepResult(
    IReadOnlyList<double> Theta,
    IReadOnlyList<Complex[]> Eigenvalues,
    IReadOnlyList<double> Polarity,
    IReadOnlyList<double> SubspaceRotation,
    double MaxEigenvalueDriftAcrossTheta,
    IReadOnlyList<ComplexMatrix> SlowBasis,
    IReadOnlyList<double> CumulativeRotation);

public static class DimensionSweep
{
    /// <summary>Sweep <paramref name="axis"/>, reading the marks and the in-between at every θ.
    /// The slow subspace is spanned by the <paramref name="slowCount"/> eigenvalues with the
    /// smallest |Re λ| (the longest-lived Liouvillian modes).</summary>
    public static DimensionSweepResult Compute(DimensionAxis axis, int slowCount)
    {
        if (axis is null) throw new ArgumentNullException(nameof(axis));
        if (slowCount < 1) throw new ArgumentOutOfRangeException(nameof(slowCount), $"need at least one slow mode; got {slowCount}");

        int points = axis.Theta.Count;
        var sortedEigenvalues = new Complex[points][];   // by (Re, then Im), for the marks-drift read
        var projectors = new ComplexMatrix[points];      // orthonormal slow-subspace projectors, for the in-between
        var slowBasis = new ComplexMatrix[points];        // per-θ orthonormal Q, for the cumulative principal-angle read
        var polarity = new double[points];

        for (int p = 0; p < points; p++)
        {
            double theta = axis.Theta[p];
            var L = PauliDephasingDissipator.BuildZ(axis.Hamiltonian(theta), axis.GammaPerSite);
            var evd = L.Evd();
            var evals = evd.EigenValues;
            int dim = evals.Count;

            // Marks: the eigenvalue multiset, sorted by (Re, then Im) so θ-to-θ comparison is
            // index-aligned (the multiset is similarity-invariant; sorting fixes the labelling).
            // The keys are rounded before comparison: the slow modes come in degenerate clusters
            // with equal real parts up to the Evd floor (~1e-15), and a raw (Re, Im) sort would
            // let that noise flip the primary order inside a cluster, mis-pairing a +iν mode with
            // its −iν twin and inflating the drift by ~2ν. Rounding collapses the noise so the
            // cluster keeps a stable, index-aligned order.
            var sorted = new Complex[dim];
            for (int i = 0; i < dim; i++) sorted[i] = evals[i];
            Array.Sort(sorted, (a, b) => CompareEigenvalues(a, b));
            sortedEigenvalues[p] = sorted;

            // In-between: the slowCount eigenvalues with smallest |Re λ|, their right eigenvectors
            // orthonormalized (modified Gram-Schmidt) into Q, then the projector P = Q·Qᴴ.
            int k = Math.Min(slowCount, dim);
            var slowIdx = Enumerable.Range(0, dim)
                .OrderBy(i => Math.Abs(evals[i].Real))
                .Take(k)
                .ToArray();
            var slowVecs = Matrix<Complex>.Build.Dense(dim, k);
            for (int j = 0; j < k; j++)
                slowVecs.SetColumn(j, evd.EigenVectors.Column(slowIdx[j]));
            var Q = OrthonormalizeColumns(slowVecs);
            slowBasis[p] = Q;
            projectors[p] = Q * Q.ConjugateTranspose();

            // Polarity ladder reading (F99): sin²(θ)/2, ¼ at the symmetric crossover, ½ at pure YZ.
            double s = Math.Sin(theta);
            polarity[p] = s * s / 2.0;
        }

        // Marks drift: max over θ and over sorted index k of |λ_k(θ) − λ_k(θ_0)|.
        double maxDrift = 0.0;
        var reference = sortedEigenvalues[0];
        for (int p = 1; p < points; p++)
        {
            var sorted = sortedEigenvalues[p];
            for (int i = 0; i < reference.Length; i++)
            {
                double d = (sorted[i] - reference[i]).Magnitude;
                if (d > maxDrift) maxDrift = d;
            }
        }

        // In-between metric (coarse): Frobenius distance between consecutive slow-subspace
        // projectors. Saturates near √(2k) once the subspace has turned a finite amount.
        var subspaceRotation = new double[points - 1];
        for (int p = 0; p < points - 1; p++)
            subspaceRotation[p] = (projectors[p] - projectors[p + 1]).FrobeniusNorm();

        // In-between metric (resolving): the largest principal angle (radians) between the slow
        // subspace at θ₀ and at θ[p], cumulative from θ₀. The principal angles are arccos of the
        // singular values of Q(θ₀)ᴴ·Q(θ[p]); the largest angle is arccos of the smallest singular
        // value. CumulativeRotation[0] = 0 by construction (Q(θ₀)ᴴ·Q(θ₀) has every σ = 1), and
        // every angle is bounded above by π/2.
        //
        // Caveat measured on the N=3 crossover axis: the LARGEST principal angle of a highly
        // degenerate slow manifold saturates near π/2 already at the first nonzero θ-step (some
        // direction in the k-mode set rotates out almost immediately), so here it is no finer an
        // eyepiece than the chordal SubspaceRotation. The smallest (first) principal angle is the
        // one that grows ~linearly for small rotation; the largest is the most pessimistic and
        // saturates first. The exposed SlowBasis lets a downstream reader form whichever principal
        // angle (or the geodesic √Σθ_i²) the question needs.
        var cumulativeRotation = new double[points];
        cumulativeRotation[0] = 0.0;
        for (int p = 1; p < points; p++)
            cumulativeRotation[p] = LargestPrincipalAngle(slowBasis[0], slowBasis[p]);

        return new DimensionSweepResult(
            Theta: axis.Theta,
            Eigenvalues: sortedEigenvalues,
            Polarity: polarity,
            SubspaceRotation: subspaceRotation,
            MaxEigenvalueDriftAcrossTheta: maxDrift,
            SlowBasis: slowBasis,
            CumulativeRotation: cumulativeRotation);
    }

    /// <summary>Order eigenvalues by (Re, then Im) with both keys rounded to
    /// <paramref name="decimals"/> places first. The rounding makes the order stable inside a
    /// degenerate cluster (equal real parts up to the Evd floor), which is what keeps the sorted
    /// θ-to-θ comparison index-aligned and the marks-drift metric honest.</summary>
    private static int CompareEigenvalues(Complex a, Complex b, int decimals = 9)
    {
        int r = Math.Round(a.Real, decimals).CompareTo(Math.Round(b.Real, decimals));
        if (r != 0) return r;
        return Math.Round(a.Imaginary, decimals).CompareTo(Math.Round(b.Imaginary, decimals));
    }

    /// <summary>The largest principal angle (radians) between the column spans of two orthonormal
    /// bases <paramref name="qa"/> and <paramref name="qb"/>. Principal angles θ_i are arccos of the
    /// singular values σ_i of qaᴴ·qb (each σ_i ∈ [0, 1]); the largest angle is arccos of the
    /// smallest σ, clamped into [0, 1] against round-off. Returns a value in [0, π/2].
    ///
    /// <para>Guard: if the two bases have different column counts (slow-mode membership changed
    /// across θ, so the subspaces have different dimensions), the overlap matrix qaᴴ·qb is
    /// rectangular and only min(cols) singular values exist. We compare on that common
    /// min-dimension subspace rather than crash; the largest principal angle is then read from the
    /// min(cols) singular values that the SVD returns.</para></summary>
    private static double LargestPrincipalAngle(ComplexMatrix qa, ComplexMatrix qb)
    {
        int ka = qa.ColumnCount;
        int kb = qb.ColumnCount;
        if (ka == 0 || kb == 0) return 0.0; // degenerate: no shared subspace to rotate

        // Overlap matrix M = qaᴴ·qb (ka × kb). Its singular values are the cosines of the principal
        // angles between the two spans; there are min(ka, kb) of them.
        var overlap = qa.ConjugateTranspose() * qb;
        var singular = overlap.Svd().S; // descending real singular values, length min(ka, kb)

        double smallest = double.PositiveInfinity;
        for (int i = 0; i < singular.Count; i++)
        {
            double s = singular[i].Real;
            if (s < smallest) smallest = s;
        }
        if (double.IsPositiveInfinity(smallest)) return 0.0;

        double cosTheta = Math.Clamp(smallest, 0.0, 1.0);
        return Math.Acos(cosTheta);
    }

    /// <summary>Modified Gram-Schmidt: returns an orthonormal basis (in the standard Hermitian
    /// inner product) for the column span of <paramref name="m"/>. Columns whose residual norm
    /// falls below <paramref name="tol"/> (linearly dependent within numerical precision) are
    /// dropped, so the result may have fewer columns than the input; this keeps the projector
    /// P = Q·Qᴴ well-defined even when the slow eigenvectors are near-degenerate.</summary>
    private static ComplexMatrix OrthonormalizeColumns(ComplexMatrix m, double tol = 1e-10)
    {
        int rows = m.RowCount;
        int cols = m.ColumnCount;
        var kept = new List<ComplexVector>(cols);
        for (int j = 0; j < cols; j++)
        {
            var v = m.Column(j);
            foreach (var q in kept)
            {
                // projection coefficient ⟨q, v⟩ = qᴴ v
                var coeff = q.ConjugateDotProduct(v);
                v = v - coeff * q;
            }
            double norm = v.L2Norm();
            if (norm > tol)
                kept.Add(v / norm);
        }
        var result = Matrix<Complex>.Build.Dense(rows, kept.Count);
        for (int j = 0; j < kept.Count; j++)
            result.SetColumn(j, kept[j]);
        return result;
    }
}
