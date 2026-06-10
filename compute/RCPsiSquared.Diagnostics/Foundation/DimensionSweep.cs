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
/// <para>The eyepiece progression. <see cref="DimensionSweepResult.SubspaceRotation"/> is the
/// coarsest: the chordal (Frobenius) projector distance saturates near its ceiling √(2k) once the
/// subspace has turned a finite amount. <see cref="DimensionSweepResult.CumulativeRotation"/>, the
/// largest principal angle between the slow subspace at θ₀ and at θ[p], saturates almost as fast on
/// a degenerate manifold, because the largest angle reports only the single fastest-rotating
/// direction. The resolving eyepiece is the full
/// <see cref="DimensionSweepResult.PrincipalAngleSpectrum"/>: all k principal angles per θ, sorted
/// ascending. It shows the fan, the directions whose angle stays ≈ 0 are the invariant core (itself
/// more contract, on the crossover axis the {I, Z} shadow that R_z fixes), the directions whose
/// angle grows are the true in-between (the {X, Y} lit part that R_z turns). The slow manifold
/// carries its own marks-and-in-between split, and the spectrum makes it visible. All readings are
/// read from the singular values of Q(θ₀)ᴴ·Q(θ[p]) via <see cref="DimensionSweepResult.SlowBasis"/>,
/// the per-θ orthonormal Q.</para>
///
/// <para>The light coordinate. <see cref="DimensionSweepResult.PerSiteLight"/> is the slow
/// manifold's per-site light profile at each θ: light_l = Tr(Π_V·Δ_l)/k, the trace of the
/// orthonormal slow projector against the diagonal mask Δ_l(x) = 1 iff bit l differs between bra
/// and ket of the coherence index x = a·d + b (the same Δ_l as the Absorption Theorem and
/// <see cref="Ptf.SlowLightDistribution"/>). The coordinate is gauge-free three times over: no
/// eigenvector phase or intra-degenerate basis enters (it is a projector trace), no vec
/// convention enters (the popcount of a⊕b is bra/ket symmetric), and no window-membership gauge
/// enters, because the profile is read on <see cref="DimensionSweepResult.ClusterClosedBasis"/>,
/// the slowCount window extended to the rate-cluster boundary. The closure is what makes the
/// coordinate honest: a cut INSIDE a degenerate |Re λ| cluster selects a solver-arbitrary slice
/// of that cluster (pinned 2026-06-10 at the N = 3 crossover, slowCount = 16: the raw-window
/// profile drifts 4.9·10⁻² across θ, while the closed window, 28 modes, is flat to 1.7·10⁻¹⁶,
/// profile exactly [11/28, 11/28, 2/28]). On the crossover axis the flatness is a theorem, not a
/// numerical accident: Ad_{R_z(θ)} is diagonal on coherence space, so it commutes with every
/// Δ_l, and the slow manifold's per-site light profile is EXACTLY θ-invariant even as the
/// manifold itself rotates; the turn moves the light's carriers, never its per-site profile.
/// <see cref="DimensionSweepResult.MaxLightDriftAcrossTheta"/> is the flatness gate. On a
/// non-similarity axis (the J-defect) the same number reads the REAL light redistribution.</para>
/// </summary>
public sealed record DimensionSweepResult(
    IReadOnlyList<double> Theta,
    IReadOnlyList<Complex[]> Eigenvalues,
    IReadOnlyList<double> Polarity,
    IReadOnlyList<double> SubspaceRotation,
    double MaxEigenvalueDriftAcrossTheta,
    IReadOnlyList<ComplexMatrix> SlowBasis,
    IReadOnlyList<double> CumulativeRotation,
    IReadOnlyList<double[]> PrincipalAngleSpectrum,
    IReadOnlyList<ComplexMatrix> ClusterClosedBasis,
    IReadOnlyList<double[]> PerSiteLight,
    double MaxLightDriftAcrossTheta);

public static class DimensionSweep
{
    /// <summary>Eigenvalues whose |Re λ| differ by less than this belong to the same rate cluster.
    /// The cluster-closed window extends the slowCount cut up to the next boundary where the gap
    /// exceeds it. Same standing as <see cref="Ptf.SlowLightDistribution.DefaultClusterTolerance"/>:
    /// the dense-Evd floor sits at ~10⁻¹⁴ and the real rate gaps at N ≤ 6 are O(γ), so anything in
    /// between works.</summary>
    public const double ClusterTolerance = 1e-6;

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
        var closedBasis = new ComplexMatrix[points];      // the cluster-closed window, for the gauge-free light read
        var perSiteLight = new double[points][];
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
            var byRate = Enumerable.Range(0, dim)
                .OrderBy(i => Math.Abs(evals[i].Real))
                .ToArray();
            var slowVecs = Matrix<Complex>.Build.Dense(dim, k);
            for (int j = 0; j < k; j++)
                slowVecs.SetColumn(j, evd.EigenVectors.Column(byRate[j]));
            var Q = OrthonormalizeColumns(slowVecs);
            slowBasis[p] = Q;
            projectors[p] = Q * Q.ConjugateTranspose();

            // The cluster-closed window: extend the cut to the rate-cluster boundary, so the
            // window is a union of full |Re λ| clusters. A cut INSIDE a degenerate cluster picks a
            // solver-arbitrary slice (which members survive the cut depends on the Evd's internal
            // ordering, a membership gauge); the closed window is the smallest well-defined slow
            // manifold containing the request, and on a similarity axis it is the exact rotated
            // image of its θ₀ self. The light coordinate is read here, never on the raw slice.
            int kk = k;
            while (kk < dim && Math.Abs(Math.Abs(evals[byRate[kk]].Real) - Math.Abs(evals[byRate[kk - 1]].Real)) < ClusterTolerance)
                kk++;
            ComplexMatrix closed;
            if (kk == k)
            {
                closed = Q;
            }
            else
            {
                var closedVecs = Matrix<Complex>.Build.Dense(dim, kk);
                for (int j = 0; j < kk; j++)
                    closedVecs.SetColumn(j, evd.EigenVectors.Column(byRate[j]));
                closed = OrthonormalizeColumns(closedVecs);
            }
            closedBasis[p] = closed;
            perSiteLight[p] = PerSiteLightOf(closed, axis.N);

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

        // In-between (the full eyepiece): the principal-angle SPECTRUM between the slow subspace at
        // θ₀ and at θ[p], all min(k) angles (radians) sorted ascending. The angles are arccos of the
        // singular values of Q(θ₀)ᴴ·Q(θ[p]). The angles staying near 0 are the invariant core (the
        // part of the slow manifold the rotation fixes; on the crossover axis the {I,Z} shadow R_z
        // leaves alone); the angles growing are the rotating directions (the {X,Y} lit in-between
        // R_z turns). The slow manifold carries its own marks-and-in-between split; the largest
        // angle alone hides it (it saturates at the first step), the spectrum resolves it.
        // CumulativeRotation keeps the largest angle (the last entry, ascending) for the summary.
        var principalSpectrum = new double[points][];
        var cumulativeRotation = new double[points];
        principalSpectrum[0] = new double[slowBasis[0].ColumnCount]; // θ₀ vs itself: every angle 0
        cumulativeRotation[0] = 0.0;
        for (int p = 1; p < points; p++)
        {
            var angles = PrincipalAngles(slowBasis[0], slowBasis[p]); // ascending, each in [0, π/2]
            principalSpectrum[p] = angles;
            cumulativeRotation[p] = angles.Length == 0 ? 0.0 : angles[^1]; // largest = last
        }

        // Light flatness gate: max over θ and over site l of |light_l(θ) − light_l(θ₀)|. On the
        // crossover axis Ad_{R_z(θ)} commutes with every Δ_l, so this is exactly zero up to the
        // Evd floor (pinned 1.7·10⁻¹⁶ at N = 3); on a real perturbation axis it reads the real
        // light redistribution.
        double maxLightDrift = 0.0;
        var lightRef = perSiteLight[0];
        for (int p = 1; p < points; p++)
            for (int l = 0; l < lightRef.Length; l++)
            {
                double dl = Math.Abs(perSiteLight[p][l] - lightRef[l]);
                if (dl > maxLightDrift) maxLightDrift = dl;
            }

        return new DimensionSweepResult(
            Theta: axis.Theta,
            Eigenvalues: sortedEigenvalues,
            Polarity: polarity,
            SubspaceRotation: subspaceRotation,
            MaxEigenvalueDriftAcrossTheta: maxDrift,
            SlowBasis: slowBasis,
            CumulativeRotation: cumulativeRotation,
            PrincipalAngleSpectrum: principalSpectrum,
            ClusterClosedBasis: closedBasis,
            PerSiteLight: perSiteLight,
            MaxLightDriftAcrossTheta: maxLightDrift);
    }

    /// <summary>The per-site light profile of the subspace spanned by the orthonormal columns of
    /// <paramref name="q"/>: light_l = Tr(Π_V·Δ_l)/k = Σ_x Δ_l(x)·Σ_j |q[x,j]|² / k, where
    /// Δ_l(x) = 1 iff bit l differs between bra a = x/d and ket b = x mod d of the coherence
    /// index (site l ↔ bit n−1−l, the <see cref="RCPsiSquared.Core.Lindblad.PauliDephasingDissipator"/>
    /// convention). Gauge-free: a projector trace, so no eigenvector phase or intra-degenerate
    /// basis choice enters, and a⊕b is bra/ket symmetric, so no vec convention enters either.
    /// Each light_l ∈ [0, 1].</summary>
    private static double[] PerSiteLightOf(ComplexMatrix q, int n)
    {
        var light = new double[n];
        int k = q.ColumnCount;
        if (k == 0) return light;
        int d = 1 << n;
        int dim = q.RowCount; // d²
        for (int x = 0; x < dim; x++)
        {
            double w = 0.0;
            for (int j = 0; j < k; j++)
            {
                var c = q[x, j];
                w += c.Real * c.Real + c.Imaginary * c.Imaginary;
            }
            if (w == 0.0) continue;
            int diff = (x / d) ^ (x % d);
            for (int l = 0; l < n; l++)
                if (((diff >> (n - 1 - l)) & 1) != 0) light[l] += w;
        }
        for (int l = 0; l < n; l++) light[l] /= k;
        return light;
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

    /// <summary>All principal angles (radians, sorted ascending) between the column spans of two
    /// orthonormal bases <paramref name="qa"/> and <paramref name="qb"/>. Principal angles θ_i are
    /// arccos of the singular values σ_i of qaᴴ·qb (each σ_i ∈ [0, 1]), clamped against round-off;
    /// there are min(cols) of them. Sorted ascending, the leading entries near 0 are the invariant
    /// directions (the core the rotation fixes) and the trailing entries near π/2 are the rotated-out
    /// directions (the in-between).
    ///
    /// <para>Guard: if the two bases have different column counts (slow-mode membership changed
    /// across θ), the overlap qaᴴ·qb is rectangular and only min(cols) singular values exist; we read
    /// those rather than crash.</para></summary>
    private static double[] PrincipalAngles(ComplexMatrix qa, ComplexMatrix qb)
    {
        int ka = qa.ColumnCount;
        int kb = qb.ColumnCount;
        if (ka == 0 || kb == 0) return Array.Empty<double>(); // no shared subspace to rotate

        // Overlap M = qaᴴ·qb (ka × kb); its singular values are the cosines of the principal angles.
        var overlap = qa.ConjugateTranspose() * qb;
        var singular = overlap.Svd().S; // descending real singular values, length min(ka, kb)

        var angles = new double[singular.Count];
        for (int i = 0; i < singular.Count; i++)
            angles[i] = Math.Acos(Math.Clamp(singular[i].Real, 0.0, 1.0));
        Array.Sort(angles); // ascending: invariant core first, rotated-out last
        return angles;
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
