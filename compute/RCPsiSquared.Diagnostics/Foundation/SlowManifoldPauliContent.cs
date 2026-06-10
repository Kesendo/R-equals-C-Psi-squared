using System.Numerics;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The Pauli reading of the fan split: it takes the two slow-subspace frames Q(θ₀) and Q(θ)
/// the <see cref="DimensionSweep"/> exposes, separates them into the invariant core (principal angle
/// ≈ 0) and the rotating remainder, and projects each onto the 4^N Pauli basis to name what they are.
///
/// <para>On the crossover axis the prediction is sharp. The rotation Ad_{R_z} on the lit sites fixes
/// exactly the operators that read {I, Z} on those sites and turns the ones carrying an X or a Y
/// there. So the core should sit entirely on {I, Z}-on-lit Pauli strings (lit-XY weight ≈ 0) and the
/// rotating part should carry lit-XY content. This is the operator-algebra face of the PTF banks: the
/// {I, Z} shadow is the settled near bank the turn cannot move, the {X, Y} light is the far bank it
/// turns.</para>
///
/// <para>The projection is done in the row-major vec convention the dissipator and SlowBasis use
/// (vec(σ)[a·d + b] = σ[a, b]), built string by string here rather than through
/// <see cref="PauliBasis.VecToPauliBasisTransform"/> (whose columns are column-major vec_F), so the
/// inner products line up with the SlowBasis vectors without a transpose twist.</para>
///
/// <para>Two core reads live here: <see cref="Compute"/> is the angle-thresholded two-θ read (the
/// fan's eyepiece, with its slowCount-window caveat), and <see cref="IntersectionCoreRank"/> is the
/// threshold-free single-θ read, dim(slow ∩ fixed cell), which needs neither an angle tolerance nor
/// a θ-pair and retires the window artifact.</para></summary>
public static class SlowManifoldPauliContent
{
    /// <summary>One Pauli string's share of a subspace's mass, with the lit-XY flag (does it carry an
    /// X or a Y on a lit site, the directions the turn rotates).</summary>
    public sealed record StringWeight(string Label, double Weight, bool LitXY);

    /// <summary>The reading: the dimensions of the split, the lit-XY mass fraction of each half (the
    /// headline test of the {I, Z} / {X, Y} attribution), and the dominant Pauli strings of each.</summary>
    public sealed record Reading(
        int CoreDim,
        int RotatingDim,
        double CoreLitXYWeight,
        double RotatingLitXYWeight,
        IReadOnlyList<StringWeight> CoreTop,
        IReadOnlyList<StringWeight> RotatingTop);

    /// <summary>The threshold-free core: dim(slow ∩ fixed cell), read from ONE subspace at ONE θ.
    ///
    /// <para>The fixed cell C is the joint-fixed cell of the per-site rotation characters on the lit
    /// sites: coherences x = a·d + b with Δ_l(x) = 0 for every lit l (bra and ket agree on every lit
    /// bit), equivalently the span of Pauli strings reading {I, Z} on every lit site. Ad_{R_z(θ)} is
    /// DIAGONAL in the coherence basis and acts on x as the phase exp(iθ·Σ_lit(bit_l(a) − bit_l(b))),
    /// so it is the IDENTITY exactly on C: every direction in C is genuinely rotation-fixed, no
    /// principal-angle threshold and no θ-pair needed. The projector Π_C is diagonal (Π_C(x) = 1 iff
    /// x ∈ C), and the core is rank(Π_C·P_slow·Π_C), counted through the projected Gram matrix
    /// (Π_C·Q)ᴴ·(Π_C·Q): its eigenvalues are the squared singular values of the cell-restricted slow
    /// basis, and the count at 1 is the rank (the pinned N = 3 spectrum has the near-1 cluster
    /// separated from the next singular value by ≥ 0.05, so <paramref name="rankTol"/> is a numerical
    /// rank cut, not a lens).</para>
    ///
    /// <para>This retires the angle-core's slowCount-window artifact (see the dated note at
    /// <see cref="DimensionField"/>'s core split): the angle count inflated 4 → 18 → 27 at slowCount
    /// 16 → 24 → 32 because the θ₀ and θ windows re-include each other's rotated images; the
    /// intersection core counts only genuinely fixed directions. Pass the cluster-closed basis
    /// (<see cref="DimensionSweepResult.ClusterClosedBasis"/>): on it the read is exact, θ-stable,
    /// and window-free in the meaningful sense, every requested slowCount that closes to the same
    /// manifold reads the same rank (N = 3 crossover: requested 8, 16, 24 all close to the 28-dim
    /// manifold and read 10; requested 32 closes to 36 and honestly reads 12, real fixed content
    /// that became slow, not a lens artifact; the historical core 4 is the closed 6-dim manifold's
    /// read). On a RAW window that slices a degenerate rate cluster the rank inherits the
    /// membership gauge and can flicker (pinned 2026-06-10: raw slowCount = 16 reads 4 at twelve of
    /// thirteen θ and 5 at θ = π/4).</para>
    ///
    /// <para>With empty <paramref name="litSites"/> the cell is the whole coherence space and the
    /// rank is just the basis dimension (the J-defect axis has no lit/shadow split).</para>
    ///
    /// <para><b>Two quarters, two genealogies (numerological guard, 2026-06-10).</b> The window
    /// ratios of the ranks (4/6, 10/28, 12/36 at N = 3) are NOT a law: at N = 4 the corresponding
    /// reads are 6/10, 8/12, 10/24, 12/32, no pattern recurs (the tempting 12/36 = 1/3 was a
    /// cluster-size coincidence). What IS a law is the fixed CELL's dimension fraction:
    /// (1/2)^(#lit), the structural half per lit site ({I,Z} is half of the four letters), which
    /// lands on ¼ at exactly two lit sites (N = 3) and ⅛ at three (N = 4). That ¼ is a COUSIN of
    /// the cusp horizon's CΨ = ¼ through the shared root ½ (the half-of-half node,
    /// HalfAsStructuralFixedPointClaim / QuarterAsBilinearMaxvalClaim), not the same object: the
    /// horizon's quarter is the bilinear apex p(1−p) ≤ ¼ and N-free; the cell fraction halves per
    /// lit site. Same number at N = 3, two family trees.</para></summary>
    public static int IntersectionCoreRank(ComplexMatrix slowBasis, int N, IReadOnlyList<int> litSites,
        double rankTol = 1e-6)
    {
        int k = slowBasis.ColumnCount;
        if (k == 0) return 0;
        int d = 1 << N;
        int dim = slowBasis.RowCount; // d²
        if (dim != d * d)
            throw new ArgumentException($"slowBasis must have 4^N = {d * d} rows for N={N}, got {dim}", nameof(slowBasis));

        // The fixed cell's rows: bra and ket agree on every lit bit (site l ↔ bit N−1−l).
        var rows = new List<int>(dim);
        for (int x = 0; x < dim; x++)
        {
            int diff = (x / d) ^ (x % d);
            bool inCell = true;
            foreach (int l in litSites)
                if (((diff >> (N - 1 - l)) & 1) != 0) { inCell = false; break; }
            if (inCell) rows.Add(x);
        }

        // The projected Gram matrix G = (Π_C·Q)ᴴ·(Π_C·Q), k × k Hermitian PSD; its eigenvalues are
        // the σ² of the cell-restricted basis, and an eigenvalue at 1 ⟺ a slow direction lies
        // entirely in the cell. Counting them is rank(Π_C·P_slow·Π_C). (Gram + Hermitian Evd rather
        // than an SVD of the restricted matrix: the restriction is wide, |C| < k, on the closed
        // windows, and the Gram route is shape-independent.)
        var restricted = ComplexMatrix.Build.Dense(rows.Count, k, (r, j) => slowBasis[rows[r], j]);
        var gram = restricted.ConjugateTranspose() * restricted;
        var eigenvalues = gram.Evd(MathNet.Numerics.LinearAlgebra.Symmetricity.Hermitian).EigenValues;
        double cut = (1.0 - rankTol) * (1.0 - rankTol);
        int rank = 0;
        for (int i = 0; i < eigenvalues.Count; i++)
            if (eigenvalues[i].Real >= cut) rank++;
        return rank;
    }

    /// <summary>Split Q(θ₀) versus Q(θ) into the core (principal angle &lt; <paramref name="angleTolDeg"/>)
    /// and the rotating remainder, and project both onto the Pauli basis. <paramref name="litSites"/>
    /// are the sites the axis rotates (for the crossover, the X/Y carriers). Returns the lit-XY mass
    /// fraction of each half and their top strings by weight.</summary>
    public static Reading Compute(ComplexMatrix q0, ComplexMatrix qTheta, int N,
        IReadOnlyList<int> litSites, double angleTolDeg = 1.0, int topCount = 6)
    {
        int d = 1 << N;
        long d2 = 1L << (2 * N);

        // Principal directions in Q(θ₀): the SVD of the overlap gives the cosines of the principal
        // angles (descending), so the leading columns of q0·U are the least-rotated (the core), the
        // trailing the most-rotated (the in-between). q0 and U are orthonormal/unitary, so the columns
        // of `principal` are an orthonormal basis aligned with the rotation.
        var svd = (q0.ConjugateTranspose() * qTheta).Svd();
        var principal = q0 * svd.U; // d² × k, columns ordered by singular value (cos angle) descending
        int k = principal.ColumnCount;

        double cosTol = Math.Cos(angleTolDeg * Math.PI / 180.0);
        int coreDim = 0;
        for (int j = 0; j < k; j++) if (svd.S[j].Real >= cosTol) coreDim++;
        int rotDim = k - coreDim;

        // The Pauli basis in the row-major vec convention: e_α = vec(σ_α)/√d (orthonormal, since
        // Tr(σ_α† σ_β) = d·δ), plus the lit-XY flag (σ_α carries an X or a Y on a lit site).
        var e = new ComplexVector[d2];
        var label = new string[d2];
        var litXY = new bool[d2];
        double invSqrtD = 1.0 / Math.Sqrt(d);
        var litSet = new HashSet<int>(litSites);
        for (long a = 0; a < d2; a++)
        {
            var letters = PauliIndex.FromFlat(a, N);
            var sigma = PauliString.Build(letters);
            var vec = ComplexVector.Build.Dense((int)d2);
            for (int r = 0; r < d; r++)
                for (int c = 0; c < d; c++)
                    vec[r * d + c] = sigma[r, c] * invSqrtD;
            e[a] = vec;
            label[a] = PauliLabel.Format(letters);

            bool lit = false;
            for (int site = 0; site < N; site++)
                if (litSet.Contains(site) && letters[site].BitA() == 1) { lit = true; break; }
            litXY[a] = lit;
        }

        var coreTop = ProjectSubspace(principal, 0, coreDim, e, label, litXY, d2, out double coreLitXY, topCount);
        var rotTop = ProjectSubspace(principal, coreDim, rotDim, e, label, litXY, d2, out double rotLitXY, topCount);

        return new Reading(coreDim, rotDim, coreLitXY, rotLitXY, coreTop, rotTop);
    }

    /// <summary>Project the subspace spanned by principal columns [<paramref name="start"/>,
    /// start + <paramref name="dim"/>) onto the Pauli basis: per string, the summed |⟨e_α, col⟩|² over
    /// the columns, divided by dim so the fractions sum to 1 (the columns are orthonormal and the e_α
    /// an orthonormal basis). <paramref name="litXYWeight"/> is set to the fraction on lit-XY strings;
    /// returns the top strings by weight.</summary>
    private static IReadOnlyList<StringWeight> ProjectSubspace(ComplexMatrix principal, int start, int dim,
        ComplexVector[] e, string[] label, bool[] litXY, long d2, out double litXYWeight, int topCount)
    {
        var cols = new ComplexVector[dim];
        for (int j = 0; j < dim; j++) cols[j] = principal.Column(start + j);

        var weight = new double[d2];
        double litSum = 0.0;
        for (long a = 0; a < d2; a++)
        {
            double w = 0.0;
            for (int j = 0; j < dim; j++)
            {
                Complex ip = e[a].ConjugateDotProduct(cols[j]); // ⟨e_α, col⟩
                w += ip.Real * ip.Real + ip.Imaginary * ip.Imaginary;
            }
            weight[a] = w;
            if (litXY[a]) litSum += w;
        }
        litXYWeight = dim > 0 ? litSum / dim : 0.0;

        double inv = dim > 0 ? 1.0 / dim : 0.0;
        var top = Enumerable.Range(0, (int)d2)
            .Where(a => weight[a] > 1e-9)
            .OrderByDescending(a => weight[a])
            .Take(topCount)
            .Select(a => new StringWeight(label[a], weight[a] * inv, litXY[a]))
            .ToList();
        return top;
    }
}
