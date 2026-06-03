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
/// inner products line up with the SlowBasis vectors without a transpose twist.</para></summary>
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
