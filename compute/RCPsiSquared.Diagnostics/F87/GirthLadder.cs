using System.Numerics;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>The branch of the girth ladder a windowed diagonal-cell pair sits on
/// (<see cref="GirthLadder.Forecast"/>).</summary>
public enum GirthLadderBranch
{
    /// <summary>The hopping graph has no odd cycle and H has no diagonal lift: every odd power sum
    /// vanishes exactly, the spectrum is symmetric, the pair is SOFT.</summary>
    Bipartite,

    /// <summary>t_ℓ ≠ 0 at some site: the first asymmetric moment is the deg-1 rung,
    /// m* = 2ℓ+1, with the positive monomial p_{m*} = P_{2ℓ+1,1}·γ where
    /// P_{2ℓ+1,1} = (2ℓ+1)·C(2ℓ,ℓ)·Σ_l (t_ℓ^(l))² &gt; 0. A positive monomial has no positive
    /// root, so the pair is HARD at EVERY γ &gt; 0 outright, no residual.</summary>
    Deg1Outright,

    /// <summary>t_ℓ ≡ 0 at every site: the γ¹ class is dead at 2ℓ+1 (and at 2ℓ+3), so the ladder
    /// continues upward: m* = 2ℓ+deg with deg odd ≥ 3 (γ³ when that rung fires, higher when not;
    /// first γ⁵ witness IIXY+ZXZY at m* = 11). The reported m* = 2ℓ+3 is a LOWER BOUND. Pascal-Gram
    /// positivity (F117) still closes hardness at whatever rung fires; only the exact rung needs
    /// the exact engine.</summary>
    HigherRung,
}

/// <summary>The girth-ladder reading of one windowed diagonal-cell pair: where the first
/// asymmetric moment of the recentred Liouvillian M = A + γQ sits, and whether positivity
/// is already explicit at the deg-1 rung.</summary>
/// <param name="Ell">The effective odd-cycle order ℓ: 1 for a diagonal lift, the hopping graph's
/// unsigned odd-girth for a pure cycle, 0 for bipartite.</param>
/// <param name="EllKind">"diagonal-lift", "{ℓ}-cycle" (e.g. "3-cycle"), or "bipartite";
/// mirrors the Python anchor's <c>effective_ell</c> kind strings.</param>
/// <param name="Branch">The ladder branch (see <see cref="GirthLadderBranch"/>).</param>
/// <param name="MStar">The first nonvanishing odd power sum: exact 2ℓ+1 for
/// <see cref="GirthLadderBranch.Deg1Outright"/>; the lower bound 2ℓ+3 for
/// <see cref="GirthLadderBranch.HigherRung"/>; 0 for bipartite (no odd moment ever fires).</param>
/// <param name="MStarIsLowerBound">True only on the <see cref="GirthLadderBranch.HigherRung"/>
/// branch, where the actual rung may sit higher (deg odd ≥ 3).</param>
/// <param name="Coefficient">P_{2ℓ+1,1} = (2ℓ+1)·C(2ℓ,ℓ)·Σ_l (t_ℓ^(l))² for
/// <see cref="GirthLadderBranch.Deg1Outright"/> (strictly positive there); 0 otherwise.</param>
/// <param name="GirthMoments">The per-site girth moments t_ℓ^(l) = Tr(Z_l·H^ℓ), l = 0..N−1
/// (real for Hermitian H); empty for bipartite (ℓ = 0 has no girth level).</param>
public sealed record GirthLadderForecast(
    int Ell,
    string EllKind,
    GirthLadderBranch Branch,
    int MStar,
    bool MStarIsLowerBound,
    double Coefficient,
    IReadOnlyList<double> GirthMoments);

/// <summary>The girth-ladder predictor: for THIS windowed diagonal-cell pair, where does the first
/// asymmetric moment sit, and is it hard at all γ outright? Answered WITHOUT any Liouvillian:
/// only the dense 2^N chain Hamiltonian (built via
/// <see cref="PauliKBodyChainExtensions.ChainKBody"/>, the same path
/// <see cref="PauliPairTrichotomy"/>'s k-body overload uses), its hopping graph, and ℓ dense
/// multiplies. Like the soft certifier, it never touches the 4^N space.
///
/// <para><b>The girth dichotomy</b> (the ladder that retired R-deg, 2026-06-10;
/// <c>simulations/f87_girth_dichotomy.py</c>, anchored on
/// <c>simulations/f87_windowed_monomial_converse.py</c>; proof
/// <c>docs/proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md</c> §4–§5): let ℓ be the effective
/// odd-cycle order of H (a nonzero diagonal ⟹ ℓ = 1 diagonal lift; else the unsigned odd-girth of
/// the basis-state hopping graph; no odd cycle ⟹ bipartite ⟹ soft). The girth moments
/// t_ℓ^(l) = Tr(Z_l·H^ℓ) decide the branch:</para>
///
/// <list type="bullet">
///   <item>t_ℓ ≠ 0 somewhere ⟹ m* = 2ℓ+1, deg 1, and the first odd power sum is the positive
///   monomial P_{2ℓ+1,1}·γ with P_{2ℓ+1,1} = (2ℓ+1)·C(2ℓ,ℓ)·Σ_l (t_ℓ^(l))² &gt; 0. A positive
///   monomial has no positive root ⟹ HARD at EVERY γ &gt; 0 outright, no residual. The
///   single-site-Z lift is the ℓ = 1 face: t_1^(l) = 2^N·c_l reproduces
///   P_{3,1} = 6·4^N·Σ_l c_l² bit-exactly.</item>
///   <item>t_ℓ ≡ 0 everywhere ⟹ the deg-1 rung is dead and m* = 2ℓ+deg with deg odd ≥ 3 (the
///   forecast reports the lower bound 2ℓ+3). Hardness at the firing rung is closed by the
///   Pascal-Gram positivity theorem (F117, <c>simulations/f87_pascal_gram_positivity.py</c>):
///   every #Q class at m* is a sum of squares or exactly zero, so hard at one γ ⟹ hard at ALL γ.
///   The typed carrier of the closed converse is <see cref="WindowedConverseAllGammaClaim"/>
///   (Tier1Derived, no residual).</item>
/// </list>
/// </summary>
public static class GirthLadder
{
    /// <summary>Entry tolerance for the dense H: a diagonal entry or hopping-graph edge exists
    /// where the magnitude exceeds this. The windowed chain builders produce exact small integers
    /// (times the template coefficients), so this only absorbs float round-off.</summary>
    public const double EntryTolerance = 1e-12;

    /// <summary>Branch tolerance on the girth moments: t_ℓ counts closed walks (integers for
    /// integer-coefficient templates), so |t_ℓ| &gt; 1e-9 cleanly separates fired from dead.</summary>
    public const double MomentTolerance = 1e-9;

    /// <summary>Forecast the girth-ladder branch for the k-body template set placed on an N-site
    /// chain (sliding-window sum, the <see cref="PauliKBodyChainExtensions.ChainKBody"/> input
    /// shape shared with <see cref="PauliPairTrichotomy"/>'s k-body overload).</summary>
    /// <param name="termTemplates">The templates (each <see cref="PauliTerm.Letters"/> of length
    /// k ≤ n plus a coefficient); mixed spans are allowed, exactly as in ChainKBody.</param>
    /// <param name="n">The chain length N; the dense H is 2^N × 2^N.</param>
    public static GirthLadderForecast Forecast(IReadOnlyList<PauliTerm> termTemplates, int n)
    {
        if (termTemplates is null) throw new ArgumentNullException(nameof(termTemplates));

        var H = termTemplates.ChainKBody(n);

        // Effective ℓ: a diagonal lift short-circuits to ℓ = 1 (no hopping cycle needed);
        // otherwise the unsigned odd-girth of the basis-state hopping graph.
        int ell;
        string kind;
        if (HasNonzeroDiagonal(H))
        {
            ell = 1;
            kind = "diagonal-lift";
        }
        else
        {
            ell = ShortestOddCycle(BuildHoppingAdjacency(H));
            if (ell == 0)
                return new GirthLadderForecast(0, "bipartite", GirthLadderBranch.Bipartite,
                    MStar: 0, MStarIsLowerBound: false, Coefficient: 0.0,
                    GirthMoments: Array.Empty<double>());
            kind = $"{ell}-cycle";
        }

        double[] moments = ComputeGirthMoments(H, n, ell);
        bool fires = moments.Any(t => Math.Abs(t) > MomentTolerance);

        if (fires)
        {
            double sumSq = moments.Sum(t => t * t);
            double coefficient = (2 * ell + 1) * CentralBinomial(ell) * sumSq;
            return new GirthLadderForecast(ell, kind, GirthLadderBranch.Deg1Outright,
                MStar: 2 * ell + 1, MStarIsLowerBound: false, coefficient, moments);
        }

        return new GirthLadderForecast(ell, kind, GirthLadderBranch.HigherRung,
            MStar: 2 * ell + 3, MStarIsLowerBound: true, Coefficient: 0.0, moments);
    }

    private static bool HasNonzeroDiagonal(ComplexMatrix H)
    {
        for (int i = 0; i < H.RowCount; i++)
            if (H[i, i].Magnitude > EntryTolerance)
                return true;
        return false;
    }

    /// <summary>Adjacency lists of the basis-state hopping graph: an undirected edge wherever
    /// |H[i,j]| &gt; tol with i ≠ j (H Hermitian, so the dense scan is symmetric).</summary>
    private static List<int>[] BuildHoppingAdjacency(ComplexMatrix H)
    {
        int d = H.RowCount;
        var adj = new List<int>[d];
        for (int i = 0; i < d; i++) adj[i] = new List<int>();
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                if (i != j && H[i, j].Magnitude > EntryTolerance)
                    adj[i].Add(j);
        return adj;
    }

    /// <summary>The unsigned odd-girth: BFS from every start vertex s; every edge (u, v) whose
    /// endpoints share a BFS layer (dist_s(u) == dist_s(v)) closes an odd walk of length
    /// 2·dist_s(u) + 1 through s, and the global minimum over s and edges is exactly the shortest
    /// odd cycle (for any shortest odd cycle C and any s on C, BFS from s exposes such an edge at
    /// length |C|). Returns 0 if no odd cycle exists (bipartite).</summary>
    private static int ShortestOddCycle(List<int>[] adj)
    {
        int d = adj.Length;
        int best = int.MaxValue;
        var dist = new int[d];
        var queue = new Queue<int>();
        for (int s = 0; s < d; s++)
        {
            Array.Fill(dist, -1);
            dist[s] = 0;
            queue.Clear();
            queue.Enqueue(s);
            while (queue.Count > 0)
            {
                int u = queue.Dequeue();
                foreach (int v in adj[u])
                {
                    if (dist[v] < 0)
                    {
                        dist[v] = dist[u] + 1;
                        queue.Enqueue(v);
                    }
                    else if (dist[v] == dist[u])
                    {
                        int candidate = dist[u] + dist[v] + 1;
                        if (candidate < best) best = candidate;
                    }
                }
            }
        }
        return best == int.MaxValue ? 0 : best;
    }

    /// <summary>The per-site girth moments t_ℓ^(l) = Tr(Z_l·H^ℓ). Z_l is diagonal ±1, so each
    /// trace is a sign-weighted diagonal sum of H^ℓ (ℓ dense multiplies, cheap at 2^N). The signs
    /// are read from <see cref="PauliString.SiteOp"/>'s diagonal to stay bit-convention-safe.
    /// For Hermitian H every t_ℓ is real (trace of a product of two Hermitian matrices); the
    /// imaginary part is asserted below round-off.</summary>
    private static double[] ComputeGirthMoments(ComplexMatrix H, int n, int ell)
    {
        var hPow = H;
        for (int p = 1; p < ell; p++) hPow *= H;

        int d = hPow.RowCount;
        var moments = new double[n];
        for (int l = 0; l < n; l++)
        {
            var zl = PauliString.SiteOp(n, l, PauliLetter.Z);
            Complex t = Complex.Zero;
            for (int i = 0; i < d; i++)
                t += zl[i, i] * hPow[i, i];
            if (Math.Abs(t.Imaginary) > 1e-9)
                throw new InvalidOperationException(
                    $"girth moment t_{ell}^({l}) has a nonzero imaginary part {t.Imaginary:E2}; " +
                    "is H Hermitian (real template coefficients)?");
            moments[l] = t.Real;
        }
        return moments;
    }

    /// <summary>C(2ℓ, ℓ) by the multiplicative formula (every prefix is the integer C(ℓ+i, i),
    /// so the integer division is exact); exact in long for every ℓ this predictor meets.</summary>
    private static double CentralBinomial(int ell)
    {
        long c = 1;
        for (int i = 1; i <= ell; i++)
            c = c * (ell + i) / i;
        return c;
    }
}
