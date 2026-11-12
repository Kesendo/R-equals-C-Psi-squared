using System.Numerics;
using RCPsiSquared.Core.Numerics;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The lemma of F173's Theorem H (<c>docs/proofs/PROOF_SIGNED_MASK_BALANCE.md</c>): L = −i[J·H, ·] +
/// Σ_l γ_l (Z_l·Z_l − ·) at second order in J on its two edges. On the populations |a⟩⟨a| the 2^N eigenvalues near 0
/// are −J²·spec(𝓛) + O(J³), and on the far cells |a⟩⟨ā| the 2^N near −2σ are −2σ + J²·spec(𝓜) + O(J³), with γ(m) the
/// sum of the rates over the bits of m and
/// <code>
///   𝓛_ba = −|H_ba|²/γ(a ⊕ b),        𝓛_aa = Σ_b |H_ba|²/γ(a ⊕ b),
///   𝓜_ba = −H_ba·H_āb̄/γ(a ⊕ b),      𝓜_aa = Σ_b (|H_ba|² + |H_b̄ā|²)/(2γ(a ⊕ b)).
/// </code>
/// 𝓛 is the Zeno generator of PROOF_WEIGHT1_DEGENERACY's Theorem E (a) on every coset; 𝓜 is, edge by edge, |h||h′|/γ times a
/// magnetic Laplacian's block carrying Theorem 2's ratio as a phase, plus (|h| − |h′|)²/(2γ) on both cells, h = H_ba and
/// h′ = H_āb̄, the potential present exactly where condition (i) fails. Both are Hermitian and positive semidefinite, and their
/// kernels have Theorem 2's two dimensions, the components of Γ_H and the good components of Γ. The diagonal of H must
/// satisfy Theorem 2's condition (ii), H_aa = H_āā, at every bitstring; otherwise the far edge moves at first order.</summary>
public static class StrongDephasingEdges
{
    /// <summary>𝓛 and 𝓜 from H's columns (<see cref="ComplementConnectionGraph.Columns"/>, bit l = site l) and the rates.</summary>
    public static (ComplexMatrix Near, ComplexMatrix Far) Build(int n, Dictionary<int, GaussianInteger>[] h, IReadOnlyList<double> gamma)
    {
        if (gamma.Count != n) throw new ArgumentException($"{gamma.Count} rates for {n} sites", nameof(gamma));
        if (gamma.Any(g => !(g > 0))) throw new ArgumentException("every rate strictly positive", nameof(gamma));
        int dim = 1 << n, full = dim - 1;
        Complex At(int r, int c) { var e = ComplementConnectionGraph.Entry(h, r, c); return new Complex((double)e.Re, (double)e.Im); }
        double Rate(int m) { double s = 0; for (int l = 0; l < n; l++) if (((m >> l) & 1) == 1) s += gamma[l]; return s; }

        var near = ComplexMatrix.Build.Dense(dim, dim);
        var far = ComplexMatrix.Build.Dense(dim, dim);
        for (int a = 0; a < dim; a++)
        {
            if (ComplementConnectionGraph.Entry(h, a, a) != ComplementConnectionGraph.Entry(h, a ^ full, a ^ full))
                throw new ArgumentException($"H_aa ≠ H_āā at a = {a}: the far edge moves at first order", nameof(h));
            // the neighbours of a in Γ: the entries of column a, and the complements of column ā's
            var neighbours = new HashSet<int>(h[a].Keys);
            foreach (int k in h[a ^ full].Keys) neighbours.Add(k ^ full);
            neighbours.Remove(a);
            foreach (int b in neighbours)
            {
                double s = Rate(a ^ b);
                Complex hba = At(b, a), hfar = At(a ^ full, b ^ full);
                double m2 = hba.Magnitude * hba.Magnitude, f2 = hfar.Magnitude * hfar.Magnitude;
                near[b, a] -= m2 / s;
                near[a, a] += m2 / s;
                far[b, a] -= hba * hfar / s;
                far[a, a] += (m2 + f2) / (2 * s);
            }
        }
        return (near, far);
    }

    /// <summary>Corollary H2's cycle: the walk along the flip sets <paramref name="walk"/> (their XOR zero) from every
    /// bitstring x, each a k-cycle of the far edge's graph, read on 𝓜's own edge blocks. Returns the least eigenvalue of
    /// the cycle's part of 𝓜 over all x, the least modulus of a hop on any of them, and the least |arg| of a cycle's
    /// holonomy, the product of −𝓜_ba/|𝓜_ba| along it; null where a walk revisits a bitstring or meets a hop of zero
    /// modulus (an edge where H_ba or H_āb̄ vanishes). Moduli that differ leave the hop |H_ba||H_āb̄|/γ and add the
    /// lemma's potential, which the cycle's part keeps.</summary>
    public static (double LambdaStar, double HopMin, double ThetaMin)? CycleReading(int n, Dictionary<int, GaussianInteger>[] h,
        IReadOnlyList<double> gamma, IReadOnlyList<int> walk)
    {
        int dim = 1 << n, full = dim - 1, k = walk.Count;
        if (k < 3) throw new ArgumentException("a walk of at least three flip sets", nameof(walk));
        if (walk.Aggregate(0, (acc, m) => acc ^ m) != 0) throw new ArgumentException("the walk's flip sets do not close", nameof(walk));
        Complex At(int r, int c) { var e = ComplementConnectionGraph.Entry(h, r, c); return new Complex((double)e.Re, (double)e.Im); }
        double Rate(int m) { double s = 0; for (int l = 0; l < n; l++) if (((m >> l) & 1) == 1) s += gamma[l]; return s; }
        double lambdaStar = double.PositiveInfinity, hopMin = double.PositiveInfinity, thetaMin = double.PositiveInfinity;
        for (int x = 0; x < dim; x++)
        {
            var verts = new int[k];
            verts[0] = x;
            for (int t = 1; t < k; t++) verts[t] = verts[t - 1] ^ walk[t - 1];
            if (verts.Distinct().Count() < k) return null;
            var cycle = ComplexMatrix.Build.Dense(k, k);
            Complex holonomy = Complex.One;
            for (int t = 0; t < k; t++)
            {
                int a = verts[t], b = verts[(t + 1) % k];
                double s = Rate(a ^ b);
                Complex hba = At(b, a), hfar = At(a ^ full, b ^ full);
                Complex off = -hba * hfar / s;
                if (off.Magnitude == 0) return null;
                double diagonal = (hba.Magnitude * hba.Magnitude + hfar.Magnitude * hfar.Magnitude) / (2 * s);
                cycle[(t + 1) % k, t] += off;
                cycle[t, (t + 1) % k] += Complex.Conjugate(off);
                cycle[t, t] += diagonal;
                cycle[(t + 1) % k, (t + 1) % k] += diagonal;
                hopMin = Math.Min(hopMin, off.Magnitude);
                holonomy *= -off / off.Magnitude;
            }
            lambdaStar = Math.Min(lambdaStar, cycle.Evd(MathNet.Numerics.LinearAlgebra.Symmetricity.Hermitian).EigenValues.Min(z => z.Real));
            thetaMin = Math.Min(thetaMin, Math.Abs(holonomy.Phase));
        }
        return (lambdaStar, hopMin, thetaMin);
    }

    /// <summary>Corollary H3's reading: where every flip set's entry H_(x⊕m, x) vanishes nowhere and has one modulus at every
    /// x, on the far side as well (|H_ba| = |H_āb̄|, no potential), 𝓜 = Σ_m μ_m·(1 − U_m) with μ_m = |H_(x⊕m, x)|²/γ(m) and U_m
    /// the hop along m carrying the phase H_ba·H_āb̄/|H_ba|², a Hermitian involution. Where the U_m commute, on each coset x + W
    /// of the flip sets' span the products of the U_m over the relations are signs e_u(x), and the far spectrum there is
    /// {2·Σ_m μ_m·y_m : y with Σ_(m∈u) y_m = e_u(x) on a basis of the relations}, read in Gaussian integers up to the
    /// rates. Where every two of them anticommute, λ_min = Σ μ − (Σ μ²)^(1/2). Null where an entry vanishes or a modulus
    /// differs; <c>Spectrum</c> null where the U_m do not commute, <c>LambdaMin</c> null where they neither commute nor
    /// all anticommute. <c>Split</c>: the signs e_u differ between cosets (under a character ratio, read sums that do not
    /// vanish). The lemma's hypotheses as in <see cref="Build"/>: H_aa = H_āā at every a, every rate strictly positive. At
    /// most 20 flip sets.</summary>
    public static HopReading? ReadHops(int n, Dictionary<int, GaussianInteger>[] h, IReadOnlyList<double> gamma)
    {
        if (gamma.Count != n) throw new ArgumentException($"{gamma.Count} rates for {n} sites", nameof(gamma));
        if (gamma.Any(g => !(g > 0))) throw new ArgumentException("every rate strictly positive", nameof(gamma));
        int dim = 1 << n, full = dim - 1;
        for (int a = 0; a < dim; a++)
            if (ComplementConnectionGraph.Entry(h, a, a) != ComplementConnectionGraph.Entry(h, a ^ full, a ^ full))
                throw new ArgumentException($"H_aa ≠ H_āā at a = {a}: the far edge moves at first order", nameof(h));
        double Rate(int m) { double s = 0; for (int l = 0; l < n; l++) if (((m >> l) & 1) == 1) s += gamma[l]; return s; }
        static BigInteger Norm(GaussianInteger z) => z.Re * z.Re + z.Im * z.Im;
        var flips = new SortedSet<int>();
        for (int a = 0; a < dim; a++) foreach (var (r, e) in h[a]) if (r != a && e != GaussianInteger.Zero) flips.Add(r ^ a);
        var masks = flips.ToArray();
        int f = masks.Length;
        if (f > 20) throw new ArgumentException($"{f} flip sets, at most 20 read", nameof(h));
        // the hop's phase numerator p_m(a) = H_ba·H_āb̄ (b = a ⊕ m), of norm K_m² where K_m = |H_ba|²
        var p = new GaussianInteger[f][];
        var k = new BigInteger[f];
        var mu = new double[f];
        for (int i = 0; i < f; i++)
        {
            int m = masks[i];
            p[i] = new GaussianInteger[dim];
            for (int a = 0; a < dim; a++)
            {
                var hba = ComplementConnectionGraph.Entry(h, a ^ m, a);
                var hfar = ComplementConnectionGraph.Entry(h, a ^ full, a ^ m ^ full);
                BigInteger nb = Norm(hba);
                if (nb.IsZero || Norm(hfar) != nb) return null;
                if (a == 0) k[i] = nb; else if (nb != k[i]) return null;
                p[i][a] = hba * hfar;
            }
            mu[i] = (double)k[i] / Rate(m);
        }
        // U_i U_j = ± U_j U_i at every y: p_i(y)·p_j(y ⊕ m_i) against p_j(y)·p_i(y ⊕ m_j)
        var anti = new bool[f, f];
        bool commuting = true, allAnti = f >= 2;
        for (int i = 0; i < f; i++)
            for (int j = i + 1; j < f; j++)
            {
                bool comm = true, ac = true;
                for (int y = 0; y < dim && (comm || ac); y++)
                {
                    var left = p[i][y] * p[j][y ^ masks[i]];
                    var right = p[j][y] * p[i][y ^ masks[j]];
                    if (left != right) comm = false;
                    if (left != -right) ac = false;
                }
                anti[i, j] = anti[j, i] = ac;
                if (!comm) commuting = false;
                if (!ac) allAnti = false;
            }
        if (!commuting)
        {
            double? floor = allAnti ? mu.Sum() - Math.Sqrt(mu.Sum(x => x * x)) : null;
            return new HopReading(masks, mu, anti, null, floor, false);
        }
        // a basis of the relations among the flip sets (GF(2) elimination with the combinations recorded)
        var relations = new List<int>();
        var pivots = new List<(int V, int Tag)>();
        for (int i = 0; i < f; i++)
        {
            int v = masks[i], tag = 1 << i;
            foreach (var (pv, pt) in pivots) if ((v ^ pv) < v) { v ^= pv; tag ^= pt; }
            if (v == 0) relations.Add(tag);
            else { pivots.Add((v, tag)); pivots.Sort((x, y) => y.V.CompareTo(x.V)); }
        }
        var span = new HashSet<int> { 0 };
        foreach (int m in masks) foreach (int x in span.ToArray()) span.Add(x ^ m);
        var reps = Enumerable.Range(0, dim).Where(x => span.All(s => (x ^ s) >= x)).ToArray();
        var spectrum = new List<double>(dim);
        var signsSeen = new HashSet<int>();
        foreach (int x0 in reps)
        {
            int signs = 0;
            for (int r = 0; r < relations.Count; r++)
            {
                // the walk of the relation from x0, its flip sets in ascending order: Π p over Π K is ±1
                GaussianInteger prod = GaussianInteger.One;
                BigInteger norm = BigInteger.One;
                int y = x0;
                for (int i = 0; i < f; i++)
                    if (((relations[r] >> i) & 1) == 1) { prod *= p[i][y]; norm *= k[i]; y ^= masks[i]; }
                if (prod == new GaussianInteger(-norm, BigInteger.Zero)) signs |= 1 << r;
                else if (prod != new GaussianInteger(norm, BigInteger.Zero))
                    throw new InvalidOperationException($"commuting hops with a relation product {prod} over {norm} at {x0}");
            }
            signsSeen.Add(signs);
            for (int yv = 0; yv < 1 << f; yv++)
            {
                bool ok = true;
                for (int r = 0; r < relations.Count && ok; r++)
                    if ((BitOperations.PopCount((uint)(yv & relations[r])) & 1) != ((signs >> r) & 1)) ok = false;
                if (!ok) continue;
                double s = 0;
                for (int i = 0; i < f; i++) if (((yv >> i) & 1) == 1) s += mu[i];
                spectrum.Add(2 * s);
            }
        }
        if (spectrum.Count != dim) throw new InvalidOperationException($"{spectrum.Count} levels read for {dim} cells");
        spectrum.Sort();
        return new HopReading(masks, mu, anti, spectrum.ToArray(), spectrum[0], signsSeen.Count > 1);
    }
}

/// <summary>Corollary H3's reading (<see cref="StrongDephasingEdges.ReadHops"/>): the flip sets, their hop weights μ_m,
/// which pairs of hops anticommute, the far edge's second-order spectrum where the hops commute, its least eigenvalue
/// where they commute or every two anticommute, and whether the relations' signs differ between cosets.</summary>
public sealed record HopReading(int[] Masks, double[] Mu, bool[,] Anticommute, double[]? Spectrum, double? LambdaMin, bool Split);
