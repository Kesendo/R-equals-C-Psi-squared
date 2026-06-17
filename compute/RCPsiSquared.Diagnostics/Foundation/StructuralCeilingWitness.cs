using System;
using System.Collections.Generic;
using System.Globalization;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The structural ceiling computed live (typed home: <c>StructuralCeilingClaim</c>). The high-Q gap
/// rate <c>g2 = strict_gap / 2γ</c> of an XY network under uniform Z-dephasing saturates at a topology-specific
/// ceiling below 1 when a dark <c>[H,A]=0</c> coherence undercuts the single-excitation band edge. This witness
/// recomputes the ceiling from the representation structure (not by building the 4^N Liouvillian): in the
/// |a⟩⟨b| coherence basis N_XY is diagonal with entry hamming(a,b), and the high-Q ceiling is the smallest
/// nonzero eigenvalue of N_XY restricted to the H-commutant (the ad_H kernel) of a small excitation sector.
/// J-independent (the high-Q limit), so the witness takes no J/γ.
///
/// <para>It confirms, live and exact: complete <c>g2(K_N) = 4/N</c> (N≥5, the (1,1) S_N-standard-rep sector);
/// star <c>g2(star_N) = 4/(N−1)</c> (N≥6, the (1,1) leaf manifold); and the N=4 outlier, the (2,2)
/// half-filling sector shared by K_4 (= 2 − 2/√3, below the floor) and ring-4 (= 1, co-occupied).</para></summary>
public sealed class StructuralCeilingWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private const double NullTol = 1e-7;   // |ad_H eigenvalue| < tol  ⇒  in the H-commutant (Ω=0)
    private const double NzTol = 1e-7;     // N_XY eigenvalue > tol     ⇒  non-steady (n_XY ≠ 0)

    // ---- topology + sector linear algebra (a faithful port of topology_ceiling_rep_derivation.py) ----

    private static List<(int, int)> Bonds(string topo, int n)
    {
        var b = new List<(int, int)>();
        switch (topo)
        {
            case "chain": for (int i = 0; i < n - 1; i++) b.Add((i, i + 1)); break;
            case "ring": for (int i = 0; i < n; i++) b.Add((i, (i + 1) % n)); break;
            case "star": for (int i = 1; i < n; i++) b.Add((0, i)); break;
            case "complete": for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) b.Add((i, j)); break;
            default: throw new ArgumentException($"unknown topology '{topo}'", nameof(topo));
        }
        return b;
    }

    private static int[] SectorStates(int n, int p)
    {
        var list = new List<int>();
        for (int x = 0; x < (1 << n); x++)
            if (BitOperations.PopCount((uint)x) == p) list.Add(x);
        return list.ToArray();
    }

    /// <summary>The XY hopping Hamiltonian (amplitude 1) restricted to the popcount-p Hilbert sector,
    /// real symmetric: an excitation hops across a bond between two sites that disagree.</summary>
    private static Matrix<double> SectorH(string topo, int n, int[] states)
    {
        int m = states.Length;
        var idx = new Dictionary<int, int>(m);
        for (int a = 0; a < m; a++) idx[states[a]] = a;
        var H = Matrix<double>.Build.Dense(m, m);
        var bonds = Bonds(topo, n);
        for (int a = 0; a < m; a++)
        {
            int s = states[a];
            foreach (var (i, j) in bonds)
            {
                int bi = (s >> i) & 1, bj = (s >> j) & 1;
                if (bi != bj)
                {
                    int s2 = s ^ (1 << i) ^ (1 << j);   // hop across the bond
                    H[idx[s2], a] += 1.0;               // symmetric: the reverse hop fills the transpose entry
                }
            }
        }
        return H;
    }

    /// <summary>The darkest H-commutant coherence in the (p,q) excitation sector: the smallest nonzero
    /// eigenvalue of N_XY restricted to the ad_H kernel. This is the high-Q g2 contribution of that sector.
    /// Returns null if the sector is empty or carries no non-steady commutant mode.</summary>
    public static double? CommutantDarkest(string topo, int n, int p, int q)
    {
        var A = SectorStates(n, p);
        var B = SectorStates(n, q);
        int na = A.Length, nb = B.Length;
        if (na == 0 || nb == 0) return null;
        var Hp = SectorH(topo, n, A);
        var Hq = SectorH(topo, n, B);
        int dim = na * nb;

        // ad_H = kron(Hp, I_nb) - kron(I_na, Hq)   (Hq real symmetric, so Hq^T = Hq); index a*nb + b
        var adH = Matrix<double>.Build.Dense(dim, dim);
        for (int a = 0; a < na; a++)
            for (int c = 0; c < na; c++)
            {
                double h = Hp[a, c];
                if (h != 0.0) for (int b = 0; b < nb; b++) adH[a * nb + b, c * nb + b] += h;
            }
        for (int a = 0; a < na; a++)
            for (int b = 0; b < nb; b++)
                for (int d = 0; d < nb; d++)
                {
                    double h = Hq[b, d];
                    if (h != 0.0) adH[a * nb + b, a * nb + d] -= h;
                }

        var evd = adH.Evd(Symmetricity.Symmetric);
        var omega = evd.EigenValues;
        var V = evd.EigenVectors;

        // N_XY diagonal: entry at (a,b) = hamming(A[a], B[b])
        var diag = new double[dim];
        for (int a = 0; a < na; a++)
            for (int b = 0; b < nb; b++)
                diag[a * nb + b] = BitOperations.PopCount((uint)(A[a] ^ B[b]));

        // commutant = the Ω=0 (ad_H kernel) columns of V
        var cols = new List<int>();
        for (int k = 0; k < dim; k++) if (Math.Abs(omega[k].Real) < NullTol) cols.Add(k);
        if (cols.Count == 0) return null;

        // project: Ntil = U^T diag(N_XY) U, U = the commutant columns
        int mc = cols.Count;
        var Ntil = Matrix<double>.Build.Dense(mc, mc);
        for (int x = 0; x < mc; x++)
            for (int y = x; y < mc; y++)
            {
                double acc = 0.0;
                int cx = cols[x], cy = cols[y];
                for (int r = 0; r < dim; r++) acc += V[r, cx] * diag[r] * V[r, cy];
                Ntil[x, y] = acc;
                Ntil[y, x] = acc;
            }

        var w = Ntil.Evd(Symmetricity.Symmetric).EigenValues;
        double min = double.PositiveInfinity;
        foreach (var e in w) { double v = e.Real; if (v > NzTol && v < min) min = v; }
        return double.IsInfinity(min) ? (double?)null : min;
    }

    // ---- IInspectable ----

    public string DisplayName => "StructuralCeilingWitness (g2 closed forms from the commutant rep structure)";

    public string Summary =>
        "the structural ceiling computed live (typed home: StructuralCeilingClaim): g2 = the darkest [H,A]=0 " +
        "coherence in the largest degenerate single-particle level. Complete 4/N, star 4/(N−1), the N=4 (2,2) " +
        "outlier 2−2/√3; recomputed from the rep structure, not the 4^N Liouvillian. " +
        "Sector overview: inspect --root blockspectrum (this zooms the (1,1) commutant sector).";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return LadderNode("complete", "the complete graph K_N: g2 = 4/N (the (1,1) S_N standard rep)",
                new[] { 4, 5, 6, 7 }, n => 4.0 / n, "4/N");
            yield return LadderNode("star", "the star K_{1,N−1}: g2 = 4/(N−1) (the (1,1) leaf manifold)",
                new[] { 5, 6, 7, 8 }, n => 4.0 / (n - 1), "4/(N−1)");
            yield return OutlierNode();
        }
    }

    private static InspectableNode LadderNode(string topo, string title, int[] ns, Func<int, double> closed, string label)
    {
        var rows = new List<IInspectable>();
        foreach (var n in ns)
        {
            double? g = CommutantDarkest(topo, n, 1, 1);
            double cf = closed(n);
            string match = g.HasValue && Math.Abs(g.Value - cf) < 1e-9 ? "match" : "MISMATCH";
            string state = g.HasValue && g.Value < 1.0 - 1e-9 ? "ceiling" : "= 1 (band edge protects)";
            rows.Add(new InspectableNode($"{topo} N={n}",
                summary: $"g2 = {(g?.ToString("0.#########", Inv) ?? "—")} ({label} = {cf.ToString("0.#########", Inv)}, " +
                         $"{match}); {state}"));
        }
        return new InspectableNode(title,
            summary: $"the (1,1) commutant darkest coherence = {label}, the high-Q ceiling (recomputed live from the rep structure).",
            children: rows);
    }

    private static InspectableNode OutlierNode()
    {
        double k4Ladder = CommutantDarkest("complete", 4, 1, 1) ?? double.NaN;   // expect 1.0 (= 4/4)
        double k4Ceiling = CommutantDarkest("complete", 4, 2, 2) ?? double.NaN;  // expect 2 − 2/√3
        double ring4Ceiling = CommutantDarkest("ring", 4, 2, 2) ?? double.NaN;   // expect 1.0 (co-occupied)
        double k4Closed = 2.0 - 2.0 / Math.Sqrt(3.0);
        var rows = new List<IInspectable>
        {
            new InspectableNode("K_4 (1,1) ladder",
                summary: $"g2 = {k4Ladder.ToString("0.#########", Inv)} = 4/4 = 1.0 — the ladder reaches the band edge, no ceiling here"),
            new InspectableNode("K_4 (2,2) half-filling",
                summary: $"g2 = {k4Ceiling.ToString("0.#########", Inv)} (2 − 2/√3 = {k4Closed.ToString("0.#########", Inv)}, " +
                         $"{(Math.Abs(k4Ceiling - k4Closed) < 1e-7 ? "match" : "MISMATCH")}) — below the floor: the K_4 ceiling"),
            new InspectableNode("ring-4 (2,2) half-filling",
                summary: $"g2 = {ring4Ceiling.ToString("0.#########", Inv)} = 1.0 — co-occupies the floor (the known ring-4 (2,2) mode, Im = 2√2·J)"),
        };
        return new InspectableNode("the N=4 outlier (the (2,2) half-filling sector, shared by K_4 and ring-4)",
            summary: "the 4/N ladder hits 1.0 at N=4, so the ceiling moves to the (2,2) two-excitation sector: " +
                     "K_4 = 2 − 2/√3 dips below the floor, ring-4 = 1 co-occupies it. One mechanism, two topologies.",
            children: rows);
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
