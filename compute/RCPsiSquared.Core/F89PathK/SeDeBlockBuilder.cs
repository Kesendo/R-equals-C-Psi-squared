using System.Collections.Generic;

namespace RCPsiSquared.Core.F89PathK;

/// <summary>The ONE (SE,DE) Liouvillian sub-block builder, shared by every consumer (the topology
/// Galois-writability witness and the Galois-vs-spectral-chaos witness) so the block is built once,
/// never re-implemented. Basis = (marked point i, 2-subset {j,k}) over N sites, dim N·C(N,2). At
/// J=γ=1 the block is L = re + i·im (integer matrices): the diagonal AT rate is −2 on overlap
/// (i ∈ {j,k}, dephasing 2γ) / −6 otherwise (no-overlap, 6γ), and im is the XY hopping (SE ket hop
/// −2, DE bra hops +2) governed by the graph adjacency. Scale to a physical knob q=J/γ by L(q) =
/// re + i·q·im. Extracted verbatim from TopologyGaloisWritabilityWitness.Block, with the ring (C_N)
/// added to the chain/star/complete the witness carried.</summary>
public static class SeDeBlockBuilder
{
    /// <summary>The (marked point i, 2-subset {j,k}, j&lt;k) basis over n sites; i outer, pair inner.</summary>
    public static IReadOnlyList<(int i, int j, int k)> Basis(int n)
    {
        var b = new List<(int, int, int)>();
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                for (int k = j + 1; k < n; k++)
                    b.Add((i, j, k));
        return b;
    }

    private static (int, int) SortPair(int x, int y) => x < y ? (x, y) : (y, x);

    /// <summary>Undirected adjacency of the wiring graph on n sites: complete K_N, star K_{1,N−1}
    /// (hub = site 0), ring C_N, or chain P_N (the default).</summary>
    public static bool[,] Adjacency(int n, string topo)
    {
        var a = new bool[n, n];
        if (topo == "complete")
        {
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    a[i, j] = i != j;
        }
        else if (topo == "star")
        {
            for (int l = 1; l < n; l++) { a[0, l] = true; a[l, 0] = true; }
        }
        else if (topo == "ring")
        {
            for (int i = 0; i < n; i++) { int j = (i + 1) % n; a[i, j] = true; a[j, i] = true; }
        }
        else   // chain
        {
            for (int i = 0; i + 1 < n; i++) { a[i, i + 1] = true; a[i + 1, i] = true; }
        }
        return a;
    }

    /// <summary>L = Re (diagonal AT rates −2/−6) + i·Im (the ±2 XY hopping), integers at J=γ=1.</summary>
    public static (int[,] re, int[,] im) Build(int n, bool[,] adj)
    {
        var basis = Basis(n);
        int d = basis.Count;
        var idx = new Dictionary<(int, int, int), int>();
        for (int t = 0; t < d; t++) idx[basis[t]] = t;
        var re = new int[d, d];
        var im = new int[d, d];
        for (int col = 0; col < d; col++)
        {
            var (i, j, k) = basis[col];
            re[col, col] = (i == j || i == k) ? -2 : -6;                  // overlap 2γ / no-overlap 6γ
            for (int i2 = 0; i2 < n; i2++) if (adj[i, i2]) im[idx[(i2, j, k)], col] += -2;   // SE hop (ket)
            for (int nj = 0; nj < n; nj++) if (adj[j, nj] && nj != k) { var (a, b) = SortPair(nj, k); im[idx[(i, a, b)], col] += 2; }
            for (int nk = 0; nk < n; nk++) if (adj[k, nk] && nk != j) { var (a, b) = SortPair(j, nk); im[idx[(i, a, b)], col] += 2; }
        }
        return (re, im);
    }

    /// <summary>L = Re + i·Im for a named topology (complete / star / ring / chain).</summary>
    public static (int[,] re, int[,] im) Build(int n, string topo) => Build(n, Adjacency(n, topo));
}
