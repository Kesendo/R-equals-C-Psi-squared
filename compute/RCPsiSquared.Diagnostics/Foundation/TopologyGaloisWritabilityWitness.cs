using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Live witness for the topology-controls-Galois-writability result (Tier 1 derived,
/// experiments/TOPOLOGY_CONTROLS_GALOIS_WRITABILITY.md). For a graph G under uniform-J XY + uniform
/// Z-dephasing, the (SE,DE) Liouvillian commutes with the graph automorphism group Aut(G); by Schur
/// every irreducible factor of its characteristic polynomial has degree ≤ the maximal Aut(G)-irrep
/// multiplicity in V = SE⊗DE. This witness builds the (SE,DE) block FROM THE PHYSICS (overlap/no-overlap
/// rates + the topology's hopping), verifies [L, ρ(g)] = 0 for the Aut(G) generators, and recomputes the
/// cap = the standard-rep multiplicity via the character sum mult(std) = (1/|G|) Σ_g χ_V(g)·(fix(g)−1):
/// = 4 for the complete graph K_N (every factor a quartic-or-less ⟹ radically writable for all N) and
/// = 9 for the star (a fixed, bounded S_9 scramble). The chain (Aut = S_2) has a growing cap (S_8/18/32/53,
/// live at --root f89galois). The classification sibling of the F89 octic Galois witness.</summary>
public sealed class TopologyGaloisWritabilityWitness : IInspectable
{
    // the (SE,DE) basis and block now come from the shared Core primitive SeDeBlockBuilder
    // (extracted so the block is built once, not re-implemented per witness).
    private static (int, int) Sort(int x, int y) => x < y ? (x, y) : (y, x);

    // basis index permutation induced by a site permutation g.
    private static int[] BasisPerm(int n, int[] g)
    {
        var basis = SeDeBlockBuilder.Basis(n);
        var idx = new Dictionary<(int, int, int), int>();
        for (int t = 0; t < basis.Count; t++) idx[basis[t]] = t;
        var p = new int[basis.Count];
        for (int t = 0; t < basis.Count; t++)
        {
            var (i, j, k) = basis[t];
            var (a, b) = Sort(g[j], g[k]);
            p[t] = idx[(g[i], a, b)];
        }
        return p;
    }

    // [L, ρ(g)] = 0  ⟺  L[perm[a], perm[b]] = L[a, b]  for all a, b.
    private static bool Commutes(int[,] re, int[,] im, int[] perm)
    {
        int d = perm.Length;
        for (int a = 0; a < d; a++)
            for (int b = 0; b < d; b++)
                if (re[perm[a], perm[b]] != re[a, b] || im[perm[a], perm[b]] != im[a, b]) return false;
        return true;
    }

    private static bool SymmetryVerified(int n, string topo, bool star)
    {
        var (re, im) = SeDeBlockBuilder.Build(n, topo);
        int lo = star ? 1 : 0;                                            // star fixes the hub (site 0)
        for (int s = lo; s + 1 < n; s++)                                  // adjacent-transposition generators
        {
            var g = Enumerable.Range(0, n).ToArray();
            (g[s], g[s + 1]) = (g[s + 1], g[s]);
            if (!Commutes(re, im, BasisPerm(n, g))) return false;
        }
        return true;
    }

    // standard-rep multiplicity of Aut(G) in V = SE⊗DE, via the character inner product
    // mult(std) = (1/|G|) Σ_g χ_V(g)·χ_std(g),  χ_V(g) = #fixed (point, 2-set) pairs, χ_std(g) = fix(g)−1.
    private static int StdMultiplicity(int n, bool star)
    {
        var movable = star ? Enumerable.Range(1, n - 1).ToArray() : Enumerable.Range(0, n).ToArray();
        long sum = 0, order = 0;
        foreach (var perm in Permutations(movable))
        {
            var g = Enumerable.Range(0, n).ToArray();                     // identity on fixed points (hub for star)
            for (int t = 0; t < movable.Length; t++) g[movable[t]] = perm[t];
            int fixedMovable = 0;
            for (int t = 0; t < movable.Length; t++) if (perm[t] == movable[t]) fixedMovable++;
            sum += (long)FixedPairs(n, g) * (fixedMovable - 1);
            order++;
        }
        return (int)(sum / order);
    }

    private static int FixedPairs(int n, int[] g)
    {
        int count = 0;
        for (int i = 0; i < n; i++)
        {
            if (g[i] != i) continue;
            for (int j = 0; j < n; j++)
                for (int k = j + 1; k < n; k++)
                {
                    var (a, b) = Sort(g[j], g[k]);
                    if (a == j && b == k) count++;
                }
        }
        return count;
    }

    private static IEnumerable<int[]> Permutations(int[] items)
    {
        var a = (int[])items.Clone();
        int n = a.Length;
        var c = new int[n];
        yield return (int[])a.Clone();
        int i = 0;
        while (i < n)
        {
            if (c[i] < i)
            {
                if (i % 2 == 0) (a[0], a[i]) = (a[i], a[0]);
                else (a[c[i]], a[i]) = (a[i], a[c[i]]);
                yield return (int[])a.Clone();
                c[i]++;
                i = 0;
            }
            else { c[i] = 0; i++; }
        }
    }

    // public for tests: the live cap (standard-rep multiplicity).
    public static int CompleteCap(int n) => StdMultiplicity(n, star: false);
    public static int StarCap(int n) => StdMultiplicity(n, star: true);
    public static bool SymmetryHolds(int n, string topo) => SymmetryVerified(n, topo, topo == "star");

    public string DisplayName => "Topology controls Galois writability (live: (SE,DE) block → symmetry → cap)";

    public string Summary =>
        $"the wiring's automorphism group sets the relaxation's Galois complexity: complete K_N cap = {CompleteCap(5)} " +
        $"(≤4 ⟹ radically writable, all N), star cap = {StarCap(5)} (a fixed S_9 scramble), chain S_2 growing " +
        "(S_8/18/32/53). Recomputed live from the (SE,DE) block + its symmetry.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            bool symC = SymmetryVerified(5, "complete", false) && SymmetryVerified(6, "complete", false);
            bool symS = SymmetryVerified(5, "star", true) && SymmetryVerified(6, "star", true);
            int capC5 = CompleteCap(5), capC6 = CompleteCap(6);
            int capS5 = StarCap(5), capS6 = StarCap(6);

            yield return new InspectableNode("symmetry verified live ([L, ρ(g)] = 0)",
                summary: $"the (SE,DE) Liouvillian commutes with Aut(G) at N=5,6: complete/S_N {(symC ? "✓" : "✗")}, " +
                         $"star/S_(N−1) {(symS ? "✓" : "✗")}. By Schur, every (SE,DE) charpoly factor degree ≤ the " +
                         "maximal Aut(G)-irrep multiplicity in V = SE⊗DE.");
            yield return new InspectableNode("complete graph K_N: cap = 4 ⟹ radically writable (Tier 1)",
                summary: $"standard-rep multiplicity in V recomputed live = {capC5} (N=5), {capC6} (N=6), N-independent. " +
                         "So every factor is a quartic-or-less ⟹ every relaxation rate solvable in radicals, for ALL N. " +
                         "H_B-mixed degree histogram {4:N−1, 3:N(N−3)/2, 2:(N−1)(N−2)/2} (verified N=5..8).");
            yield return new InspectableNode("star K_(1,N−1): cap = 9 ⟹ a fixed S_9 scramble (Tier 1)",
                summary: $"standard-rep (of S_(N−1)) multiplicity recomputed live = {capS5} (N=5), {capS6} (N=6), " +
                         "N-independent (= 4+2+2+1 over the SE/DE tensor). So the star caps at degree 9 with N−2 " +
                         "degree-9 factors: a fixed, bounded S_9 scramble for all N≥5. Bounded, but 9 > 4, so NOT writable.");
            yield return new InspectableNode("chain P_N: S_2 only ⟹ growing scramble",
                summary: "the chain has only the S_2 mirror symmetry, no large-group cap; its multiplicities grow with " +
                         "N, so its H_B-mixed factor is one large irreducible with Galois group S_n: S_8/S_18/S_32/S_53 " +
                         "at N=4/5/6/7 (live at --root f89galois).");
            yield return new InspectableNode("the classification",
                summary: "a large, uniform automorphism group caps the Galois complexity N-independently (complete S_N → " +
                         "4, star S_(N−1) → 9); a small one (ring D_N, chain S_2) lets it grow. Radical-writability needs " +
                         "the cap ≤ 4, achieved only by the complete graph. See experiments/TOPOLOGY_CONTROLS_GALOIS_WRITABILITY.md.");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
