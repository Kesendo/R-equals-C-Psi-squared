using System;
using System.Collections.Generic;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>The shape of the windowed F87 hardness obstruction, by (k, N). The §7.5/§7.6 derivation
/// proved soft ⟺ bipartite for any k (so the edge-mask set S having no odd 𝔽₂-relation IS the
/// soft class); this scan settles the size of that obstruction in pure GF(2) bit arithmetic (no
/// Hamiltonian, no Liouvillian). The answer (PROOF_F103 §7.7) is the law
/// max minimal-odd-cycle = min(2W − 1, 2k − 3), W = N − k + 1 windows: the K3 triangle at k = 3 is
/// just the smallest body bound (2k − 3 = 3), not a universal shape.
///
/// <para>For a k-body diagonal-cell (Klein (0,1)) Mixed term, the X/Y positions form a k-bit
/// window-mask; placed on the N-site chain by the sliding-window builder it contributes
/// {mask &lt;&lt; w : 0 ≤ w ≤ N−k}. The pair's edge set S is the union of both terms' window-masks,
/// and bipartite ⟺ S carries no odd subset XOR-ing to 0. The minimal such odd subset is the
/// obstruction; this scan reports its size distribution over the hard pairs. Reading a mask as a
/// GF(2)[x] polynomial (bit j = x^j) and a window shift as ×x^w turns the obstruction into a gcd
/// question (see <see cref="GcdFormulaSize"/>, <see cref="ValuationAtOnePlusX"/>); pure bit ops then
/// scale far past the Liouvillian route. The k=3 result reproduces the F103 §6 anchor (42 pairs,
/// 16 hard, 26 soft, every obstruction a triangle).</para></summary>
public static class WindowedObstructionScan
{
    /// <summary>A diagonal-cell (0,1) Mixed term reduced to what the GF(2) scan needs: its X/Y
    /// window-mask (k bits) and its #Y parity (pairs must be y_par-homogeneous).</summary>
    public readonly record struct CellTerm(ulong WindowMask, int YParity);

    /// <summary>Enumerate the k-body diagonal-cell (0,1) Mixed terms for Z-dephasing: #(X+Y) even
    /// and ≥ 2 (Mixed, nonzero window-mask), #(Y+Z) odd. Letters encoded base-4 (I=0, X=1, Y=2, Z=3).</summary>
    public static List<CellTerm> CellTerms(int k)
    {
        if (k < 1 || k > 30) throw new ArgumentOutOfRangeException(nameof(k));
        var terms = new List<CellTerm>();
        long total = 1L << (2 * k); // 4^k
        for (long code = 0; code < total; code++)
        {
            int na = 0, nb = 0, ny = 0;
            ulong mask = 0;
            for (int i = 0; i < k; i++)
            {
                int letter = (int)((code >> (2 * i)) & 3); // I=0, X=1, Y=2, Z=3
                bool bitA = letter == 1 || letter == 2;     // X or Y
                bool bitB = letter == 2 || letter == 3;     // Y or Z
                if (bitA) { na++; mask |= 1UL << i; }
                if (bitB) nb++;
                if (letter == 2) ny++;
            }
            if (na % 2 == 0 && nb % 2 == 1 && na >= 2)       // Klein (0,1), Mixed
                terms.Add(new CellTerm(mask, ny & 1));
        }
        return terms;
    }

    /// <summary>Distinct chain edge-masks S of a pair: each term's window-mask shifted across all
    /// windows [0, N−k].</summary>
    private static List<ulong> EdgeMasks(CellTerm t1, CellTerm t2, int k, int n, HashSet<ulong> buf)
    {
        buf.Clear();
        for (int w = 0; w <= n - k; w++) { buf.Add(t1.WindowMask << w); buf.Add(t2.WindowMask << w); }
        return new List<ulong>(buf);
    }

    /// <summary>Size of the smallest odd-cardinality subset of S that XORs to 0 (the minimal odd
    /// 𝔽₂-relation = the shortest odd cycle in the Cayley graph of S). 0 if none (bipartite / soft).</summary>
    public static int MinOddCycle(IReadOnlyList<ulong> s)
    {
        int n = s.Count;
        for (int size = 3; size <= n; size += 2)
        {
            var idx = new int[size];
            for (int i = 0; i < size; i++) idx[i] = i;
            while (true)
            {
                ulong x = 0;
                for (int i = 0; i < size; i++) x ^= s[idx[i]];
                if (x == 0) return size;
                int p = size - 1;
                while (p >= 0 && idx[p] == n - size + p) p--;
                if (p < 0) break;
                idx[p]++;
                for (int q = p + 1; q < size; q++) idx[q] = idx[q - 1] + 1;
            }
        }
        return 0;
    }

    /// <summary>Result of a (k, N) scan: pair counts and the size distribution of the minimal odd
    /// obstruction over the hard (non-bipartite) pairs.</summary>
    public readonly record struct ScanResult(
        int K, int N, int Pairs, int Hard, int Soft, IReadOnlyDictionary<int, int> MinOddCycleSizes)
    {
        /// <summary>True iff every hard pair's obstruction is a triangle (size 3) and there is at
        /// least one hard pair.</summary>
        public bool ObstructionIsAlwaysTriangle =>
            Hard > 0 && MinOddCycleSizes.Count == 1 && MinOddCycleSizes.ContainsKey(3);
    }

    // ---- GF(2)[x] view: a k-bit X/Y mask is a polynomial (bit j = coeff of x^j); a window shift is
    // multiplication by x^w. An odd relation is q_A p1 = q_B p2, and via g = gcd(p1,p2) the minimal
    // candidate is popcount(p1/g)+popcount(p2/g). This underlies the size law max = min(2W-1, 2k-3).

    /// <summary>GF(2)[x] degree (highest set bit) of a polynomial; -1 for 0.</summary>
    public static int PolyDegree(ulong p) => p == 0 ? -1 : 63 - System.Numerics.BitOperations.LeadingZeroCount(p);

    /// <summary>GF(2)[x] quotient a / b (b != 0), discarding the remainder.</summary>
    public static ulong PolyDivQuotient(ulong a, ulong b)
    {
        if (b == 0) throw new DivideByZeroException();
        ulong q = 0;
        int db = PolyDegree(b);
        while (a != 0 && PolyDegree(a) >= db)
        {
            int shift = PolyDegree(a) - db;
            q ^= 1UL << shift;
            a ^= b << shift;
        }
        return q;
    }

    /// <summary>GF(2)[x] gcd via the Euclidean algorithm.</summary>
    public static ulong PolyGcd(ulong a, ulong b)
    {
        while (b != 0)
        {
            int db = PolyDegree(b);
            ulong r = a;
            while (r != 0 && PolyDegree(r) >= db) r ^= b << (PolyDegree(r) - db);
            a = b; b = r;
        }
        return a;
    }

    /// <summary>The (1+x)-adic valuation of a polynomial: how many times (1+x) = 0b11 divides it.
    /// A diagonal-cell mask has even popcount, so p(1) = 0 and the valuation is ≥ 1.</summary>
    public static int ValuationAtOnePlusX(ulong p)
    {
        int v = 0;
        while (p != 0 && (System.Numerics.BitOperations.PopCount(p) & 1) == 0) // p(1) == 0
        {
            p = PolyDivQuotient(p, 0b11);
            v++;
        }
        return v;
    }

    /// <summary>The gcd-formula candidate odd-relation size popcount(p1/g)+popcount(p2/g), an upper
    /// bound on the minimal odd cycle (the s=1 relation; a shorter one can exist via cancellation).</summary>
    public static int GcdFormulaSize(ulong p1, ulong p2)
    {
        ulong g = PolyGcd(p1, p2);
        return System.Numerics.BitOperations.PopCount(PolyDivQuotient(p1, g))
             + System.Numerics.BitOperations.PopCount(PolyDivQuotient(p2, g));
    }

    /// <summary>Scan all y_par-homogeneous diagonal-cell Mixed pairs at body count k on an N-site
    /// chain (k &lt; N is the windowed regime), classifying each by its minimal odd obstruction.</summary>
    public static ScanResult Scan(int k, int n)
    {
        if (n <= k) throw new ArgumentOutOfRangeException(nameof(n), "windowed regime needs n > k");
        var terms = CellTerms(k);
        int pairs = 0, hard = 0, soft = 0;
        var dist = new Dictionary<int, int>();
        var buf = new HashSet<ulong>();
        for (int a = 0; a < terms.Count; a++)
            for (int b = a; b < terms.Count; b++)
            {
                if (terms[a].YParity != terms[b].YParity) continue; // y_par-homogeneous
                pairs++;
                int cyc = MinOddCycle(EdgeMasks(terms[a], terms[b], k, n, buf));
                if (cyc == 0) soft++;
                else { hard++; dist[cyc] = dist.GetValueOrDefault(cyc) + 1; }
            }
        return new ScanResult(k, n, pairs, hard, soft, dist);
    }
}
