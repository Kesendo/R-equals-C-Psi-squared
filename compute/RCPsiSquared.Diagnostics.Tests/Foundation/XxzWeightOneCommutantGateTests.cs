using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The weight-1 commutant of the open XXZ chain H = Σ_b (X_bX_{b+1} + Y_bY_{b+1} + Δ·Z_bZ_{b+1}),
/// block by block. Every weight-1 Pauli string (one X or Y letter, I or Z elsewhere) lives in a block (p, p+1) or
/// its transpose, on the cells |a⟩⟨b| with b = a plus one occupied site, and by the theorem of
/// docs/proofs/PROOF_WEIGHT1_DEGENERACY.md ("The count at exceptional couplings") the real modes at Re λ = −2γ
/// under uniform Z dephasing are, off a finite set of γ/J, exactly this commutant.
///
/// <para>The closed form. With Δ = cos θ, put on the cell that adds site l the coefficient T_|m|(Δ) = cos(m·θ), m the
/// number of occupied minus the number of empty sites left of l: v_p = Σ_l cos(θ·Σ_{k&lt;l}(2n̂_k − 1))·σ_l restricted
/// to the block, the real part of Σ_l q^(Σ_{k&lt;l}(2n̂_k − 1))·σ_l with q = e^{iθ}, a Jordan-Wigner string whose
/// sign has been turned into a phase (q = i at Δ = 0, q = 1 at Δ = 1). Its residual
/// [H, v_p] is (integer)·(1 − Δ²)·U_{|z_p|−1}(Δ) on the cells adding the LAST site and zero elsewhere, with
/// z_p = N − 1 − 2p the Z sum of the other N − 1 sites; so v_p commutes with H exactly when z_p = 0 (the middle
/// block at odd N), Δ = ±1, or, for |Δ| ≤ 1, sin(z_p·θ) = 0; for |Δ| > 1 the factor is sinh φ·sinh(|z_p|·φ),
/// |Δ| = cosh φ, and only the middle block qualifies.</para>
///
/// <para>Two exact routes, meeting. (1) The closed form's residual, in integers: Δ = a/c, H scaled by c, v scaled by
/// c^(N−1), every entry a BigInteger compared exactly against ±2 times its predicted factor, with nonzero blocks
/// asserted nonzero, and v_p itself nonzero (a lower bound, kernel ≥ 1 where the law holds). (2) dim ker(ad_H on the block's weight-1 cells) by a
/// rank modulo a prime, which can only overcount the rational kernel (an upper bound). The count reads exactly
/// 1 where the law holds and 0 elsewhere, so the two bounds close. Inputs: Δ ∈ {0, ±1/2, 3/10, 4/5, ±1, 2, 7/3},
/// including the root-of-unity points that open extra blocks (Δ = 0: every block at odd N; Δ = ±1/2: 3 | z_p),
/// irrational θ (Niven) and |Δ| > 1, which open only the middle block.</para></summary>
public class XxzWeightOneCommutantGateTests
{
    private readonly ITestOutputHelper _out;
    public XxzWeightOneCommutantGateTests(ITestOutputHelper o) => _out = o;

    private const long Prime = 1_000_000_007;

    private static int Pop(int x) => BitOperations.PopCount((uint)x);

    /// <summary>c·H|s⟩ for H = Σ_b (XX + YY) + (a/c)·ZZ on the open chain: hop amplitude 2c, diagonal a·Σ z_i z_j.</summary>
    private static List<(int State, BigInteger Coeff)> ScaledH(int n, int a, int c, int s)
    {
        var res = new List<(int, BigInteger)>();
        BigInteger diag = 0;
        for (int i = 0; i + 1 < n; i++)
        {
            int bi = (s >> i) & 1, bj = (s >> (i + 1)) & 1;
            diag += bi == bj ? a : -a;
            if (bi != bj) res.Add((s ^ (1 << i) ^ (1 << (i + 1)), 2 * c));
        }
        res.Add((s, diag));
        return res;
    }

    /// <summary>Numerators of T_k(a/c) and U_k(a/c) over c^k.</summary>
    private static BigInteger[] ChebNum(int a, int c, int kmax, bool second)
    {
        var t = new BigInteger[Math.Max(kmax + 1, 2)];
        t[0] = 1; t[1] = second ? 2 * a : a;
        for (int k = 1; k + 1 < t.Length; k++) t[k + 1] = 2 * a * t[k] - (BigInteger)c * c * t[k - 1];
        return t;
    }

    private static IEnumerable<(int Ket, int Bra, int Site)> Cells(int n, int p)
    {
        for (int ket = 0; ket < 1 << n; ket++)
        {
            if (Pop(ket) != p) continue;
            for (int l = 0; l < n; l++)
                if (((ket >> l) & 1) == 0) yield return (ket, ket | (1 << l), l);
        }
    }

    private static int LeftImbalance(int ket, int l)
    {
        int occ = Pop(ket & ((1 << l) - 1));
        return occ - (l - occ);
    }

    public static IEnumerable<object[]> Grid()
    {
        var deltas = new[] { (0, 1), (1, 2), (-1, 2), (3, 10), (4, 5), (1, 1), (-1, 1), (2, 1), (7, 3) };
        foreach (int n in new[] { 3, 4, 5, 6, 7, 8 })
            foreach (var (a, c) in deltas)
                yield return new object[] { n, a, c };
    }

    /// <summary>The law's factor for block p, scaled: (c² − a²)·Un_{|m|−1}·c^(N−1−|m|); zero iff the closed form commutes.</summary>
    private static BigInteger PredictedFactor(int n, int a, int c, int p)
    {
        int m = Math.Abs(n - 1 - 2 * p);
        if (m == 0) return 0;
        var u = ChebNum(a, c, n, second: true);
        return ((BigInteger)c * c - (BigInteger)a * a) * u[m - 1] * BigInteger.Pow(c, n - 1 - m);
    }

    [Theory]
    [MemberData(nameof(Grid))]
    public void ClosedForm_Residual_IsTheLawsFactor_OnTheLastSiteOnly(int n, int a, int c)
    {
        var tn = ChebNum(a, c, n, second: false);
        for (int p = 0; p < n; p++)
        {
            var residual = new Dictionary<(int, int), BigInteger>();
            void Add(int ket, int bra, BigInteger x)
            {
                residual.TryGetValue((ket, bra), out var y);
                residual[(ket, bra)] = y + x;
            }
            bool anyCoefficient = false;
            foreach (var (ket, bra, site) in Cells(n, p))
            {
                int m = Math.Abs(LeftImbalance(ket, site));
                BigInteger coeff = tn[m] * BigInteger.Pow(c, n - 1 - m);
                anyCoefficient |= !coeff.IsZero;
                foreach (var (k2, h) in ScaledH(n, a, c, ket)) Add(k2, bra, h * coeff);
                foreach (var (b2, h) in ScaledH(n, a, c, bra)) Add(ket, b2, -h * coeff);
            }
            Assert.True(anyCoefficient, $"N={n} Δ={a}/{c} p={p}: the closed form is the zero operator");
            var nonzero = residual.Where(kv => !kv.Value.IsZero).ToList();
            BigInteger factor = PredictedFactor(n, a, c, p);
            if (factor.IsZero)
            {
                Assert.True(nonzero.Count == 0,
                    $"N={n} Δ={a}/{c} p={p}: the law says the closed form commutes, residual has {nonzero.Count} nonzero entries");
            }
            else
            {
                Assert.True(nonzero.Count > 0, $"N={n} Δ={a}/{c} p={p}: predicted nonzero residual, got zero");
                foreach (var ((ket, bra), x) in nonzero)
                {
                    Assert.True((ket ^ bra) == 1 << (n - 1),
                        $"N={n} Δ={a}/{c} p={p}: residual on a cell not adding the last site");
                    Assert.True(BigInteger.Abs(x) == 2 * BigInteger.Abs(factor),
                        $"N={n} Δ={a}/{c} p={p}: residual {x} is not ±2 times the law's factor {factor}");
                }
            }
        }
    }

    [Theory]
    [MemberData(nameof(Grid))]
    public void WeightOneCommutant_IsOneWhereTheLawHolds_AndZeroElsewhere(int n, int a, int c)
    {
        var counts = new int[n];
        for (int p = 0; p < n; p++)
        {
            var cells = Cells(n, p).ToList();
            var rowIndex = new Dictionary<(int, int), int>();
            var vectors = new List<Dictionary<int, long>>();
            foreach (var (ket, bra, _) in cells)
            {
                var v = new Dictionary<int, long>();
                void Add(int k, int b, BigInteger x)
                {
                    long xm = (long)(((x % Prime) + Prime) % Prime);
                    if (xm == 0) return;
                    if (!rowIndex.TryGetValue((k, b), out int r)) rowIndex[(k, b)] = r = rowIndex.Count;
                    v.TryGetValue(r, out long y);
                    long z = (y + xm) % Prime;
                    if (z == 0) v.Remove(r); else v[r] = z;
                }
                foreach (var (k2, h) in ScaledH(n, a, c, ket)) Add(k2, bra, h);
                foreach (var (b2, h) in ScaledH(n, a, c, bra)) Add(ket, b2, -h);
                vectors.Add(v);
            }
            counts[p] = cells.Count - RankModPrime(vectors);
            int expected = PredictedFactor(n, a, c, p).IsZero ? 1 : 0;
            Assert.True(counts[p] == expected,
                $"N={n} Δ={a}/{c} block ({p},{p + 1}): weight-1 commutant {counts[p]}, law predicts {expected}");
        }
        _out.WriteLine($"N={n} Δ={a}/{c}: per block (p,p+1) = [{string.Join(", ", counts)}]");
    }

    private static int RankModPrime(List<Dictionary<int, long>> vectors)
    {
        var pivots = new Dictionary<int, Dictionary<int, long>>();
        foreach (var src in vectors)
        {
            var r = new Dictionary<int, long>(src);
            while (r.Count > 0)
            {
                int col = r.Keys.Min();
                if (pivots.TryGetValue(col, out var pr))
                {
                    long f = r[col];
                    foreach (var (k, val) in pr)
                    {
                        r.TryGetValue(k, out long y);
                        long z = ((y - f * val) % Prime + Prime) % Prime;
                        if (z == 0) r.Remove(k); else r[k] = z;
                    }
                }
                else
                {
                    long inv = ModPow(r[col], Prime - 2);
                    pivots[col] = r.ToDictionary(kv => kv.Key, kv => kv.Value * inv % Prime);
                    break;
                }
            }
        }
        return pivots.Count;
    }

    private static long ModPow(long b, long e)
    {
        long result = 1; b %= Prime;
        while (e > 0) { if ((e & 1) == 1) result = result * b % Prime; b = b * b % Prime; e >>= 1; }
        return result;
    }
}
