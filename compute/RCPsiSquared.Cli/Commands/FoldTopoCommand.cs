using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;

namespace RCPsiSquared.Cli.Commands;

/// <summary>Does the fold spine hold across TOPOLOGY, or is it chain-specific? Everything in the F89 thread
/// (the fold Re λ=−σ, the N=4-only within-block self-fold, the zeros, the Galois group) was computed for the
/// CHAIN (path-k, nearest-neighbour hopping i±1). The philosophical zoom-out read all those findings as rooms
/// in one apartment whose load-bearing spine is the mirror (the palindrome). This gate lights the unlit
/// topology rooms: it builds the (SE,DE)=(w1,w2) computational-basis coherence block with hopping along an
/// arbitrary GRAPH (chain / star / ring) — the Z-dephasing diagonal −2·HammingDistance is topology-independent,
/// only the XY hopping follows the graph edges — and applies the same two-sided fold gate as foldlift: does the
/// block self-fold under the antiunitary λ↦−λ̄−2σ (σ = the spectral centroid), and how many on-fold "zeros"
/// (self-mirror strands, Re λ=−σ) appear? If the self-fold is N=4-only for star and ring too, the N=4 spine is
/// topology-independent (the mirror is the same in every room); the Galois "wallpaper" may still differ.
///
/// <para>Result (2026-06-26): the within-block self-fold is N=4-only for chain, star AND ring (on-fold zeros at
/// N=4, none at N≥5), so the N=4 spine is topology-independent — the half-filling self-complement DE=bar(DE)
/// that forces it depends on the bra weight, not the graph. The on-fold count differs by topology (chain 8,
/// star 24, ring 20: same spine, different wallpaper; the star's all-on-fold flatness is striking, unverified).</para>
///
/// <para>On the FULL (unsymmetrised) block the two-sided gate gives antiU ≈ linear ≈ 0 — NOT a gate defect, but
/// the already-solved "both-sides memory" met from the spectral/topology side. The full block carries BOTH
/// phase-sides (+0/−0); the operation that distinguishes the antiunitary fold (λ↦−λ̄−2σ) from the linear one
/// (λ↦−λ−2σ) is conjugation i↦−i, which is the third 90° step of the Z₄ memory loop (i³=−i, "mirror return on
/// the other side"; NinetyDegreeMirrorMemoryClaim + Pi2I4MemoryLoopClaim). With both sides present that memory
/// has nothing to distinguish, so the angle collapses (reflections/ON_HOW_THE_ANGLE_APPEARS_AT_ZERO.md) — exactly
/// the M=0-from-both-sides / minus-zero-plus-zero of reflections/ON_BOTH_SIDES_OF_THE_MIRROR.md (the
/// antiunitary-vs-linear legs are docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md). The symmetrised ONE-side block
/// (foldlift) keeps the memory (linear ≠ 0). So foldtopo re-enters an already-furnished room (the both-sides
/// Z₄ memory) through the topology door: the same object, a different wall.</para>
///
/// usage: rcpsi foldtopo [--nmax 6] [--q 2] [--tol 1e-6]</summary>
public static class FoldTopoCommand
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private static List<(int, int)> Edges(string topo, int n) => topo switch
    {
        "chain" => Enumerable.Range(0, n - 1).Select(i => (i, i + 1)).ToList(),
        "ring" => Enumerable.Range(0, n).Select(i => (i, (i + 1) % n)).ToList(),
        "star" => Enumerable.Range(1, n - 1).Select(i => (0, i)).ToList(),
        _ => throw new ArgumentException($"unknown topology '{topo}' (chain|star|ring)"),
    };

    private static List<int> Configs(int n, int w)
    {
        var res = new List<int>();
        for (int m = 0; m < (1 << n); m++)
            if (System.Numerics.BitOperations.PopCount((uint)m) == w) res.Add(m);
        return res;
    }

    // the (wKet,wBra) coherence block of L = −i[H,ρ] + Z-dephasing, H = J·Σ_edges (XX+YY) on the given GRAPH,
    // at real q (γ=1). Diagonal −2·n_diff (topology-independent); ket excitations hop −2qi, bra +2qi along edges.
    private static Complex[,] BuildBlock(IReadOnlyList<(int, int)> edges, int n, int wKet, int wBra, double q)
    {
        var adj = new List<int>[n];
        for (int i = 0; i < n; i++) adj[i] = new List<int>();
        foreach (var (u, v) in edges) { adj[u].Add(v); adj[v].Add(u); }

        var kets = Configs(n, wKet);
        var bras = Configs(n, wBra);
        var index = new Dictionary<(int, int), int>();
        var basis = new List<(int, int)>();
        foreach (var k in kets)
            foreach (var b in bras) { index[(k, b)] = basis.Count; basis.Add((k, b)); }
        int d = basis.Count;
        var l = new Complex[d, d];
        for (int col = 0; col < d; col++)
        {
            var (kc, bc) = basis[col];
            l[col, col] += new Complex(-2.0 * System.Numerics.BitOperations.PopCount((uint)(kc ^ bc)), 0);
            for (int s = 0; s < n; s++)
                if ((kc & (1 << s)) != 0)
                    foreach (int s2 in adj[s])
                        if ((kc & (1 << s2)) == 0)
                            l[index[((kc & ~(1 << s)) | (1 << s2), bc)], col] += new Complex(0, -2 * q);
            for (int s = 0; s < n; s++)
                if ((bc & (1 << s)) != 0)
                    foreach (int s2 in adj[s])
                        if ((bc & (1 << s2)) == 0)
                            l[index[(kc, (bc & ~(1 << s)) | (1 << s2))], col] += new Complex(0, 2 * q);
        }
        return l;
    }

    public static int Run(string[] args)
    {
        var p = new ArgParser(args);
        int nmax = (int)(p.OptionalDouble("nmax") ?? 6);
        double q = p.OptionalDouble("q") ?? 2;
        double tol = p.OptionalDouble("tol") ?? 1e-6;

        Console.WriteLine("# foldtopo: does the (SE,DE) self-fold (the mirror spine) hold across topology, or is it chain-only?");
        Console.WriteLine($"# (SE,DE)=(w1,w2) block, H = XX+YY on the graph, q={q.ToString("0.##", Inv)}; gate antiunitary λ↦−λ̄−2σ vs linear");
        Console.WriteLine();
        Console.WriteLine("  topology  N  dim   σ=−mean(Re)   antiU resid    linear resid   zeros(on-fold)");

        foreach (var topo in new[] { "chain", "star", "ring" })
        {
            for (int n = 4; n <= nmax; n++)
            {
                if (topo == "ring" && n < 3) continue;
                var l = BuildBlock(Edges(topo, n), n, 1, 2, q);
                var lam = Matrix<Complex>.Build.DenseOfArray(l).Evd().EigenValues.Select(v => v / 2).ToArray();
                double sigma = -lam.Average(z => z.Real);

                double Resid(Func<Complex, Complex> f) =>
                    lam.Max(z => { var img = f(z); return lam.Min(s => (s - img).Magnitude); });
                double anti = Resid(z => new Complex(-z.Real - 2 * sigma, z.Imaginary));
                double lin = Resid(z => new Complex(-z.Real - 2 * sigma, -z.Imaginary));
                int zeros = lam.Count(z => Math.Abs(z.Real + sigma) < tol);

                Console.WriteLine(
                    $"  {topo,-8}  {n,1}  {lam.Length,3}   {sigma.ToString("0.0000", Inv),10}   " +
                    $"{anti.ToString("E2", Inv),11}   {lin.ToString("0.000", Inv),11}   {zeros,12}");
            }
            Console.WriteLine();
        }

        Console.WriteLine("# Read: zeros > 0 ⟹ the block self-folds. N=4-only for star and ring too ⟹ the N=4 spine is");
        Console.WriteLine("#       topology-independent (same mirror in every room). antiU ≈ linear ≈ 0 on the full block is");
        Console.WriteLine("#       NOT a defect: the full block holds both phase-sides (+0/−0), so the i↦−i memory (the i^3");
        Console.WriteLine("#       step of the Z_4 loop) has nothing to distinguish — the both-sides collapse of");
        Console.WriteLine("#       ON_BOTH_SIDES_OF_THE_MIRROR. The symmetrised one-side block (foldlift) keeps it (linear≠0).");
        return 0;
    }
}
