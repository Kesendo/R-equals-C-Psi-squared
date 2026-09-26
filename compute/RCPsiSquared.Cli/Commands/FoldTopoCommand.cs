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
/// arbitrary GRAPH (chain / star / ring); the Z-dephasing diagonal −2·HammingDistance is topology-independent,
/// only the XY hopping follows the graph edges. It applies the same two-sided fold gate as foldlift: does the
/// block self-fold under the antiunitary λ↦−λ̄−2σ (σ = the spectral centroid), and how many on-fold "zeros"
/// (self-mirror strands, Re λ=−σ) appear? If the self-fold is N=4-only for star and ring too, the N=4 spine is
/// topology-independent (the mirror is the same in every room); the Galois "wallpaper" may still differ.
///
/// <para>Result (2026-06-26): the within-block self-fold is N=4-only for chain, star AND ring (on-fold zeros at
/// N=4, none at N≥5), so the N=4 spine is topology-independent: the half-filling self-complement DE=bar(DE)
/// that forces it depends on the bra weight, not the graph. The on-fold count differs by topology and by q (chain 8,
/// star 24, ring 20 at q = 2; 4, 8, 4 at q = 0.7): same spine, different wallpaper.</para>
///
/// <para>On the FULL (unsymmetrised) block the antiunitary fold (λ↦−λ̄−2σ) and the linear one (λ↦−λ−2σ)
/// differ only by complex conjugation, so on a block spectrum S that is closed under conjugation (S = S̄) the
/// two folds are one test and their residuals coincide. The run shows exactly that: antiU ≈ linear ≈ 0 at N=4
/// on all three graphs, and equal nonzero residuals at N=5 and 6 on the chain, the star and the even ring,
/// while the odd ring (N=5) is the one row where they differ (0.471 against 1.402 at q=2), its block spectrum
/// not being conjugation-closed. docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md §3.4 closes the spectrum of a
/// Pauli-real generator under conjugation; that covers the full L and does not reach inside this block: for the
/// real hopping H here, conjugation flips L_H and the transpose ρ ↦ ρᵀ flips it back while keeping the
/// dephasing diagonal, so the conjugate of the (w1, w2) block's spectrum is the (w2, w1) block's (checked
/// numerically on the N=4 chain and the N=5 ring, to 4e-14). Closure within the one block is extra, and the
/// odd ring shows it can fail. The symmetrised ONE-side block (<c>foldlift</c>) keeps a nonzero linear
/// residual.</para>
///
/// <para>A reading: conjugation i↦−i is the third 90° step of the Z₄ memory loop (i³ = −i, the mirror's return
/// on the other side; NinetyDegreeMirrorMemoryClaim + Pi2I4MemoryLoopClaim). Where the full block carries both
/// phase sides, that memory has nothing to distinguish and the angle collapses
/// (reflections/ON_HOW_THE_ANGLE_APPEARS_AT_ZERO.md, the minus-zero-plus-zero of
/// reflections/ON_BOTH_SIDES_OF_THE_MIRROR.md), while the one-side block keeps it: the both-sides memory met
/// again through the topology door.</para>
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
        Console.WriteLine("#       topology-independent (same mirror in every room). antiU = linear wherever the full block's spectrum is");
        Console.WriteLine("#       conjugation-closed (≈ 0 at N=4; the odd ring is the exception). That agreement is");
        Console.WriteLine("#       NOT a defect: where the full block holds both phase-sides (+0/−0), the i↦−i memory (the i^3");
        Console.WriteLine("#       step of the Z_4 loop) has nothing to distinguish: the both-sides collapse of");
        Console.WriteLine("#       ON_BOTH_SIDES_OF_THE_MIRROR. The symmetrised one-side block (foldlift) keeps it (linear≠0).");
        return 0;
    }
}
