using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;

namespace RCPsiSquared.Cli.Commands;

/// <summary>Where do the N=4 "zeros" (on-line self-mirror strands) GO for N≥5? The branch-locus self-fold of
/// the (SE,DE)=(w1,w2) block is N=4-only (see `foldlift`), because at N=4 the DE sector is its own Hamming
/// complement. The F89c structural lemma (experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md, Tier 1) gives the right
/// partner: the palindrome's column bit-flip ρ[a,b] → ρ[a,bar(b)] (flip the BRA index) carries n_diff(a,b) ↔
/// n_diff(a,bar(b)) = N − n_diff, so it maps the (SE,DE) block to the (SE, w_{N−2}) block (ket stays w1, bra
/// w2 → w_{N−2}). At N=4, w_{N−2} = w2 = DE → partner is itself (the self-fold). For N≥5 the partner is a
/// DIFFERENT block: (SE,TE)=(w1,w3) at N=5, (w1,w4) at N=6, ...
///
/// <para>This probe builds the (SE,DE) block and its partner (SE,w_{N−2}) from scratch (a general
/// computational-basis coherence-block builder; diagonal −2·HammingDistance = the Absorption Theorem rate,
/// ket excitations hop −2qi, bra excitations +2qi), and tests the CROSS-fold: does spec(SE,DE) fold onto
/// spec(SE,w_{N−2}) under the antiunitary λ ↦ −λ̄ − 2σ about the global centre −σ = −N? If yes, the global
/// palindrome LIFTS to all N as a cross-block mirror; only the N=4 self-fold (partner = self) is special, and
/// the "zeros" become cross-block mirror partners rather than on-line self-mirrors. No monodromy, no rebuild
/// of the path-3 witness: a spectrum check on two exact blocks.</para>
///
/// usage: rcpsi foldcross [--nmax 6] [--q 2]</summary>
public static class FoldCrossCommand
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    // all bitmasks on N sites with exactly w set bits.
    private static List<int> Configs(int n, int w)
    {
        var res = new List<int>();
        for (int m = 0; m < (1 << n); m++)
            if (System.Numerics.BitOperations.PopCount((uint)m) == w) res.Add(m);
        return res;
    }

    // the (wKet, wBra) computational-basis coherence block of L = −i[H,ρ] + Z-dephasing on an N-site chain at
    // real q (γ=1). Diagonal −2·n_diff (Absorption Theorem 2γ·n_diff); ket excitations hop −2qi, bra +2qi
    // (the L = −iHρ + iρH split), nearest-neighbour, Pauli exclusion respected.
    private static Complex[,] BuildBlock(int n, int wKet, int wBra, double q)
    {
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
                if ((kc & (1 << s)) != 0)                                   // ket excitation hops (−2qi)
                    foreach (int s2 in new[] { s - 1, s + 1 })
                        if (s2 >= 0 && s2 < n && (kc & (1 << s2)) == 0)
                            l[index[((kc & ~(1 << s)) | (1 << s2), bc)], col] += new Complex(0, -2 * q);
            for (int s = 0; s < n; s++)
                if ((bc & (1 << s)) != 0)                                   // bra excitation hops (+2qi)
                    foreach (int s2 in new[] { s - 1, s + 1 })
                        if (s2 >= 0 && s2 < n && (bc & (1 << s2)) == 0)
                            l[index[(kc, (bc & ~(1 << s)) | (1 << s2))], col] += new Complex(0, 2 * q);
        }
        return l;
    }

    private static Complex[] Spectrum(Complex[,] l)
        => Matrix<Complex>.Build.DenseOfArray(l).Evd().EigenValues.ToArray();

    public static int Run(string[] args)
    {
        var p = new ArgParser(args);
        int nmax = (int)(p.OptionalDouble("nmax") ?? 6);
        double q = p.OptionalDouble("q") ?? 2;

        Console.WriteLine("# foldcross: the (SE,DE)=(w1,w2) block and its palindrome partner (SE,w_{N-2}) (column bit-flip");
        Console.WriteLine("# rho[a,b] -> rho[a,bar b]); does the GLOBAL palindrome lift as a CROSS-block fold for N>=5?");
        Console.WriteLine($"# q={q.ToString("0.##", Inv)}; fold lambda -> -conj(lambda) - 2*sigma, sigma = midpoint of the two block centroids");
        Console.WriteLine();
        Console.WriteLine("  N   partner   dim(1,2)  centroid(1,2)  centroid(1,N-2)   sigma   self-fold(1,2)   CROSS-fold resid");
        for (int n = 4; n <= nmax; n++)
        {
            int wp = n - 2;
            var s12 = Spectrum(BuildBlock(n, 1, 2, q));
            var sp = Spectrum(BuildBlock(n, 1, wp, q));
            double m12 = s12.Average(z => z.Real), mp = sp.Average(z => z.Real);
            double sigma = -(m12 + mp) / 2;                                 // global fold centre −σ = centroid midpoint

            double Resid(Complex[] a, Complex[] b, double sig) =>
                a.Max(l => { var img = new Complex(-l.Real - 2 * sig, l.Imaginary); return b.Min(s => (s - img).Magnitude); });

            double cross = Resid(s12, sp, sigma);
            double self = Resid(s12, s12, -m12);                            // self-fold about (1,2)'s own centroid

            Console.WriteLine(
                $"  {n}   (1,{wp})    {s12.Length,5}   {m12.ToString("0.0000", Inv),12}   {mp.ToString("0.0000", Inv),13}   " +
                $"{sigma.ToString("0.00", Inv),5}   {self.ToString("E2", Inv),12}   {cross.ToString("E2", Inv),14}");
        }

        Console.WriteLine();
        Console.WriteLine("# Read: CROSS resid ~ 0 with sigma = N  =>  the global palindrome lifts as the");
        Console.WriteLine("#       (SE,DE) <-> (SE,w_{N-2}) cross-fold; the N=4 self-fold is the degenerate w_{N-2}=w2=DE");
        Console.WriteLine("#       case (partner=self). The N=4 on-line 'zeros' become cross-block mirror partners for N>=5.");
        return 0;
    }
}
