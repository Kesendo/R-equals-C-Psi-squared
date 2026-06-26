using System;
using System.Globalization;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Cli.Commands;

/// <summary>Does the F89 branch-locus palindrome LIFT beyond path-3? The fold Re λ = −σ on the path-3
/// (SE,DE) octic is the antiunitary symmetry T·L(q)·T⁻¹ = −L(q̄)−2σ, claimed topology/N-general in
/// experiments/F89_BRANCH_LOCUS_PALINDROME.md (the F1 palindrome + weight-complement P are topology-general).
/// This probe FEEDS the existing exact block builder F89PathKSeDeBlock.BuildTwoTimesSymBlock(q, nBlock) for
/// path-k = nBlock−1, computes the block spectrum at a real integer q, and applies the same two-sided gate as
/// BranchLocusPalindromeWitness: the spectrum must close under the ANTIUNITARY mirror λ↦−λ̄−2σ (residual ≈ 0)
/// and NOT under the bare LINEAR one λ↦−λ−2σ (Im flipped; residual large). It reports the centre σ = −mean(Re λ),
/// the two residuals, and the count of on-line "zeros" (σ_T-fixed strands, Re λ = −σ) vs mirror pairs vs orphans.
/// No monodromy / no rebuild of the path-3 witness: a spectrum check on the trusted Z[i] block.
///
/// usage: rcpsi foldlift [--nmax 7] [--q 2] [--tol 1e-6]</summary>
public static class FoldLiftCommand
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    public static int Run(string[] args)
    {
        var p = new ArgParser(args);
        int nmax = (int)(p.OptionalDouble("nmax") ?? 7);
        int q = (int)(p.OptionalDouble("q") ?? 2);
        double tol = p.OptionalDouble("tol") ?? 1e-6;

        Console.WriteLine($"# foldlift: does the F89 branch-locus palindrome (fold Re λ = −σ) lift beyond path-3?");
        Console.WriteLine($"# (SE,DE) block at q={q}, λ = eig(2M)/2; two-sided gate antiunitary λ↦−λ̄−2σ vs linear λ↦−λ−2σ");
        Console.WriteLine();
        Console.WriteLine($"  path-k   N  dim   σ=−mean(Re)   antiU resid    linear resid   zeros(on-line)  mirror-pairs  orphans");

        for (int nBlock = 4; nBlock <= nmax; nBlock++)
        {
            var g = F89PathKSeDeBlock.BuildTwoTimesSymBlock(q0: q, nBlock: nBlock);
            int d = g.GetLength(0);
            var m = Matrix<Complex>.Build.Dense(d, d, (r, c) =>
                new Complex((double)g[r, c].Re, (double)g[r, c].Im));
            var lam = m.Evd().EigenValues.Select(v => v / 2).ToArray();   // physical λ_k

            double sigma = -lam.Average(z => z.Real);                     // centroid: mean(Re) = −σ if the fold holds

            Complex AntiU(Complex z) => new(-z.Real - 2 * sigma, z.Imaginary);
            Complex Lin(Complex z) => new(-z.Real - 2 * sigma, -z.Imaginary);
            double Residual(Func<Complex, Complex> f) =>
                lam.Max(l => { var img = f(l); return lam.Min(s => (s - img).Magnitude); });

            double anti = Residual(AntiU);
            double lin = Residual(Lin);

            int zeros = lam.Count(l => Math.Abs(l.Real + sigma) < tol);

            // mirror-pair / orphan accounting on the off-line strands (greedy match under AntiU).
            var off = lam.Where(l => Math.Abs(l.Real + sigma) >= tol).ToList();
            var used = new bool[off.Count];
            int pairs = 0, orphans = 0;
            for (int i = 0; i < off.Count; i++)
            {
                if (used[i]) continue;
                var img = AntiU(off[i]);
                int best = -1; double bd = double.PositiveInfinity;
                for (int j = 0; j < off.Count; j++)
                    if (j != i && !used[j]) { double dd = (off[j] - img).Magnitude; if (dd < bd) { bd = dd; best = j; } }
                if (best >= 0 && bd < 1e-6) { used[i] = used[best] = true; pairs++; }
                else { used[i] = true; orphans++; }
            }

            Console.WriteLine(
                $"  path-{nBlock - 1}   {nBlock,1}  {d,3}   {sigma.ToString("0.0000", Inv),10}   " +
                $"{anti.ToString("E2", Inv),11}   {lin.ToString("0.000", Inv),11}   {zeros,12}   {pairs,12}   {orphans,8}");
        }

        Console.WriteLine();
        Console.WriteLine("# Read: antiU resid ≈ 0 with linear resid large  ⟹  the fold lifts (spectrum mirror about −σ,");
        Console.WriteLine("#       Re reflected, Im preserved). zeros = σ_T-fixed strands on the fold; the candidate '0's");
        Console.WriteLine("#       whose connection-through-twins is the next, harder (monodromy) question.");
        return 0;
    }
}
