using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>F173 against the Liouvillian itself, at random rates. The battery checks the signed mask balance against
/// Theorem 2 of the complement connection on the bitstrings or against the balance; here the
/// spectrum of L = −i[H, ·] + Σ_l γ_l (Z_l·Z_l − ·) is computed with random, non-uniform rates γ_l and read for the
/// pairing about −σ, σ = Σ_l γ_l, and for the two ends (the eigenvalues at 0 and at −2σ), beside the page's closed
/// form and Theorem 2's two counts. Theorem 2 sees no rate (F158's ends are the dark and the lit strings at every
/// strictly positive profile), so a rate profile is a free input of the check: the rates are drawn per row, and which
/// site gets which rate does not matter.
///
/// <para>The spectrum is cut into exact blocks, not approximated: |x⟩⟨y| with x ∈ a + W and y ∈ b + W, W the GF(2)
/// span of the flip sets, is invariant under L, so each block is built by <see cref="PerBlockLiouvillianBuilder"/> and
/// diagonalised alone, and all eigenvalues are pooled before the pairing is read (the reflection moves blocks).</para>
///
/// <para>No exact route here (an eigensolver), so the tolerance is an error model: a backward-stable solver on
/// non-defective blocks places each eigenvalue within about n·eps·‖L‖ of the true one (n = 4^N, the whole spectrum,
/// ‖L‖ read as the largest eigenvalue modulus): a generous bound, not a law, whose failures would be false fails. What
/// carries the reading is the separation: a paired row pairs as a multiset (counts within 10⁻⁶·‖L‖ of μ and of −μ
/// equal at every eigenvalue) and to within that bound, a broken row misses by more than 10⁻³·σ, and the gap between
/// the two is read off; the ends are counted within 10⁻⁸·‖L‖ and every other eigenvalue lies more than 10⁻³·σ from
/// them, so the count does not depend on that number. That a broken row is broken as a set, not only in its
/// multiplicities, is measured here and not promised by F158. The control that must fail: the same paired spectra
/// read about a centre moved by 1 % of σ.</para></summary>
public class SignedMaskBalanceLiouvillianTests
{
    enum Reader { TheoremB, TheoremCD, CorollaryE2, TheoremE }

    sealed record Row(string Name, (string T, long C)[] Templates, int N, bool Ring, bool RandomMagnitudes, Reader By, bool? Page);

    static readonly Row[] Rows =
    {
        new("XY ring (Theorem C, breaks at odd N)", new[] { ("XY", 1L) }, 5, true, true, Reader.TheoremCD, false),
        new("XY ring (Theorem C, pairs at even N)", new[] { ("XY", 1L) }, 6, true, true, Reader.TheoremCD, true),
        new("+XIY + XXX chain (Theorem B, N* = 5)", new[] { ("XIY", 1L), ("XXX", 1L) }, 4, false, true, Reader.TheoremB, true),
        new("+XIY + XXX chain (Theorem B, N* = 5)", new[] { ("XIY", 1L), ("XXX", 1L) }, 5, false, true, Reader.TheoremB, false),
        new("XX + XIY ring (Theorem D, breaks on every ring)", new[] { ("XX", 1L), ("XIY", 1L) }, 6, true, true, Reader.TheoremCD, false),
        new("XIY + XZX ring (below Corollary E2's guard, Theorem 2 alone)", new[] { ("XIY", 1L), ("XZX", 1L) }, 4, true, true, Reader.CorollaryE2, null),
        new("XIY + XZX ring (Corollary E2, 4 | N)", new[] { ("XIY", 1L), ("XZX", 1L) }, 6, true, true, Reader.CorollaryE2, false),
        new("XZ + ZY ring (Corollary E2, every ring)", new[] { ("XZ", 1L), ("ZY", 1L) }, 4, true, false, Reader.CorollaryE2, true),
        new("XZ + ZY ring (Corollary E2, every ring)", new[] { ("XZ", 1L), ("ZY", 1L) }, 5, true, false, Reader.CorollaryE2, true),
        new("XZ + ZY chain (the end breaks it)", new[] { ("XZ", 1L), ("ZY", 1L) }, 4, false, false, Reader.TheoremE, false),
        new("X·(X + Y)·X ring (Corollary E2, breaks at 3 | N)", new[] { ("XXX", 1L), ("XYX", 1L) }, 5, true, false, Reader.CorollaryE2, true),
        new("X·(X + Y)·X ring (Corollary E2, breaks at 3 | N)", new[] { ("XXX", 1L), ("XYX", 1L) }, 6, true, false, Reader.CorollaryE2, false),
        new("XXX + YYY ring (Corollary E2, every ring)", new[] { ("XXX", 1L), ("YYY", 1L) }, 6, true, false, Reader.CorollaryE2, true),
        new("(X − Y)⊗(X + Y) ring (Corollary E2, even N)", new[] { ("XX", 1L), ("XY", 1L), ("YX", -1L), ("YY", -1L) }, 5, true, false, Reader.CorollaryE2, false),
        new("(X − Y)⊗(X + Y) ring (Corollary E2, even N)", new[] { ("XX", 1L), ("XY", 1L), ("YX", -1L), ("YY", -1L) }, 6, true, false, Reader.CorollaryE2, true),
    };

    /// <summary>The typed claim's own reading of a row (pairs = true), null where its reader refuses the row.</summary>
    static bool? Typed(Row row)
    {
        try
        {
            switch (row.By)
            {
                case Reader.TheoremB:
                    var (breaks, from) = SignedMaskBalanceClaim.Chain(row.Templates[0].T, row.Templates[1].T);
                    return !breaks || row.N < from;
                case Reader.TheoremCD:
                    return !SignedMaskBalanceClaim.RingBreaks(row.Templates.Select(t => t.T).ToArray(), row.N);
                case Reader.CorollaryE2:
                    return !SignedMaskBalanceClaim.RingOneShape(row.Templates, row.N).Breaks;
                default:
                    return SignedMaskBalanceClaim.DecideCharacter(SignedMaskBalanceClaim.PlaceStrings(row.Templates, row.N, row.Ring)) switch
                    {
                        SignedMaskBalanceClaim.CharacterVerdict.Pairs => true,
                        SignedMaskBalanceClaim.CharacterVerdict.Breaks => false,
                        _ => null,
                    };
            }
        }
        catch (ArgumentException) { return null; }
    }

    sealed record Reading(Row Row, int Near, int Far, int End0, int End2, double EndGap, double Pairing, int MultiplicityMisses, double Shifted, double Scale, double Sigma, int Dim);

    static Reading Read(Row row, Random rng)
    {
        int n = row.N, dim = 1 << n;
        // one magnitude per placed string where the row asks for it (Theorem A reads the signs alone)
        var placed = SignedMaskBalanceClaim.PlaceStrings(row.Templates, n, row.Ring)
            .Select(s => (new PauliMask(s.X, s.Z), (BigInteger)(s.C * (row.RandomMagnitudes ? rng.Next(1, 10) : 1)))).ToList();
        var columns = ComplementConnectionGraph.Columns(n, placed);
        var (near, far, _) = ComplementConnectionGraph.EverySiteDephased(n, columns);

        var H = ComplexMatrix.Build.Dense(dim, dim);
        var flips = new List<int>();
        for (int c = 0; c < dim; c++)
            foreach (var (r, e) in columns[c])
            {
                H[r, c] = new Complex((double)e.Re, (double)e.Im);
                if (r != c) flips.Add(r ^ c);
            }
        var gamma = Enumerable.Range(0, n).Select(_ => 0.3 + 1.4 * rng.NextDouble()).ToArray();
        double sigma = gamma.Sum();

        // the cosets of W, the GF(2) span of the flip sets
        var basis = new List<int>();
        foreach (int f in flips)
        {
            int v = f;
            foreach (int b in basis) v = Math.Min(v, v ^ b);
            if (v != 0) { basis.Add(v); basis.Sort((p, q) => q.CompareTo(p)); }
        }
        int Reduce(int v) { foreach (int b in basis) v = Math.Min(v, v ^ b); return v; }
        var cosets = Enumerable.Range(0, dim).GroupBy(Reduce).Select(g => g.ToArray()).ToArray();

        var eigs = new List<Complex>();
        foreach (var a in cosets)
            foreach (var b in cosets)
            {
                var flat = (from x in a from y in b select x * dim + y).ToArray();
                var block = PerBlockLiouvillianBuilder.BuildBlockZ(H, gamma, flat);
                eigs.AddRange(block.Evd().EigenValues);
            }

        double scale = eigs.Max(l => l.Magnitude);
        double tol = 1e-8 * scale;
        int end0 = eigs.Count(l => l.Magnitude < tol), end2 = eigs.Count(l => (l + 2 * sigma).Magnitude < tol);
        // how far the nearest eigenvalue off the ends lies from them: the ends' tolerance reads nothing else
        double endGap = eigs.Where(l => l.Magnitude >= tol && (l + 2 * sigma).Magnitude >= tol)
                            .Select(l => Math.Min(l.Magnitude, (l + 2 * sigma).Magnitude)).DefaultIfEmpty(double.PositiveInfinity).Min();
        return new Reading(row, near, far, end0, end2, endGap, PairingDistance(eigs, sigma), MultiplicityMisses(eigs, sigma, 1e-6 * scale),
                           PairingDistance(eigs, 1.01 * sigma), scale, sigma, eigs.Count);
    }

    /// <summary>max over the eigenvalues λ of the distance from −λ − 2c to the nearest eigenvalue: zero exactly when
    /// the set (not the multiset) is symmetric about −c; the multiplicities are read by <see cref="MultiplicityMisses"/>.</summary>
    static double PairingDistance(List<Complex> eigs, double centre)
    {
        var mu = eigs.Select(l => l + centre).OrderBy(z => z.Real).ToArray();
        var re = mu.Select(z => z.Real).ToArray();
        double worst = 0;
        foreach (var z in mu)
        {
            var target = -z;
            double best = double.PositiveInfinity;
            int i = Array.BinarySearch(re, target.Real);
            if (i < 0) i = ~i;
            for (int j = i; j < mu.Length && mu[j].Real - target.Real < best; j++) best = Math.Min(best, (mu[j] - target).Magnitude);
            for (int j = i - 1; j >= 0 && target.Real - mu[j].Real < best; j--) best = Math.Min(best, (mu[j] - target).Magnitude);
            worst = Math.Max(worst, best);
        }
        return worst;
    }

    /// <summary>The multiset: for every eigenvalue μ (centred, μ = λ + c), the count of eigenvalues within r of μ against
    /// the count within r of −μ; the number of μ where the two differ. Zero for a multiset symmetric about −c whenever r
    /// lies above the solver's error and below the distances between distinct eigenvalues.</summary>
    static int MultiplicityMisses(List<Complex> eigs, double centre, double r)
    {
        var mu = eigs.Select(l => l + centre).OrderBy(z => z.Real).ToArray();
        var re = mu.Select(z => z.Real).ToArray();
        int Count(Complex at)
        {
            int i = Array.BinarySearch(re, at.Real - r);
            if (i < 0) i = ~i;
            int k = 0;
            for (int j = i; j < mu.Length && mu[j].Real <= at.Real + r; j++) if ((mu[j] - at).Magnitude < r) k++;
            return k;
        }
        return mu.Count(z => Count(z) != Count(-z));
    }

    [Fact]
    public void TheLiouvillianPairsExactlyWhereTheoremTwoAndThePageSay_AtRandomRates()
    {
        var rng = new Random(20261112);
        var readings = Rows.Select(r => Read(r, rng)).ToList();
        double eps = Math.Pow(2, -52);
        var paired = new List<double>();
        var broken = new List<double>();
        foreach (var g in readings)
        {
            string at = $"{g.Row.Name} at N = {g.Row.N}";
            bool thm2 = g.Near == g.Far;
            Assert.True(g.Row.Page == Typed(g.Row), $"{at}: the page says {g.Row.Page}, the typed claim {Typed(g.Row)}");
            if (g.Row.Page is bool page) Assert.True(page == thm2, $"{at}: the page says {page}, Theorem 2 {thm2}");
            Assert.True((g.End0, g.End2) == (g.Near, g.Far), $"{at}: spectral ends {g.End0}/{g.End2}, Theorem 2 {g.Near}/{g.Far}");
            Assert.True(g.EndGap > 1e-3 * g.Sigma, $"{at}: an eigenvalue {g.EndGap:E2} from an end");
            double ratio = g.Pairing / (g.Dim * eps * g.Scale);
            if (thm2)
            {
                Assert.True(ratio < 1, $"{at}: pairing distance {g.Pairing:E2} is {ratio:F3} of n·eps·‖L‖");
                Assert.True(g.MultiplicityMisses == 0, $"{at}: {g.MultiplicityMisses} eigenvalues without a partner of equal multiplicity");
                Assert.True(g.Shifted > 1e-3 * g.Sigma, $"{at}: the moved centre still pairs ({g.Shifted:E2})");   // the control
                paired.Add(g.Pairing);
            }
            else
            {
                Assert.True(g.MultiplicityMisses > 0, $"{at}: broken by Theorem 2 but the multiset pairs");
                Assert.True(g.Pairing > 1e-3 * g.Sigma, $"{at}: broken by Theorem 2 but pairs to {g.Pairing:E2}");
                broken.Add(g.Pairing);
            }
        }
        Assert.Equal((8, 7), (paired.Count, broken.Count));
        Assert.True(broken.Min() / paired.Max() > 1e8, $"gap {broken.Min():E2} / {paired.Max():E2}");
        Assert.True(broken.Min() > 0.31, $"{broken.Min():E3}");                                                // measured, quoted on the page
    }
}
