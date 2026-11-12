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
    enum Reader { TheoremB, TheoremCD, CorollaryE2, TheoremE, TheoremG }

    // Patch: the open L × M patch, N = L·M, the templates written as for the torus (rows split by '/')
    sealed record Row(string Name, (string T, long C)[] Templates, int N, bool Ring, bool RandomMagnitudes, Reader By, bool? Page, (int L, int M)? Patch = null);

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
        new("XY across + XX/X on the 3 × 2 patch (Theorem G, breaks from 3 × 2)", new[] { ("XY", 1L), ("XX/X", 1L) }, 6, false, true, Reader.TheoremG, false, (3, 2)),
        new("XY across + XX/X on the 2 × 3 patch (Theorem G, pairs below L* = 3)", new[] { ("XY", 1L), ("XX/X", 1L) }, 6, false, true, Reader.TheoremG, true, (2, 3)),
        new("XX across + XY down on the 3 × 2 patch (Theorem G, c = 0)", new[] { ("XX", 1L), ("X/Y", 1L) }, 6, false, true, Reader.TheoremG, true, (3, 2)),
    };

    /// <summary>The strings of a row placed on its sites: the chain or the ring, or the open patch at every translate that fits.</summary>
    static IEnumerable<(ulong X, ulong Z, long C)> Place(Row row)
    {
        if (row.Patch is not { } p) return SignedMaskBalanceClaim.PlaceStrings(row.Templates, row.N, row.Ring).Select(s => (s.X, s.Z, s.C));
        if (row.N != p.L * p.M) throw new ArgumentException($"{row.Name}: N = {row.N} is not {p.L} × {p.M}");
        var placed = new List<(ulong, ulong, long)>();
        foreach (var (t, c) in row.Templates)
        {
            var (wx, wy) = SignedMaskBalanceClaim.Extent(t);
            for (int b = 0; b + wy <= p.M; b++)
                for (int a = 0; a + wx <= p.L; a++)
                {
                    ulong x = 0, z = 0;
                    foreach (var (i, j, letter) in SignedMaskBalanceClaim.TorusLetters(t))
                    {
                        int site = i + a + p.L * (j + b);
                        if (letter is 'X' or 'Y') x |= 1UL << site;
                        if (letter is 'Y' or 'Z') z |= 1UL << site;
                    }
                    placed.Add((x, z, c));
                }
        }
        return placed;
    }

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
                case Reader.TheoremG:
                    return !SignedMaskBalanceClaim.PatchBreaks(row.Templates[0].T, row.Templates[1].T, row.Patch!.Value.L, row.Patch.Value.M);
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
        var placed = Place(row)
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
            string at = g.Row.Patch is not null ? g.Row.Name : $"{g.Row.Name} at N = {g.Row.N}";
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
        Assert.Equal((10, 8), (paired.Count, broken.Count));
        Assert.True(broken.Min() / paired.Max() > 1e8, $"gap {broken.Min():E2} / {paired.Max():E2}");
        Assert.True(broken.Min() > 0.31, $"{broken.Min():E3}");                                                // measured, quoted on the page
    }

    /// <summary>Theorem H's lemma against the Liouvillian, on the eight broken rows: H scaled to J·H with J = j/max|c|,
    /// j = 0.08, 0.04, 0.02, rates drawn in [0.5, 1.5]. The 2^N eigenvalues of L within 0.5 of 0 (in the blocks |x⟩⟨y|
    /// with x, y in one coset of the flip sets' span) and within 0.5 of −2σ (in the blocks holding the far cells |x⟩⟨x̄|),
    /// each sorted by its real part, against −J²·spec(𝓛) and −2σ + J²·spec(𝓜), 𝓛 and 𝓜 built from H's entries
    /// (<see cref="StrongDephasingEdges"/>); the clusters must hold exactly 2^N eigenvalues each. The deviation over J²
    /// is O(J) by the expansion and is gated to fall at least linearly from one J to the next; measured, it falls as J on
    /// XIY + XZX on the ring of six, where third-order terms move its levels, and as J² on the other seven. The second-order spectra
    /// differ on every row (the break at that order), and their kernels are Theorem 2's two counts. Where no component of
    /// Γ is good, the far edge's nearest eigenvalue d(J) = min |λ + 2σ| is J²·λ_min(𝓜) to a remainder r(J) = d/J² − λ_min
    /// whose fall is measured, not proven: the ratio r(J)/r(J/2) gated near 4, J² (a remainder linear in J would read 2);
    /// and λ_min against the page's closed forms, the weighted code distance on the rows with one string per flip set and,
    /// on XZ + ZY on the chain of four with |c| equal, the two end squares q(1/γ₀, 2/γ₁) + q(1/γ₃, 2/γ₂),
    /// q(u, v) = u + v − √(u² + v²). Where some component is good (X·(X + Y)·X on the ring of six) λ_min(𝓜) = 0 and the far
    /// edge keeps eigenvalues at −2σ, while the clusters still differ. The eigensolver's part, n·eps·‖L‖/J² (n = 4^N),
    /// lies two orders below every deviation it is read against.</summary>
    [Fact]
    public void AtStrongDephasingTheTwoEdgesAreTheSecondOrderGenerators()
    {
        var rng = new Random(20261113);
        double eps = Math.Pow(2, -52);
        int dRows = 0, partialRows = 0, codeRows = 0, endRows = 0;
        var devRatios = new List<(string Row, double Ratio)>();
        var rRatios = new List<double>();
        var linear = new List<string>();
        foreach (var row in Rows.Where(r => r.Page == false))
        {
            int n = row.N, dim = 1 << n, full = dim - 1;
            var placed = Place(row).Select(s => (new PauliMask(s.X, s.Z), (BigInteger)(s.C * (row.RandomMagnitudes ? rng.Next(1, 10) : 1)))).ToList();
            var columns = ComplementConnectionGraph.Columns(n, placed);
            var gamma = Enumerable.Range(0, n).Select(_ => 0.5 + rng.NextDouble()).ToArray();
            double sigma = gamma.Sum(), maxC = placed.Max(p => (double)BigInteger.Abs(p.Item2));
            var (nearOp, farOp) = StrongDephasingEdges.Build(n, columns, gamma);
            double[] Spec(ComplexMatrix m) => m.Evd(MathNet.Numerics.LinearAlgebra.Symmetricity.Hermitian).EigenValues.Select(z => z.Real).OrderBy(v => v).ToArray();
            var nu = Spec(nearOp);
            var mu = Spec(farOp);
            double bound = dim * eps * Math.Max(nu[^1], mu[^1]);
            var (nearCount, farCount, _) = ComplementConnectionGraph.EverySiteDephased(n, columns);
            Assert.True(nu[0] > -bound && mu[0] > -bound, $"{row.Name}: an eigenvalue of 𝓛 or 𝓜 below zero ({nu[0]:E2}, {mu[0]:E2})");
            Assert.True(nu.Count(v => v < 1e-9 * nu[^1]) == nearCount && mu.Count(v => v < 1e-9 * mu[^1]) == farCount,
                        $"{row.Name}: the kernels against Theorem 2's {nearCount}/{farCount}");
            double secondOrderBreak = nu.Zip(mu, (x, y) => Math.Abs(x - y)).Max();
            Assert.True(secondOrderBreak > 0.05, $"{row.Name}: the second-order spectra differ by {secondOrderBreak:E2} only");

            double lam = mu[0];
            if (farCount == 0)
            {
                dRows++;
                Assert.True(lam > 0.05, $"{row.Name}: λ_min(𝓜) = {lam:E3} with no good component");
                if (row.By is Reader.TheoremB or Reader.TheoremCD or Reader.TheoremG)
                {
                    double code = SignedMaskBalanceClaim.StrongDephasingBreak(
                        placed.Select(p => new SignedMaskBalanceClaim.PlacedString(p.Item1.X, p.Item1.Z, (long)p.Item2)).ToList(), gamma);
                    Assert.True(Math.Abs(lam - code) <= bound, $"{row.Name}: λ_min(𝓜) = {lam:R}, code distance {code:R}");
                    codeRows++;
                }
                if (row.Name.StartsWith("XZ + ZY chain"))
                {
                    static double Q(double u, double v) => u + v - Math.Sqrt(u * u + v * v);
                    double ends = Q(1 / gamma[0], 2 / gamma[1]) + Q(1 / gamma[3], 2 / gamma[2]);
                    Assert.True(Math.Abs(lam - ends) <= bound, $"{row.Name}: λ_min(𝓜) = {lam:R}, the two end squares {ends:R}");
                    endRows++;
                }
            }
            else
            {
                partialRows++;
                Assert.True(Math.Abs(lam) <= bound, $"{row.Name}: λ_min(𝓜) = {lam:E2} with {farCount} good components");
            }

            // H itself, the cosets of the flip sets' span, and the blocks holding the two edges
            var H = ComplexMatrix.Build.Dense(dim, dim);
            var flips = new List<int>();
            for (int c = 0; c < dim; c++)
                foreach (var (r, e) in columns[c])
                {
                    H[r, c] = new Complex((double)e.Re, (double)e.Im);
                    if (r != c) flips.Add(r ^ c);
                }
            var basis = new List<int>();
            foreach (int f in flips)
            {
                int v = f;
                foreach (int b in basis) v = Math.Min(v, v ^ b);
                if (v != 0) { basis.Add(v); basis.Sort((p, q) => q.CompareTo(p)); }
            }
            int Reduce(int v) { foreach (int b in basis) v = Math.Min(v, v ^ b); return v; }
            var cosets = Enumerable.Range(0, dim).GroupBy(Reduce).ToDictionary(g => g.Key, g => g.ToArray());
            // the block builder reads site l as bit N − 1 − l, the columns above site l as bit l: here the rates matter
            var bySite = gamma.Reverse().ToArray();

            var devs = new List<double>();
            var rs = new List<double>();
            foreach (double j in new[] { 0.08, 0.04, 0.02 })
            {
                double J = j / maxC, scale = 0;
                var nearCluster = new List<Complex>();
                var farCluster = new List<Complex>();
                foreach (var a in cosets.Values)
                {
                    var partner = cosets[Reduce(a[0] ^ full)];
                    // the block (a, a) holds the populations of a, the block (a, ā's coset) the far cells x ∈ a; where the two
                    // cosets coincide it is one block, read once
                    var bras = partner == a ? new[] { a } : new[] { a, partner };
                    foreach (var bra in bras)
                    {
                        var flat = (from x in a from y in bra select x * dim + y).ToArray();
                        foreach (var l in PerBlockLiouvillianBuilder.BuildBlockZ(J * H, bySite, flat).Evd().EigenValues)
                        {
                            scale = Math.Max(scale, l.Magnitude);
                            if (bra == a && l.Magnitude < 0.5) nearCluster.Add(l);
                            if (bra == partner && (l + 2 * sigma).Magnitude < 0.5) farCluster.Add(l);
                        }
                    }
                }
                Assert.True(nearCluster.Count == dim && farCluster.Count == dim, $"{row.Name} at j = {j}: clusters of {nearCluster.Count} and {farCluster.Count}, not {dim}");
                var nearSorted = nearCluster.OrderByDescending(z => z.Real).ToArray();
                var farSorted = farCluster.OrderBy(z => z.Real).ToArray();
                double dev = 0;
                for (int k = 0; k < dim; k++)
                {
                    dev = Math.Max(dev, (nearSorted[k] + J * J * nu[k]).Magnitude / (J * J));
                    dev = Math.Max(dev, (farSorted[k] + 2 * sigma - J * J * mu[k]).Magnitude / (J * J));
                }
                double solver = dim * dim * eps * scale / (J * J);
                Assert.True(solver < 1e-2 * dev, $"{row.Name} at j = {j}: the solver's part {solver:E2} against the deviation {dev:E2}");
                devs.Add(dev);
                if (farCount == 0)
                {
                    double r = farCluster.Min(l => (l + 2 * sigma).Magnitude) / (J * J) - lam;
                    Assert.True(solver < 1e-2 * Math.Abs(r), $"{row.Name} at j = {j}: the solver's part {solver:E2} against r = {r:E2}");
                    rs.Add(r);
                }
            }
            if (devs[0] / devs[1] < 3 && devs[1] / devs[2] < 3) linear.Add(row.Name);
            for (int k = 0; k < 2; k++)
            {
                double dr = devs[k] / devs[k + 1];
                devRatios.Add((row.Name, dr));
                Assert.True(dr > 1.8, $"{row.Name}: the clusters' deviation falls by {dr:F3} only (devs {string.Join(", ", devs.Select(x => x.ToString("E3")))})");
                if (farCount != 0) continue;
                double rr = rs[k] / rs[k + 1];
                rRatios.Add(rr);
                Assert.True(rr > 3.5 && rr < 4.5, $"{row.Name}: r(J)/r(J/2) = {rr:F4} (r = {string.Join(", ", rs.Select(x => x.ToString("E3")))})");
            }
        }
        Assert.Equal((7, 1, 4, 1), (dRows, partialRows, codeRows, endRows));
        Assert.Equal((16, 14), (devRatios.Count, rRatios.Count));
        // measured: the clusters' deviation over J² falls as J² on seven rows and as J on one, where third-order terms move
        // its levels; the far edge's nearest eigenvalue falls as J² beyond J²·λ_min on all seven rows that have one
        Assert.Equal(new[] { "XIY + XZX ring (Corollary E2, 4 | N)" }, linear);
        Assert.All(devRatios.Where(d => !linear.Contains(d.Row)), d => Assert.True(d.Ratio > 3.5 && d.Ratio < 4.5, $"{d.Row}: {d.Ratio:F3}"));
    }
}
