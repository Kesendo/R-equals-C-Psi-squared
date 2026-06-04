using System.Collections.Generic;
using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F87;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>The windowed F87 hardness obstruction is uniformly a triangle. The §7.5/§7.6 proof made
/// soft ⟺ bipartite a theorem for any k, leaving only the SHAPE of the minimal odd 𝔽₂-relation
/// k-specific. <see cref="WindowedObstructionScan"/> computes that shape in pure GF(2) bit
/// arithmetic (no Hamiltonian), so it reaches body counts the Liouvillian route cannot. These tests
/// (1) reproduce the F103 §6 k=3 anchor, (2) cross-validate the pure-mask scan against the canonical
/// H-based bipartite classifier <see cref="BipartiteChirality"/>, and (3) confirm the obstruction is
/// a size-3 triangle for every hard pair at k = 3..8.</summary>
public class WindowedObstructionScanTests
{
    [Fact]
    public void Scan_K3N4_ReproducesF103Anchor()
    {
        var r = WindowedObstructionScan.Scan(k: 3, n: 4);
        Assert.Equal(42, r.Pairs);
        Assert.Equal(16, r.Hard);
        Assert.Equal(26, r.Soft);
        Assert.True(r.ObstructionIsAlwaysTriangle);   // every hard pair a triangle
        Assert.Equal(16, r.MinOddCycleSizes[3]);
    }

    // Counts cross-checked against the Python scout simulations/_f87_oddcycle_kscaling.py.
    [Theory]
    [InlineData(3, 4, 42, 16)]
    [InlineData(4, 5, 828, 192)]
    [InlineData(5, 6, 14520, 1792)]
    public void Scan_PairAndHardCounts_MatchPythonScout(int k, int n, int pairs, int hard)
    {
        var r = WindowedObstructionScan.Scan(k, n);
        Assert.Equal(pairs, r.Pairs);
        Assert.Equal(hard, r.Hard);
    }

    /// <summary>The pure-GF(2) mask scan agrees with the canonical H-based bipartite classifier
    /// (<see cref="BipartiteChirality"/>, which builds the rotated hopping graph and 2-colours it)
    /// on every windowed k=3, N=4 diagonal-cell pair: MinOddCycle == 0 ⟺ bipartite ⟺ soft.</summary>
    [Fact]
    public void PureGF2Scan_AgreesWithHBasedBipartite_K3N4()
    {
        const int k = 3, n = 4;
        var chain = new ChainSystem(N: n, J: 1.0, GammaZero: 0.05);
        var withLetters = EnumerateCellTermsWithLetters(k);

        int checked_ = 0;
        for (int a = 0; a < withLetters.Count; a++)
            for (int b = a; b < withLetters.Count; b++)
            {
                var (letters1, mask1, yp1) = withLetters[a];
                var (letters2, mask2, yp2) = withLetters[b];
                if (yp1 != yp2) continue; // y_par-homogeneous

                var s = new HashSet<ulong>();
                for (int w = 0; w <= n - k; w++) { s.Add(mask1 << w); s.Add(mask2 << w); }
                bool gf2Bipartite = WindowedObstructionScan.MinOddCycle(new List<ulong>(s)) == 0;

                var terms = new List<PauliTerm> { new(letters1, Complex.One), new(letters2, Complex.One) };
                bool hBipartite = BipartiteChirality.Classify(chain, terms, PauliLetter.Z).IsBipartite;

                Assert.Equal(hBipartite, gf2Bipartite);
                checked_++;
            }
        Assert.Equal(42, checked_);
    }

    /// <summary>The minimal odd obstruction is a triangle (size 3) for every hard pair, across the
    /// windowed body counts k = 3..8 (N = k+1, the smallest windowed chain). The k ≥ 6 cases are the
    /// new data the Liouvillian route could not reach; pure GF(2) bit arithmetic makes them cheap.</summary>
    [Theory]
    [InlineData(3, 4)]
    [InlineData(4, 5)]
    [InlineData(5, 6)]
    [InlineData(6, 7)]
    [InlineData(7, 8)]
    [InlineData(8, 9)]
    public void Obstruction_IsAlwaysATriangle_Windowed(int k, int n)
    {
        var r = WindowedObstructionScan.Scan(k, n);
        Assert.True(r.Hard > 0, $"expected hard pairs at k={k}, N={n}");
        Assert.True(r.ObstructionIsAlwaysTriangle,
            $"k={k}, N={n}: obstruction sizes {string.Join(",", r.MinOddCycleSizes)}");
    }

    // Parallel to WindowedObstructionScan.CellTerms, retaining the PauliLetter[] for the H-based
    // cross-check. Letters base-4: I=0, X=1, Y=2, Z=3.
    private static List<(PauliLetter[] letters, ulong mask, int yPar)> EnumerateCellTermsWithLetters(int k)
    {
        var outp = new List<(PauliLetter[], ulong, int)>();
        var lut = new[] { PauliLetter.I, PauliLetter.X, PauliLetter.Y, PauliLetter.Z };
        long total = 1L << (2 * k);
        for (long code = 0; code < total; code++)
        {
            int na = 0, nb = 0, ny = 0;
            ulong mask = 0;
            var letters = new PauliLetter[k];
            for (int i = 0; i < k; i++)
            {
                int letter = (int)((code >> (2 * i)) & 3);
                letters[i] = lut[letter];
                bool bitA = letter == 1 || letter == 2;
                bool bitB = letter == 2 || letter == 3;
                if (bitA) { na++; mask |= 1UL << i; }
                if (bitB) nb++;
                if (letter == 2) ny++;
            }
            if (na % 2 == 0 && nb % 2 == 1 && na >= 2)
                outp.Add((letters, mask, ny & 1));
        }
        return outp;
    }
}
