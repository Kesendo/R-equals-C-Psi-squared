using System.Numerics;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F87;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>The F87 diagonal Klein cell (0,1) under Z dephasing, decided pair by pair by Theorem 2 of
/// <c>docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md</c> (<see cref="PalindromeSoftCertifier.ComplementConnectionAtN"/>:
/// a graph on 2^N bitstrings, exact over the Gaussian integers, no eigensolver), against the counts F106, F103,
/// F110, F111 and F115 state.
///
/// <para>Hardness: at full support (N = k) every hard pair contains a pure-Z template (F111's rule), at k = 4
/// (228, the F106 anchor, read there by the spectral classifier) and at k = 5 (2056). Below full support hard
/// pairs of two mixed templates appear (above the window floor N ≥ 2k − 2 exactly F115's valuation criterion), and "a pure-Z template makes the
/// pair hard" holds at every N: a pure-Z template has an odd Z count, so its windows put a diagonal on H that the
/// complement flips, H_xx ≠ H_x̄x̄ at some x, and condition (ii) of Theorem 2 fails; a mixed template of the cell
/// has no diagonal and cannot compensate. Every count below is a literal, so a change in the graph test or in the
/// enumeration moves it.</para>
///
/// <para>F115's y-parity gate: over every diagonal-cell Mixed string pair, y-parity homogeneous or not, the
/// (1+x)-valuation verdict equals the graph verdict at N = 2k − 2 and 2k − 1, so the criterion needs no y-parity
/// condition; the mismatched pairs are hard as often as the homogeneous ones.</para></summary>
public class DiagonalCellComplementConnectionTests
{
    private readonly ITestOutputHelper _out;
    public DiagonalCellComplementConnectionTests(ITestOutputHelper output) => _out = output;

    private static bool PureZ(PauliTerm t) =>
        t.Letters.All(l => l == PauliLetter.I || l == PauliLetter.Z) && t.Letters.Any(l => l == PauliLetter.Z);

    [Theory]
    [InlineData(2, 2, 3, 0, 3)]
    [InlineData(2, 5, 3, 0, 3)]
    [InlineData(3, 4, 42, 8, 34)]
    [InlineData(3, 6, 42, 8, 34)]
    [InlineData(4, 4, 228, 0, 228)]
    [InlineData(4, 5, 292, 128, 228)]
    [InlineData(4, 6, 356, 320, 228)]
    [InlineData(4, 7, 356, 320, 228)]
    [InlineData(4, 8, 356, 320, 228)]
    [InlineData(5, 5, 2056, 0, 2056)]
    [InlineData(5, 6, 2952, 896, 2056)]
    [InlineData(5, 7, 4488, 2432, 2056)]
    [InlineData(5, 8, 6536, 4480, 2056)]
    public void DiagonalCell_HardPairs_ByYParity_AndThePureZRule(int k, int n, int hardY0, int hardY1, int pureZPairs)
    {
        var items = Z2HomogeneousKBodyEnumeration.Enumerate(k).Where(x => x.Klein == (0, 1)).ToList();
        var rows = items.AsParallel().Select(x =>
        {
            var (palindrome, _, _) = PalindromeSoftCertifier.ComplementConnectionAtN(new[] { x.Term1, x.Term2 }, n);
            return (x.YPar, hard: !palindrome, pure: PureZ(x.Term1) || PureZ(x.Term2));
        }).ToList();
        int h0 = rows.Count(r => r.hard && r.YPar == 0), h1 = rows.Count(r => r.hard && r.YPar == 1);
        _out.WriteLine($"k={k} N={n}: hard (y_par 0, 1) = ({h0}, {h1}); pairs with a pure-Z template {rows.Count(r => r.pure)}");
        Assert.Equal(hardY0, h0);
        Assert.Equal(hardY1, h1);
        // a pure-Z template makes the pair hard, at every N
        Assert.Equal(pureZPairs, rows.Count(r => r.pure));
        Assert.True(rows.Where(r => r.pure).All(r => r.hard), "a pair with a pure-Z template came out palindromic");
        // at full support the converse holds too: every hard pair has a pure-Z template (F111)
        if (n == k) Assert.True(rows.Where(r => r.hard).All(r => r.pure), "a hard pair without a pure-Z template at N = k");
        // below full support pairs of two mixed templates turn hard from k = 3 on; at k <= 2 the mixed masks share one
        // valuation and the split stays pure at every N
        else if (k >= 3) Assert.True(rows.Any(r => r.hard && !r.pure), "below full support no mixed+mixed pair came out hard");
        else Assert.True(rows.Where(r => r.hard).All(r => r.pure), "a mixed+mixed pair came out hard at k <= 2");
    }

    [Theory]
    [InlineData(3, 4, 30, 36, 16, 16)]
    [InlineData(3, 5, 30, 36, 16, 16)]
    [InlineData(4, 6, 772, 768, 448, 448)]
    [InlineData(4, 7, 772, 768, 448, 448)]
    [InlineData(5, 8, 14280, 14400, 8960, 8960)]
    [InlineData(5, 9, 14280, 14400, 8960, 8960)]
    public void F115Valuation_AgreesWithTheGraph_InBothYParityHalves(int k, int n, int pairsSame, int pairsMismatched,
        int hardSame, int hardMismatched)
    {
        var alphabet = new[] { PauliLetter.I, PauliLetter.X, PauliLetter.Y, PauliLetter.Z };
        var strings = new List<(PauliLetter[] L, ulong Mask, int YPar)>();
        for (long code = 0; code < 1L << (2 * k); code++)
        {
            var letters = new PauliLetter[k]; int na = 0, nb = 0, ny = 0; ulong mask = 0;
            for (int i = 0; i < k; i++)
            {
                int c = (int)((code >> (2 * i)) & 3); letters[i] = alphabet[c];
                if (c == 1 || c == 2) { na++; mask |= 1UL << i; }
                if (c == 2 || c == 3) nb++;
                if (c == 2) ny++;
            }
            if (na % 2 == 0 && nb % 2 == 1 && na >= 2) strings.Add((letters, mask, ny & 1));
        }
        var pairs = new List<(int A, int B)>();
        for (int a = 0; a < strings.Count; a++) for (int b = a + 1; b < strings.Count; b++) pairs.Add((a, b));
        var rows = pairs.AsParallel().Select(p =>
        {
            var s1 = strings[p.A]; var s2 = strings[p.B];
            bool valuation = WindowedObstructionScan.IsHardPair(s1.Mask, s2.Mask);
            var (palindrome, _, _) = PalindromeSoftCertifier.ComplementConnectionAtN(
                new[] { new PauliTerm(s1.L, Complex.One), new PauliTerm(s2.L, Complex.One) }, n);
            return (same: s1.YPar == s2.YPar, valuation, hard: !palindrome);
        }).ToList();
        var same = rows.Where(r => r.same).ToList(); var mismatched = rows.Where(r => !r.same).ToList();
        _out.WriteLine($"k={k} N={n}: homogeneous {same.Count} pairs, {same.Count(r => r.hard)} hard; " +
                       $"mismatched {mismatched.Count} pairs, {mismatched.Count(r => r.hard)} hard");
        Assert.Equal(pairsSame, same.Count);
        Assert.Equal(pairsMismatched, mismatched.Count);
        Assert.Equal(hardSame, same.Count(r => r.hard));
        Assert.Equal(hardMismatched, mismatched.Count(r => r.hard));
        Assert.True(rows.All(r => r.valuation == r.hard), "the (1+x)-valuation and the graph disagree on a pair");
    }
}
