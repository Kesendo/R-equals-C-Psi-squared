using System;
using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F87;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>The GF(2)[x] syzygy production criterion (the girth-ladder deg-1 channel test),
/// mirror of simulations/f87_girth_dichotomy.py Block 5. A k-letter window template carries two
/// mask polynomials over the slot index j: p (bit j ⟺ letter_j ∈ {X, Y}) and q (bit j ⟺
/// letter_j ∈ {Y, Z}); the pair's syzygy B = (p2/g)·q1 + (p1/g)·q2 (g = gcd(p1, p2), all
/// carry-less) is a MONOMIAL ⟺ the single-site-Z production channel exists at the mask level.
/// The deg-1 verdict additionally needs the y_par = 0 face: at y_par = 1 word reversal is a
/// sign-reversing involution that kills every odd girth moment regardless of production.
///
/// <para>Pinned anchors (canonical pairs from simulations/f87_windowed_monomial_converse.py):
/// the K3 triangle XXZ+XZX has B = x + x² + x³ (no production, deg-3 branch, m* = 9);
/// IXXZ+XIXZ has B = x³ and fires p₇ = 573440·γ; the flux pair IXY+XIY and the γ⁵ rung
/// IIXY+ZXZY have monomial B but sit on the y_par = 1 face, so the channel never converts.
/// Exhaustive cross-checks pin the Python oracle counts at k = 3 (0 of 76) and k = 4 (28).</para></summary>
public class SyzygyBTests
{
    private static PauliLetter[] L(string s) =>
        s.Select(PauliLetterExtensions.FromSymbol).ToArray();

    [Fact]
    public void PolyMul_IsTheCarrylessGF2Product()
    {
        Assert.Equal(0b101UL, WindowedObstructionScan.PolyMul(0b11, 0b11));      // (1+x)² = 1+x²
        Assert.Equal(0b1100UL, WindowedObstructionScan.PolyMul(0b11, 0b100));    // (1+x)·x² = x²+x³
        Assert.Equal(0UL, WindowedObstructionScan.PolyMul(0b1011, 0));
        Assert.Equal(0b1011UL, WindowedObstructionScan.PolyMul(0b1011, 1));

        // Round trip against the existing quotient: (a·b)/b = a for exact products.
        foreach (var (a, b) in new (ulong, ulong)[] { (0b101, 0b11), (0b1101, 0b111), (0b1, 0b1011) })
            Assert.Equal(a, WindowedObstructionScan.PolyDivQuotient(
                WindowedObstructionScan.PolyMul(a, b), b));
    }

    [Fact]
    public void IsMonomial_TrueIffExactlyOneBitSet()
    {
        Assert.True(WindowedObstructionScan.IsMonomial(1));
        Assert.True(WindowedObstructionScan.IsMonomial(0b1000));
        Assert.True(WindowedObstructionScan.IsMonomial(1UL << 40));
        Assert.False(WindowedObstructionScan.IsMonomial(0));
        Assert.False(WindowedObstructionScan.IsMonomial(0b11));
        Assert.False(WindowedObstructionScan.IsMonomial(0b1110));
    }

    /// <summary>The K3 triangle (k = 3 XXZ-family cell triangle, K3_EVEN of the Python anchor):
    /// B = x + x² + x³, NOT monomial. No production channel, which is why the k = 3 cell had no
    /// deg-1 pure cycles at all (the old "deg = 1 only for a single-site-Z lift" taxonomy was a
    /// k = 3 cell fact). K3 fires on the deg-3 branch: m* = 9 = 2·3+3, pure γ³.</summary>
    [Fact]
    public void SyzygyB_K3Triangle_IsNotMonomial()
    {
        var (p1, q1, y1) = WindowedObstructionScan.TemplateMasks(L("XXZ"));
        var (p2, q2, y2) = WindowedObstructionScan.TemplateMasks(L("XZX"));
        Assert.Equal((0b011UL, 0b100UL, 0), (p1, q1, y1));
        Assert.Equal((0b101UL, 0b010UL, 0), (p2, q2, y2));

        ulong b = WindowedObstructionScan.SyzygyB(p1, q1, p2, q2);
        Assert.Equal(0b1110UL, b);                                   // x + x² + x³
        Assert.False(WindowedObstructionScan.IsMonomial(b));
        Assert.False(WindowedObstructionScan.Deg1ChannelOpen(L("XXZ"), L("XZX")));
    }

    /// <summary>IXXZ+XIXZ (k = 4): B = x³ IS monomial and the pair sits on the y_par = 0 face, so
    /// the deg-1 channel is open and fires: m* = 7 = 2·3+1, p₇ = 573440·γ (girth dichotomy
    /// Block 3, the first of the 20 k = 4 deg-1 pure cycles).</summary>
    [Fact]
    public void Deg1ChannelOpen_IXXZplusXIXZ_IsOpen()
    {
        var (p1, q1, _) = WindowedObstructionScan.TemplateMasks(L("IXXZ"));
        var (p2, q2, _) = WindowedObstructionScan.TemplateMasks(L("XIXZ"));
        ulong b = WindowedObstructionScan.SyzygyB(p1, q1, p2, q2);
        Assert.Equal(0b1000UL, b);                                   // x³
        Assert.True(WindowedObstructionScan.IsMonomial(b));
        Assert.True(WindowedObstructionScan.Deg1ChannelOpen(L("IXXZ"), L("XIXZ")));
    }

    /// <summary>The y_par = 1 face is closed regardless of production: word reversal is a
    /// sign-reversing involution killing every odd t_j. IIXY+ZXZY has monomial B yet fires only
    /// at γ⁵ (m* = 11, coefficient 86507520); the flux pair IXY+XIY has monomial B yet t_3 = 0
    /// (production open but sign-cancelled, the Cartier-Foata kill).</summary>
    [Fact]
    public void Deg1ChannelOpen_YParityOneFace_IsClosedRegardlessOfProduction()
    {
        var (p1, q1, y1) = WindowedObstructionScan.TemplateMasks(L("IIXY"));
        var (p2, q2, y2) = WindowedObstructionScan.TemplateMasks(L("ZXZY"));
        Assert.Equal((1, 1), (y1, y2));
        Assert.True(WindowedObstructionScan.IsMonomial(
            WindowedObstructionScan.SyzygyB(p1, q1, p2, q2)));       // B = x: production is open
        Assert.False(WindowedObstructionScan.Deg1ChannelOpen(L("IIXY"), L("ZXZY")));

        var (f1, g1, _) = WindowedObstructionScan.TemplateMasks(L("IXY"));
        var (f2, g2, _) = WindowedObstructionScan.TemplateMasks(L("XIY"));
        Assert.True(WindowedObstructionScan.IsMonomial(
            WindowedObstructionScan.SyzygyB(f1, g1, f2, g2)));       // B = x²
        Assert.False(WindowedObstructionScan.Deg1ChannelOpen(L("IXY"), L("XIY")));
    }

    [Fact]
    public void Deg1ChannelOpen_DifferentYParity_IsClosed()
    {
        // XXZ is y_par = 0, XYI is y_par = 1: not y_par-homogeneous, channel closed by definition.
        Assert.False(WindowedObstructionScan.Deg1ChannelOpen(L("XXZ"), L("XYI")));
    }

    /// <summary>The degenerate faces. One diagonal template (p = 0): B reduces to that template's
    /// q-mask, the single-site-Z lift channel; IIZ+XXZ is open (single-site Z, the lift fires
    /// deg-1 at m* = 3) while ZZZ+XXZ is not (the MULTIZ multi-Z lift: t_1 = 0, deg 3, m* = 5).
    /// The fully diagonal pair (p1 = p2 = 0) has no syzygy and throws, mirroring the Python
    /// anchor's syzygy_B returning None.</summary>
    [Fact]
    public void SyzygyB_DiagonalFaces_ReduceToTheQLiftChannel()
    {
        var (pd, qd, _) = WindowedObstructionScan.TemplateMasks(L("IIZ"));
        var (pm, qm, _) = WindowedObstructionScan.TemplateMasks(L("XXZ"));
        Assert.Equal(0UL, pd);
        Assert.Equal(qd, WindowedObstructionScan.SyzygyB(pd, qd, pm, qm));   // p1 = 0 face: B = q1
        Assert.Equal(qd, WindowedObstructionScan.SyzygyB(pm, qm, pd, qd));   // p2 = 0 face: B = q2
        Assert.True(WindowedObstructionScan.Deg1ChannelOpen(L("IIZ"), L("XXZ")));
        Assert.True(WindowedObstructionScan.Deg1ChannelOpen(L("XXZ"), L("IIZ")));

        var (pz, qz, _) = WindowedObstructionScan.TemplateMasks(L("ZZZ"));
        Assert.Equal(qz, WindowedObstructionScan.SyzygyB(pz, qz, pm, qm));
        Assert.False(WindowedObstructionScan.IsMonomial(qz));                // 0b111: multi-Z lift
        Assert.False(WindowedObstructionScan.Deg1ChannelOpen(L("ZZZ"), L("XXZ")));

        Assert.Throws<ArgumentException>(() => WindowedObstructionScan.SyzygyB(0, qd, 0, qz));
        Assert.Throws<ArgumentException>(
            () => WindowedObstructionScan.Deg1ChannelOpen(L("IIZ"), L("IZI")));
    }

    /// <summary>Exhaustive k = 3 cross-check against the Python oracle (one-shot run 2026-06-10
    /// over simulations/f87_girth_dichotomy.py syzygy_B/_is_monomial): 16 diagonal-cell templates,
    /// 76 same-y-par pairs, 42 both-nondiagonal; NO y_par = 0 pure-cycle (both-nondiagonal) pair
    /// has B monomial (0 of 21; the k = 3 cell fact behind the old deg-1 taxonomy), while the
    /// y_par = 1 face carries 4 monomial-B pairs (all reversal-killed, flux among them). Soft
    /// pairs (equal (1+x)-valuations) always have B(1) = 0, i.e. even popcount.</summary>
    [Fact]
    public void Exhaustive_K3_NoYParZeroPureCyclePairHasMonomialB()
    {
        var counts = EnumerateAndCount(k: 3);
        Assert.Equal(16, counts.Cells);
        Assert.Equal(76, counts.SameYParPairs);
        Assert.Equal(42, counts.NondiagonalPairs);
        Assert.Equal(21, counts.YPar0Pairs);
        Assert.Equal(0, counts.YPar0MonomialB);
        Assert.Equal(21, counts.YPar1Pairs);
        Assert.Equal(4, counts.YPar1MonomialB);
        Assert.Equal(0, counts.ChannelOpenPairs);   // no k=3 pure-cycle pair opens the channel
    }

    /// <summary>Exhaustive k = 4 cross-check against the Python oracle (one-shot run 2026-06-10
    /// over simulations/f87_girth_dichotomy.py syzygy_B/_is_monomial on all same-y-par cell
    /// pairs): 64 templates, 1056 same-y-par pairs, 828 both-nondiagonal; exactly 28 of the 300
    /// y_par = 0 nondiagonal pairs have B monomial (= Deg1ChannelOpen). The same oracle run
    /// classified the 28: all are saturated-hard, 20 fire deg-1 at N = 5 (the girth-dichotomy
    /// Block 3 pure cycles) and 8 are window-starved (IXZX+XIZX-type, firing from N = 6). The
    /// y_par = 1 face carries 80 monomial-B pairs, every one reversal-killed.</summary>
    [Fact]
    public void Exhaustive_K4_MonomialBCounts_MatchPythonOracle()
    {
        var counts = EnumerateAndCount(k: 4);
        Assert.Equal(64, counts.Cells);
        Assert.Equal(1056, counts.SameYParPairs);
        Assert.Equal(828, counts.NondiagonalPairs);
        Assert.Equal(300, counts.YPar0Pairs);
        Assert.Equal(28, counts.YPar0MonomialB);
        Assert.Equal(528, counts.YPar1Pairs);
        Assert.Equal(80, counts.YPar1MonomialB);
        Assert.Equal(28, counts.ChannelOpenPairs);  // the wrapper opens exactly the y_par=0 monomials
    }

    private readonly record struct SyzygyCensus(
        int Cells, int SameYParPairs, int NondiagonalPairs,
        int YPar0Pairs, int YPar0MonomialB, int YPar1Pairs, int YPar1MonomialB,
        int ChannelOpenPairs);

    /// <summary>Enumerate all diagonal-cell (Klein (0,1): #(X+Y) even, #(Y+Z) odd) k-letter
    /// templates including the diagonal (p = 0) ones, pair them y_par-homogeneously (unordered,
    /// self-pairs included, exactly the Python oracle's combinations_with_replacement), and count
    /// the monomial-B faces (B and the channel verdict only over both-nondiagonal pairs, matching
    /// the oracle). Also asserts the two structural facts on the way: a monomial B implies
    /// <see cref="WindowedObstructionScan.IsHardPair"/>, and a soft pair's B has even popcount
    /// (B(1) = 0, never monomial).</summary>
    private static SyzygyCensus EnumerateAndCount(int k)
    {
        var lut = new[] { PauliLetter.I, PauliLetter.X, PauliLetter.Y, PauliLetter.Z };
        var cells = new List<PauliLetter[]>();
        for (long code = 0; code < 1L << (2 * k); code++)
        {
            var letters = new PauliLetter[k];
            int na = 0, nb = 0;
            for (int i = 0; i < k; i++)
            {
                letters[i] = lut[(int)((code >> (2 * i)) & 3)];
                na += letters[i].BitA();
                nb += letters[i].BitB();
            }
            if (na % 2 == 0 && nb % 2 == 1) cells.Add(letters);      // diagonal cell, incl. p = 0
        }

        int sameYPar = 0, nondiag = 0, y0 = 0, y0Mono = 0, y1 = 0, y1Mono = 0, open = 0;
        for (int a = 0; a < cells.Count; a++)
            for (int b = a; b < cells.Count; b++)
            {
                var (p1, q1, yp1) = WindowedObstructionScan.TemplateMasks(cells[a]);
                var (p2, q2, yp2) = WindowedObstructionScan.TemplateMasks(cells[b]);
                if (yp1 != yp2) continue;
                sameYPar++;
                if (p1 == 0 || p2 == 0) continue;                    // oracle skips diagonal faces
                nondiag++;
                if (WindowedObstructionScan.Deg1ChannelOpen(cells[a], cells[b])) open++;
                ulong bPoly = WindowedObstructionScan.SyzygyB(p1, q1, p2, q2);
                bool mono = WindowedObstructionScan.IsMonomial(bPoly);
                if (mono)
                    Assert.True(WindowedObstructionScan.IsHardPair(p1, p2),
                        "monomial B on a pair with equal (1+x)-valuations");
                else if (!WindowedObstructionScan.IsHardPair(p1, p2))
                    Assert.True((System.Numerics.BitOperations.PopCount(bPoly) & 1) == 0,
                        "soft pair with odd-popcount B (B(1) != 0)");
                if (yp1 == 0) { y0++; if (mono) y0Mono++; }
                else { y1++; if (mono) y1Mono++; }
            }
        return new SyzygyCensus(cells.Count, sameYPar, nondiag, y0, y0Mono, y1, y1Mono, open);
    }
}
