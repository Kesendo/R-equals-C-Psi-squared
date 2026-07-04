using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the F87-hardness bloc (adopted 2026-07-04 from PROOF_F103_F87_Z2_CUBED_
// REFINEMENT.md sections 6-7 (F103/F110/F111 + section 7.7 = F115), PROOF_F87_WINDOWED_MONOMIAL_
// CONVERSE.md section 5 (F117), PROOF_F107_TRULY_Y_PARITY_ZERO_PURITY.md and PROOF_F109_MOTHER_
// SOFT_Y_PARITY_ONE_PURITY.md, via the registry entries). MirrorWorld adopts the bloc's
// EIGENSOLVER-FREE certificates, not the spectral classifier: hardness as a GF(2)[x] valuation
// difference, purity as letter parity, and the all-gamma converse's trace face (odd power-sums of
// M = A + gamma Q -- traces, never eigenvalues). The crown pin: the two independent certificates
// (valuation vs traces) agree on the K3 trio.
public class HardnessTests
{
    static readonly World W = new();

    // F102: the cube's third axis. y_par = #Y mod 2 relates to the Klein axes as
    // y_par = (k + (bit_a XOR bit_b)) mod 2 on {X,Y,Z}^k -- equal at even k, flipped at odd k,
    // hence an independent third Z2 once k ranges over both parities (XYZ vs III: same Klein
    // cell (0,0), different y_par).
    [Fact]
    public void The_Cube_Has_An_Independent_Third_Axis()
    {
        foreach (int k in new[] { 2, 3, 4 })
            foreach (var letters in Tuples("XYZ", k))
            {
                var (a, b, y) = Hardness.Cube(letters.ToCharArray());
                Assert.Equal((k + (a ^ b)) % 2, y);
            }
        Assert.Equal((0, 0, 1), Hardness.Cube("XYZ".ToCharArray()));
        Assert.Equal((0, 0, 0), Hardness.Cube("III".ToCharArray()));
    }

    // F107: truly forces y_par = 0 under every dephase letter. The per-letter truly criteria
    // (Z: #Y,#Z even; X: #X,#Y even; Y: #Y,#Z even) all contain "#Y even". All 64 k=3 and
    // 256 k=4 letter sequences, three dephase letters.
    [Fact]
    public void Truly_Forces_YPar_Zero_Under_Every_Dephase_Letter()
    {
        foreach (int k in new[] { 3, 4 })
            foreach (var letters in Tuples("IXYZ", k))
                foreach (char d in "ZXY")
                    if (Hardness.Truly(letters.ToCharArray(), d))
                        Assert.Equal(0, letters.Count(c => c == 'Y') % 2);
    }

    // F109: the mother sector's soft side is y_par = 1 pure. In Klein (0,0) the three letter
    // counts share one parity, so truly = all even (y_par 0) and non-truly = all odd (y_par 1).
    // The non-truly Klein-(0,0) census over {X,Y,Z}^k: 6 sequences at k=3 (the XYZ permutations,
    // 21 pairs-with-self) and 24 at k=4 (300 pairs) -- the F103/F106 mother-soft counts.
    [Fact]
    public void The_Mother_Sector_Splits_All_Even_Truly_Against_All_Odd_YParOne()
    {
        foreach (int k in new[] { 3, 4 })
        {
            int nonTruly = 0;
            foreach (var s in Tuples("IXYZ", k))
            {
                var letters = s.ToCharArray();
                var (a, b, y) = Hardness.Cube(letters);
                if ((a, b) != (0, 0)) continue;
                if (Hardness.Truly(letters, 'Z'))
                    Assert.Equal(0, y);
                else { Assert.Equal(1, y); nonTruly++; }
            }
            int expected = k == 3 ? 6 : 24;
            Assert.Equal(expected, nonTruly);
            Assert.Equal(k == 3 ? 21 : 300, nonTruly * (nonTruly + 1) / 2);   // pairs incl. self
        }
    }

    // F115: hardness is a valuation difference. The (1+x)-adic valuations of the K3 trio's
    // window masks (XXZ -> 1+x, XZX -> 1+x^2 = (1+x)^2, ZXX -> x+x^2 = x(1+x)) are {1, 2, 1};
    // a diagonal-cell pair is hard iff the valuations differ.
    [Fact]
    public void The_Valuation_Difference_Decides_The_K3_Trio()
    {
        Assert.Equal(1, Hardness.Valuation(0b011));                          // 1 + x
        Assert.Equal(2, Hardness.Valuation(0b101));                          // (1+x)^2
        Assert.Equal(1, Hardness.Valuation(0b110));                          // x(1+x)
        Assert.False(Hardness.IsHardPair(0b011, 0b110));                     // equal valuation: soft
        Assert.True(Hardness.IsHardPair(0b011, 0b101));                      // hard
        Assert.True(Hardness.IsHardPair(0b110, 0b101));                      // hard
    }

    // F115: the valuation classes and the closed-form hard count. Even-popcount nonzero k-bit
    // masks split by valuation v = 1..k-1 into classes of size 2^(k-1-v); the cross-class pair
    // count matches (4^(k-1) - 3 2^(k-1) + 2)/3 (OEIS A203241), and the Klein/y-parity dressing
    // 2^(2k-3) gives 448 / 8960 / 158720 at k = 4, 5, 6.
    [Fact]
    public void The_Valuation_Classes_And_The_Hard_Count_Close()
    {
        foreach (int k in new[] { 4, 5, 6 })
        {
            var masks = new List<ulong>();
            for (ulong m = 1; m < (1UL << k); m++)
                if (System.Numerics.BitOperations.PopCount(m) % 2 == 0) masks.Add(m);
            for (int v = 1; v <= k - 1; v++)
                Assert.Equal(1L << (k - 1 - v), masks.Count(m => Hardness.Valuation(m) == v));
            long hardPairs = 0;
            for (int i = 0; i < masks.Count; i++)
                for (int j = i + 1; j < masks.Count; j++)
                    if (Hardness.IsHardPair(masks[i], masks[j])) hardPairs++;
            Assert.Equal(Hardness.HardMaskPairCount(k), hardPairs);
        }
        Assert.Equal(14, Hardness.HardMaskPairCount(4));
        Assert.Equal(448, Hardness.DressedHardCount(4));
        Assert.Equal(8960, Hardness.DressedHardCount(5));
        Assert.Equal(158720, Hardness.DressedHardCount(6));
        Assert.Equal(3, Hardness.ObstructionCeiling(3, 2));                  // the K3 always-triangle
        Assert.Equal(Math.Min(2 * 4 - 1, 2 * 6 - 3), Hardness.ObstructionCeiling(6, 4));
    }

    // F110/F111 (with F105's N-stability): hard pairs live only in the diagonal Klein cell of
    // the dephase letter, with the Y-inversion; the adopted splits are 42:8 at k=3 (identical at
    // N=4 and N=5) and 228:0 at k=N=4 via the pure-D template rule, whose decomposition
    // 36 + 192 + 300 = 528 is the pairs arithmetic of 8 pure and 24 mixed cell members.
    [Fact]
    public void Hard_Sits_On_The_Diagonal_Cell_With_The_Y_Inversion()
    {
        Assert.Equal((0, 1), Hardness.DiagonalCell('Z'));
        Assert.Equal((1, 0), Hardness.DiagonalCell('X'));
        Assert.Equal((1, 1), Hardness.DiagonalCell('Y'));
        Assert.Equal((42, 8), Hardness.HardSplitK3);
        Assert.Equal((36, 192, 300), Hardness.TemplateDecompositionK4);
        Assert.Equal(228, 36 + 192);                                         // the F106 pure split
        Assert.Equal(528, 36 + 192 + 300);
        Assert.Equal(192, 8 * 24);                                           // 8 pure-D, 24 mixed members
        Assert.Equal(36, 8 * 9 / 2);
        Assert.Equal(300, 24 * 25 / 2);
        // pure-D templates carry y_par(D) by construction (the Y-inversion's k=4 mechanism)
        foreach (char d in "ZXY")
            foreach (var s in Tuples("IXYZ", 4))
            {
                var letters = s.ToCharArray();
                if (!Hardness.IsPureTemplate(letters, d)) continue;
                if (Hardness.Cube(letters) is var (a, b, y) && (a, b) == Hardness.DiagonalCell(d))
                    Assert.Equal(d == 'Y' ? 1 : 0, y);
            }
    }

    // F117, the cell-free m = 3 face: for H = Z_0 (a bare single-site-Z component, c_0 = 1) the
    // first odd power-sum of M = A + gamma Q is p_3 = 6 4^N c^2 gamma, exactly -- computed here
    // from raw matrix TRACES, no eigensolver, and pinned against the closed form at N = 3.
    [Fact]
    public void The_CellFree_Face_Reads_Six_FourToTheN_From_Traces()
    {
        var hardness = new Hardness(W);
        foreach (double gamma in new[] { 0.5, 1.0 })
        {
            double p3 = hardness.OddPowerSums(new[] { "ZII".ToCharArray() }, n: 3, j: 1.0, gamma, upToOdd: 3)[1];
            Assert.Equal(6.0 * 64.0 * gamma, p3, 6);                         // P_{3,1} = 6 4^N sum c_l^2
        }
    }

    // The crown: the two eigensolver-free certificates agree on the K3 trio at N = 4, and the
    // trace face lands on the registry's exact CRT integer. The soft pair (XXZ, ZXX) has every
    // odd power-sum vanishing through m = 9; the two hard pairs are SILENT on the deg-1 girth
    // ladder (p_7 = 0 -- silence is not softness) and fire at m* = 2*girth + 3 = 9 through the
    // d = 3 class: p_9 = 2064384 gamma^3 exactly, the F117 K3 coefficient, a single positive
    // monomial (the derived monomiality), so hard at EVERY gamma. Traces against GF(2), no
    // spectrum anywhere, and the verdicts match the valuation criterion pair for pair.
    [Fact]
    public void The_Traces_And_The_Valuation_Agree_On_The_K3_Trio()
    {
        var hardness = new Hardness(W);
        var xxz = "XXZ".ToCharArray(); var xzx = "XZX".ToCharArray(); var zxx = "ZXX".ToCharArray();

        var soft = hardness.OddPowerSums(new[] { xxz, zxx }, n: 4, j: 1.0, gamma: 0.5, upToOdd: 9);
        for (int i = 0; i < soft.Length; i++)
            Assert.True(Math.Abs(soft[i]) < 1e-4, $"soft pair: p_{2 * i + 1} = {soft[i]:E2} must vanish");

        foreach (double gamma in new[] { 0.5, 1.0 })
            foreach (var pair in new[] { new[] { xxz, xzx }, new[] { xzx, zxx } })
            {
                var p = hardness.OddPowerSums(pair, n: 4, j: 1.0, gamma, upToOdd: 9);
                for (int i = 0; i < 4; i++)
                    Assert.True(Math.Abs(p[i]) < 1e-4, $"hard pair: p_{2 * i + 1} = {p[i]:E2} must vanish below m* = 9");
                Assert.Equal(2064384.0 * gamma * gamma * gamma, p[4], 3);    // the F117 K3 integer, exactly
            }
    }

    static IEnumerable<string> Tuples(string alphabet, int k)
        => k == 0 ? new[] { "" } : Tuples(alphabet, k - 1).SelectMany(t => alphabet.Select(c => t + c));
}
