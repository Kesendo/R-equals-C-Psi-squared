using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>From-below pins for <see cref="ComplementConnectionWitness"/> and
/// <see cref="PalindromeComplementConnectionClaim"/> (proof
/// <c>docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md</c>). The graph route is judged against
/// <see cref="PalindromeStringSpanWitness"/>, which ranks commutators on string spans and never forms a
/// basis state, and against Theorems 1 and 3's rules, which carry no graph at all. The five N = 3 rows
/// with every site dephased are the string witness's own test rows, with its counts.</summary>
public class ComplementConnectionWitnessTests
{
    [Theory]
    [InlineData("ZZZ", "...", 4, 4, true)]
    [InlineData("ZZZ", "Z..", 4, 0, false)]
    [InlineData("ZZZ", "X.X", 1, 1, true)]
    [InlineData("XZZ", "..Y", 1, 1, true)]
    [InlineData("XYZ", "...", 1, 0, false)]
    public void The_Graph_Meets_The_String_Route_On_Its_Own_Rows(string deph, string field, int near, int far, bool pal)
    {
        var r = new ComplementConnectionWitness(3, deph, field).Read();
        var s = new PalindromeStringSpanWitness(3, deph, field).Read();
        Assert.Equal((near, far, pal), (r.HoppingComponents, r.Good, r.Palindrome));
        Assert.Equal((s.NearUpper, s.FarUpper), (r.HoppingComponents, r.Good));
        Assert.Equal(pal, r.Predicted);
    }

    // Every axis assignment and every letter field pattern at N = 3, on the chain and the ring: the graph's
    // verdict is the theorem's rule on every row, the near count is 1 wherever the axes are mixed (Lemma A),
    // and every seventh row is read beside the string route, counts and verdict both.
    [Theory]
    [InlineData("chain")]
    [InlineData("ring")]
    public void Every_N3_Row_Follows_Theorem_One_Or_Three(string topology)
    {
        int rows = 0, pal = 0, mixedBroken = 0;
        foreach (int a in Enumerable.Range(0, 27))
        foreach (int f in Enumerable.Range(0, 64))
        {
            string deph = new(Enumerable.Range(0, 3).Select(l => "XYZ"[a / Pow3(l) % 3]).ToArray());
            string field = new(Enumerable.Range(0, 3).Select(l => ".XYZ"[f >> (2 * l) & 3]).ToArray());
            var r = new ComplementConnectionWitness(3, deph, field, topology).Read();
            Assert.True(r.Predicted == r.Palindrome, $"{deph} {field} {topology}: {r}");
            bool mixed = deph.Distinct().Count() > 1;
            if (mixed)
            {
                Assert.Equal(1, r.HoppingComponents);
                Assert.Equal(PalindromeComplementConnectionClaim.MixedAxesPalindrome(deph, field.Where(c => c != '.')),
                    r.Palindrome);
                if (!r.Palindrome) { mixedBroken++; Assert.Equal(0, r.Good); }
            }
            if (rows++ % 7 == 0)
            {
                var s = new PalindromeStringSpanWitness(3, deph, field, topology).Read();
                Assert.Equal((s.NearUpper, s.FarUpper), (r.HoppingComponents, r.Good));
                Assert.Equal(PalindromeStringSpanWitness.IsPalindrome(s.Verdict), r.Palindrome);
            }
            if (r.Palindrome) pal++;
        }
        Assert.True(pal > 100 && mixedBroken > 100, $"palindromic {pal}, mixed broken {mixedBroken}");
    }

    private static int Pow3(int l) => l == 0 ? 1 : 3 * Pow3(l - 1);

    // Past the string route's comfortable size, a few rows at N = 6 on the ring, read beside it.
    [Theory]
    [InlineData("ZXZXZX", "YYY.Y.", true)]
    [InlineData("ZXZXZX", "YYX...", false)]
    [InlineData("XYZXYZ", "......", false)]
    [InlineData("YYYYYY", "X..X..", true)]
    [InlineData("YYYYYY", "X..Z..", false)]
    [InlineData("YYYYYY", "X..Y..", false)]
    public void Six_Sites_On_The_Ring_Meet_The_String_Route(string deph, string field, bool pal)
    {
        var r = new ComplementConnectionWitness(6, deph, field, "ring").Read();
        var s = new PalindromeStringSpanWitness(6, deph, field, "ring").Read();
        Assert.Equal(pal, r.Palindrome);
        Assert.Equal(pal, r.Predicted);
        Assert.Equal((s.NearUpper, s.FarUpper), (r.HoppingComponents, r.Good));
    }

    // Two sites (the ring doubles the bond, as twoend and twoendstrings do) and the complete graph at N = 4.
    [Theory]
    [InlineData(2, "ring", "XY", "ZZ", true)]
    [InlineData(2, "ring", "XY", "Z.", true)]
    [InlineData(2, "chain", "XY", "X.", false)]
    [InlineData(2, "chain", "ZZ", "..", true)]
    [InlineData(4, "complete", "XXYY", "ZZ.Z", true)]
    [InlineData(4, "complete", "XXYZ", "....", false)]
    [InlineData(4, "complete", "ZZZZ", "X..Y", false)]
    public void Small_And_Complete_Rows_Meet_The_String_Route(int n, string topology, string deph, string field, bool pal)
    {
        var r = new ComplementConnectionWitness(n, deph, field, topology).Read();
        var s = new PalindromeStringSpanWitness(n, deph, field, topology).Read();
        Assert.Equal((pal, pal), (r.Palindrome, r.Predicted));
        Assert.Equal((s.NearUpper, s.FarUpper), (r.HoppingComponents, r.Good));
        if (deph.Distinct().Count() > 1) Assert.Equal(1, r.HoppingComponents);
    }

    // The turn to Z is a proper rotation of each site's letters: it sends the jump to Z and keeps the
    // product rule X·Y = i·Z, which an improper map (a partial transpose) would break by a sign.
    [Theory]
    [InlineData('X')]
    [InlineData('Y')]
    [InlineData('Z')]
    public void The_Turn_Sends_The_Jump_To_Z_And_Keeps_XY_Equal_To_iZ(char jump)
    {
        Assert.Equal(('Z', 1), ComplementConnectionWitness.Turn(jump, jump));
        var (x, sx) = ComplementConnectionWitness.Turn(jump, 'X');
        var (y, sy) = ComplementConnectionWitness.Turn(jump, 'Y');
        var (z, sz) = ComplementConnectionWitness.Turn(jump, 'Z');
        var (product, phase) = PauliMask.Multiply(PauliMask.Parse(x.ToString()), PauliMask.Parse(y.ToString()));
        Assert.Equal(PauliMask.Parse(z.ToString()), product);
        // (sx·x)(sy·y) = sx·sy·i^phase·z must equal i·(sz·z): phase 1 with sx·sy = sz, or phase 3 with
        // sx·sy = −sz (the Hadamard: Z·Y = −i·X, and Y turns to −Y)
        Assert.True(phase is 1 or 3, $"phase {phase}");
        Assert.Equal(sz, phase == 1 ? sx * sy : -sx * sy);
    }

    [Fact]
    public void An_Undephased_Site_Is_Outside_And_Refused()
    {
        Assert.Throws<ArgumentException>(() => new ComplementConnectionWitness(3, "Z.Z"));
        Assert.Throws<ArgumentOutOfRangeException>(() => new ComplementConnectionWitness(ComplementConnectionWitness.MaxN + 1));
    }

    [Fact]
    public void The_Claim_Stands_On_F158_And_Reads_Theorem_Three()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<PalindromeComplementConnectionClaim>();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
        Assert.Same(registry.Get<PalindromeTwoEndCountClaim>(), claim.TwoEndCount);
        Assert.Contains("PROOF_PALINDROME_COMPLEMENT_CONNECTION.md", claim.Anchor);
        Assert.True(PalindromeComplementConnectionClaim.MixedAxesPalindrome("ZXZ", "YY"));
        Assert.False(PalindromeComplementConnectionClaim.MixedAxesPalindrome("ZXZ", "XY"));
        Assert.False(PalindromeComplementConnectionClaim.MixedAxesPalindrome("XYZ", ""));
        Assert.Throws<ArgumentException>(() => PalindromeComplementConnectionClaim.MixedAxesPalindrome("ZZ", ""));
    }

    [Fact]
    public void The_Children_Carry_The_Meetings()
    {
        var w = new ComplementConnectionWitness(4, "ZXZX", "YY..", "ring");
        var text = string.Join("\n", w.Children.Select(c => c.Summary));
        Assert.DoesNotContain("DO NOT MEET", text);
        Assert.Contains("They meet.", text);
    }
}
