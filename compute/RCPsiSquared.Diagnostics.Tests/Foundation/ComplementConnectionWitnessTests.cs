using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>From-below pins for <see cref="ComplementConnectionWitness"/> and
/// <see cref="PalindromeComplementConnectionClaim"/> (proof
/// <c>docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md</c>). The graph route (Theorems 2 and 4) is
/// judged against <see cref="PalindromeStringSpanWitness"/>, which ranks commutators on string spans and
/// never forms a basis state, and against Theorems 1 and 3's rules, which carry no graph at all; Theorem 5's
/// rule, where no graph reads the row, against the string route's exact verdict. The five N = 3 rows with
/// every site dephased are the string witness's own test rows, with its counts.</summary>
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
        Assert.Equal((near, far, pal), (r.Near, r.Far, r.Palindrome));
        Assert.Equal((s.NearUpper, s.FarUpper), (r.Near, r.Far));
        Assert.Equal(pal, r.Predicted!.Value);
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
                Assert.Equal(1, r.Near);
                Assert.Equal(PalindromeComplementConnectionClaim.MixedAxesPalindrome(deph, field.Where(c => c != '.')),
                    r.Palindrome);
                if (!r.Palindrome) { mixedBroken++; Assert.Equal(0, r.Far); }
            }
            if (rows++ % 7 == 0)
            {
                var s = new PalindromeStringSpanWitness(3, deph, field, topology).Read();
                Assert.Equal((s.NearUpper, s.FarUpper), (r.Near, r.Far));
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
        Assert.Equal(pal, r.Predicted!.Value);
        Assert.Equal((s.NearUpper, s.FarUpper), (r.Near, r.Far));
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
        Assert.Equal((pal, pal), (r.Palindrome, r.Predicted!.Value));
        Assert.Equal((s.NearUpper, s.FarUpper), (r.Near, r.Far));
        if (deph.Distinct().Count() > 1) Assert.Equal(1, r.Near);
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
    public void A_Row_Without_A_Jump_And_An_Unknown_Bond_Are_Refused()
    {
        Assert.Throws<ArgumentException>(() => new ComplementConnectionWitness(3, "..."));
        Assert.Throws<ArgumentException>(() => new ComplementConnectionWitness(3, "ZZZ", bonds: "XZ"));
        Assert.Throws<ArgumentOutOfRangeException>(() => new ComplementConnectionWitness(ComplementConnectionWitness.MaxN + 1));
    }

    // Theorem 4 on its own graph: ZZ bonds, every pattern at N = 3 on the chain and the ring (and at N = 4 on
    // the chain, every fifth row) with at least one site dephased, each undephased site keeping Z (its field absent or Z). The
    // sector sums equal the string route's counts on every row, and the verdict its verdict wherever that one
    // is exact; enough rows carry a nonzero cross-sector term that dropping those terms would fail the counts.
    [Theory]
    [InlineData(3, "chain")]
    [InlineData(3, "ring")]
    [InlineData(4, "chain")]
    public void Undephased_Sites_That_Keep_A_Letter_Follow_Theorem_Four(int n, string topology)
    {
        int rows = 0, cross = 0, pal = 0, exact = 0, seen = 0;
        foreach (string deph in Patterns(n, ".XYZ"))
        foreach (string field in Patterns(n, ".XYZ"))
        {
            if (!deph.Contains('.') || deph.All(c => c == '.')) continue;
            if (Enumerable.Range(0, n).Any(s => deph[s] == '.' && field[s] is 'X' or 'Y')) continue;
            if (n == 4 && seen++ % 5 != 0) continue;
            var w = new ComplementConnectionWitness(n, deph, field, topology, "ZZ");
            Assert.True(w.GraphReads, $"{deph} {field}");
            var r = w.Read();
            var s = new PalindromeStringSpanWitness(n, deph, field, topology, bonds: "ZZ").Read();
            Assert.True((s.NearUpper, s.FarUpper) == (r.Near!.Value, r.Far!.Value),
                $"{deph} {field} {topology}: graph {r.Near}, {r.Far}; strings {s.NearUpper}, {s.FarUpper}");
            Assert.True(s.NearLower <= r.Near && s.FarLower <= r.Far, $"{deph} {field}: below a lower bound");
            if (PalindromeStringSpanWitness.IsExact(s.Verdict))
            {
                exact++;
                Assert.Equal(PalindromeStringSpanWitness.IsPalindrome(s.Verdict), r.Palindrome);
            }
            if (r.CrossNear > 0 || r.CrossFar > 0) cross++;
            if (r.Palindrome) pal++;
            rows++;
        }
        Assert.True(cross > 10 && pal > 10 && exact > 10, $"cross {cross}, palindromic {pal}, exact {exact}");
        if (n == 3) Assert.Equal(864 + 144, rows); // one undephased site (3·9·2·16) and two (3·3·4·4), as the census
    }

    // Theorem 5 beside the string route: Heisenberg bonds, exactly one undephased site, every axis and field
    // pattern at N = 3 on the chain and the ring. No graph reads these rows; the rule meets the string route's
    // verdict on every row where that verdict is exact, and every row is exact here.
    [Theory]
    [InlineData("chain")]
    [InlineData("ring")]
    public void One_Undephased_Site_Follows_Theorem_Five(string topology)
    {
        int rows = 0, pal = 0;
        foreach (string deph in Patterns(3, ".XYZ"))
        foreach (string field in Patterns(3, ".XYZ"))
        {
            if (deph.Count(c => c == '.') != 1) continue;
            var w = new ComplementConnectionWitness(3, deph, field, topology);
            Assert.False(w.GraphReads);
            var r = w.Read();
            Assert.StartsWith("Theorem 5", r.TheoremClass);
            Assert.True(r.Exact, $"{deph} {field} {topology}: a rank reading");
            Assert.True(r.Predicted == r.Palindrome, $"{deph} {field} {topology}: rule {r.Predicted}, strings {r.Palindrome}");
            if (r.Palindrome) pal++;
            rows++;
        }
        Assert.Equal(3 * 9 * 64, rows); // the undephased site, two axes, a field pattern
        Assert.Equal(279, pal);          // N·(3·2^(N−1)·2^N − 3), the count the MirrorWorld gate derives
    }

    // XX + YY bonds with every site dephased, every pattern at N = 3 on the chain: under uniform Z Theorem 1's
    // rule meets the graph, every other assignment has no closed rule; the graph meets the string route on
    // every seventh row either way. ZZ bonds with every site dephased have no rule either.
    [Fact]
    public void XY_Bonds_Read_Theorem_One_Under_Z_And_Theorem_Two_Alone_Otherwise()
    {
        int rows = 0, ruled = 0;
        foreach (string deph in Patterns(3, "XYZ"))
        foreach (string field in Patterns(3, ".XYZ"))
        {
            var r = new ComplementConnectionWitness(3, deph, field, "chain", "XY").Read();
            if (deph == "ZZZ")
            {
                ruled++;
                Assert.StartsWith("Theorem 1", r.TheoremClass);
                Assert.True(r.Predicted == r.Palindrome, $"{deph} {field}: {r}");
            }
            else
            {
                Assert.Null(r.Predicted);
                Assert.StartsWith("Theorem 2 alone", r.TheoremClass);
            }
            if (rows++ % 7 == 0)
            {
                var s = new PalindromeStringSpanWitness(3, deph, field, "chain", bonds: "XY").Read();
                Assert.Equal((s.NearUpper, s.FarUpper), (r.Near!.Value, r.Far!.Value));
            }
        }
        Assert.Equal(64, ruled);
        var zz = new ComplementConnectionWitness(3, "ZZZ", "X..", "chain", "ZZ").Read();
        Assert.Null(zz.Predicted);
        Assert.StartsWith("Theorem 2 alone", zz.TheoremClass);
    }

    // Past the string route's span (2^(N + undephased) strings) the graph still reads Theorem 4, and the
    // inspect says the string route was not run instead of failing.
    [Fact]
    public void A_Wide_Theorem_Four_Row_Is_Read_By_The_Graph_Alone()
    {
        var w = new ComplementConnectionWitness(12, "ZZZZZZZ.....", null, "chain", "ZZ");
        Assert.True(w.GraphReads);
        Assert.Equal(1 << 5, w.Read().Sectors);
        Assert.Contains(w.Children, c => c.Summary.StartsWith("Not run", StringComparison.Ordinal));
    }

    // Two undephased sites under Heisenberg bonds: no theorem gives a rule, and the witness says so.
    [Fact]
    public void Two_Undephased_Sites_Under_Heisenberg_Bonds_Have_No_Rule()
    {
        var r = new ComplementConnectionWitness(3, ".X.", "X.X").Read();
        Assert.Null(r.Predicted);
        Assert.Null(r.Near);
        Assert.StartsWith("outside", r.TheoremClass);
    }

    private static IEnumerable<string> Patterns(int n, string alphabet)
    {
        int k = alphabet.Length, total = 1;
        for (int i = 0; i < n; i++) total *= k;
        for (int code = 0; code < total; code++)
        {
            var c = new char[n];
            for (int i = 0, v = code; i < n; i++, v /= k) c[i] = alphabet[v % k];
            yield return new string(c);
        }
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
