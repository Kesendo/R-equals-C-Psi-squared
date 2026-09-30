using System.Numerics;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using D = RCPsiSquared.Diagnostics.Foundation.PalindromeStringSpanWitness.Decision;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>From-below pins for <see cref="PalindromeStringSpanWitness"/>, F158's string route
/// (proof <c>docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md</c>, second corollary to Lemma 1 and
/// section (f5)).
///
/// <para>Every count is judged against <see cref="PalindromeTwoEndCountWitness"/>, which ranks the
/// Liouvillian itself and never forms a span, so each meeting is two constructions agreeing. The
/// seven N = 3 rows are that witness's own test rows. The decision each row takes is derived, not
/// read off: a Heisenberg bond forces one colour per component, which the fields must accept; the
/// Z-field row's word is Tr(H Z_0) = 2^3 · 3 at fields 3; a global pi rotation about an axis
/// orthogonal to a jump's letter fixes the Heisenberg chain and flips the jump, so every one-jump
/// word vanishes on the XYZ row and the far bound 0 decides; on (.X., X.X) H^1 A has no X_1 to pick up
/// and H^2 A picks up 2 · (10 · 3) from each bond, 120, times 2^3.</para>
///
/// <para>The two rank rows come out of a scan of all 4032 N = 3 chain rows of this family (numpy
/// ranks of the integer ad matrices and fw.odd_word_obstruction, 2026-09-30, and an independent review
/// scan): twelve rows carry no colouring, fire no one-jump word up to H^4 and have a far bound not
/// below the near lower bound, so only the ranks decide. Six are broken at (2, 1), such as (.X., YZY),
/// where the dense witness's characteristic polynomial separates the two sides at a point, a proof over
/// Q(i); six pair at (1, 1), such as (.X., Y.Z), where the polynomial identity holds at the sampled
/// points. These rows carry fields of EQUAL magnitude, which the colouring page's census (distinct
/// magnitudes) did not cover. The review scan took exact rational ranks and found the same twelve;
/// on the broken six no one-jump word fires up to H^12 either, nor a three-jump word up to power 4,
/// so on them the rank is the only reading at every budget tried.</para>
///
/// <para>One property is stated rather than gated: keeping the smaller nullity over TWO primes.
/// On every row tried here the two primes agree, and the Hamiltonian's coefficients are fixed (bonds
/// 10, fields 3), so no row was found that reduces badly at one prime, and a witness reading only the
/// first prime passes these tests, as MirrorWorld's CollisionGap records for its own second prime.</para></summary>
public class PalindromeStringSpanWitnessTests
{
    [Theory]
    [InlineData("ZZZ", "...", 4, 4, D.PalindromeByColour)]
    [InlineData("ZZZ", "Z..", 4, 0, D.BrokenByWord)]
    [InlineData("ZZZ", "X.X", 1, 1, D.PalindromeByColour)]
    [InlineData("XZZ", "..Y", 1, 1, D.PalindromeByColour)]
    [InlineData("Z.Z", "X.X", 1, 1, D.PalindromeByColour)]
    [InlineData("XYZ", "...", 1, 0, D.BrokenByCount)]
    [InlineData(".X.", "X.X", 6, 2, D.BrokenByWord)]
    [InlineData(".X.", "YZY", 2, 1, D.BrokenByRank)]
    [InlineData(".X.", "Y.Z", 1, 1, D.PalindromeByRank)]
    public void The_String_Route_Meets_The_Dense_Witness(string deph, string field, int near, int far, D decision)
    {
        var r = new PalindromeStringSpanWitness(3, deph, field).Read();
        var dense = new PalindromeTwoEndCountWitness(3, deph, field).Read();
        Assert.Equal((dense.NearCount, dense.FarCount), (r.NearUpper, r.FarUpper));
        Assert.Equal((near, far), (r.NearUpper, r.FarUpper));
        Assert.Equal(decision, r.Verdict);
        Assert.Equal(dense.PolynomialSaysPalindrome, PalindromeStringSpanWitness.IsPalindrome(r.Verdict));
        Assert.InRange(r.NearLower, 1, r.NearUpper);
        Assert.InRange(r.FarLower, 0, r.FarUpper);
        if (r.FarLower > 0) Assert.Equal(r.NearLower, r.FarLower);
    }

    [Fact]
    public void The_Word_Traces_Are_The_Derived_Integers()
    {
        var z = new PalindromeStringSpanWitness(3, "ZZZ", "Z..").Read();
        Assert.Equal("H^1·ZII", z.Word);
        Assert.Equal(new GaussianInteger(24, 0), z.WordTrace);
        var x = new PalindromeStringSpanWitness(3, ".X.", "X.X").Read();
        Assert.Equal("H^2·IXI", x.Word);
        Assert.Equal(new GaussianInteger(960, 0), x.WordTrace);
        Assert.Null(new PalindromeStringSpanWitness(3, "XYZ", "...").Read().Word);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(6)]
    [InlineData(8)]
    [InlineData(10)]
    public void The_Canonical_Chain_Holds_N_Plus_One_At_Both_Ends_Past_The_Dense_Guard(int n)
    {
        var r = new PalindromeStringSpanWitness(n).Read();
        int expected = PalindromeTwoEndCountClaim.CanonicalChainCount(n);
        Assert.Equal((expected, expected), (r.NearUpper, r.FarUpper));
        Assert.Equal(1 << n, r.DarkStrings);
        Assert.Equal(1 << n, r.LitStrings);
        Assert.Equal(new[] { new string('X', n), new string('Y', n) }, r.Colourings.OrderBy(s => s));
        Assert.Equal(D.PalindromeByColour, r.Verdict);
    }

    [Theory]
    [InlineData("ring")]
    [InlineData("complete")]
    public void Topology_Is_Carried_And_Meets_The_Dense_Witness(string topology)
    {
        foreach (var (deph, field) in new[] { ("ZZZZ", "...."), ("ZZZZ", "X..."), (".X.Z", "Y.Y.") })
        {
            var r = new PalindromeStringSpanWitness(4, deph, field, topology).Read();
            var dense = new PalindromeTwoEndCountWitness(4, deph, field, topology).Read();
            Assert.Equal((dense.NearCount, dense.FarCount), (r.NearUpper, r.FarUpper));
            Assert.Equal(dense.PolynomialSaysPalindrome, PalindromeStringSpanWitness.IsPalindrome(r.Verdict));
        }
    }

    [Fact]
    public void Exactness_Is_Labelled_Per_Decision()
    {
        Assert.True(PalindromeStringSpanWitness.IsExact(D.PalindromeByColour));
        Assert.True(PalindromeStringSpanWitness.IsExact(D.BrokenByWord));
        Assert.True(PalindromeStringSpanWitness.IsExact(D.BrokenByCount));
        Assert.False(PalindromeStringSpanWitness.IsExact(D.PalindromeByRank));
        Assert.False(PalindromeStringSpanWitness.IsExact(D.BrokenByRank));
        Assert.Contains("a rank reading", new PalindromeStringSpanWitness(3, ".X.", "YZY").Summary, StringComparison.Ordinal);
        Assert.Contains("an exact reading", new PalindromeStringSpanWitness(3).Summary, StringComparison.Ordinal);
    }

    [Fact]
    public void The_Guards_Are_Guards()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new PalindromeStringSpanWitness(1));
        Assert.Throws<ArgumentOutOfRangeException>(() => new PalindromeStringSpanWitness(PalindromeStringSpanWitness.MaxN + 1));
        Assert.Throws<ArgumentException>(() => new PalindromeStringSpanWitness(3, deph: "..."));
        Assert.Throws<ArgumentException>(() => new PalindromeStringSpanWitness(3, deph: "ZQZ"));
        Assert.Throws<ArgumentException>(() => new PalindromeStringSpanWitness(3, deph: "ZZ"));
        Assert.Throws<ArgumentException>(() => new PalindromeStringSpanWitness(3, topology: "star"));
        // one jump on twelve sites leaves 2^23 dark strings: refused, not enumerated
        Assert.Throws<InvalidOperationException>(() =>
            new PalindromeStringSpanWitness(12, deph: "Z" + new string('.', 11)).Read());
    }

    [Fact]
    public void The_Children_Recompute_And_The_Dense_Meeting_Is_Printed()
    {
        var w = new PalindromeStringSpanWitness(3, ".X.", "YZY");
        var children = w.Children.ToList();
        Assert.Equal(7, children.Count);
        Assert.All(children, c => Assert.False(string.IsNullOrWhiteSpace(c.Summary)));
        Assert.Contains("They meet.", children[4].Summary, StringComparison.Ordinal);
        Assert.Contains("Not run", new PalindromeStringSpanWitness(6).Children.ElementAt(4).Summary, StringComparison.Ordinal);
        // the falsifier row beside the canonical one, at a size the dense witness cannot reach
        var beside = new PalindromeStringSpanWitness(6).Children.ElementAt(5).Summary;
        Assert.Contains("near 7, far 7", beside, StringComparison.Ordinal);
        Assert.Contains("near 7, far 0, BrokenByWord by the word H^1·ZIIIII", beside, StringComparison.Ordinal);
    }
}
