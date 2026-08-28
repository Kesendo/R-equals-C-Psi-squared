using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>From-below pins for <see cref="PalindromeTwoEndCountWitness"/>, the F158 live lab
/// (proof <c>docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md</c>, gate
/// <c>simulations/f138_rank_criterion.py</c>).
///
/// <para>The witness is an independent construction and not a port: it rebuilds the Liouvillian in
/// C# over GF(p), so agreeing with the committed Python run is two computations meeting rather than
/// one repeated. The expected values below are read from that run and from the proof file, never
/// recomputed from the criterion the witness is being checked against.</para>
///
/// <para><b>Two-sided throughout, and deliberately so.</b> Every row where the criterion says
/// PALINDROME sits beside one where it says BROKEN, because a witness that had lost the ability to
/// report a break would pass a one-sided suite unchanged. The sharpest pair is the canonical chain
/// against the SAME dephasing with one Z field added: the near count does not move and the far one
/// empties, which is the criterion doing the only thing it can do.</para></summary>
public class PalindromeTwoEndCountWitnessTests
{
    [Fact]
    public void CanonicalChain_BothCountsAreNPlusOne_AndThePalindromeHolds()
    {
        foreach (int n in new[] { 2, 3, 4 })
        {
            var r = new PalindromeTwoEndCountWitness(n).Read();
            Assert.Equal(n + 1, r.NearCount);
            Assert.Equal(n + 1, r.FarCount);
            Assert.Equal(PalindromeTwoEndCountClaim.CanonicalChainCount(n), r.NearCount);
            Assert.True(r.CriterionSaysPalindrome);
            Assert.True(r.PolynomialSaysPalindrome);
        }
    }

    [Fact]
    public void OneZFieldOnTheDephasingAxis_EmptiesTheFarEndAndBreaksThePalindrome()
    {
        // The same dephasing as the canonical row: the near count must NOT move, and the far one
        // must go to zero. A witness that had broken both ends together would pass a test that
        // only checked the verdict.
        var r = new PalindromeTwoEndCountWitness(3, deph: "ZZZ", field: "Z..").Read();
        Assert.Equal(4, r.NearCount);
        Assert.Equal(0, r.FarCount);
        Assert.False(r.CriterionSaysPalindrome);
        Assert.False(r.PolynomialSaysPalindrome);
    }

    [Fact]
    public void TheTwoRoutesToTheSameTwoDimensionsMeet_OnEveryRowTried()
    {
        // Lemma 1 live: the nullities of L and L + 2 sigma against the nullities of the operator
        // conditions, computed by a route that never forms L at all.
        var rows = new (string Deph, string Field)[]
        {
            ("ZZZ", "..."), ("ZZZ", "Z.."), ("ZZZ", "X.X"), ("XZZ", "..Y"),
            ("Z.Z", "X.X"), ("XYZ", "..."), (".X.", "X.X"),
        };
        foreach (var (deph, field) in rows)
        {
            var r = new PalindromeTwoEndCountWitness(3, deph, field).Read();
            Assert.Equal(r.NearCount, r.NearFromOperatorConditions);
            Assert.Equal(r.FarCount, r.FarFromOperatorConditions);
        }
    }

    [Fact]
    public void TheCriterionAndThePolynomialAgree_OnRowsOfBothVerdicts()
    {
        var rows = new (string Deph, string Field)[]
        {
            ("ZZZ", "..."), ("ZZZ", "Z.."), ("ZZZ", "X.X"), ("ZZZ", "XY."),
            ("XZZ", "..Y"), ("Z.Z", "X.X"), ("XYZ", "..."), (".X.", "X.X"),
            ("ZZZ", "XXZ"), ("Z..", "XYZ"),
        };
        int holds = 0, breaks = 0;
        foreach (var (deph, field) in rows)
        {
            var r = new PalindromeTwoEndCountWitness(3, deph, field).Read();
            Assert.Equal(r.PolynomialSaysPalindrome, r.CriterionSaysPalindrome);
            if (r.PolynomialSaysPalindrome) holds++; else breaks++;
        }
        // Both verdicts must appear, or the agreement above is one-sided and proves nothing.
        Assert.True(holds > 0, "no palindromic row in the suite: the agreement would be one-sided");
        Assert.True(breaks > 0, "no broken row in the suite: the agreement would be one-sided");
    }

    [Fact]
    public void TheKernelDominates_OnEveryRowTried()
    {
        // Lemma 3's inequality, which holds with no palindrome anywhere in sight.
        foreach (var deph in new[] { "ZZZ", "XZZ", "XYZ", "Z.Z", ".X." })
            foreach (var field in new[] { "...", "Z..", "X.X", "XY." })
            {
                var r = new PalindromeTwoEndCountWitness(3, deph, field).Read();
                Assert.True(r.FarCount <= r.NearCount,
                    $"dim ker(L + 2 sigma) = {r.FarCount} exceeded dim ker L = {r.NearCount} " +
                    $"at deph {deph}, field {field}");
            }
    }

    [Fact]
    public void TopologyIsCarried_AndTheRingIsNotTheChain()
    {
        var chain = new PalindromeTwoEndCountWitness(4, topology: "chain").Read();
        var ring = new PalindromeTwoEndCountWitness(4, topology: "ring").Read();
        Assert.Equal(5, chain.NearCount);
        Assert.Equal(5, ring.NearCount);
        Assert.True(chain.CriterionSaysPalindrome);
        Assert.True(ring.CriterionSaysPalindrome);
    }

    [Fact]
    public void TheGuardsAreGuards_AndSayWhatTheyAre()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new PalindromeTwoEndCountWitness(5));
        Assert.Throws<ArgumentOutOfRangeException>(() => new PalindromeTwoEndCountWitness(1));
        Assert.Throws<ArgumentException>(() => new PalindromeTwoEndCountWitness(3, deph: "ZZ"));
        Assert.Throws<ArgumentException>(() => new PalindromeTwoEndCountWitness(3, deph: "ZQZ"));
        Assert.Throws<ArgumentException>(() => new PalindromeTwoEndCountWitness(3, topology: "star"));
    }

    [Fact]
    public void ChildrenRecompute_AndTheSummaryCarriesBothVerdicts()
    {
        var holds = new PalindromeTwoEndCountWitness(3);
        var breaks = new PalindromeTwoEndCountWitness(3, deph: "ZZZ", field: "Z..");
        Assert.Contains("PALINDROME", holds.Summary, StringComparison.Ordinal);
        Assert.Contains("BROKEN", breaks.Summary, StringComparison.Ordinal);
        Assert.Equal(5, holds.Children.Count());
        Assert.All(holds.Children, c => Assert.False(string.IsNullOrWhiteSpace(c.Summary)));
    }
}
