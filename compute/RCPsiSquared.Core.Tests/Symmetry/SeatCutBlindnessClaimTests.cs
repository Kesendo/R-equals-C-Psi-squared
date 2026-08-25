using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>From-below pins for <see cref="SeatCutBlindnessClaim"/> (F157, proof
/// <c>docs/proofs/PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md</c>). The closed forms are pinned
/// against values printed by the committed run
/// <c>simulations/results/seat_cut_blindness/seat_cut_blindness_run.txt</c> and by the run of the
/// proof's own gate script, never against a second evaluation of the same expression.
///
/// <para>Two-sided throughout: every nonzero blind count sits beside a zero the same expression
/// has to produce, because a law that only ever returns positive numbers would not distinguish a
/// working gcd from one that returns its first argument.</para></summary>
public class SeatCutBlindnessClaimTests
{
    private static SeatCutBlindnessClaim BuildClaim() =>
        new(new F4KernelDimensionByComponentsClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void Anchor_NamesTheProofFileAndBothSourcePages()
    {
        var anchor = BuildClaim().Anchor;
        Assert.Contains("docs/proofs/PROOF_BLIND_SEAT_SPAN_AND_NODE_LEMMA.md", anchor);
        Assert.Contains("experiments/THE_SEAT_THAT_CUTS.md", anchor);
        Assert.Contains("experiments/THE_BLIND_SITE.md", anchor);
        Assert.Contains("docs/ANALYTICAL_FORMULAS.md", anchor);
    }

    [Fact]
    public void Parent_IsHeldAsATypedEdge()
    {
        var parent = new F4KernelDimensionByComponentsClaim();
        Assert.Same(parent, new SeatCutBlindnessClaim(parent).KernelByComponents);
        Assert.Throws<ArgumentNullException>(() => new SeatCutBlindnessClaim(null!));
    }

    // The Heisenberg law at N = 7 and N = 11: the run file's "kernel" part prints dim ker = 1 + blind
    // at every seat, so blind = dim ker - 1. N = 11 centre 5 is the value the arc opened on.
    [Theory]
    [InlineData(7, 0, 0)]
    [InlineData(7, 1, 0)]
    [InlineData(7, 2, 0)]
    [InlineData(7, 3, 3)]
    [InlineData(7, 4, 0)]
    [InlineData(7, 5, 0)]
    [InlineData(7, 6, 0)]
    [InlineData(11, 5, 5)]
    [InlineData(9, 4, 4)]
    public void BlindHeisenberg_MatchesTheCommittedRun(int n, int seat, int expected)
    {
        Assert.Equal(expected, SeatCutBlindnessClaim.BlindHeisenberg(n, seat));
    }

    // The XY law at N = 7: the run's "xy" part prints Heisenberg 3:3 against XY 1:1 3:3 5:1, the
    // disagreement that made the two books two objects.
    [Theory]
    [InlineData(7, 0, 0)]
    [InlineData(7, 1, 1)]
    [InlineData(7, 2, 0)]
    [InlineData(7, 3, 3)]
    [InlineData(7, 5, 1)]
    [InlineData(11, 7, 3)]
    [InlineData(11, 8, 2)]
    [InlineData(11, 9, 1)]
    public void BlindXy_MatchesTheCommittedRun(int n, int seat, int expected)
    {
        Assert.Equal(expected, SeatCutBlindnessClaim.BlindXy(n, seat));
    }

    [Fact]
    public void TheTwoBooksDisagreeLoudly_AndAgreeAtTheFixedCentre()
    {
        // N = 12: Heisenberg blinds four seats, XY none.
        var heisenberg = Enumerable.Range(0, 12)
            .Where(s => SeatCutBlindnessClaim.BlindHeisenberg(12, s) > 0).ToArray();
        var xy = Enumerable.Range(0, 12)
            .Where(s => SeatCutBlindnessClaim.BlindXy(12, s) > 0).ToArray();
        Assert.Equal(4, heisenberg.Length);
        Assert.Empty(xy);

        // At the reflection-fixed centre of an odd chain the two laws agree at (N-1)/2.
        foreach (int n in new[] { 5, 7, 9, 11, 13 })
        {
            int centre = (n - 1) / 2;
            Assert.Equal(centre, SeatCutBlindnessClaim.BlindHeisenberg(n, centre));
            Assert.Equal(centre, SeatCutBlindnessClaim.BlindXy(n, centre));
        }
    }

    [Fact]
    public void TheEndSeatsAreNeverBlind_TheLemmaJ1Zero()
    {
        // Lemma J1: no eigenvector of an unreduced Jacobi matrix vanishes at an end, so an end seat
        // sees every eigenspace. This is the zero that the positive counts above sit beside.
        for (int n = 2; n <= 30; n++)
        {
            Assert.Equal(0, SeatCutBlindnessClaim.BlindHeisenberg(n, 0));
            Assert.Equal(0, SeatCutBlindnessClaim.BlindHeisenberg(n, n - 1));
            Assert.Equal(0, SeatCutBlindnessClaim.BlindXy(n, 0));
            Assert.Equal(0, SeatCutBlindnessClaim.BlindXy(n, n - 1));
        }
    }

    [Fact]
    public void APrimeChainHidesTheHeisenbergLaw_WhichIsWhyItStayedInvisible()
    {
        // Against a reflection-parity reading the Heisenberg law coincides at every seat of a PRIME
        // chain, and the committed data lived at N = 5 and N = 11, both prime: the only blind seat
        // of a prime chain is its centre. The first discriminating case is N = 6.
        foreach (int n in new[] { 5, 7, 11, 13 })
        {
            var blind = Enumerable.Range(0, n)
                .Where(s => SeatCutBlindnessClaim.BlindHeisenberg(n, s) > 0).ToArray();
            Assert.Equal(new[] { (n - 1) / 2 }, blind);
        }
        // N = 6 is the first case that discriminates the two readings, so it is pinned exactly
        // rather than as "at least one": a law that blinded every seat would pass the weaker form.
        Assert.Equal(new[] { 0, 1, 0, 0, 1, 0 },
            Enumerable.Range(0, 6).Select(s => SeatCutBlindnessClaim.BlindHeisenberg(6, s)).ToArray());
    }

    [Fact]
    public void TheInspectNodes_IllustrateTheSentencesTheyCarry()
    {
        // The claim's ExtraChildren interpolate live values into prose. Nothing else exercises them,
        // and an example that renders as "0 against 0" would silently illustrate the opposite of its
        // own sentence, which is what happened before this test existed.
        var nodes = BuildClaim().Children.ToList();
        var books = nodes.Single(c => c.DisplayName.Contains("two uniform chain laws"));
        Assert.Contains("N = 12: Heisenberg blinds seats 1, 4, 7 and 10, XY none", books.Summary);
        Assert.Contains("seat 4: 1 against 0", books.Summary);
        Assert.Equal(4, Enumerable.Range(0, 12).Count(s => SeatCutBlindnessClaim.BlindHeisenberg(12, s) > 0));

        var parity = nodes.Single(c => c.DisplayName.Contains("parity-forced"));
        Assert.Contains("parity-forced = True", parity.Summary);
        Assert.Contains("(False)", parity.Summary);
    }

    [Fact]
    public void ParityForcing_IsXyOnly_AndOddSeatOfOddChainOnly()
    {
        Assert.True(SeatCutBlindnessClaim.ParityForcedXy(9, 3, SeatCutBook.Xy));
        Assert.True(SeatCutBlindnessClaim.ParityForcedXy(9, 1, SeatCutBook.Xy));
        // The falsifiers: the same call on the other book, on an even seat, and on an even chain.
        Assert.False(SeatCutBlindnessClaim.ParityForcedXy(9, 3, SeatCutBook.Heisenberg));
        Assert.False(SeatCutBlindnessClaim.ParityForcedXy(9, 2, SeatCutBook.Xy));
        Assert.False(SeatCutBlindnessClaim.ParityForcedXy(8, 3, SeatCutBook.Xy));

        // And where parity forces, the count is genuinely positive, which is the content.
        foreach (int n in new[] { 5, 7, 9, 11 })
            foreach (int seat in Enumerable.Range(0, n).Where(s => s % 2 == 1))
            {
                Assert.True(SeatCutBlindnessClaim.ParityForcedXy(n, seat, SeatCutBook.Xy));
                Assert.True(SeatCutBlindnessClaim.BlindXy(n, seat) > 0);
            }
    }

    [Fact]
    public void SpanOnTheZeroFreeUniformChain_IsOnePlusBlind()
    {
        // Corollary B, written as literals rather than as a second evaluation of the same gcd: a
        // gate whose expectation is computed from the thing under test cannot fail.
        // N = 7 Heisenberg, seats 0..6, from the run file's "kernel" table:
        Assert.Equal(new[] { 1, 1, 1, 4, 1, 1, 1 },
            Enumerable.Range(0, 7)
                .Select(s => SeatCutBlindnessClaim.SpanOnZeroFreeUniformChain(7, s, SeatCutBook.Heisenberg))
                .ToArray());
        // N = 7 XY, where the identity 1 + (gcd(j+1, N+1) - 1) collapses to the bare gcd(j+1, 8):
        Assert.Equal(new[] { 1, 2, 1, 4, 1, 2, 1 },
            Enumerable.Range(0, 7)
                .Select(s => SeatCutBlindnessClaim.SpanOnZeroFreeUniformChain(7, s, SeatCutBook.Xy))
                .ToArray());
    }

    [Theory]
    [InlineData(1, 0)]
    [InlineData(0, 0)]
    [InlineData(-3, 0)]
    public void TooFewSites_Throws(int n, int seat)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => SeatCutBlindnessClaim.BlindHeisenberg(n, seat));
        Assert.Throws<ArgumentOutOfRangeException>(() => SeatCutBlindnessClaim.BlindXy(n, seat));
    }

    [Theory]
    [InlineData(7, -1)]
    [InlineData(7, 7)]
    [InlineData(7, 99)]
    public void SeatOutOfRange_Throws(int n, int seat)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => SeatCutBlindnessClaim.BlindHeisenberg(n, seat));
        Assert.Throws<ArgumentOutOfRangeException>(() => SeatCutBlindnessClaim.BlindXy(n, seat));
        Assert.Throws<ArgumentOutOfRangeException>(() => SeatCutBlindnessClaim.ParityForcedXy(n, seat, SeatCutBook.Xy));
    }
}
