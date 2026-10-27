using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>From-below pins for <see cref="HeisenbergColourLawWitness"/> and <see cref="HeisenbergColourLawClaim"/>
/// (F169, proof <c>docs/proofs/PROOF_HEISENBERG_COLOURS_OFF_THE_LETTERS.md</c>).</summary>
public class HeisenbergColourLawWitnessTests
{
    // Every default row with a prediction: the far kernel from the commutator in the string basis (fields included, no
    // split) equals the law's prediction, proven rows and measured ones alike; the predictions are not all equal (0, 1
    // and N + 1 occur). The law predicts no dimension under a field along the common axis, only a broken palindrome.
    [Fact]
    public void The_Commutator_Meets_The_Law_On_Every_Row()
    {
        var readings = new HeisenbergColourLawWitness().Read();
        Assert.All(readings, r => Assert.True(r.Agrees, $"{r.Row.Name}: law {r.Predicted}, commutator {r.Measured}"));
        Assert.Equal(new int?[] { 0, 1, 5 }, readings.Where(r => r.Compared).Select(r => r.Predicted).Distinct().OrderBy(x => x));
        // a uniform field along the common axis at even N: the palindrome breaks, yet the far kernel keeps one
        // (non-invertible) element, which is why the law predicts no dimension there
        var along = readings.Single(r => !r.Compared);
        Assert.Equal(1, along.Measured);
        Assert.Contains(readings, r => r.Proven);
        Assert.Contains(readings, r => !r.Proven);
    }

    // The witness can disagree: off the theorem's hypotheses (two disconnected bonds, generic axes) the far kernel is
    // the product of the two bonds' own, one each, while the law, written for a connected graph, would say none.
    [Fact]
    public void Off_A_Connected_Graph_The_Commutator_Leaves_The_Law()
    {
        var row = new HeisenbergColourLawWitness.Row("two bonds apart", 4, new[] { (0, 1), (2, 3) }, new[] { 1L, 2 },
            new (long, long, long)[] { (1, 2, 0), (0, 1, 2), (2, 0, 1), (1, -1, 2) }, new (long, long, long)[4]);
        var r = HeisenbergColourLawWitness.ReadRow(row);
        Assert.False(r.Proven);
        Assert.Equal(0, r.Predicted);
        Assert.Equal(1, r.Measured);
    }

    [Fact]
    public void The_Claim_Stands_On_F158_And_F165_And_Reads_Its_Class()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<HeisenbergColourLawClaim>();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
        Assert.Same(registry.Get<PalindromeTwoEndCountClaim>(), claim.TwoEndCount);
        Assert.Same(registry.Get<PalindromeComplementConnectionClaim>(), claim.ComplementConnection);
        Assert.Contains("PROOF_HEISENBERG_COLOURS_OFF_THE_LETTERS.md", claim.Anchor);

        var ring = new[] { (0, 1), (1, 2), (2, 3), (3, 0) };
        var k4 = new[] { (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3) };
        var generic = new (long, long, long)[] { (1, 2, 0), (0, 1, 2), (2, 0, 1), (1, -1, 2) };
        Assert.Equal(1, HeisenbergColourLawClaim.ReducedCoreCycles(4, ring, generic));
        Assert.Equal(3, HeisenbergColourLawClaim.ReducedCoreCycles(4, k4, generic));
        Assert.Equal(1, HeisenbergColourLawClaim.ReducedCoreCycles(4, k4, new[] { generic[0], generic[0], generic[1], generic[1] }));
        Assert.True(HeisenbergColourLawClaim.Proven(4, ring, generic));
        Assert.False(HeisenbergColourLawClaim.Proven(4, k4, generic));
        Assert.True(HeisenbergColourLawClaim.Proven(4, k4, new (long, long, long)[] { (1, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0) }));
        Assert.Equal(1, HeisenbergColourLawClaim.PredictedFarDimension(new (long, long, long)[] { (2, 2, 1), (1, 2, -1), (1, 0, 2) }, new (long, long, long)[3]));
        Assert.Equal(0, HeisenbergColourLawClaim.PredictedFarDimension(generic, new (long, long, long)[4]));
        Assert.Equal(4, HeisenbergColourLawClaim.PredictedFarDimension(Enumerable.Repeat((0L, 0L, 1L), 3).ToArray(), new (long, long, long)[3]));
        Assert.Null(HeisenbergColourLawClaim.PredictedFarDimension(Enumerable.Repeat((0L, 0L, 1L), 3).ToArray(), new (long, long, long)[] { (0, 0, 1), (0, 0, 0), (0, 0, 0) }));
    }
}
