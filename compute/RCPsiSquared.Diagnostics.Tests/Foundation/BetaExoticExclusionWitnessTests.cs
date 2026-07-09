using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live β-exotic witness (<c>inspect --root betaexotic</c>): it must recompute the
/// certified disc-layer reading at inspect time and report the exclusion only when the certificate
/// certifies. The number it reads was an unvalidated diagnostic before 2026-07-09; these tests pin
/// that the witness never reports an exclusion from an uncertified reading.</summary>
public class BetaExoticExclusionWitnessTests
{
    [Fact]
    [Trait("Category", "FOLDRESULTANT")]
    public void Witness_ExcludesTheBetaExotic_AtBothParities()
    {
        var w = new BetaExoticExclusionWitness();

        Assert.Equal(5, w.N);
        // a beta-exotic (Puiseux exponent 3/2) forces a disc zero of order 3; we must be strictly below it
        Assert.True(w.MaxDiscMultiplicityREven < BetaExoticExclusionWitness.BetaExoticDiscOrder);
        Assert.True(w.MaxDiscMultiplicityROdd < BetaExoticExclusionWitness.BetaExoticDiscOrder);
        Assert.True(w.BetaExoticExcluded);
        Assert.Contains("EXCLUDED", w.Summary);
    }

    /// <summary>The exclusion must be CONJUNCTIVE over both parity sectors and both certification
    /// flags: the physics review found the four N=5 forced seeds split 2 R-even + 2 R-odd, so a
    /// single-sector reading would be silently partial.</summary>
    [Fact]
    [Trait("Category", "FOLDRESULTANT")]
    public void Witness_ReportsPerParity_AndNamesTheOpenScope()
    {
        var w = new BetaExoticExclusionWitness();
        var children = w.Children.ToList();

        Assert.Contains(children, c => c.DisplayName.Contains("R-odd"));
        Assert.Contains(children, c => c.DisplayName.Contains("R-even"));
        // the scope boundary must be stated where a reader will meet it
        Assert.Contains(children, c => c.Summary.Contains("s₆"));
        Assert.Contains("Per-N", w.Summary);
    }

    [Fact]
    [Trait("Category", "FOLDRESULTANT")]
    public void Witness_RefusesEveryUncertifiedN()
    {
        // the witness must never certify an N the certificate has not been run at. N = 9 is out of reach
        // by this route (a 324-dimensional block); N = 3 was never in scope.
        Assert.Throws<ArgumentOutOfRangeException>(() => new BetaExoticExclusionWitness(9));
        Assert.Throws<ArgumentOutOfRangeException>(() => new BetaExoticExclusionWitness(3));
        Assert.Equal(new[] { 5, 7 }, BetaExoticExclusionWitness.CertifiedN);
    }
}
