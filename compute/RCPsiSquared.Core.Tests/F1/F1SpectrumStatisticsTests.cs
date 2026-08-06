using System.Numerics;
using RCPsiSquared.Core.F1;
using Xunit;

namespace RCPsiSquared.Core.Tests.F1;

/// <summary>The F1 symmetry distance (<see cref="F1SpectrumStatistics.MaxF1PairingDistance"/>): it
/// must test MULTISET-closure under λ ↦ −2σ − λ, not set-closure. A set/Hausdorff distance is blind
/// to a dropped/duplicated eigenvalue that still has a same-valued neighbour — the exact failure a
/// spectrum-reconstruction bug produces. These pin the multiplicity-aware behaviour at its Core home.</summary>
public class F1SpectrumStatisticsTests
{
    [Fact]
    public void MaxF1PairingDistance_IsZero_ForAnF1SymmetricMultiset()
    {
        // {0, -2, -1+i, -1-i} is closed under λ ↦ -2σ - λ with σ=1 (center -1).
        var spectrum = new[]
        {
            new Complex(0, 0), new Complex(-2, 0),
            new Complex(-1, 1), new Complex(-1, -1),
        };
        Assert.True(F1SpectrumStatistics.MaxF1PairingDistance(spectrum, sigma: 1.0) < 1e-12);
    }

    [Fact]
    public void MaxF1PairingDistance_DetectsAMultiplicityDefect()
    {
        // {+1 x3, -1 x1} is closed under λ ↦ -λ (σ=0) as a SET but NOT as a multiset
        // (reflected = {-1 x3, +1 x1}). A set distance reports ~0; the greedy NN with removal must
        // see the surplus +1 forced onto a -1 at distance 2.
        var defective = new[]
        {
            new Complex(1, 0), new Complex(1, 0), new Complex(1, 0), new Complex(-1, 0),
        };
        Assert.Equal(2.0, F1SpectrumStatistics.MaxF1PairingDistance(defective, sigma: 0.0), 9);
    }

    [Fact]
    public void MaxF1PairingDistance_RejectsAnEmptySpectrum_RatherThanScoringItPerfect()
    {
        // Until 2026-08-06 this returned 0.0, the best possible reading, for a spectrum that was
        // never scored: the matching loop never runs and the seed survives. It is pinned here
        // rather than in the rate sibling alone, because the sibling's guard was written first
        // and its stated reason was true of this method too.
        Assert.Throws<ArgumentException>(
            () => F1SpectrumStatistics.MaxF1PairingDistance(System.Array.Empty<Complex>(), sigma: 1.0));
    }

    // ----- the rate-space sibling: d ↦ 2σ − d via the embedding λ = −d -----
    //
    // The multiplicity behaviour itself is pinned above and is NOT re-pinned here; what these
    // add is the part the eigenvalue tests cannot see: that the rate entry point carries the
    // factor of two as d ↦ 2σ − d, that it routes through the same matcher rather than a
    // private one, and that the embedding is exact rather than close.
    //
    // Mind the centre, which this file got wrong once: the reflection d ↦ 2σ − d has its FIXED
    // POINT at σ, so the palindrome is centred on σ and it is the PARTNERS that sum to 2σ.
    //
    // The first four tests were checked against mutants: flipping the embedding sign kills four
    // of them, disabling the taken[] removal kills three. The last two pin POLICY rather than
    // arithmetic (a blind spot and an input contract), so there is no wrong implementation for
    // them to fall over; they exist to stop a future edit from quietly reversing a decision.

    [Fact]
    public void MaxF1RatePairingDistance_IsZero_ForRatesPalindromicAboutSigma()
    {
        // Centred on σ = 0.5: {0.2, 0.8} and {0.4, 0.6} each sum to 2σ = 1.0 and so pair.
        // Against the off-by-two implementation, d ↦ σ − d, this reads 1.1 rather than 0.
        var rates = new[] { 0.2, 0.4, 0.6, 0.8 };
        Assert.True(F1SpectrumStatistics.MaxF1RatePairingDistance(rates, sigma: 0.5) < 1e-15);
    }

    [Fact]
    public void MaxF1RatePairingDistance_UsesTwoSigmaNotSigma()
    {
        // The discriminating pair for the factor of two, stated as a positive value rather than
        // as "not zero". The reflection is d ↦ 2σ − d, so partners SUM to 2σ while the fixed
        // point sits at σ. {0.25, 0.75} sums to 1, so it pairs exactly at σ = 0.5. Feed it
        // σ = 0.25 instead, the value an off-by-two caller would pass, and the reflection becomes
        // d ↦ 0.5 − d: 0.25 maps to itself and consumes its own partner, leaving 0.75 to match
        // the surviving −0.25 at distance 1.0.
        var rates = new[] { 0.25, 0.75 };
        Assert.True(F1SpectrumStatistics.MaxF1RatePairingDistance(rates, sigma: 0.5) < 1e-15);
        Assert.Equal(1.0, F1SpectrumStatistics.MaxF1RatePairingDistance(rates, sigma: 0.25), 12);
    }

    [Fact]
    public void MaxF1RatePairingDistance_DetectsAMultiplicityDefect_InRateSpace()
    {
        // The rate-space image of the eigenvalue test above: {0.5 ×3, 1.5 ×1} at σ = 1 is
        // set-closed under d ↦ 2 − d and multiset-broken. A matcher without removal reports ~0;
        // with removal the surplus 0.5 is forced onto a 1.5 at distance 1.
        var defective = new[] { 0.5, 0.5, 0.5, 1.5 };
        Assert.Equal(1.0, F1SpectrumStatistics.MaxF1RatePairingDistance(defective, sigma: 1.0), 12);
    }

    [Fact]
    public void MaxF1RatePairingDistance_EqualsTheEmbeddedEigenvalueRoute_Exactly()
    {
        // The identity the rate entry point exists to express. Deliberately an EXACT comparison:
        // d ↦ −d is an isometry and preserves the index order the greedy walk consumes, so there
        // is no rounding to tolerate. A tolerance here would hide a changed matcher.
        var rates = new[] { 0.13, 0.29, 0.31, 0.47, 0.53, 0.69, 0.71, 0.87 };
        const double sigma = 0.5;

        var embedded = new Complex[rates.Length];
        for (int i = 0; i < rates.Length; i++) embedded[i] = new Complex(-rates[i], 0.0);

        Assert.Equal(
            F1SpectrumStatistics.MaxF1PairingDistance(embedded, sigma),
            F1SpectrumStatistics.MaxF1RatePairingDistance(rates, sigma));
    }

    [Fact]
    public void MaxF1RatePairingDistance_IsBlind_ToABreakThatMovesOnlyAFrequency()
    {
        // The projection's blind spot, pinned as a fact rather than left as a caveat, because a
        // measured ratio between the two readings does NOT establish it: on an N=6 chain at
        // γ=0.5 the two agree to 1.0x, and reading that as "the projection is as good" is the
        // error this test exists to forbid.
        //
        // {−0.1±i, −0.3±2i} at σ=0.2 is NOT F1-closed: −0.1+i reflects to −0.3−i, which is not
        // in the multiset. But its rates {0.1, 0.1, 0.3, 0.3} pair exactly, centred on σ = 0.2
        // with partners summing to 2σ = 0.4, because the projection has thrown away the
        // frequencies that carry the break.
        var spectrum = new[]
        {
            new Complex(-0.1, 1.0), new Complex(-0.1, -1.0),
            new Complex(-0.3, 2.0), new Complex(-0.3, -2.0),
        };
        var rates = new[] { 0.1, 0.1, 0.3, 0.3 };
        const double sigma = 0.2;

        Assert.True(F1SpectrumStatistics.MaxF1RatePairingDistance(rates, sigma) < 1e-15,
            "the rate projection must read this spectrum as clean, which is the point");
        Assert.True(F1SpectrumStatistics.MaxF1PairingDistance(spectrum, sigma) > 0.1,
            "the full complex check must see the break the projection cannot");
    }

    [Fact]
    public void MaxF1RatePairingDistance_RejectsAnEmptyMultiset_RatherThanScoringItPerfect()
    {
        // An empty input has no distance. Returning 0.0 would publish the BEST possible reading
        // for a spectrum that was never scored, which is how a metric fails quiet.
        Assert.Throws<ArgumentException>(
            () => F1SpectrumStatistics.MaxF1RatePairingDistance(System.Array.Empty<double>(), sigma: 1.0));
        Assert.Throws<ArgumentNullException>(
            () => F1SpectrumStatistics.MaxF1RatePairingDistance(null!, sigma: 1.0));
    }
}
