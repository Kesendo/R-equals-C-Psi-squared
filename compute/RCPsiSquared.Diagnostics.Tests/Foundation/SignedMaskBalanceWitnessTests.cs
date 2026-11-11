using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live battery of F173: the claim's GF(2) readings agree with Theorem 2 of the complement connection,
/// walked on the bitstrings (Theorem D also against the balance, to N = 32), on every family the witness reads, and the
/// five controls (the pieces of Theorem B swapped, the power dropped from Theorem C, Theorem D with the parity of N or
/// with i₀ at the greatest valuation) disagree, so the battery can fail.</summary>
public class SignedMaskBalanceWitnessTests
{
    [Fact]
    public void EveryFamilyAgreesWithTheoremTwo_AndTheControlsDisagree()
    {
        var w = new SignedMaskBalanceWitness();
        Assert.All(w.Families, f => Assert.True(f.Agrees, $"{f.Name}: {f.FirstMismatch}"));
        Assert.True(w.AllAgree);
        Assert.True(w.ControlsBreak, $"controls {w.SwappedWidthsMismatches}, {w.SwappedQuotientsMismatches}, {w.DroppedRingPowerMismatches}, " +
                                     $"{w.ParityOnlyMismatches}, {w.MostValuedMismatches}");
        Assert.Equal(new[] { (864, 3024), (368, 1036), (19, 32), (49, 290), (1, 8), (1694, 9140), (1046, 28232), (416, 10418), (20, 100), (14, 14) },
            w.Families.Select(f => (f.Cells, f.Readings)).ToArray());                                // measured
        Assert.Equal(new Dictionary<int, int> { [5] = 6, [6] = 8 }, w.F115FloorsAtK4);
    }

    // The readings themselves, by hand: the page's worked examples.
    [Fact]
    public void TheClosedFormsReadThePagesExamples()
    {
        Assert.Equal((true, (int?)6), SignedMaskBalanceClaim.Chain("IXZX", "XIZX"));
        Assert.Equal((true, (int?)5), SignedMaskBalanceClaim.Chain("XIY", "XXX"));
        Assert.Equal((true, (int?)4), SignedMaskBalanceClaim.Chain("XIY", "XXZ"));
        Assert.Equal((false, (int?)null), SignedMaskBalanceClaim.Chain("XIX", "XXX"));
        Assert.True(SignedMaskBalanceClaim.RingBreaks("XY", 5));
        Assert.False(SignedMaskBalanceClaim.RingBreaks("XY", 6));
        Assert.True(SignedMaskBalanceClaim.RingBreaks("XIY", 6));
        Assert.False(SignedMaskBalanceClaim.RingBreaks("XIY", 8));
        Assert.True(SignedMaskBalanceClaim.RingBreaks("XIIIY", 12));
        Assert.False(SignedMaskBalanceClaim.RingBreaks("XIIIY", 16));
        Assert.False(SignedMaskBalanceClaim.RingBreaks("XYX", 7));
        Assert.False(SignedMaskBalanceClaim.RingBreaks("YY", 5));
        Assert.True(SignedMaskBalanceClaim.RingBreaks(new[] { "XXZ", "XZY" }, 7));     // v = 1, 2 at T = 1: the odd sign
        Assert.False(SignedMaskBalanceClaim.RingBreaks(new[] { "XXZ", "XZY" }, 6));    // T = 2: s_j + s_i₀·0 = 0
        Assert.True(SignedMaskBalanceClaim.RingBreaks(new[] { "XX", "XIY" }, 16));     // the larger valuation's sign is odd
        Assert.True(SignedMaskBalanceClaim.RingBreaks(new[] { "XX", "XIIY" }, 8));     // v = 1, 1: s₁ ≠ s₂
        Assert.False(SignedMaskBalanceClaim.RingBreaks(new[] { "XY", "XIIY" }, 8));    // v = 1, 1, both odd
        Assert.True(SignedMaskBalanceClaim.RingBreaks(new[] { "XY", "XIIY" }, 7));     // T = 1: every mask even, a sign odd
        Assert.True(SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnChain(new[] { "XXI", "XIX" }, 8)));
        Assert.False(SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnChain(new[] { "XXI", "XIY" }, 4)));
        Assert.Throws<System.ArgumentOutOfRangeException>(() => SignedMaskBalanceClaim.OnRing(new[] { "XIY" }, 4).ToList());
        Assert.Throws<System.ArgumentException>(() => SignedMaskBalanceClaim.OnChain(new[] { "XX", "XY" }, 3).ToList());
        Assert.Throws<System.ArgumentException>(() => SignedMaskBalanceClaim.RingBreaks(System.Array.Empty<string>(), 8));
        Assert.Throws<System.ArgumentOutOfRangeException>(() => SignedMaskBalanceClaim.OnRing(new[] { "XY" }, 65).ToList());
    }
}
