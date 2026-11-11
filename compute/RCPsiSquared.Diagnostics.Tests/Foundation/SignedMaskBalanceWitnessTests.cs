using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live battery of F173: the claim's GF(2) readings agree with Theorem 2 of the complement connection,
/// walked on the bitstrings (Theorem D also against the balance, to N = 32), on every family the witness reads, and the
/// twelve controls (the pieces of Theorem B swapped, the power dropped from Theorem C, Theorem D with the parity of N or
/// with i₀ at the greatest valuation, Theorem E's conditions each dropped, the tests of Corollaries E2 and E3 each
/// dropped) disagree, so the battery can fail.</summary>
public class SignedMaskBalanceWitnessTests
{
    [Fact]
    public void EveryFamilyAgreesWithTheoremTwo_AndTheControlsDisagree()
    {
        var w = new SignedMaskBalanceWitness();
        Assert.All(w.Families, f => Assert.True(f.Agrees, $"{f.Name}: {f.FirstMismatch}"));
        Assert.True(w.AllAgree);
        Assert.True(w.ControlsBreak, $"controls {w.SwappedWidthsMismatches}, {w.SwappedQuotientsMismatches}, {w.DroppedRingPowerMismatches}, " +
                                     $"{w.ParityOnlyMismatches}, {w.MostValuedMismatches}, {string.Join(",", w.DroppedConditionMismatches)}, {w.DroppedReadSetMismatches}, {w.DroppedOrderingMismatches}, {w.DroppedReadSumMismatches}, {w.DroppedE3OrderingMismatches}");
        Assert.Equal(new[] { (864, 3024), (368, 1036), (19, 32), (49, 290), (1, 8), (1694, 9140), (1046, 28232), (416, 10418), (20, 100), (428, 3936), (16, 128), (352, 2176), (352, 10560), (128, 768), (928, 1856), (32, 224), (21, 21) },
            w.Families.Select(f => (f.Cells, f.Readings)).ToArray());                                // measured
        Assert.Equal(new Dictionary<int, int> { [5] = 6, [6] = 8 }, w.F115FloorsAtK4);
        Assert.Equal(new Dictionary<SignedMaskBalanceClaim.ChainShapeBranch, int>                     // measured, quoted on the page
            {
                [SignedMaskBalanceClaim.ChainShapeBranch.SquaresFail] = 928, [SignedMaskBalanceClaim.ChainShapeBranch.NoRelation] = 192,
                [SignedMaskBalanceClaim.ChainShapeBranch.ReadSumNonzero] = 928, [SignedMaskBalanceClaim.ChainShapeBranch.HolonomyPlus] = 400,
                [SignedMaskBalanceClaim.ChainShapeBranch.HolonomyMinus] = 400,
            }, w.E3Branches);
        Assert.Equal((464, 576), (w.DroppedReadSumMismatches, w.DroppedE3OrderingMismatches));
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

    // Theorem E by hand: the read sets of two-string flip sets and the branches of Corollary E2.
    [Fact]
    public void TheCharacterConnectionReadsThePagesExamples()
    {
        var mid = SignedMaskBalanceClaim.ReadRatio(0b111, SignedMaskBalanceClaim.PlaceStrings(new[] { ("XXX", 1L), ("XYX", 1L) }, 3, ring: false));
        Assert.Equal(SignedMaskBalanceClaim.RatioKind.Character, mid.Kind);
        Assert.Equal(0b010UL, mid.Read);                       // X·(X + Y)·X reads its middle site
        Assert.Equal(1, mid.Quarter % 2);                      // ω = ±i, the middle site inside the mask
        var unequal = SignedMaskBalanceClaim.ReadRatio(0b111, SignedMaskBalanceClaim.PlaceStrings(new[] { ("XXX", 1L), ("XYX", 2L) }, 3, ring: false));
        Assert.Equal(SignedMaskBalanceClaim.RatioKind.NotCharacter, unequal.Kind);   // c₁² ≠ c₂²
        var flipFlop = SignedMaskBalanceClaim.ReadRatio(0b11, SignedMaskBalanceClaim.PlaceStrings(new[] { ("XX", 1L), ("YY", 1L) }, 2, ring: false));
        Assert.Equal(SignedMaskBalanceClaim.RatioKind.AbsentEdge, flipFlop.Kind);    // XX + YY vanishes where the bits agree
        Assert.Equal(SignedMaskBalanceClaim.RingShapeBranch.ReadSetNotDivisible, SignedMaskBalanceClaim.RingOneShape(new[] { ("XXX", 1L), ("XYX", 1L) }, 9).Branch);
        Assert.Equal(SignedMaskBalanceClaim.RingShapeBranch.NoRelation, SignedMaskBalanceClaim.RingOneShape(new[] { ("XXX", 1L), ("XYX", 1L) }, 8).Branch);
        Assert.Equal(SignedMaskBalanceClaim.RingShapeBranch.OnePlusXDoesNotDivide, SignedMaskBalanceClaim.RingOneShape(new[] { ("XXX", 1L), ("YYY", 1L) }, 6).Branch);
        Assert.Equal(SignedMaskBalanceClaim.RingShapeBranch.HolonomyMinus, SignedMaskBalanceClaim.RingOneShape(new[] { ("XIY", 1L), ("XZX", 1L) }, 6).Branch);
        Assert.Equal(SignedMaskBalanceClaim.RingShapeBranch.HolonomyPlus, SignedMaskBalanceClaim.RingOneShape(new[] { ("XIY", 1L), ("XZX", 1L) }, 8).Branch);
        var product = new[] { ("XX", 1L), ("XY", 1L), ("YX", -1L), ("YY", -1L) };          // (X − Y)⊗(X + Y): hol(g) = (−1)^(N − 2)
        Assert.True(SignedMaskBalanceClaim.RingOneShape(product, 7).Breaks);
        Assert.False(SignedMaskBalanceClaim.RingOneShape(product, 8).Breaks);
        Assert.False(SignedMaskBalanceClaim.RingOneShape(product, 7, orderingTerm: false).Breaks);
        Assert.Equal(SignedMaskBalanceClaim.CharacterVerdict.Breaks,
            SignedMaskBalanceClaim.DecideCharacter(SignedMaskBalanceClaim.PlaceStrings(new[] { ("XZ", 1L), ("ZY", 1L) }, 6, ring: false)));   // the chain's end
        Assert.Equal(SignedMaskBalanceClaim.CharacterVerdict.Outside,
            SignedMaskBalanceClaim.DecideCharacter(SignedMaskBalanceClaim.PlaceStrings(new[] { ("XX", 1L), ("YY", 1L) }, 4, ring: false)));
    }

    // Corollary E3 by hand: two flip shapes on the chain, the read sum D and hol(u₀).
    [Fact]
    public void TheTwoShapeChainReadsThePagesExamples()
    {
        var field = new[] { ("X", 1L), ("Y", 1L) };
        var plus = new[] { ("XXX", 1L), ("YYY", 1L) };
        var minus = new[] { ("XXX", 1L), ("YYY", -1L) };
        Assert.Equal((3, 0UL), SignedMaskBalanceClaim.ChainTwoShapeData(field, plus));                  // q₂R₁ + q₁R₂ = f·1 + 1·f
        Assert.Equal((false, SignedMaskBalanceClaim.ChainShapeBranch.HolonomyPlus), SignedMaskBalanceClaim.ChainTwoShapes(field, plus, 5));
        Assert.Equal((true, SignedMaskBalanceClaim.ChainShapeBranch.HolonomyMinus), SignedMaskBalanceClaim.ChainTwoShapes(field, minus, 5));
        Assert.True(SignedMaskBalanceClaim.ChainTwoShapes(field, plus, 5, orderingTerm: false).Breaks);   // the ordering term is −1 here
        var middle = new[] { ("XXX", 1L), ("XYX", 1L) };
        var wide = new[] { ("XIIX", 1L) };
        Assert.Equal((4, 0b110UL), SignedMaskBalanceClaim.ChainTwoShapeData(middle, wide));              // D = (1 + x)·x
        Assert.Equal(SignedMaskBalanceClaim.ChainShapeBranch.ReadSumNonzero, SignedMaskBalanceClaim.ChainTwoShapes(middle, wide, 4).Branch);
        Assert.Equal(SignedMaskBalanceClaim.ChainShapeBranch.SquaresFail, SignedMaskBalanceClaim.ChainTwoShapes(middle, wide, 5).Branch);
        Assert.Equal((false, SignedMaskBalanceClaim.ChainShapeBranch.NoRelation), SignedMaskBalanceClaim.ChainTwoShapes(new[] { ("XIY", 1L) }, new[] { ("XXX", 1L) }, 4));
        Assert.Equal((true, SignedMaskBalanceClaim.ChainShapeBranch.HolonomyMinus), SignedMaskBalanceClaim.ChainTwoShapes(new[] { ("XIY", 1L) }, new[] { ("XXX", 1L) }, 5));   // Theorem B's c = 1 from N* = 5
        var fieldBesideProduct = SignedMaskBalanceClaim.ChainTwoShapes(field, new[] { ("XX", 1L), ("XY", 1L), ("YX", -1L), ("YY", -1L) }, 4);
        Assert.Equal((true, SignedMaskBalanceClaim.ChainShapeBranch.HolonomyMinus), fieldBesideProduct);  // ω_f²·ω_p = −1, the ordering term +1
        Assert.Equal(fieldBesideProduct, SignedMaskBalanceClaim.ChainTwoShapes(field, new[] { ("XX", 1L), ("XY", 1L), ("YX", -1L), ("YY", -1L) }, 4, orderingTerm: false));
        Assert.Equal((true, SignedMaskBalanceClaim.ChainShapeBranch.ModuliDiffer),
            SignedMaskBalanceClaim.ChainTwoShapes(new[] { ("IX", 1L), ("ZX", 1L) }, new[] { ("XXX", 1L) }, 3));   // 1 + (−1)^x₀ vanishes on one side only
        Assert.Equal((false, SignedMaskBalanceClaim.ChainShapeBranch.Outside),
            SignedMaskBalanceClaim.ChainTwoShapes(new[] { ("XX", 1L), ("YY", 1L) }, new[] { ("X", 1L) }, 3));     // an absent edge: undecided
        Assert.Throws<System.ArgumentException>(() => SignedMaskBalanceClaim.ChainTwoShapeData(new[] { ("XX", 1L), ("YY", 1L) }, new[] { ("X", 1L) }));
        Assert.Throws<System.ArgumentException>(() => SignedMaskBalanceClaim.ChainTwoShapes(new[] { ("XX", 1L) }, new[] { ("XY", 1L) }, 4));
        Assert.Throws<System.ArgumentException>(() => SignedMaskBalanceClaim.ChainTwoShapes(new[] { ("XX", 1L), ("XIY", 1L) }, wide, 5));
    }
}
