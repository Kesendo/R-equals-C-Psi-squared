using System;
using System.Linq;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86;

/// <summary>F86 Stage E1: top-level c=2 universal-shape synthesis Claim. Verifies
/// integration into <see cref="F86KnowledgeBase"/>, the empirical witness collection
/// surfacing, the directional Endpoint &gt; Interior split, and the c=2 guard.
///
/// <para>The Tier honesty is load-bearing: <see cref="C2UniversalShapeDerivation.IsClosedFormDerived"/>
/// must be false in the current session (Tier1Candidate, not Tier1Derived). Future sessions
/// land the closed form via cross-block perturbation, projector-overlap lift, or char-poly
/// factorisation; promotion flips both flags together.</para>
/// </summary>
public class C2UniversalShapeDerivationTests
{
    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    public void KnowledgeBase_C2_HasUniversalShapeDerivation(int N)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var kb = new F86KnowledgeBase(block);
        var claim = kb.C2UniversalShape;
        Assert.NotNull(claim);
        Assert.Equal(Tier.Tier1Candidate, claim!.Tier);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    public void C2UniversalShapeDerivation_ExposesEmpiricalWitnesses(int N)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var derivation = C2UniversalShapeDerivation.Build(block);
        Assert.Equal(block.NumBonds, derivation.Witnesses.Count);
        // Directional structure pinned: Endpoint > Interior at the class-mean level.
        Assert.True(derivation.EndpointMean > derivation.InteriorMean,
            $"At N={N}: EndpointMean={derivation.EndpointMean:F4} should exceed " +
            $"InteriorMean={derivation.InteriorMean:F4}");
        // Closed form not yet derived this session — Tier1Candidate, not Tier1Derived.
        Assert.False(derivation.IsClosedFormDerived);
        Assert.NotNull(derivation.PendingDerivationNote);
    }

    [Fact]
    public void C2UniversalShapeDerivation_DirectionalGap_IsEmpiricallyConsistent()
    {
        // Empirical anchor: the Endpoint − Interior gap is ≈ 0.022 ± 0.003 across N=5..8.
        var block = new CoherenceBlock(N: 7, n: 1, gammaZero: 0.05);
        var derivation = C2UniversalShapeDerivation.Build(block);
        Assert.InRange(derivation.DirectionalGap, 0.015, 0.030);
    }

    [Fact]
    public void C2UniversalShapeDerivation_ThrowsIfNotC2()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05); // c=3
        Assert.Throws<ArgumentException>(() => C2UniversalShapeDerivation.Build(block));
    }

    [Fact]
    public void F86OpenQuestions_Item1Prime_ReflectsC2Closure()
    {
        var item1Prime = F86OpenQuestions.Standard
            .FirstOrDefault(q => q.Name.Contains("Item 1") || q.Description.Contains("Item 1"));
        Assert.NotNull(item1Prime);
        // The updated description should mention c=2 closure status.
        Assert.Contains("c=2", item1Prime!.Description);
        Assert.Contains("Tier1Candidate", item1Prime.Description, StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public void F86KnowledgeBase_ForC1Block_LeavesC2UniversalShapeNull()
    {
        // c=1 has only one HD channel — no inter-channel SVD, no c=2 universal-shape claim.
        var block = new CoherenceBlock(N: 3, n: 0, gammaZero: 0.05);
        var kb = new F86KnowledgeBase(block);
        Assert.Equal(1, block.C);
        Assert.Null(kb.C2UniversalShape);
    }

    [Fact]
    public void F86KnowledgeBase_ForC3Block_LeavesC2UniversalShapeNull()
    {
        // c=3 falls under Item 4' (multi-k extension), not the c=2 derivation plan.
        var block = new CoherenceBlock(N: 7, n: 2, gammaZero: 0.05);
        var kb = new F86KnowledgeBase(block);
        Assert.Equal(3, block.C);
        Assert.Null(kb.C2UniversalShape);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    public void C2UniversalShapeDerivation_HwhmRatioComposition_MirrorsStageD2Outcome(int N)
    {
        // The C2UniversalShapeDerivation must wrap the Stage D2 C2HwhmRatio: its means
        // and witnesses must agree with the underlying primitive.
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var derivation = C2UniversalShapeDerivation.Build(block);

        Assert.Equal(
            derivation.HwhmRatio.HwhmLeftOverQPeakMean(BondClass.Endpoint),
            derivation.EndpointMean,
            12);
        Assert.Equal(
            derivation.HwhmRatio.HwhmLeftOverQPeakMean(BondClass.Interior),
            derivation.InteriorMean,
            12);
        Assert.Equal(
            derivation.HwhmRatio.IsAnalyticallyDerived,
            derivation.IsClosedFormDerived);
        Assert.Same(derivation.HwhmRatio.Witnesses, derivation.Witnesses);
    }
}
