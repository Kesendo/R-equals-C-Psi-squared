using System;
using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Wiring audit for <see cref="SidewaysSpinLadderClaim"/> (the <c>sideways_spin_ladder</c> open
/// arc, typed 2026-08-09): S⁺(ρ) = Σ_l (−1)^l c_l†ρc_l† intertwines L for Σ-odd real-symmetric pure
/// hopping and any γ_l, and the F125 fold family is the two S⁺ chain interiors, one spin-(N−3)/2 multiplet
/// per chain with measured Clebsch-Gordan S⁺ transport norms (the η chains inferred from the count).
/// One typed parent: <see cref="SpectatorIntertwinerClaim"/> (F125,
/// Tier1Derived); this claim Tier1Candidate (the multiplet half is measured at N=5,7,9, not derived).</summary>
public class SidewaysSpinLadderClaimRegistrationTests
{
    [Fact]
    public void BuildDefault_ContainsClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<SidewaysSpinLadderClaim>());
    }

    [Fact]
    public void Claim_IsTier1Candidate()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.Equal(Tier.Tier1Candidate, registry.Get<SidewaysSpinLadderClaim>().Tier);
    }

    [Fact]
    public void Claim_Ancestors_ContainTypedParent()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var ancestors = registry.AncestorsOf<SidewaysSpinLadderClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(SpectatorIntertwinerClaim), ancestors);
    }

    [Fact]
    public void ClosedForms_MatchTheMeasuredNumbers()
    {
        // ℓ = (N−3)/2: spin 1 at N=5 (F125's stated case), spin 2 at N=7.
        Assert.Equal(1.0, SidewaysSpinLadderClaim.ChainSpin(5), 15);
        Assert.Equal(2.0, SidewaysSpinLadderClaim.ChainSpin(7), 15);
        // The measured N=7 norms 2, √6, √6, 2 are the CG coefficients at ℓ = 2.
        Assert.Equal(2.0, SidewaysSpinLadderClaim.CgTransportNorm(2.0, -2.0), 15);
        Assert.Equal(Math.Sqrt(6.0), SidewaysSpinLadderClaim.CgTransportNorm(2.0, -1.0), 15);
        Assert.Equal(Math.Sqrt(6.0), SidewaysSpinLadderClaim.CgTransportNorm(2.0, 0.0), 15);
        Assert.Equal(2.0, SidewaysSpinLadderClaim.CgTransportNorm(2.0, 1.0), 15);
        // Odd-N orbit size 4N−8 = 4(N−2), the multiplet dimension: 12 at N=5, 20 at N=7, 28 at N=9.
        // Even N refuses (the proof gives 4N−12 there and the chain reading is untested).
        Assert.Equal(12, SidewaysSpinLadderClaim.OrbitSizeOddN(5));
        Assert.Equal(20, SidewaysSpinLadderClaim.OrbitSizeOddN(7));
        Assert.Equal(28, SidewaysSpinLadderClaim.OrbitSizeOddN(9));
        Assert.Throws<ArgumentOutOfRangeException>(() => SidewaysSpinLadderClaim.OrbitSizeOddN(6));
    }
}
