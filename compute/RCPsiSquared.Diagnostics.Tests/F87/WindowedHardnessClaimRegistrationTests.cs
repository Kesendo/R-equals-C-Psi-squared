using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>Wiring audit for <see cref="WindowedHardnessClaim"/>: it surfaces F115 (the windowed
/// F87 hardness combinatorial theory) as a registered, registry-queryable Claim. F115 is the
/// mask-combinatorial reading of the same diagonal-cell hard/soft line that
/// <see cref="F87DiagonalCellBipartiteWitnessSet"/> reads through the chiral-K bipartite criterion,
/// so it is wired as that set's child. Built from the full
/// <see cref="KnowledgeRegistryFactory.BuildDefault"/> registry, so the §7.5/§7.6 converse parent
/// (Tier1Derived) is present and the strength-inheritance check (5 ≥ 5) is exercised in
/// production.</summary>
public class WindowedHardnessClaimRegistrationTests
{
    [Fact]
    public void BuildDefault_ContainsWindowedHardnessClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<WindowedHardnessClaim>());
    }

    [Fact]
    public void WindowedHardnessClaim_TypedParent_IsPresentInRegistry()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<F87DiagonalCellBipartiteWitnessSet>(),
            "typed parent F87DiagonalCellBipartiteWitnessSet must be registered");
    }

    [Fact]
    public void WindowedHardnessClaim_Ancestors_ContainBipartiteWitnessSet()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();

        var ancestors = registry.AncestorsOf<WindowedHardnessClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F87DiagonalCellBipartiteWitnessSet), ancestors);
    }

    [Fact]
    public void WindowedHardnessClaim_TierIsTier1Derived()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.Equal(Tier.Tier1Derived, registry.Get<WindowedHardnessClaim>().Tier);
    }

    [Fact]
    public void WindowedHardnessClaim_SelfCheck_AllCasesPass()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<WindowedHardnessClaim>();

        Assert.True(claim.Cases.Count >= 3,
            "battery should hold the hard-differing-valuation, soft-equal-valuation, and size-law cases");
        foreach (var c in claim.Cases)
            Assert.True(c.Passes, $"battery case '{c.Name}': expected {c.Expected}, got {c.Actual}");
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }

    /// <summary>The criterion the Claim wraps is exactly <see cref="WindowedObstructionScan"/>'s
    /// one-number rule: hard iff the two X/Y masks have different (1+x)-adic valuations. A hard
    /// witness (the extremal family p1 = x + x^{k-1}, p2 = 1 + x^{k-1} at k=4) has differing
    /// valuations; a soft witness (two distinct masks of equal valuation) has equal valuations;
    /// the size-law helper saturates min(2W-1, 2k-3).</summary>
    [Fact]
    public void WindowedHardnessClaim_Battery_MatchesWindowedObstructionScan()
    {
        // Hard: differing (1+x)-valuations -> IsHardPair true (extremal body-bound family, k=4).
        ulong hard1 = (1UL << 1) | (1UL << 3); // x + x^3
        ulong hard2 = 1UL | (1UL << 3);        // 1 + x^3
        Assert.True(WindowedObstructionScan.IsHardPair(hard1, hard2));
        Assert.NotEqual(WindowedObstructionScan.ValuationAtOnePlusX(hard1),
                        WindowedObstructionScan.ValuationAtOnePlusX(hard2));

        // Soft: equal (1+x)-valuations -> IsHardPair false (two distinct v=1 masks).
        ulong soft1 = 0b11;   // 1 + x        (v = 1)
        ulong soft2 = 0b110;  // x + x^2      (v = 1)
        Assert.False(WindowedObstructionScan.IsHardPair(soft1, soft2));
        Assert.Equal(WindowedObstructionScan.ValuationAtOnePlusX(soft1),
                     WindowedObstructionScan.ValuationAtOnePlusX(soft2));

        // Size law: the extremal family saturates the body bound 2k-3 (k=4 -> 5).
        Assert.Equal(2 * 4 - 3, WindowedObstructionScan.GcdFormulaSize(hard1, hard2));

        // Size law (windowed regime): max minimal-odd-cycle == min(2W-1, 2k-3) at k=3, N=4.
        var scan = WindowedObstructionScan.Scan(k: 3, n: 4);
        int w = 4 - 3 + 1;
        Assert.Equal(System.Math.Min(2 * w - 1, 2 * 3 - 3), scan.MinOddCycleSizes.Keys.Max());
    }
}
