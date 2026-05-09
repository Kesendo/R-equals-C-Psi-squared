using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F87;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

public class F87Pi2InheritanceTests
{
    private static F87Pi2Inheritance Build()
    {
        var ladder = new Pi2DyadicLadderClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f1 = new F1Pi2Inheritance(ladder, memoryLoop);
        return new F87Pi2Inheritance(f1);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Fact]
    public void KleinClassCount_IsFour()
    {
        Assert.Equal(4, F87Pi2Inheritance.KleinClassCount);
    }

    [Fact]
    public void TransitivelyInheritedTwoFactor_IsExactlyTwo()
    {
        // F87's discriminator M is F1's residual; F1's "2" coefficient is a_0 = 2 on
        // the Pi2 dyadic ladder. F87 inherits this number-level constant transitively
        // through F1Pi2Inheritance.
        Assert.Equal(2.0, Build().TransitivelyInheritedTwoFactor, precision: 14);
    }

    [Fact]
    public void KleinAlignedClassNames_HasFourEntries_MatchingKleinFourCells()
    {
        // F87's 4-way Π²-class refinement IS the KleinFour decomposition; the names
        // align Pp / Pm / Mp / Mm = truly / Π²-even non-truly / Π²-odd A / Π²-odd B.
        var f = Build();
        Assert.Equal(4, f.KleinAlignedClassNames.Count);
        Assert.Contains("truly", f.KleinAlignedClassNames);
        Assert.Contains("pi2_even_nontruly", f.KleinAlignedClassNames);
        Assert.Contains("pi2_odd_subgroup_A", f.KleinAlignedClassNames);
        Assert.Contains("pi2_odd_subgroup_B", f.KleinAlignedClassNames);
    }

    [Fact]
    public void Constructor_RejectsNullF1Inheritance()
    {
        Assert.Throws<ArgumentNullException>(() => new F87Pi2Inheritance(null!));
    }

    [Fact]
    public void Anchor_References_F87_AndF1Pi2_AndKleinFour()
    {
        var f = Build();
        Assert.Contains("F87TrichotomyClassification.cs", f.Anchor);
        Assert.Contains("F1Pi2Inheritance.cs", f.Anchor);
        Assert.Contains("Pi2KnowledgeBaseClaims.cs", f.Anchor);
    }
}
