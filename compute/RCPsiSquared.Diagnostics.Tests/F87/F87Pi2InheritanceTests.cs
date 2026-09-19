using System.Linq;
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
        var f1 = new F1Pi2Inheritance(new RCPsiSquared.Core.F1.F1PalindromeIdentity(), ladder, memoryLoop);
        return new F87Pi2Inheritance(f1);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
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
    public void Constructor_RejectsNullF1Inheritance()
    {
        Assert.Throws<ArgumentNullException>(() => new F87Pi2Inheritance(null!));
    }

    [Fact]
    public void Anchor_References_F87_AndF1Pi2()
    {
        // F87's only typed parent edges are F87TrichotomyClassification (the F87
        // closed form itself) and F1Pi2Inheritance (the transitive "2" = a_0 anchor).
        // No KleinFour edge: F87 entstand vor Klein cells, sie sind unterschiedlich.
        var f = Build();
        Assert.Contains("F87TrichotomyClassification.cs", f.Anchor);
        Assert.Contains("F1Pi2Inheritance.cs", f.Anchor);
        Assert.DoesNotContain("Pi2KnowledgeBaseClaims.cs", f.Anchor);
    }

    [Fact]
    public void LiveChildren_KeepTheExactF87TrichotomyDefinition()
    {
        var definition = Assert.Single(Build().Children, c => c.DisplayName == "F87 closed form");

        Assert.Equal(
            "trichotomy via F1 residual: truly iff ‖M‖ < ε; soft iff M ≠ 0 but spectrum pairs; hard iff no pairing",
            definition.Summary);
    }

    [Fact]
    public void LiveChildren_ExposeOnlyTheFiniteMarrakeshHardwareRow()
    {
        var children = Build().Children.ToArray();
        var marrakesh = Assert.Single(children, c => c.DisplayName == "Marrakesh finite Δ row");
        string rendered = string.Join("\n", children.Select(c => $"{c.DisplayName}\n{c.Summary}"));

        Assert.Equal(
            "hardware job d7mjnjjaq2pc73a1pk4g (2026-04-26): Δ(soft − truly) = −0.722",
            marrakesh.Summary);
        Assert.DoesNotContain("F87 hardware confirmation", children.Select(c => c.DisplayName));
        Assert.DoesNotContain("Kingston", rendered);
        Assert.DoesNotContain("regime-uniformity", rendered);
    }
}
