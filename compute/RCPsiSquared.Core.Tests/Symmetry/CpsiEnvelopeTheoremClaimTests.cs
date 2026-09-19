using System;
using System.Linq;
using System.Reflection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class CpsiEnvelopeTheoremClaimTests
{
    [Fact]
    public void Claim_TierIsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, CpsiEnvelopeTheoremClaim.Shared.Tier);
    }

    [Fact]
    public void Claim_StatementCarriesTheNegativeBoundaryAndRetainedPositiveContent()
    {
        var name = CpsiEnvelopeTheoremClaim.Shared.Name;
        Assert.Contains("historical Envelope package is false", name);
        Assert.Contains("pointwise", name);
        Assert.Contains("absorbing", name);
        Assert.Contains("local-control", name);
        Assert.Contains("autonomous N=2", name);
        Assert.Contains("remains unproved", name);
        Assert.Contains("conditional convergence implication", name);
        Assert.Contains("Finite N/Q/K atlas rows are numerical readings only", name);
        Assert.DoesNotContain("proven Tier-1 for N=2", name);
    }

    [Fact]
    public void Claim_PublicDirectParentProperties_AreExactlyF25AndQuarter()
    {
        PropertyInfo[] directParents = typeof(CpsiEnvelopeTheoremClaim)
            .GetProperties(BindingFlags.Instance | BindingFlags.Public | BindingFlags.DeclaredOnly)
            .Where(property => typeof(Claim).IsAssignableFrom(property.PropertyType))
            .ToArray();

        Assert.Equal(2, directParents.Length);

        var actual = directParents
            .Select(property => (Name: property.Name, PropertyType: property.PropertyType))
            .OrderBy(parent => parent.Name, StringComparer.Ordinal)
            .ToArray();
        var expected = new[]
        {
            (Name: "F25", PropertyType: typeof(F25CPsiBellPlusPi2Inheritance)),
            (Name: "Quarter", PropertyType: typeof(QuarterAsBilinearMaxvalClaim)),
        };

        Assert.Equal(expected, actual);
    }

    [Fact]
    public void Build_SharesOneQuarterInstance_AcrossF25AndTheDirectEdge()
    {
        var c = CpsiEnvelopeTheoremClaim.Shared;
        // Build() threads ONE algebraic QuarterAsBilinearMaxvalClaim into both F25 and the direct edge.
        Assert.Same(c.Quarter, c.F25.Quarter);
    }
}
