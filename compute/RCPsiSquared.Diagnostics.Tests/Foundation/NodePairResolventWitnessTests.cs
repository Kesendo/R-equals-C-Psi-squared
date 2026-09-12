using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class NodePairResolventWitnessTests
{
    [Fact]
    public void CanonicalReducedResolvent_HasNodePairZeroAndNonNodeControl()
    {
        var reading = new NodePairResolventWitness().Reading;

        Assert.Equal(BigRational.Zero, reading.NodePairEntry);
        Assert.Equal(new BigRational(-1, 8), reading.NodeNonNodeEntry);
    }

    [Fact]
    public void DirectCharacteristicPolynomial_MeetsTheFactorizedFormula()
    {
        var reading = new NodePairResolventWitness().Reading;

        Assert.True(reading.FactorizationMatches);
        Assert.NotEmpty(reading.DirectPolynomial.CoefficientsLowToHigh);
        Assert.Equal(reading.DirectPolynomial, reading.FactorizedPolynomial);
    }

    [Fact]
    public void UniformCentrePointwiseCases_MeetEndpointNodeCountAwayFromExceptionalSet()
    {
        var reading = new NodePairResolventWitness().Reading;

        Assert.NotEmpty(reading.CriterionCases);
        Assert.All(reading.CriterionCases, item =>
        {
            Assert.True(item.R != BigRational.Zero);
            Assert.True(item.R != BigRational.One);
            Assert.True(item.R != -BigRational.One);
            Assert.Equal(item.EndpointNodeCount, item.DetunedBlindCount);
        });
        Assert.Contains(reading.CriterionCases, item => item.EndpointNodeCount == 0);
        Assert.Contains(reading.CriterionCases, item => item.EndpointNodeCount > 0);
    }

    [Fact]
    public void ExceptionalSet_IsComputedAndFencedRatherThanFedToTheIff()
    {
        var reading = new NodePairResolventWitness().Reading;

        Assert.Equal(new[] { -BigRational.One, BigRational.Zero, BigRational.One },
            reading.ExceptionalCases.Select(item => item.R).OrderBy(item => item.Sign).ToArray());
        Assert.All(reading.ExceptionalCases, item => Assert.False(item.InTheoremDomain));
        Assert.Contains(reading.ExceptionalCases,
            item => item.R == BigRational.One && item.DetunedBlindCount > item.EndpointNodeCount);
        Assert.Contains(reading.ExceptionalCases,
            item => item.R == -BigRational.One && item.DetunedBlindCount > item.EndpointNodeCount);
        Assert.Contains(reading.ExceptionalCases,
            item => item.R == BigRational.Zero && item.Fence.Contains("zero bond"));
    }

    [Fact]
    public void InspectionTree_SeparatesLiveRecomputationsFromStoredScopeBreadcrumb()
    {
        IInspectable witness = new NodePairResolventWitness();
        var children = witness.Children.ToArray();

        Assert.Contains("exact", witness.Summary, StringComparison.OrdinalIgnoreCase);
        Assert.Contains(children, child => child.DisplayName.Contains("nonzero control")
                                          && child.Provenance == NodeProvenance.Live);
        Assert.Contains(children, child => child.DisplayName.Contains("uniform-centre iff")
                                          && child.Provenance == NodeProvenance.Live);
        Assert.Contains(children, child => child.DisplayName.Contains("exceptional set")
                                          && child.Provenance == NodeProvenance.Live);
        Assert.Contains(children, child => child.DisplayName.Contains("scope and proof")
                                          && child.Provenance == NodeProvenance.Stored);
    }
}
