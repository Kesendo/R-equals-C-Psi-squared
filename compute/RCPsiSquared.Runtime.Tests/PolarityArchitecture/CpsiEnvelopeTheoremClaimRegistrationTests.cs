using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class CpsiEnvelopeTheoremClaimRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterAbsorptionTheoremClaim()
            .RegisterF25CPsiBellPlusPi2Inheritance();

    [Fact]
    public void Register_AddsClaim()
    {
        var registry = BuildBaseRegistry().RegisterCpsiEnvelopeTheoremClaim().Build();
        Assert.True(registry.Contains<CpsiEnvelopeTheoremClaim>());
    }

    [Fact]
    public void Register_TierIsOpenQuestion()
    {
        var registry = BuildBaseRegistry().RegisterCpsiEnvelopeTheoremClaim().Build();
        Assert.Equal(Tier.OpenQuestion, registry.Get<CpsiEnvelopeTheoremClaim>().Tier);
    }

    [Fact]
    public void Register_DirectGetDependencies_AreExactlyF25AndQuarter()
    {
        var registry = BuildBaseRegistry().RegisterCpsiEnvelopeTheoremClaim().Build();
        var directEdges = registry.EdgesInto<CpsiEnvelopeTheoremClaim>().ToList();
        var directParentTypes = directEdges.Select(edge => edge.Parent).ToHashSet();

        Assert.Equal(2, directEdges.Count);
        Assert.Equal(2, directParentTypes.Count);
        Assert.Equal(
            new[]
            {
                typeof(F25CPsiBellPlusPi2Inheritance),
                typeof(QuarterAsBilinearMaxvalClaim),
            }.OrderBy(type => type.FullName),
            directParentTypes.OrderBy(type => type.FullName));
    }

    [Fact]
    public void Register_ResolvedClaimShowsTheRisesItRecomputes()
    {
        // The registry-resolved claim carries both counterexample classes as live nodes: the maxima they
        // report are the ones the closed forms yield now, and in each pair the later one is higher.
        var claim = BuildBaseRegistry().RegisterCpsiEnvelopeTheoremClaim().Build()
            .Get<CpsiEnvelopeTheoremClaim>();
        var children = ((IInspectable)claim).Children.ToList();

        var (first, second) = CpsiEnvelopeTheoremClaim.LocalFieldFirstTwoMaxima();
        Assert.True(second.Cpsi > first.Cpsi);
        var fields = children.Single(c => c.DisplayName.StartsWith("local fields"));
        Assert.Equal(NodeProvenance.Live, fields.Provenance);
        Assert.Contains(first.Cpsi.ToString("0.000000", CultureInfo.InvariantCulture), fields.Summary);
        Assert.Contains(second.Cpsi.ToString("0.000000", CultureInfo.InvariantCulture), fields.Summary);

        var micro = CpsiEnvelopeTheoremClaim.MicroMaximumExample();
        Assert.True(micro.NextMainPeak.Cpsi > micro.Micro.Cpsi);
        var split = children.Single(c => c.DisplayName.StartsWith("number-conserving H"));
        Assert.Equal(NodeProvenance.Live, split.Provenance);
        Assert.Contains(micro.Micro.Cpsi.ToString("0.000000000", CultureInfo.InvariantCulture), split.Summary);
        Assert.Contains(micro.NextMainPeak.Cpsi.ToString("0.0000000", CultureInfo.InvariantCulture), split.Summary);
    }
}
