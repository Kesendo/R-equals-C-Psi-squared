using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.Tests.F1Family;

public class F1OpenQuestionsRegistrationTests
{
    private static ClaimRegistry BuildWithF1OpenQuestions() =>
        new ClaimRegistryBuilder()
            .RegisterF1Family(new RCPsiSquared.Core.ChainSystems.ChainSystem(
                N: 5, J: 1.0, GammaZero: 0.05,
                HType: RCPsiSquared.Core.ChainSystems.HamiltonianType.XY,
                Topology: RCPsiSquared.Core.ChainSystems.TopologyKind.Chain))
            .RegisterF1OpenQuestions()
            .Build();

    [Fact]
    public void RegisterF1OpenQuestions_AddsCollection()
    {
        var registry = BuildWithF1OpenQuestions();
        Assert.True(registry.Contains<OpenQuestionCollection<F1Marker>>());
    }

    [Fact]
    public void RegisterF1OpenQuestions_TierIsOpenQuestion()
    {
        var registry = BuildWithF1OpenQuestions();
        var collection = registry.Get<OpenQuestionCollection<F1Marker>>();
        Assert.Equal(Tier.OpenQuestion, collection.Tier);
    }

    [Fact]
    public void RegisterF1OpenQuestions_ItemsMatchStandard()
    {
        var registry = BuildWithF1OpenQuestions();
        var collection = registry.Get<OpenQuestionCollection<F1Marker>>();
        Assert.Equal(F1OpenQuestions.Standard.Count, collection.Items.Count);
        Assert.Same(F1OpenQuestions.Standard[0], collection.Items[0]);
    }

    [Fact]
    public void RegisterF1OpenQuestions_AncestorContainsF1Palindrome()
    {
        var registry = BuildWithF1OpenQuestions();
        var ancestors = registry.AncestorsOf<OpenQuestionCollection<F1Marker>>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F1PalindromeIdentity), ancestors);
    }
}
