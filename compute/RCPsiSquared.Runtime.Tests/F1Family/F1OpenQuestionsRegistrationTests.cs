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
        // F1OpenQuestions.Standard is empty as of 2026-05-18 (all four F1 items
        // closed in the May 2026 sprint); the registered collection inherits the
        // empty list. The Same-reference assertion at index 0 from the pre-closure
        // version is no longer applicable.
    }

    [Fact]
    public void RegisterF1OpenQuestions_StandardIsEmpty()
    {
        // F1 family is open-question-free as of 2026-05-18 — first time the
        // collection is empty. See F1OpenQuestions XML doc for the per-item
        // closure references (T1 closed form, depol closed form, non-uniform γ
        // negative result, general topology synthesis proof + verification record).
        var registry = BuildWithF1OpenQuestions();
        var collection = registry.Get<OpenQuestionCollection<F1Marker>>();
        Assert.Empty(collection.Items);
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
