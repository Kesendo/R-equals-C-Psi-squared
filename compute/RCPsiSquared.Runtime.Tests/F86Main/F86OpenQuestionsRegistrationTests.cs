using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F86Main;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.F86Main;

public class F86OpenQuestionsRegistrationTests
{
    private static ClaimRegistry BuildRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterAbsorptionTheoremClaim()
            .RegisterF1PalindromeIdentity()
            .RegisterF86Main(gammaZero: 0.05, gEff: 1.0)
            .RegisterF86OpenQuestions()
            .Build();

    [Fact]
    public void FullSrpQuestion_IsExplicitAndInheritsFromSecuredPTypeClaim()
    {
        var registry = BuildRegistry();
        Assert.Equal(Tier.OpenQuestion, registry.Get<FullIrreducibleSrpClassQuestion>().Tier);
        var ancestors = registry.AncestorsOf<FullIrreducibleSrpClassQuestion>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(ShiftedGeneratorSectorwisePClaim), ancestors);
    }

    [Fact]
    public void FamilyOpenQuestionCollection_InheritsFromExplicitSrpQuestion()
    {
        var registry = BuildRegistry();
        var ancestors = registry.AncestorsOf<OpenQuestionCollection<F86Marker>>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(FullIrreducibleSrpClassQuestion), ancestors);
        Assert.Contains(typeof(ShiftedGeneratorSectorwisePClaim), ancestors);
    }
}
