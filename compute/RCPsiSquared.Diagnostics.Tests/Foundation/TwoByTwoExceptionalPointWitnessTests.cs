using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class TwoByTwoExceptionalPointWitnessTests
{
    private static readonly TwoByTwoExceptionalPointWitness Default = new();

    [Fact]
    public void EveryRow_JordanOnItsPoint_NoneOneStepOff()
    {
        Assert.Equal(13, Default.Rows.Count);
        Assert.Equal(9, Default.DistinctOnPointMatrices);   // three companions alike, ring = proton, edge = toy
        Assert.All(Default.Rows, r => Assert.True(r.Holds, r.Setting));
        Assert.Equal(3, Default.Rows.Count(r => r.Companion));
    }

    [Fact]
    public void SelfOrthogonality_FollowsTheWriting_AndDiagOneICarriesIt()
    {
        Assert.True(Default.SelfOrthogonalityFollowsTheWriting);
        Assert.True(Default.WritingCarriesSelfOrthogonality);   // (1, −1), vᵀv = 2 → (1, i), vᵀv = 0
        Assert.Contains(Default.Rows, r => r.SelfOrthogonal == false);
        Assert.Contains(Default.Rows, r => r.SelfOrthogonal == true);
    }

    [Fact]
    public void Dimer_SingleExcitationJordan_NoJumpGeneratorSizeThree()
    {
        Assert.True(Default.DimerBlockIsJordan);
        Assert.Equal(new[] { 10, 15, 16, 16 }, Default.DimerNoJumpNullities);
        Assert.True(Default.DimerNoJumpSizeThree);
        Assert.True(Default.Holds);
    }

    [Fact]
    public void Owners_TheToyModelAndQEpLaw_AgreeExactly()
    {
        Assert.True(Default.OwnerToyAgrees);    // TwoLevelEpModel at γ₀ = 1, J = 2, g_eff = 1: discriminant 0.0, centre −4
        Assert.True(Default.QEpLawSitsOnTheToyPoint);   // QEpLaw's 2/g_eff as J: d = 0, a quarter step off not
        Assert.Contains(Default.Rows, r => r.Setting.Contains("owner's point"));
    }

    [Fact]
    public void Law_ScalarAtTheDoubleRoot_IsNoJordanBlock()
    {
        // a scalar 2 × 2 has d = 0 and M₀ = 0: semisimple, the control for the M₀ ≠ 0 half of the rule
        var scalar = new GaussianRational[,] { { new(3, 0), GaussianRational.Zero }, { GaussianRational.Zero, new(3, 0) } };
        Assert.True(TwoByTwoExceptionalPointClaim.D(scalar).IsZero);
        Assert.False(TwoByTwoExceptionalPointClaim.IsJordanBlock(scalar));
    }
}
