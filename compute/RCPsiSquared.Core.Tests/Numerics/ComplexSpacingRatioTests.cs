using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.Numerics;

/// <summary>The complex spacing ratio (Sá, Ribeiro, Prosen, PRX 2020), the non-Hermitian RMT
/// diagnostic: z_k = (NN−λ_k)/(NNN−λ_k) in C. The two reference classes are CALCULATED here (not
/// hardcoded) via the same diagnostic and must separate: a 2D-Poisson cloud → ⟨|z|⟩≈0.658,
/// ⟨cos θ⟩≈0 (no repulsion); a GinUE spectrum → ⟨|z|⟩≈0.738, ⟨cos θ⟩≈−0.241 (level repulsion +
/// angular avoidance). The references live in the production helper so the witness reuses them.</summary>
public class ComplexSpacingRatioTests
{
    [Fact]
    public void Poisson2D_Reference_Reads_AbsZ_Near0p658_AndCos_Near0()
    {
        var (meanAbs, meanCos) = ComplexSpacingRatio.PoissonDiskReference(count: 4000, seed: 1);
        Assert.InRange(meanAbs, 0.63, 0.69);
        Assert.InRange(meanCos, -0.04, 0.04);
    }

    [Fact]
    public void Ginue_Reference_Reads_AbsZ_Near0p738_AndCos_Negative()
    {
        var (meanAbs, meanCos) = ComplexSpacingRatio.GinueReference(n: 500, seed: 2);
        Assert.InRange(meanAbs, 0.71, 0.77);
        Assert.InRange(meanCos, -0.30, -0.16);
    }

    [Fact]
    public void TheTwoClasses_Separate_GinueRepelsPoissonDoesNot()
    {
        var pois = ComplexSpacingRatio.PoissonDiskReference(count: 4000, seed: 3);
        var gin = ComplexSpacingRatio.GinueReference(n: 500, seed: 4);
        Assert.True(gin.meanAbs > pois.meanAbs, "GinUE repels ⟹ larger ⟨|z|⟩ than Poisson");
        Assert.True(gin.meanCos < pois.meanCos - 0.1, "GinUE has angular avoidance (⟨cos θ⟩ < 0)");
    }
}
