using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
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

    // ---- per-point z values (Reviewer B: pool the dimensionless z's across spectra, NEVER raw
    //      eigenvalues; needed for the Door-C CSR sweep's bootstrap + finite-size-matched references) ----

    /// <summary>A deterministic, distinct, low-discrepancy point cloud (golden-ratio Halton-like), so
    /// the per-point/aggregate consistency check does not depend on any RNG.</summary>
    private static List<Complex> DeterministicCloud(int n)
    {
        var pts = new List<Complex>(n);
        for (int k = 1; k <= n; k++)
            pts.Add(new Complex((k * 0.6180339887498949) % 1.0, (k * 0.7548776662466927) % 1.0));
        return pts;
    }

    /// <summary>The new per-point ZValues must reproduce the trusted aggregate Of() exactly: ⟨|z|⟩ and
    /// ⟨cos arg z⟩ pooled over the returned z's equal Of()'s meanAbs / meanCos to machine precision.</summary>
    [Fact]
    public void ZValues_PooledMean_MatchesOfAggregate_Exactly()
    {
        var pts = DeterministicCloud(300);
        var (ofAbs, ofCos, _) = ComplexSpacingRatio.Of(pts);
        var zs = ComplexSpacingRatio.ZValues(pts);

        Assert.True(zs.Count > 100, $"expected many z's, got {zs.Count}");
        Assert.Equal(ofAbs, zs.Average(z => z.Magnitude), precision: 12);
        Assert.Equal(ofCos, zs.Average(z => Math.Cos(z.Phase)), precision: 12);
    }

    /// <summary>Below the distinct-point gate (10) ZValues returns an empty list, mirroring Of()'s NaN
    /// gate, so a too-small spectrum contributes nothing to a pool.</summary>
    [Fact]
    public void ZValues_BelowGate_IsEmpty()
    {
        var pts = DeterministicCloud(7);
        Assert.Empty(ComplexSpacingRatio.ZValues(pts));
    }

    /// <summary>The finite-size reference z-value draws must agree with the existing aggregate
    /// references when pooled over the SAME seed (same cloud), so a pool of small finite-size draws is
    /// a faithful, bias-matched reference rather than the asymptotic 0.658/0.738.</summary>
    [Fact]
    public void PoissonAndGinueZValues_PooledMean_MatchTheirAggregateReferences()
    {
        var poisZ = ComplexSpacingRatio.PoissonDiskZValues(count: 4000, seed: 5);
        var poisAgg = ComplexSpacingRatio.PoissonDiskReference(count: 4000, seed: 5);
        Assert.Equal(poisAgg.meanAbs, poisZ.Average(z => z.Magnitude), precision: 12);

        var ginZ = ComplexSpacingRatio.GinueZValues(n: 300, seed: 6);
        var ginAgg = ComplexSpacingRatio.GinueReference(n: 300, seed: 6);
        Assert.Equal(ginAgg.meanAbs, ginZ.Average(z => z.Magnitude), precision: 12);
    }
}
