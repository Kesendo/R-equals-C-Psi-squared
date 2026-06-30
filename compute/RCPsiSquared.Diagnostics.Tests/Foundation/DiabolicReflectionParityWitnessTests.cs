using System;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The from-below grounding of the odd-N real-q diabolic onset (the dimension-mismatch / sector-swap
/// mechanism). The reflection R splits the full (SE,DE) block into R-even and R-odd; the realness
/// antiunitarity (Σ L Σ = L†) maps the two sectors onto each other at EVEN N (equal dims, σ_even = conj σ_odd)
/// so neither is self-conjugate and a real-axis collision is the generic pseudo-Hermitian defective EP, while
/// at ODD N the reflection-fixed central site makes dim(even) − dim(odd) = (N−1)/2 ≠ 0, forbidding the
/// cross-pairing, so each sector is self-conjugate and carries real eigenvalues that can cross semisimply
/// (the real-q diabolic). These tests pin that structure live, including the reproduction of the known
/// N=7/N=9 diabolics in the R-even sector. The witness for experiments/F89_PATH_K_DIABOLIC.md.</summary>
public class DiabolicReflectionParityWitnessTests
{
    [Theory]
    [InlineData(7, 3)]      // odd: dim(even) − dim(odd) = (N−1)/2 = #reflection-fixed singletons
    [InlineData(9, 4)]
    [InlineData(6, 0)]      // even: balanced sectors, no fixed central site
    [InlineData(8, 0)]
    public void DimMismatch_IsTheCenterSiteFingerprint(int nBlock, int expectedDiff)
    {
        var r = new DiabolicReflectionParityWitness().Read(nBlock, 1.0);
        Assert.Equal(expectedDiff, r.EvenDim - r.OddDim);
        Assert.Equal(expectedDiff, r.FixedSingletons);
    }

    [Theory]
    [InlineData(7)]
    [InlineData(9)]
    public void OddN_BothSectorsSelfConjugate(int nBlock)
    {
        var r = new DiabolicReflectionParityWitness().Read(nBlock, 1.0);
        Assert.True(r.SelfConjEven < 1e-7, $"R-even not self-conjugate at odd N={nBlock}: {r.SelfConjEven:E2}");
        Assert.True(r.SelfConjOdd < 1e-7, $"R-odd not self-conjugate at odd N={nBlock}: {r.SelfConjOdd:E2}");
    }

    [Theory]
    [InlineData(6)]
    [InlineData(8)]
    public void EvenN_SectorsAreConjugatePaired_NeitherSelfConjugate(int nBlock)
    {
        var r = new DiabolicReflectionParityWitness().Read(nBlock, 1.0);
        Assert.True(r.CrossConj < 1e-7, $"σ_even != conj σ_odd at even N={nBlock}: {r.CrossConj:E2}");
        Assert.True(r.SelfConjEven > 1.0, $"R-even unexpectedly self-conjugate at even N={nBlock}: {r.SelfConjEven:E2}");
    }

    [Theory]
    [InlineData(7, 1.1264, -4.942)]   // the C# scout diabolics, reproduced from below in the R-even sector
    [InlineData(9, 0.4755, -5.424)]
    public void REvenSector_ReproducesKnownDiabolic(int nBlock, double qRe, double lambda)
    {
        var (l1, _, gap) = new DiabolicReflectionParityWitness().ReproduceDiabolic(nBlock, qRe, lambda);
        Assert.True(Math.Abs(l1 - lambda) < 0.01, $"nearest R-even eigenvalue {l1:F4} far from {lambda}");
        Assert.True(gap < 1e-3, $"not a coalescence at q={qRe}: gap={gap:E2}");
    }

    [Fact]
    public void Summary_StatesTheParityVerdict()
    {
        var s = new DiabolicReflectionParityWitness().Summary;
        Assert.Contains("odd", s, StringComparison.OrdinalIgnoreCase);
        Assert.False(string.IsNullOrWhiteSpace(s));
    }
}
