using System;
using System.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Guards <see cref="ReferenceDefectiveLoci"/>: the Multi-Sector Monodromy Census reference loci,
/// the (1,2) octic's DEFECTIVE (braid-carrying P₁₀ √-branch) exceptional-point q-values at which the census
/// probes every other sector. The N=4 (path-3) octic has four REAL defective EPs at q≈0.460/0.854/0.857/1.738
/// (the 0.854/0.857 pair is a ~0.003-split near-degenerate twin the 0.05-cell monodromy lasso encloses
/// together; <see cref="PathKMonodromyScout.FindDiabolicsExact"/> at the tuned fine cell resolves them),
/// plus complex EPs; the silent SEMISIMPLE DIABOLIC at q≈0.659 (a (3q⁴+q²−1) root) is EXCLUDED. The loci
/// are conjugation-closed. Companion: <see cref="SectorEpProbe"/> reads the defective-vs-diabolic character
/// of one such coalescence.</summary>
public class ReferenceDefectiveLociTests
{
    [Fact]
    public void For_N4_HasFourRealDefectiveLoci_ExcludesDiabolic_ConjugationClosed()
    {
        var loci = ReferenceDefectiveLoci.For(4);
        double[] expectedReal = { 0.460, 0.854, 0.857, 1.738 }; // real-q P10 EPs (review from-below, ~3 digits)
        foreach (var e in expectedReal)
            Assert.Contains(loci, q => Math.Abs(q.Real - e) < 5e-3 && Math.Abs(q.Imaginary) < 1e-6);
        Assert.DoesNotContain(loci, q => Math.Abs(q.Real - 0.659) < 5e-3 && Math.Abs(q.Imaginary) < 1e-6); // the diabolic, excluded
        foreach (var q in loci) // conjugation-closed
            Assert.Contains(loci, r => (r - Complex.Conjugate(q)).Magnitude < 1e-4);
    }
}
