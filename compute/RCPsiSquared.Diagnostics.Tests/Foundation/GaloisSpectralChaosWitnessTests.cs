using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The Galois-vs-spectral-chaos witness: the sector-resolved RMT test that asks whether the
/// H_B-mixed half (whose decay rates have Galois group S_n, no radical closure) is spectrally CHAOTIC
/// (GinUE) at fixed q. The finding is a clean NULL: at accessible N the H_B half reads Poisson-like /
/// sub-Poisson (it sits on the integrable frequency lattice), NOT GinUE. Algebraic chaos (Galois over
/// the q-field) and spectral chaos (RMT at fixed q) are distinct here — the former does not imply the
/// latter. These tests pin that null on the trusted machine (the same (SE,DE) block, MathNet EVD).</summary>
public class GaloisSpectralChaosWitnessTests
{
    private static double[] Qs() => Enumerable.Range(0, 12).Select(i => 0.3 + i * (3.7 / 11)).ToArray();

    [Fact]
    public void ChainHbMixedSector_IsNotGinibre_ButPoissonLikeOrSubPoisson()
    {
        // chain N=7: the H_B-mixed factor is the S_53 irreducible (no radical closure). If the
        // Galois-chaos conjecture held it would read GinUE (⟨|z|⟩≈0.738). It does not.
        var (meanAbs, meanCos, avgCount) = GaloisSpectralChaosWitness.HbMixedCsr(7, "chain", Qs());
        Assert.True(avgCount > 20, $"need enough distinct H_B points for a CSR (got {avgCount:F0}/q)");
        Assert.True(meanAbs < 0.66, $"chain H_B ⟨|z|⟩={meanAbs:F3} is Poisson-like/sub-Poisson, NOT GinUE (~0.738)");
        Assert.True(meanCos > -0.10, $"chain H_B ⟨cos θ⟩={meanCos:F3} shows no GinUE angular repulsion (~−0.24)");
    }

    [Fact]
    public void TheMachineCanSeeGinibre_SoTheNullIsMeaningful()
    {
        // the null is only meaningful if the same diagnostic DOES flag a genuine Ginibre spectrum.
        var (ginAbs, _) = Core.Numerics.ComplexSpacingRatio.GinueReference(n: 300, seed: 7);
        Assert.True(ginAbs > 0.71, $"GinUE reference ⟨|z|⟩={ginAbs:F3} must read chaotic (~0.738)");
    }

    [Fact]
    public void Witness_Renders_TheNullVerdict()
    {
        var children = ((IInspectable)new GaloisSpectralChaosWitness()).Children.ToList();
        Assert.Contains(children, c =>
            c.Summary.Contains("not", System.StringComparison.OrdinalIgnoreCase) &&
            (c.Summary.Contains("Ginibre") || c.Summary.Contains("GinUE")));
    }
}
