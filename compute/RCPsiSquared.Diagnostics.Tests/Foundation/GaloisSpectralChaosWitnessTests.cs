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
    public void ChainAtLockedSector_IsPicketFence_NotGinibre()
    {
        // The other side of the door's comparison: the AT-locked half (rates −2γ/−6γ carrying free-fermion
        // Bloch frequencies) is a sparse structured set, not a GinUE chaos cloud. The discriminator is ⟨|z|⟩
        // (not the angle): it reads LOW (clustering, below GinUE's ~0.74) over far fewer distinct points
        // than the H_B half. Its ⟨cos θ⟩ can run negative (lattice angular order + few-point noise), so the
        // verdict rests on the magnitude, which a genuine GinUE spectrum would push high.
        var qs = Qs();
        var (atAbs, _, atCnt) = GaloisSpectralChaosWitness.AtLockedCsr(7, "chain", qs);
        var (_, _, hbCnt) = GaloisSpectralChaosWitness.HbMixedCsr(7, "chain", qs);
        Assert.False(double.IsNaN(atAbs), "AT-locked CSR should be defined at N=7");
        Assert.True(atAbs < 0.70, $"AT-locked ⟨|z|⟩={atAbs:F3} is below GinUE (~0.74): clustering, not a chaos cloud");
        Assert.True(atCnt < hbCnt, $"AT-locked is sparser ({atCnt:F0}/q) than the H_B half ({hbCnt:F0}/q): the picket-fence");
    }

    [Fact]
    public void TheMachineCanSeeGinibre_SoTheNullIsMeaningful()
    {
        // the null is only meaningful if the same diagnostic DOES flag a genuine Ginibre spectrum.
        var (ginAbs, _) = Core.Numerics.ComplexSpacingRatio.GinueReference(n: 300, seed: 7);
        Assert.True(ginAbs > 0.71, $"GinUE reference ⟨|z|⟩={ginAbs:F3} must read chaotic (~0.738)");
    }

    [Fact]
    public void Witness_Renders_BothHalves_AsNonChaotic()
    {
        var children = ((IInspectable)new GaloisSpectralChaosWitness()).Children.ToList();
        // the H_B-mixed null
        Assert.Contains(children, c =>
            c.Summary.Contains("not", System.StringComparison.OrdinalIgnoreCase) &&
            (c.Summary.Contains("Ginibre") || c.Summary.Contains("GinUE")));
        // the AT-locked picket-fence (the other side of the door's comparison)
        Assert.Contains(children, c =>
            c.DisplayName.Contains("AT-locked") && c.DisplayName.Contains("picket-fence"));
    }
}
