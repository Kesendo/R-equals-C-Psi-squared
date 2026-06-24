using System;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The q-parametric monodromy of the F89 path-3 octic: where the Galois structure actually
/// lives spectrally. The fixed-q geometry was a null (--root galoischaos); here the octic's 8 roots are
/// tracked as q loops the complex plane. Gate G2 (this file): a loop around the diabolic point q_EP
/// returns the IDENTITY (the two coalescing roots do not braid), confirming f89octic's "semisimple, not
/// defective" by an INDEPENDENT route (monodromy, not the Riesz projector). The monodromy machine this
/// builds is the gateway to the S_8-generation gate (Galois = monodromy, from below).</summary>
public class GaloisMonodromyWitnessTests
{
    [Fact]
    public void LoopAroundTheDiabolicPoint_ReturnsIdentity_SemisimpleNotDefective()
    {
        // at an ISOLATING radius (only q_EP enclosed) the diabolic loop does not braid: the two
        // coalescing roots return to themselves (double discriminant zero, transversal crossing).
        var perm = GaloisMonodromyWitness.MonodromyAroundDiabolic(radius: 0.02, steps: 600);
        Assert.Equal(12, perm.Length);
        for (int i = 0; i < perm.Length; i++)
            Assert.Equal(i, perm[i]);
    }

    [Fact]
    public void WiderLoop_EnclosesGenuineBranchPoints_CarryingTranspositions()
    {
        // a radius-0.1 loop also encloses genuine simple branch points (zeros of P_10 = defective EPs)
        // flanking q_EP; each carries a transposition. This is the seed of the S_8-generation gate (G3):
        // the diabolic itself is silent, but real EPs nearby DO braid.
        var perm = GaloisMonodromyWitness.MonodromyAroundDiabolic(radius: 0.1, steps: 600);
        int moved = Enumerable.Range(0, perm.Length).Count(i => perm[i] != i);
        Assert.True(moved >= 2 && moved % 2 == 0,
            $"genuine branch points braid (even # of moved strands, got {moved}), unlike the silent diabolic");
    }

    [Fact]
    public void TheTwoOcticRoots_Coalesce_AtTheDiabolicPoint()
    {
        var (qMin, lamMid, gap) = GaloisMonodromyWitness.ClosestApproachOnRealAxis();
        Assert.InRange(qMin, 0.60, 0.72);              // q_EP = sqrt((-1+sqrt13)/6) ≈ 0.659
        Assert.True(gap < 0.06, $"the two octic roots nearly coincide there (gap={gap:F4})");
        Assert.InRange(lamMid.Real, -4.3, -3.7);        // λ_EP real part = −4γ
        Assert.InRange(lamMid.Imaginary, 1.0, 1.7);     // λ_EP imag part = 2J = 2·q_EP ≈ 1.32
    }

    [Fact]
    public void TheMonodromyReconstructsS8FromBelow_TheTranspositionGraphIsConnected()
    {
        // every EP lassoed from a common base (detour-routed over the q=0 super-branch), its braid read on
        // the 8 octic strands in one labelling; a lasso that nets a single k-cycle certifies those strands
        // are one monodromy orbit. The graph is CONNECTED on all 8 strands ⟹ the transpositions generate the
        // FULL symmetric group: Gal(F_8) = S_8, reconstructed purely from eigenvalue braids (monodromy =
        // Galois), an independent route to the algebraic Frobenius certificate.
        var g3 = GaloisMonodromyWitness.GeneratesS8();
        Assert.Equal(8, g3.largest);            // all 8 octic strands in one component
        Assert.Equal(1, g3.components);
        Assert.True(g3.connected, "the transposition graph is connected ⟹ Gal(F_8) = S_8 from below");
    }

    [Fact]
    public void Witness_Renders_TheMonodromyIdentityVerdict()
    {
        var children = ((IInspectable)new GaloisMonodromyWitness()).Children.ToList();
        Assert.Contains(children, c =>
            c.Summary.Contains("monodromy", StringComparison.OrdinalIgnoreCase) &&
            c.Summary.Contains("identity", StringComparison.OrdinalIgnoreCase));
    }
}
