using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.F89PathK;

public class F89Path3OcticBlockTests
{
    private static double QEp => F89Path3OcticEpClaim.QEp;          // √((−1+√13)/6)

    [Fact]
    public void SymBlock_AtQEp_HasCoalescingPairAt_MinusFourGammaPlusTwoIJ()
    {
        var l = F89Path3OcticBlock.BuildSeDeSymBlock(QEp, 1.0);     // γ=1, J=q_EP
        Assert.Equal(12, l.RowCount);
        var evals = l.Evd().EigenValues.ToArray();

        // the octic = the 8 non-AT-locked modes (rate −Re ∉ {2,6})
        var octic = evals.Where(e => System.Math.Abs(-e.Real - 2.0) > 1e-6
                                  && System.Math.Abs(-e.Real - 6.0) > 1e-6).ToArray();
        Assert.Equal(8, octic.Length);

        // exactly one coalescing pair at λ_EP = −4 + 2i·q_EP
        var lamEp = F89Path3OcticEpClaim.MergedEigenvalue(1.0, QEp);
        int nearEp = octic.Count(e => (e - lamEp).Magnitude < 1e-6);
        Assert.Equal(2, nearEp);                                   // the double root
    }

    [Fact]
    public void SymBlock_SpectrumIsASubSpectrumOfTheFullN4Liouvillian()
    {
        // The (SE,DE) coherence block is an invariant sub-block (XY conserves Sz;
        // Z-dephasing is computational-basis-diagonal), so its spectrum ⊆ the full 4^N L.
        // This certifies the hand-rolled block is the real object (review finding A.1).
        double j = QEp, g = 1.0;
        var blockEvals = F89Path3OcticBlock.BuildSeDeSymBlock(j, g).Evd().EigenValues.ToArray();
        var fullEvals = F89BlockLiouvillian.BuildBlockL(j, g, 4).Evd().EigenValues.ToArray();
        foreach (var s in blockEvals)
        {
            double nearest = fullEvals.Min(f => (f - s).Magnitude);
            Assert.True(nearest < 1e-9, $"block eigenvalue {s} is not in the full spectrum (nearest {nearest:E2})");
        }
    }

    [Fact]
    public void EpCharacter_AtQEp_ReadsDiabolic_NotDefective()
    {
        // Gate-first: the artifact-free instrument must read DIABOLIC at the exact double root.
        // Load-bearing facts (basis-independent): Kind, geo=alg=2, dep≈0. ‖P‖ is grid-sensitive
        // and is NOT the diabolic discriminator (a closed-block defective EP also reads bounded
        // ‖P‖ — see Gate0a); the verdict rests on geo=alg & dep≈0. ‖P‖<1e3 is only a loose
        // finiteness sanity bound.
        double g = 1.0, j = QEp;
        var l = F89Path3OcticBlock.BuildSeDeSymBlock(j, g);
        var center = F89Path3OcticEpClaim.MergedEigenvalue(g, j);    // λ_EP = −4γ + 2iJ
        var r = EpCharacter.Characterize(l, center, radius: 0.5);

        Assert.Equal(EpCharacter.EpKind.Diabolic, r.Kind);          // the verdict
        Assert.Equal(2, r.Algebraic);
        Assert.Equal(2, r.Geometric);                              // geo == alg ⟹ diabolic
        Assert.True(r.Departure < 1e-6, $"diabolic departure {r.Departure} should be ≈ 0");
        Assert.True(r.ProjectorNorm < 1e3, $"‖P‖ {r.ProjectorNorm} should be a bounded sanity value (not a diabolic discriminator)");
        Assert.True(r.EigenvectorMergeCos < 0.99, $"|cos| {r.EigenvectorMergeCos} must stay < 1 (independent eigenvectors)");
    }
}
