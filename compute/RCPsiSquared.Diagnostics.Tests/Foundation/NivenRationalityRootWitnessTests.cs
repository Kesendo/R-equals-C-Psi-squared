using System;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live Niven-root witness recomputes the three cyclotomic faces and the N=4 first-golden hinge.
/// Pure arithmetic; the sympy-exact minimal-polynomial proof is simulations/niven_rationality_root.py.</summary>
public class NivenRationalityRootWitnessTests
{
    [Fact]
    public void Witness_RendersFiveFaces_AllGatesPass()
    {
        var w = new NivenRationalityRootWitness();
        var kids = ((IInspectable)w).Children.ToList();
        Assert.Equal(5, kids.Count);
        // every gated face must report "gate ✓" in its display name (no GATE FIRED)
        foreach (var k in kids)
            Assert.DoesNotContain("GATE FIRED", k.DisplayName);
        Assert.Contains(kids, k => k.DisplayName.Contains("gate ✓"));   // at least the gated faces pass
    }

    [Fact]
    public void Hinge_N4_BandEdgeIsPhi_AndRatesCarrySqrt5()
    {
        // band edge at N=4 = φ = 2cos(π/5)
        double be4 = 2.0 * Math.Cos(Math.PI / 5.0);
        Assert.Equal((1.0 + Math.Sqrt(5.0)) / 2.0, be4, 12);
        // the N=4 rate α_1/γ₀ = (4/5)sin²(π/5) = (5−√5)/10 (golden family)
        double rate41 = 4.0 / 5.0 * Math.Pow(Math.Sin(Math.PI / 5.0), 2);
        Assert.Equal((5.0 - Math.Sqrt(5.0)) / 10.0, rate41, 12);
    }

    [Fact]
    public void VFace_GoldenAtN5_SilverAtN4()
    {
        // the V-Effect face: golden shifted to N=5 (angle π/N), silver at N=4
        Assert.Equal((5.0 + Math.Sqrt(5.0)) / 4.0, 1.0 + Math.Cos(Math.PI / 5.0), 12);   // golden, N=5
        Assert.Equal(1.0 + Math.Sqrt(2.0) / 2.0, 1.0 + Math.Cos(Math.PI / 4.0), 12);     // silver, N=4
    }
}
