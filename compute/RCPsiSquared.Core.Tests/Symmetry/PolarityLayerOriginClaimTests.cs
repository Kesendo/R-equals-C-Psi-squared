using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.States;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class PolarityLayerOriginClaimTests
{
    [Fact]
    public void Claim_IsTier1Derived()
    {
        var claim = new PolarityLayerOriginClaim();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void Claim_AnchorReferences_ThePolarityLayerDoc()
    {
        var claim = new PolarityLayerOriginClaim();
        Assert.Contains("THE_POLARITY_LAYER.md", claim.Anchor);
    }

    [Fact]
    public void Claim_HasNamedCascadeChildren_FromSubstrateToReadings()
    {
        var claim = new PolarityLayerOriginClaim();
        IInspectable c = claim;
        var labels = c.Children.Select(ch => ch.DisplayName).ToList();
        // tier + anchor metadata (from base) + cascade children
        Assert.Contains("layer −1 (substrate)", labels);
        Assert.Contains("layer 0 (Π involution)", labels);
        Assert.Contains("layer 1 (X-basis diagonalization)", labels);
        Assert.Contains("layer 2 (the 0.5-shift)", labels);
        Assert.Contains("layer 3 (multi-axis refinement)", labels);
        Assert.Contains("layer 4 (bridges to the trio)", labels);
        Assert.Contains("reading (a): palindromic-over-the-layer", labels);
        Assert.Contains("reading (b): ARE-the-layer", labels);
    }

    [Fact]
    public void XPlusState_HasBlochXEqualsPlusOne_ConfirmingPolarityIsPlusZero()
    {
        // Operational anchor for layer 1: |+⟩ has Bloch x = +1, the +0 polarity.
        var psi = PolarityState.Uniform(N: 1, sign: +1);
        var rho = OuterProduct(psi);
        double bx = ((PauliString.Build(new[] { PauliLetter.X }) * rho).Trace()).Real;
        Assert.Equal(1.0, bx, precision: 12);
    }

    [Fact]
    public void XMinusState_HasBlochXEqualsMinusOne_ConfirmingPolarityIsMinusZero()
    {
        // Operational anchor for layer 1: |−⟩ has Bloch x = −1, the −0 polarity.
        var psi = PolarityState.Uniform(N: 1, sign: -1);
        var rho = OuterProduct(psi);
        double bx = ((PauliString.Build(new[] { PauliLetter.X }) * rho).Trace()).Real;
        Assert.Equal(-1.0, bx, precision: 12);
    }

    [Fact]
    public void HalfShift_AtPureXEigenstates_GivesOneAndZeroDiagonals()
    {
        // Layer 2 verification: ρ = (I + r·X)/2 at r = +1 has diagonals 1 and 0.
        // This IS the 0.5-shift: 1/2 + r/2 = 1/2 + 0.5 = 1, 1/2 − 0.5 = 0.
        var psi = PolarityState.Uniform(N: 1, sign: +1);
        var rho = OuterProduct(psi);
        // |+⟩⟨+| = (1/2 1/2; 1/2 1/2) — diagonals are both 1/2 in the Z basis,
        // because |+⟩⟨+| has 1/2 + 0.5·X/2-trace = 1/2 + 1/2 along X eigenbasis.
        // To witness the 0.5-shift directly: project onto X eigenbasis.
        // Equivalently: ⟨+|+⟩⟨+|+⟩ = 1, ⟨−|+⟩⟨+|−⟩ = 0.
        var plusPsi = PolarityState.Uniform(N: 1, sign: +1);
        var minusPsi = PolarityState.Uniform(N: 1, sign: -1);
        Complex pp = (plusPsi.ToColumnMatrix().ConjugateTranspose() * rho * plusPsi.ToColumnMatrix())[0, 0];
        Complex mm = (minusPsi.ToColumnMatrix().ConjugateTranspose() * rho * minusPsi.ToColumnMatrix())[0, 0];
        Assert.Equal(1.0, pp.Real, precision: 12);
        Assert.Equal(0.0, mm.Real, precision: 12);
    }

    [Fact]
    public void PolarityPair_CenteredOnD0Root_WithHalfMagnitudeAndUnitSpan()
    {
        // Pin the polarity pair {−a_2, +a_2} = {−0.5, +0.5} to three Pi2 identities at
        // precision 15: sum = 0 (d=0 root), span = 1, magnitude = a_2. Anchored to
        // Pi2DyadicLadder.Term(2) for rename-safety, not a hardcoded 0.5.
        var ladder = new Pi2DyadicLadderClaim();
        double half = ladder.Term(2);   // a_2 = 1/2 = HalfAsStructuralFixedPoint magnitude
        double plus = +half;
        double minus = -half;

        // pair sum equals the d=0 polynomial root (mirror axis)
        Assert.Equal(0.0, plus + minus, precision: 15);

        // pair span equals the unit interval (the 0.5-shift is half this span)
        Assert.Equal(1.0, plus - minus, precision: 15);

        // each side's magnitude equals a_2 = 1/2 = HalfAsStructuralFixedPoint
        Assert.Equal(half, Math.Abs(plus), precision: 15);
        Assert.Equal(half, Math.Abs(minus), precision: 15);
    }

    private static ComplexMatrix OuterProduct(MathNet.Numerics.LinearAlgebra.Vector<Complex> psi)
    {
        var col = psi.ToColumnMatrix();
        return col * col.ConjugateTranspose();
    }
}
