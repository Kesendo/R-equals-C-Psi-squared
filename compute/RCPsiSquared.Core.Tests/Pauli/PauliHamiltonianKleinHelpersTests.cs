using System.Numerics;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Tests.Pauli;

public class PauliHamiltonianKleinHelpersTests
{
    private static PauliHamiltonian Build(int N, params (int Site1, PauliLetter L1, int Site2, PauliLetter L2)[] terms)
    {
        var ts = terms.Select(t => PauliTerm.TwoSite(N, t.Site1, t.L1, t.Site2, t.L2, Complex.One)).ToList();
        return new PauliHamiltonian(N, ts);
    }

    [Fact]
    public void XXPlusYY_KleinSetIsMother_AndIsKleinHomogeneous()
    {
        // XX + YY: both bilinears in Klein (0,0) Mother (truly).
        var H = Build(3, (0, PauliLetter.X, 1, PauliLetter.X), (0, PauliLetter.Y, 1, PauliLetter.Y));
        Assert.Single(H.KleinSet);
        Assert.Contains((0, 0), H.KleinSet);
        Assert.True(H.IsKleinHomogeneous);
        Assert.True(H.IsZ2Homogeneous);
    }

    [Fact]
    public void XYPlusYX_KleinSetIsZeroOne_AndIsKleinHomogeneous()
    {
        // XY + YX: both Klein (0,1) — Π²-odd subgroup.
        var H = Build(3, (0, PauliLetter.X, 1, PauliLetter.Y), (0, PauliLetter.Y, 1, PauliLetter.X));
        Assert.Single(H.KleinSet);
        Assert.Contains((0, 1), H.KleinSet);
        Assert.True(H.IsKleinHomogeneous);
    }

    [Fact]
    public void YZPlusZY_KleinSetIsOneZero_Pi2EvenNontruly()
    {
        // YZ + ZY: Klein (1,0) — Π²-even non-truly.
        var H = Build(3, (0, PauliLetter.Y, 1, PauliLetter.Z), (0, PauliLetter.Z, 1, PauliLetter.Y));
        Assert.Single(H.KleinSet);
        Assert.Contains((1, 0), H.KleinSet);
        Assert.True(H.IsKleinHomogeneous);
    }

    [Fact]
    public void XXPlusXY_KleinSetHasTwo_NotKleinHomogeneous()
    {
        // XX in (0,0) Mother; XY in (0,1) Π²-odd → mixed Klein.
        var H = Build(3, (0, PauliLetter.X, 1, PauliLetter.X), (0, PauliLetter.X, 1, PauliLetter.Y));
        Assert.Equal(2, H.KleinSet.Count);
        Assert.Contains((0, 0), H.KleinSet);
        Assert.Contains((0, 1), H.KleinSet);
        Assert.False(H.IsKleinHomogeneous);
    }

    [Fact]
    public void IdentityTerms_AreExcludedFromKleinSet()
    {
        // I⊗I has KBody=0 → filtered. Add it alongside XX and verify KleinSet still {(0,0)}.
        var ts = new List<PauliTerm>
        {
            PauliTerm.TwoSite(3, 0, PauliLetter.X, 1, PauliLetter.X, Complex.One),
            new PauliTerm(new[] { PauliLetter.I, PauliLetter.I, PauliLetter.I }, Complex.One),
        };
        var H = new PauliHamiltonian(3, ts);
        Assert.Single(H.KleinSet);
        Assert.Contains((0, 0), H.KleinSet);
    }

    [Fact]
    public void XYZ_AtK3_IsKleinHomogeneousButZ2HomogeneityRefinesViaYParity()
    {
        // X⊗Y⊗Z (k=3): bit_a = 1+1+0 = 0; bit_b = 0+1+1 = 0; Y-par = 1.
        // Klein = (0,0) (looks like Mother), but Y-par = 1 distinguishes it from
        // a true Mother term. Demonstrates Z₂³ > Z₂² at k=3.
        var t1 = new PauliTerm(new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.Z }, Complex.One);
        var H1 = new PauliHamiltonian(3, new[] { t1 });
        Assert.Equal((0, 0), t1.KleinIndex);
        Assert.Equal((0, 0, 1), t1.FullZ2Signature);
        Assert.True(H1.IsKleinHomogeneous);
        Assert.True(H1.IsZ2Homogeneous); // single term

        // Add XX⊗I (Klein (0,0), Y-par 0) → still Klein-homogeneous but NOT Z₂-homogeneous.
        var t2 = new PauliTerm(new[] { PauliLetter.X, PauliLetter.X, PauliLetter.I }, Complex.One);
        Assert.Equal((0, 0), t2.KleinIndex);
        Assert.Equal((0, 0, 0), t2.FullZ2Signature);
        var H2 = new PauliHamiltonian(3, new[] { t1, t2 });
        Assert.True(H2.IsKleinHomogeneous);     // both Klein (0,0)
        Assert.False(H2.IsZ2Homogeneous);       // Y-par 0 vs 1 differ
        Assert.Equal(2, H2.FullZ2SignatureSet.Count);
    }
}
