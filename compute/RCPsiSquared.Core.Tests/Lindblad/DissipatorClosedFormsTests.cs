using System.Numerics;
using RCPsiSquared.Core.Lindblad;

namespace RCPsiSquared.Core.Tests.Lindblad;

public class DissipatorClosedFormsTests
{
    [Fact]
    public void T1_PauliDecomposition_ReproducesC1C2()
    {
        // T1 σ⁻ = 0.5 X − 0.5i Y; expected c1 = 3, c2 = 4 from HARDWARE_DISSIPATORS table.
        var (c1, c2) = DissipatorClosedForms.C1C2FromPauli(0.5, new Complex(0, -0.5), 0.0);
        Assert.Equal(3.0, c1, 10);
        Assert.Equal(4.0, c2, 10);
    }

    [Fact]
    public void Tphi_PureZ_HasC1Zero_C2Sixteen()
    {
        var (c1, c2) = DissipatorClosedForms.C1C2FromPauli(0.0, 0.0, 1.0);
        Assert.Equal(0.0, c1, 10);
        Assert.Equal(16.0, c2, 10);
    }

    [Fact]
    public void XNoise_HasC1Sixteen_C2Sixteen()
    {
        var (c1, c2) = DissipatorClosedForms.C1C2FromPauli(1.0, 0.0, 0.0);
        Assert.Equal(16.0, c1, 10);
        Assert.Equal(16.0, c2, 10);
    }

    [Fact]
    public void YNoise_HasC1Zero_C2Sixteen()
    {
        var (c1, c2) = DissipatorClosedForms.C1C2FromPauli(0.0, 1.0, 0.0);
        Assert.Equal(0.0, c1, 10);
        Assert.Equal(16.0, c2, 10);
    }

    [Fact]
    public void D2_T1_T1_HasExpectedValue()
    {
        // T1 ⊗ T1: ‖c‖² = 0.25 + 0.25 + 0 = 0.5. d2 = 32 · 0.5 · 0.5 = 8.
        double d2 = DissipatorClosedForms.D2FromPauli(0.5, new Complex(0, -0.5), 0.0,
                                                     0.5, new Complex(0, -0.5), 0.0);
        Assert.Equal(8.0, d2, 10);
    }
}
