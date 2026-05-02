using RCPsiSquared.Core.Lindblad;

namespace RCPsiSquared.Core.Tests.Lindblad;

public class CpsiBellPlusTests
{
    [Fact]
    public void AtTimeZero_BellPlus_HasCpsiOneThird()
    {
        // At t=0 all u=v=w=1, so CΨ = 1·(1+1+1+1)/12 = 4/12 = 1/3.
        Assert.Equal(1.0 / 3.0, CpsiBellPlus.At(0.05, 0.03, 0.07, t: 0.0), 12);
    }

    [Fact]
    public void AtLargeTime_DecaysToZero()
    {
        // For nonzero rates, all of u, v, w → 0 as t → ∞, so CΨ → 0.
        Assert.True(CpsiBellPlus.At(0.05, 0.03, 0.07, t: 100.0) < 1e-6);
    }

    [Fact]
    public void PureZ_AtCuspK_IsOneQuarter()
    {
        // K_Z = 0.0374 corresponds to γ·t at which CΨ first crosses 1/4 under pure Z noise.
        // For pure Z (γ_x = γ_y = 0): α = 4γ_z, β = 4γ_z, δ = 0 → u = v = exp(-4γt), w = 1.
        // CΨ = u(1 + u² + v² + 1)/12 = u(1 + 2u² + 1)/12 = u(2 + 2u²)/12 = u(1+u²)/6.
        // At K = γt = 0.0374: u = exp(-4·0.0374) = exp(-0.1496) ≈ 0.86099
        // CΨ = 0.86099·(1 + 0.86099²)/6 ≈ 0.86099·1.7413/6 ≈ 0.2499 ≈ 1/4 ✓
        double cpsi = CpsiBellPlus.At(0.0, 0.0, 1.0, CpsiBellPlus.CuspK.PureZ);
        Assert.InRange(cpsi, 0.249, 0.251);
    }

    [Fact]
    public void PureZ_EqualsPureY_ByF26Symmetry()
    {
        // F26 functional symmetry: pure Y has α=4γ_y, β=0, δ=4γ_y → identical CΨ structure.
        double t = 0.5;
        double cz = CpsiBellPlus.At(0.0, 0.0, 1.0, t);
        double cy = CpsiBellPlus.At(0.0, 1.0, 0.0, t);
        Assert.Equal(cz, cy, 12);
    }
}
