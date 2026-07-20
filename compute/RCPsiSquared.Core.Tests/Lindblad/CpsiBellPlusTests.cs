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
    public void PureX_EqualsPureY_L1IsPinnedCorrelation()
    {
        // Pure Y has physical rates α=4γ, β=0 (violating the WLOG α ≤ β), so L₁ = max(u,v) = 1:
        // the YY correlation is pinned, CΨ = (1+u²)/6, functionally identical to pure X, NOT pure Z
        // (PROOF_MONOTONICITY_CPSI.md Part 2 K table; F27's 2026-06-22 revert of the K_Y=K_Z error).
        double t = 0.5;
        double cx = CpsiBellPlus.At(1.0, 0.0, 0.0, t);
        double cy = CpsiBellPlus.At(0.0, 1.0, 0.0, t);
        double cz = CpsiBellPlus.At(0.0, 0.0, 1.0, t);
        Assert.Equal(cx, cy, 12);
        Assert.NotEqual(cz, cy);
    }

    [Fact]
    public void PureY_AtCuspK_IsOneQuarter()
    {
        // K_Y = ln(2)/8 ≈ 0.0867: CΨ = (1+u²)/6 = 1/4 at u² = e^{-8γt} = 1/2.
        double cpsi = CpsiBellPlus.At(0.0, 1.0, 0.0, CpsiBellPlus.CuspK.PureY);
        Assert.InRange(cpsi, 0.249, 0.251);
    }
}
