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
        // K_Z is γ·t at which CΨ first crosses 1/4 under pure Z noise. For pure Z
        // (γ_x = γ_y = 0): α = 4γ_z, β = 4γ_z, δ = 0 → u = v = exp(-4γt), w = 1, so
        // CΨ = u(1 + u² + v² + 1)/12 = u(1+u²)/6.
        //
        // The constant is a four-place reading of −¼·ln(u*), u* the real root of u³+u = 3/2,
        // so it carries a rounding of at most 5e-5 and the crossing is displaced by that times
        // the slope: |dCΨ/dK| = (1+3u*²)/6 · 4u* = 1.8516, giving 9.3e-5. That product IS the
        // error model, and the gate is it, not a round window ten times wider.
        double cpsi = CpsiBellPlus.At(0.0, 0.0, 1.0, CpsiBellPlus.CuspK.PureZ);
        const double slope = 1.8516;           // |dCΨ/dK| at the pure-Z crossing
        const double quoteRounding = 0.5e-4;   // half-width of the four-place constant
        Assert.True(Math.Abs(cpsi - 0.25) <= slope * quoteRounding,
            $"CΨ at K_Z should sit within the quoted constant's own rounding, was {cpsi:R}");
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
        // K_Y = ln(2)/8 = 0.0866433…: CΨ = (1+v²)/6 = 1/4 at v² = e^{-8γt} = 1/2.
        // In exact arithmetic the crossing is 1/4 exactly. In floating point the only error is
        // the rounding of Math.Exp and the divide, so the model is a few ulp of 1/4 and the gate
        // is that model rather than a decimal window a wrong fourth digit would fit through.
        double cpsi = CpsiBellPlus.At(0.0, 1.0, 0.0, CpsiBellPlus.CuspK.PureY);
        double ulp = Math.BitIncrement(0.25) - 0.25;
        Assert.True(Math.Abs(cpsi - 0.25) <= 4 * ulp,
            $"CΨ at K_Y should sit within a few ulp of 1/4, was {cpsi:R} (off by {Math.Abs(cpsi - 0.25) / ulp:F1} ulp)");
    }
}
