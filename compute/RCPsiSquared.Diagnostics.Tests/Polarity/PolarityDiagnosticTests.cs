using RCPsiSquared.Core.States;
using RCPsiSquared.Diagnostics.Polarity;

namespace RCPsiSquared.Diagnostics.Tests.Polarity;

public class PolarityDiagnosticTests
{
    [Fact]
    public void UniformPlus_HasFullPositivePolarityAxis()
    {
        // |+⟩^N has ⟨X_i⟩ = 1, ⟨Y_i⟩ = ⟨Z_i⟩ = 0 per site.
        var psi = PolarityState.Uniform(N: 3, sign: +1);
        var rho = DensityMatrix.FromStateVector(psi);
        var p = PolarityDiagnostic.FromDensityMatrix(rho);
        foreach (var x in p.PolarityAxis) Assert.Equal(1.0, x, 10);
        foreach (var o in p.OffAxis) Assert.Equal(0.0, o, 10);
        Assert.Equal(1.0, p.AggregatePolarity, 10);
        Assert.Equal(1.0, p.OnAxisFraction, 10);
    }

    [Fact]
    public void UniformMinus_HasFullNegativePolarityAxis()
    {
        var psi = PolarityState.Uniform(N: 3, sign: -1);
        var rho = DensityMatrix.FromStateVector(psi);
        var p = PolarityDiagnostic.FromDensityMatrix(rho);
        foreach (var x in p.PolarityAxis) Assert.Equal(-1.0, x, 10);
        Assert.Equal(1.0, p.AggregatePolarity, 10);
    }

    [Fact]
    public void NeelPattern_AlternatesPolarity()
    {
        // |+, -, +⟩
        var psi = PolarityState.Build(N: 3, signs: new[] { +1, -1, +1 });
        var rho = DensityMatrix.FromStateVector(psi);
        var p = PolarityDiagnostic.FromDensityMatrix(rho);
        Assert.Equal(+1.0, p.PolarityAxis[0], 10);
        Assert.Equal(-1.0, p.PolarityAxis[1], 10);
        Assert.Equal(+1.0, p.PolarityAxis[2], 10);
    }
}
