using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Propagation;
using Xunit;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Propagation;

public class SpectralPropagatorTests
{
    [Fact]
    public void PureDecay_MatchesClosedFormSurvival()
    {
        // L = diag(-1, -2). vec(rho0) = [1, 1]. Two observable covectors pick each component.
        var L = ComplexMatrix.Build.DenseOfDiagonalArray(new[] { new Complex(-1, 0), new Complex(-2, 0) });
        var rho0 = ComplexVector.Build.DenseOfArray(new[] { Complex.One, Complex.One });
        var oA = ComplexVector.Build.DenseOfArray(new[] { Complex.One, Complex.Zero });
        var oB = ComplexVector.Build.DenseOfArray(new[] { Complex.Zero, Complex.One });
        var taus = new[] { 0.0, 1.0, 2.0 };

        var result = SpectralPropagator.Evolve(L, rho0, new[] { oA, oB }, taus);

        Assert.Equal(1.0, result[0][0], 12);
        Assert.Equal(Math.Exp(-1.0), result[0][1], 12);
        Assert.Equal(Math.Exp(-2.0), result[0][2], 12);
        Assert.Equal(1.0, result[1][0], 12);
        Assert.Equal(Math.Exp(-2.0), result[1][1], 12);
        Assert.Equal(Math.Exp(-4.0), result[1][2], 12);
    }

    [Fact]
    public void EvolveWithSpectrum_ReturnsEigenvaluesAndSameObservables()
    {
        var L = ComplexMatrix.Build.DenseOfDiagonalArray(new[] { new Complex(-1, 0), new Complex(-2, 0) });
        var rho0 = ComplexVector.Build.DenseOfArray(new[] { Complex.One, Complex.One });
        var oA = ComplexVector.Build.DenseOfArray(new[] { Complex.One, Complex.Zero });
        var taus = new[] { 0.0, 1.0 };

        var result = SpectralPropagator.EvolveWithSpectrum(L, rho0, new[] { oA }, taus);

        var reals = result.Eigenvalues.Select(z => z.Real).OrderBy(x => x).ToArray();
        Assert.Equal(-2.0, reals[0], 12);
        Assert.Equal(-1.0, reals[1], 12);
        var plain = SpectralPropagator.Evolve(L, rho0, new[] { oA }, taus);
        Assert.Equal(plain[0][0], result.Observables[0][0], 12);
        Assert.Equal(plain[0][1], result.Observables[0][1], 12);
    }

    [Fact]
    public void EvolveStateVectors_ReturnsFullStateConsistentWithObservables()
    {
        var L = ComplexMatrix.Build.DenseOfDiagonalArray(new[] { new Complex(-1, 0), new Complex(-2, 0) });
        var rho0 = ComplexVector.Build.DenseOfArray(new[] { Complex.One, Complex.One });
        var oA = ComplexVector.Build.DenseOfArray(new[] { Complex.One, Complex.Zero });
        var taus = new[] { 0.0, 1.0, 2.0 };

        var states = SpectralPropagator.EvolveStateVectors(L, rho0, taus);

        // at τ=0 the state is exactly vec(ρ₀)
        Assert.Equal(1.0, states[0][0].Real, 12);
        Assert.Equal(1.0, states[0][1].Real, 12);
        // closed-form decay on each component
        Assert.Equal(Math.Exp(-1.0), states[1][0].Real, 12);
        Assert.Equal(Math.Exp(-4.0), states[2][1].Real, 12);
        // the full state reproduces the linear observable: w · vec(ρ(τ)) == Evolve
        var obs = SpectralPropagator.Evolve(L, rho0, new[] { oA }, taus);
        for (int ti = 0; ti < taus.Length; ti++)
        {
            Complex acc = oA[0] * states[ti][0] + oA[1] * states[ti][1];
            Assert.Equal(obs[0][ti], acc.Real, 12);
        }
    }
}
