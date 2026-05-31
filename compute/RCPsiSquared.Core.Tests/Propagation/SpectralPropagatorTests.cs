using System;
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
}
