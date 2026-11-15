using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Quadrupolar;
using Xunit;
using M = RCPsiSquared.Core.Numerics.GaussianRationalMatrix;

namespace RCPsiSquared.Core.Tests.Quadrupolar;

/// <summary>The Kronecker sum is the generator of two independent parts on the product space: on a
/// product state ρ₁ ⊗ ρ₂ it gives L₁ρ₁ ⊗ ρ₂ + ρ₁ ⊗ L₂ρ₂, on the product space's row-major vec.</summary>
public class GeneratorCompositionTests
{
    [Fact]
    public void KroneckerSum_ActsAsTheLeibnizRuleOnProducts()
    {
        var l1 = new QuadrupolarGenerator(SpinSpace.Of(2), new QuadrupolarRates(1, 2, 3), SecularHamiltonian.Of(1, 2)).Gauged();
        var l2 = new QuadrupolarGenerator(SpinSpace.Of(1), new QuadrupolarRates(0, 0, 0), SecularHamiltonian.Of(0, 5)).Gauged();
        int d1 = 3, d2 = 2, d = 6;
        var sum = GeneratorComposition.KroneckerSum(l1, d1, l2, d2);

        var rho1 = M.Zero(d1); rho1[0, 1] = 1; rho1[2, 2] = new GaussianRational(2, 1);
        var rho2 = M.Zero(d2); rho2[1, 0] = 3; rho2[0, 0] = -1;
        var lhs = Apply(sum, Kron(rho1, rho2), d);
        var rhs = M.Add(Kron(Apply(l1, rho1, d1), rho2), Kron(rho1, Apply(l2, rho2, d2)));
        Assert.True(M.AreEqual(lhs, rhs));
    }

    private static GaussianRational[,] Kron(GaussianRational[,] a, GaussianRational[,] b)
    {
        int n = a.GetLength(0), m = b.GetLength(0);
        var r = M.Zero(n * m);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                for (int k = 0; k < m; k++)
                    for (int l = 0; l < m; l++) r[i * m + k, j * m + l] = a[i, j] * b[k, l];
        return r;
    }

    private static GaussianRational[,] Apply(GaussianRational[,] super, GaussianRational[,] rho, int d)
    {
        var y = M.Zero(d);
        for (int r = 0; r < d * d; r++)
        {
            var acc = GaussianRational.Zero;
            for (int c = 0; c < d * d; c++) acc += super[r, c] * rho[c / d, c % d];
            y[r / d, r % d] = acc;
        }
        return y;
    }
}
