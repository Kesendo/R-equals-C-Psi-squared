using RCPsiSquared.Core.Numerics;
using M = RCPsiSquared.Core.Numerics.GaussianRationalMatrix;

namespace RCPsiSquared.Core.Quadrupolar;

/// <summary>Composition of generators, the hook for nuclei and pools: the Kronecker sum
/// L₁ ⊗ 1 + 1 ⊗ L₂ written on the product space's own row-major vec, index (d₁d₂)·i + j with
/// i = d₂·i₁ + i₂ and j = d₂·j₁ + j₂, so the result is a generator like any other.</summary>
public static class GeneratorComposition
{
    public static GaussianRational[,] KroneckerSum(GaussianRational[,] l1, int d1, GaussianRational[,] l2, int d2)
    {
        int d = d1 * d2;
        var l = M.Zero(d * d);
        int Cell(int i1, int i2, int j1, int j2) => d * (d2 * i1 + i2) + (d2 * j1 + j2);
        for (int r = 0; r < d1 * d1; r++)
            for (int c = 0; c < d1 * d1; c++)
            {
                var v = l1[r, c];
                if (v.IsZero) continue;
                for (int i2 = 0; i2 < d2; i2++)
                    for (int j2 = 0; j2 < d2; j2++)
                        l[Cell(r / d1, i2, r % d1, j2), Cell(c / d1, i2, c % d1, j2)] += v;
            }
        for (int r = 0; r < d2 * d2; r++)
            for (int c = 0; c < d2 * d2; c++)
            {
                var v = l2[r, c];
                if (v.IsZero) continue;
                for (int i1 = 0; i1 < d1; i1++)
                    for (int j1 = 0; j1 < d1; j1++)
                        l[Cell(i1, r / d2, j1, r % d2), Cell(i1, c / d2, j1, c % d2)] += v;
            }
        return l;
    }
}
