using System;
using System.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The per-site Z-field overload of the (SE,DE) block, for the F89 Door-C CSR sweep Stage 2
/// (random-field disorder ensemble). A field Σ_k w_k Z_k is diagonal and Hermitian, so like the Δ·ZZ
/// term it leaves the absorption-theorem real rate untouched and only shifts the imaginary frequency by
/// −i·q·(fieldEnergy(ket) − fieldEnergy(bra)). A random field breaks integrability AND the S₂ reflection
/// AND conjugation symmetry (so Stage 2 uses BuildFull, OffReal domain).</summary>
public class XxzCoherenceBlockFieldTests
{
    /// <summary>Consistency anchor: a zero field reproduces BuildFull exactly.</summary>
    [Fact]
    public void BuildFullWithField_ZeroField_EqualsBuildFull()
    {
        int n = 5; var q = new Complex(2, 0); double delta = 0.5;
        var bare = XxzCoherenceBlock.BuildFull(n, q, delta);
        var zero = XxzCoherenceBlock.BuildFullWithField(n, q, delta, new double[n]);
        int d = bare.GetLength(0);
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                Assert.Equal(bare[i, j], zero[i, j]);
    }

    /// <summary>A nonzero field shifts the IMAGINARY diagonal (the frequency) but leaves the AT real rate
    /// and every off-diagonal entry untouched.</summary>
    [Fact]
    public void BuildFullWithField_ShiftsImaginaryDiagonal_LeavesRestUnchanged()
    {
        int n = 5; var q = new Complex(2, 0); double delta = 0.0;
        var bare = XxzCoherenceBlock.BuildFull(n, q, delta);
        var w = new[] { 0.3, -0.5, 0.1, 0.2, -0.4 };
        var fielded = XxzCoherenceBlock.BuildFullWithField(n, q, delta, w);
        int d = bare.GetLength(0);
        int shifted = 0;
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
            {
                if (i != j) Assert.Equal(bare[i, j], fielded[i, j]);                      // off-diagonal unchanged
                else
                {
                    Assert.Equal(bare[i, i].Real, fielded[i, i].Real, 12);                // AT real rate untouched
                    if (Math.Abs(fielded[i, i].Imaginary - bare[i, i].Imaginary) > 1e-9) shifted++;
                }
            }
        Assert.True(shifted > 0, "a nonzero field must shift some imaginary diagonal entries");
    }
}
