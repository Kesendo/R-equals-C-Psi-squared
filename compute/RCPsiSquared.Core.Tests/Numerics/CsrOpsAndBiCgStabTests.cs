using System;
using System.Numerics;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.Numerics;

public class CsrOpsAndBiCgStabTests
{
    private static WeightCoherenceSectorCsr.Csr SmallBlock()
        => WeightCoherenceSectorCsr.BuildFull(5, 2, 3, new Complex(0.7, 0.0));

    [Fact]
    public void Multiply_MatchesDenseMultiply()
    {
        var m = SmallBlock();
        var dense = WeightCoherenceBlock.Build(5, 2, 3, new Complex(0.7, 0.0));
        var rng = new Random(7);
        var x = new Complex[m.Dim];
        for (int i = 0; i < m.Dim; i++) x[i] = new Complex(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5);
        var y = new Complex[m.Dim];
        CsrOps.Multiply(m, x, y);
        for (int r = 0; r < m.Dim; r++)
        {
            Complex expect = Complex.Zero;
            for (int c = 0; c < m.Dim; c++) expect += dense[r, c] * x[c];
            Assert.True((expect - y[r]).Magnitude <= 1e-12 * (1 + expect.Magnitude));
        }
    }

    [Fact]
    public void HermitianTranspose_IsExactConjTranspose()
    {
        var m = SmallBlock();
        var h = CsrOps.HermitianTranspose(m);
        var back = CsrOps.HermitianTranspose(h);
        Assert.Equal(m.RowPtr, back.RowPtr);
        Assert.Equal(m.ColIdx, back.ColIdx);
        for (int k = 0; k < m.Values.Length; k++)
            Assert.True((m.Values[k] - back.Values[k]).Magnitude == 0.0);
    }

    [Fact]
    public void HermitianTranspose_MatchesHandComputed_On3x3()
    {
        // Known asymmetric 3x3 with complex entries. The double-transpose test alone would pass for
        // the identity function; this pins the transpose against a hand-computed conjugate transpose.
        //   M = [ 1+2i    0    3-1i ]
        //       [   0    4i     0   ]
        //       [   5   6+7i  -2i   ]
        var rowPtr = new[] { 0, 2, 3, 6 };
        var colIdx = new[] { 0, 2, 1, 0, 1, 2 };
        var values = new[]
        {
            new Complex(1, 2), new Complex(3, -1),
            new Complex(0, 4),
            new Complex(5, 0), new Complex(6, 7), new Complex(0, -2),
        };
        var m = new WeightCoherenceSectorCsr.Csr(3, rowPtr, colIdx, values);

        var h = CsrOps.HermitianTranspose(m);

        // H[i,j] = conj(M[j,i]):
        //   H = [ 1-2i    0     5   ]
        //       [   0   -4i   6-7i  ]
        //       [ 3+1i    0    2i   ]
        var expectedRowPtr = new[] { 0, 2, 4, 6 };
        var expectedColIdx = new[] { 0, 2, 1, 2, 0, 2 };
        var expectedValues = new[]
        {
            new Complex(1, -2), new Complex(5, 0),
            new Complex(0, -4), new Complex(6, -7),
            new Complex(3, 1), new Complex(0, 2),
        };

        Assert.Equal(m.Dim, h.Dim);
        Assert.Equal(expectedRowPtr, h.RowPtr);
        Assert.Equal(expectedColIdx, h.ColIdx);
        for (int k = 0; k < expectedValues.Length; k++)
            Assert.True((expectedValues[k] - h.Values[k]).Magnitude == 0.0,
                $"entry {k}: expected {expectedValues[k]} got {h.Values[k]}");
    }

    [Fact]
    public void BiCgStab_SolvesShiftedSystem_AtModerateConditioning()
    {
        var m = SmallBlock();
        var shift = new Complex(-3.0, 0.4);       // generic: far from the spectrum's tight spots
        var rng = new Random(11);
        var rhs = new Complex[m.Dim];
        for (int i = 0; i < m.Dim; i++) rhs[i] = new Complex(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5);
        var x = new Complex[m.Dim];
        var outcome = BiCgStabSolver.Solve(m, shift, rhs, x, relTol: 1e-10, maxIter: 20000);
        Assert.True(outcome.Converged, $"stalled at relres {outcome.RelResidual}");
        var r = new Complex[m.Dim];
        CsrOps.MultiplyShifted(m, shift, x, r);
        double num = 0, den = 0;
        for (int i = 0; i < m.Dim; i++) { num += MagSq(r[i] - rhs[i]); den += MagSq(rhs[i]); }
        Assert.True(Math.Sqrt(num / den) <= 1e-8, "true residual must match the reported one");
    }

    private static double MagSq(Complex z) => z.Real * z.Real + z.Imaginary * z.Imaginary;
}
