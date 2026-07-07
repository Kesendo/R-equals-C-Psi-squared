using System;
using System.Numerics;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.Numerics;

/// <summary>The eigenVECTOR-frame holonomy around a defective exceptional point (the companion of
/// <see cref="Monodromy"/>'s eigenVALUE swap). Encircling a defective EP2, the two coalescing right
/// eigenvectors — transported in the biorthogonal vᵀv gauge natural to a complex-symmetric matrix —
/// rotate by a generator with eigenvalues ±i, so the frame is single-valued only after 4 loops
/// (M₁ ~ 90° rotation, M₂ = −I, M₄ = I: i⁴ = 1). A Hermitian vᴴv gauge sees only the mod-2 swap
/// (eigenvalues ±1, M₂ = +I); the vᵀv self-orthogonality is what promotes it to the mod-4 (±i) frame —
/// the load-bearing piece (the static |cos| merge of EpCharacter is a different, non-holonomy object).</summary>
public class EigenvectorHolonomyTests
{
    [Fact]
    public void ComplexSymmetricEP2_FrameRotatesBy90Degrees_i4Equals1()
    {
        // M(s) = [[s,1],[1,−s]], complex-symmetric, defective EP2 at s=i (eigenvalues ±√(s²+1), λ*=0),
        // self-orthogonal eigenvector [1,−i]. Encircle s=i.
        Complex[,] Mat(Complex s) => new Complex[,] { { s, Complex.One }, { Complex.One, -s } };
        var res = EigenvectorHolonomy.FrameMonodromy(
            Mat, dim: 2, lambda0: Complex.Zero,
            center: Complex.ImaginaryOne, radius: 0.1, nLoops: 4, stepsPerLoop: 2000);

        // the 2D span is preserved on EVERY loop (odd and even) — the clean-EP signature
        for (int k = 0; k < res.SpanResidual.Length; k++)
            Assert.True(res.SpanResidual[k] < 1e-6, $"loop {k + 1} span residual {res.SpanResidual[k]:e2} too large");

        // M₁ has eigenvalues ±i: trace 0, det 1
        var ev = res.M1Eigenvalues;
        Assert.True((ev[0] + ev[1]).Magnitude < 1e-3, $"tr M₁ = {ev[0] + ev[1]} not ~0");
        Assert.True((ev[0] * ev[1] - Complex.One).Magnitude < 1e-3, $"det M₁ = {ev[0] * ev[1]} not ~1");

        // M₂ = −I  and  M₄ = +I  (the mod-4 memory loop)
        Assert.True(DistFromScalarI(res.LoopMonodromy[1], -1.0) < 1e-3, "M₂ not ~ −I");
        Assert.True(DistFromScalarI(res.LoopMonodromy[3], +1.0) < 1e-3, "M₄ not ~ +I");
    }

    [Fact]
    public void DiabolicCrossing_IsNotTheMod4Loop_FrameDoesNotRotate()
    {
        // Specificity: M(s)=diag(s,−s) — eigenvalues ±s CROSS at s=0 with INDEPENDENT eigenvectors (a
        // DIABOLIC crossing, not a defective EP). No √-branch ⟹ no eigenvalue swap ⟹ the frame does not
        // rotate: M₁ ≈ I (tr ≈ 2), not ±i. The witness cannot mistake a diabolic/spectator for the mod-4 loop.
        Complex[,] Mat(Complex s) => new Complex[,] { { s, Complex.Zero }, { Complex.Zero, -s } };
        var res = EigenvectorHolonomy.FrameMonodromy(
            Mat, dim: 2, lambda0: Complex.Zero, center: Complex.Zero, radius: 0.1, nLoops: 4, stepsPerLoop: 800);
        var ev = res.M1Eigenvalues;
        Assert.True((ev[0] + ev[1]).Magnitude > 1.0,
            $"tr M₁ = {ev[0] + ev[1]} — a diabolic crossing must NOT produce the ±i swap (expect tr≈2, M₁≈I)");
    }

    [Fact]
    public void N5LivingSeed_EigenvectorFrameIsTheMod4MemoryLoop_i4Equals1()
    {
        // The N=5 (1,2)-block defective seed: octic q* ≈ 0.643037, λ* ≈ −3.8196 (WeightCoherenceBlock's q
        // is the octic axis, hop ±2iq). The full 50-dim block tracks cleanly at N=5 (no sector isolation
        // needed here). This test is ALSO the gauge guard: replacing BNorm's vᵀv with a Hermitian magnitude
        // drops this to ±1 (mod-2) and fails, while the 2×2 toy above stays ±i either way.
        const int n = 5;
        Complex qstar = 0.643037, lam0 = -3.8196;

        int dim = WeightCoherenceBlock.Build(n, 1, 2, qstar).GetLength(0);   // C(5,1)·C(5,2) = 50
        Complex[,] Mat(Complex q) => WeightCoherenceBlock.Build(n, 1, 2, q);

        var res = EigenvectorHolonomy.FrameMonodromy(
            Mat, dim, lam0, qstar, radius: 0.001, nLoops: 4, stepsPerLoop: 400);

        for (int k = 0; k < res.SpanResidual.Length; k++)
            Assert.True(res.SpanResidual[k] < 1e-4, $"loop {k + 1} span residual {res.SpanResidual[k]:e2} too large");

        var ev = res.M1Eigenvalues;
        Assert.True((ev[0] + ev[1]).Magnitude < 1e-2, $"tr M₁ = {ev[0] + ev[1]} not ~0 (defective ±i)");
        Assert.True((ev[0] * ev[1] - Complex.One).Magnitude < 1e-2, $"det M₁ = {ev[0] * ev[1]} not ~1");
        Assert.True(DistFromScalarI(res.LoopMonodromy[1], -1.0) < 1e-2, "M₂ not ~ −I");
        Assert.True(DistFromScalarI(res.LoopMonodromy[3], +1.0) < 1e-2, "M₄ not ~ +I");
    }

    // Frobenius distance of a row-major 2×2 [a,b,c,d] from c·I
    static double DistFromScalarI(Complex[] m, double c)
    {
        var diff = new[] { m[0] - c, m[1], m[2], m[3] - c };
        double s = 0;
        foreach (var z in diff) s += z.Magnitude * z.Magnitude;
        return Math.Sqrt(s);
    }
}
