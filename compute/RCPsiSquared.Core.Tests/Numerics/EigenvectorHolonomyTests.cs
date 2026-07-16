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
    public void DisguisedPiFreeEP2_SameMod4Holonomy_TheThirdClockIsGeneric()
    {
        // The third-clock control (decides SeedHolonomyClaim's "noted correspondence"): a defective EP2
        // with NO palindrome, chain, or Clifford structure carries the identical ±i / M₂=−I / M₄=+I
        // holonomy, so the seed's mod-4 loop cannot be Π's Z₄ (a Π-free system shows it) — it is generic
        // EP2 geometry. Construction: B(s) = blockdiag([[s,1],[1,−s]], S₄) has its EP2 at s = i EXACTLY;
        // conjugating by a complex-orthogonal Q (QᵀQ = I, complex Givens rotations) preserves complex
        // symmetry and defectiveness while mixing the toy across all six dimensions — unlike the bare
        // 2×2 toy above, the coalescing pair here lives in no preferred basis and has spectators, the
        // regime where the vᵀv tracking is actually exercised. Owner: experiments/SEED_HOLONOMY_THIRD_CLOCK.md.
        const int n = 6;
        var s4 = RandomComplexSymmetricShifted(seedState: 20260716);
        var q = ComplexOrthogonal(n, rotations: 12, seedState: 424242);

        Complex[,] Mat(Complex s)
        {
            var b = new Complex[n, n];
            b[0, 0] = s; b[0, 1] = Complex.One; b[1, 0] = Complex.One; b[1, 1] = -s;
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    b[i + 2, j + 2] = s4[i, j];
            return Congruence(q, b, n);   // Q · B · Qᵀ, complex symmetric by construction
        }

        // sanity: complex symmetry survived the mixing, and the spectators sit away from λ* = 0
        var mStar = Mat(Complex.ImaginaryOne);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                Assert.True((mStar[i, j] - mStar[j, i]).Magnitude < 1e-12, "disguised matrix not complex symmetric");

        var res = EigenvectorHolonomy.FrameMonodromy(
            Mat, dim: n, lambda0: Complex.Zero,
            center: Complex.ImaginaryOne, radius: 0.05, nLoops: 4, stepsPerLoop: 1200);

        for (int k = 0; k < res.SpanResidual.Length; k++)
            Assert.True(res.SpanResidual[k] < 1e-6, $"loop {k + 1} span residual {res.SpanResidual[k]:e2} too large");

        var ev = res.M1Eigenvalues;
        Assert.True((ev[0] + ev[1]).Magnitude < 1e-3, $"tr M₁ = {ev[0] + ev[1]} not ~0 (±i)");
        Assert.True((ev[0] * ev[1] - Complex.One).Magnitude < 1e-3, $"det M₁ = {ev[0] * ev[1]} not ~1");
        Assert.True(DistFromScalarI(res.LoopMonodromy[1], -1.0) < 1e-3, "M₂ not ~ −I");
        Assert.True(DistFromScalarI(res.LoopMonodromy[3], +1.0) < 1e-3, "M₄ not ~ +I");
    }

    // deterministic LCG (no System.Random: its sequence is an implementation detail across .NET versions)
    static double NextUniform(ref ulong state)
    {
        state = state * 6364136223846793005UL + 1442695040888963407UL;
        return (state >> 11) * (1.0 / (1UL << 53)) * 2.0 - 1.0;
    }

    // random complex-symmetric 4×4, spectrum shifted by +3 so the spectators clear the EP eigenvalue 0
    static Complex[,] RandomComplexSymmetricShifted(ulong seedState)
    {
        var m = new Complex[4, 4];
        for (int i = 0; i < 4; i++)
            for (int j = i; j < 4; j++)
            {
                var z = new Complex(NextUniform(ref seedState), NextUniform(ref seedState));
                m[i, j] = z; m[j, i] = z;
            }
        for (int i = 0; i < 4; i++) m[i, i] += 3.0;
        return m;
    }

    // complex-orthogonal Q = product of complex Givens rotations: QᵀQ = I (cos²θ + sin²θ = 1 for complex θ)
    static Complex[,] ComplexOrthogonal(int n, int rotations, ulong seedState)
    {
        var q = new Complex[n, n];
        for (int i = 0; i < n; i++) q[i, i] = Complex.One;
        for (int r = 0; r < rotations; r++)
        {
            int i = (int)(Math.Abs(NextUniform(ref seedState)) * n) % n;
            int j = (int)(Math.Abs(NextUniform(ref seedState)) * n) % n;
            if (i == j) { r--; continue; }
            var theta = new Complex(NextUniform(ref seedState), 0.4 * NextUniform(ref seedState));
            Complex c = Complex.Cos(theta), s = Complex.Sin(theta);
            for (int row = 0; row < n; row++)
            {
                Complex qi = q[row, i], qj = q[row, j];
                q[row, i] = qi * c - qj * s;
                q[row, j] = qi * s + qj * c;
            }
        }
        return q;
    }

    // Q · B · Qᵀ
    static Complex[,] Congruence(Complex[,] q, Complex[,] b, int n)
    {
        var qb = new Complex[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                Complex acc = Complex.Zero;
                for (int k = 0; k < n; k++) acc += q[i, k] * b[k, j];
                qb[i, j] = acc;
            }
        var outM = new Complex[n, n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                Complex acc = Complex.Zero;
                for (int k = 0; k < n; k++) acc += qb[i, k] * q[j, k];
                outM[i, j] = acc;
            }
        return outM;
    }

    [Fact]
    public void WeightCoherenceBlock_IsComplexSymmetric_TheVTvGaugePrecondition()
    {
        // The vᵀv gauge presupposes Lᵀ = L. Established in prose (F89_BETA_EXOTIC_GENERICITY: diagonal real,
        // hop C = iK with K real symmetric) but never pinned from below until now. N=5 and N=9 (1,2)-blocks.
        foreach (int n in new[] { 5, 9 })
        {
            var m = WeightCoherenceBlock.Build(n, 1, 2, new Complex(0.7, 0.3));
            int d = m.GetLength(0);
            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                    Assert.True(m[i, j] == m[j, i], $"N={n}: L[{i},{j}] != L[{j},{i}] — block not complex symmetric");
        }
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
