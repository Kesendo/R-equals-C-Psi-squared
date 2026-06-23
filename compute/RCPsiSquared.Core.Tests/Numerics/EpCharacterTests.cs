using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Numerics;
using Xunit;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Numerics;

/// <summary>Gate-first validation of the artifact-free EP-character diagnostic on KNOWN answers:
/// a defective toy Jordan 2×2 MUST read DEFECTIVE; a diabolic repeated eigenvalue MUST read
/// DIABOLIC. Mirrors GATE 0 of <c>simulations/_review_coherence_horizon_ep.py</c> and the toy gate
/// of <c>simulations/review_f86a_diabolic_vs_defective.py</c>. None of the three measures reads an
/// <c>eig</c> eigenvector pairing, so a verdict here cannot be the F86a-misfire family.</summary>
public class EpCharacterTests
{
    private static ComplexMatrix Mat(Complex[,] a) => Matrix<Complex>.Build.DenseOfArray(a);

    // ---- GATE 0: KNOWN-ANSWER validation ----

    [Fact]
    public void Gate0a_ToyJordanBlock_ReadsDefective()
    {
        // [[-2g, i*w],[0, -2g]]: a genuine non-trivial Jordan block (double eigenvalue −2g, geo mult 1).
        double g = 1.0, w = 1.3;
        var l = Mat(new Complex[,]
        {
            { new Complex(-2 * g, 0), new Complex(0, w) },
            { Complex.Zero,           new Complex(-2 * g, 0) },
        });
        var r = EpCharacter.Characterize(l, new Complex(-2 * g, 0), radius: 0.5);

        Assert.Equal(EpCharacter.EpKind.Defective, r.Kind);
        Assert.Equal(2, r.Algebraic);
        Assert.Equal(1, r.Geometric);             // geo < alg ⟹ Jordan
        Assert.Equal(1.3, r.Departure, 3);        // dep = |w| = 1.3 (Python GATE 0(a))
        Assert.Equal(1.0, r.ProjectorNorm, 4);    // closed 2×2 ⟹ ‖P‖ ≈ 1 (supplementary, not gating)
    }

    [Fact]
    public void Gate0b_ToyLeffEp_ReadsDefective_WithDepartureFour()
    {
        // [[-2g, i*J],[i*J, -6g]] at J=2g: the F86a toy_Leff EP, double root at −4g, Jordan.
        double g = 1.0, j = 2.0 * g;
        var l = Mat(new Complex[,]
        {
            { new Complex(-2 * g, 0), new Complex(0, j) },
            { new Complex(0, j),      new Complex(-6 * g, 0) },
        });
        var r = EpCharacter.Characterize(l, new Complex(-4 * g, 0), radius: 0.5);

        Assert.Equal(EpCharacter.EpKind.Defective, r.Kind);
        Assert.Equal(2, r.Algebraic);
        Assert.Equal(1, r.Geometric);
        Assert.Equal(4.0, r.Departure, 3);        // dep = 4.0 exactly (Python GATE 0(b))
    }

    [Fact]
    public void Gate0c_DiagonalRepeated_ReadsDiabolic()
    {
        // diag(−2g, −2g): two independent equal eigenvalues, the canonical diabolic point.
        double g = 1.0;
        var l = Mat(new Complex[,]
        {
            { new Complex(-2 * g, 0), Complex.Zero },
            { Complex.Zero,           new Complex(-2 * g, 0) },
        });
        var r = EpCharacter.Characterize(l, new Complex(-2 * g, 0), radius: 0.5);

        Assert.Equal(EpCharacter.EpKind.Diabolic, r.Kind);
        Assert.Equal(2, r.Algebraic);
        Assert.Equal(2, r.Geometric);             // geo == alg ⟹ diabolic
        Assert.True(r.Departure < 1e-6, $"diabolic departure {r.Departure} should be ≈ 0");
    }

    [Fact]
    public void Gate0d_RotatedDiagonal_StillReadsDiabolic()
    {
        // A unitary rotation of diag(−2g,−2g) is still normal with a repeated eigenvalue ⟹ diabolic.
        double g = 1.0, th = 0.7;
        var q = Mat(new Complex[,]
        {
            { new Complex(System.Math.Cos(th), 0), new Complex(-System.Math.Sin(th), 0) },
            { new Complex(System.Math.Sin(th), 0), new Complex(System.Math.Cos(th), 0) },
        });
        var diag = Mat(new Complex[,]
        {
            { new Complex(-2 * g, 0), Complex.Zero },
            { Complex.Zero,           new Complex(-2 * g, 0) },
        });
        var l = q * diag * q.ConjugateTranspose();
        var r = EpCharacter.Characterize(l, new Complex(-2 * g, 0), radius: 0.5);

        Assert.Equal(EpCharacter.EpKind.Diabolic, r.Kind);
        Assert.Equal(2, r.Geometric);
        Assert.True(r.Departure < 1e-6, $"rotated-diabolic departure {r.Departure} should be ≈ 0");
    }

    [Fact]
    public void Gate0e_ObliquelyEmbeddedJordan_StillReadsDefective()
    {
        // A 2×2 Jordan block direct-summed with distinct eigenvalues, then conjugated by a
        // fixed non-orthogonal S so the eigenspace is OBLIQUE / non-normal (‖P‖>1).
        // The dep/geo signal must NOT be fooled (the failure mode that sank F86a's eig instruments).
        var lam = new Complex(-4, 1.3);
        var jb = Mat(new Complex[,]
        {
            { lam, Complex.One, 0, 0, 0, 0 },
            { 0, lam, 0, 0, 0, 0 },
            { 0, 0, new Complex(1,0), 0, 0, 0 },
            { 0, 0, 0, new Complex(2,0), 0, 0 },
            { 0, 0, 0, 0, new Complex(3,0), 0 },
            { 0, 0, 0, 0, 0, new Complex(-1,0) },
        });
        var s = Mat(FixedObliqueSimilarity());
        var l = s * jb * s.Inverse();
        var r = EpCharacter.Characterize(l, lam, radius: 0.5);
        Assert.Equal(EpCharacter.EpKind.Defective, r.Kind);   // load-bearing: still defective despite obliqueness
        Assert.Equal(1, r.Geometric);                         // load-bearing: geo=1 (Jordan)
        Assert.Equal(2, r.Algebraic);
        Assert.True(r.ProjectorNorm > 1.5, $"oblique embedding should give ‖P‖>1.5, got {r.ProjectorNorm}");
    }

    [Fact]
    public void Gate0f_F89ModelAtGammaZero_ReadsDiabolic_NotRigged()
    {
        // At γ=0 the (SE,DE) block is anti-Hermitian (normal) ⟹ all degeneracies are diabolic.
        // If the instrument read defective here, it would be rigged.
        var l = RCPsiSquared.Core.F89PathK.F89Path3OcticBlock.BuildSeDeSymBlock(0.69, 0.0);
        var (center, _) = NearestPairMidpoint(l);
        var r = EpCharacter.Characterize(l, center, radius: 0.3);
        Assert.NotEqual(EpCharacter.EpKind.Defective, r.Kind);
    }

    // ---- the eigenvector-merge corroborator ----

    [Fact]
    public void EigenvectorMerge_JordanPairCosineApproachesOne()
    {
        // A near-defective 2×2 [[λ, 1],[ε, λ]]: the two eigenvectors are nearly parallel ⟹ |cos| → 1.
        var l = Mat(new Complex[,]
        {
            { Complex.One, Complex.One },
            { new Complex(1e-6, 0), Complex.One },
        });
        double cos = EpCharacter.EigenvectorMergeCosine(l);
        Assert.True(cos > 0.999, $"near-Jordan eigenvector |cos| {cos} should approach 1");
    }

    [Fact]
    public void EigenvectorMerge_DiabolicPairStaysOrthogonalish()
    {
        // diag(1, 2): distinct normal eigenvectors are orthonormal ⟹ |cos| ≈ 0.
        var l = Mat(new Complex[,]
        {
            { Complex.One, Complex.Zero },
            { Complex.Zero, new Complex(2, 0) },
        });
        double cos = EpCharacter.EigenvectorMergeCosine(l);
        Assert.True(cos < 1e-9, $"distinct-eigenvector |cos| {cos} should be ≈ 0");
    }

    // ---- the building blocks in isolation ----

    [Fact]
    public void RieszProjector_EnclosingOneEigenvalue_HasTraceOne()
    {
        // diag(0, −5): a circle of radius 1 about 0 encloses exactly one eigenvalue ⟹ tr P = 1.
        var l = Mat(new Complex[,]
        {
            { Complex.Zero, Complex.Zero },
            { Complex.Zero, new Complex(-5, 0) },
        });
        var p = EpCharacter.RieszProjector(l, Complex.Zero, radius: 1.0);
        Assert.Equal(1.0, p.Trace().Real, 4);
    }

    [Fact]
    public void DepartureFromNormality_IsZeroForNormal_PositiveForJordan()
    {
        var normal = Mat(new Complex[,]
        {
            { new Complex(1, 0), Complex.Zero },
            { Complex.Zero, new Complex(2, 0) },
        });
        // The departure is limited by Evd precision (the eigenvalue magnitudes carry ~1e-7 error), so a
        // normal matrix reads ≈ 0 to ~1e-6, not to machine ε — exactly why the verdict uses a RELATIVE
        // departure tolerance (1e-2), never an absolute < 1e-9.
        Assert.True(EpCharacter.DepartureFromNormality(normal) < 1e-6,
            "a normal matrix should have departure-from-normality ≈ 0");

        var jordan = Mat(new Complex[,]
        {
            { Complex.One, new Complex(3, 0) },
            { Complex.Zero, Complex.One },
        });
        Assert.Equal(3.0, EpCharacter.DepartureFromNormality(jordan), 6);
    }

    // ---- deterministic helpers for the oblique-embedding / γ=0 gates ----

    private static Complex[,] FixedObliqueSimilarity()
    {
        // I + a fixed off-diagonal perturbation (deterministic, invertible: det = 1).
        // Magnitudes calibrated so the conjugated Jordan eigenspace is oblique enough that ‖P‖ ≈ 2.2 > 1.5
        // (the original 0.7/-0.5/… set was only mildly oblique at ‖P‖ ≈ 1.40); the verdict stays Defective/geo=1.
        var a = new Complex[6, 6];
        for (int i = 0; i < 6; i++) a[i, i] = Complex.One;
        a[0, 2] = 1.1; a[1, 4] = -0.8; a[3, 0] = 0.6; a[5, 1] = 0.9; a[2, 5] = -0.5;
        return a;
    }

    private static (Complex Center, double Dist) NearestPairMidpoint(
        MathNet.Numerics.LinearAlgebra.Matrix<Complex> l)
    {
        var ev = l.Evd().EigenValues.ToArray();
        double best = double.MaxValue; Complex mid = Complex.Zero;
        for (int i = 0; i < ev.Length; i++)
            for (int k = i + 1; k < ev.Length; k++)
            {
                double d = (ev[i] - ev[k]).Magnitude;
                if (d < best) { best = d; mid = 0.5 * (ev[i] + ev[k]); }
            }
        return (mid, best);
    }
}
