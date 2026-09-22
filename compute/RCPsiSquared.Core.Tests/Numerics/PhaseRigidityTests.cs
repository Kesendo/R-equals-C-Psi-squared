using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.Numerics;

public class PhaseRigidityTests
{
    [Fact]
    public void SimpleIsolatedModeOfNonNormalMatrix_CanHaveRigidityBelowOne()
    {
        // The third column is the simple lambda=1 mode of a diagonalizable non-normal matrix.
        // Isolation makes the per-mode Petermann value well-defined; it does not make r equal one.
        var rightEigenvectors = Matrix<Complex>.Build.DenseOfArray(new Complex[,]
        {
            { Complex.One, Complex.Zero, new Complex(3.0, 0.0) },
            { Complex.Zero, Complex.One, Complex.Zero },
            { Complex.Zero, Complex.Zero, Complex.One },
        });

        var rigidity = PhaseRigidity.RigiditiesFromRightEigenvectors(rightEigenvectors)[2];

        Assert.Equal(1.0 / Math.Sqrt(10.0), rigidity, 12);
        Assert.True(rigidity < 1.0);
    }

    [Fact]
    public void ExactSemisimpleDegeneracy_PerVectorPetermannValuesDependOnEigenbasis()
    {
        // The same diagonalizable matrix has a two-dimensional lambda=0 eigenspace.
        // Rotating only that eigenspace leaves L and its spectral projector unchanged,
        // but redistributes the individual Petermann factors.  Therefore an exact
        // degeneracy needs an invariant-subspace/Jordan diagnostic, not per-vector K.
        double c = 1.0 / Math.Sqrt(2.0);
        var eigenvalues = Matrix<Complex>.Build.DiagonalOfDiagonalArray(new[]
        {
            Complex.Zero,
            Complex.Zero,
            Complex.One,
        });
        var selected = Matrix<Complex>.Build.DiagonalOfDiagonalArray(new[]
        {
            Complex.One,
            Complex.One,
            Complex.Zero,
        });
        var basisA = Matrix<Complex>.Build.DenseOfArray(new Complex[,]
        {
            { Complex.One, Complex.Zero, new Complex(3.0, 0.0) },
            { Complex.Zero, Complex.One, Complex.Zero },
            { Complex.Zero, Complex.Zero, Complex.One },
        });
        var rotateDegenerateBlock = Matrix<Complex>.Build.DenseOfArray(new Complex[,]
        {
            { new Complex(c, 0.0), new Complex(-c, 0.0), Complex.Zero },
            { new Complex(c, 0.0), new Complex(c, 0.0), Complex.Zero },
            { Complex.Zero, Complex.Zero, Complex.One },
        });
        var basisB = basisA * rotateDegenerateBlock;

        var matrixA = basisA * eigenvalues * basisA.Inverse();
        var matrixB = basisB * eigenvalues * basisB.Inverse();
        Assert.True((matrixA - matrixB).FrobeniusNorm() < 1e-12);
        Assert.True(
            (matrixA * matrixA.ConjugateTranspose() - matrixA.ConjugateTranspose() * matrixA)
                .FrobeniusNorm() > 1.0,
            "the counterexample must remain non-normal");

        var projectorA = basisA * selected * basisA.Inverse();
        var projectorB = basisB * selected * basisB.Inverse();
        Assert.True((projectorA - projectorB).FrobeniusNorm() < 1e-12);

        var rigidityA = PhaseRigidity.RigiditiesFromRightEigenvectors(basisA);
        var rigidityB = PhaseRigidity.RigiditiesFromRightEigenvectors(basisB);
        var petermannA = rigidityA.Take(2).Select(r => 1.0 / (r * r)).ToArray();
        var petermannB = rigidityB.Take(2).Select(r => 1.0 / (r * r)).ToArray();

        Assert.Equal(10.0, petermannA[0], 12);
        Assert.Equal(1.0, petermannA[1], 12);
        Assert.All(petermannB, k => Assert.Equal(5.5, k, 12));
    }

    [Fact]
    public void Hermitian_AllRigiditiesAreOne()
    {
        // A Hermitian matrix is normal: left = right eigenvectors, so r = 1 for every mode.
        var H = Matrix<Complex>.Build.DenseOfArray(new Complex[,]
        {
            { new Complex(2, 0), new Complex(0, 1) },
            { new Complex(0, -1), new Complex(3, 0) },
        });
        var modes = PhaseRigidity.Compute(H);
        Assert.All(modes, m => Assert.Equal(1.0, m.Rigidity, 6));
    }

    [Fact]
    public void NearDefectiveBlock_RigidityCollapsesToZero()
    {
        // A near-defective 2×2 block [[λ,1],[ε,λ]] sits just off the exceptional point [[λ,1],[0,λ]]:
        // its two eigenvectors are nearly parallel, so the phase rigidity collapses toward 0.
        var lam = new Complex(1, 0);
        var J = Matrix<Complex>.Build.DenseOfArray(new Complex[,]
        {
            { lam, Complex.One },
            { new Complex(1e-8, 0), lam },
        });
        var minR = PhaseRigidity.Compute(J).Min(m => m.Rigidity);
        Assert.True(minR < 1e-3, $"near-defective min rigidity {minR} should be ~0 (an approached EP)");
    }

    [Fact]
    public void N2HorizonLiouvillianAtQ1_HasACoalescingGapMode()
    {
        // The N=2 coherence horizon is at Q=1 (γ=J): the {0,2}-coherence gap mode is a genuine EP.
        var L = new ChainSystem(2, 1.0, 1.0).BuildLiouvillian();
        var nz = PhaseRigidity.Compute(L).Where(m => m.Lambda.Real < -1e-6).ToList();
        double gap = nz.Max(m => m.Lambda.Real);
        var gapModes = nz.Where(m => m.Lambda.Real > gap - 0.1).ToList();
        double minR = gapModes.Min(m => m.Rigidity);
        Assert.True(minR < 0.05, $"N=2 Q=1 gap-mode min rigidity {minR:F4} should signal the EP");
    }
}
