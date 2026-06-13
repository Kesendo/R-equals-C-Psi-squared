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
