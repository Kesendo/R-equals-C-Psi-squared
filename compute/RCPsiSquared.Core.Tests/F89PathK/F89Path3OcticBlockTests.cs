using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.F89PathK;

public class F89Path3OcticBlockTests
{
    private static double QEp => F89Path3OcticEpClaim.QEp;          // √((−1+√13)/6)

    [Fact]
    public void SymBlock_AtQEp_HasCoalescingPairAt_MinusFourGammaPlusTwoIJ()
    {
        var l = F89Path3OcticBlock.BuildSeDeSymBlock(QEp, 1.0);     // γ=1, J=q_EP
        Assert.Equal(12, l.RowCount);
        var evals = l.Evd().EigenValues.ToArray();

        // the octic = the 8 non-AT-locked modes (rate −Re ∉ {2,6})
        var octic = evals.Where(e => System.Math.Abs(-e.Real - 2.0) > 1e-6
                                  && System.Math.Abs(-e.Real - 6.0) > 1e-6).ToArray();
        Assert.Equal(8, octic.Length);

        // exactly one coalescing pair at λ_EP = −4 + 2i·q_EP
        var lamEp = new Complex(-4.0, 2.0 * QEp);
        int nearEp = octic.Count(e => (e - lamEp).Magnitude < 1e-6);
        Assert.Equal(2, nearEp);                                   // the double root
    }
}
