using System.Linq;
using System.Reflection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F87;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>Parallel architecture-enforcement test for the Diagnostics assembly: scope-extends
/// the IZ2AxisClaim requirement so every <c>*Pi2Inheritance</c> Claim that lives in
/// Diagnostics (currently <see cref="F87Pi2Inheritance"/>) also carries a Z₂-axis tag.
///
/// <para>The Core-side equivalent
/// (<c>RCPsiSquared.Core.Tests.Symmetry.Pi2InheritanceZ2AxisTests.EveryPi2InheritanceClaim_ImplementsIZ2AxisClaim</c>)
/// scans the Core assembly only and would miss F87Pi2Inheritance because Diagnostics is a
/// downstream assembly. This test closes that gap.</para></summary>
public class F87Pi2InheritanceZ2AxisTests
{
    private static readonly Assembly DiagnosticsAssembly = typeof(F87Pi2Inheritance).Assembly;

    [Fact]
    public void EveryPi2InheritanceClaim_InDiagnostics_ImplementsIZ2AxisClaim()
    {
        var pi2InheritanceClaims = DiagnosticsAssembly.GetTypes()
            .Where(t => t.IsClass && !t.IsAbstract)
            .Where(t => typeof(Claim).IsAssignableFrom(t))
            .Where(t => t.Name.EndsWith("Pi2Inheritance"))
            .ToList();

        Assert.NotEmpty(pi2InheritanceClaims);

        var missing = pi2InheritanceClaims
            .Where(t => !typeof(IZ2AxisClaim).IsAssignableFrom(t))
            .Select(t => t.Name)
            .OrderBy(n => n, StringComparer.Ordinal)
            .ToList();

        Assert.True(missing.Count == 0,
            "Pi²-Inheritance Claims in Diagnostics that do not implement IZ2AxisClaim: "
            + string.Join(", ", missing));
    }

    [Fact]
    public void F87Pi2Inheritance_HasZ2AxisBitB()
    {
        Assert.True(typeof(IZ2AxisClaim).IsAssignableFrom(typeof(F87Pi2Inheritance)));

        var ladder = new Pi2DyadicLadderClaim();
        var memoryLoop = new Pi2I4MemoryLoopClaim();
        var f1 = new F1Pi2Inheritance(new RCPsiSquared.Core.F1.F1PalindromeIdentity(), ladder, memoryLoop);
        var f87 = new F87Pi2Inheritance(f1);

        Assert.Equal(Z2Axis.BitB, f87.Z2Axis);
        Assert.Null(f87.BitATwin);
    }
}
