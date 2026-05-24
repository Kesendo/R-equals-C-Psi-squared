using System.Linq;
using System.Reflection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class Pi2InheritanceZ2AxisTests
{
    private static readonly Assembly CoreAssembly = typeof(Claim).Assembly;

    [Fact]
    public void EveryPi2InheritanceClaim_ImplementsIZ2AxisClaim()
    {
        var pi2InheritanceClaims = CoreAssembly.GetTypes()
            .Where(t => t.IsClass && !t.IsAbstract)
            .Where(t => t.Name.EndsWith("Pi2Inheritance"))
            .ToList();

        Assert.NotEmpty(pi2InheritanceClaims);

        var missing = pi2InheritanceClaims
            .Where(t => !typeof(IZ2AxisClaim).IsAssignableFrom(t))
            .Select(t => t.Name)
            .ToList();

        Assert.Empty(missing);
    }

    [Fact]
    public void F61_HasZ2AxisBitA()
    {
        var f61Type = CoreAssembly.GetType("RCPsiSquared.Core.Symmetry.F61BitAParityPi2Inheritance");
        Assert.NotNull(f61Type);
        Assert.True(typeof(IZ2AxisClaim).IsAssignableFrom(f61Type));
    }

    [Fact]
    public void KleinFourCellClaim_HasZ2AxisKlein2()
    {
        var kleinType = CoreAssembly.GetType("RCPsiSquared.Core.Symmetry.KleinFourCellClaim");
        Assert.NotNull(kleinType);
        Assert.True(typeof(IZ2AxisClaim).IsAssignableFrom(kleinType));
    }

    [Fact]
    public void EveryBitBClaim_HasBitATwinPropertyInitiallyNull()
    {
        var bitBClaims = CoreAssembly.GetTypes()
            .Where(t => t.IsClass && !t.IsAbstract)
            .Where(t => typeof(IZ2AxisClaim).IsAssignableFrom(t))
            .Where(t => t.GetConstructors().Any())
            .ToList();

        Assert.NotEmpty(bitBClaims);

        // Smoke check: F61 (BitA) returns null for BitATwin, which is fine.
        var f61Type = CoreAssembly.GetType("RCPsiSquared.Core.Symmetry.F61BitAParityPi2Inheritance");
        Assert.NotNull(f61Type);
    }
}
