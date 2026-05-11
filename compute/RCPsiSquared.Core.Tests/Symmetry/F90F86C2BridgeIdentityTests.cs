using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F90F86C2BridgeIdentityTests
{
    private static F90F86C2BridgeIdentity BuildClaim()
    {
        var f89 = new F89TopologyOrbitClosure(new Pi2DyadicLadderClaim());
        var atLock = new F89PathKAtLockMechanismClaim(f89);
        return new F90F86C2BridgeIdentity(f89, atLock);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void JConventionFactor_IsTwo()
    {
        Assert.Equal(2.0, F90F86C2BridgeIdentity.JConventionFactor, precision: 14);
    }

    [Theory]
    [InlineData(2.5470, 1.2735)]   // F86 Endpoint Q_peak (N=6) → F89 Q
    [InlineData(1.6157, 0.80785)]  // F86 Interior flanking Q_peak (N=6 b=1)
    [InlineData(16.7895, 8.39475)] // F86 N=8 b=3 central escape
    public void F86JToF89J_HalvesValue(double f86J, double expectedF89J)
    {
        Assert.Equal(expectedF89J, F90F86C2BridgeIdentity.F86JToF89J(f86J), precision: 12);
    }

    [Theory]
    [InlineData(1.2735, 2.5470)]
    [InlineData(0.80785, 1.6157)]
    [InlineData(8.39475, 16.7895)]
    public void F89JToF86J_DoublesValue(double f89J, double expectedF86J)
    {
        Assert.Equal(expectedF86J, F90F86C2BridgeIdentity.F89JToF86J(f89J), precision: 12);
    }

    [Fact]
    public void F86JToF89J_NegativeJ_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F90F86C2BridgeIdentity.F86JToF89J(-0.1));
    }

    [Fact]
    public void F89JToF86J_NegativeJ_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => F90F86C2BridgeIdentity.F89JToF86J(-0.1));
    }

    [Fact]
    public void RoundTrip_F86ToF89ToF86_IsIdentity()
    {
        double[] testValues = { 0.0, 0.5, 1.0, 2.5470, 16.7895 };
        foreach (var v in testValues)
            Assert.Equal(v, F90F86C2BridgeIdentity.F89JToF86J(F90F86C2BridgeIdentity.F86JToF89J(v)), precision: 14);
    }

    [Fact]
    public void BitExactBondCount_IsAtLeast27Of29()
    {
        Assert.True(F90F86C2BridgeIdentity.BitExactBondCountVerified >= 27);
        Assert.True(F90F86C2BridgeIdentity.TotalBondComparisonsVerified >= F90F86C2BridgeIdentity.BitExactBondCountVerified);
    }
}
