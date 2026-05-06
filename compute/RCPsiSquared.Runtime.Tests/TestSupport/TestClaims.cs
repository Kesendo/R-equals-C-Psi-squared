using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Runtime.Tests.TestSupport;

/// <summary>Synthetic claims used across registry tests. They reference no real anchor
/// file; tests that need an anchor-existence check provide their own.</summary>
public sealed class FooT1D : Claim
{
    public FooT1D() : base("synthetic Foo (T1D)", Tier.Tier1Derived,
        "compute/RCPsiSquared.Runtime.Tests/TestSupport/TestClaims.cs") { }
    public override string DisplayName => "Foo";
    public override string Summary => "synthetic Tier1Derived";
}

public sealed class BarT2E : Claim
{
    public BarT2E() : base("synthetic Bar (T2E)", Tier.Tier2Empirical,
        "compute/RCPsiSquared.Runtime.Tests/TestSupport/TestClaims.cs") { }
    public override string DisplayName => "Bar";
    public override string Summary => "synthetic Tier2Empirical";
}
