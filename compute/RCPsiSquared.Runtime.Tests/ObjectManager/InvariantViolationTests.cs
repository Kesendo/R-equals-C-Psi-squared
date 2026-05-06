using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.Tests.TestSupport;

namespace RCPsiSquared.Runtime.Tests.ObjectManager;

public class InvariantViolationTests
{
    private sealed class NeedsBaz : Claim
    {
        public Claim Baz { get; }
        public NeedsBaz(Claim baz) : base("synthetic NeedsBaz", Tier.Tier1Derived,
            "compute/RCPsiSquared.Runtime.Tests/TestSupport/TestClaims.cs") { Baz = baz; }
        public override string DisplayName => "NeedsBaz";
        public override string Summary => "synthetic claim with unresolved Baz dependency";
    }

    private sealed class Baz : Claim
    {
        public Baz() : base("Baz", Tier.Tier1Derived, "x") { }
        public override string DisplayName => "Baz";
        public override string Summary => "synthetic Baz";
    }

    [Fact]
    public void Build_MissingParent_Throws_WithDependencyName()
    {
        var ex = Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .Register<NeedsBaz>(b => new NeedsBaz(b.Get<Baz>()))
                .Build());

        Assert.Equal("MissingParent", ex.Rule);
        Assert.Contains("Baz", ex.Message);
        Assert.Contains("NeedsBaz", ex.Message);
    }
}
