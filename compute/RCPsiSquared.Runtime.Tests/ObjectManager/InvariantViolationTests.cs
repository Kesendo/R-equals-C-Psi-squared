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

    private sealed class CycleA : Claim
    {
        public Claim Other { get; }
        public CycleA(Claim other) : base("CycleA", Tier.Tier1Derived, "x") { Other = other; }
        public override string DisplayName => "CycleA";
        public override string Summary => "synthetic A in A→B→A cycle";
    }

    private sealed class CycleB : Claim
    {
        public Claim Other { get; }
        public CycleB(Claim other) : base("CycleB", Tier.Tier1Derived, "x") { Other = other; }
        public override string DisplayName => "CycleB";
        public override string Summary => "synthetic B in A→B→A cycle";
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

    [Fact]
    public void Build_Cycle_Throws_WithCyclePath()
    {
        var ex = Assert.Throws<InvariantViolationException>(() =>
            new ClaimRegistryBuilder()
                .Register<CycleA>(b => new CycleA(b.Get<CycleB>()))
                .Register<CycleB>(b => new CycleB(b.Get<CycleA>()))
                .Build());

        Assert.Equal("Cycle", ex.Rule);
        Assert.Contains(typeof(CycleA), ex.Path);
        Assert.Contains(typeof(CycleB), ex.Path);
    }
}
