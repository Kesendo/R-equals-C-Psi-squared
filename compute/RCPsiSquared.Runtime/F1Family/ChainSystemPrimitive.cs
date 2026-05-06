using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Runtime.F1Family;

/// <summary>Schicht 1 wrapper: exposes a <see cref="ChainSystem"/> as a registered
/// <see cref="Claim"/> so F1 / F73 Claims can declare it as a parent. The Claim itself is
/// trivially Tier1Derived (a parameterisation, no hypothesis); the anchor points at the
/// Core record.</summary>
public sealed class ChainSystemPrimitive : Claim
{
    public ChainSystem System { get; }

    public ChainSystemPrimitive(ChainSystem system)
        : base($"ChainSystem(N={system.N}, J={system.J}, gamma={system.GammaZero})",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/ChainSystems/ChainSystem.cs")
    {
        System = system;
    }

    public override string DisplayName =>
        $"ChainSystemPrimitive (N={System.N}, {System.HType}, {System.Topology})";

    public override string Summary =>
        $"N={System.N} qubits, J={System.J}, gamma={System.GammaZero}, H={System.HType}, topology={System.Topology}";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("N", System.N);
            yield return InspectableNode.RealScalar("J", System.J);
            yield return InspectableNode.RealScalar("gamma_zero", System.GammaZero);
            yield return new InspectableNode("HType", summary: System.HType.ToString());
            yield return new InspectableNode("Topology", summary: System.Topology.ToString());
        }
    }
}
