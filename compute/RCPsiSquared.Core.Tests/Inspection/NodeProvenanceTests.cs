using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Core.Tests.Inspection;

/// <summary>The [live]/[stored] provenance badge (stranger_door's fifth door). The rule is
/// value-origin honesty: a frozen carrier is Stored, a computational witness is Live, and a
/// node opts into the other when its displayed value's true origin differs from its type.</summary>
public class NodeProvenanceTests
{
    [Fact]
    public void FrozenCarrier_DefaultsToStored()
    {
        // A bare InspectableNode holds a string frozen at construction -> written down, not recomputed.
        var node = new InspectableNode("n", "s");
        Assert.Equal(NodeProvenance.Stored, node.Provenance);
    }

    [Fact]
    public void Carrier_CanOptIntoLive()
    {
        // A live producer minting a child for a value it just computed stamps it Live.
        var node = new InspectableNode("n", "s", provenance: NodeProvenance.Live);
        Assert.Equal(NodeProvenance.Live, node.Provenance);
    }

    [Fact]
    public void Witness_WithoutOverride_IsLive()
    {
        // The interface contract is "a computational unit, not a data carrier": default Live.
        IInspectable witness = new BareWitness();
        Assert.Equal(NodeProvenance.Live, witness.Provenance);
    }

    /// <summary>A minimal real IInspectable that declares only the four required members and
    /// no Provenance, so it exercises the interface default.</summary>
    private sealed class BareWitness : IInspectable
    {
        public string DisplayName => "bare";
        public string Summary => "";
        public IEnumerable<IInspectable> Children => Array.Empty<IInspectable>();
        public InspectablePayload Payload => InspectablePayload.Empty;
    }
}
