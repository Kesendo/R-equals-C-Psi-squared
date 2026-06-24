using System.Collections.Generic;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The F89 path-k H_B-mixed Galois groups, k = 3..6: Gal(F_d/Q(i)(q)) = S_d (non-solvable)
/// ⟹ the Liouvillian decay rates λ_k(q) admit no radical closure in q = J/γ (Abel-Ruffini). The
/// witness-first capstone for the landed path-3..6 = S_8/18/32/53 result: path-3 (octic) is
/// recomputed fully live (block → Berkowitz → isolate F_d by dividing out the AT factor); path-4/5/6
/// are semi-live (the committed oracle F_d is verified against the live block charpoly, then read).
/// Breadcrumbed from F89PathKHbMixedDegreesClaim (degree table Tier1Candidate / Galois verdict
/// Tier1Derived).</summary>
public sealed class F89PathKGaloisWitness : IInspectable
{
    public string DisplayName => "F89 path-k H_B-mixed Galois groups (k=3..6): S_8/18/32/53, non-solvable";

    public string Summary =>
        "Gal(F_d/Q(i)(q)) = S_d for d = 8,18,32,53 (the path-3..6 H_B-mixed Liouvillian factors) — non-solvable, " +
        "so the decay rates λ_k(q) admit no radical closure in q = J/γ. path-3 fully live; path-4/5/6 semi-live (the committed F_d verified against the live block charpoly).";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new F89OcticGaloisWitness();            // path-3, fully live (octic isolated from the block)
            yield return new F89PathSemiLiveGaloisWitness(4);    // path-4, semi-live (oracle verified live)
            yield return new F89PathSemiLiveGaloisWitness(5);    // path-5, semi-live
            yield return new F89PathSemiLiveGaloisWitness(6);    // path-6, semi-live (heavy: degree-53 / 75×75 Berkowitz)
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
