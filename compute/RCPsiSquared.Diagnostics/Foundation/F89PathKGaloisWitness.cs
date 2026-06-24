using System.Collections.Generic;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The F89 path-k H_B-mixed Galois groups, k = 3..6: Gal(F_d/Q(i)(q)) = S_d (non-solvable)
/// ⟹ the Liouvillian decay rates λ_k(q) admit no radical closure in q = J/γ (Abel-Ruffini). The
/// witness-first capstone for the landed path-3..6 = S_8/18/32/53 result: ALL paths are recomputed
/// fully live (full Option D) — the block is built live, the AT-locked factor is reconstructed (path-3
/// from its F_a/F_b quadratics; path-4/5/6 from the rate-confined invariant subspace), divided out to
/// isolate F_d, and read by Frobenius. F_d is never imported. Breadcrumbed from
/// F89PathKHbMixedDegreesClaim (degree table Tier1Candidate / Galois verdict Tier1Derived).</summary>
public sealed class F89PathKGaloisWitness : IInspectable
{
    public string DisplayName => "F89 path-k H_B-mixed Galois groups (k=3..6): S_8/18/32/53, non-solvable";

    public string Summary =>
        "Gal(F_d/Q(i)(q)) = S_d for d = 8,18,32,53 (the path-3..6 H_B-mixed Liouvillian factors) — non-solvable, " +
        "so the decay rates λ_k(q) admit no radical closure in q = J/γ. All four paths fully live: F_d isolated from the live block, never imported.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new F89OcticGaloisWitness();            // path-3, fully live (octic from F_a/F_b)
            yield return new F89PathKLiveGaloisWitness(4);       // path-4, fully live (reconstructed AT)
            yield return new F89PathKLiveGaloisWitness(5);       // path-5, fully live
            yield return new F89PathKLiveGaloisWitness(6);       // path-6, fully live (degree-53 / 75×75 Berkowitz)
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
