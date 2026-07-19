using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wiring of <see cref="RecordParityLawClaim"/> (F135; Tier1Derived). One typed
/// parent: <see cref="AbsorptionTheoremClaim"/> — Proposition 1 (the pair reduction) is the
/// absorption substrate |A⟩⟨B|(t) read on a pair page, and the record radius β_j is the
/// surviving coherence magnitude. Resolution is topological, so registration order among
/// siblings is immaterial.</summary>
public static class RecordParityLawClaimRegistration
{
    public static ClaimRegistryBuilder RegisterRecordParityLawClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<RecordParityLawClaim>(b =>
            new RecordParityLawClaim(b.Get<AbsorptionTheoremClaim>()));
}
