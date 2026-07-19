using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wiring of <see cref="RecordLetterLawClaim"/> (F136; Tier1Derived). One typed
/// parent: <see cref="RecordParityLawClaim"/> (F135) — every channel of the letter law is
/// Proposition 1 read on a channel class, and Law A is the D = ∅ column of the trichotomy.
/// The MirrorWorld adoption (Witness.cs) is a CONSUMER of F136, not a parent, and MirrorWorld
/// is standalone: no edge points there.</summary>
public static class RecordLetterLawClaimRegistration
{
    public static ClaimRegistryBuilder RegisterRecordLetterLawClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<RecordLetterLawClaim>(b =>
            new RecordLetterLawClaim(b.Get<RecordParityLawClaim>()));
}
