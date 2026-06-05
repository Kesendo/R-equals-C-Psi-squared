using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>Schicht-1 wiring of <see cref="PalindromeSoftCertifierClaim"/> (the §7.12 Liouvillian-free
/// soft-certifier). Parameterised by <see cref="ChainSystemPrimitive"/> (the certifier's soundness and
/// the ceiling witness are spectrally checked at the registered chain's N) and depends on two typed
/// parent claims:
/// <see cref="F87DiagonalCellBipartiteWitnessSet"/> (the §7 diagonal-K bipartite criterion the
/// certifier's linear strategy scales) and <see cref="F87TrichotomyClassification"/> (the spectral
/// authority the soundness is checked against).
///
/// <para>Requires: <see cref="Runtime.F1Family.F1FamilyRegistration.RegisterF1Family"/> for
/// <see cref="ChainSystemPrimitive"/>,
/// <see cref="F87DiagonalCellBipartiteWitnessSetRegistration.RegisterF87DiagonalCellBipartiteWitnessSet"/>
/// for <see cref="F87DiagonalCellBipartiteWitnessSet"/>, and
/// <see cref="F87FamilyRegistration.RegisterF87Family"/> for
/// <see cref="F87TrichotomyClassification"/>; the builder errors with <c>MissingParent</c> if any are
/// absent.</para>
///
/// <para>Tier consistency: PalindromeSoftCertifierClaim is Tier1Candidate ←
/// F87DiagonalCellBipartiteWitnessSet (Tier1Candidate) and ← F87TrichotomyClassification
/// (Tier1Derived). The strength-inheritance check is parent ≥ child, i.e. 4 ≥ 4 and 5 ≥ 4, both
/// pass.</para></summary>
public static class PalindromeSoftCertifierClaimRegistration
{
    public static ClaimRegistryBuilder RegisterPalindromeSoftCertifierClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<PalindromeSoftCertifierClaim>(b =>
        {
            _ = b.Get<F87DiagonalCellBipartiteWitnessSet>();
            _ = b.Get<F87TrichotomyClassification>();
            var chain = b.Get<ChainSystemPrimitive>();
            return new PalindromeSoftCertifierClaim(chain.System);
        });
}
