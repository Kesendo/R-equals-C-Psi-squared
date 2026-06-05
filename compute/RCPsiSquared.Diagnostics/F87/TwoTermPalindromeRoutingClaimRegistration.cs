using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>Schicht-1 wiring of <see cref="TwoTermPalindromeRoutingClaim"/> (the Liouvillian-free
/// two-term Q-pair router). Parameterised by <see cref="ChainSystemPrimitive"/> (the registry chain
/// supplies the Claim's identity; the authority cross-check itself runs on an internal N=4 chain) and
/// depends on three typed parent claims:
/// <see cref="F87TrichotomyClassification"/> (the spectral authority the router is verified against),
/// <see cref="F87DiagonalCellBipartiteWitnessSet"/> (the diagonal P1-family special case this
/// generalises to the full hidden-Q routing), and <see cref="CrossoverMirrorSqrtNinetyClaim"/> (the
/// Continuous-family crossover mirror's existing C# realisation).
///
/// <para>Requires: <see cref="Runtime.F1Family.F1FamilyRegistration.RegisterF1Family"/> for
/// <see cref="ChainSystemPrimitive"/>,
/// <see cref="F87FamilyRegistration.RegisterF87Family"/> for
/// <see cref="F87TrichotomyClassification"/>,
/// <see cref="F87DiagonalCellBipartiteWitnessSetRegistration.RegisterF87DiagonalCellBipartiteWitnessSet"/>
/// for <see cref="F87DiagonalCellBipartiteWitnessSet"/>, and the CrossoverMirror registration for
/// <see cref="CrossoverMirrorSqrtNinetyClaim"/>; the builder errors with <c>MissingParent</c> if any
/// are absent.</para>
///
/// <para>Tier consistency: TwoTermPalindromeRoutingClaim is Tier2Empirical ←
/// F87TrichotomyClassification (Tier1Derived), ← F87DiagonalCellBipartiteWitnessSet (Tier1Candidate),
/// and ← CrossoverMirrorSqrtNinetyClaim (Tier1Derived). All three parents are at least as strong as
/// the child (parent ≥ child), so the strength-inheritance check passes.</para></summary>
public static class TwoTermPalindromeRoutingClaimRegistration
{
    public static ClaimRegistryBuilder RegisterTwoTermPalindromeRoutingClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<TwoTermPalindromeRoutingClaim>(b =>
        {
            _ = b.Get<F87TrichotomyClassification>();
            _ = b.Get<F87DiagonalCellBipartiteWitnessSet>();
            _ = b.Get<CrossoverMirrorSqrtNinetyClaim>();
            var chain = b.Get<ChainSystemPrimitive>();
            return new TwoTermPalindromeRoutingClaim(chain.System);
        });
}
