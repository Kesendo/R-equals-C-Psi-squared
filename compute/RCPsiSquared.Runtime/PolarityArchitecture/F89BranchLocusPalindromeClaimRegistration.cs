using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F89BranchLocusPalindromeClaim"/>: the path-3 octic branch
/// locus is mirror-symmetric about Re λ = −4, forced by the F1 palindrome carried antiunitarily on the
/// (SE,DE) block. Two typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="F1PalindromeIdentity"/>: the F1 palindrome the branch locus inherits (the mirror
///         λ ↦ −λ − 2σ realised on the block as the antiunitary T).</item>
///   <item><see cref="F89Path3OcticEpClaim"/>: the octic and the diabolic on the line, with the AT-midpoint
///         centre −σ = −4 (the absorption rungs −2γ, −6γ).</item>
/// </list>
///
/// <para>Tier consistency: Tier 1 derived (forced by the antiunitary palindrome + verified to machine
/// precision). Both parents Tier 1 derived.</para>
///
/// <para>Requires <see cref="F1PalindromeIdentity"/> and <see cref="F89Path3OcticEpClaim"/> to be
/// registered (the builder topo-resolves, so the order of registration is free).</para></summary>
public static class F89BranchLocusPalindromeClaimRegistration
{
    public static ClaimRegistryBuilder RegisterF89BranchLocusPalindromeClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F89BranchLocusPalindromeClaim>(b =>
            new F89BranchLocusPalindromeClaim(
                b.Get<F1PalindromeIdentity>(),
                b.Get<F89Path3OcticEpClaim>()));
}
