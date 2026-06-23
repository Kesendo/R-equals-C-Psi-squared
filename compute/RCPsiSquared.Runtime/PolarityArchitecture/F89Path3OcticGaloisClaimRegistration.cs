using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F89Path3OcticGaloisClaim"/>:
/// Gal(F_8 / Q(i)(q)) = S_8 (Tier 1 derived; non-solvable, no radical closure in q),
/// via disc-non-square (⊄ A_8) plus a (5,2,1) Frobenius certificate at 𝔭|5, q0=2
/// (5-cycle ⇒ primitive ⇒ ⊇A_8 by Jordan). Two typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="F89TopologyOrbitClosure"/>: orbit-closure framework parent.</item>
///   <item><see cref="F89PathKAtLockMechanismClaim"/>: the H_B-mixed sub-factor
///         (where the octic lives) is defined relative to the AT-lock subspace.</item>
/// </list>
///
/// <para>Tier consistency: Tier 1 derived for Gal = S_8 (and the non-solvability /
/// no-radical-closure consequence). The certificate (specialization + Dedekind +
/// Jordan) lives in the docstring + the cited script; not promoted to a separate Claim.</para>
///
/// <para>Requires: <see cref="F89TopologyOrbitClosureRegistration.RegisterF89TopologyOrbitClosure"/>
/// and <see cref="F89PathKAtLockMechanismClaimRegistration.RegisterF89PathKAtLockMechanismClaim"/>.</para></summary>
public static class F89Path3OcticGaloisClaimRegistration
{
    public static ClaimRegistryBuilder RegisterF89Path3OcticGaloisClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F89Path3OcticGaloisClaim>(b =>
        {
            var f89 = b.Get<F89TopologyOrbitClosure>();
            var atLock = b.Get<F89PathKAtLockMechanismClaim>();
            return new F89Path3OcticGaloisClaim(f89, atLock);
        });
}
