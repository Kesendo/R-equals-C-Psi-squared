using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F89Path3OcticGaloisClaim"/>: Gal(F_8)
/// ⊄ A_8 from disc-non-square (Tier 1); non-solvability conjecture Tier 2 open.
/// Two typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="F89TopologyOrbitClosure"/>: orbit-closure framework parent.</item>
///   <item><see cref="F89PathKAtLockMechanismClaim"/>: the H_B-mixed sub-factor
///         (where the octic lives) is defined relative to the AT-lock subspace.</item>
/// </list>
///
/// <para>Tier consistency: Tier 1 derived for the Gal ⊄ A_8 statement. The
/// stronger non-solvability conjecture is documented as Tier 2 open in the
/// docstring and is NOT promoted to a separate Claim.</para>
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
