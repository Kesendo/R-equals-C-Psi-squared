using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F89OcticMonodromyClaim"/> (monodromy = S_8, the geometric
/// route). Two typed parent edges: <see cref="F89Path3OcticGaloisClaim"/> (the algebraic S_8 this braid
/// route reproduces) and <see cref="F89Path3OcticEpClaim"/> (the EP locations the lassos enclose; the
/// silent diabolic). Both parents Tier 1 derived. The builder topo-resolves, so registration order is
/// free, but both parents must be registered.</summary>
public static class F89OcticMonodromyClaimRegistration
{
    public static ClaimRegistryBuilder RegisterF89OcticMonodromyClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F89OcticMonodromyClaim>(b =>
            new F89OcticMonodromyClaim(
                b.Get<F89Path3OcticGaloisClaim>(),
                b.Get<F89Path3OcticEpClaim>()));
}
