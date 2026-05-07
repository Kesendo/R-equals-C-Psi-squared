using RCPsiSquared.Core.F71;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.F71Family;

/// <summary>Registers the F71 mirror-symmetry family: <see cref="C1MirrorIdentity"/> as
/// the analytic identity, two N-parameterised consequences (<see cref="F71MirrorOperator"/>,
/// <see cref="F71BondOrbitDecomposition"/>), and the cross-F86 generalisation
/// (<see cref="F86MirrorGeneralisationLink"/>) that lifts c₁(b) = c₁(N−2−b) to
/// Q_peak(b) = Q_peak(N−2−b).
///
/// <code>
///   C1MirrorIdentity (analytic identity)
///     ├── F71MirrorOperator (N-parameterised)
///     ├── F71BondOrbitDecomposition (N-parameterised)
///     └── F86MirrorGeneralisationLink (cross-KB lift to F86)
/// </code>
///
/// <para>All four claims are Tier1Derived. F71MirrorOperator and F71BondOrbitDecomposition
/// take an N argument; the registration takes the same N so the Runtime instance matches
/// the F1Family <see cref="ChainSystemPrimitive"/> N value used elsewhere.</para></summary>
public static class F71FamilyRegistration
{
    public static ClaimRegistryBuilder RegisterF71Family(
        this ClaimRegistryBuilder builder,
        int N) =>
        builder
            .Register<C1MirrorIdentity>(_ => new C1MirrorIdentity())
            .Register<F71MirrorOperator>(b =>
            {
                _ = b.Get<C1MirrorIdentity>();
                return new F71MirrorOperator(N);
            })
            .Register<F71BondOrbitDecomposition>(b =>
            {
                _ = b.Get<C1MirrorIdentity>();
                return new F71BondOrbitDecomposition(N);
            })
            .Register<F86MirrorGeneralisationLink>(b =>
            {
                _ = b.Get<C1MirrorIdentity>();
                return new F86MirrorGeneralisationLink();
            });
}
