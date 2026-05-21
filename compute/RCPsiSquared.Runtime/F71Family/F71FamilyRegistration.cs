using RCPsiSquared.Core.F71;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.F71Family;

/// <summary>Registers the F71 mirror-symmetry family: <see cref="C1MirrorIdentity"/> as
/// the analytic identity, two N-parameterised consequences (<see cref="F71MirrorOperator"/>,
/// <see cref="F71BondOrbitDecomposition"/>), the cross-F86 generalisation
/// (<see cref="F86MirrorGeneralisationLink"/>) that lifts c₁(b) = c₁(N−2−b) to
/// Q_peak(b) = Q_peak(N−2−b), the F100 non-uniform-J extension
/// (<see cref="C1QPeakMirrorJParity"/>), and the F101 non-uniform-γ extension
/// (<see cref="C1MirrorGammaParity"/>).
///
/// <code>
///   C1MirrorIdentity (analytic identity)
///     ├── F71MirrorOperator (N-parameterised)
///     ├── F71BondOrbitDecomposition (N-parameterised)
///     ├── F86MirrorGeneralisationLink (cross-KB lift to F86)
///     ├── C1QPeakMirrorJParity (F100: graceful breakdown under non-uniform J)
///     └── C1MirrorGammaParity (F101: graceful breakdown under non-uniform γ)
/// </code>
///
/// <para>All six claims are Tier1Derived. F71MirrorOperator and F71BondOrbitDecomposition
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
            })
            .Register<C1QPeakMirrorJParity>(b =>
            {
                _ = b.Get<C1MirrorIdentity>();
                return new C1QPeakMirrorJParity();
            })
            .Register<C1MirrorGammaParity>(b =>
            {
                _ = b.Get<C1MirrorIdentity>();
                return new C1MirrorGammaParity();
            });
}
