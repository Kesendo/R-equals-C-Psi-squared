using RCPsiSquared.Core.F71;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.F86Main;

/// <summary>Schicht-1 wiring for <see cref="PerF71OrbitObservation"/> — the Tier2Empirical
/// refinement of F86 Statement 3 (F71 mirror invariance) that pins the per-F71-orbit
/// Q_peak substructure within the Interior bond class.
///
/// <para>This is the typed-knowledge surface for the empirical observation that
/// motivated the JW track's truly-vs-flanking-innermost question: the c=2 N=6 witness
/// records central Q_peak ≈ 1.43 vs flanking Q_peak ≈ 1.63 — i.e., Interior is NOT a
/// single F71 orbit at N ≥ 6. Wiring it as a Schicht-1 claim makes the inversion +
/// high-Q-plateau facts queryable as Edge data instead of leaving them as floating
/// memory notes.</para>
///
/// <para>Edge declared: <see cref="PerF71OrbitObservation"/> ←
/// <see cref="F71BondOrbitDecomposition"/>. The orbit table refines the F71 bond-orbit
/// partitioning with (c, N)-specific Q_peak measurements; the parent provides the orbit
/// structure itself.</para>
///
/// <para>Tier consistency: PerF71OrbitObservation is Tier2Empirical (frozen pinned
/// witnesses, no closed form);F71BondOrbitDecomposition is Tier1Derived. TierStrength
/// inheritance passes (5 ≥ 2). The closed-form classification of orbit pattern + plateau
/// enhancement is the open Tier1 promotion path documented in PROOF_F86_QPEAK Open
/// Element 5.</para>
///
/// <para>Requires: <see cref="F71Family.F71FamilyRegistration.RegisterF71Family"/> must
/// have been called first so <see cref="F71BondOrbitDecomposition"/> is available as a
/// parent. The N parameter passed to that registration is independent of the witnesses
/// here — the witnesses pin specific (c, N) cases as a frozen table.</para></summary>
public static class F86PerF71OrbitObservationRegistration
{
    public static ClaimRegistryBuilder RegisterF86PerF71OrbitObservation(
        this ClaimRegistryBuilder builder) =>
        builder.Register<PerF71OrbitObservation>(b =>
        {
            _ = b.Get<F71BondOrbitDecomposition>();
            return new PerF71OrbitObservation();
        });
}
