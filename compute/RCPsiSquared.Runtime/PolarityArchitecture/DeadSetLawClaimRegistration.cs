using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="DeadSetLawClaim"/> (F132, 2026-07-16): the
/// dead-set law. The mirror-composition antiunitaries V_g = (U_g·X^N)∘conj are exact
/// symmetries of the XY+field Lindblad flow at every fixed h; their kill sign is the
/// N-free mod-4 function of the Majorana degree (ε_odd = (−1)^(d(d−1)/2),
/// ε_even = (−1)^(d(d+1)/2)); with the popcount blocks and the conserved degree the
/// identically-zero-at-every-h readout set closes into one line:
/// alive ⟺ (K_pop ∧ d ≡ 0 mod 4) ∨ (coherence ∧ K_coh ∧ d = N). Necessity derived,
/// sufficiency gated; the zz knob is the free-world fence.
///
/// <para>Tier1Derived (the necessity face; self-check battery at N = 3 in the ctor,
/// matrix + moment level, no eigensolver). Typed parents: <see cref="ChiralKClaim"/>
/// (the sublattice gauge, the hopping flip), <see cref="AntilinearTriangleClaim"/>
/// (conj as the antilinear leg), <see cref="MirrorOrderSortingClaim"/> (F131, whose
/// third-mirror sighting opened the doubly-mirrored zeros this law closes). Anchor:
/// <c>experiments/LATTICE_DEAD_SET_RULE.md</c> +
/// <c>simulations/lattice_dead_set_rule.py</c> (gate, 21 checks) +
/// <c>docs/ANALYTICAL_FORMULAS.md</c> (F132).</para></summary>
public static class DeadSetLawClaimRegistration
{
    public static ClaimRegistryBuilder RegisterDeadSetLawClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<DeadSetLawClaim>(b =>
            new DeadSetLawClaim(
                b.Get<ChiralKClaim>(),
                b.Get<AntilinearTriangleClaim>(),
                b.Get<MirrorOrderSortingClaim>()));
}
