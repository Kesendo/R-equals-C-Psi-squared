using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.F1Family;

/// <summary>Builder extension that registers the F1 family: <see cref="ChainSystemPrimitive"/>,
/// then <see cref="F1PalindromeIdentity"/>, then five Tier-1-derived F1 children that the
/// <see cref="F1KnowledgeBase"/> exposes as top-level properties (the KB's sixth Tier-1
/// child <c>SingleBodyScaling</c> is the omitted one; see the SingleBody paragraph below),
/// plus the Tier-2-verified <see cref="F1GeneralTopologyVerifiedClaim"/> general-topology
/// verification record, plus three F1-anchored topology-bound bridges
/// (<see cref="F4KernelDimensionByComponentsClaim"/>, <see cref="RingN4DihedralLockClaim"/>,
/// <see cref="StarImMaxBoundClaim"/>):
///
/// <list type="bullet">
///   <item><see cref="PalindromeResidualScalingClaim"/> (<see cref="HamiltonianClass.Main"/>):
///         residual ‖M‖² scaling on the main class. The <paramref name="hClass"/> parameter
///         on <c>RegisterF1Family</c> selects which class this single registration carries.</item>
///   <item><see cref="F1T1ResidualClosedForm"/>: ‖M(T1)‖² = 4^(N−1)·[3·Σγ² + 4·(Σγ)²].</item>
///   <item><see cref="F1T1ResidualPi2Decomposition"/>: Π²-orthogonal Pythagorean split of
///         the T1 residual into (anti, sym) parts; depends on F1T1ResidualClosedForm.</item>
///   <item><see cref="F1DepolResidualClosedForm"/>: ‖M(depol)‖² = 4^(N−1)·[(16/9)·Σγ² + 16·(Σγ)²].</item>
///   <item><see cref="F49NonUniformCrossTermClaim"/>: ‖{L_H, L_Dc}‖² = 4·Σ_b ‖L_H^bond‖²·Σ_{m∉bond}γ_m² +
///         Σ_b G(bond, H)·(γ_i−γ_j)²; F49's non-uniform γ extension.</item>
///   <item><see cref="F1GeneralTopologyVerifiedClaim"/> (Tier 2): verification record
///         that the (B, D2) parameterisation of ‖M(N, G)‖² extends bit-exactly to
///         disconnected, weighted, and random connected graphs at N=5..7. Depends on
///         <see cref="PalindromeResidualScalingClaim"/> (the closed form whose
///         universality is verified) and <see cref="F1PalindromeIdentity"/> (parent F1).</item>
///   <item><see cref="F4KernelDimensionByComponentsClaim"/> (Tier 1 derived,
///         promoted 2026-05-19; originally landed Tier 1 candidate 2026-05-18): F4
///         disconnected-graph extension dim ker L_H = Π_c (|c|+1) surfaced by the F1
///         SLOW_N8 sweep. The formula owner is the F4 family (parent F4 =
///         <see cref="F4StationaryModeCountPi2Inheritance"/>) but the claim sits
///         inside the F1 registration because the four bit-exact N=8 anchors live in
///         the same JSON metric files written by the F1 sweep. The connected-case
///         upper bound is closed by DEGENERACY_PALINDROME Result 2 (magnetization
///         conservation); the multi-component product follows from standard
///         tensor-sum kernel factorisation. Parent edge: depend on
///         <see cref="F1PalindromeIdentity"/> (Tier 1 derived, strength 5) which is
///         the topology-order anchor of the entire F1 family. The same edge that
///         served the Tier 1 candidate (strength 4) child satisfies the now
///         Tier 1 derived (strength 5) child (<see cref="TierStrength"/> check
///         5 ≥ 5). <see cref="F1GeneralTopologyVerifiedClaim"/> (Tier 2 verified,
///         strength 3) cannot serve as direct parent; the sister-relationship is
///         documented in the claim's XML doc and proof file instead.</item>
///   <item><see cref="RingN4DihedralLockClaim"/> (Tier 1 derived, 2026-05-19):
///         <c>Im_max(ring, N=4, J) = (3/4)·J·N = 3·J</c> Q-universal Im-spectral
///         saturation surfaced by the Q-sweep extension of the SLOW_N8 + N=9 chain
///         sprint. Closed form via the K_{2,2} = C_4 bipartite-complete graph
///         isomorphism: <c>H = J·S⃗_A·S⃗_B</c> with sublattice total spins, Casimir
///         spectrum <c>{−2J, −J, 0³, +J}</c>, max gap 3J realised by the
///         <c>|Ψ_+⟩⟨Ψ_−|</c> Liouvillian eigenmode between the S_tot=2 ferromagnet
///         and the (S_A=1, S_B=1, S_tot=0) singlet. Pure-dephasing dissipator only
///         adds real decay so no L-mode exceeds the H-spread bound. Parent edge:
///         <see cref="F1PalindromeIdentity"/> (Tier 1 derived, strength 5); the
///         Im-max bound lives in the L-spectrum the F1 palindrome partitions, and
///         the eigenmode-construction machinery is shared with the F4 kernel-dim
///         sister bridge above. See <c>docs/proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md</c>.</item>
///   <item><see cref="StarImMaxBoundClaim"/> (Tier 1 derived, 2026-05-19):
///         <c>Im_max(star, N, J) = J·N/2</c> Q-universal saturation surfaced by
///         the same Q-sweep (24 anchors at γ₀=0.05) plus the SLOW_N8 anchor at
///         <c>star_N8.json</c>. Closed form via the SU(2)/Schur-Weyl hub-leaf
///         Casimir factorisation <c>H_star = J·S⃗_0·S⃗_L</c> on the star bipartite
///         split A = {hub}, B = {N-1 leaves}; maximum-leaf-spin S_L = (N-1)/2
///         ferromagnetic sector gives the ΔE_max = J·N/2 realised by the
///         <c>|Ψ_+⟩⟨Ψ_−|</c> eigenmode between the S_tot = N/2 fully-aligned and
///         the S_tot = (N-2)/2 hub-anti-aligned states. Same proof skeleton as
///         <see cref="RingN4DihedralLockClaim"/> with sublattice sizes (1, N-1)
///         instead of (2, 2); the Marrakesh-convention <c>Im/σ = 1 ↔ J = 2γ</c>
///         reading is the Q=2 column of the universal <c>Im/σ = Q/2</c> lock.
///         Parent edge: <see cref="F1PalindromeIdentity"/> (Tier 1 derived,
///         strength 5). See <c>docs/proofs/PROOF_STAR_OPTICAL_CONFOCAL_SATURATION.md</c>.</item>
/// </list>
///
/// <para>Dependency edges (parent → child):</para>
/// <list type="bullet">
///   <item>ChainSystemPrimitive → F1PalindromeIdentity (topological-order hint;
///         F1PalindromeIdentity has no Core-side parameter but the chain parameterises
///         the surrounding KB context).</item>
///   <item>ChainSystemPrimitive → PalindromeResidualScalingClaim (the N from the chain
///         is what the scaling claim is built against).</item>
///   <item>F1PalindromeIdentity → PalindromeResidualScalingClaim.</item>
///   <item>F1PalindromeIdentity → F1T1ResidualClosedForm.</item>
///   <item>F1PalindromeIdentity → F1T1ResidualPi2Decomposition.</item>
///   <item>F1T1ResidualClosedForm → F1T1ResidualPi2Decomposition (the Π²-decomposition
///         closes the parent total Pythagorically: (anti) + (sym) = (3·Σγ² + 4·(Σγ)²)).</item>
///   <item>F1PalindromeIdentity → F1DepolResidualClosedForm.</item>
///   <item>F1PalindromeIdentity → F49NonUniformCrossTermClaim (extends F49 by relaxing
///         the uniform-γ assumption; the F1 σ-shift L_Dc = L_D + σ·I that frames the
///         cross-term is the F1 identity's centering convention).</item>
///   <item>F1PalindromeIdentity → F1GeneralTopologyVerifiedClaim (the master F1
///         identity grounds the verification record).</item>
///   <item>PalindromeResidualScalingClaim → F1GeneralTopologyVerifiedClaim (the
///         verification record asserts universality of the (B, D2) parameterisation
///         carried by the scaling claim).</item>
/// </list>
///
/// <para><b>SingleBody scaling: deliberately omitted.</b> The
/// <see cref="F1KnowledgeBase"/> also exposes a <c>SingleBodyScaling</c> Tier-1 property
/// (the same <see cref="PalindromeResidualScalingClaim"/> type with
/// <see cref="HamiltonianClass.SingleBody"/>). The runtime
/// <see cref="ClaimRegistryBuilder"/> is type-keyed (one factory per
/// <see cref="System.Type"/>); attempting to register a second
/// <see cref="PalindromeResidualScalingClaim"/> with a different class throws
/// <c>InvariantViolationException(rule: "DuplicateRegistration")</c>. No native
/// discriminator parameter or generic wrapper key exists on the builder, so we surface
/// only one scaling claim per registry build. The <paramref name="hClass"/> parameter
/// lets the caller choose which class this is (Main by default, matching the original
/// behaviour). The KB still exposes both classes for inspection.</para></summary>
public static class F1FamilyRegistration
{
    /// <summary>Register the eleven F1-family Claims (ChainSystemPrimitive +
    /// F1PalindromeIdentity + PalindromeResidualScalingClaim + F1T1ResidualClosedForm +
    /// F1T1ResidualPi2Decomposition + F1DepolResidualClosedForm +
    /// F49NonUniformCrossTermClaim + F1GeneralTopologyVerifiedClaim +
    /// F4KernelDimensionByComponentsClaim + RingN4DihedralLockClaim +
    /// StarImMaxBoundClaim) for a given <paramref name="chain"/>. Default
    /// Hamiltonian class for the scaling claim is <see cref="HamiltonianClass.Main"/>;
    /// chain bond count and degree-squared sum default to <c>null</c> (use
    /// <c>FactorChain</c>).</summary>
    public static ClaimRegistryBuilder RegisterF1Family(
        this ClaimRegistryBuilder builder,
        ChainSystem chain,
        HamiltonianClass hClass = HamiltonianClass.Main)
    {
        return builder
            .Register<ChainSystemPrimitive>(_ => new ChainSystemPrimitive(chain))
            .Register<F1PalindromeIdentity>(b =>
            {
                _ = b.Get<ChainSystemPrimitive>();
                return new F1PalindromeIdentity();
            })
            .Register<PalindromeResidualScalingClaim>(b =>
            {
                var primitive = b.Get<ChainSystemPrimitive>();
                _ = b.Get<F1PalindromeIdentity>();
                return new PalindromeResidualScalingClaim(
                    N: primitive.System.N,
                    hClass: hClass);
            })
            .Register<F1T1ResidualClosedForm>(b =>
            {
                _ = b.Get<F1PalindromeIdentity>();
                return new F1T1ResidualClosedForm();
            })
            .Register<F1T1ResidualPi2Decomposition>(b =>
            {
                _ = b.Get<F1PalindromeIdentity>();
                _ = b.Get<F1T1ResidualClosedForm>();
                return new F1T1ResidualPi2Decomposition();
            })
            .Register<F1DepolResidualClosedForm>(b =>
            {
                _ = b.Get<F1PalindromeIdentity>();
                return new F1DepolResidualClosedForm();
            })
            .Register<F49NonUniformCrossTermClaim>(b =>
            {
                _ = b.Get<F1PalindromeIdentity>();
                return new F49NonUniformCrossTermClaim();
            })
            .Register<F1GeneralTopologyVerifiedClaim>(b =>
            {
                // The general-topology verification record is anchored to the (B, D2)
                // closed form of PalindromeResidualScalingClaim and the master F1
                // identity. Both parents are pulled into the dependency graph so the
                // Tier-2 verified record sits downstream of its Tier-1 anchors.
                _ = b.Get<F1PalindromeIdentity>();
                _ = b.Get<PalindromeResidualScalingClaim>();
                return new F1GeneralTopologyVerifiedClaim();
            })
            .Register<F4KernelDimensionByComponentsClaim>(b =>
            {
                // F4 disconnected-graph extension surfaced by the F1 SLOW_N8 sweep
                // (2026-05-18). Promoted Tier 1 candidate → Tier 1 derived 2026-05-19
                // after DEGENERACY_PALINDROME Result 2 was identified as the closure
                // of the connected-case upper bound. The claim's formula owner is F4
                // (parent), but the anchors live in the F1 SLOW_N8 JSON metric files.
                // We pull F1PalindromeIdentity (Tier 1 derived, strength 5) as the
                // explicit parent: a Tier 1 derived (strength 5) child needs a parent
                // of strength ≥ 5, and 5 ≥ 5 holds.
                // F1GeneralTopologyVerifiedClaim (Tier 2 verified, strength 3) cannot
                // serve as parent directly; the sister-relationship is documented in
                // the claim's XML doc and proof file instead.
                _ = b.Get<F1PalindromeIdentity>();
                return new F4KernelDimensionByComponentsClaim();
            })
            .Register<RingN4DihedralLockClaim>(b =>
            {
                // Ring N=4 dihedral lock surfaced by the 2026-05-19 Q-sweep extension
                // of the F1 SLOW_N8 + N=9 chain bridge sprint: Im_max(ring, N=4, J) =
                // (3/4)·J·N = 3·J Q-universal, derived from the K_{2,2} = C_4 bipartite-
                // complete graph isomorphism and the resulting Casimir spectrum {−2J,
                // −J, 0³, +J} with max gap 3J. The Liouvillian eigenmode realising the
                // bound is the standard transition between S_tot=2 ferromagnet and
                // (S_A=1, S_B=1, S_tot=0) singlet. Pure-dephasing dissipator only adds
                // real decay so no L-mode can exceed the H-spread bound. Tier 1
                // derived; see PROOF_RING_N4_DIHEDRAL_LOCK.md.
                //
                // Parent edge: F1PalindromeIdentity (Tier 1 derived, strength 5)
                // serves; the Im-max bound lives in the same L-spectrum the F1
                // palindrome partitions, and the analytical machinery (eigenmode
                // construction, dissipator decoupling) is the same as for the F4
                // kernel-dim sister claim above.
                _ = b.Get<F1PalindromeIdentity>();
                return new RingN4DihedralLockClaim();
            })
            .Register<StarImMaxBoundClaim>(b =>
            {
                // Star Im-max saturation surfaced by the 2026-05-19 Q-sweep extension
                // of the F1 SLOW_N8 + N=9 chain bridge sprint: Im_max(star, N, J) =
                // J·N/2 Q-universal for any N ≥ 3, derived from the SU(2)/Schur-Weyl
                // hub-leaf Casimir factorisation H_star = J·S⃗_0·S⃗_L. The maximum-leaf-
                // spin S_L = (N-1)/2 ferromagnetic sector gives ΔE_max = J·N/2,
                // realised by the Liouvillian eigenmode between the S_tot = N/2
                // fully-aligned state and the S_tot = (N-2)/2 hub-anti-aligned state.
                // Pure-dephasing dissipator only adds real decay so no L-mode can
                // exceed the H-spread bound. Tier 1 derived; see
                // PROOF_STAR_OPTICAL_CONFOCAL_SATURATION.md.
                //
                // Parent edge: F1PalindromeIdentity (Tier 1 derived, strength 5)
                // serves, same as for the Ring N=4 sister bound above. The two
                // claims share the same proof skeleton (bipartite split → all-pairs
                // bonding → Casimir → maximum-S_tot ferromagnet eigenmode).
                _ = b.Get<F1PalindromeIdentity>();
                return new StarImMaxBoundClaim();
            });
    }
}
