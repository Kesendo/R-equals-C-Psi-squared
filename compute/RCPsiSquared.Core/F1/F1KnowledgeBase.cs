using RCPsiSquared.Core.Confirmations;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.F1;

/// <summary>"Everything we know about F1" as a typed OOP knowledge graph. The root
/// <see cref="IInspectable"/> for <c>--root f1</c>.
///
/// <para>F1 is the Liouvillian palindrome identity Π·L·Π⁻¹ = −L − 2Σγ·I, the namesake
/// theorem of the project (R=CΨ²). The KB is parametrised by chain length N, with optional
/// graph parameters (B, D2) for non-chain residual scaling.</para>
///
/// <para>What is in here:</para>
/// <list type="bullet">
///   <item>Tier-1 derived: <see cref="PalindromeIdentity"/> (F1 itself), main and
///         single-body residual scaling claims (<see cref="MainScaling"/>,
///         <see cref="SingleBodyScaling"/>), the T1 amplitude-damping closed form
///         (<see cref="T1ResidualClosedForm"/>) and its Π²-orthogonal Pythagorean
///         split into M_anti (F82/F84 amplitude-damping content) and M_sym
///         (<see cref="T1ResidualPi2Decomposition"/>), the depolarizing-noise
///         closed form (<see cref="DepolResidualClosedForm"/>), and the F49
///         non-uniform γ cross-term closed form
///         (<see cref="F49NonUniformCrossTerm"/>).</item>
///   <item>Tier-2 verified: hardware confirmations from
///         <see cref="ConfirmationsRegistry"/> that exercise F1 / palindrome trichotomy
///         on Marrakesh and related machines; plus the general-topology verification
///         record (<see cref="GeneralTopologyVerification"/>) that extends the
///         (B, D2) parameterisation to disconnected, weighted, and random connected
///         graphs at N=5, 6, 7, 8, 9 (N=8 via the opt-in SLOW_N8 block-spectrum dogfood,
///         N=9 chain via the SLOW_N9 dogfood routed through the <c>MklDirect</c> ILP64
///         bridge that landed 2026-05-19, both with full <see cref="F1SpectrumStatistics"/>
///         metric capture under <c>simulations/results/f1_n8_n9_metrics/</c>). New
///         frontier at N=10 is memory-pressure rather than the LP64 ceiling; see
///         <see cref="F1GeneralTopologyVerifiedClaim.ScaleFrontierBlockedAtN"/>.</item>
///   <item>Open: <see cref="OpenQuestions"/> is EMPTY as of 2026-05-18: the first
///         time the F1 family is open-question-free. All four items from the May 2026
///         sprint closed: T1 closed form (Tier-1 derived), depol closed form (Tier-1
///         derived), non-uniform γ (negative-result closure), and general topology
///         (synthesis proof + verification record). See <see cref="F1OpenQuestions"/>
///         XML doc for the per-item closure references.</item>
///   <item>Tier-1 derived bridge from the F1 N=8 SLOW_N8 sweep
///         (<see cref="KernelDimensionByComponents"/>): F4 disconnected-graph extension
///         dim ker L_H = Π_c (|c|+1) over connected components c. Promoted from
///         Tier 1 candidate to Tier 1 derived 2026-05-19 after DEGENERACY_PALINDROME
///         Result 2 (magnetization conservation) was identified as the closure of
///         the per-component upper bound. Bit-exact across the four N=8 topologies
///         captured in <see cref="F1GeneralTopologyVerifiedClaim.SpectrumMetricsDataFiles"/>;
///         the bonus discovery surfaced by that sweep. Lives on the F1 KB as a
///         bridge because the data came from F1 work, even though the formula itself
///         is structurally a child of F4 (parent
///         <see cref="F4StationaryModeCountPi2Inheritance"/>).</item>
///   <item>Tier-1 derived topology Im-max bound: Ring N=4
///         (<see cref="RingN4DihedralLock"/>) saturates Im_max(ring, N=4, J) =
///         (3/4)·J·N = 3·J Q-universally via the C_4 = K_{2,2} bipartite-complete
///         Casimir factorisation; 6 Q-sweep anchors bit-exact. See
///         <c>docs/proofs/PROOF_RING_N4_DIHEDRAL_LOCK.md</c>.</item>
///   <item>Tier-1 derived topology Im-max bound: Star
///         (<see cref="StarImMaxBound"/>) saturates Im_max(star, N, J) = J·N/2
///         Q-universally for any N ≥ 3 via the SU(2)/Schur-Weyl hub-leaf Casimir
///         factorisation H_star = J·S⃗_0·S⃗_L; 24 Q-sweep + 1 N=8 + 4 Python anchors
///         all bit-exact (29 anchors total). The Marrakesh-convention "Im/σ = 1
///         ↔ J = 2γ" reading is the Q = 2 specialization of the universal
///         Im/σ = Q/2 lock. See <c>docs/proofs/PROOF_STAR_OPTICAL_CONFOCAL_SATURATION.md</c>.</item>
/// </list>
/// </summary>
public sealed class F1KnowledgeBase : IInspectable
{
    private static readonly string[] _f1ConfirmationNames = new[]
    {
        "palindrome_trichotomy",
        "lebensader_skeleton_trace_decoupling",
        "f83_pi2_class_signature_marrakesh",
    };

    public int N { get; }
    public int? BondCount { get; }
    public int? DegreeSquaredSum { get; }
    public F1PalindromeIdentity PalindromeIdentity { get; }
    public PalindromeResidualScalingClaim MainScaling { get; }
    public PalindromeResidualScalingClaim SingleBodyScaling { get; }
    public F1T1ResidualClosedForm T1ResidualClosedForm { get; }
    public F1T1ResidualPi2Decomposition T1ResidualPi2Decomposition { get; }
    public F1DepolResidualClosedForm DepolResidualClosedForm { get; }
    public F49NonUniformCrossTermClaim F49NonUniformCrossTerm { get; }
    public F1GeneralTopologyVerifiedClaim GeneralTopologyVerification { get; }

    /// <summary>F4 disconnected-graph extension surfaced by the F1 SLOW_N8 sweep
    /// (2026-05-18): dim ker L_H = Π_c (|c|+1) over connected components c. The
    /// formula itself is structurally a child of F4 (parent
    /// <see cref="F4StationaryModeCountPi2Inheritance"/>), but it surfaces on the
    /// F1 KB as a bridge because the four bit-exact N=8 anchors live in the same
    /// JSON files written by the F1 sweep (see
    /// <see cref="F1GeneralTopologyVerifiedClaim.SpectrumMetricsDataFiles"/>).
    /// Tier 1 derived (promoted 2026-05-19): connected-case upper bound closed by
    /// <c>experiments/DEGENERACY_PALINDROME.md</c> Result 2 (magnetization
    /// conservation), multi-component product follows from standard tensor-sum
    /// kernel factorisation. See <c>PROOF_F4_KERNEL_DIMENSION_BY_COMPONENTS</c>
    /// § "Upper-bound closure (resolved 2026-05-18)".</summary>
    public F4KernelDimensionByComponentsClaim KernelDimensionByComponents { get; }

    /// <summary>Ring N=4 dihedral lock surfaced by the 2026-05-19 Q-sweep:
    /// Im_max(ring, N=4, J) = (3/4)·J·N = 3·J Q-universal. The 4-cycle is
    /// graph-isomorphic to the bipartite-complete graph K_{2,2}, so the Heisenberg
    /// Hamiltonian factors through total sublattice spins (H = J·S⃗_A·S⃗_B). The
    /// Casimir spectrum {−2J, −J, 0³, +J} pins the maximum H eigenvalue gap to 3J,
    /// which the Liouvillian eigenmode |Ψ_+⟩⟨Ψ_−| between the S_tot=2 ferromagnet
    /// and the (S_A=1, S_B=1, S_tot=0) singlet realises exactly. Pure-dephasing
    /// dissipator only adds real decay so no L-mode can exceed the H-spread bound.
    /// Tier 1 derived; see <c>PROOF_RING_N4_DIHEDRAL_LOCK.md</c>. Surfaces on the
    /// F1 KB alongside <see cref="KernelDimensionByComponents"/> because the
    /// empirical anchors come from the same Q-sweep that the F1 family scaffolded;
    /// the formula itself is an Im-max bound, conceptually a sister to the star
    /// saturation in <see cref="StarImMaxBound"/>.</summary>
    public RingN4DihedralLockClaim RingN4DihedralLock { get; }

    /// <summary>Star Im-max saturation surfaced by the 2026-05-19 Q-sweep and the
    /// 2026-05-18 SLOW_N8 sweep: Im_max(star, N, J) = (1/2)·J·N Q-universally for
    /// any N ≥ 3. The star Hamiltonian H_star = J·S⃗_0·S⃗_L (hub spin · total leaf
    /// spin) factors through SU(2)/Schur-Weyl Casimir; the maximum-leaf-spin
    /// S_L = (N-1)/2 ferromagnetic sector gives ΔE_max = J·N/2 realised by the
    /// Liouvillian eigenmode between the S_tot = N/2 fully-aligned state and the
    /// S_tot = (N-2)/2 hub-anti-aligned state. Pure-dephasing dissipator only adds
    /// real decay so no L-mode can exceed the H-spread bound. Tier 1 derived; see
    /// <c>PROOF_STAR_OPTICAL_CONFOCAL_SATURATION.md</c>. The Marrakesh-convention
    /// reading <c>Im/σ = 1 ↔ J = 2γ</c> is the Q=2 specialization of the
    /// universal <c>Im/σ = Q/2</c> lock. Sister bound to <see cref="RingN4DihedralLock"/>
    /// via the same Casimir technique (star uses sublattice sizes |A|=1, |B|=N-1;
    /// ring N=4 uses |A|=|B|=2).</summary>
    public StarImMaxBoundClaim StarImMaxBound { get; }

    public IReadOnlyList<HardwareConfirmationClaim> HardwareConfirmations { get; }
    public IReadOnlyList<OpenQuestion> OpenQuestions { get; }

    public F1KnowledgeBase(int N, int? bondCount = null, int? degreeSquaredSum = null)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), $"N must be ≥ 2; got {N}");
        if (bondCount is null != degreeSquaredSum is null)
            throw new ArgumentException(
                "BondCount and DegreeSquaredSum must both be provided for graph-aware scaling, or both omitted for the chain default.");

        this.N = N;
        BondCount = bondCount;
        DegreeSquaredSum = degreeSquaredSum;

        PalindromeIdentity = new F1PalindromeIdentity();
        MainScaling = new PalindromeResidualScalingClaim(N, HamiltonianClass.Main, bondCount, degreeSquaredSum);
        SingleBodyScaling = new PalindromeResidualScalingClaim(N, HamiltonianClass.SingleBody, bondCount, degreeSquaredSum);
        T1ResidualClosedForm = new F1T1ResidualClosedForm();
        T1ResidualPi2Decomposition = new F1T1ResidualPi2Decomposition();
        DepolResidualClosedForm = new F1DepolResidualClosedForm();
        F49NonUniformCrossTerm = new F49NonUniformCrossTermClaim();
        GeneralTopologyVerification = new F1GeneralTopologyVerifiedClaim();
        KernelDimensionByComponents = new F4KernelDimensionByComponentsClaim();
        RingN4DihedralLock = new RingN4DihedralLockClaim();
        StarImMaxBound = new StarImMaxBoundClaim();

        HardwareConfirmations = HardwareConfirmationClaim.LookupAll(_f1ConfirmationNames);

        OpenQuestions = F1OpenQuestions.Standard;
    }

    public string DisplayName =>
        BondCount is { } b
            ? $"F1 knowledge base (N={N}, B={b}, D2={DegreeSquaredSum})"
            : $"F1 knowledge base (N={N}, chain)";

    public string Summary =>
        $"Π·L·Π⁻¹ = −L − 2Σγ·I; main F = {MainScaling.Factor:G6}, single-body F = {SingleBodyScaling.Factor:G6}; " +
        $"{HardwareConfirmations.Count} hardware confirmations, {OpenQuestions.Count} open items";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("N (chain length)",
                summary: $"{N} qubits, {1L << (2 * N)}-dim Pauli-string operator space");

            // Tier 1 (derived) group: F1 master + scaling + closed forms + the F4
            // disconnected-graph bridge (promoted from Tier 1 candidate to Tier 1
            // derived 2026-05-19; closure via DEGENERACY_PALINDROME Result 2 for the
            // connected-case upper bound plus standard tensor-sum kernel factorisation)
            // + the Ring N=4 dihedral lock (Tier 1 derived 2026-05-19; closure via
            // K_{2,2} = C_4 bipartite-complete + Casimir spectrum, see
            // PROOF_RING_N4_DIHEDRAL_LOCK.md) + the Star Im-max saturation
            // (Tier 1 derived 2026-05-19; closure via SU(2)/Schur-Weyl hub-leaf
            // Casimir + maximum-S_L ferromagnet eigenmode, see
            // PROOF_STAR_OPTICAL_CONFOCAL_SATURATION.md).
            yield return InspectableNode.Group("Tier 1 (derived)",
                PalindromeIdentity, MainScaling, SingleBodyScaling,
                T1ResidualClosedForm, T1ResidualPi2Decomposition, DepolResidualClosedForm,
                F49NonUniformCrossTerm, KernelDimensionByComponents, RingN4DihedralLock,
                StarImMaxBound);

            // Tier 2 group includes the hardware confirmations + the
            // general-topology verification record (numerical sweep across N=5..8).
            var tier2Children = new List<IInspectable>(HardwareConfirmations.Count + 1)
            {
                GeneralTopologyVerification,
            };
            tier2Children.AddRange(HardwareConfirmations);
            yield return InspectableNode.Group("Tier 2 (verified)", tier2Children.ToArray());

            yield return InspectableNode.Group("open questions",
                OpenQuestions.Cast<IInspectable>().ToArray());
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
