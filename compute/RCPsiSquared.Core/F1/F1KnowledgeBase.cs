using RCPsiSquared.Core.Confirmations;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;

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
///         graphs at N=5, 6, 7, 8 (N=8 via the opt-in SLOW_N8 block-spectrum dogfood
///         with full <see cref="F1SpectrumStatistics"/> metric capture under
///         <c>simulations/results/f1_n8_n9_metrics/</c>). N=9 is wired as a
///         <see cref="Xunit.SkippableFactAttribute"/> but blocked at the LP64 MKL
///         P/Invoke marshalling ceiling; see
///         <see cref="F1GeneralTopologyVerifiedClaim.ScaleFrontierBlockedAtN"/>.</item>
///   <item>Open: <see cref="OpenQuestions"/> is EMPTY as of 2026-05-18: the first
///         time the F1 family is open-question-free. All four items from the May 2026
///         sprint closed: T1 closed form (Tier-1 derived), depol closed form (Tier-1
///         derived), non-uniform γ (negative-result closure), and general topology
///         (synthesis proof + verification record). See <see cref="F1OpenQuestions"/>
///         XML doc for the per-item closure references.</item>
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

            yield return InspectableNode.Group("Tier 1 (derived)",
                PalindromeIdentity, MainScaling, SingleBodyScaling,
                T1ResidualClosedForm, T1ResidualPi2Decomposition, DepolResidualClosedForm,
                F49NonUniformCrossTerm);

            // Tier 2 group includes the hardware confirmations + the
            // general-topology verification record (numerical sweep across N=5..7).
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
