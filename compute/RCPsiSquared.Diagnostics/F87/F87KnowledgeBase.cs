using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Confirmations;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>"Everything we know about F87" as a typed OOP knowledge graph parametrised by
/// a <see cref="ChainSystem"/>. The root <see cref="IInspectable"/> for <c>--root f87</c>.
///
/// <para>F87 is the Pauli-pair trichotomy classifier (truly / soft / hard). It uses F1's
/// palindrome residual M as its discriminator; see <see cref="F87TrichotomyClassification"/>
/// for the full F-chain layout that hangs off this entry point.</para>
///
/// <para>What is in here:</para>
/// <list type="bullet">
///   <item>Tier-1 derived: <see cref="Trichotomy"/> (the F87 statement itself);
///         <see cref="DissipatorResonance"/> (F87-hardness aligns with dephase-letter
///         Klein index, SU(2)-symmetric); <see cref="DissipatorAxisSelectsPolarity"/>
///         (typed bridge: dissipator letter = polarity-axis selector).</item>
///   <item>Tier-2 empirical: <see cref="CanonicalWitnesses"/>: five canonical Pauli-pair
///         Hamiltonians (XX+YY, Heisenberg, YZ+ZY, XX+XY, XY+YX) classified live on the
///         chain.</item>
///   <item>Tier-2 verified: <see cref="HardwareConfirmations"/>: four Marrakesh entries
///         from <see cref="Core.Confirmations.ConfirmationsRegistry"/>: palindrome trichotomy,
///         F83 4-class signature, Π-protection (YZ+ZY soft), Lebensader skeleton/trace.
///         (Two more, γ_Z calibration and d=0 sector trichotomy, exist only in the Python
///         <c>simulations/framework/confirmations.py</c> and have not been ported to the C#
///         registry yet.)</item>
///   <item>Open: <see cref="OpenQuestions"/>: F80 k≥3, F81/F82/F84 2-qubit dissipators,
///         F83 topology at higher body, per-backend amplification predictor, 4-way Π²-class
///         enum.</item>
/// </list>
///
/// <para>Historical note: pre-2026-05-03 this lived informally as "F77" in the codebase; the
/// formal F77 entry in the registry is unrelated (Multi-drop MM(0) saturation). See memory
/// entry <c>project_f77_f87_rename</c>.</para>
/// </summary>
public sealed class F87KnowledgeBase : IInspectable
{
    private static readonly string[] _f87ConfirmationNames = new[]
    {
        "palindrome_trichotomy",
        "f83_pi2_class_signature_marrakesh",
        "pi_protected_xiz_yzzy",
        "lebensader_skeleton_trace_decoupling",
    };

    public ChainSystem Chain { get; }
    public F87TrichotomyClassification Trichotomy { get; }
    public DissipatorResonanceLaw DissipatorResonance { get; }
    public DissipatorAxisSelectsPolarityClaim DissipatorAxisSelectsPolarity { get; }
    public IReadOnlyList<F87CanonicalWitness> CanonicalWitnesses { get; }
    public IReadOnlyList<HardwareConfirmationClaim> HardwareConfirmations { get; }
    public IReadOnlyList<OpenQuestion> OpenQuestions { get; }

    public F87KnowledgeBase(ChainSystem chain)
    {
        Chain = chain;
        Trichotomy = new F87TrichotomyClassification();
        DissipatorResonance = new DissipatorResonanceLaw();
        DissipatorAxisSelectsPolarity = new DissipatorAxisSelectsPolarityClaim();
        CanonicalWitnesses = F87CanonicalWitness.StandardSet(chain);
        HardwareConfirmations = HardwareConfirmationClaim.LookupAll(_f87ConfirmationNames);
        OpenQuestions = F87OpenQuestions.Standard;
    }

    public string DisplayName => $"F87 knowledge base (N={Chain.N}, J={Chain.J:G3}, γ₀={Chain.GammaZero:G3}, {Chain.Topology})";

    public string Summary =>
        $"Pauli-pair trichotomy (truly/soft/hard) via F1 residual + dissipator-resonance law (F87-hardness ∈ matched Klein cell) + polarity-axis-selector bridge to PolarityLayerOrigin; " +
        $"{CanonicalWitnesses.Count} canonical witnesses, {HardwareConfirmations.Count} hardware confirmations, {OpenQuestions.Count} open items";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("Chain (ChainSystem)",
                summary: $"N={Chain.N}, J={Chain.J:G3}, γ₀={Chain.GammaZero:G3}, topology={Chain.Topology}, H={Chain.HType}");

            yield return InspectableNode.Group("Tier 1 (derived)",
                Trichotomy, DissipatorResonance, DissipatorAxisSelectsPolarity);

            yield return InspectableNode.Group("Tier 2 (empirical canonical witnesses)",
                CanonicalWitnesses.Cast<IInspectable>().ToArray());

            yield return InspectableNode.Group("Tier 2 (hardware-verified)",
                HardwareConfirmations.Cast<IInspectable>().ToArray());

            yield return InspectableNode.Group("open questions",
                OpenQuestions.Cast<IInspectable>().ToArray());
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
