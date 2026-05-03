using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Confirmations;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>"Everything we know about the Π² involution and its Klein decomposition" as a
/// typed OOP knowledge graph parametrised by a <see cref="ChainSystem"/>. The root
/// <see cref="IInspectable"/> for <c>--root pi2</c>.
///
/// <para>The Π² Klein layer is the foundational structural primitive that underlies the
/// F-chain (F1 palindrome → F87 trichotomy → F80–F85 dissipator analysis): Π² acts
/// diagonally on Pauli strings as ±1, and its (Π²_Z, Π²_X) two-axis decomposition
/// classifies the 4^N Pauli operator space into 4 Klein cells (F88).</para>
///
/// <para>What's in here:</para>
/// <list type="bullet">
///   <item>Tier-1 derived: <see cref="PolynomialFoundation"/> (d²−2d=0 ↔ R=CΨ², the trunk),
///         <see cref="RootAnchor"/> (1/d = 1/2 lineage), <see cref="Involution"/>,
///         <see cref="KleinDecomposition"/>, <see cref="BilinearApex"/>, <see cref="MirrorRegime"/>,
///         <see cref="HalfFixedPoint"/> (three faces close), <see cref="MirrorMemory"/>
///         (90° back to the mirror, F80's i)</item>
///   <item>Tier-2 empirical: <see cref="BilinearTable"/> (9 Pauli-pair × 4 cells)</item>
///   <item>Tier-2 hardware-verified: <see cref="HardwareConfirmations"/> (Marrakesh f83
///         X-axis-flip pattern)</item>
///   <item>Open: <see cref="OpenQuestions"/> (X-flip mechanism, 2:2 truly-kernel,
///         N≥4 transition, k-body extension, mirror-regime relation)</item>
/// </list>
///
/// <para>Schicht-Hierarchie der Π²-Primitive (für Implementierung):</para>
/// <list type="number">
///   <item>Schicht 1 (one axis): <see cref="Pi2Projection.SymmetricProjector"/> +
///         <see cref="Pi2Projection.AntisymmetricProjector"/> + <see cref="Pi2Projection.Split"/></item>
///   <item>Schicht 1.5 (typed enum): <see cref="Pauli.Pi2Class"/> +
///         <see cref="Pauli.PauliPairBondTermExtensions.Pi2ClassOf"/></item>
///   <item>Schicht 2 (two axes): <see cref="Pi2Projection.KleinSplit"/> +
///         <see cref="KleinDecomposition"/> record</item>
///   <item>Schicht 3 (Klein × spectrum): <see cref="Pi2KleinSpectralView.ComputeFor"/></item>
/// </list>
/// </summary>
public sealed class Pi2KnowledgeBase : IInspectable
{
    private static readonly string[] _hardwareConfirmationNames = new[]
    {
        "f83_pi2_class_signature_marrakesh",
        "palindrome_trichotomy",
    };

    public ChainSystem Chain { get; }
    public PolynomialFoundationClaim PolynomialFoundation { get; }
    public QubitDimensionalAnchorClaim RootAnchor { get; }
    public Pi2InvolutionClaim Involution { get; }
    public KleinFourCellClaim KleinDecomposition { get; }
    public BilinearApexClaim BilinearApex { get; }
    public HalfIntegerMirrorClaim MirrorRegime { get; }
    public HalfAsStructuralFixedPointClaim HalfFixedPoint { get; }
    public NinetyDegreeMirrorMemoryClaim MirrorMemory { get; }
    public Pi2KleinBilinearTable BilinearTable { get; }
    public IReadOnlyList<HardwareConfirmationClaim> HardwareConfirmations { get; }
    public IReadOnlyList<OpenQuestion> OpenQuestions { get; }

    public Pi2KnowledgeBase(ChainSystem chain)
    {
        Chain = chain;
        PolynomialFoundation = new PolynomialFoundationClaim();
        RootAnchor = new QubitDimensionalAnchorClaim();
        Involution = new Pi2InvolutionClaim();
        KleinDecomposition = new KleinFourCellClaim();
        BilinearApex = new BilinearApexClaim();
        MirrorRegime = new HalfIntegerMirrorClaim(chain.N);
        HalfFixedPoint = new HalfAsStructuralFixedPointClaim();
        MirrorMemory = new NinetyDegreeMirrorMemoryClaim();
        BilinearTable = new Pi2KleinBilinearTable();
        HardwareConfirmations = HardwareConfirmationClaim.LookupAll(_hardwareConfirmationNames);
        OpenQuestions = Pi2OpenQuestions.Standard;
    }

    public string DisplayName =>
        $"Π² knowledge base (N={Chain.N}, w_XY={Chain.N / 2.0:F1}, " +
        $"{(Chain.N % 2 == 1 ? "half-integer" : "integer")}-mirror)";

    public string Summary =>
        $"Π² involution + Klein 4-cell decomposition (F88); apex 1/2 anchor; " +
        $"{HardwareConfirmations.Count} hardware confirmations; {OpenQuestions.Count} open items";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("Chain (ChainSystem)",
                summary: $"N={Chain.N}, J={Chain.J:G3}, γ₀={Chain.GammaZero:G3}, " +
                         $"topology={Chain.Topology}, H={Chain.HType}");

            yield return InspectableNode.Group("Tier 1 (derived)",
                PolynomialFoundation, RootAnchor, Involution, KleinDecomposition, BilinearApex, MirrorRegime, HalfFixedPoint, MirrorMemory);

            yield return InspectableNode.Group("Tier 2 (empirical)",
                BilinearTable);

            yield return InspectableNode.Group("Tier 2 (hardware-verified)",
                HardwareConfirmations.Cast<IInspectable>().ToArray());

            yield return InspectableNode.Group("open questions",
                OpenQuestions.Cast<IInspectable>().ToArray());
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
