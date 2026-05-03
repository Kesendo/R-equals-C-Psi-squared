using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F71;

/// <summary>"Everything we know about F71" as a typed OOP knowledge graph parametrised by
/// chain length N. The root <see cref="IInspectable"/> for <c>--root f71</c>.
///
/// <para>F71 is purely kinematic (dimension-only), so the KB depends on N alone; no
/// chromaticity / γ₀ required. The chain-mirror operator R, bond-orbit decomposition, and
/// the c₁ / Q_peak mirror identities are all N-parametrised.</para>
///
/// <para>What is in here:</para>
/// <list type="bullet">
///   <item>Tier-1 derived: <see cref="MirrorOperator"/> (R, R²=I), <see cref="BondOrbits"/>
///         (orbit decomposition), <see cref="C1Identity"/> (c₁ closed-form mirror identity),
///         <see cref="F86Generalisation"/> (Q_peak mirror cross-reference to F86 KB).</item>
///   <item>Open: <see cref="OpenQuestions"/>: non-uniform J/γ extension, asymmetric ρ₀,
///         per-orbit Q_peak derivation.</item>
/// </list>
///
/// <para>The F71 → F86 link is a Tier-1 derived cross-reference, not a duplicate of the
/// F86 KB's own F71MirrorInvariance claim. The F86 KB carries the empirical F86 substructure
/// (witness data, per-orbit observation, etc.); the F71 KB carries the kinematic argument
/// and points to F86 as one of its consequences.</para>
/// </summary>
public sealed class F71KnowledgeBase : IInspectable
{
    public int N { get; }
    public F71MirrorOperator MirrorOperator { get; }
    public F71BondOrbitDecomposition BondOrbits { get; }
    public C1MirrorIdentity C1Identity { get; }
    public F86MirrorGeneralisationLink F86Generalisation { get; }
    public IReadOnlyList<OpenQuestion> OpenQuestions { get; }

    public F71KnowledgeBase(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), $"N must be ≥ 2; got {N}");
        this.N = N;
        MirrorOperator = new F71MirrorOperator(N);
        BondOrbits = new F71BondOrbitDecomposition(N);
        C1Identity = new C1MirrorIdentity();
        F86Generalisation = new F86MirrorGeneralisationLink();
        OpenQuestions = F71OpenQuestions.Standard;
    }

    public string DisplayName => $"F71 knowledge base (N={N})";

    public string Summary =>
        $"chain-mirror R on {1 << N}-dim Hilbert space; {BondOrbits.NumOrbits} bond orbits " +
        $"({(BondOrbits.HasSelfPairedCentralOrbit ? "self-paired central present" : "no self-paired central")}); " +
        $"c₁ identity + F86 generalisation, {OpenQuestions.Count} open items";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("N (chain length)",
                summary: $"{N} qubits, {N - 1} bonds, {1 << N}-dim Hilbert space");

            yield return InspectableNode.Group("Tier 1 (derived)",
                MirrorOperator, BondOrbits, C1Identity, F86Generalisation);

            yield return InspectableNode.Group("open questions",
                OpenQuestions.Cast<IInspectable>().ToArray());
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
