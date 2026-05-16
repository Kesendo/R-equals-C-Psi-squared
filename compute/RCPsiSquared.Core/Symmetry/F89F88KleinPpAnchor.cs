using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F89 ↔ F88 bridge: F89's bond Hamiltonian (XX and YY per bond) sits
/// entirely in F88's Klein <b>Pp</b> cell (Π²_Z = Π²_X = +1), and F89's initial state
/// ρ_cc = (|S_1⟩⟨S_2| + |S_2⟩⟨S_1|)/2 is the F88 popcount-pair coherence configuration
/// (n_p=1, n_q=2). Both anchors place F89 within F88's typed taxonomy without ambiguity.
///
/// <para>The two F89-side ingredients (Hamiltonian + initial state) each land in
/// specific F88 cells; together they pin F89's evolution to the Pp/Pm sub-algebra
/// under uniform Z-dephasing:</para>
///
/// <list type="bullet">
///   <item><b>Bond Hamiltonian.</b> XX has bit_a sum 1+1=2 (even → Π²_X=+1) and
///   bit_b sum 0+0=0 (even → Π²_Z=+1). YY has bit_a sum 1+1=2 (even → Π²_X=+1)
///   and bit_b sum 1+1=2 (even → Π²_Z=+1). Both bilinears are <see cref="KleinFourCellClaim"/>
///   Pp-resident; verified via <see cref="KleinBilinearEntry"/>.</item>
///   <item><b>Initial state.</b> ρ_cc is a popcount-(1, 2) coherence pair on N qubits.
///   <see cref="F88PopcountPairLens"/>(N, 1, 2) classifies it:
///   PopcountMirror at N=3 (1+2=3=N → α=0); KIntermediate at N=4 (n_q=2=N/2 →
///   closed-form α); Generic at N≥5 (α=1/2).</item>
/// </list>
///
/// <para>This is the typed-knowledge counterpart to the F88 cross-reference in
/// the F89 entry of <c>docs/ANALYTICAL_FORMULAS.md</c>: F89's structural identity
/// with F88's Pp-Klein-cell + popcount-pair-lens framework, no separate proof needed
/// (both Pi2 and F88 are upstream Tier-1-Derived; this bridge propagates their
/// classifications onto F89's specific operators and state).</para>
///
/// <para><b>Why this matters for cross-claim work</b>: F89's evolution preserves
/// the Pp ⊕ Pm sub-algebra (Z-dephasing closes within Π²_X=+1 cells when the
/// Hamiltonian is Pp-resident); F88's popcount-pair α anchors then directly feed
/// the Π²-odd-in-memory fraction <c>(1/2 − α·s)/(1 − s)</c> for F89's (1, 2)
/// coherence block at each Hamming distance.</para>
///
/// <para>Anchors: <see cref="KleinFourCellClaim"/> (operator-level Pp cell),
/// <see cref="F88PopcountPairLens"/> (state-level α anchors),
/// <see cref="F89TopologyOrbitClosure"/> (F89's bond-graph orbit closure),
/// <c>docs/ANALYTICAL_FORMULAS.md</c> F88 + F89.</para></summary>
public sealed class F89F88KleinPpAnchor : Claim
{
    public const string F89BondKleinCell = "Pp";
    public const int F89BondPi2ZEigenvalue = +1;
    public const int F89BondPi2XEigenvalue = +1;

    public const int F89InitialStateNp = 1;
    public const int F89InitialStateNq = 2;

    private readonly KleinFourCellClaim _kleinFourCell;
    private readonly F89TopologyOrbitClosure _f89;

    /// <summary>Cell label for F89's XX bond term (always Pp; computed from
    /// <see cref="KleinBilinearEntry"/>).</summary>
    public static KleinBilinearEntry F89XxEntry { get; } =
        new KleinBilinearEntry(PauliLetter.X, PauliLetter.X);

    /// <summary>Cell label for F89's YY bond term (always Pp; computed from
    /// <see cref="KleinBilinearEntry"/>).</summary>
    public static KleinBilinearEntry F89YyEntry { get; } =
        new KleinBilinearEntry(PauliLetter.Y, PauliLetter.Y);

    /// <summary>F88 popcount-pair lens for F89's (S_1, S_2) Dicke-coherence initial
    /// state on an N-qubit chain. Returns the typed Tier-1-Derived F88 lens with
    /// its closed-form α anchor.</summary>
    public static F88PopcountPairLens InitialStateLens(int N) =>
        new F88PopcountPairLens(N, F89InitialStateNp, F89InitialStateNq);

    public F89F88KleinPpAnchor(KleinFourCellClaim kleinFourCell, F89TopologyOrbitClosure f89)
        : base("F89 ↔ F88 bridge: F89 bond Hamiltonian (XX, YY) is F88-Pp-Klein-cell-resident; ρ_cc initial state is the F88 popcount-(1, 2) coherence pair, classified by F88PopcountPairLens(N, 1, 2) as PopcountMirror at N=3, KIntermediate at N=4, Generic at N≥5",
               Tier.Tier1Derived,
               "KleinFourCellClaim Pp cell (operator-level), F88PopcountPairLens (state-level α anchors), F89TopologyOrbitClosure (F89 bond-graph orbit closure), docs/ANALYTICAL_FORMULAS.md F88 + F89")
    {
        _kleinFourCell = kleinFourCell ?? throw new ArgumentNullException(nameof(kleinFourCell));
        _f89 = f89 ?? throw new ArgumentNullException(nameof(f89));
    }

    public override string DisplayName => "F89 bond (Pp) + ρ_cc (popcount-(1, 2)) anchor in F88";

    public override string Summary =>
        $"XX, YY both F88-Pp (Π²_Z = Π²_X = +1); ρ_cc = popcount-(1, 2) coherence; F88PopcountPairLens(N, 1, 2) " +
        $"is PopcountMirror at N=3, KIntermediate at N=4, Generic at N≥5 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F89 XX bond",
                summary: $"{F89XxEntry.DisplayName} — {F89XxEntry.Summary}");
            yield return new InspectableNode("F89 YY bond",
                summary: $"{F89YyEntry.DisplayName} — {F89YyEntry.Summary}");
            yield return new InspectableNode("F89 initial state ρ_cc",
                summary: $"(|S_{F89InitialStateNp}⟩⟨S_{F89InitialStateNq}| + |S_{F89InitialStateNq}⟩⟨S_{F89InitialStateNp}|)/2; F88 popcount-pair coherence with (n_p, n_q) = ({F89InitialStateNp}, {F89InitialStateNq})");
            yield return new InspectableNode("F88 lens classification per N",
                summary: "N=3 → PopcountMirror (1+2=3=N, α=0); N=4 → KIntermediate (n_q=2=N/2, closed-form α); N≥5 → Generic (α=1/2)");
            yield return new InspectableNode("Sub-algebra closure",
                summary: "Pp-resident H + Z-dephasing closes within Pp ⊕ Pm; F89 evolution stays in the Π²_X=+1 axis throughout");
            yield return new InspectableNode("Cross-bridge to F89 AT-lock",
                summary: "Pp Klein-cell residency = subset of F87-Truly; combined with F89F87TrulyInheritance the F89 AT-lock Re(λ_n) = −2γ₀ follows");
        }
    }
}
