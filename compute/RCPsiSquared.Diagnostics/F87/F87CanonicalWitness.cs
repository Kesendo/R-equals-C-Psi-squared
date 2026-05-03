using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>An F87 trichotomy witness: a named canonical Pauli-pair Hamiltonian on a chain
/// that exercises one branch of the trichotomy. Self-computes its classification on demand
/// via <see cref="PauliPairTrichotomy.Classify"/>, lazily, so listing many witnesses doesn't
/// trigger N L-builds upfront.
///
/// <para>The expected <see cref="ExpectedClass"/> is the analytical / empirical anchor; the
/// <see cref="ActualClass"/> is computed on the chain. Disagreement signals either a chain
/// configuration that breaks the witness (e.g. wrong N, T1 added, non-uniform γ) or a real
/// regression in the classifier.</para>
/// </summary>
public sealed class F87CanonicalWitness : Claim
{
    private readonly Lazy<TrichotomyClass> _actual;

    public ChainSystem Chain { get; }
    public string Description { get; }
    public IReadOnlyList<PauliPairBondTerm> Terms { get; }
    public TrichotomyClass ExpectedClass { get; }

    public TrichotomyClass ActualClass => _actual.Value;
    public bool Matches => ActualClass == ExpectedClass;

    public F87CanonicalWitness(string name, string description, ChainSystem chain,
        IReadOnlyList<PauliPairBondTerm> terms, TrichotomyClass expectedClass, string anchor)
        : base($"F87 witness: {name}", Tier.Tier2Empirical, anchor)
    {
        Chain = chain;
        Description = description;
        Terms = terms;
        ExpectedClass = expectedClass;
        _actual = new Lazy<TrichotomyClass>(() => PauliPairTrichotomy.Classify(chain, terms));
    }

    public override string DisplayName => $"F87 witness: {Name.Replace("F87 witness: ", "")} (expect {ExpectedClass})";

    public override string Summary =>
        $"chain N={Chain.N}, terms [{TermsLabel()}] → expected {ExpectedClass}, actual {ActualClass} ({(Matches ? "match" : "MISMATCH")})";

    private string TermsLabel() =>
        string.Join(" + ", Terms.Select(t => $"{t.LetterA}{t.LetterB}"));

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("description", summary: Description);
            yield return new InspectableNode("terms", summary: TermsLabel());
            yield return new InspectableNode("expected", summary: ExpectedClass.ToString());
            yield return new InspectableNode("actual", summary: ActualClass.ToString());
            yield return new InspectableNode("agreement", summary: Matches ? "MATCH" : "MISMATCH");
        }
    }

    /// <summary>Five canonical witnesses spanning all three trichotomy branches plus the
    /// EQ-030 Marrakesh-confirmed soft anchor.</summary>
    public static IReadOnlyList<F87CanonicalWitness> StandardSet(ChainSystem chain) => new[]
    {
        new F87CanonicalWitness(
            name: "XX+YY",
            description: "Canonical XY chain. Π·L·Π⁻¹ + L + 2σ·I = 0 holds bit-exactly.",
            chain: chain,
            terms: new[]
            {
                new PauliPairBondTerm(PauliLetter.X, PauliLetter.X),
                new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Y),
            },
            expectedClass: TrichotomyClass.Truly,
            anchor: "docs/ANALYTICAL_FORMULAS.md F1 + F87"),
        new F87CanonicalWitness(
            name: "XX+YY+ZZ (Heisenberg)",
            description: "Heisenberg chain: both XX+YY and ZZ are Π²-even truly; Z-dephasing leaves M=0.",
            chain: chain,
            terms: new[]
            {
                new PauliPairBondTerm(PauliLetter.X, PauliLetter.X),
                new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Y),
                new PauliPairBondTerm(PauliLetter.Z, PauliLetter.Z),
            },
            expectedClass: TrichotomyClass.Truly,
            anchor: "docs/ANALYTICAL_FORMULAS.md F1 + F87"),
        new F87CanonicalWitness(
            name: "YZ+ZY (EQ-030 soft)",
            description: "Bond-flipped YZ: bit_b sum is even but the Hamiltonian is not truly. Marrakesh-confirmed soft Hamiltonian (drop=28, EQ-030).",
            chain: chain,
            terms: new[]
            {
                new PauliPairBondTerm(PauliLetter.Y, PauliLetter.Z),
                new PauliPairBondTerm(PauliLetter.Z, PauliLetter.Y),
            },
            expectedClass: TrichotomyClass.Soft,
            anchor: "review/EMERGING_QUESTIONS.md EQ-030 + Marrakesh d7n3013aq2pc73a2a18g"),
        new F87CanonicalWitness(
            name: "XX+XY (mixed hard)",
            description: "XX (Π²-even, truly) + XY (Π²-odd): mixed Π²-class, no spectrum pairing. The F87 hard anchor.",
            chain: chain,
            terms: new[]
            {
                new PauliPairBondTerm(PauliLetter.X, PauliLetter.X),
                new PauliPairBondTerm(PauliLetter.X, PauliLetter.Y),
            },
            expectedClass: TrichotomyClass.Hard,
            anchor: "docs/ANALYTICAL_FORMULAS.md F87 + experiments/V_EFFECT_FINE_STRUCTURE.md"),
        new F87CanonicalWitness(
            name: "XY+YX (bond-flip soft)",
            description: "Bond-flipped Z-free pair: soft per Marrakesh skeleton-trace test. Both Π²-odd.",
            chain: chain,
            terms: new[]
            {
                new PauliPairBondTerm(PauliLetter.X, PauliLetter.Y),
                new PauliPairBondTerm(PauliLetter.Y, PauliLetter.X),
            },
            expectedClass: TrichotomyClass.Soft,
            anchor: "Marrakesh lebensader_skeleton_trace_decoupling"),
    };
}
