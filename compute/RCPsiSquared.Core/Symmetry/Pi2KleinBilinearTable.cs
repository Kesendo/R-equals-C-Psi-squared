using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The 9-Pauli-pair-bilinear × 4-Klein-cell assignment table (Tier 2 empirical).
/// Self-computing: for each of the 9 non-identity 2-body bilinears (P, Q) ∈ {X, Y, Z}²,
/// the cell is determined by (Π²_Z, Π²_X) eigenvalues from <see cref="PiOperator.SquaredEigenvalue"/>
/// applied to the letter pair.</summary>
public sealed class Pi2KleinBilinearTable : Claim
{
    public IReadOnlyList<KleinBilinearEntry> Entries { get; }

    public Pi2KleinBilinearTable()
        : base("9 Pauli-pair bilinears × 4 Klein cells",
               Tier.Tier2Empirical,
               "docs/ANALYTICAL_FORMULAS.md F88 + Pi2KleinViewTests")
    {
        var letters = new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.Z };
        var entries = new List<KleinBilinearEntry>();
        foreach (var p in letters)
            foreach (var q in letters)
                entries.Add(new KleinBilinearEntry(p, q));
        Entries = entries;
    }

    public override string DisplayName => "Klein bilinear table (9 × 4 cells)";

    public override string Summary
    {
        get
        {
            var counts = Entries.GroupBy(e => e.Cell).ToDictionary(g => g.Key, g => g.Count());
            return $"Pp={counts.GetValueOrDefault("Pp", 0)} Pm={counts.GetValueOrDefault("Pm", 0)} " +
                   $"Mp={counts.GetValueOrDefault("Mp", 0)} Mm={counts.GetValueOrDefault("Mm", 0)}";
        }
    }

    protected override IEnumerable<IInspectable> ExtraChildren =>
        Entries.Cast<IInspectable>();
}

public sealed class KleinBilinearEntry : IInspectable
{
    public PauliLetter LetterA { get; }
    public PauliLetter LetterB { get; }
    public int Pi2ZEigenvalue { get; }
    public int Pi2XEigenvalue { get; }

    public KleinBilinearEntry(PauliLetter a, PauliLetter b)
    {
        LetterA = a;
        LetterB = b;
        var pair = new[] { a, b };
        Pi2ZEigenvalue = PiOperator.SquaredEigenvalue(pair, PauliLetter.Z);
        Pi2XEigenvalue = PiOperator.SquaredEigenvalue(pair, PauliLetter.X);
    }

    public string Cell => (Pi2ZEigenvalue, Pi2XEigenvalue) switch
    {
        (+1, +1) => "Pp",
        (+1, -1) => "Pm",
        (-1, +1) => "Mp",
        (-1, -1) => "Mm",
        _ => "??",
    };

    public string DisplayName => $"{LetterA}{LetterB} → {Cell} ({Pi2ZSign}{Pi2XSign})";

    private string Pi2ZSign => Pi2ZEigenvalue == +1 ? "+" : "−";
    private string Pi2XSign => Pi2XEigenvalue == +1 ? "+" : "−";

    public string Summary
    {
        get
        {
            bool truly = (LetterA, LetterB) is (PauliLetter.X, PauliLetter.X)
                or (PauliLetter.Y, PauliLetter.Y)
                or (PauliLetter.Z, PauliLetter.Z);
            return $"Π²_Z = {Pi2ZSign}1, Π²_X = {Pi2XSign}1{(truly ? "; truly" : "")}";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("Π²_Z eigenvalue", summary: Pi2ZEigenvalue.ToString("+#;−#"));
            yield return new InspectableNode("Π²_X eigenvalue", summary: Pi2XEigenvalue.ToString("+#;−#"));
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
