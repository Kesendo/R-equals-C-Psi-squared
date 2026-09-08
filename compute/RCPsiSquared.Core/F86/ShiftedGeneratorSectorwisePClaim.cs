using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86;

/// <summary>The secured F86 symmetry statement for the shifted generator
/// <c>L' = L + Σγ I</c>.  The F1 relation gives <c>Π L' Π^-1 = -L'</c>.
/// On a resolved Π² sector with character <c>p_x ∈ {+1,-1}</c>, the phase-normalized
/// generator <c>P_px = sqrt(p_x) Π</c> is involutive and still anticommutes with L'.
/// This is a sectorwise P-type generator statement, not a global SRP-class assignment.</summary>
public sealed class ShiftedGeneratorSectorwisePClaim : Claim
{
    public ShiftedGeneratorSectorwisePClaim()
        : base("shifted generator has sectorwise phase-normalized P-type anticommutation",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F86A_EP_MECHANISM.md + experiments/PT_SYMMETRY_ANALYSIS.md")
    { }

    public static Complex PhaseForPi2Character(int px) => px switch
    {
        +1 => Complex.One,
        -1 => Complex.ImaginaryOne,
        _ => throw new ArgumentOutOfRangeException(nameof(px), px, "Π² character must be +1 or -1."),
    };

    public static Complex NormalizedGeneratorSquare(int px)
    {
        Complex phase = PhaseForPi2Character(px);
        return phase * phase * px;
    }

    public override string DisplayName =>
        "F86 shifted-generator symmetry: sectorwise involutive P generator";

    public override string Summary =>
        "For L'=L+ΣγI, ΠL'Π^-1=-L'. After resolving Π² character p_x=±1, " +
        "P_px=sqrt(p_x)Π obeys P_px²=I and anticommutes with L'. This secures a " +
        "sectorwise P-type generator; it does not assign the full irreducible SRP class.";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("shifted generator",
                summary: "L' = L + Σγ I; the F1 palindrome becomes Π L' Π^-1 = -L'.");
            yield return new InspectableNode("Π²-sector normalization",
                summary: "p_x=+1 uses P=Π; p_x=-1 uses P=iΠ. In either sector P²=I and {P,L'}=0.");
            yield return new InspectableNode("classification boundary",
                summary: "The full irreducible strong-symmetry-sector algebra has not been resolved; no global SRP class is assigned.");
        }
    }
}

/// <summary>Open classification task downstream of the secured sectorwise P generator.
/// A complete SRP label requires reduction by all strong symmetries and determination of
/// the remaining unitary and antiunitary relations in each irreducible block.</summary>
public sealed class FullIrreducibleSrpClassQuestion : Claim
{
    public FullIrreducibleSrpClassQuestion()
        : base("full irreducible-sector SRP class — OPEN",
               Tier.OpenQuestion,
               "docs/proofs/PROOF_F86A_EP_MECHANISM.md + experiments/PT_SYMMETRY_ANALYSIS.md")
    { }

    public override string DisplayName => "F86 full irreducible-sector SRP class — OPEN";

    public override string Summary =>
        "OpenQuestion: reduce by the complete strong-symmetry algebra, then determine all " +
        "surviving unitary and antiunitary relations per irreducible block. The secured " +
        "sectorwise P-type anticommutation alone does not fix an SRP class.";
}
