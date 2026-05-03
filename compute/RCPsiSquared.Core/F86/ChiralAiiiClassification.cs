using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 algebraic-classification claim: the same-sign-imaginary 2-level form is
/// "PT-phenomenology-like" (EP at finite coupling, spectral flow) but algebraically lives
/// in <b>class AIII chiral</b> — distinct from Bender-Boettcher PT (Π is linear; classical
/// PT requires anti-linear operators).
///
/// <para>The local 2-level EP at Q_EP = 2/g_eff is the rate-channel instance of the chiral
/// classification established for the full Liouvillian; the Hopf bifurcation in
/// <c>hypotheses/FRAGILE_BRIDGE.md</c> is the global instance with Petermann factor K=403
/// signaling an EP in the complex γ plane. Whether the local-2-level-EP and the
/// global-complex-γ-EP are connected algebraically is itself open.</para>
/// </summary>
public sealed class ChiralAiiiClassification : F86Claim
{
    public ChiralAiiiClassification()
        : base("class AIII chiral (NOT Bender-Boettcher PT)",
               Tier.Tier1Derived,
               "experiments/PT_SYMMETRY_ANALYSIS.md + hypotheses/FRAGILE_BRIDGE.md + docs/ANALYTICAL_FORMULAS.md F86")
    { }

    public override string DisplayName => "F86 algebraic class: AIII chiral, NOT Bender-Boettcher PT";

    public override string Summary =>
        $"same-sign +iJ·g_eff form admits EP at finite coupling (Tier 1); algebraically chiral (Π linear), NOT classical PT (anti-linear required)";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("local instance",
                summary: "2-level rate-channel EP at Q_EP = 2/g_eff (F86 Statement 1)");
            yield return new InspectableNode("global instance",
                summary: "Hopf bifurcation in FRAGILE_BRIDGE, Petermann K=403 in complex γ plane");
            yield return new InspectableNode("classification anchor",
                summary: "experiments/PT_SYMMETRY_ANALYSIS.md — Π is class AIII chiral (linear)");
            yield return new InspectableNode("PTF c=1 sibling",
                summary: "F86's c≥2 EP machinery generalises PTF's c=1 (vac-SE block) α_i closure law");
        }
    }
}
