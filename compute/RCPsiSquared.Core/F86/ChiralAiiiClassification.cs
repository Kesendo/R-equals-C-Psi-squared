using RCPsiSquared.Core.Inspection;

using RCPsiSquared.Core.Knowledge;
namespace RCPsiSquared.Core.F86;

/// <summary>F86 algebraic-classification claim: the conjugation operator Π places the
/// <b>centered</b> Liouvillian L_c (the F1 palindrome Π·L·Π⁻¹ + L + 2Σγ·I = 0, so the spectrum
/// is symmetric about −Σγ) in <b>class AIII chiral</b>, distinct from Bender-Boettcher PT
/// (Π is linear, order 4; classical PT requires anti-linear operators). This is a general
/// property of the full operator structure, established in
/// <c>experiments/PT_SYMMETRY_ANALYSIS.md</c>.
///
/// <para>The chiral defining property is λ ↔ −λ pairing. It is <b>exact</b> in the Σγ=0
/// gain-loss system of <c>hypotheses/FRAGILE_BRIDGE.md</c> (centre at 0): the Hopf bifurcation
/// there IS the chiral-symmetry-breaking event, with Petermann factor K=403 signalling a
/// genuine defective EP in the complex γ plane. The toy 2-level EP at Q_EP = 2/g_eff is a
/// <b>separate object</b>: its eigenvalues −4γ₀ ± √(4γ₀²−J²·g_eff²) coalesce on the real axis
/// at Re=−4γ₀ (centre −Nγ₀, NOT 0), so in isolation it carries no λ ↔ −λ pairing; it is
/// PT-like phenomenology (EP at finite coupling, spectral flow), not an instance OF the chiral
/// classification. These are the two firmly-established genuine defective EPs (the Σγ=0
/// gain-loss EP and the toy 2×2), distinct mechanisms that happen to share the same-sign-
/// imaginary 2×2 algebraic form.</para>
///
/// <para>Whether the toy 2×2's coalescence reflects a defective-EP / exact-chiral-pairing
/// structure shared by the FULL Σγ=N·γ₀ coherence block (and hence linked to FRAGILE_BRIDGE)
/// is an OPEN question (<see cref="LocalGlobalEpLink"/>, demoted Tier2Verified → OpenQuestion
/// by the F86a-retraction 2026-06-21): the full block is genuinely non-normal on the real Q
/// axis but has NO eigenvalue coalescence there (eigenvalues simple). This is open at the
/// EP-link level only; the AIII chiral classification of the centered L_c itself is general
/// (via F1, for every Σγ) and is not in question. The surviving firmly-shared content is
/// exactly that AIII chiral algebraic label.</para>
/// </summary>
public sealed class ChiralAiiiClassification : Claim
{
    public ChiralAiiiClassification()
        : base("class AIII chiral (NOT Bender-Boettcher PT)",
               Tier.Tier1Derived,
               "experiments/PT_SYMMETRY_ANALYSIS.md + hypotheses/FRAGILE_BRIDGE.md + docs/ANALYTICAL_FORMULAS.md F86")
    { }

    public override string DisplayName => "F86 algebraic class: AIII chiral, NOT Bender-Boettcher PT";

    public override string Summary =>
        $"Π places the centered L_c (F1 palindrome) in class AIII chiral (Π linear, order 4), NOT classical PT (anti-linear required); the same-sign +iJ·g_eff toy EP is PT-like phenomenology, a separate object";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("toy 2×2 EP (separate object, NOT an AIII instance)",
                summary: "2-level rate-channel EP at Q_EP = 2/g_eff (F86 Statement 1); coalesces on the real axis at Re=−4γ₀ (centre −Nγ₀), no λ↔−λ pairing in isolation; PT-like phenomenology");
            yield return new InspectableNode("Σγ=0 chiral-exact EP",
                summary: "FRAGILE_BRIDGE gain-loss (centre 0): exact λ↔−λ pairing, the Hopf bifurcation IS chiral-symmetry breaking; Petermann K=403, genuine defective EP in complex γ plane");
            yield return new InspectableNode("classification anchor",
                summary: "experiments/PT_SYMMETRY_ANALYSIS.md: Π places the centered L_c in class AIII chiral (linear, order 4)");
            yield return new InspectableNode("PTF c=1 sibling",
                summary: "F86's c≥2 EP machinery generalises PTF's c=1 (vac-SE block) α_i closure law");
        }
    }
}
