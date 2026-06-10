using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>The Tier1Derived two-reflection spine of the F87 windowed-converse monomial theorem
/// (Phase B, 2026-06-09). Recenter the Liouvillian as M = A + γQ on the coherence space, A = −i[H,·]
/// = A_L + A_R (left/right hop), Q = Σ_l Z_l⊗Z_l diagonal. Two involutions, 𝓕 = F⊗F and R = I⊗F
/// (F = X^⊗N), give the sign table 𝓕A𝓕 = −A, 𝓕Q𝓕 = +Q, R A_L R = +A_L, R A_R R = −A_R, R Q R = −Q
/// (holding for complex H via FH^TF = (FHF)^T = −H^T). Hence every word contributing to an ODD
/// power-sum has #A_L, #A_R, #Q all odd. Since a nonzero diagonal trace needs the left/right hops to
/// be closed walks on the hopping graph, and an odd closed walk exists iff non-bipartite (minimal
/// length the unsigned odd-girth ℓ):
///
/// <list type="bullet">
///   <item>bipartite ⟹ no odd closed walk ⟹ all odd power-sums ≡ 0 (SOFT), a second proof of
///         bipartite ⟹ soft, complementary to <see cref="Core.Symmetry.ChiralKClaim"/>;</item>
///   <item>non-bipartite ⟹ #A_L, #A_R ≥ ℓ ⟹ #A ≥ 2ℓ (THRESHOLD): the first asymmetric odd moment
///         sits at the odd-girth out-and-back;</item>
///   <item>deg-1 positivity: P_{3,1} = 6·4^N·Σ_l c_l² over H's single-site-Z Pauli coefficients
///         (tensor traces: Tr(A²Q) = 2·Σ_l Tr(HZ_l)²), manifestly ≥ 0 and > 0 exactly for a
///         single-site-Z lift.</item>
/// </list>
///
/// <para>Tier1Derived: this spine is proven general-N, including (since 2026-06-10) the girth
/// dichotomy: the deg-1 class in closed form at every m, with hard-at-all-γ OUTRIGHT on the
/// t_ℓ ≠ 0 branch. On the t_ℓ = 0 branch the monomial form is proven and only R-sign is open,
/// carried by <see cref="WindowedConverseAllGammaClaim"/>, so nothing open is typed as Derived
/// here. Anchor: <c>docs/proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md</c> +
/// <c>simulations/f87_windowed_monomial_converse.py</c> +
/// <c>simulations/f87_girth_dichotomy.py</c>.</para></summary>
public sealed class WindowedConverseThresholdClaim : Claim
{
    /// <summary>One self-check tying the claim to the GF(2) odd-girth / m* arithmetic.</summary>
    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    public IReadOnlyList<BatteryCase> Cases { get; }
    public int PassCount => Cases.Count(c => c.Passes);

    public WindowedConverseThresholdClaim()
        : base("F87 windowed-converse threshold (two-reflection spine): all-odd (#A_L,#A_R,#Q) parity ⟹ bipartite soft + non-bipartite #A≥2ℓ; the deg-1 class in closed form (girth dichotomy). Tier1Derived; hard-at-all-γ outright on the t_ℓ≠0 branch, modulo R-sign on the t_ℓ=0 branch",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md + " +
               "simulations/f87_windowed_monomial_converse.py")
    {
        Cases = BuildBattery();
    }

    public string TwoReflections =>
        "𝓕 = F⊗F, R = I⊗F (F = X^⊗N): 𝓕A𝓕 = −A, 𝓕Q𝓕 = +Q, R·A_L·R = +A_L, R·A_R·R = −A_R, R·Q·R = −Q " +
        "(complex H included, via FH^TF = (FHF)^T = −H^T). For odd power-sums this forces #A_L, #A_R, #Q all odd.";

    public string Threshold =>
        "A nonzero diagonal trace needs the left/right hops to be closed walks on the hopping graph; an odd " +
        "closed walk exists iff non-bipartite (minimal length the unsigned odd-girth ℓ). So bipartite ⟹ all odd " +
        "power-sums ≡ 0 (soft); non-bipartite ⟹ #A_L,#A_R ≥ ℓ ⟹ #A ≥ 2ℓ (first asymmetry at the odd-girth).";

    public string Deg1Positivity =>
        "deg-1 in closed form at every m (the girth dichotomy, 2026-06-10): the supertrace factorization " +
        "Tr(Q·A^{2k}) = (−1)^k Σ_l Σ_j (−1)^j C(2k,j)·t_j·t_{2k−j} with t_j = Tr(Z_l H^j); F-chirality kills " +
        "even j, the unsigned odd-girth kills j < ℓ, and P_{2ℓ+1,1} = (2ℓ+1)·C(2ℓ,ℓ)·Σ_l t_ℓ² ≥ 0 (a sum of " +
        "squares). The ℓ=1 face is P_{3,1} = 6·4^N·Σ_l c_l² (tensor traces: Tr(A²Q) = 2·Σ_l Tr(HZ_l)²). " +
        "Closed form, RIGOROUS-GENERAL; anchor simulations/f87_girth_dichotomy.py.";

    public string FollowsModulo =>
        "The girth ladder (2026-06-10): t_ℓ ≠ 0 ⟹ m* = 2ℓ+1, deg 1, positive outright (hard ∀γ>0, no " +
        "residual); t_ℓ ≡ 0 ⟹ the γ¹ class is proven dead at 2ℓ+1 and 2ℓ+3 and the first firing odd rung " +
        "sets m* = 2ℓ+deg (γ³ when it fires; higher when not, first γ⁵ witness IIXY+ZXZY), with R-sign in " +
        "ladder form (first surviving class single + positive) open, carried by WindowedConverseAllGammaClaim.";

    public override string DisplayName =>
        "F87 windowed-converse threshold (two-reflection spine, Tier1Derived)";

    public override string Summary =>
        "two reflections force #A_L,#A_R,#Q all odd in every odd power-sum word; bipartite ⟹ all odd power-sums " +
        "≡ 0 (soft), non-bipartite ⟹ #A ≥ 2ℓ, deg-1 positivity closed-form; monomial/all-γ-hard follows once " +
        $"R-deg + R-sign hold; {PassCount}/{Cases.Count} PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Two reflections (𝓕, R) ⟹ all-odd parity", summary: TwoReflections);
            yield return new InspectableNode("Threshold #A ≥ 2ℓ + soft re-proof", summary: Threshold);
            yield return new InspectableNode("deg-1 positivity (closed form)", summary: Deg1Positivity);
            yield return new InspectableNode("Follows once R-deg + R-sign hold", summary: FollowsModulo);
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
        }
    }

    /// <summary>Light pure-GF(2) self-check (no Liouvillian): the canonical triangle's odd-girth ℓ
    /// (an odd 3-relation among the masks {0b011, 0b110, 0b101}) is 3, and m* = 2ℓ + deg with deg = 3
    /// is 9, the structural numbers the proof and the verification script certify.</summary>
    private static IReadOnlyList<BatteryCase> BuildBattery()
    {
        var cases = new List<BatteryCase>();
        int ell = WindowedObstructionScan.MinOddCycle(new ulong[] { 0b011, 0b110, 0b101 });
        cases.Add(new BatteryCase(
            Name: "triangle odd-girth ℓ = 3",
            Detail: "MinOddCycle({011,110,101})",
            Expected: "3",
            Actual: ell.ToString(CultureInfo.InvariantCulture)));
        int mStar = 2 * ell + 3; // m* = 2ℓ + deg, deg = 3 for a genuine cycle
        cases.Add(new BatteryCase(
            Name: "first asymmetric moment m* = 2ℓ + deg (deg = 3)",
            Detail: $"2·{ell} + 3",
            Expected: "9",
            Actual: mStar.ToString(CultureInfo.InvariantCulture)));
        return cases;
    }
}
