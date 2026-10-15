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
/// t_ℓ ≠ 0 branch. The t_ℓ = 0 branch closed the same day via the Pascal-Gram positivity
/// theorem, carried by <see cref="WindowedConverseAllGammaClaim"/> (now Tier1Derived, no
/// residual). Anchor: <c>docs/proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md</c> +
/// <c>simulations/f87_windowed_monomial_converse.py</c> +
/// <c>simulations/f87_girth_dichotomy.py</c> +
/// <c>simulations/f87_pascal_gram_positivity.py</c>.</para></summary>
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
        : base("F87 windowed-converse threshold (two-reflection spine): all-odd (#A_L,#A_R,#Q) parity ⟹ bipartite soft + non-bipartite #A≥2ℓ; the deg-1 class in closed form (girth dichotomy). Tier1Derived; the all-γ closure on both branches is carried by WindowedConverseAllGammaClaim (Pascal-Gram positivity, no residual since 2026-06-10)",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md + " +
               "simulations/f87_windowed_monomial_converse.py")
    {
        Cases = BuildBattery();
    }

    public string TwoReflections =>
        "𝓕 = F⊗F, R = I⊗F (F = X^⊗N): 𝓕A𝓕 = −A, 𝓕Q𝓕 = +Q, R·A_L·R = +A_L, R·A_R·R = −A_R, R·Q·R = −Q " +
        "(complex H included, via FH^TF = (FHF)^T = −H^T). For odd power-sums this forces #A_L, #A_R, #Q all odd.";

    public string AnyLitChiral =>
        "F enters only through F H F = −H and F Z_l F = −Z_l. Any K with K H K⁻¹ = −H that anticommutes with every " +
        "Z_l (a lit chiral) gives the same sign table through ρ ↦ KρK⁻¹ and ρ ↦ ρK⁻¹, and with it the parity, the " +
        "threshold, the girth dichotomy and the Pascal-Gram positivity: an H with a lit chiral pairs its spectrum at " +
        "every γ > 0 or at none. The golden ceiling has K = a^⊗N and the linear colouring ⊗[I,I,Z,Z], its soft branch.";

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
        "sets m* = 2ℓ+deg (γ³ when it fires; higher when not, first γ⁵ witness IIXY+ZXZY). R-sign in ladder " +
        "form was resolved the same day by the Pascal-Gram positivity theorem (every class at m* a sum of " +
        "squares or exactly zero), carried by WindowedConverseAllGammaClaim: hard ∀γ>0 on BOTH branches.";

    public override string DisplayName =>
        "F87 windowed-converse threshold (two-reflection spine, Tier1Derived)";

    public override string Summary =>
        "two reflections force #A_L,#A_R,#Q all odd in every odd power-sum word; bipartite ⟹ all odd power-sums " +
        "≡ 0 (soft), non-bipartite ⟹ #A ≥ 2ℓ, deg-1 positivity closed-form; the all-γ closure is carried by " +
        $"WindowedConverseAllGammaClaim (Pascal-Gram positivity, no residual); {PassCount}/{Cases.Count} PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Two reflections (𝓕, R) ⟹ all-odd parity", summary: TwoReflections);
            yield return new InspectableNode("Any lit chiral serves (§7)", summary: AnyLitChiral);
            yield return new InspectableNode("Threshold #A ≥ 2ℓ + soft re-proof", summary: Threshold);
            yield return new InspectableNode("deg-1 positivity (closed form)", summary: Deg1Positivity);
            yield return new InspectableNode("The girth ladder + Pascal-Gram closure (2026-06-10)", summary: FollowsModulo);
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
        int mStar = 2 * ell + 3; // m* = 2ℓ + deg; deg = 3 here because the K3 triangle sits on the t_ℓ ≡ 0 branch of the girth dichotomy (pure cycles with t_ℓ ≠ 0 fire deg 1 at 2ℓ+1)
        cases.Add(new BatteryCase(
            Name: "first asymmetric moment m* = 2ℓ + deg (deg = 3)",
            Detail: $"2·{ell} + 3",
            Expected: "9",
            Actual: mStar.ToString(CultureInfo.InvariantCulture)));
        // Any lit chiral serves (PROOF §7): the instance N = 3, H = 3·XZX + 2·XZY + 2·YZX + b₀ + b₂ with
        // a = 2X + Y, b = X − 2Y. K = a^⊗3 negates H and every Z_l (a² = 5, so K⁻¹ = K/125), while F = X^⊗3
        // does not negate H, so the F-based converse does not reach this H.
        var lit = ReadLitChiralInstance();
        cases.Add(new BatteryCase(
            Name: "any lit chiral serves: K = (2X+Y)^⊗3 carries the sign table where F = X^⊗3 does not",
            Detail: "K·H·K = −125·H and K·Z_l·K = −125·Z_l for l = 0, 1, 2; F·H·F ≠ −H (Gaussian integers)",
            Expected: "K lit chiral, F not chiral",
            Actual: (lit.KLitChiral ? "K lit chiral" : "K FAILS") + ", " + (lit.FChiral ? "F chiral" : "F not chiral")));
        // The moments do not see a transpose or the sign of −i (H → −H, H → Hᵀ leave them fixed), so the two
        // helpers that build M are pinned directly.
        var y = GMat.Letter('Y');
        bool helpersRight = GMat.Equal(GMat.Transpose(y), GMat.Combine(y, -1, y, 0))
                         && GMat.Equal(GMat.TimesMinusI(GMat.Letter('I')), new GMat(new long[2, 2], new long[,] { { -1, 0 }, { 0, -1 } }));
        cases.Add(new BatteryCase(
            Name: "the coherence-space builders: Yᵀ = −Y and −i·I",
            Detail: "GMat.Transpose and GMat.TimesMinusI, which the moments below cannot detect",
            Expected: "right",
            Actual: helpersRight ? "right" : "WRONG"));
        cases.Add(new BatteryCase(
            Name: "the lit-chiral instance fires at m* = 7 with deg 1, positive and odd in γ",
            Detail: "odd p_m = Tr(M^m), m = 1, 3, 5, 7, at γ = 1 and γ = 2 (exact)",
            Expected: "0, 0, 0, 8064000 | 0, 0, 0, 16128000",
            Actual: string.Join(", ", lit.MomentsAtOne) + " | " + string.Join(", ", lit.MomentsAtTwo)));
        cases.Add(new BatteryCase(
            Name: "the girth moment gives it: p₇ = 7·C(6,3)·Σ_l t₃²·γ, t₃ = Tr(Z_l H³)",
            Detail: "t_j for j = 1, 2, 3 at l = 0, 1, 2; Σ t₃² and 140·Σ t₃²",
            Expected: "t₁ = t₂ = 0, t₃ = (0, −240, 0), 140·57600 = 8064000",
            Actual: lit.GirthMoments));
        return cases;
    }

    private readonly record struct LitChiralInstance(
        bool KLitChiral, bool FChiral, IReadOnlyList<string> MomentsAtOne, IReadOnlyList<string> MomentsAtTwo, string GirthMoments);

    private static LitChiralInstance ReadLitChiralInstance()
    {
        const int n = 3, d = 1 << n;
        var h = GMat.Zero(d);
        foreach (var (label, weight) in new[] { ("XZX", 3L), ("XZY", 2L), ("YZX", 2L), ("bII", 1L), ("IIb", 1L) })
            h = GMat.Combine(h, 1, GMat.Product(label), weight);
        var zs = Enumerable.Range(0, n)
            .Select(l => GMat.Product(string.Concat(Enumerable.Range(0, n).Select(s => s == l ? 'Z' : 'I'))))
            .ToList();

        var k = GMat.Product("aaa");
        long scale = 125;
        bool kChiral = GMat.Equal(GMat.Mul(GMat.Mul(k, h), k), GMat.Combine(h, -scale, h, 0))
                    && zs.All(z => GMat.Equal(GMat.Mul(GMat.Mul(k, z), k), GMat.Combine(z, -scale, z, 0)));
        var f = GMat.Product("XXX");
        bool fChiral = GMat.Equal(GMat.Mul(GMat.Mul(f, h), f), GMat.Combine(h, -1, h, 0));

        // M = −i(H ⊗ I − I ⊗ Hᵀ) + γ·Σ_l Z_l ⊗ Z_l on the row-stacked coherence space.
        var id = GMat.Identity(d);
        var a = GMat.TimesMinusI(GMat.Combine(GMat.Kron(h, id), 1, GMat.Kron(id, GMat.Transpose(h)), -1));
        var q = zs.Aggregate(GMat.Zero(d * d), (acc, z) => GMat.Combine(acc, 1, GMat.Kron(z, z), 1));
        IReadOnlyList<string> OddMoments(long gamma)
        {
            var m = GMat.Combine(a, 1, q, gamma);
            var power = GMat.Identity(d * d);
            var list = new List<string>();
            for (int p = 1; p <= 7; p++)
            {
                power = GMat.Mul(power, m);
                if (p % 2 == 0) continue;
                var (re, im) = GMat.Trace(power);
                list.Add(im == 0 ? re.ToString(CultureInfo.InvariantCulture) : $"{re}{(im < 0 ? "−" : "+")}{Math.Abs(im)}i");
            }
            return list;
        }

        var hPowers = new[] { h, GMat.Mul(h, h), GMat.Mul(GMat.Mul(h, h), h) };
        var t = zs.Select(z => hPowers.Select(hp => GMat.Trace(GMat.Mul(z, hp))).ToArray()).ToArray();
        bool lowVanish = t.All(row => row[0] == (0, 0) && row[1] == (0, 0));
        string t3 = string.Join(", ", t.Select(row => row[2].Im == 0
            ? row[2].Re.ToString(CultureInfo.InvariantCulture).Replace('-', '−') : "complex"));
        long sumSquares = t.Sum(row => row[2].Re * row[2].Re + row[2].Im * row[2].Im);
        string girth = (lowVanish ? "t₁ = t₂ = 0" : "t₁/t₂ NONZERO") + $", t₃ = ({t3}), 140·{sumSquares} = {140 * sumSquares}";

        return new LitChiralInstance(kChiral, fChiral, OddMoments(1), OddMoments(2), girth);
    }
}
