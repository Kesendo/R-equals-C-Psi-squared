using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>F116, the golden/metallic ceiling router (PROOF_CEILING_GOLDEN_ROUTER.md, reviewed bit-exact
/// 2026-06-22) as a standalone Tier1Derived Claim. This is the "C# witness first" inverse gap: the live
/// witness <see cref="Foundation.GoldenRouterWitness"/> (inspect --root router) and the reused helper
/// <see cref="KBodyPalindromeRouting"/> already carry the from-below truth; only the tiered Claim wrapper
/// was missing (open arc f116_golden_router_typed_claim).
///
/// <para><b>The result.</b> For the open-chain sliding-window Hamiltonian H = Σ_windows (XZX+XZY+YZX) (or
/// its X↔Y sibling YZY+XZY+YZX) at any N ≥ 3, under Z-dephasing at arbitrary site rates γ_l ≥ 0, there is
/// an invertible period-4 per-site product W = ⊗_l q_{l mod 4} with W L W⁻¹ = −L − 2σ, σ = Σ_l γ_l. The
/// per-site maps follow the [a, a, b, b] rhythm a = φX + Y, b = X − φY on the golden locus
/// α² − αβ − β² = 0, with h_l = q_l(Z) = (−1)^(l+1)·i·R(g_l) and q_l² = −(2+φ)·I (so cond(W) = 1). The
/// two-sided form W(ρ) = (2+φ)^(N/2)·P ρ Q reduces {W, [H,·]} = 0 to the chiral conditions
/// P H P⁻¹ = −H and Q H Q⁻¹ = −H. The palindromizer is a per-site PRODUCT: the Z-middle ceiling cases are
/// continuous-periodic-LOCAL, not non-local.</para>
///
/// <para><b>Mechanism (the window lemma).</b> The three-site, template-summed anticommutator
/// {q_w⊗q_{w+1}⊗q_{w+2}, [XZX+XZY+YZX,·]₃} = 0 EXACTLY (in the ring ℤ[φ]+iℤ[φ]) at every window offset,
/// while no single template's anticommutator vanishes: the cancellation is CROSS-TEMPLATE inside one
/// window. Additivity over windows gives {W, A} = 0 at every N ≥ 3. (This is why the per-term certifier
/// could not see the router; the per-term <see cref="KBodyPalindromeRouting.Routes"/> declines these
/// cases.)</para>
///
/// <para><b>Exclusion (deductive).</b> {W,A}=0 on the identity column forces G = ⊗_l g_l to commute with
/// H, reducing per window to two bilinears K₁, K₂. From them: no uniform router (g = 0 only), no period-2,
/// and the golden locus α²−αβ−β²=0 is the only gate; period-3 is impossible for N ≥ 5 (an alternating
/// assignment cannot close an odd cycle; the N=4 period-3 solution is the now-explained 2026-06-07
/// small-N artifact); the discrete Klein candidates P1/P4/M2/M are all off-locus.</para>
///
/// <para><b>Metallic family.</b> On the soft line t₂ = t₃, the router transports verbatim with a = (r, 1),
/// b = (1, −r), r(c) = (c+√(c²+4))/2 the metallic mean of c = t₁/t₂ (golden c=1 → φ, silver c=2 → 1+√2,
/// bronze c=3, c=0 the 45° frame), q_l² = −(1+r²)·I, and W_c L_c W_c⁻¹ = −L_c − 2σ for every real c. The
/// window lemma is a polynomial identity in r (degree ≤ 5, verified exactly zero at 8 rational nodes), so
/// it holds for all real (indeed complex) c by derivation.</para>
///
/// <para><b>Scope fences (the soft part, NOT Tier1Derived).</b> Rigidity (zero continuous physical moduli,
/// the router is unique modulo gauge) is PROVEN only at the golden and silver stations (Jacobian nullity =
/// gauge, N=5,6). The c=0 station's "8 physical moduli" is a finite-difference Jacobian count at N=5 ONLY
/// (spectral gap 3.5e9), NOT a closed-form/analytic result; 4 of 8 moduli are catalogued. Uniqueness is
/// scoped to INVERTIBLE W (the anticommutation equation alone carries abundant singular strata); within
/// the closed-form family at N=5 the router is essentially unique modulo per-site gauge, the four pattern
/// shifts, and an explicit order-32 sign group (128 gauge classes, every one a genuine palindromizer).</para>
///
/// <para>Tier1Derived: the existence + closed form + identity + exclusion + metallic-family carry no float
/// (exact ring / exact-Fraction interpolation); the end-to-end Lindbladian check is ~2e-15. Two typed
/// parents, both Tier1Derived (5 ≥ 5): <see cref="F1PalindromeIdentity"/> (the global palindrome
/// Π·L·Π⁻¹ = −L − 2Σγ·I that W realizes locally for the ceiling class) and
/// <see cref="WindowedConverseThresholdClaim"/> (the F87 two-reflection chiral spine 𝓕=F⊗F whose chiral
/// driving F H F = −H the two-sided form distributes into P ≠ Q). This router is the construction that
/// CLOSES <see cref="PalindromeSoftCertifierClaim"/>'s k=3 windowed soft-certifier ceiling at zero (the
/// 2 → 0 step); the certifier USES it as a helper, so the certifier is logically downstream, expressed
/// here as a see-cref, not a parent edge. Anchor: docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md +
/// ANALYTICAL_FORMULAS.md F116 + simulations/ceiling_golden_router.py + simulations/metallic_router_family.py
/// (both self-validating, exact ring arithmetic) + the live witness GoldenRouterWitness
/// (inspect --root router).</para></summary>
public sealed class GoldenRouterClaim : Claim
{
    /// <summary>The chain length used for the live window-summed routing re-check. The golden router is
    /// N-free; soft is established at N ≥ 4 for the Z-middle cases, so 4 suffices and is cheap.</summary>
    public const int CheckN = 4;

    /// <summary>Tolerance for the live-vs-closed-form metallic-mean MATCH (matches the witness/helper
    /// floor; the live ratio is a bisection root, the closed form the metallic quadratic).</summary>
    public const double MatchTolerance = 1e-9;

    /// <summary>The two typed parents, exposed for inspection.</summary>
    public F1PalindromeIdentity F1Palindrome { get; }
    public WindowedConverseThresholdClaim ChiralSpine { get; }

    /// <summary>One self-check case tying the Claim to the from-below router machinery
    /// (<see cref="KBodyPalindromeRouting"/>), in the style of <see cref="WindowedConverseThresholdClaim"/>.</summary>
    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    /// <summary>One metallic-family reading: the weight c, the live router-derived frame ratio (root of
    /// ‖{W,S}‖_F = 0), the closed-form metallic mean, and whether they agree within
    /// <see cref="MatchTolerance"/>.</summary>
    public readonly record struct MetallicReading(double C, double LiveRatio, double ClosedForm)
    {
        public double Difference => Math.Abs(LiveRatio - ClosedForm);
        public bool Matches => Difference < MatchTolerance;
    }

    // --- live readings, recomputed once at construct time from the reused router helper ---

    /// <summary>The window-summed router certificate for the golden case XZX+XZY+YZX, or null. Expected
    /// <see cref="KBodyPalindromeRouting.GoldenDescription"/>.</summary>
    public string? GoldenRoutes { get; }

    /// <summary>True iff the per-term lens DECLINES the golden case (the documented cross-template coverage
    /// gap, NOT non-locality): the window-summed router certifies what the per-term lens cannot.</summary>
    public bool GoldenPerTermDeclines { get; }

    /// <summary>The window-summed router certificate for the X↔Y sibling YZY+XZY+YZX, or null. Expected
    /// <see cref="KBodyPalindromeRouting.GoldenMirrorDescription"/>.</summary>
    public string? SiblingRoutes { get; }

    /// <summary>The metallic family evaluated live at c ∈ {0, 1, 2, 3}: live root vs closed-form mean.</summary>
    public IReadOnlyList<MetallicReading> MetallicReadings { get; }

    public IReadOnlyList<BatteryCase> Battery { get; }
    public int PassCount => Battery.Count(c => c.Passes);
    public bool SelfCheckPasses => PassCount == Battery.Count;

    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private static PauliTerm T(string label) => new(PauliLabel.Parse(label), Complex.One);
    private static List<PauliTerm> Terms(params string[] labels) => labels.Select(T).ToList();

    public GoldenRouterClaim(F1PalindromeIdentity f1Palindrome, WindowedConverseThresholdClaim chiralSpine)
        : base(
            "F116 golden/metallic router: an invertible period-4 per-site product W = ⊗_l q_{l mod 4} " +
            "([a,a,b,b], a=φX+Y, b=X−φY on the golden locus α²−αβ−β²=0, q_l²=−(2+φ)I) palindromizes the " +
            "Z-middle sliding-window ceiling H=Σ(XZX+XZY+YZX) (and its X↔Y sibling) under Z-dephasing: " +
            "W L W⁻¹ = −L − 2σ at every N≥3 for arbitrary site rates γ_l≥0; the Z-middle ceiling cases are " +
            "continuous-periodic-LOCAL, not non-local. Generalises to the metallic family r(c)=(c+√(c²+4))/2 " +
            "as a polynomial identity in r (golden c=1, silver c=2, bronze c=3). Exclusion (no uniform/no " +
            "period-2/period-3-impossible N≥5/Klein off-locus) is deductive from the K₁,K₂ identity-column " +
            "functional",
            Tier.Tier1Derived,
            "docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md + docs/ANALYTICAL_FORMULAS.md F116 + " +
            "simulations/ceiling_golden_router.py + simulations/metallic_router_family.py")
    {
        F1Palindrome = f1Palindrome ?? throw new ArgumentNullException(nameof(f1Palindrome));
        ChiralSpine = chiralSpine ?? throw new ArgumentNullException(nameof(chiralSpine));

        // Re-run the reused router machinery from below, once.
        var golden = Terms("XZX", "XZY", "YZX");
        var sibling = Terms("YZY", "XZY", "YZX");
        GoldenRoutes = KBodyPalindromeRouting.RoutesWindowSummed(golden, CheckN);
        GoldenPerTermDeclines = !KBodyPalindromeRouting.Routes(golden, CheckN);
        SiblingRoutes = KBodyPalindromeRouting.RoutesWindowSummed(sibling, CheckN);

        var readings = new List<MetallicReading>();
        foreach (double c in new[] { 0.0, 1.0, 2.0, 3.0 })
            readings.Add(new MetallicReading(c,
                KBodyPalindromeRouting.LiveMetallicRatio(c),
                KBodyPalindromeRouting.MetallicMean(c)));
        MetallicReadings = readings;

        Battery = BuildBattery();
    }

    private IReadOnlyList<BatteryCase> BuildBattery()
    {
        var cases = new List<BatteryCase>();

        cases.Add(new BatteryCase(
            Name: "golden case routes (window-summed)",
            Detail: "RoutesWindowSummed(XZX+XZY+YZX, N=4)",
            Expected: KBodyPalindromeRouting.GoldenDescription,
            Actual: GoldenRoutes ?? "null"));

        cases.Add(new BatteryCase(
            Name: "golden case declines per-term (cross-template gap)",
            Detail: "!Routes(XZX+XZY+YZX, N=4) (per-term lens cannot see cross-template cancellation)",
            Expected: bool.TrueString,
            Actual: GoldenPerTermDeclines.ToString()));

        cases.Add(new BatteryCase(
            Name: "X↔Y sibling routes (window-summed mirror)",
            Detail: "RoutesWindowSummed(YZY+XZY+YZX, N=4)",
            Expected: KBodyPalindromeRouting.GoldenMirrorDescription,
            Actual: SiblingRoutes ?? "null"));

        // The genuine metallic check: live router root = closed-form mean at the golden + silver stations.
        var golden = MetallicReadings.Single(r => r.C == 1.0);
        cases.Add(new BatteryCase(
            Name: "golden mean live = closed form (c=1 → φ)",
            Detail: $"|live {golden.LiveRatio.ToString("0.######", Inv)} − closed {golden.ClosedForm.ToString("0.######", Inv)}| < {MatchTolerance:0.0e+00}",
            Expected: "MATCH",
            Actual: golden.Matches ? "MATCH" : "MISMATCH"));

        var silver = MetallicReadings.Single(r => r.C == 2.0);
        cases.Add(new BatteryCase(
            Name: "silver mean live = closed form (c=2 → 1+√2)",
            Detail: $"|live {silver.LiveRatio.ToString("0.######", Inv)} − closed {silver.ClosedForm.ToString("0.######", Inv)}| < {MatchTolerance:0.0e+00}",
            Expected: "MATCH",
            Actual: silver.Matches ? "MATCH" : "MISMATCH"));

        // The residual ‖{W,S}‖_F vanishes at the golden mean and is nonzero off it (the dip the witness draws).
        double phi = KBodyPalindromeRouting.MetallicMean(1.0);
        double atMean = KBodyPalindromeRouting.WindowSummedAnticommutatorNorm(1.0, phi);
        double offMean = KBodyPalindromeRouting.WindowSummedAnticommutatorNorm(1.0, phi + 0.25);
        cases.Add(new BatteryCase(
            Name: "residual ‖{W,S}‖_F vanishes at r = φ, nonzero off it",
            Detail: $"‖{{W,S}}‖_F(c=1, r=φ) = {atMean.ToString("0.##e+00", Inv)} (≈0), at r=φ+0.25 = {offMean.ToString("0.###", Inv)} (>0)",
            Expected: "vanishes",
            Actual: (atMean < MatchTolerance && offMean > MatchTolerance) ? "vanishes" : "FAIL"));

        return cases;
    }

    public override string DisplayName =>
        "F116: W L W⁻¹ = −L − 2σ (golden period-4 local palindromizer for the Z-middle ceiling, Tier1Derived)";

    public override string Summary
    {
        get
        {
            int matched = MetallicReadings.Count(r => r.Matches);
            return "an invertible period-4 per-site product W ([a,a,b,b] on the golden locus) palindromizes " +
                   "the Z-middle sliding-window ceiling H=Σ(XZX+XZY+YZX): W L W⁻¹ = −L − 2σ at every N≥3, " +
                   "arbitrary site γ; LOCAL, not non-local. Live: golden→" + (GoldenRoutes ?? "null") +
                   ", sibling→" + (SiblingRoutes ?? "null") + ", per-term declines=" + GoldenPerTermDeclines +
                   $", metallic r(c) live=closed at {matched}/{MetallicReadings.Count} c; {PassCount}/{Battery.Count} PASS ({Tier.Label()})";
        }
    }

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the closed form (W, golden locus, two-sided form)",
                summary: "W = ⊗_l q_{l mod 4}, period-4 [a,a,b,b] with a = φX+Y, b = X−φY on the golden locus " +
                         "α²−αβ−β²=0 (slopes 1/φ and −φ, frame angle tan 2θ = 2), h_l = (−1)^(l+1)·i·R(g_l); " +
                         "q_l² = −(2+φ)·I so cond(W)=1. W L W⁻¹ = −L − 2σ, σ=Σ_l γ_l, at every N≥3 for arbitrary " +
                         "site-dependent γ_l≥0 (H = Σ_windows XZX+XZY+YZX and its X↔Y sibling). Two-sided form " +
                         "W(ρ) = (2+φ)^(N/2)·P ρ Q reduces {W,[H,·]}=0 to the chiral conditions P H P⁻¹ = −H, " +
                         "Q H Q⁻¹ = −H; G = PQ is a Hermitian product involution, an exact eigenmode at the floor " +
                         "L(G) = −2σ·G.");

            yield return new InspectableNode("the window lemma (cross-template, exact in the ring)",
                summary: "{q_w⊗q_{w+1}⊗q_{w+2}, [XZX+XZY+YZX,·]₃} = 0 exactly (ℤ[φ]+iℤ[φ]) at every window " +
                         "offset, while no single template's anticommutator vanishes: the cancellation is " +
                         "cross-template inside one window. Additivity over windows ⟹ {W,A}=0 at every N≥3 " +
                         "(single window included). This is why the per-term certifier cannot see the router; " +
                         "the per-term lens (KBodyPalindromeRouting.Routes) correctly returns false here.");

            yield return new InspectableNode("the exclusion theorems (deductive from K₁, K₂)",
                summary: "{W,A}=0 on the identity column I^⊗N forces G = ⊗_l g_l to commute with H, reducing per " +
                         "window to two bilinears K₁ = α_w α_{w+2}+α_w β_{w+2}+β_w α_{w+2}, K₂ = α_w β_{w+2}+β_w " +
                         "α_{w+2}−β_w β_{w+2}. From them, deductively: no uniform router (g=0 only), no period-2; " +
                         "the golden locus α²−αβ−β²=0 is the only gate (kernel map exchanges the two locus " +
                         "directions, so each parity chain alternates a↔b); period-3 impossible for N≥5 (an " +
                         "alternating assignment cannot close an odd cycle; the N=4 period-3 solution is the " +
                         "now-explained 2026-06-07 small-N artifact); the discrete Klein candidates P1/P4/M2/M " +
                         "are all off-locus (the one-line reason the discrete-periodic search came back empty).");

            yield return new InspectableNode("the metallic family (polynomial identity in r)",
                summary: "on the soft line t₂=t₃, the router transports verbatim with a=(r,1), b=(1,−r), " +
                         "r(c)=(c+√(c²+4))/2 the metallic mean of c=t₁/t₂ (golden c=1→φ, silver c=2→1+√2, " +
                         "bronze c=3→(3+√13)/2, c=0→r=1 the 45° frame), q_l²=−(1+r²)·I, W_c L_c W_c⁻¹ = −L_c − 2σ " +
                         "for every real c. The window lemma is a polynomial identity in r (degree ≤ 5, verified " +
                         "exactly zero at 8 rational nodes via exact Fraction arithmetic), so it holds for all " +
                         "real, indeed all complex, c by derivation. The c≠0 identity-column determinant factors " +
                         "as c·(α²−cαβ−β²): the metallic locus is the gate.");

            yield return new InspectableNode("SCOPE FENCE: the c=0 station is soft (NOT Tier1Derived)",
                summary: "Rigidity (zero continuous physical moduli; the router unique modulo gauge) is PROVEN " +
                         "only at the golden and silver stations: finite-difference Jacobian nullity = gauge " +
                         "dimension at N=5,6. The c=0 station's '8 physical moduli' (nullity 16 = 8 gauge + 8 " +
                         "physical, spectral gap 3.5e9) is a FINITE-DIFFERENCE JACOBIAN COUNT AT N=5 ONLY, NOT a " +
                         "closed-form/analytic result; 4 of 8 moduli are catalogued (the [v,v,v̄,v̄] X-axis-mirror " +
                         "continuum + the period-2 Pauli-axis routers, the only place period-2 exists in the " +
                         "family). The hard side off the line rests on the girth-ladder witness at (1,2,1) and " +
                         "the immediate-tilt scan, not on a closed-form p_{m*} for general weights. The existence " +
                         "side carries no float anywhere. This soft sub-result is held below the Tier1Derived line.");

            yield return new InspectableNode("uniqueness, scoped to invertible W",
                summary: "Uniqueness statements are scoped to INVERTIBLE W: the anticommutation equation alone " +
                         "also carries abundant singular strata. Within the closed-form family at N=5, probed " +
                         "exhaustively (48 unbiased searches, every found zero matched), the invertible solution " +
                         "set is 4 cyclic shifts × an explicit order-32 sign-decoration group (128 gauge classes), " +
                         "every one a genuine palindromizer; the sign group is the diagonal ±1 (anti)commutant of " +
                         "A. Honest phrasing: essentially unique, modulo per-site gauge, the four pattern shifts, " +
                         "and the explicit structural sign group.");

            // Live battery (from-below re-run of the reused router machinery).
            foreach (var c in Battery)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));

            // The relationship to the soft-certifier (a see-cref, deliberately NOT a parent edge).
            yield return new InspectableNode("closes the soft-certifier ceiling (downstream, not a parent)",
                summary: "This router is the construction that closes PalindromeSoftCertifierClaim's k=3 windowed " +
                         "soft-certifier ceiling at zero (the 2→0 step of the 6→4→2→0 arc): the two Z-middle " +
                         "cases move from the ceiling to certified-positive witnesses (RoutingWindowSummed). The " +
                         "certifier USES this router as a helper, so it is logically DOWNSTREAM of this claim; it " +
                         "is referenced here, not wired as a parent (and the certifier is Tier1Candidate < " +
                         "Tier1Derived, so parenting on it would also violate strength inheritance). The live lab " +
                         "is GoldenRouterWitness (inspect --root router).");

            yield return new InspectableNode("parent: F1PalindromeIdentity",
                summary: "the global palindrome Π·L·Π⁻¹ = −L − 2Σγ·I; the golden router W realizes that same " +
                         "palindrome form LOCALLY (as a per-site product) for the Z-middle ceiling Hamiltonian " +
                         "class, where the standard Π argument does not directly give a per-site product.");

            yield return new InspectableNode("parent: WindowedConverseThresholdClaim",
                summary: "the F87 two-reflection chiral spine (𝓕=F⊗F, R=I⊗F) and its chiral driving lemma " +
                         "F H F = −H (PROOF_F87_WINDOWED_MONOMIAL_CONVERSE §2); the router's two-sided form " +
                         "distributes that chiral action into DIFFERENT left and right factors P ≠ Q, the move " +
                         "that routes the Z-middle case where the single-reflection F H F = −H fails on the " +
                         "cross-terms.");
        }
    }
}
