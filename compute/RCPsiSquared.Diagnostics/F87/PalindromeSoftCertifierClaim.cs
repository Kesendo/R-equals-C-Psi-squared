using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using Strategy = RCPsiSquared.Diagnostics.F87.PalindromeSoftCertifier.SoftStrategy;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>Schicht-1 surface for the §7.12 Liouvillian-free SOFT-certifier
/// (<see cref="PalindromeSoftCertifier"/>) as a single registry-queryable Claim. The certifier carries the
/// DIAGONAL chiral K (the structured basis-state colourings, plus the site-swap reflection) AND the
/// non-diagonal hidden-Q routing mechanisms: the 2-body family-intersection routing (Stufe A, one uniform
/// per-site product Q palindromizes a sum of 2-body bilinears sharing a Q-family), the derived k-body
/// per-term routing (Stufe B, the per-term k-site anticommutator {Q_k, [T,·]_k} = 0 reaches the routable
/// k-body cases the 2-body table misses), the single-site-field products, and the window-summed golden
/// routing (Stufe B′, F116: the period-4 golden router whose window lemma the per-term lens cannot see).
/// It certifies "soft" iff one strategy applies; it never claims hard (NotCertified carries no claim). The
/// wrapped certifier now also exposes the two-sided <see cref="PalindromeSoftCertifier.Decide"/> (soft via
/// <see cref="PalindromeSoftCertifier.Certify"/> plus the N-free hard verdict from the F115 diagonal-cell
/// valuation); this Claim asserts ONLY the two settled facts:
///
/// <list type="number">
///   <item><b>Soundness (one-sided)</b>: every case in the soundness battery is both Certified by
///     <see cref="PalindromeSoftCertifier.Certify"/> AND not Hard by the spectral authority
///     <see cref="PauliPairTrichotomy.Classify"/>. So a certificate never lies: certified ⟹ not hard.
///     The battery spans the strategies: XY+YX (ExcitationPairing), XZ+ZX (ExcitationParity),
///     XY−YX (LinearSiteColoring, the chain chiral-K reading the §7 diagonal-K criterion scales),
///     XX+XZ (Routing, the 2-body non-diagonal hidden-Q family {P4}), XIX+XXY+YXX (RoutingKBody, the
///     derived k-body per-term routing, routed by the P4 pattern), IXI+IIY+YII (SingleSiteField), and the
///     two Z-middle cases XZX+XZY+YZX, YZY+XZY+YZX (RoutingWindowSummed, the golden period-4 router and
///     its X↔Y conjugate).</item>
///   <item><b>The structural ceiling is CLOSED at zero</b> (PROOF_F103 §7.12 +
///     docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md, F116): the arc ran 6 → 4 → 2 → 0. The 2 Z-middle cases,
///     once read as "palindromized only by a non-local Π", are per-site routable after all: the golden
///     period-4 router palindromizes each under the WINDOW-SUMMED anticommutator condition (the
///     cancellation is cross-template inside one window, which the per-term Stufe B correctly cannot see),
///     and the RoutingWindowSummed strategy certifies them. The recomputed ceiling (the soft-yet-NotCertified
///     members of the k=3 sliding-window family) is EMPTY; the claim asserts Count == 0. The remaining
///     incompleteness of the certifier is coverage/scalability (soft cases whose routers sit outside the
///     scalable strategies), not locality. NotCertified still does not imply not-soft.</item>
/// </list>
///
/// <para>The earlier ceiling members are all LOCAL: XIX+XIY+YIX and YIY+XIY+YIX route via a
/// continuous-uniform per-site Q (continuous-sum, the 6 to 4 step), IXI+IIY+YII, IYI+IIX+XII route via a
/// site-varying per-site product of single-site crossover maps, certified by SingleSiteField (the 4 to 2
/// step), and the Z-middle pair routes via the golden window-summed router (the 2 to 0 step), see
/// experiments/CEILING_FOUR_NONLOCAL_CASES.md.</para>
///
/// <para>Tier: Tier1Candidate. The routing mechanisms are reused as HELPERS (like
/// <see cref="PalindromeMaskClassifier"/>; the k-body and window-summed legs via
/// <see cref="KBodyPalindromeRouting"/>), so no new typed parent and no tier change. Typed parents:
/// <see cref="F87DiagonalCellBipartiteWitnessSet"/> (the §7 diagonal-K bipartite criterion the certifier's
/// linear strategy scales, Tier1Derived) and <see cref="F87TrichotomyClassification"/> (the spectral
/// authority the soundness is checked against, Tier1Derived). The strength-inheritance check is
/// parent ≥ child, i.e. 5 ≥ 4 and 5 ≥ 4, both pass.</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md</c> §7.12 +
/// <c>docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md</c> (F116, the Stufe B′ leg) +
/// <see cref="PalindromeSoftCertifier"/> + <see cref="KBodyPalindromeRouting"/> (the Stufe B + B′ legs) +
/// <c>PalindromeSoftCertifierCeilingTests</c>.</para></summary>
public sealed class PalindromeSoftCertifierClaim : Claim
{
    /// <summary>One soundness-battery case: a named soft Hamiltonian with the certifier's verdict
    /// (Certified + which Strategy) and the spectral authority's verdict (NotHard). A case
    /// <see cref="Passes"/> iff the certifier certified it AND the spectrum confirms it is not hard:
    /// the one-sided soundness property "certificate ⟹ not hard".</summary>
    public readonly record struct SoundnessCase(
        string Name, string Detail, bool Certified, Strategy Strategy, bool NotHard)
    {
        public bool Passes => Certified && NotHard;
    }

    /// <summary>A §7.12 structural-ceiling member: a case that is soft (by the spectral authority) yet
    /// NotCertified by every scalable strategy. <see cref="Holds"/> iff the member genuinely sits on the
    /// ceiling. Since the window-summed golden routing landed (Stufe B′, F116) the recomputed ceiling is
    /// EMPTY: the former witnesses XZX+XZY+YZX, YZY+XZY+YZX are certified-positive battery cases now, and
    /// the claim asserts that no member remains.</summary>
    public readonly record struct CeilingWitness(string Name, bool IsSoft, bool Certified)
    {
        public bool Holds => IsSoft && !Certified;
    }

    /// <summary>Chain providing N for the spectral checks.</summary>
    public ChainSystem Chain { get; }

    /// <summary>Soundness battery: each case certified-and-not-hard (the one-sided soundness property).</summary>
    public IReadOnlyList<SoundnessCase> SoundnessBattery { get; }

    /// <summary>The recomputed §7.12 structural ceiling: the soft-yet-NotCertified members of the k=3
    /// sliding-window family, EMPTY since F116 (the golden window-summed router certifies the last two
    /// Z-middle cases, the 2 to 0 step). The claim asserts Count == 0; a non-empty list would be a
    /// regression of a certification strategy.</summary>
    public IReadOnlyList<CeilingWitness> Ceiling { get; }

    public PalindromeSoftCertifierClaim(ChainSystem chain)
        : base("§7.12 soft certifier: sound, the k=3 windowed structural ceiling closed at zero (F116)",
               Tier.Tier1Candidate,
               "docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md §7.12 + " +
               "docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md (F116) + " +
               "compute/RCPsiSquared.Diagnostics/F87/PalindromeSoftCertifier.cs + " +
               "compute/RCPsiSquared.Diagnostics/F87/KBodyPalindromeRouting.cs + " +
               "compute/RCPsiSquared.Diagnostics.Tests/F87/PalindromeSoftCertifierCeilingTests.cs")
    {
        Chain = chain ?? throw new ArgumentNullException(nameof(chain));
        SoundnessBattery = BuildSoundnessBattery(chain);
        Ceiling = BuildCeiling(chain);
    }

    /// <summary>How many soundness cases pass (certified-and-not-hard).</summary>
    public int SoundnessPassCount => SoundnessBattery.Count(c => c.Passes);

    /// <summary>The full §7.12 self-check: every soundness case certified-and-not-hard (including the two
    /// Z-middle golden cases via RoutingWindowSummed), and the recomputed structural ceiling EMPTY.</summary>
    public bool SelfCheckPasses =>
        SoundnessBattery.All(c => c.Passes) && Ceiling.Count == 0;

    public override string DisplayName =>
        $"§7.12 soft certifier (sound, ceiling closed at zero, N={Chain.N}, {SoundnessBattery.Count} soundness cases)";

    public override string Summary =>
        $"sound one-sided certifier (chiral K + routing Stufe A/B/B′ + single-site fields): " +
        $"{SoundnessPassCount}/{SoundnessBattery.Count} battery cases certified-and-not-hard; " +
        $"structural ceiling closed at zero (F116 golden router): {Ceiling.Count} remaining " +
        $"({(Ceiling.Count == 0 ? "PASS" : "FAIL")}) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("chain",
                summary: $"N={Chain.N}, {Chain.HType}, {Chain.Topology}");
            yield return new InspectableNode("soundness (one-sided)",
                summary: "Certify(terms, N).Certified AND PauliPairTrichotomy.Classify(chain, terms) != Hard");
            foreach (var c in SoundnessBattery)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; certified={c.Certified} (strategy {c.Strategy}), not-hard={c.NotHard}, " +
                             (c.Passes ? "PASS" : "FAIL"));
            yield return new InspectableNode("ceiling (closed at zero, F116)",
                summary: $"6 → 4 → 2 → 0: the golden window-summed router certifies the Z-middle pair; " +
                         $"{Ceiling.Count} soft-yet-NotCertified members remain " +
                         "(" + (Ceiling.Count == 0 ? "PASS" : "FAIL") + ")");
            foreach (var c in Ceiling)
                yield return new InspectableNode($"ceiling REGRESSION: {c.Name}",
                    summary: $"soft={c.IsSoft}, certified={c.Certified}: a certification strategy regressed " +
                             "(the ceiling must be empty)");
        }
    }

    private static PauliTerm T(string label, Complex coefficient) =>
        new(PauliLabel.Parse(label), coefficient);

    private static List<PauliTerm> H(params string[] labels) =>
        labels.Select(l => T(l, Complex.One)).ToList();

    /// <summary>The soundness battery, one case per certification strategy, each self-checked against
    /// the spectral authority. All are soft, so all must be certified AND not hard.
    /// <list type="bullet">
    ///   <item><b>XY+YX</b>: pure pairing, certified by ExcitationPairing (topology-independent).</item>
    ///   <item><b>XZ+ZX</b>: all-odd-flip, bit_b-homogeneous, certified by ExcitationParity.</item>
    ///   <item><b>XY−YX</b>: pure hopping, chain-bipartite flip-masks, certified by LinearSiteColoring
    ///     (the chain chiral-K reading the §7 diagonal-K criterion scales).</item>
    ///   <item><b>XX+XZ</b>: non-diagonal 2-body hidden-Q routing (Stufe A), the shared uniform family {P4},
    ///     certified by Routing (a non-bipartite basis-state graph the colourings cannot reach).</item>
    ///   <item><b>XIX+XXY+YXX</b>: derived k-body per-term hidden-Q routing (Stufe B), certified by
    ///     RoutingKBody (a routable 3-body case the 2-body family table misses, routed by the P4 pattern).</item>
    ///   <item><b>IXI+IIY+YII</b>: a sum of single-site transverse fields, certified by SingleSiteField (the
    ///     per-site crossover product palindromizes the chain; the 4 to 2 step, an I-heavy case now local).</item>
    ///   <item><b>XZX+XZY+YZX</b> and <b>YZY+XZY+YZX</b>: the former Z-middle ceiling pair, certified by
    ///     RoutingWindowSummed (Stufe B′, F116: the golden period-4 router and its X↔Y conjugate, the
    ///     window-summed anticommutator condition the per-term lens cannot see; the 2 to 0 step).</item>
    /// </list></summary>
    private static IReadOnlyList<SoundnessCase> BuildSoundnessBattery(ChainSystem chain)
    {
        var battery = new (string Name, string Detail, List<PauliTerm> Terms)[]
        {
            ("XY+YX (ExcitationPairing)", "pure pairing, σ± colouring ⌊n/2⌋ mod 2", H("XY", "YX")),
            ("XZ+ZX (ExcitationParity)", "all-odd-flip, bit_b-homogeneous, n mod 2 colouring", H("XZ", "ZX")),
            ("XY-YX (LinearSiteColoring)", "pure hopping, chain-bipartite flip-masks (chiral K)",
                new List<PauliTerm> { T("XY", Complex.One), T("YX", -Complex.One) }),
            ("XX+XZ (Routing)", "hidden-Q uniform family {P4}, non-bipartite basis-state graph",
                H("XX", "XZ")),
            ("XIX+XXY+YXX (RoutingKBody)", "derived per-term k-site hidden-Q routing",
                H("XIX", "XXY", "YXX")),
            ("IXI+IIY+YII (SingleSiteField)", "sum of single-site transverse fields, per-site crossover product",
                H("IXI", "IIY", "YII")),
            ("XZX+XZY+YZX (RoutingWindowSummed)", "the golden period-4 window-summed router (Stufe B′, F116)",
                H("XZX", "XZY", "YZX")),
            ("YZY+XZY+YZX (RoutingWindowSummed)", "the X↔Y golden sibling, window-summed router (Stufe B′, F116)",
                H("YZY", "XZY", "YZX")),
        };

        var cases = new List<SoundnessCase>(battery.Length);
        foreach (var (name, detail, terms) in battery)
        {
            var cert = PalindromeSoftCertifier.Certify(terms, chain.N);
            var spectral = PauliPairTrichotomy.Classify(chain, terms);
            cases.Add(new SoundnessCase(
                Name: name,
                Detail: detail,
                Certified: cert.Certified,
                Strategy: cert.Strategy,
                NotHard: spectral != TrichotomyClass.Hard));
        }
        return cases;
    }

    /// <summary>Recompute the §7.12 structural ceiling: of the last two candidate cases (the Z-middle
    /// XZX+XZY+YZX, YZY+XZY+YZX, the 6 → 4 → 2 arc's survivors), keep any that is STILL soft-yet-NotCertified.
    /// With the window-summed golden routing strategy (Stufe B′, F116, PROOF_CEILING_GOLDEN_ROUTER.md) both
    /// are certified, so the returned list is EMPTY: the ceiling closed at zero, and the two cases live in
    /// the soundness battery as certified-positive witnesses (strategy RoutingWindowSummed) instead. A
    /// non-empty result is a certification regression, surfaced by <see cref="SelfCheckPasses"/>. Soft is
    /// established at N≥4, so a fixed N=4 chain is used when the claim's chain has N&lt;4.</summary>
    private static IReadOnlyList<CeilingWitness> BuildCeiling(ChainSystem chain)
    {
        var soundChain = chain.N >= 4 ? chain : new ChainSystem(N: 4, J: 1.0, GammaZero: 0.05);
        var labels = new[]
        {
            ("XZX+XZY+YZX", new[] { "XZX", "XZY", "YZX" }),
            ("YZY+XZY+YZX", new[] { "YZY", "XZY", "YZX" }),
        };
        var ceiling = new List<CeilingWitness>();
        foreach (var (name, ls) in labels)
        {
            var terms = H(ls);
            var spectral = PauliPairTrichotomy.Classify(soundChain, terms);
            var cert = PalindromeSoftCertifier.Certify(terms, soundChain.N);
            if (spectral == TrichotomyClass.Soft && !cert.Certified)
                ceiling.Add(new CeilingWitness(name, IsSoft: true, Certified: false));
        }
        return ceiling;
    }
}
