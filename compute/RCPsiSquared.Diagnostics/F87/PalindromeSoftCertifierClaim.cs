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
/// DIAGONAL chiral K (the structured basis-state colourings, plus the site-swap reflection) AND BOTH
/// non-diagonal hidden-Q routing mechanisms: the 2-body family-intersection routing (Stufe A, one uniform
/// per-site product Q palindromizes a sum of 2-body bilinears sharing a Q-family) AND the derived k-body
/// per-term routing (Stufe B, the per-term k-site anticommutator {Q_k, [T,·]_k} = 0 reaches the routable
/// k-body cases the 2-body table misses). It certifies "soft" iff one strategy applies; it never claims
/// hard (NotCertified carries no claim). This Claim asserts ONLY the two settled facts:
///
/// <list type="number">
///   <item><b>Soundness (one-sided)</b>: every case in the soundness battery is both Certified by
///     <see cref="PalindromeSoftCertifier.Certify"/> AND not Hard by the spectral authority
///     <see cref="PauliPairTrichotomy.Classify"/>. So a certificate never lies: certified ⟹ not hard.
///     The battery spans the strategies: XY+YX (ExcitationPairing), XZ+ZX (ExcitationParity),
///     XY−YX (LinearSiteColoring, the chain chiral-K reading the §7 diagonal-K criterion scales),
///     XX+XZ (Routing, the 2-body non-diagonal hidden-Q family {P4}), and XIX+XXY+YXX (RoutingKBody, the
///     derived k-body per-term routing, routed by the P4 pattern).</item>
///   <item><b>The structural ceiling</b> (the receded incompleteness, PROOF_F103 §7.12): with both routing
///     mechanisms added, the non-bipartite-soft 2-body class (XX+XZ, Stufe A) and the routable k-body cases
///     (Stufe B) are CERTIFIED, no longer the ceiling. The remaining ceiling is the NON-LOCAL k-body
///     routed-soft frontier, the 4 cases XZX+XZY+YZX, YZY+XZY+YZX, IXI+IIY+YII, IYI+IIX+XII: each is soft
///     (<see cref="PauliPairTrichotomy.Classify"/> == Soft, the k-body overload, verified at N=4,5,6) yet
///     NotCertified, because each admits NO per-site product Q at all (palindromized only by a non-local Π),
///     so the derived per-term routing declines it. NotCertified therefore does not imply not-soft.</item>
/// </list>
///
/// <para>Certifying the 4 non-local ceiling cases (XZX+XZY+YZX, YZY+XZY+YZX, IXI+IIY+YII, IYI+IIX+XII), which
/// admit no per-site product Q, is out of scope and is NOT asserted here. The 2 cases once counted with them,
/// XIX+XIY+YIX and YIY+XIY+YIX, are in fact LOCAL (a continuous-uniform per-site Q palindromizes them, verified;
/// NotCertified only because that Q routes via continuous-sum, outside the scalable strategies),
/// see experiments/CEILING_FOUR_NONLOCAL_CASES.md.</para>
///
/// <para>Tier: Tier1Candidate. Both routing mechanisms are reused as HELPERS (like
/// <see cref="PalindromeMaskClassifier"/>; the k-body leg via <see cref="KBodyPalindromeRouting"/>), so no
/// new typed parent and no tier change. Typed parents: <see cref="F87DiagonalCellBipartiteWitnessSet"/>
/// (the §7 diagonal-K bipartite criterion the certifier's linear strategy scales, Tier1Candidate) and
/// <see cref="F87TrichotomyClassification"/> (the spectral authority the soundness is checked against,
/// Tier1Derived). The strength-inheritance check is parent ≥ child, i.e. 4 ≥ 4 and 5 ≥ 4, both pass.</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md</c> §7.12 +
/// <see cref="PalindromeSoftCertifier"/> + <see cref="KBodyPalindromeRouting"/> (the Stufe B leg) +
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

    /// <summary>The §7.12 structural ceiling witness (XZX+XZY+YZX): a NON-LOCAL k-body routed-soft case
    /// that admits NO per-site product Q (palindromized only by a non-local Π), so the derived k-body
    /// per-term routing (Stufe B) declines it and the certifier cannot (and must not) certify it. The pair
    /// <see cref="Holds"/> iff it is soft (by the spectral authority) and NotCertified.</summary>
    public readonly record struct CeilingWitness(string Name, bool IsSoft, bool Certified)
    {
        public bool Holds => IsSoft && !Certified;
    }

    /// <summary>Chain providing N for the spectral checks.</summary>
    public ChainSystem Chain { get; }

    /// <summary>Soundness battery: each case certified-and-not-hard (the one-sided soundness property).</summary>
    public IReadOnlyList<SoundnessCase> SoundnessBattery { get; }

    /// <summary>The 4 k-body routed-soft non-local ceiling cases (XZX+XZY+YZX, YZY+XZY+YZX, IXI+IIY+YII,
    /// IYI+IIX+XII): each soft (spectral authority) yet NotCertified, admitting NO per-site product Q.</summary>
    public IReadOnlyList<CeilingWitness> Ceiling { get; }

    public PalindromeSoftCertifierClaim(ChainSystem chain)
        : base("§7.12 soft certifier: sound, structurally incomplete (the k-body routed-soft ceiling)",
               Tier.Tier1Candidate,
               "docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md §7.12 + " +
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

    /// <summary>The full §7.12 self-check: every soundness case certified-and-not-hard, and the
    /// XZX+XZY+YZX k-body ceiling pair (soft, NotCertified) holding.</summary>
    public bool SelfCheckPasses =>
        SoundnessBattery.All(c => c.Passes) && Ceiling.All(c => c.Holds);

    public override string DisplayName =>
        $"§7.12 soft certifier (sound + k-body routed-soft ceiling, N={Chain.N}, {SoundnessBattery.Count} soundness cases)";

    public override string Summary =>
        $"sound one-sided certifier (chiral K + 2-body routing Stufe A + derived k-body routing Stufe B): " +
        $"{SoundnessPassCount}/{SoundnessBattery.Count} battery cases certified-and-not-hard; " +
        $"non-local ceiling: {Ceiling.Count(c => c.Holds)}/{Ceiling.Count} cases soft+NotCertified " +
        $"({(Ceiling.All(c => c.Holds) ? "PASS" : "FAIL")}) ({Tier.Label()})";

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
            foreach (var c in Ceiling)
                yield return new InspectableNode($"ceiling: {c.Name}",
                    summary: $"non-local (no per-site product Q): soft={c.IsSoft}, certified={c.Certified} " +
                             "(" + (c.Holds ? "PASS" : "FAIL") + ")");
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

    /// <summary>The 4 §7.12 non-local ceiling witnesses. Each is a NON-LOCAL 3-body routed-soft case: soft
    /// by the spectral authority (the k-body overload, verified at N=4,5,6) yet NotCertified, admitting NO
    /// per-site product Q (palindromized only by a non-local Π). The 2 formerly-counted cases XIX+XIY+YIX,
    /// YIY+XIY+YIX are NOT here: they are LOCAL (continuous-uniform routable, verified N=3,4,5; NotCertified
    /// only because their router is an arbitrary continuous map outside the scalable strategies); see
    /// experiments/CEILING_FOUR_NONLOCAL_CASES.md. Soft is established at N≥4, so a fixed N=4 chain is used
    /// when the claim's chain has N<4.</summary>
    private static IReadOnlyList<CeilingWitness> BuildCeiling(ChainSystem chain)
    {
        var soundChain = chain.N >= 4 ? chain : new ChainSystem(N: 4, J: 1.0, GammaZero: 0.05);
        var labels = new[]
        {
            ("XZX+XZY+YZX", new[] { "XZX", "XZY", "YZX" }),
            ("YZY+XZY+YZX", new[] { "YZY", "XZY", "YZX" }),
            ("IXI+IIY+YII", new[] { "IXI", "IIY", "YII" }),
            ("IYI+IIX+XII", new[] { "IYI", "IIX", "XII" }),
        };
        var ceiling = new List<CeilingWitness>(labels.Length);
        foreach (var (name, ls) in labels)
        {
            var terms = H(ls);
            var spectral = PauliPairTrichotomy.Classify(soundChain, terms);
            var cert = PalindromeSoftCertifier.Certify(terms, soundChain.N);
            ceiling.Add(new CeilingWitness(name, IsSoft: spectral == TrichotomyClass.Soft, Certified: cert.Certified));
        }
        return ceiling;
    }
}
