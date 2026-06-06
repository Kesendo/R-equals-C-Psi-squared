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
/// (<see cref="PalindromeSoftCertifier"/>) as a single registry-queryable Claim. The certifier carries BOTH
/// soft mechanisms: the DIAGONAL chiral K (the structured basis-state colourings, plus the site-swap
/// reflection) AND the NON-DIAGONAL hidden-Q routing (one uniform per-site product Q palindromizes a sum of
/// 2-body bilinears sharing a Q-family). It certifies "soft" iff one strategy applies; it never claims hard
/// (NotCertified carries no claim). This Claim asserts ONLY the two settled facts:
///
/// <list type="number">
///   <item><b>Soundness (one-sided)</b>: every case in the soundness battery is both Certified by
///     <see cref="PalindromeSoftCertifier.Certify"/> AND not Hard by the spectral authority
///     <see cref="PauliPairTrichotomy.Classify"/>. So a certificate never lies: certified ⟹ not hard.
///     The battery spans the strategies: XY+YX (ExcitationPairing), XZ+ZX (ExcitationParity),
///     XY−YX (LinearSiteColoring, the chain chiral-K reading the §7 diagonal-K criterion scales), and
///     XX+XZ (Routing, the non-diagonal hidden-Q family {P4}).</item>
///   <item><b>The structural ceiling</b> (the receded incompleteness, PROOF_F103 §7.12): with routing
///     added, the non-bipartite-soft 2-body class (XX+XZ) is now CERTIFIED, no longer the ceiling. The
///     remaining ceiling is the k-body routed-soft frontier (Stufe B): XZX+XZY+YZX is soft
///     (<see cref="PauliPairTrichotomy.Classify"/> == Soft, the k-body overload, verified at N=4,5,6) yet
///     NotCertified, because the routing family table is 2-body and cannot reach a 3-body routed case.
///     NotCertified therefore does not imply not-soft.</item>
/// </list>
///
/// <para>The k-body routing mechanism (Stufe B, the soft-criterion for the routed-soft frontier beyond the
/// 2-body family table) is out of scope and is NOT asserted here.</para>
///
/// <para>Tier: Tier1Candidate. Routing is reused as a HELPER (like <see cref="PalindromeMaskClassifier"/>),
/// so no new typed parent and no tier change. Typed parents: <see cref="F87DiagonalCellBipartiteWitnessSet"/>
/// (the §7 diagonal-K bipartite criterion the certifier's linear strategy scales, Tier1Candidate) and
/// <see cref="F87TrichotomyClassification"/> (the spectral authority the soundness is checked against,
/// Tier1Derived). The strength-inheritance check is parent ≥ child, i.e. 4 ≥ 4 and 5 ≥ 4, both pass.</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md</c> §7.12 +
/// <see cref="PalindromeSoftCertifier"/> + <c>PalindromeSoftCertifierCeilingTests</c>.</para></summary>
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

    /// <summary>The §7.12 structural ceiling witness (XZX+XZY+YZX): a k-body routed-soft case the
    /// 2-body-gated certifier cannot reach, so it cannot (and must not) certify it. The pair
    /// <see cref="Holds"/> iff it is soft (by the spectral authority) and NotCertified.</summary>
    public readonly record struct CeilingWitness(string Name, bool IsSoft, bool Certified)
    {
        public bool Holds => IsSoft && !Certified;
    }

    /// <summary>Chain providing N for the spectral checks.</summary>
    public ChainSystem Chain { get; }

    /// <summary>Soundness battery: each case certified-and-not-hard (the one-sided soundness property).</summary>
    public IReadOnlyList<SoundnessCase> SoundnessBattery { get; }

    /// <summary>The XZX+XZY+YZX k-body routed-soft ceiling witness: soft, NotCertified.</summary>
    public CeilingWitness Ceiling { get; }

    public PalindromeSoftCertifierClaim(ChainSystem chain)
        : base("§7.12 soft certifier: sound, structurally incomplete (the k-body routed-soft ceiling)",
               Tier.Tier1Candidate,
               "docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md §7.12 + " +
               "compute/RCPsiSquared.Diagnostics/F87/PalindromeSoftCertifier.cs + " +
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
        SoundnessBattery.All(c => c.Passes) && Ceiling.Holds;

    public override string DisplayName =>
        $"§7.12 soft certifier (sound + k-body routed-soft ceiling, N={Chain.N}, {SoundnessBattery.Count} soundness cases)";

    public override string Summary =>
        $"sound one-sided certifier: {SoundnessPassCount}/{SoundnessBattery.Count} battery cases certified-and-not-hard; " +
        $"ceiling: {Ceiling.Name} soft, NotCertified ({(Ceiling.Holds ? "PASS" : "FAIL")}) ({Tier.Label()})";

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
            yield return new InspectableNode("ceiling (the k-body routed-soft frontier)",
                summary: $"{Ceiling.Name}: soft={Ceiling.IsSoft}, certified={Ceiling.Certified} " +
                         "(" + (Ceiling.Holds ? "PASS" : "FAIL") + ")");
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
    ///   <item><b>XX+XZ</b>: non-diagonal hidden-Q routing, the shared uniform family {P4}, certified by
    ///     Routing (a non-bipartite basis-state graph the colourings cannot reach).</item>
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

    /// <summary>The §7.12 ceiling witness XZX+XZY+YZX: a 3-body routed-soft case (Stufe B). It is soft by
    /// the spectral authority (the k-body <see cref="PauliPairTrichotomy.Classify(ChainSystem,
    /// IReadOnlyList{PauliTerm}, double, double, PauliLetter)"/> overload) yet NotCertified, because the
    /// routing family table is 2-body and cannot reach a 3-body routed case. Its soft verdict is
    /// established at N = 4, 5, 6 (NOT at N = 3, where k = 3 fills the whole chain, a different regime), so
    /// it is checked on the claim's chain when N ≥ 4, else on a fixed N = 4 chain.
    /// Mirrors <c>PalindromeSoftCertifierCeilingTests.KBodyRoutedSoft_IsRealAndBeyondTheRoutingTable</c>.</summary>
    private static CeilingWitness BuildCeiling(ChainSystem chain)
    {
        var terms = H("XZX", "XZY", "YZX");
        // The 3-body witness needs room: its soft verdict is established at N ≥ 4. At N = 3 (k fills the
        // whole chain) it is a different regime, so fall back to a fixed N = 4 chain there.
        var soundChain = chain.N >= 4 ? chain : new ChainSystem(N: 4, J: 1.0, GammaZero: 0.05);
        var spectral = PauliPairTrichotomy.Classify(soundChain, terms);
        var cert = PalindromeSoftCertifier.Certify(terms, soundChain.N);
        return new CeilingWitness(
            Name: "XZX+XZY+YZX",
            IsSoft: spectral == TrichotomyClass.Soft,
            Certified: cert.Certified);
    }
}
