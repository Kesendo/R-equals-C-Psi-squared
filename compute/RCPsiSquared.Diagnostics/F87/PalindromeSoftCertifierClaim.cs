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
/// (<see cref="PalindromeSoftCertifier"/>) as a single registry-queryable Claim. The certifier tries
/// three scalable structured 2-colourings of the basis-state graph (linear chiral K, excitation
/// pairing, excitation parity) and certifies "soft" iff one applies; it never claims hard
/// (NotCertified carries no claim). This Claim asserts ONLY the two settled facts:
///
/// <list type="number">
///   <item><b>Soundness (one-sided)</b>: every case in the soundness battery is both Certified by
///     <see cref="PalindromeSoftCertifier.Certify"/> AND not Hard by the spectral authority
///     <see cref="PauliPairTrichotomy.Classify"/>. So a certificate never lies: certified ⟹ not hard.
///     The battery spans all three strategies: XY+YX (ExcitationPairing), XZ+ZX (ExcitationParity),
///     and XY−YX (LinearSiteColoring, the chain chiral-K reading the §7 diagonal-K criterion scales).</item>
///   <item><b>The structural ceiling</b> (the proven incompleteness, PROOF_F103 §7.12): XX+XZ is soft
///     (<see cref="PauliPairTrichotomy.Classify"/> == Soft) yet its basis-state graph is NON-bipartite
///     (<see cref="BipartiteChirality"/> reports IsBipartite == false), so NO colouring at any degree can
///     reach it and the certifier returns NotCertified. NotCertified therefore does not imply not-soft.</item>
/// </list>
///
/// <para>The OPEN non-diagonal mechanism behind the non-bipartite-soft class (what soft-criterion lies
/// beyond basis-state bipartiteness) is out of scope and is NOT asserted here.</para>
///
/// <para>Tier: Tier1Candidate. Typed parents: <see cref="F87DiagonalCellBipartiteWitnessSet"/> (the §7
/// diagonal-K bipartite criterion the certifier's linear strategy scales, Tier1Candidate) and
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

    /// <summary>The §7.12 structural ceiling witness (XX+XZ): a soft pair whose basis-state graph is
    /// non-bipartite, so the certifier cannot (and must not) certify it. The triple
    /// <see cref="Holds"/> iff it is soft, non-bipartite, and NotCertified.</summary>
    public readonly record struct CeilingWitness(
        string Name, bool IsSoft, bool IsBipartite, bool Certified)
    {
        public bool Holds => IsSoft && !IsBipartite && !Certified;
    }

    /// <summary>Chain providing N for the spectral checks.</summary>
    public ChainSystem Chain { get; }

    /// <summary>Soundness battery: each case certified-and-not-hard (the one-sided soundness property).</summary>
    public IReadOnlyList<SoundnessCase> SoundnessBattery { get; }

    /// <summary>The XX+XZ structural ceiling witness: soft, non-bipartite, NotCertified.</summary>
    public CeilingWitness Ceiling { get; }

    public PalindromeSoftCertifierClaim(ChainSystem chain)
        : base("§7.12 soft certifier: sound, structurally incomplete (the non-bipartite-soft ceiling)",
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

    /// <summary>The full §7.12 self-check: every soundness case certified-and-not-hard, and the XX+XZ
    /// ceiling triple (soft, non-bipartite, NotCertified) holding.</summary>
    public bool SelfCheckPasses =>
        SoundnessBattery.All(c => c.Passes) && Ceiling.Holds;

    public override string DisplayName =>
        $"§7.12 soft certifier (sound + non-bipartite-soft ceiling, N={Chain.N}, {SoundnessBattery.Count} soundness cases)";

    public override string Summary =>
        $"sound one-sided certifier: {SoundnessPassCount}/{SoundnessBattery.Count} battery cases certified-and-not-hard; " +
        $"ceiling: XX+XZ soft, non-bipartite, NotCertified ({(Ceiling.Holds ? "PASS" : "FAIL")}) ({Tier.Label()})";

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
            yield return new InspectableNode("ceiling (the proven incompleteness)",
                summary: $"{Ceiling.Name}: soft={Ceiling.IsSoft}, bipartite={Ceiling.IsBipartite}, " +
                         $"certified={Ceiling.Certified} (" + (Ceiling.Holds ? "PASS" : "FAIL") + ")");
        }
    }

    private static PauliTerm T(string label, Complex coefficient) =>
        new(PauliLabel.Parse(label), coefficient);

    private static List<PauliTerm> H(params string[] labels) =>
        labels.Select(l => T(l, Complex.One)).ToList();

    /// <summary>The soundness battery, one case per certification strategy, each self-checked against
    /// the spectral authority. All three are soft, so all three must be certified AND not hard.
    /// <list type="bullet">
    ///   <item><b>XY+YX</b>: pure pairing, certified by ExcitationPairing (topology-independent).</item>
    ///   <item><b>XZ+ZX</b>: all-odd-flip, bit_b-homogeneous, certified by ExcitationParity.</item>
    ///   <item><b>XY−YX</b>: pure hopping, chain-bipartite flip-masks, certified by LinearSiteColoring
    ///     (the chain chiral-K reading the §7 diagonal-K criterion scales).</item>
    /// </list></summary>
    private static IReadOnlyList<SoundnessCase> BuildSoundnessBattery(ChainSystem chain)
    {
        var battery = new (string Name, string Detail, List<PauliTerm> Terms)[]
        {
            ("XY+YX (ExcitationPairing)", "pure pairing, σ± colouring ⌊n/2⌋ mod 2", H("XY", "YX")),
            ("XZ+ZX (ExcitationParity)", "all-odd-flip, bit_b-homogeneous, n mod 2 colouring", H("XZ", "ZX")),
            ("XY-YX (LinearSiteColoring)", "pure hopping, chain-bipartite flip-masks (chiral K)",
                new List<PauliTerm> { T("XY", Complex.One), T("YX", -Complex.One) }),
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

    /// <summary>The §7.12 ceiling witness XX+XZ: soft by the spectral authority, but with a
    /// non-bipartite basis-state graph (so no chiral K exists and the certifier returns NotCertified).
    /// Mirrors <c>PalindromeSoftCertifierCeilingTests.NonBipartiteSoft_IsRealAndBeyondAnyColouring</c>.</summary>
    private static CeilingWitness BuildCeiling(ChainSystem chain)
    {
        var terms = H("XX", "XZ");
        var bc = BipartiteChirality.Classify(chain, terms);
        var cert = PalindromeSoftCertifier.Certify(terms, chain.N);
        return new CeilingWitness(
            Name: "XX+XZ",
            IsSoft: bc.ActualClass == TrichotomyClass.Soft,
            IsBipartite: bc.IsBipartite,
            Certified: cert.Certified);
    }
}
