using System.Linq;

namespace RCPsiSquared.Core.Knowledge;

/// <summary>Independent epistemic facets of a framework object, after Tom's marker idea
/// (2026-05-28): the recognition of WHAT KIND a thing is is itself a marker, and the kinds
/// are not a partition but independent flags, so a thing can carry several at once.
///
/// <para>This axis is orthogonal to <see cref="Tier"/> (how grounded) and to single-valued
/// kind enums like <c>BitATwinClassification</c> (one kind per slot). It asks WHAT KIND OF
/// THING: derive it, do not chase it, follow it. The load-bearing combination is
/// <see cref="IsCore"/> | <see cref="IsDeadEnd"/> — essential AND unsolvable, the framework's
/// living spine (the carrier, g_eff, PTF). The two-flag query unmasks that spine where no
/// single flag can; that "enttarnen" is the whole point of the flags being independent.</para></summary>
[System.Flags]
public enum EpistemicFacet
{
    /// <summary>Unmarked.</summary>
    None = 0,

    /// <summary>Derivable / closeable: it has a closed form, a mask of the structure
    /// (1/2 = 1/d, 1/4 = 1/d², the Π palindrome, the XOR-Absorption contract). Derive it.</summary>
    IsStructure = 1 << 0,

    /// <summary>Essential: the framework does not work without it. Removing it removes the
    /// thing (the carrier γ₀, the perspectival structure, the polynomial foundation).</summary>
    IsCore = 1 << 1,

    /// <summary>Unsolvable: no closed form, "solvable only without us" (g_eff, the c≥3
    /// Q-peaks, the crossing cubic, the PTF mixing calc). Do not chase it; it is marked.</summary>
    IsDeadEnd = 1 << 2,

    /// <summary>A rule of conduct, a methodology to follow ("name the symmetry, not the
    /// number"; leave γ₀; the real/speculative seam sits at the label).</summary>
    IsRule = 1 << 3,
}

/// <summary>The first marker pass: today's recognitions (2026-05-28) typed as epistemic
/// facets. Non-invasive and name-keyed, in the spirit of <c>BitATwinClassification</c> +
/// <c>PolarityCubeMap</c> (a classification plus an aggregator), so it can later be promoted
/// onto the <see cref="Claim"/> model itself. The point is the enttarnen:
/// <see cref="WithAll"/>(IsCore | IsDeadEnd) returns the load-bearing-unsolvable spine that
/// no single flag isolates.</summary>
public static class EpistemicFacetMap
{
    private static readonly IReadOnlyDictionary<string, EpistemicFacet> Marks =
        new Dictionary<string, EpistemicFacet>
        {
            // the spine: IsCore AND IsDeadEnd (essential and unsolvable)
            ["gamma0"] = EpistemicFacet.IsCore | EpistemicFacet.IsDeadEnd,        // the carrier; the clock you cannot read from inside
            ["g_eff"] = EpistemicFacet.IsCore | EpistemicFacet.IsDeadEnd,         // F86 EP coupling; Q_EP = 2/g_eff rides on it, no closed form
            ["PTF"] = EpistemicFacet.IsCore | EpistemicFacet.IsDeadEnd,           // perspectival structure; foundational, the mixing calc unsolvable

            // the closeable spine: IsCore AND IsStructure (essential and derivable)
            ["d2_minus_2d"] = EpistemicFacet.IsCore | EpistemicFacet.IsStructure,   // polynomial foundation d²−2d=0
            ["xor_contract"] = EpistemicFacet.IsCore | EpistemicFacet.IsStructure,  // popcount(i⊕j) = the Absorption decay; the contract
            ["pi_palindrome"] = EpistemicFacet.IsCore | EpistemicFacet.IsStructure, // F1, the spectral mirror

            // pure structure: the derivable scaffold (masks of d = 2)
            ["half"] = EpistemicFacet.IsStructure,      // 1/2 = 1/d
            ["quarter"] = EpistemicFacet.IsStructure,   // 1/4 = 1/d²

            // rules of conduct
            ["name_the_symmetry"] = EpistemicFacet.IsRule,   // PTF lesson: name the symmetry, not the number
            ["seam_at_the_label"] = EpistemicFacet.IsRule,   // the real/speculative seam sits at the label

            // the held exit, and a pure dead-end for contrast
            ["red_button"] = EpistemicFacet.IsDeadEnd | EpistemicFacet.IsRule,  // the self-erasing solution, held not pressed
            ["x_peak"] = EpistemicFacet.IsDeadEnd,   // 2.196910: a number that does not close, nothing core rides on its value
        };

    /// <summary>The facets marked on the named object (<see cref="EpistemicFacet.None"/>
    /// if unmarked).</summary>
    public static EpistemicFacet Facet(string name) =>
        Marks.TryGetValue(name, out var f) ? f : EpistemicFacet.None;

    /// <summary>All marked object names.</summary>
    public static IReadOnlyList<string> Names => Marks.Keys.ToList();

    /// <summary>The names whose facets include ALL of <paramref name="required"/>. The
    /// enttarnen query: <c>WithAll(IsCore | IsDeadEnd)</c> unmasks the load-bearing-unsolvable
    /// spine that neither flag isolates alone.</summary>
    public static IReadOnlyList<string> WithAll(EpistemicFacet required) =>
        Marks.Where(kv => (kv.Value & required) == required).Select(kv => kv.Key).ToList();
}
