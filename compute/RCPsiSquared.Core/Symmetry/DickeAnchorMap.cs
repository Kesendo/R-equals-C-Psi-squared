using System.Text;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Aggregates <see cref="IDickeAnchorBearing"/> claims into a single
/// tabulation of the <see cref="DickeAnchor"/> 3-anchor set
/// {Mirror, KIntermediate, Generic}. For each anchor lists the Direct claims
/// (if any), the Parent claims that feed inheritance into it, and the Covers
/// claims that include it in a multi-anchor theorem. Gaps (anchors with no
/// Direct claim) surface explicitly via <see cref="GapAnchors"/>.
///
/// <para>Mirror of <see cref="F99AnchorMap"/> for the Dicke 3-anchor subset.
/// The two maps cover overlapping but not identical territories:
/// DickeAnchor's three positions {Mirror, KIntermediate, Generic} correspond
/// to F99's first / fourth / fifth anchors {0, 3/8, 1/2}; the F99 second and
/// third anchors {1/8, 1/4} (depth-3 and silver-Dicke) are non-uniform Dicke
/// and live only in F99's lens.</para>
///
/// <para>This is the second I*Bearing + *Map pair built tonight (after
/// F99AnchorMap). The pattern is now established: tiny interface + map +
/// test confirms what we have and surfaces what's missing, for any closed
/// finite anchor set.</para>
/// </summary>
public sealed class DickeAnchorMap
{
    /// <summary>The canonical DickeAnchor set, in order of decreasing γ
    /// (Mirror γ=1 → KIntermediate γ=1/2 → Generic γ=0, equivalently
    /// increasing α: 0 → 3/8 → 1/2).</summary>
    public static IReadOnlyList<DickeAnchor> CanonicalAnchors { get; } =
        new[] { DickeAnchor.Mirror, DickeAnchor.KIntermediate, DickeAnchor.Generic };

    /// <summary>Human-readable label per anchor (γ + α + state description).</summary>
    public static IReadOnlyDictionary<DickeAnchor, string> AnchorLabels { get; } =
        new Dictionary<DickeAnchor, string>
        {
            [DickeAnchor.Mirror]        = "γ=1, α=0 (X⊗N eigenstate; c²=∞ uniform-Dicke limit)",
            [DickeAnchor.KIntermediate] = "γ=1/2, α=3/8 (uniform Dicke c=1, K-intermediate)",
            [DickeAnchor.Generic]       = "γ=0, α=1/2 (any Π²-odd state; c²=0 limit, W/GHZ/Bell)",
        };

    private readonly IReadOnlyList<(Claim Claim, IDickeAnchorBearing Bearing)> _claims;

    public DickeAnchorMap(params Claim[] claims)
    {
        if (claims is null) throw new ArgumentNullException(nameof(claims));
        var bearing = new List<(Claim, IDickeAnchorBearing)>();
        foreach (var c in claims)
        {
            if (c is null)
                throw new ArgumentException("DickeAnchorMap claims must not contain null.", nameof(claims));
            if (c is not IDickeAnchorBearing ib)
                throw new ArgumentException(
                    $"Claim '{c.Name}' does not implement IDickeAnchorBearing.", nameof(claims));
            bearing.Add((c, ib));
        }
        _claims = bearing;
    }

    /// <summary>All claims in insertion order.</summary>
    public IReadOnlyList<Claim> Claims => _claims.Select(t => t.Claim).ToList();

    /// <summary>Direct claims at the given DickeAnchor.</summary>
    public IReadOnlyList<Claim> DirectClaimsAt(DickeAnchor anchor) =>
        _claims.Where(t => t.Bearing.DickeRole == DickeAnchorRole.Direct
                           && t.Bearing.DickeAnchors.Contains(anchor))
               .Select(t => t.Claim)
               .ToList();

    /// <summary>Covers claims that include the given DickeAnchor.</summary>
    public IReadOnlyList<Claim> CoversClaimsAt(DickeAnchor anchor) =>
        _claims.Where(t => t.Bearing.DickeRole == DickeAnchorRole.Covers
                           && t.Bearing.DickeAnchors.Contains(anchor))
               .Select(t => t.Claim)
               .ToList();

    /// <summary>All Parent claims.</summary>
    public IReadOnlyList<Claim> ParentClaims =>
        _claims.Where(t => t.Bearing.DickeRole == DickeAnchorRole.Parent)
               .Select(t => t.Claim)
               .ToList();

    /// <summary>DickeAnchors with no Direct claim. These are the open
    /// vantage-lens lücken — anchors only covered by F99 (via Covers role)
    /// with no dedicated Dicke-specific theorem class.</summary>
    public IReadOnlyList<DickeAnchor> GapAnchors =>
        CanonicalAnchors.Where(a => DirectClaimsAt(a).Count == 0).ToList();

    /// <summary>DickeAnchors with at least one Direct claim.</summary>
    public IReadOnlyList<DickeAnchor> CoveredAnchors =>
        CanonicalAnchors.Where(a => DirectClaimsAt(a).Count > 0).ToList();

    /// <summary>Render the table as a string suitable for stdout / test
    /// output. Three rows (one per Dicke anchor) plus a footer listing
    /// gaps and Parent claims.</summary>
    public string Render()
    {
        var sb = new StringBuilder();
        sb.AppendLine();
        sb.AppendLine("DickeAnchor 3-Set Coverage Map");
        sb.AppendLine(new string('=', 88));
        sb.AppendLine();
        sb.AppendLine($"  {"DickeAnchor",-16} {"γ / α / state",-50} {"Direct claim",-30}");
        sb.AppendLine($"  {new string('-', 16)} {new string('-', 50)} {new string('-', 30)}");

        foreach (var anchor in CanonicalAnchors)
        {
            string name = anchor.ToString();
            string label = AnchorLabels[anchor];
            var direct = DirectClaimsAt(anchor);
            var covers = CoversClaimsAt(anchor);
            string directStr = direct.Count == 0 ? "(GAP)" : string.Join(", ", direct.Select(ShortName));
            string coversStr = covers.Count == 0 ? "—" : string.Join(", ", covers.Select(ShortName));
            sb.AppendLine($"  {name,-16} {label,-50} {directStr,-30}");
            sb.AppendLine($"  {"",-16} {"  covered by:",-50} {coversStr}");
        }

        sb.AppendLine();
        sb.AppendLine("Coverage summary");
        sb.AppendLine(new string('-', 88));
        sb.AppendLine($"  Direct claims:   {CoveredAnchors.Count} of {CanonicalAnchors.Count} anchors  " +
                      $"({string.Join(", ", CoveredAnchors)})");
        sb.AppendLine($"  Gap anchors:     {GapAnchors.Count} of {CanonicalAnchors.Count} anchors  " +
                      $"({string.Join(", ", GapAnchors)})");
        sb.AppendLine($"  Parent claims:   {ParentClaims.Count} feeding the inheritance graph " +
                      $"({string.Join(", ", ParentClaims.Select(ShortName))})");
        sb.AppendLine($"  Covers claims:   {_claims.Count(t => t.Bearing.DickeRole == DickeAnchorRole.Covers)} multi-anchor theorem(s)");

        if (GapAnchors.Count > 0)
        {
            sb.AppendLine();
            sb.AppendLine("Open lücken (DickeAnchors with no Direct claim of their own):");
            foreach (var gap in GapAnchors)
            {
                sb.AppendLine($"  {gap}: {AnchorLabels[gap]}");
            }
            sb.AppendLine();
            sb.AppendLine("Each gap is a Dicke vantage not yet dedicated as a typed theorem class.");
            sb.AppendLine("F99 covers them as a subset of its 5-anchor sin²(θ)/2 mechanism, but");
            sb.AppendLine("dedicated Dicke-Direct claims would let us anchor analysis at each");
            sb.AppendLine("anchor specifically — particularly Mirror (γ=1, X⊗N eigenstate facts");
            sb.AppendLine("as a standalone theorem) and Generic (γ=0, the W/GHZ/Bell family).");
        }
        return sb.ToString();
    }

    private static string ShortName(Claim c)
    {
        var t = c.GetType().Name;
        if (t.EndsWith("Pi2Inheritance"))
            t = t[..^"Pi2Inheritance".Length];
        if (t.EndsWith("Claim"))
            t = t[..^"Claim".Length];
        return t;
    }
}
