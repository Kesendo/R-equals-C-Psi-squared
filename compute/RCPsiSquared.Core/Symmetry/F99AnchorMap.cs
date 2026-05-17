using System.Text;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Aggregates <see cref="IF99AnchorBearing"/> claims into a single
/// tabulation of the F99 canonical-trig-angle anchor set
/// {0, 1/8, 1/4, 3/8, 1/2}. For each anchor lists the Direct claims (if any),
/// the Parent claims that feed inheritance into it, and the Covers claims that
/// include it in a multi-anchor theorem. Gaps (anchors with no Direct claim)
/// are surfaced explicitly via <see cref="GapAnchors"/>.
///
/// <para>This is the operational answer to Tom's question (2026-05-17 night):
/// "wenn man sich dann den Print der Vererbung anschaut, sieht man die
/// offenen Lücken?" The map walks all IF99AnchorBearing claims passed in,
/// tabulates by anchor value, and the rendered print shows immediately which
/// F99 anchors have dedicated typed claims and which are only covered by F99
/// itself (= the open vantage-lens lücken the framework has not yet built
/// dedicated tooling for).</para>
///
/// <para>Construction is explicit (claims passed in via constructor): there is
/// no global registry walk. This keeps the map auditable — the test/caller
/// names which claims it considers, so the gap list cannot silently widen
/// when an unrelated claim is added elsewhere.</para>
/// </summary>
public sealed class F99AnchorMap
{
    /// <summary>The canonical F99 anchor set, in ascending order of α.</summary>
    public static IReadOnlyList<double> CanonicalAnchors { get; } =
        new[] { 0.0, 1.0 / 8.0, 1.0 / 4.0, 3.0 / 8.0, 1.0 / 2.0 };

    /// <summary>Human-readable label per anchor (γ value + brief structural role).</summary>
    public static IReadOnlyDictionary<double, string> AnchorLabels { get; } =
        new Dictionary<double, string>
        {
            [0.0]       = "γ=1 (Mirror endpoint, |+...+⟩ X⊗N-eigenstate)",
            [1.0 / 8.0] = "γ=√3/2 (depth-3, non-uniform Dicke c²=2√3+3)",
            [1.0 / 4.0] = "γ=√2/2 (silver-Dicke c²=1+√2, Bloch equator)",
            [3.0 / 8.0] = "γ=1/2 (K-intermediate, uniform Dicke c=1)",
            [1.0 / 2.0] = "γ=0 (Generic, any Π²-odd: W, GHZ, Bell)",
        };

    private const double Tol = 1e-12;

    private readonly IReadOnlyList<(Claim Claim, IF99AnchorBearing Anchor)> _claims;

    public F99AnchorMap(params Claim[] claims)
    {
        if (claims is null) throw new ArgumentNullException(nameof(claims));
        var bearing = new List<(Claim, IF99AnchorBearing)>();
        foreach (var c in claims)
        {
            if (c is null)
                throw new ArgumentException("F99AnchorMap claims must not contain null.", nameof(claims));
            if (c is not IF99AnchorBearing ib)
                throw new ArgumentException(
                    $"Claim '{c.Name}' does not implement IF99AnchorBearing.", nameof(claims));
            bearing.Add((c, ib));
        }
        _claims = bearing;
    }

    /// <summary>All claims that participate in the map, in insertion order.</summary>
    public IReadOnlyList<Claim> Claims => _claims.Select(t => t.Claim).ToList();

    /// <summary>Claims that have <see cref="F99AnchorRole.Direct"/> coverage of
    /// the given anchor value (within tolerance 1e-12).</summary>
    public IReadOnlyList<Claim> DirectClaimsAt(double anchor) =>
        _claims.Where(t => t.Anchor.F99Role == F99AnchorRole.Direct
                           && t.Anchor.F99AnchorValues.Any(a => Math.Abs(a - anchor) < Tol))
               .Select(t => t.Claim)
               .ToList();

    /// <summary>Claims that have <see cref="F99AnchorRole.Covers"/> coverage of
    /// the given anchor value.</summary>
    public IReadOnlyList<Claim> CoversClaimsAt(double anchor) =>
        _claims.Where(t => t.Anchor.F99Role == F99AnchorRole.Covers
                           && t.Anchor.F99AnchorValues.Any(a => Math.Abs(a - anchor) < Tol))
               .Select(t => t.Claim)
               .ToList();

    /// <summary>All Parent claims passed in (Parent role is anchor-agnostic;
    /// each Parent feeds the inheritance graph at all anchors it relates to
    /// structurally, without picking out specific α values).</summary>
    public IReadOnlyList<Claim> ParentClaims =>
        _claims.Where(t => t.Anchor.F99Role == F99AnchorRole.Parent)
               .Select(t => t.Claim)
               .ToList();

    /// <summary>F99 anchors with no Direct claim. These are the open vantage-
    /// lens lücken — anchors only covered by F99 itself, with no dedicated
    /// theorem class.</summary>
    public IReadOnlyList<double> GapAnchors =>
        CanonicalAnchors.Where(a => DirectClaimsAt(a).Count == 0).ToList();

    /// <summary>F99 anchors with at least one Direct claim.</summary>
    public IReadOnlyList<double> CoveredAnchors =>
        CanonicalAnchors.Where(a => DirectClaimsAt(a).Count > 0).ToList();

    /// <summary>Render the table as a string suitable for stdout / test output.
    /// Five rows (one per canonical anchor), four columns
    /// (α / γ-label / Direct claim / Covers claim), with a footer listing the
    /// gaps and Parent claims.</summary>
    public string Render()
    {
        var sb = new StringBuilder();
        sb.AppendLine();
        sb.AppendLine("F99 Canonical-Trig-Angle Anchor Coverage Map");
        sb.AppendLine(new string('=', 88));
        sb.AppendLine();
        sb.AppendLine($"  {"α",-8} {"γ / role",-50} {"Direct claim",-30} {"Covered by"}");
        sb.AppendLine($"  {new string('-', 8)} {new string('-', 50)} " +
                      $"{new string('-', 30)} {new string('-', 12)}");

        foreach (var anchor in CanonicalAnchors)
        {
            string alpha = AlphaLabel(anchor);
            string label = AnchorLabels[anchor];
            var direct = DirectClaimsAt(anchor);
            var covers = CoversClaimsAt(anchor);
            string directStr = direct.Count == 0 ? "(GAP)" : string.Join(", ", direct.Select(c => ShortName(c)));
            string coversStr = covers.Count == 0 ? "—"     : string.Join(", ", covers.Select(c => ShortName(c)));
            sb.AppendLine($"  {alpha,-8} {label,-50} {directStr,-30} {coversStr}");
        }

        sb.AppendLine();
        sb.AppendLine("Coverage summary");
        sb.AppendLine(new string('-', 88));
        sb.AppendLine($"  Direct claims:   {CoveredAnchors.Count} of {CanonicalAnchors.Count} anchors  " +
                      $"({string.Join(", ", CoveredAnchors.Select(AlphaLabel))})");
        sb.AppendLine($"  Gap anchors:     {GapAnchors.Count} of {CanonicalAnchors.Count} anchors  " +
                      $"({string.Join(", ", GapAnchors.Select(AlphaLabel))})");
        sb.AppendLine($"  Parent claims:   {ParentClaims.Count} feeding the inheritance graph " +
                      $"({string.Join(", ", ParentClaims.Select(ShortName))})");
        sb.AppendLine($"  Covers claims:   {_claims.Count(t => t.Anchor.F99Role == F99AnchorRole.Covers)} multi-anchor theorem(s)");

        if (GapAnchors.Count > 0)
        {
            sb.AppendLine();
            sb.AppendLine("Open lücken (F99 anchors with no Direct claim of their own):");
            foreach (var gap in GapAnchors)
            {
                sb.AppendLine($"  α = {AlphaLabel(gap)}: {AnchorLabels[gap]}");
            }
            sb.AppendLine();
            sb.AppendLine("Each gap is a vantage-lens not yet built. F99 covers all five anchors");
            sb.AppendLine("as one theorem (α = sin²(θ)/2 at the five canonical trig angles), but");
            sb.AppendLine("dedicated Direct claims would let us anchor analysis tooling at each");
            sb.AppendLine("vantage individually (per the 'anchors as vantages' reading from");
            sb.AppendLine("2026-05-17 night: each anchor IS a sichtweise; choosing one selects what");
            sb.AppendLine("counts as foreground/stable vs background/wave-breaking).");
        }
        return sb.ToString();
    }

    private static string AlphaLabel(double anchor)
    {
        if (Math.Abs(anchor)             < Tol) return "0";
        if (Math.Abs(anchor - 1.0 / 8.0) < Tol) return "1/8";
        if (Math.Abs(anchor - 1.0 / 4.0) < Tol) return "1/4";
        if (Math.Abs(anchor - 3.0 / 8.0) < Tol) return "3/8";
        if (Math.Abs(anchor - 1.0 / 2.0) < Tol) return "1/2";
        return anchor.ToString("G6");
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
