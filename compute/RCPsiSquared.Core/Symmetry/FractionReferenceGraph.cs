using System.Text;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The static catalog of fraction-to-fraction references in the
/// framework's algebra (α-axis coordinates), plus query methods for walking
/// the graph.
///
/// <para>The "Level" structure isn't a stack — it's the WEB OF REFERENCES
/// between fractions. Each anker fraction in {0, 1/8, 1/4, 3/8, 1/2} is a
/// node on the α-axis; each operation (squaring, F86b α(γ), F98 long-time
/// asymptote, Π²-parity complement, ...) is an edge.</para>
///
/// <para><b>Important: this graph lives on the α-axis [0, 1/2], which is
/// ALREADY the folded picture under F86b α(γ) = (1−γ²)/2.</b> Two γ-values
/// ±|γ| map to the same α, so the α-axis cannot itself express polarity
/// mirror partners. 0 is therefore NOT a back-reference root — it is the
/// convergence point where the ±γ-polarity sides meet under folding.
/// The proper polarity-mirror structure lives on the γ-axis in a
/// separate map (see <c>PolarityMirrorMap</c>), and the existing typed
/// <see cref="PolarityLayerOriginClaim"/> already documents 0 as the active
/// substrate axis, not a passive midpoint or tree root.</para>
///
/// <para><b>Operationally</b>: the basis Q = J/γ₀ then evaluates these
/// references at specific Q-anchors. At Q=1 (Balance), F99 anker fractions
/// land directly in the Liouvillian Im=0 spectrum at small N (verified in
/// <c>simulations/wave_breaking_q_anchor_scan.py</c>). The static reference
/// graph IS algebraic; the Q-basis chooses WHICH references become
/// dynamically visible.</para>
///
/// <para><b>Edges curated here</b> (not exhaustive): the references encoded
/// in existing typed claims. New edges are added as new claims surface new
/// operations between fractions.</para>
/// </summary>
public sealed class FractionReferenceGraph
{
    /// <summary>The F99 canonical anker set (rational sin²(θ)/2 hits).</summary>
    public static IReadOnlyList<double> F99Ankers { get; } =
        new[] { 0.0, 1.0 / 8.0, 1.0 / 4.0, 3.0 / 8.0, 1.0 / 2.0 };

    /// <summary>The 3 wave-breaking-scan Q-anker basis values {1.0, 1.5, 2.0} — the
    /// subset of the canonical 9 Q-anchors that have explicit dynamical-visibility
    /// documented in <c>simulations/wave_breaking_q_anchor_scan.py</c>. The full
    /// 9-anchor structure (including onset edges 0.2/0.35, peak start 1.2, c=3 Q_peak
    /// 1.6, peak end / c=4,5 Q_peak 1.8, Endpoint orbit candidate 2.5) is typed in
    /// <see cref="QAnchorMap"/>; this property derives from
    /// <see cref="QAnchorMap.CanonicalAnchors"/> by filtering for the wave-breaking-scan
    /// source.</summary>
    public static IReadOnlyList<double> QBasisAnkers { get; } =
        QAnchorMap.CanonicalAnchors
            .Where(a => a.DocumentingSource.Contains("wave_breaking_q_anchor_scan"))
            .Select(a => a.Q)
            .ToList();

    /// <summary>All 9 named Q-anchors per <see cref="QAnchorMap"/>: onset edges 0.2 /
    /// 0.35, Balance 1.0, peak band {1.2, 1.5, 1.6, 1.8} (incl. c-specific Q_peak),
    /// Q_EP idealized 2.0, Endpoint orbit candidate 2.5.</summary>
    public static IReadOnlyList<double> NamedQAnchors { get; } =
        QAnchorMap.CanonicalAnchors.Select(a => a.Q).ToList();

    /// <summary>The standard catalog of fraction references, hard-coded from
    /// the typed claims that establish each operation. Adding a new claim
    /// that establishes a new fraction-to-fraction operation = add a new
    /// edge here. Order doesn't matter; queries don't depend on it.</summary>
    public static IReadOnlyList<FractionReference> StandardReferences { get; } =
        new FractionReference[]
        {
            // === 1/2 ↔ 1/4 (the most-connected pair, multi-viewpoint) ===
            new(1.0 / 2.0, 1.0 / 4.0, "squaring: (1/2)² = 1/4",
                FractionReferenceDirection.Backward, "QuarterAsBilinearMaxvalClaim"),
            new(1.0 / 2.0, 1.0 / 4.0, "argmax/maxval pair: p(1−p) peaks at p=1/2 with maxval 1/4",
                FractionReferenceDirection.Backward, "BilinearApexClaim + QuarterAsBilinearMaxvalClaim"),
            new(1.0 / 2.0, 1.0 / 4.0, "cardioid invariants: |z*|=1/2 and |z*|²=1/4 both invariant around F97 cardioid",
                FractionReferenceDirection.Backward, "F97CardioidHalfFixedPointPi2Inheritance"),
            new(1.0 / 2.0, 1.0 / 4.0, "dyadic ladder: a_2 = 1/2, a_3 = 1/4 (adjacent rungs)",
                FractionReferenceDirection.Backward, "Pi2DyadicLadderClaim"),

            // === 3/8 → 1/4 (F98 long-time asymptote) ===
            new(3.0 / 8.0, 1.0 / 4.0, "F98 long-time: (N+2)/[4(N+1)] from 3/8 at N=2 to 1/4 at N→∞",
                FractionReferenceDirection.Backward, "KIntermediateAsymptoteQuarterInheritance"),

            // === F86b α(γ) endpoints ===
            new(0.0, 1.0 / 2.0, "F86b α(γ) = (1−γ²)/2 from γ=1 (α=0) to γ=0 (α=1/2)",
                FractionReferenceDirection.Forward, "CanonicalTrigAnchorPi2Inheritance"),
            new(1.0 / 2.0, 0.0, "F86b α(γ) backward: γ=0 to γ=1 traverses non-anker continuum down to 0",
                FractionReferenceDirection.Backward, "CanonicalTrigAnchorPi2Inheritance"),

            // === α=0 as polarity-mirror convergence (NOT root) ===
            new(0.0, 0.0, "Mirror convergence: X⊗N|ψ⟩ = ±|ψ⟩ at γ=+1 OR γ=-1 both fold to α=0 via (1−γ²)/2; the ±γ-polarity sides meet here, NOT a tree root (see PolarityLayerOriginClaim, PolarityMirrorMap on γ-axis)",
                FractionReferenceDirection.Polarity, "XGlobalEigenstateMirrorPi2Inheritance + PolarityLayerOriginClaim"),

            // === F99 canonical-trig adjacency (anker chain via cos at canonical angles) ===
            new(1.0 / 2.0, 3.0 / 8.0, "F99 canonical-trig adjacent: γ=cos(60°)=1/2 ↔ γ=cos(90°)=0, α 3/8 ↔ 1/2",
                FractionReferenceDirection.Backward, "CanonicalTrigAnchorPi2Inheritance"),
            new(3.0 / 8.0, 1.0 / 4.0, "F99 canonical-trig adjacent: γ=cos(45°)=√2/2 → γ=cos(60°)=1/2",
                FractionReferenceDirection.Backward, "CanonicalTrigAnchorPi2Inheritance"),
            new(1.0 / 4.0, 1.0 / 8.0, "F99 canonical-trig adjacent: γ=cos(30°)=√3/2 → γ=cos(45°)=√2/2",
                FractionReferenceDirection.Backward, "CanonicalTrigAnchorPi2Inheritance"),
            new(1.0 / 8.0, 0.0, "F99 canonical-trig adjacent: γ=cos(0°)=1 → γ=cos(30°)=√3/2",
                FractionReferenceDirection.Backward, "CanonicalTrigAnchorPi2Inheritance"),

            // === Π²-parity complements (periodic-table fractions n/8 ↔ (8−n)/8) ===
            // The full n/8 ↔ (8−n)/8 complement family at valence n: three non-trivial
            // pairs + one self-mirror at n=4 (carbon). Pair-sum-constant is the F1
            // periodic palindrome over periods 2-6 (project_periodic_palindrome).
            new(1.0 / 8.0, 7.0 / 8.0, "Π²-parity complement: Li/Na (1/8) ↔ F/Cl (7/8) periodic-table",
                FractionReferenceDirection.Mirror, "project_v_effect_combinatorial memory"),
            new(1.0 / 4.0, 3.0 / 4.0, "Π²-parity complement: Be/Mg (1/4) ↔ O/S (3/4) periodic-table (alkaline-earth ↔ chalcogen, F1 palindrome)",
                FractionReferenceDirection.Mirror, "project_v_effect_combinatorial memory + project_periodic_palindrome"),
            new(3.0 / 8.0, 5.0 / 8.0, "Π²-parity complement: B/Al (3/8) ↔ N/P (5/8) periodic-table",
                FractionReferenceDirection.Mirror, "project_v_effect_combinatorial memory"),
            new(1.0 / 2.0, 1.0 / 2.0, "Π²-parity self-mirror at n=4: carbon (C/Si/Ge) is its own complement (4/8 ↔ 4/8); structural parallel to PolarityMirrorMap γ=0 self-mirror",
                FractionReferenceDirection.Mirror, "project_v_effect_combinatorial memory + project_periodic_palindrome"),

            // === Q-basis observation edges (Q=1 Balance hits F99 ankers in Liouvillian Im=0) ===
            // These are observation/dynamic edges, not algebraic. They show what the basis Q=J/γ₀
            // reveals when set to specific Q-anchors. Verified in wave_breaking_q_anchor_scan.py.
            new(1.0 / 2.0, 3.0 / 8.0, "Q=1 basis: N=2 chain Liouvillian Im=0 = 3/8 (lands on F99 KIntermediate anker)",
                FractionReferenceDirection.Backward, "wave_breaking_q_anchor_scan.py"),
            new(1.0 / 2.0, 1.0 / 2.0, "Q=1 basis: N=3 chain Liouvillian Im=0 = 1/2 (lands on F99 Generic anker)",
                FractionReferenceDirection.Backward, "wave_breaking_q_anchor_scan.py"),
            new(1.0 / 2.0, 1.0 / 4.0, "Q=1 basis: N=5 chain Liouvillian Im=0 = 1/4 (lands on F99 silver anker)",
                FractionReferenceDirection.Backward, "wave_breaking_q_anchor_scan.py"),
        };

    public IReadOnlyList<FractionReference> References { get; }

    public FractionReferenceGraph() : this(StandardReferences) { }

    public FractionReferenceGraph(IReadOnlyList<FractionReference> references)
    {
        References = references ?? throw new ArgumentNullException(nameof(references));
    }

    /// <summary>All references whose source is the given fraction (within
    /// tolerance 1e-12).</summary>
    public IReadOnlyList<FractionReference> ReferencesFrom(double fromFraction) =>
        References.Where(r => Math.Abs(r.FromFraction - fromFraction) < AnchorConstants.Tol).ToList();

    /// <summary>All references whose target is the given fraction.</summary>
    public IReadOnlyList<FractionReference> ReferencesTo(double toFraction) =>
        References.Where(r => Math.Abs(r.ToFraction - toFraction) < AnchorConstants.Tol).ToList();

    /// <summary>References from the given fraction in a specific direction.</summary>
    public IReadOnlyList<FractionReference> BackwardFrom(double fromFraction) =>
        ReferencesFrom(fromFraction)
            .Where(r => r.Direction == FractionReferenceDirection.Backward)
            .ToList();

    /// <summary>Walk all backward references from <paramref name="start"/>
    /// recursively, collecting reached fractions. Returns the set of
    /// fractions reachable by following backward arrows on the α-axis.</summary>
    public IReadOnlySet<double> BackwardClosure(double start)
    {
        var visited = new HashSet<double>();
        var queue = new Queue<double>();
        queue.Enqueue(start);
        visited.Add(start);
        while (queue.Count > 0)
        {
            double current = queue.Dequeue();
            foreach (var r in BackwardFrom(current))
            {
                // Skip self-loops at α=0 (the Polarity convergence edge is not
                // a backward step; it documents that ±γ both fold to α=0).
                if (Math.Abs(r.ToFraction - r.FromFraction) < AnchorConstants.Tol) continue;
                if (!visited.Any(v => Math.Abs(v - r.ToFraction) < AnchorConstants.Tol))
                {
                    visited.Add(r.ToFraction);
                    queue.Enqueue(r.ToFraction);
                }
            }
        }
        return visited;
    }

    /// <summary>Verify the structural claim: every F99 anker can reach α=0
    /// via backward references on the α-axis. α=0 is the Mirror convergence
    /// point where ±γ-polarity sides meet under (1−γ²)/2 folding — NOT a
    /// tree root (per <see cref="PolarityLayerOriginClaim"/> which types 0
    /// as the active substrate axis, not a passive midpoint).</summary>
    public bool AllAnkersConvergeToMirrorAxis()
    {
        foreach (var anker in F99Ankers)
        {
            if (Math.Abs(anker) < AnchorConstants.Tol) continue; // α=0 trivially is itself
            var closure = BackwardClosure(anker);
            if (!closure.Any(v => Math.Abs(v) < AnchorConstants.Tol)) return false;
        }
        return true;
    }

    /// <summary>Count how many distinct references (edges) exist between
    /// each fraction pair. Multi-edges = multiple viewpoints on the same
    /// fraction relationship.</summary>
    public IReadOnlyDictionary<(double, double), int> EdgeCounts()
    {
        var counts = new Dictionary<(double, double), int>();
        foreach (var r in References)
        {
            // Round to 6 decimals for grouping
            var key = (Math.Round(r.FromFraction, 6), Math.Round(r.ToFraction, 6));
            counts.TryGetValue(key, out int c);
            counts[key] = c + 1;
        }
        return counts;
    }

    /// <summary>Render the graph as a printable table grouped by source
    /// fraction, then by direction.</summary>
    public string Render()
    {
        var sb = new StringBuilder();
        sb.AppendLine();
        sb.AppendLine("Fraction Reference Graph");
        sb.AppendLine(new string('=', 88));
        sb.AppendLine();
        sb.AppendLine($"  F99 ankers: {string.Join(", ", F99Ankers.Select(AnchorConstants.FormatEighthFraction))}");
        sb.AppendLine($"  Q basis ankers: {string.Join(", ", QBasisAnkers)}");
        sb.AppendLine($"  Total edges: {References.Count}");
        sb.AppendLine();

        // Group by source fraction
        var bySource = References
            .GroupBy(r => Math.Round(r.FromFraction, 6))
            .OrderBy(g => g.Key);

        foreach (var group in bySource)
        {
            sb.AppendLine($"  From {AnchorConstants.FormatEighthFraction(group.Key)}:");
            foreach (var r in group)
            {
                string arrow = r.Direction switch
                {
                    FractionReferenceDirection.Forward  => "→",
                    FractionReferenceDirection.Backward => "↞",
                    FractionReferenceDirection.Mirror   => "↔",
                    FractionReferenceDirection.Polarity => "⇌",
                    _ => "?",
                };
                sb.AppendLine($"    {arrow} {AnchorConstants.FormatEighthFraction(r.ToFraction)}  [{r.Operation}]");
                sb.AppendLine($"      (from {r.DocumentingClaim})");
            }
            sb.AppendLine();
        }

        // Multi-edge summary
        sb.AppendLine("Multi-edge fraction pairs (more than one viewpoint):");
        var counts = EdgeCounts().OrderByDescending(kv => kv.Value);
        foreach (var (pair, count) in counts.Where(kv => kv.Value > 1))
        {
            sb.AppendLine($"  {AnchorConstants.FormatEighthFraction(pair.Item1)} → {AnchorConstants.FormatEighthFraction(pair.Item2)}: {count} viewpoints");
        }
        sb.AppendLine();

        sb.AppendLine($"Mirror-axis convergence: AllAnkersConvergeToMirrorAxis = {AllAnkersConvergeToMirrorAxis()}");
        sb.AppendLine("  (α=0 is the convergence point where ±γ-polarity sides meet under");
        sb.AppendLine("  (1−γ²)/2 folding, NOT a tree root — see PolarityLayerOriginClaim");
        sb.AppendLine("  and PolarityMirrorMap on the γ-axis for the unfolded polarity structure.)");
        return sb.ToString();
    }

}
