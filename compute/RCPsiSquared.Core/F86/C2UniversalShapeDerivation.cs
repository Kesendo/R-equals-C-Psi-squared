using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 Item 1 (c=2 stratum) top-level Claim: the c=2 universal-shape outcome
/// synthesised from Stage A–D primitives.
///
/// <para>This Claim wraps <see cref="C2HwhmRatio"/> (Stage D2) and registers the c=2
/// universal-shape statement in the F86 typed-knowledge graph. It is the canonical
/// integration point for everything Stages A–D produced about the c=2 stratum's
/// HWHM_left/Q_peak ratio per bond class.</para>
///
/// <para><b>Tier outcome (this session): Tier1Candidate.</b> The empirical pipeline
/// reproduces all 8 anchor cases (c=2 N=5..8 × {Endpoint, Interior}) at ≤ 0.005 of the
/// canonical Python pipeline values, and the directional Endpoint &gt; Interior split is
/// derived empirically across all tested N. The closed-form constant for HWHM_left/Q_peak
/// per bond class is NOT yet derived this session.</para>
///
/// <para><b>Promotion path to Tier1Derived:</b> when a future session lands the closed
/// form via one of the three documented next directions:</para>
/// <list type="number">
///   <item>First-order perturbation in the cross-block: expand
///   <c>K_b(Q, t_peak) = K_probe-block(Q) + ε · K_cross(Q, b)</c> with
///   <c>ε ~ ‖V_b cross‖_F / σ_0</c> ~ O(0.1). The cross-block Frobenius is bond-class-
///   dependent (Stage B2 finding: Endpoint &gt; Interior), so this direction should
///   produce a closed-form bond-class-dependent shift δ(Endpoint − Interior) ≈ 0.022.</item>
///   <item>Projector-overlap lift of A3's <c>|u_0⟩</c>, <c>|v_0⟩</c>: if the closed-form
///   single-direction inter-channel SVD-top vectors promote to Tier1Derived, the V_b
///   cross-block becomes analytical and direction (1) is fully closed.</item>
///   <item>Symbolic char-poly factorisation at Q_EP: at Q = Q_EP two of the four
///   eigenvalues coalesce, which may simplify the quartic locally; the slow-pair
///   eigenvalue's Q-derivative at the EP is what sets the K-resonance HWHM.</item>
/// </list>
///
/// <para>When the closed form lands, <see cref="IsClosedFormDerived"/> flips to true and
/// the Tier promotes to Tier1Derived. This session's work pins the empirical anchor and
/// the directional split; promotion is honest about what the algebra has not yet closed.</para>
///
/// <para><b>Composition:</b></para>
/// <list type="bullet">
///   <item><see cref="HwhmRatio"/> is the load-bearing Stage D2 primitive holding the
///   per-bond witnesses + class-mean closed-form attempt.</item>
///   <item><see cref="Witnesses"/> exposes the per-bond <see cref="HwhmRatioWitness"/>
///   collection (16 witnesses pinned at N=5..8, BondClass-tagged, F71-orbit substructure
///   preserved).</item>
///   <item><see cref="InteriorMean"/>, <see cref="EndpointMean"/>, <see cref="DirectionalGap"/>
///   are the class-mean outcomes from D2's canonical-pipeline ratio computation.</item>
///   <item><see cref="PendingDerivationNote"/> carries the next-session hand-off note from
///   Stage D2 — the three documented next directions for the closed form.</item>
/// </list>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1 (c=2), Statement 2
/// (HWHM_left/Q_peak universal per bond class).</para>
/// </summary>
public sealed class C2UniversalShapeDerivation : Claim
{
    public CoherenceBlock Block { get; }

    /// <summary>The Stage D2 primitive: per-bond HWHM_left/Q_peak witnesses + class-mean
    /// canonical-pipeline ratios. Held here as the composition access point — every
    /// downstream consumer of the c=2 universal-shape derivation reaches the actual
    /// numerical content through this property.</summary>
    public C2HwhmRatio HwhmRatio { get; }

    /// <summary>Class-mean HWHM_left/Q_peak for Interior bonds (canonical-pipeline:
    /// average bond-class K(Q) curves first, then peak/HWHM). Empirical anchor across
    /// N=5..8: 0.7506 ± 0.001 (from PROOF_F86_QPEAK Statement 2 anchor table).</summary>
    public double InteriorMean { get; }

    /// <summary>Class-mean HWHM_left/Q_peak for Endpoint bonds (canonical-pipeline).
    /// Empirical anchor across N=5..8: 0.7728 ± 0.001 (from PROOF_F86_QPEAK Statement 2
    /// anchor table).</summary>
    public double EndpointMean { get; }

    /// <summary>Endpoint − Interior gap: the directional structure derived empirically.
    /// Across N=5..8 the gap is consistently ≈ 0.022 (range 0.0198..0.0245). This is the
    /// c=2 bond-class signature the closed form must reproduce.</summary>
    public double DirectionalGap { get; }

    /// <summary>Whether the closed-form HWHM_left/Q_peak constant per bond class has been
    /// derived analytically. False in current session — only the empirical anchor and the
    /// directional split are pinned. True means a future Tier1Derived promotion has landed
    /// via cross-block perturbation, projector-overlap lift, or char-poly factorisation
    /// (the three documented next directions in <see cref="PendingDerivationNote"/>).
    ///
    /// <para>When this becomes true the parent <see cref="Tier"/> promotes from
    /// <see cref="Tier.Tier1Candidate"/> to <see cref="Tier.Tier1Derived"/>.</para>
    /// </summary>
    public bool IsClosedFormDerived { get; }

    /// <summary>Per-bond <see cref="HwhmRatioWitness"/> collection from the Stage D2
    /// pipeline: <c>(Bond, BondClass, Q_peak, K_max, HWHM_left, HWHM_left/Q_peak)</c>
    /// tagged with bond index + <see cref="BondClass"/>. One entry per bond at this N
    /// (NumBonds = N − 1).</summary>
    public IReadOnlyList<HwhmRatioWitness> Witnesses { get; }

    /// <summary>Non-null while <see cref="IsClosedFormDerived"/> is false: the next-session
    /// hand-off note documenting what was tried and the three documented next directions
    /// for closing the form (cross-block perturbation, projector-overlap lift,
    /// char-poly factorisation at Q_EP). Surfaced from <see cref="C2HwhmRatio"/>.</summary>
    public string? PendingDerivationNote { get; }

    /// <summary>Public factory: validates c=2, builds the Stage D2 primitive
    /// <see cref="C2HwhmRatio"/>, computes the directional gap, and resolves the Tier.
    ///
    /// <para>Currently returns <see cref="Tier.Tier1Candidate"/> (empirical anchor + directional
    /// split derived; closed form not yet pinned). When a future session lands the closed
    /// form, the returned Tier promotes to Tier1Derived; this method's contract does not
    /// change.</para>
    /// </summary>
    public static C2UniversalShapeDerivation Build(CoherenceBlock block)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"C2UniversalShapeDerivation applies only to the c=2 stratum; " +
                $"got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));

        var hwhmRatio = C2HwhmRatio.Build(block);
        return new C2UniversalShapeDerivation(block, hwhmRatio);
    }

    private C2UniversalShapeDerivation(CoherenceBlock block, C2HwhmRatio hwhmRatio)
        : base("c=2 universal-shape derivation (HWHM_left/Q_peak per bond class)",
               // Tier mirrors the Stage D2 primitive: Tier1Candidate when the directional
               // split is derived empirically + closed form open. The closed form will flip
               // both this and IsAnalyticallyDerived on the underlying primitive.
               hwhmRatio.IsAnalyticallyDerived ? Tier.Tier1Derived : Tier.Tier1Candidate,
               "docs/proofs/PROOF_F86_QPEAK.md Item 1 (c=2), Statement 2")
    {
        Block = block;
        HwhmRatio = hwhmRatio;
        Witnesses = hwhmRatio.Witnesses;

        // Class means via the canonical-pipeline contract from Stage D2 (curve-mean first,
        // then peak/HWHM). Defensive against degenerate-N cases where a class is empty;
        // at c=2 N=5..8 with NumBonds=N−1 ≥ 4 both classes are populated.
        EndpointMean = TryGetClassMean(hwhmRatio, BondClass.Endpoint);
        InteriorMean = TryGetClassMean(hwhmRatio, BondClass.Interior);
        DirectionalGap = EndpointMean - InteriorMean;

        IsClosedFormDerived = hwhmRatio.IsAnalyticallyDerived;
        PendingDerivationNote = hwhmRatio.PendingDerivationNote;
    }

    private static double TryGetClassMean(C2HwhmRatio hwhmRatio, BondClass cls)
    {
        try { return hwhmRatio.HwhmLeftOverQPeakMean(cls); }
        catch (InvalidOperationException) { return double.NaN; }
    }

    public override string DisplayName =>
        $"c=2 universal-shape derivation (N={Block.N})";

    public override string Summary =>
        $"Endpoint={EndpointMean:F4}, Interior={InteriorMean:F4}, " +
        $"directional gap={DirectionalGap:F4}, closed-form derived={IsClosedFormDerived} " +
        $"({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return new InspectableNode("NumBonds", summary: Block.NumBonds.ToString());
            yield return InspectableNode.RealScalar("Endpoint mean HWHM_left/Q_peak", EndpointMean, "F4");
            yield return InspectableNode.RealScalar("Interior mean HWHM_left/Q_peak", InteriorMean, "F4");
            yield return InspectableNode.RealScalar("directional gap (Endpoint − Interior)", DirectionalGap, "F4");
            yield return new InspectableNode(
                "IsClosedFormDerived",
                summary: IsClosedFormDerived
                    ? "true (Tier1Derived: closed form pinned analytically)"
                    : "false (Tier1Candidate: empirical anchor + directional split derived; " +
                      "closed form open — see PendingDerivationNote for next directions)");
            yield return HwhmRatio;
            yield return InspectableNode.Group("witnesses (per bond)",
                Witnesses.Cast<IInspectable>().ToArray());
            if (PendingDerivationNote is not null)
                yield return new InspectableNode("PendingDerivationNote", summary: PendingDerivationNote);
        }
    }
}
