using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 meta-claim: the c=2 local EP at real Q_EP and the FRAGILE_BRIDGE global
/// EP in the complex-γ plane are the same exceptional-point structure under shared algebra,
/// validated empirically by the c=2 N=5..8 Petermann-K sweep on the real Q axis.
///
/// <para><b>Tier outcome: Tier2Verified.</b> The shared algebra (same-sign-imaginary 2×2
/// form, AIII chiral classification) is a Tier1-style structural fact already pinned in
/// <see cref="ChiralAiiiClassification"/>. The empirical real-axis K-divergence + parity
/// asymmetry (encoded in the four <see cref="PetermannSpikeWitness"/> entries from the
/// 2026-05-06 sweep) validates the connection: K grows monotonically with N within each
/// parity class (odd: ~1.79× per step; even: ~2.36× per step), odd-N peaks dominate even-N
/// peaks by 2–4× at every step (matches A3's σ_0 R-even/R-odd-degeneracy prediction),
/// and by N=7 the c=2 K=2384 sits ~6× above FRAGILE_BRIDGE's K=403 ballpark on the real
/// Q axis. Verdict: F86's local EP at real Q_EP is a real-axis hit of the same EP whose
/// near-singularity FRAGILE_BRIDGE detects at K=403 in the complex-γ direction. They are
/// not off-axis siblings; they are connected.</para>
///
/// <para><b>Why Tier2Verified, not Tier1Derived:</b> the explicit complex-γ analytic
/// continuation linking the two real-axis K-divergence regimes to a single off-axis EP
/// has not been carried out in code. That promotion requires either complex-γ
/// infrastructure in <c>LindbladPropagator</c> or a closed-form K(N) at the EP. The
/// <see cref="PendingDerivationNote"/> carries the gap explicitly.</para>
///
/// <para>Block-independent meta-claim: registered at the F86 KB root for any block, not
/// just c=2 (similar to <see cref="ChiralAiiiClassification"/>). The witnesses pin the
/// concrete c=2 sweep that motivated the connection; the algebraic claim itself is shared
/// across all c.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F86_QPEAK.md</c> (F86 local instance, Statement 1),
/// <c>hypotheses/FRAGILE_BRIDGE.md</c> (global instance, K=403), and
/// <c>experiments/PT_SYMMETRY_ANALYSIS.md</c> (AIII chiral classification, the shared
/// algebraic substrate).</para>
/// </summary>
public sealed class LocalGlobalEpLink : Claim
{
    /// <summary>Same-sign-imaginary 2×2 algebra, the shared algebraic object:
    /// L_eff − (trace/2)·I = [[−Δ/2, +iJ·g_eff], [+iJ·g_eff, +Δ/2]] (PROOF_F86_QPEAK
    /// Statement 1). AIII chiral classification (PT_SYMMETRY_ANALYSIS): linear Π, Π⁴=I,
    /// anti-commutation {Π, L_c}=0.</summary>
    public string SharedAlgebra => "same-sign-imaginary 2×2; AIII chiral";

    /// <summary>F86 local instance (Σγ = N·γ₀ ≠ 0, real Q_EP = 2/g_eff).</summary>
    public string LocalInstanceAnchor => "F86 Statement 1, PROOF_F86_QPEAK.md";

    /// <summary>FRAGILE_BRIDGE global instance (Σγ = 0, complex-γ-plane EP, Petermann K=403).</summary>
    public string GlobalInstanceAnchor => "hypotheses/FRAGILE_BRIDGE.md";

    /// <summary>Empirical witnesses from the c=2 N=5..8 Petermann-K sweep
    /// (γ₀=0.05, c=2 block (1,2), Q-grid 121 pts on [0.50, 4.00], 2026-05-06).
    /// One witness per N. Hard-coded as a static pinned table — this is a frozen
    /// empirical anchor, not a live-recomputed Claim.</summary>
    public IReadOnlyList<PetermannSpikeWitness> Witnesses => _petermannWitnesses;

    /// <summary>The not-yet-executed bridge: explicit complex-γ analytic continuation.
    /// Names the surviving symmetry, the structural divergence, and the promotion path
    /// to Tier1Derived.</summary>
    public string PendingDerivationNote =>
        "Surviving symmetry: AIII chiral classification — the same linear Π with Π⁴=I and " +
        "{Π, L_c}=0 anti-commutation governs both instances. Structural divergence: the " +
        "F86 local EP sits at real Q_EP = 2/g_eff with Σγ = N·γ₀ ≠ 0 (rate-channel axis), " +
        "while FRAGILE_BRIDGE's global EP sits in the complex-γ plane with Σγ = 0 " +
        "(off-real-axis Hopf-bifurcation continuation). The c=2 N=5..8 sweep witnesses " +
        "establish the real-axis K-divergence; the explicit complex-γ analytic continuation " +
        "linking these two regimes to a single off-axis EP is NOT executed in code. " +
        "Promotion to Tier1Derived requires either (a) complex-γ infrastructure in " +
        "LindbladPropagator (so the K-divergence trajectory can be traced into the complex " +
        "γ plane), or (b) a closed-form K(N) at the EP that matches FRAGILE_BRIDGE's K=403 " +
        "ballpark via known asymptotics. Until then, Tier2Verified is honest: shared " +
        "algebra + empirically validated real-axis hit, with the analytic continuation as " +
        "the documented gap.";

    private LocalGlobalEpLink()
        : base("local–global EP link (F86 ↔ FRAGILE_BRIDGE)",
               Tier.Tier2Verified,
               "docs/proofs/PROOF_F86_QPEAK.md (F86 Statement 1) + " +
               "hypotheses/FRAGILE_BRIDGE.md (Petermann K=403) + " +
               "experiments/PT_SYMMETRY_ANALYSIS.md (AIII chiral) + " +
               "docs/superpowers/syntheses/2026-05-06-petermann-sweep-c2.md (empirical witnesses)")
    { }

    /// <summary>Public factory: builds the meta-claim with the four pinned Petermann
    /// spike witnesses. No <see cref="CoherenceBlocks.CoherenceBlock"/> required — the
    /// claim is block-independent (the witnesses pin the concrete c=2 sweep, but the
    /// algebraic statement is shared across all c).</summary>
    public static LocalGlobalEpLink Build() => new();

    public override string DisplayName =>
        "F86 ↔ FRAGILE_BRIDGE: shared EP under AIII chiral algebra";

    public override string Summary =>
        $"Tier2Verified: real-axis hit empirically validated (c=2 N=5..8 Petermann K up to {_petermannWitnesses[2].MaxKGlobal:F0} at N=7, ~6× FRAGILE_BRIDGE K=403); analytic continuation pending";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("shared algebra", summary: SharedAlgebra);
            yield return new InspectableNode("local instance anchor",
                summary: LocalInstanceAnchor);
            yield return new InspectableNode("global instance anchor",
                summary: GlobalInstanceAnchor);
            yield return InspectableNode.Group(
                "Petermann-K witnesses (c=2 N=5..8 real-axis sweep)",
                _petermannWitnesses.Cast<IInspectable>().ToArray());
            yield return new InspectableNode("PendingDerivationNote",
                summary: PendingDerivationNote);
        }
    }

    /// <summary>Pinned empirical table from the 2026-05-06 c=2 N=5..8 Petermann-K sweep.
    /// γ₀=0.05, c=2 block (1,2), Q-grid 121 pts on [0.50, 4.00].
    /// Order: N=5, 6, 7, 8 (parity alternates odd/even/odd/even).</summary>
    private static readonly PetermannSpikeWitness[] _petermannWitnesses = new[]
    {
        new PetermannSpikeWitness(N: 5, Parity: "odd",  BlockDim: 50,
            MaxKGlobal: 1333.6, ArgMaxQ: 1.288,
            MaxKInterior: 1333.6, MaxKEndpoint:  252.5),
        new PetermannSpikeWitness(N: 6, Parity: "even", BlockDim: 90,
            MaxKGlobal:  337.9, ArgMaxQ: 0.938,
            MaxKInterior:  337.9, MaxKEndpoint:  131.7),
        new PetermannSpikeWitness(N: 7, Parity: "odd",  BlockDim: 147,
            MaxKGlobal: 2384.7, ArgMaxQ: 1.842,
            MaxKInterior: 2384.7, MaxKEndpoint: 1341.2),
        new PetermannSpikeWitness(N: 8, Parity: "even", BlockDim: 224,
            MaxKGlobal:  795.4, ArgMaxQ: 2.046,
            MaxKInterior:  580.0, MaxKEndpoint:  795.4),
    };
}

/// <summary>One row of the c=2 Petermann-K sweep witness table — a frozen empirical
/// data point linking the F86 local EP at real Q to FRAGILE_BRIDGE's K=403 ballpark.
///
/// <para>Hosts the per-N spike data: the global Petermann factor maximum across the
/// Q-grid, the Q at which it peaked, and the per-bond-class breakdown. The parity
/// label (odd/even N) carries A3's σ_0 R-even/R-odd-degeneracy prediction: odd-N
/// peaks empirically dominate even-N peaks by 2–4× at every step.</para>
///
/// <para>Block dimension is the dimension of the c=2 block (n=1 → n=2 coherence
/// subspace) — verifies the sweep ran on the correct block.</para>
/// </summary>
/// <param name="N">Chain length.</param>
/// <param name="Parity">"odd" or "even" — A3 parity class.</param>
/// <param name="BlockDim">Dimension of the c=2 block at this N.</param>
/// <param name="MaxKGlobal">max Petermann K across the Q-grid (max over both bond
/// classes — the global maximum).</param>
/// <param name="ArgMaxQ">Q value at which <paramref name="MaxKGlobal"/> peaks.</param>
/// <param name="MaxKInterior">max Petermann K restricted to Interior bonds.</param>
/// <param name="MaxKEndpoint">max Petermann K restricted to Endpoint bonds.</param>
public sealed record PetermannSpikeWitness(
    int N,
    string Parity,
    int BlockDim,
    double MaxKGlobal,
    double ArgMaxQ,
    double MaxKInterior,
    double MaxKEndpoint
) : IInspectable
{
    public string DisplayName => $"Petermann-K spike at c=2 N={N} ({Parity})";

    public string Summary => $"max K = {MaxKGlobal:F1} at Q = {ArgMaxQ:F3}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return InspectableNode.RealScalar("N", N);
            yield return new InspectableNode("parity", summary: Parity);
            yield return InspectableNode.RealScalar("block dim (c=2)", BlockDim);
            yield return InspectableNode.RealScalar("max K (global)", MaxKGlobal, "F1");
            yield return InspectableNode.RealScalar("argmax Q", ArgMaxQ, "F3");
            yield return InspectableNode.RealScalar("max K (Interior)", MaxKInterior, "F1");
            yield return InspectableNode.RealScalar("max K (Endpoint)", MaxKEndpoint, "F1");
        }
    }

    public InspectablePayload Payload =>
        new InspectablePayload.Real($"max K (N={N}, {Parity})", MaxKGlobal, "F1");
}
