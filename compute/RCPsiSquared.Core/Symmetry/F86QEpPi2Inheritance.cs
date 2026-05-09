using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F86a's <c>Q_EP = 2/g_eff</c> has its "2" numerator sitting on
/// the Pi2-Foundation: exactly <c>a_0</c> on the dyadic halving ladder, the qubit
/// dimension d. So Q_EP = d/g_eff — the EP-position is the qubit dimension divided
/// by the effective coupling.
///
/// <para>The g_eff is the SVD-top inter-channel singular value σ_0 (parameterised
/// per (c, N)). The "2" in the numerator is universal and constant — it is the
/// Pi2-Foundation's <c>a_0</c>. F86 inherits the numerator from the Pi2-Foundation;
/// the per-(c, N) variation lives entirely in g_eff.</para>
///
/// <para>Tier1Derived: pure composition. F86a is Tier1Derived in
/// <c>docs/proofs/PROOF_F86_QPEAK.md</c>; <see cref="QEpLaw"/> wraps it as a typed
/// claim parameterised by g_eff. This claim makes the numerator's
/// Pi2-Foundation inheritance explicit.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Statement 1 +
/// <c>docs/ANALYTICAL_FORMULAS.md</c> F86a +
/// <c>compute/RCPsiSquared.Core/F86/QEpLaw.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F86QEpPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly QEpLaw _qEp;

    /// <summary>The "2" numerator in F86's <c>Q_EP = 2/g_eff</c>. Exactly equal to
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(0) = <c>a_0</c> = d.</summary>
    public double TwoFactor => _ladder.Term(0);

    /// <summary>The g_eff used in the parent <see cref="QEpLaw"/> instance.</summary>
    public double GEff => _qEp.GEff;

    /// <summary>The composed Q_EP value: <c>TwoFactor / GEff</c>; bit-exactly equal
    /// to <see cref="QEpLaw.Value"/>.</summary>
    public double LiveQEp => TwoFactor / GEff;

    public F86QEpPi2Inheritance(Pi2DyadicLadderClaim ladder, QEpLaw qEp)
        : base("F86 Q_EP's 2 numerator inherits from Pi2-Foundation (2 = a_0)",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F86_QPEAK.md (Statement 1) + " +
               "docs/ANALYTICAL_FORMULAS.md F86 + " +
               "compute/RCPsiSquared.Core/F86/QEpLaw.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _qEp = qEp ?? throw new ArgumentNullException(nameof(qEp));
    }

    public override string DisplayName =>
        "F86 Q_EP's 2 numerator as Pi2-Foundation a_0";

    public override string Summary =>
        $"Q_EP = 2/g_eff = a_0/g_eff = {TwoFactor}/{GEff:G6} = {LiveQEp:G6} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F86a",
                summary: "Q_EP = 2/g_eff (Tier1Derived in PROOF_F86_QPEAK)");
            yield return InspectableNode.RealScalar("TwoFactor (= a_0 = d)", TwoFactor);
            yield return InspectableNode.RealScalar("g_eff", GEff);
            yield return InspectableNode.RealScalar("LiveQEp (= TwoFactor / g_eff)", LiveQEp);
            yield return InspectableNode.RealScalar("QEpLaw.Value (parent's pinned)", _qEp.Value);
            yield return new InspectableNode("inheritance",
                summary: "the 2 in the numerator is not a free parameter; it is the qubit dimension on the Pi2 ladder");
        }
    }
}
