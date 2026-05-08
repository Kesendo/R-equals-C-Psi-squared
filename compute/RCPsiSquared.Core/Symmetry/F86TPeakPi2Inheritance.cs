using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F86 Statement 1's <c>t_peak = 1/(4γ₀)</c> has its "4" denominator sitting
/// on the Pi2-Foundation: exactly <c>a_{−1}</c> on the dyadic halving ladder — the
/// upper-side N=1 anchor (<c>4 = 4^1 = d²</c> for one qubit, the cardinality of the
/// single-qubit Pauli basis {I, X, Y, Z}).
///
/// <para>So <c>t_peak = 1/(d²·γ₀)</c> — the peak time is the inverse of (operator
/// space dimension for one qubit) times the dephasing rate. The "4" is not a free
/// parameter; it is the operator-space side of the Pi2 ladder for N=1.</para>
///
/// <para>Mirror reading: 1/4 = <c>a_3</c> on the lower (memory) side of the ladder
/// (= <see cref="QuarterAsBilinearMaxvalClaim"/>), so <c>t_peak = a_3 · (1/γ₀)</c>
/// is the alternate composition. Both readings are consistent — they are mirror
/// partners of one ladder fact.</para>
///
/// <para>Tier1Derived: pure composition. F86 Statement 1 is Tier1Derived in
/// <c>docs/proofs/PROOF_F86_QPEAK.md</c>; <see cref="TPeakLaw"/> wraps it as a typed
/// claim parameterised by γ₀.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Statement 1 +
/// <c>docs/ANALYTICAL_FORMULAS.md</c> F86 +
/// <c>compute/RCPsiSquared.Core/F86/TPeakLaw.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F86TPeakPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly TPeakLaw _tPeak;

    /// <summary>The "4" denominator in F86's <c>t_peak = 1/(4γ₀)</c>. Exactly equal to
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(−1) = <c>a_{−1}</c> = 4 = d² for N=1
    /// (operator-space side of the ladder).</summary>
    public double FourFactor => _ladder.Term(-1);

    /// <summary>The mirror partner: <c>1/4 = a_3</c> on the memory side. Equivalent
    /// reading via the inversion symmetry <c>a_n · a_{2−n} = 1</c>.</summary>
    public double OneOverFourFactor => _ladder.Term(3);

    /// <summary>The γ₀ used in the parent <see cref="TPeakLaw"/> instance.</summary>
    public double GammaZero => _tPeak.GammaZero;

    /// <summary>The composed t_peak value: <c>1 / (FourFactor · GammaZero)</c>;
    /// bit-exactly equal to <see cref="TPeakLaw.Value"/>.</summary>
    public double LiveTPeak => 1.0 / (FourFactor * GammaZero);

    public F86TPeakPi2Inheritance(Pi2DyadicLadderClaim ladder, TPeakLaw tPeak)
        : base("F86 t_peak's 4 denominator inherits from Pi2-Foundation (4 = a_{-1} = d² for N=1)",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F86_QPEAK.md (Statement 1) + " +
               "docs/ANALYTICAL_FORMULAS.md F86 + " +
               "compute/RCPsiSquared.Core/F86/TPeakLaw.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _tPeak = tPeak ?? throw new ArgumentNullException(nameof(tPeak));
    }

    public override string DisplayName =>
        "F86 t_peak's 4 denominator as Pi2-Foundation a_{-1}";

    public override string Summary =>
        $"t_peak = 1/(4γ₀): 4 = a_{{-1}} = {FourFactor} (= d² for N=1, operator-space side); mirror 1/4 = a_3 = {OneOverFourFactor} (memory side); live t_peak = {LiveTPeak:G6} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F86 Statement 1",
                summary: "t_peak = 1/(4γ₀) (Tier1Derived in PROOF_F86_QPEAK; universal across c, N, n, bond)");
            yield return InspectableNode.RealScalar("FourFactor (= a_{-1} = d² for N=1)", FourFactor);
            yield return InspectableNode.RealScalar("OneOverFourFactor (= a_3 mirror partner)", OneOverFourFactor);
            yield return InspectableNode.RealScalar("γ₀", GammaZero);
            yield return InspectableNode.RealScalar("LiveTPeak", LiveTPeak);
            yield return InspectableNode.RealScalar("TPeakLaw.Value (parent's pinned)", _tPeak.Value);
            yield return new InspectableNode("inheritance",
                summary: "the 4 in the denominator is not a free parameter; it is the operator-space dimension d² for N=1, on the upper side of the Pi2 ladder");
            yield return new InspectableNode("mirror reading",
                summary: "via inversion a_n · a_{2-n} = 1: t_peak = a_3 / γ₀ = (1/4)/γ₀; same fact, memory side");
        }
    }
}
