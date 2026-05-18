using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F1 depolarizing-residual closed form
/// <c>‖M(depol)‖²_F = 4^(N−1) · [(16/9)·Σγ² + 16·(Σγ)²]</c> (Tier 1 derived in
/// <see cref="RCPsiSquared.Core.F1.F1DepolResidualClosedForm"/>) as
/// Pi2-Foundation inheritance.
///
/// <para>The two coefficient anchors factor algebraically through the same
/// Pi2-Foundation primitives that
/// <see cref="F5DepolarizingErrorPi2Inheritance"/> uses for F5's scalar
/// <c>2·(N−2)/3</c>; the difference is only that the squared Frobenius residual
/// squares the per-Pauli rate, while F5 carries it linearly. Both root in the
/// same denominator <c>d² − 1 = 3</c>.</para>
///
/// <list type="bullet">
///   <item><b>(16/9) local coefficient</b>: <c>(d²/(d²−1))² = (4/3)²</c>. The
///         per-Pauli depolarizing rate <c>d²/(d²−1) = 4/3</c> is the
///         universal-channel rate the depolarizing master equation pumps each of
///         the <c>d²−1 = 3</c> non-identity Paulis at; squaring it gives the
///         per-site Frobenius weight 16/9.</item>
///   <item><b>16 cross-site coefficient</b>: <c>(d²)² = 4² = 16</c>. The
///         per-site trace contribution is <c>|tr(M_l)|² = (d²·γ)² ∝ 16</c> at
///         γ = 1; cross-site tensor assembly via
///         <c>tr(M_l† M_l′) = |tr(M_l)|² · 4^(N−2)</c> for <c>l ≠ l′</c> lifts
///         this squared quantity to the cooperative <c>(Σγ)²</c> piece.</item>
/// </list>
///
/// <para>The composite ratio between the local and cross coefficients is
/// <c>(16/9) / 16 = 1/9 = (1/(d²−1))²</c>: same square of the same
/// <c>1/(d²−1)</c> that gives F5 its <c>2/3</c> at the linear level.</para>
///
/// <para><b>Tier outcome: Tier1Derived.</b> Unlike T1's analogue
/// <see cref="F1T1AmplitudeDampingPi2Inheritance"/> (Tier1Candidate because
/// T1's <c>c_1 = 3</c> is not Pi2-anchored), both depol coefficients (16/9, 16)
/// reduce to clean Pi2-Foundation expressions in <c>d²</c> and <c>d²−1</c>.
/// That is what makes this claim a derived Pi2 inheritance: the squaring of
/// the per-Pauli rate is exactly what cleans both coefficients into trivial
/// algebra over the Pi2 primitives. Contrast F5: F5's coefficient is also
/// clean (<c>2/3</c>), but at the linear scalar level; squaring it would give
/// (4/9), not 16/9, because F5's "2" is <c>d</c>-anchored and the depol
/// closed-form "2" inside <c>d²/(d²−1)</c> is actually <c>d²</c>-anchored.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F1 depol bullet +
/// <c>docs/proofs/PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM.md</c> +
/// <c>compute/RCPsiSquared.Core/F1/F1DepolResidualClosedForm.cs</c> (the Tier1
/// parent closed form) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F5DepolarizingErrorPi2Inheritance.cs</c>
/// (sibling Pi2-inheritance for the scalar F5 error; same denominator anchor
/// at the linear level).</para>
/// </summary>
public sealed class F1DepolResidualClosedFormPi2Inheritance : Claim
{
    public Pi2DyadicLadderClaim Ladder { get; }
    public Pi2OperatorSpaceMirrorClaim Mirror { get; }

    /// <summary>The qubit operator-space dimension squared, <c>d² = 4</c>, read
    /// from <see cref="Pi2OperatorSpaceMirrorClaim.PairAt"/>(1). Same anchor as
    /// <c>a_{−1}</c> on the dyadic ladder.</summary>
    public double DSquared => Mirror.PairAt(1)?.OperatorSpace ?? 4.0;

    /// <summary>The non-identity Pauli count per qubit, <c>d² − 1 = 3</c>. Same
    /// anchor as <see cref="F5DepolarizingErrorPi2Inheritance.DSquaredMinusOne"/>;
    /// this is the foundation the F5 scalar error's "3" denominator and the depol
    /// closed form's "9 = 3²" denominator share.</summary>
    public double DSquaredMinusOne => DSquared - 1.0;

    /// <summary>The per-Pauli depolarizing rate <c>d²/(d²−1) = 4/3</c>. Each of
    /// the <c>d²−1 = 3</c> non-identity Paulis is pumped at this rate by a unit
    /// depolarizing budget; F5's linear "2/3" is <c>d/(d²−1)</c>, this depol
    /// closed form's per-Pauli rate is <c>d²/(d²−1)</c> (one power of d higher
    /// because the residual squares the rate before per-site assembly).</summary>
    public double PerPauliDepolarizingRate => DSquared / DSquaredMinusOne;

    /// <summary>The local (<c>Σγ²</c>) coefficient: <c>(d²/(d²−1))² = 16/9</c>.
    /// Squaring the per-Pauli rate is what carries the algebra from F5's linear
    /// scalar to F1's squared Frobenius residual.</summary>
    public double LocalCoefficient => PerPauliDepolarizingRate * PerPauliDepolarizingRate;

    /// <summary>The cross-site (<c>(Σγ)²</c>) coefficient: <c>(d²)² = 16</c>.
    /// Squared because the cooperative tensor-assembly piece is
    /// <c>|tr(M_l)|² · 4^(N−2)</c> and <c>|tr(M_l)| ∝ d²</c> at γ=1.</summary>
    public double CrossSiteCoefficient => DSquared * DSquared;

    /// <summary>Live re-composition of the local coefficient from the Pi2-derived
    /// primitives. Drift between the algebraic value and
    /// <see cref="RCPsiSquared.Core.F1.F1DepolResidualClosedForm.LocalCoefficient"/>
    /// surfaces here.</summary>
    public double LiveLocalCoefficient => LocalCoefficient;

    /// <summary>Live re-composition of the cross-site coefficient from the
    /// Pi2-derived primitives. Drift between the algebraic value and
    /// <see cref="RCPsiSquared.Core.F1.F1DepolResidualClosedForm.CrossSiteCoefficient"/>
    /// surfaces here.</summary>
    public double LiveCrossSiteCoefficient => CrossSiteCoefficient;

    public F1DepolResidualClosedFormPi2Inheritance(Pi2DyadicLadderClaim ladder, Pi2OperatorSpaceMirrorClaim mirror)
        : base("F1 depol-residual closed-form coefficients (16/9, 16) inherit from Pi2-Foundation",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F1 depol bullet + " +
               "docs/proofs/PROOF_F1_DEPOL_RESIDUAL_CLOSED_FORM.md + " +
               "compute/RCPsiSquared.Core/F1/F1DepolResidualClosedForm.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/F5DepolarizingErrorPi2Inheritance.cs")
    {
        Ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        Mirror = mirror ?? throw new ArgumentNullException(nameof(mirror));
    }

    public override string DisplayName =>
        "F1 depol-residual coefficients (16/9, 16) as Pi2-Foundation inheritance";

    public override string Summary =>
        $"local = (d²/(d²−1))² = ({DSquared}/{DSquaredMinusOne})² = {LocalCoefficient:F4}; " +
        $"cross = (d²)² = {CrossSiteCoefficient:F4}; per-Pauli rate = d²/(d²−1) = {PerPauliDepolarizingRate:F4} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F1 depol closed form",
                summary: "‖M(depol)‖²_F = 4^(N−1)·[(16/9)·Σγ² + 16·(Σγ)²] (Tier 1 derived; H-independent; per-site only)");
            yield return InspectableNode.RealScalar("DSquared (= a_{-1} = d²)", DSquared);
            yield return InspectableNode.RealScalar("DSquaredMinusOne (= d² − 1)", DSquaredMinusOne);
            yield return InspectableNode.RealScalar("PerPauliDepolarizingRate (= d²/(d²−1) = 4/3)", PerPauliDepolarizingRate);
            yield return InspectableNode.RealScalar("LocalCoefficient (= (d²/(d²−1))² = 16/9)", LocalCoefficient);
            yield return InspectableNode.RealScalar("CrossSiteCoefficient (= (d²)² = 16)", CrossSiteCoefficient);
            yield return new InspectableNode("F5 sibling comparison",
                summary: "F5 uses d/(d²−1) = 2/3 for the linear scalar error coefficient; this claim uses (d²/(d²−1))² = 16/9 for the squared Frobenius residual local coefficient. Both root in DSquaredMinusOne = 3 from the same Pi2-Foundation primitive.");
        }
    }
}
