using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F1 depolarizing-residual closed form
/// <c>‖M_F1(depol)‖²_F = 4^(N−1) · (16/9)·Σγ²</c> (Tier 1 derived in
/// <see cref="RCPsiSquared.Core.F1.F1DepolResidualClosedForm"/>) as
/// Pi2-Foundation inheritance.
///
/// <para>The centered local coefficient factors algebraically through the same
/// Pi2-Foundation primitives that
/// <see cref="F5DepolarizingErrorPi2Inheritance"/> uses for F5's scalar
/// <c>2N/3</c>; the squared Frobenius residual squares the per-Pauli rate, while
/// F5 carries the extreme pair-sum shortfall linearly. Both root in the
/// same denominator <c>d² − 1 = 3</c>.</para>
///
/// <list type="bullet">
///   <item><b>(16/9) local coefficient</b>: <c>(d²/(d²−1))² = (4/3)²</c>. The
///         per-Pauli depolarizing rate <c>d²/(d²−1) = 4/3</c> is the
///         universal-channel rate the depolarizing master equation pumps each of
///         the <c>d²−1 = 3</c> non-identity Paulis at; squaring it gives the
///         per-site Frobenius weight 16/9.</item>
///   <item><b>Zero cross-site coefficient after F1 centering.</b> The bare local
///         kernel has trace −8, but adding the F1 shift 2γ_l I makes it
///         traceless. Hence <c>tr(M_l† M_l′)=0</c> for l ≠ l′. Computed, not
///         asserted: <see cref="ResidualCoefficients"/> builds the per-site kernel
///         M_l = Π·L_l·Π⁻¹ + L_l + 2σ_l·I in integers from the same primitives
///         (L_l = −d²/(d²−1) on each non-identity Pauli, Π the per-site I↔X, Y↔Z
///         relabelling, σ_l the kernel's own spectral mean) and reads both
///         coefficients off its norm and trace. The zero is the centering at the
///         spectral mean: Π only permutes a Pauli-diagonal kernel, so
///         tr M_l = 2·tr L_l + 2d²σ_l vanishes exactly there, and the uncentered
///         kernel (σ_l = 0) gives the bare 16.</item>
/// </list>
///
/// <para><b>Tier outcome: Tier1Derived.</b> Both depol coefficients reduce to clean
/// Pi2-Foundation expressions in <c>d²</c> and <c>d²−1</c>; squaring the per-Pauli rate
/// is what cleans the algebra into Pi2 primitives. Contrast T1's analogue
/// <see cref="F1T1AmplitudeDampingPi2Inheritance"/> (Tier1Candidate, T1's <c>c_1 = 3</c>
/// is not Pi2-anchored). Contrast F5: F5's coefficient is also clean (<c>2/3</c>) but at
/// the linear scalar level; squaring it would give (4/9), not 16/9, because F5's "2" is
/// <c>d</c>-anchored and the depol closed-form "2" inside <c>d²/(d²−1)</c> is actually
/// <c>d²</c>-anchored.</para>
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
public sealed class F1DepolResidualClosedFormPi2Inheritance : Claim, IZ2AxisClaim
{

    /// <summary>The F1² / Π²_Z axis (bit_b parity, n_Y + n_Z mod 2). The
    /// canonical Pi²-Inheritance axis. The bit_a-twin (Π²_X / F61 axis) is
    /// currently not typed for this Claim.</summary>
    public Z2Axis Z2Axis => Z2Axis.BitB;

    /// <summary>The typed bit_a-twin sibling, if one exists. Currently null
    /// (no bit_a twin is typed for this Claim; this is an open slot in the
    /// cubic-architecture coverage).</summary>
    public Claim? BitATwin => null;
    public Pi2DyadicLadderClaim Ladder { get; }
    public Pi2OperatorSpaceMirrorClaim Mirror { get; }

    /// <summary>The qubit operator-space dimension squared, <c>d² = 4</c>, read as
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(−1) = <c>a_{−1}</c> on the Pi2 dyadic
    /// ladder. Matches the F1T1 sibling at
    /// <see cref="F1T1AmplitudeDampingPi2Inheritance.FourMultiplier"/> (identical
    /// <c>a_{−1} = 4</c> role). Cross-check available via
    /// <see cref="Pi2OperatorSpaceMirrorClaim.PairAt"/>(1).OperatorSpace, which feeds
    /// <see cref="DSquaredMinusOne"/>.</summary>
    public double DSquared => Ladder.Term(-1);

    /// <summary>The non-identity Pauli count per qubit, <c>d² − 1 = 3</c>. Read from
    /// <see cref="Pi2OperatorSpaceMirrorClaim.PairAt"/>(1).OperatorSpace − 1 so the
    /// operator-space-pair semantic (number of non-identity Paulis) is anchored on the
    /// mirror, parallel to <see cref="F5DepolarizingErrorPi2Inheritance.DSquaredMinusOne"/>;
    /// this is the foundation the F5 scalar error's "3" denominator and the depol
    /// closed form's "9 = 3²" denominator share.</summary>
    public double DSquaredMinusOne => (Mirror.PairAt(1)?.OperatorSpace ?? 4.0) - 1.0;

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

    /// <summary>The centered cross-site (<c>(Σγ)²</c>) coefficient, computed by
    /// <see cref="ResidualCoefficients"/>(centered: true): zero because the
    /// F1-centered local kernel is traceless.</summary>
    public double CrossSiteCoefficient => ResidualCoefficients(centered: true).Cross;

    /// <summary>The (Σγ², (Σγ)²) coefficients of ‖Σ_l M_l ⊗ I‖²_F / 4^(N−1) for the per-site
    /// depolarizing F1 kernel, computed from the Pi2 primitives in exact integer arithmetic
    /// (every entry scaled by d² − 1 = 3). The kernel is Pauli-diagonal, so
    /// ‖Σ_l M_l ⊗ I‖² = 4^(N−1)·Σ_l ‖M_l‖² + 4^(N−2)·Σ_{l≠l′} tr M_l·tr M_l′, which with
    /// M_l = γ_l·m gives Cross = (tr m)²/d² and Local = ‖m‖² − (tr m)²/d².
    /// <paramref name="centered"/> = true shifts by 2σ_l, σ_l = the kernel's spectral
    /// mean (γ_l for depolarizing); false is the bare residual (σ = 0).</summary>
    public (double Local, double Cross) ResidualCoefficients(bool centered)
    {
        int dSq = checked((int)DSquared), nonId = checked((int)DSquaredMinusOne);
        // L_l on (I, X, Y, Z) in units γ_l/(d²−1): 0 on I, −d² on each non-identity Pauli
        var l = new[] { 0, -dSq, -dSq, -dSq };
        // Π per site: I ↔ X, Y ↔ Z; on a Pauli-diagonal kernel conjugation only permutes entries
        var pi = new[] { 1, 0, 3, 2 };
        // σ_l·(d²−1) = −tr L_l / d² ·(d²−1) = the mean rate, integral here (= d²−1)
        int traceL = l.Sum();
        if (traceL % dSq != 0) throw new InvalidOperationException("spectral mean not integral in these units");
        int twoSigma = centered ? -2 * traceL / dSq : 0;
        var m = Enumerable.Range(0, 4).Select(p => l[pi[p]] + l[p] + twoSigma).ToArray();
        long norm = m.Sum(x => (long)x * x), trace = m.Sum();
        double scale = (double)nonId * nonId;
        double cross = trace * trace / (double)dSq / scale;
        double local = (norm - trace * trace / (double)dSq) / scale;
        return (local, cross);
    }

    /// <summary>Drift-check alias for <see cref="LocalCoefficient"/>. Tests compare this
    /// against the parent's hardcoded
    /// <see cref="RCPsiSquared.Core.F1.F1DepolResidualClosedForm.LocalCoefficient"/> constant
    /// to catch divergence between the algebraic Pi2-derived value and the parent's pinned
    /// number. The forward is intentional: the Pi2 derivation IS the live recomposition,
    /// the parent's constant is the independent reference.</summary>
    public double LiveLocalCoefficient => LocalCoefficient;

    /// <summary>Drift-check alias for <see cref="CrossSiteCoefficient"/>; compared against
    /// <see cref="RCPsiSquared.Core.F1.F1DepolResidualClosedForm.CrossSiteCoefficient"/>
    /// in tests. Same Pi2-vs-parent split as <see cref="LiveLocalCoefficient"/>.</summary>
    public double LiveCrossSiteCoefficient => CrossSiteCoefficient;

    public F1DepolResidualClosedFormPi2Inheritance(Pi2DyadicLadderClaim ladder, Pi2OperatorSpaceMirrorClaim mirror)
        : base("F1 centered depol-residual coefficient 16/9 and zero cross-site term inherit from Pi2-Foundation",
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
        "F1 centered depol-residual coefficients (16/9, 0) as Pi2-Foundation inheritance";

    public override string Summary =>
        $"local = (d²/(d²−1))² = ({DSquared}/{DSquaredMinusOne})² = {LocalCoefficient:F4}; " +
        $"cross = {CrossSiteCoefficient} after F1 centering (bare kernel: {ResidualCoefficients(centered: false).Cross}), " +
        $"both computed from the per-site kernel; per-Pauli rate = d²/(d²−1) = {PerPauliDepolarizingRate:F4} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F1 depol closed form",
                summary: "‖M_F1(depol)‖²_F = 4^(N−1)·(16/9)·Σγ² (Tier 1 derived; σ=Σγ; H-independent; per-site only)");
            yield return InspectableNode.RealScalar("DSquared (= a_{-1} = d²)", DSquared);
            yield return InspectableNode.RealScalar("DSquaredMinusOne (= d² − 1)", DSquaredMinusOne);
            yield return InspectableNode.RealScalar("PerPauliDepolarizingRate (= d²/(d²−1) = 4/3)", PerPauliDepolarizingRate);
            yield return InspectableNode.RealScalar("LocalCoefficient (= (d²/(d²−1))² = 16/9)", LocalCoefficient);
            yield return InspectableNode.RealScalar("CrossSiteCoefficient (= 0 after F1 centering)", CrossSiteCoefficient);
            yield return new InspectableNode("F5 sibling comparison",
                summary: "F5 uses d/(d²−1) = 2/3 for the linear scalar error coefficient; this claim uses (d²/(d²−1))² = 16/9 for the squared Frobenius residual local coefficient. Both root in DSquaredMinusOne = 3 from the same Pi2-Foundation primitive.");
        }
    }
}
