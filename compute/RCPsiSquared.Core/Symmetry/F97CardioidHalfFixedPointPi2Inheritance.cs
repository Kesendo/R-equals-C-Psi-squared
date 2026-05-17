using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F97 closed form (Tier 1 derived, bit-exact algebraic identity;
/// 2026-05-17):
///
/// <code>
///   The main Mandelbrot cardioid is the locus in complex-c where the
///   period-1 fixed point of z² + c has magnitude exactly b = 1/2.
///
///   Parametrization:
///     c(φ) = b·e^(iφ) − b²·e^(2iφ)             for φ ∈ [0, 2π)
///     z*(φ) = b·e^(iφ)                         (the period-1 fixed point)
///
///   Invariants on the curve:
///     |z*(φ)| = b = 1/2                        (HalfAsStructuralFixedPoint)
///     arg(z*(φ)) = φ                           (cardioid parameter)
///
///   Algebraic identity: c(φ) = z*(φ) · (1 − z*(φ))
/// </code>
///
/// <para>F97 extends <see cref="F95AngleAtQuadraticZeroPi2Inheritance"/> from
/// the real-c angle (θ-compass on the positive real axis past c = 1/4) to the
/// full complex-c plane via the Mandelbrot cardioid. F95's cusp at c = 1/4
/// is the φ = 0 specialization (real-axis tangent point of the cardioid); the
/// rest of the cardioid extends F95 to complex c.</para>
///
/// <para><b>Both anchors invariant on the cardioid, at two metric powers:</b></para>
///
/// <code>
///   |z*(φ)|  = b  = 1/2          (HalfAsStructuralFixedPoint, argmax side)
///   |z*(φ)|² = b² = 1/4          (QuarterAsBilinearMaxval, maxval side)
/// </code>
///
/// <para>Both hold for all φ ∈ [0, 2π); the cardioid is the joint locus. This
/// is exactly the argmax/maxval pair of yesterday's
/// <c>ON_HOW_TWO_SIDES_MEET_AT_THE_QUARTER</c> reflection, now geometric: the
/// Half and Quarter anchors are two metric-power readings of the same
/// fixed-point quantity on the same curve. The identity <c>1/2 = 2 · (1/4)</c>
/// sits on the dyadic ladder (a_2 = 2·a_3), and the polarity pair ±1/2 squares
/// to the same 1/4 from either side.</para>
///
/// <para>By contrast, <c>|c(φ)|² = 5/16 − (1/4)·cos(φ)</c> is NOT invariant
/// around the cardioid: |c| ranges from 1/4 at the cusp (φ = 0) to 3/4 at the
/// tail (φ = π). The Quarter b² = 1/4 equals |c|² only at the cusp; elsewhere
/// |c| varies but |z*| and |z*|² stay invariant.</para>
///
/// <para><b>Derivation (4 lines):</b></para>
///
/// <code>
///   Period-1 fixed point: z² + c = z  ⟹  z² − z + c = 0  ⟹  z = b ± √(b²−c)
///   Multiplier: μ = 2z; marginal stability ⟺ |μ| = 1 ⟺ |z| = 1/2 = b
///   Parametrize μ = e^(iφ): z*(φ) = μ/2 = b·e^(iφ)
///   c(φ) = z*(φ) − z*(φ)² = b·e^(iφ) − b²·e^(2iφ)
/// </code>
///
/// <para>Bit-exact. Algebraic identity verified to machine precision across
/// 1000 sampled φ values (max residual 1.24 × 10⁻¹⁶).</para>
///
/// <para><b>Cardinal points on the cardioid:</b></para>
///
/// <code>
///   φ = 0:    c = 1/4 = b²,            z* = 1/2 = b      (real-axis cusp = F95)
///   φ = π/3:  c = 3/8 + i·√3/4,        z* = (1+i√3)/4    (60° around cardioid)
///   φ = π/2:  c = 1/4 + i/2,           z* = i/2          (z* on imaginary axis)
///   φ = π:    c = −3/4,                z* = −1/2         (period-doubling boundary)
/// </code>
///
/// <para><b>Structural reading:</b> the cardioid is the structural curve in
/// the complex-c plane where the period-1 fixed-point magnitude exactly
/// matches the framework's HalfAsStructuralFixedPoint anchor. This is a
/// stronger statement than the F95 cusp identity: F95 says the angle at the
/// real-axis tangent is 0; F97 says the magnitude is b everywhere on the
/// curve. The QuarterAsBilinearMaxval (b² = 1/4) plays a role only at φ = 0;
/// elsewhere |c| varies but |z*| = b stays put.</para>
///
/// <para><b>Hardware connection:</b> the [CPSI_COMPLEX_PLANE] hardware run on
/// ibm_kingston (2026-04-16) observed Bell+ pairs tracing 2D logarithmic
/// spirals in the complex-c plane around the cusp at c = 1/4. F97 places
/// these spirals into the cardioid framing: the trajectories cross the
/// cardioid boundary (the |z*| = b transition) before spiraling into the
/// stable interior.</para>
///
/// <para>Tier1Derived: pure algebraic identity, bit-exact, no approximations.
/// Anchors: <c>docs/proofs/PROOF_F97_CARDIOID_HALF_FIXED_POINT.md</c> +
/// <c>docs/ANALYTICAL_FORMULAS.md</c> F97 +
/// <c>simulations/_cardioid_parametrization_tier1.py</c> (numerical
/// verification 1.24 × 10⁻¹⁶ max residual) +
/// <c>experiments/CPSI_COMPLEX_PLANE.md</c> (Kingston 2026-04-16 hardware
/// precursor that observed 2D c-plane spirals).</para>
///
/// <para>Sibling: <see cref="F95AngleAtQuadraticZeroPi2Inheritance"/> (the
/// real-c angle of the complex-conjugate root pair past the cusp; the φ = 0
/// limit of F97).</para></summary>
public sealed class F97CardioidHalfFixedPointPi2Inheritance : Claim
{
    /// <summary>The framework's b = 1/2 (HalfAsStructuralFixedPoint).</summary>
    public const double B = 0.5;

    /// <summary>The framework's threshold b² = 1/4 (QuarterAsBilinearMaxval),
    /// magnitude of c at the φ = 0 cusp of the cardioid.</summary>
    public const double Threshold = 0.25;

    /// <summary>1/2 magnitude pinned to the cardioid via HalfAsStructuralFixedPoint.</summary>
    public HalfAsStructuralFixedPointClaim Half { get; }

    /// <summary>1/4 = b² magnitude of c at the real-axis tangent point only.</summary>
    public QuarterAsBilinearMaxvalClaim Quarter { get; }

    /// <summary>The 90° rotation that lifts c from the real axis to the
    /// complex plane (generator of the e^(iφ) parametrization).</summary>
    public NinetyDegreeMirrorMemoryClaim NinetyDegree { get; }

    /// <summary>Polynomial-foundation typed parent: the c = 0 case where the
    /// fixed point is z* = 0 (degenerate).</summary>
    public PolynomialFoundationClaim Polynomial { get; }

    /// <summary>Compute c(φ) = b·e^(iφ) − b²·e^(2iφ) on the cardioid at angle φ.</summary>
    public Complex CardioidC(double phi)
    {
        Complex e_iphi = Complex.Exp(Complex.ImaginaryOne * phi);
        return B * e_iphi - Threshold * (e_iphi * e_iphi);
    }

    /// <summary>Compute the period-1 fixed point z*(φ) = b·e^(iφ) at the
    /// cardioid parameter φ. By construction |z*| = b and arg(z*) = φ.</summary>
    public Complex CardioidFixedPoint(double phi)
    {
        Complex e_iphi = Complex.Exp(Complex.ImaginaryOne * phi);
        return B * e_iphi;
    }

    /// <summary>Magnitude invariance: |z*(φ)| should equal b = 1/2 for ANY φ.
    /// Returns the actual computed magnitude (should be bit-exact b).
    /// This is the "argmax side" of yesterday's argmax/maxval pair: the
    /// magnitude at which the structural anchor lives.</summary>
    public double FixedPointMagnitude(double phi) =>
        Complex.Abs(CardioidFixedPoint(phi));

    /// <summary>Squared-magnitude invariance: |z*(φ)|² should equal b² = 1/4
    /// for ANY φ. This is the "maxval side" of the argmax/maxval pair: the
    /// projection value under squaring. Both <see cref="FixedPointMagnitude"/>
    /// (= b = 1/2) and this method (= b² = 1/4) are invariants of the
    /// cardioid; the Half and Quarter anchors are two metric-power readings
    /// of the same fixed-point quantity.</summary>
    public double FixedPointMagnitudeSquared(double phi)
    {
        double mag = FixedPointMagnitude(phi);
        return mag * mag;
    }

    /// <summary>Phase parameter: arg(z*(φ)) should equal φ (mod 2π).</summary>
    public double FixedPointArgument(double phi)
    {
        Complex z = CardioidFixedPoint(phi);
        return Math.Atan2(z.Imaginary, z.Real);
    }

    /// <summary>Drift indicator: |z*(φ)| = b for all φ. Tests at canonical
    /// points φ ∈ {0, π/3, π/2, π}.</summary>
    public bool MagnitudeInvariantAroundCardioid()
    {
        double[] phis = { 0.0, Math.PI / 3, Math.PI / 2, Math.PI, 4 * Math.PI / 3, 5 * Math.PI / 3 };
        foreach (double phi in phis)
        {
            if (Math.Abs(FixedPointMagnitude(phi) - B) > 1e-14) return false;
        }
        return true;
    }

    /// <summary>Drift indicator: |z*(φ)|² = b² = 1/4 for all φ. Tests at
    /// canonical points; this is the Quarter anchor's invariance reading,
    /// complementary to <see cref="MagnitudeInvariantAroundCardioid"/>.</summary>
    public bool SquaredMagnitudeInvariantAroundCardioid()
    {
        double[] phis = { 0.0, Math.PI / 3, Math.PI / 2, Math.PI, 4 * Math.PI / 3, 5 * Math.PI / 3 };
        foreach (double phi in phis)
        {
            if (Math.Abs(FixedPointMagnitudeSquared(phi) - Threshold) > 1e-14) return false;
        }
        return true;
    }

    /// <summary>Drift indicator: c(φ) = z*(1 − z*) algebraic identity.</summary>
    public bool AlgebraicIdentityHolds(double phi)
    {
        Complex c = CardioidC(phi);
        Complex z = CardioidFixedPoint(phi);
        Complex c_from_identity = z * (Complex.One - z);
        return Complex.Abs(c - c_from_identity) < 1e-14;
    }

    /// <summary>Drift indicator: at φ = 0, c = 1/4 = Threshold exactly
    /// (real-axis cusp recovery, where F95 also lives).</summary>
    public bool CuspAgreesWithF95Threshold()
    {
        Complex c_cusp = CardioidC(0.0);
        return Math.Abs(c_cusp.Real - Threshold) < 1e-15
               && Math.Abs(c_cusp.Imaginary) < 1e-15;
    }

    /// <summary>Drift indicator: at φ = π, c = −3/4 (the real-axis "tail",
    /// boundary of the period-doubling region).</summary>
    public bool TailAtMinusThreeQuarters()
    {
        Complex c_tail = CardioidC(Math.PI);
        return Math.Abs(c_tail.Real - (-0.75)) < 1e-14
               && Math.Abs(c_tail.Imaginary) < 1e-14;
    }

    public F97CardioidHalfFixedPointPi2Inheritance(
        HalfAsStructuralFixedPointClaim half,
        QuarterAsBilinearMaxvalClaim quarter,
        NinetyDegreeMirrorMemoryClaim ninetyDegree,
        PolynomialFoundationClaim polynomial)
        : base("F97 Mandelbrot cardioid parametrization at b = 1/2: " +
               "c(φ) = b·e^(iφ) − b²·e^(2iφ), z*(φ) = b·e^(iφ); " +
               "|z*| = b = 1/2 invariant around the curve; arg(z*) = φ cardioid parameter; " +
               "bit-exact algebraic identity, machine-precision numerical verification",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F97_CARDIOID_HALF_FIXED_POINT.md + " +
               "docs/ANALYTICAL_FORMULAS.md F97 + " +
               "simulations/_cardioid_parametrization_tier1.py + " +
               "experiments/CPSI_COMPLEX_PLANE.md (Kingston 2026-04-16 hardware precursor) + " +
               "experiments/BOUNDARY_NAVIGATION.md (real-c precursor, F95 origin) + " +
               "compute/RCPsiSquared.Core/Symmetry/F95AngleAtQuadraticZeroPi2Inheritance.cs (real-c angle sibling) + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (HalfAsStructuralFixedPoint, QuarterAsBilinearMaxval, NinetyDegreeMirrorMemory, PolynomialFoundation typed parents)")
    {
        Half = half ?? throw new ArgumentNullException(nameof(half));
        Quarter = quarter ?? throw new ArgumentNullException(nameof(quarter));
        NinetyDegree = ninetyDegree ?? throw new ArgumentNullException(nameof(ninetyDegree));
        Polynomial = polynomial ?? throw new ArgumentNullException(nameof(polynomial));
    }

    public override string DisplayName =>
        "F97 Mandelbrot cardioid parametrization at b = 1/2: |z*| = b invariant around the curve";

    public override string Summary =>
        $"c(φ) = b·e^(iφ) − b²·e^(2iφ) at b = {B}; z*(φ) = b·e^(iφ); " +
        $"|z*| = {B} invariant for all φ ∈ [0, 2π); cusp at φ = 0 gives c = {Threshold} (F95 limit); " +
        $"tail at φ = π gives c = −3/4 (period-doubling boundary); " +
        $"bit-exact algebraic identity ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("B (= 1/2 = HalfAsStructuralFixedPoint)", B);
            yield return InspectableNode.RealScalar("Threshold = B² (= 1/4, magnitude of c at cusp only)", Threshold);
            double[] sampleAngles = { 0.0, Math.PI / 3, Math.PI / 2, Math.PI };
            string[] sampleLabels = { "φ=0 (cusp)", "φ=π/3 (60°)", "φ=π/2 (top)", "φ=π (tail)" };
            for (int i = 0; i < sampleAngles.Length; i++)
            {
                Complex c = CardioidC(sampleAngles[i]);
                Complex z = CardioidFixedPoint(sampleAngles[i]);
                yield return new InspectableNode($"Cardioid sample: {sampleLabels[i]}",
                    summary: $"c = {c.Real:G4} + {c.Imaginary:G4}i, |c| = {Complex.Abs(c):G4}, " +
                             $"z* = {z.Real:G4} + {z.Imaginary:G4}i, |z*| = {Complex.Abs(z):G4} (= b)");
            }
            yield return new InspectableNode("Drift checks",
                summary: $"MagnitudeInvariantAroundCardioid (|z*| = b = 1/2): {MagnitudeInvariantAroundCardioid()}; " +
                         $"SquaredMagnitudeInvariantAroundCardioid (|z*|² = b² = 1/4): {SquaredMagnitudeInvariantAroundCardioid()}; " +
                         $"AlgebraicIdentityHolds (at π/3): {AlgebraicIdentityHolds(Math.PI / 3)}; " +
                         $"CuspAgreesWithF95Threshold: {CuspAgreesWithF95Threshold()}; " +
                         $"TailAtMinusThreeQuarters: {TailAtMinusThreeQuarters()}");
            yield return new InspectableNode("Structural reading (both anchors invariant)",
                summary: "The Mandelbrot cardioid carries BOTH typed anchors as invariants of the same z*, at two metric powers: " +
                         "|z*| = b = 1/2 (HalfAsStructuralFixedPoint, argmax side) and " +
                         "|z*|² = b² = 1/4 (QuarterAsBilinearMaxval, maxval side). " +
                         "This is exactly the argmax/maxval pair of ON_HOW_TWO_SIDES_MEET_AT_THE_QUARTER " +
                         "(2026-05-16), now geometric on the cardioid. 1/2 = 2·(1/4) sits on the dyadic ladder; " +
                         "the polarity sides ±1/2 each square to 1/4. " +
                         "By contrast |c|² = 5/16 − (1/4)·cos(φ) is NOT invariant; |c| = 1/4 only at the cusp.");
            yield return new InspectableNode("Sibling at real-axis cusp",
                summary: "F95AngleAtQuadraticZeroPi2Inheritance: θ(c; b) = arctan(√(c/b² − 1)) for real c > b². " +
                         "F95 covers the angle of the complex root pair past the cusp; F97 covers the full complex-c " +
                         "extension via the cardioid parametrization. F95 and F97 cover complementary regions of the same z² − 2bz + c algebra.");
        }
    }
}
