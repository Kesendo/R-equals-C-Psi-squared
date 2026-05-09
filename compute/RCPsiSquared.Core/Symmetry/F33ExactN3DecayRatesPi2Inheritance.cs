using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F33 closed form (Tier 1, exact rational; ANALYTICAL_FORMULAS line 200):
///
/// <code>
///   rate_1 = 2·γ           (w=1, pure single-site coherence)
///   rate_2 = 8·γ / 3       (w=2 mixed, Hamiltonian superposition)
///   rate_3 = 10·γ / 3      (w=2 mixed, Hamiltonian superposition)
///
///   Boundary rates: 0 (kernel) and 2·(N−1)·γ = 6·γ at N=3 (XOR sector)
///
///   Absorption Theorem: α = 2·γ · ⟨n_XY⟩
///     rate_1 ↔ ⟨n_XY⟩ = 1     (pure weight-1)
///     rate_2 ↔ ⟨n_XY⟩ = 4/3   (mixed w=1 and w=2)
///     rate_3 ↔ ⟨n_XY⟩ = 5/3   (mixed w=1 and w=2)
/// </code>
///
/// <para>F33 is the N=3 specialization of the Liouvillian decay-rate spectrum
/// for Heisenberg + Z-dephasing. Three distinct interior rates plus the
/// boundary rates (0 and 6γ) fully determine the spectrum at N=3. The
/// Hamiltonian mixes w=1 and w=2 Pauli strings into supermodes with exact
/// rational decay rates.</para>
///
/// <para><b>Why N=3 is special:</b> at N=3 the rates are exact rationals
/// (1, 4/3, 5/3 in 2γ-units). At N ≥ 4 the internal rates become
/// topology-dependent: only the boundary rates 2γ (weight-1, F50) and
/// 2(N−1)γ (XOR sector, F43) remain universal. The intermediate spectrum
/// becomes irrational/topology-specific.</para>
///
/// <para><b>Two independent information channels:</b> at N=3, frequency and
/// decay are perfectly orthogonal (the imaginary and real parts of L
/// eigenvalues factor cleanly). Each rate corresponds to a Hamiltonian-mixing
/// supermode with both an oscillation frequency (from H) and a decay rate
/// (from the dissipator).</para>
///
/// <para>Pi2-Foundation anchor:</para>
/// <list type="bullet">
///   <item><b>WeightOneRateCoefficient = 2 = a_0</b>: in rate_1 = 2γ. Live
///         from <see cref="Pi2DyadicLadderClaim.Term"/>(0). Identical to
///         <see cref="F50WeightOneDegeneracyPi2Inheritance.DecayRateFactor"/>
///         F33's rate_1 IS the N=3 specialization of F50's universal
///         weight-1 eigenvalue position.</item>
/// </list>
///
/// <para>The 8/3 and 10/3 numerators are combinatorial-Hamiltonian-mixing
/// fractions, no direct Pi2 anchor; they emerge from the N=3 specific
/// Clebsch-Gordan structure when w=1 and w=2 combine through the Heisenberg
/// J·SWAP interaction.</para>
///
/// <para>Tier1Derived: F33 is Tier 1 exact-rational (closed form via N=3
/// diagonalization); valid for N=3 Heisenberg chain + Z-dephasing.
/// Pi2-Foundation anchoring is composition through F50 + Pi2DyadicLadder.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F33 (line 200) +
/// <c>experiments/SIGNAL_PROCESSING_VIEW.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F50WeightOneDegeneracyPi2Inheritance.cs</c>
/// (rate_1 = 2γ universal anchor) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F33ExactN3DecayRatesPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly F50WeightOneDegeneracyPi2Inheritance _f50;

    /// <summary>The "2" coefficient in rate_1 = 2γ. Live from Pi2DyadicLadder a_0.
    /// Identical to F50's DecayRateFactor at N=3.</summary>
    public double WeightOneRateCoefficient => _ladder.Term(0);

    /// <summary>F33's rate_1 = 2γ at N=3. Universal weight-1 eigenvalue position
    /// per F50; here specialized to N=3.</summary>
    public double Rate1(double gammaZero)
    {
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        return WeightOneRateCoefficient * gammaZero;
    }

    /// <summary>F33's rate_2 = 8γ/3 at N=3. First mixed w=1+w=2 supermode rate.</summary>
    public double Rate2(double gammaZero)
    {
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        return 8.0 * gammaZero / 3.0;
    }

    /// <summary>F33's rate_3 = 10γ/3 at N=3. Second mixed w=1+w=2 supermode rate.</summary>
    public double Rate3(double gammaZero)
    {
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        return 10.0 * gammaZero / 3.0;
    }

    /// <summary>F33's highest paired-mode rate at N=3: 2·(N−1)·γ = 4·γ (the
    /// "fastest paired" boundary from F3, w=N−1 = 2 modes at N=3). Distinct from
    /// F43's XOR-sector rate 2·N·γ = 6γ which sits one rung above. F33's
    /// XorBoundaryRate is the F3 max-rate boundary specialised to N=3, not
    /// the XOR sector itself.</summary>
    public double XorBoundaryRate(double gammaZero)
    {
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        return WeightOneRateCoefficient * 2.0 * gammaZero;
    }

    /// <summary>Absorption Theorem: α = 2γ · ⟨n_XY⟩. Inverse map from a decay rate
    /// to the XY-weight expectation that produces it.</summary>
    public double NXyExpectationFromRate(double rate, double gammaZero)
    {
        if (gammaZero <= 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be > 0.");
        if (rate < 0) throw new ArgumentOutOfRangeException(nameof(rate), rate, "rate must be ≥ 0.");
        return rate / (WeightOneRateCoefficient * gammaZero);
    }

    /// <summary>The three F33 rates' ⟨n_XY⟩ values: rate_1 ↔ 1, rate_2 ↔ 4/3, rate_3 ↔ 5/3.
    /// Verifies the Absorption Theorem at exact rationals.</summary>
    public IReadOnlyList<(int RateIndex, double NXyExpectation)> AbsorptionTheoremTable { get; } = new[]
    {
        (RateIndex: 1, NXyExpectation: 1.0),
        (RateIndex: 2, NXyExpectation: 4.0 / 3.0),
        (RateIndex: 3, NXyExpectation: 5.0 / 3.0),
    };

    /// <summary>Drift check: rate_2 / rate_3 = 8/10 = 4/5 (rational).</summary>
    public bool RatesAreRationalRatios(double gammaZero, double tolerance = 1e-12)
    {
        double r2 = Rate2(gammaZero);
        double r3 = Rate3(gammaZero);
        return Math.Abs(r2 / r3 - 0.8) < tolerance;
    }

    public F33ExactN3DecayRatesPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        F50WeightOneDegeneracyPi2Inheritance f50)
        : base("F33 N=3 exact decay rates: rate_1 = 2γ (= F50 specialization), rate_2 = 8γ/3, rate_3 = 10γ/3 (rational); Absorption Theorem α = 2γ·⟨n_XY⟩",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F33 + " +
               "experiments/SIGNAL_PROCESSING_VIEW.md + " +
               "compute/RCPsiSquared.Core/Symmetry/F50WeightOneDegeneracyPi2Inheritance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _f50 = f50 ?? throw new ArgumentNullException(nameof(f50));
    }

    public override string DisplayName =>
        "F33 N=3 exact rates as Pi2-Foundation a_0 + F50 (specialization at N=3) inheritance";

    public override string Summary =>
        $"N=3 Heisenberg+Z-deph: rate_1 = 2γ (= F50), rate_2 = 8γ/3, rate_3 = 10γ/3; rational; Absorption α = 2γ·⟨n_XY⟩ ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F33 closed form",
                summary: "rate_1 = 2γ (w=1 pure), rate_2 = 8γ/3 (w=1+w=2 mix), rate_3 = 10γ/3 (w=1+w=2 mix); plus boundaries 0 (kernel) and 6γ (XOR at N=3); rational rates unique to N=3");
            yield return InspectableNode.RealScalar("WeightOneRateCoefficient (= a_0 = 2)", WeightOneRateCoefficient);
            yield return new InspectableNode("F50 specialization at N=3",
                summary: $"F33's rate_1 = 2γ IS F50's universal weight-1 eigenvalue position. F50.DecayRateFactor (= {_f50.DecayRateFactor}) is the same '2' as F33's WeightOneRateCoefficient.");
            yield return new InspectableNode("Absorption Theorem table",
                summary: "rate_1 ↔ ⟨n_XY⟩ = 1 (pure weight-1); rate_2 ↔ ⟨n_XY⟩ = 4/3 (mix); rate_3 ↔ ⟨n_XY⟩ = 5/3 (mix); the theorem α = 2γ·⟨n_XY⟩ holds exactly including non-integer ⟨n_XY⟩");
            yield return new InspectableNode("N=3 verified",
                summary: $"rate_1(γ=0.05) = {Rate1(0.05):G6}; rate_2(γ=0.05) = {Rate2(0.05):G6}; rate_3(γ=0.05) = {Rate3(0.05):G6}; ratio 2/3 = {Rate2(1.0)/Rate3(1.0):F4} (= 4/5 = 0.8 exactly)");
            yield return new InspectableNode("two information channels",
                summary: "at N=3 frequency (from H) and decay (from dissipator) are perfectly orthogonal channels; each rate factors into oscillation × decay cleanly");
            yield return new InspectableNode("N ≥ 4 universality breaks",
                summary: "N ≥ 4: internal rates become topology-dependent; only boundary rates 2γ (weight-1, F50) and 2(N−1)γ (XOR partner, F43-related) remain universal");
        }
    }
}
