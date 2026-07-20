using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F33, the N=3 rate ladder (ANALYTICAL_FORMULAS F33):
///
/// <code>
///   Exact at every J (pure-weight rungs, multiplicities out of 64):
///     rate = 0·γ    ⟨n_XY⟩ = 0   mult 4    (kernel)
///     rate = 2·γ    ⟨n_XY⟩ = 1   mult 14   (pure weight-1)
///     rate = 4·γ    ⟨n_XY⟩ = 2   mult 14   (pure weight-2)
///     rate = 6·γ    ⟨n_XY⟩ = 3   mult 4    (XOR sector)
///
///   A J/γ → ∞ limit only (each is a TRIPLE at finite J):
///     rate_2 →  8·γ / 3      ⟨n_XY⟩ → 4/3   (cross-sector mix)
///     rate_3 → 10·γ / 3      ⟨n_XY⟩ → 5/3   (cross-sector mix)
///
///   Absorption Theorem: α = 2·γ · ⟨n_XY⟩, exact in both cases.
/// </code>
///
/// <para>F33 is the N=3 specialization of the Liouvillian decay-rate spectrum
/// for Heisenberg + Z-dephasing. The pure-weight rungs are J-independent: a
/// mode sitting on one weight sector has nothing left to mix. The two
/// fractional rates are NOT exact rationals at finite coupling. Each is a
/// triple of distinct levels whose lowest member approaches the limit as
/// 0.46·(γ/J)² while the band's own spread closes as 0.53·(γ/J)²; in the
/// canonical regime Q = 1.5 the "8γ/3" band is 2.4607, 2.6040, 2.6980 γ, up to
/// 7.7% away from 8/3. <see cref="Rate2"/> and <see cref="Rate3"/> return the
/// limit values and must be cited as limits.</para>
///
/// <para><b>Why N=3 is special:</b> the pure-weight rungs 2γ (weight-1, F50)
/// and 2(N−1)γ (the F3 band top; F43 is the XOR sector at 2Nγ, one rung above)
/// remain universal at all N. Topology-dependence
/// of the internal spectrum starts already at N=3, not at N ≥ 4: on the triangle
/// the 2.4607 level is absent and the pure-rung multiplicities are 4,16,16,4
/// against the chain's 4,14,14,4.</para>
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
/// Clebsch-Gordan structure when weight sectors combine through the Heisenberg
/// J·SWAP interaction.</para>
///
/// <para>Tier1Derived: the pure-weight rungs are Tier 1 (exact at every J via
/// N=3 diagonalization). The two fractional rates are Tier 1 as a J/γ → ∞
/// limit, not as closed-form values at finite coupling. Valid for the N=3
/// Heisenberg chain + Z-dephasing. Pi2-Foundation anchoring is composition
/// through F50 + Pi2DyadicLadder.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F33 +
/// <c>experiments/SIGNAL_PROCESSING_VIEW.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F50WeightOneDegeneracyPi2Inheritance.cs</c>
/// (rate_1 = 2γ universal anchor) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F33ExactN3DecayRatesPi2Inheritance : Claim, IZ2AxisClaim
{

    /// <summary>The F1² / Π²_Z axis (bit_b parity, n_Y + n_Z mod 2). The
    /// canonical Pi²-Inheritance axis. The bit_a-twin (Π²_X / F61 axis) is
    /// currently not typed for this Claim.</summary>
    public Z2Axis Z2Axis => Z2Axis.BitB;

    /// <summary>The typed bit_a-twin sibling, if one exists. Currently null
    /// (no bit_a twin is typed for this Claim; this is an open slot in the
    /// cubic-architecture coverage).</summary>
    public Claim? BitATwin => null;

    // Absorption-Theorem descendant; bit_a twin holds by the Hadamard X↔Z duality (PROOF_BIT_A_TWIN_VIA_HADAMARD.md).
    public BitATwinClassification BitATwinStatus => BitATwinClassification.CoveredByHadamardDuality;
    public Pi2DyadicLadderClaim Ladder { get; }
    public F50WeightOneDegeneracyPi2Inheritance F50 { get; }

    /// <summary>The "2" coefficient in rate_1 = 2γ. Live from Pi2DyadicLadder a_0.
    /// Identical to F50's DecayRateFactor at N=3.</summary>
    public double WeightOneRateCoefficient => Ladder.Term(0);

    /// <summary>F33's rate_1 = 2γ at N=3. Universal weight-1 eigenvalue position
    /// per F50; here specialized to N=3.</summary>
    public double Rate1(double gammaZero)
    {
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        return WeightOneRateCoefficient * gammaZero;
    }

    /// <summary>F33's rate_2 → 8γ/3 at N=3, the J/γ → ∞ limit of the first mixed
    /// cross-sector supermode band. At finite J this is a triple, not one level.</summary>
    public double Rate2(double gammaZero)
    {
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        return 8.0 * gammaZero / 3.0;
    }

    /// <summary>F33's rate_3 → 10γ/3 at N=3, the J/γ → ∞ limit of the second mixed
    /// cross-sector supermode band. At finite J this is a triple, not one level.</summary>
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

    /// <summary>The three F33 rates' ⟨n_XY⟩ values: rate_1 ↔ 1 (exact at every J),
    /// rate_2 ↔ 4/3 and rate_3 ↔ 5/3 (limit values). Verifies the Absorption
    /// Theorem, which holds exactly in both cases.</summary>
    public IReadOnlyList<(int RateIndex, double NXyExpectation)> AbsorptionTheoremTable { get; } = new[]
    {
        (RateIndex: 1, NXyExpectation: 1.0),
        (RateIndex: 2, NXyExpectation: 4.0 / 3.0),
        (RateIndex: 3, NXyExpectation: 5.0 / 3.0),
    };

    /// <summary>Drift check on the limit values: rate_2 / rate_3 = 8/10 = 4/5.</summary>
    public bool RatesAreRationalRatios(double gammaZero, double tolerance = 1e-12)
    {
        double r2 = Rate2(gammaZero);
        double r3 = Rate3(gammaZero);
        return Math.Abs(r2 / r3 - 0.8) < tolerance;
    }

    public F33ExactN3DecayRatesPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        F50WeightOneDegeneracyPi2Inheritance f50)
        : base("F33 N=3 rate ladder: pure-weight rungs 0, 2γ (= F50 specialization), 4γ, 6γ exact at every J; rate_2 → 8γ/3 and rate_3 → 10γ/3 as J/γ → ∞ only; Absorption Theorem α = 2γ·⟨n_XY⟩",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F33 + " +
               "experiments/SIGNAL_PROCESSING_VIEW.md + " +
               "compute/RCPsiSquared.Core/Symmetry/F50WeightOneDegeneracyPi2Inheritance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        Ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        F50 = f50 ?? throw new ArgumentNullException(nameof(f50));
    }

    public override string DisplayName =>
        "F33 N=3 rate ladder as Pi2-Foundation a_0 + F50 (specialization at N=3) inheritance";

    public override string Summary =>
        $"N=3 Heisenberg+Z-deph: rungs 0, 2γ (= F50), 4γ, 6γ exact at every J; rate_2 → 8γ/3, rate_3 → 10γ/3 as J/γ → ∞; Absorption α = 2γ·⟨n_XY⟩ ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F33 closed form",
                summary: "exact at every J: 0 (kernel, mult 4), 2γ (w=1 pure, mult 14), 4γ (w=2 pure, mult 14), 6γ (XOR, mult 4); limit only: rate_2 → 8γ/3 and rate_3 → 10γ/3 (cross-sector mix, parity-preserving: w=0/w=2 and w=1/w=3, never w=1/w=2), each a triple at finite J whose lowest member approaches as 0.46·(γ/J)² and whose spread closes as 0.53·(γ/J)²");
            yield return InspectableNode.RealScalar("WeightOneRateCoefficient (= a_0 = 2)", WeightOneRateCoefficient);
            yield return new InspectableNode("F50 specialization at N=3",
                summary: $"F33's rate_1 = 2γ IS F50's universal weight-1 eigenvalue position. F50.DecayRateFactor (= {F50.DecayRateFactor}) is the same '2' as F33's WeightOneRateCoefficient.");
            yield return new InspectableNode("Absorption Theorem table",
                summary: "rate_1 ↔ ⟨n_XY⟩ = 1 (pure weight-1, exact at every J); rate_2 ↔ ⟨n_XY⟩ → 4/3 and rate_3 ↔ ⟨n_XY⟩ → 5/3 (mixes, limit values, since the mixing weights themselves depend on J); the theorem α = 2γ·⟨n_XY⟩ holds exactly in both cases");
            yield return new InspectableNode("N=3 verified",
                summary: $"rate_1(γ=0.05) = {Rate1(0.05):G6}; rate_2(γ=0.05) = {Rate2(0.05):G6}; rate_3(γ=0.05) = {Rate3(0.05):G6}; ratio 2/3 = {Rate2(1.0)/Rate3(1.0):F4} (= 4/5 = 0.8 exactly)");
            yield return new InspectableNode("two information channels",
                summary: "at N=3 frequency (from H) and decay (from dissipator) are perfectly orthogonal channels; each rate factors into oscillation × decay cleanly");
            yield return new InspectableNode("universality breaks",
                summary: "internal rates are topology-dependent already at N=3 (triangle rungs 4,16,16,4 vs chain 4,14,14,4); only the rung positions 2γ (weight-1, F50) and 2(N−1)γ (F3 band top) remain universal");
        }
    }
}
