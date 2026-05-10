using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F89 topology orbit closure (Tier 1 derived, verified bit-identical at
/// N=7 across all 14 topology classes and N=4 across all C(4,2)=6 site pairs).
///
/// <code>
///   For ρ_cc = (|S_1⟩⟨S_2| + h.c.) / 2, uniform-J multi-bond XY
///   H_B = J · Σ_{(p,q) ∈ B} (X_p X_q + Y_p Y_q), and uniform Z-dephasing γ₀:
///
///     S(t) = Σ_l 2 · |(ρ_l(t))_{0,1}|²
///
///   depends only on the S_N-orbit of B under σ · B = {(σ(p), σ(q)) : (p,q) ∈ B}.
///
///   For chain-restricted B (NN-bonds only), the orbit equals the bond-graph
///   topology class: the sorted multiset of connected-path-lengths.
///
///   S(0) = (N − 1) / N,  closed form (probe-only, independent of B).
/// </code>
///
/// <para>Proof: S_N-orbit transitivity. ρ_cc, the uniform Z-dissipator, and the
/// spatial-sum kernel are all S_N-symmetric; H_B transforms as
/// H_{σ·B} = U_σ H_B U_σ^†; the Lindblad solution is therefore covariant under
/// conjugation by U_σ, and the S_N-symmetric kernel projects out everything
/// except the orbit label.</para>
///
/// <para>Cited scaffolding (registered as ancestor edges in the runtime via
/// <c>F89TopologyOrbitClosureRegistration</c>):
/// <see cref="F73SpatialSumPurityClosurePi2Inheritance"/> is the closely related
/// closure for the (vac, SE) coherence block (same orbit-style argument; F73
/// gets a closed exponential thanks to a uniform per-element 2γ₀ rate, while
/// F89's (S_1, S_2) block has non-uniform per-element decay 2γ₀ on overlap and
/// 6γ₀ off, hence orbit-only with no closed exponential).
/// <see cref="F71MirrorSymmetryPi2Inheritance"/> is the spatial-Z₂ subgroup
/// b ↔ N−2−b sitting inside the full S_N argument used here.</para>
///
/// <para>All-isolated subclass closed form (Tier 1 derived sub-claim, see
/// <see cref="AllIsolatedClosedForm"/>): for the (1)^m topology classes
/// (m disjoint NN-bonds), S(t) = [(N−1)/N + 4m(N−2)(cos(4Jt)−1)/(N²(N−1))]·exp(−4γ₀ t),
/// with universal asymptotic rate 4γ₀ matching F73's vac-SE rate. Mixed-topology
/// and pure-path closed forms remain Tier 2 empirical (derivation open).</para>
///
/// <para>Pi2-Foundation inheritance (two roles, same a_{−1}): both the time-decay
/// coefficient 4γ₀ and the time-oscillation coefficient 4J in the all-isolated
/// closed form trace to <see cref="Pi2DyadicLadderClaim.Term"/>(−1) = 4. The
/// γ-axis inheritance is identical to F73's
/// <see cref="F73SpatialSumPurityClosurePi2Inheritance.DecayRateCoefficient"/>
/// (per-coherence rate 2γ₀ doubles to S-decay rate 4γ₀ on |·|²). The J-axis
/// inheritance is the same a_{0}=2 → a_{−1}=4 doubling: H_B-eigenstate
/// frequency 2J doubles to S-oscillation frequency 4J on |·|². The Pi2 ladder
/// thus anchors the time coefficients on both energy axes; the (N−1)/N baseline
/// and the 4m(N−2)/(N²(N−1)) correction prefactor are combinatorial (S_N-orbit
/// and 2-qubit-block algebra), not Pi2-anchored.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F89 (line 2436) +
/// <c>experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md</c> +
/// <c>simulations/_bond_isolate_compare_n7.py</c> +
/// <c>simulations/_bond_isolate_long_range_verify.py</c> +
/// <c>simulations/_bond_isolate_topology_classes_n7.py</c> +
/// <c>compute/RCPsiSquared.Propagate</c> <c>bond-isolate --bonds i,j,...</c>
/// CLI mode.</para></summary>
public sealed class F89TopologyOrbitClosure : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;

    /// <summary>The "4" decay-rate coefficient in <c>exp(−4γ₀ t)</c>, the asymptotic
    /// rate of S(t) for any all-isolated topology class. Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(−1) = a_{−1} = 4. Identical anchor
    /// to F73's <see cref="F73SpatialSumPurityClosurePi2Inheritance.DecayRateCoefficient"/>:
    /// the per-coherence Z-deph rate 2γ₀ (= a_{0}·γ₀) doubles to S-decay rate 4γ₀
    /// (= a_{−1}·γ₀) on |·|².</summary>
    public double DecayRateCoefficient => _ladder.Term(-1);

    /// <summary>The "4" oscillation-frequency coefficient in <c>cos(4J t)</c>, the
    /// frequency of the m-correction in S(t) for all-isolated topology classes.
    /// Live from <see cref="Pi2DyadicLadderClaim.Term"/>(−1) = a_{−1} = 4, the
    /// J-axis mirror of the same Pi2 ladder anchor that gives F73's decay rate:
    /// the H_B-eigenstate frequency 2J (= a_{0}·J) doubles to S-oscillation
    /// frequency 4J (= a_{−1}·J) on |·|².</summary>
    public double OscillationFrequencyCoefficient => _ladder.Term(-1);

    /// <summary>Live drift check: both the γ-axis decay coefficient and the J-axis
    /// oscillation coefficient resolve to <c>a_{−1} = 4</c> from the Pi2 ladder, and
    /// match the literal <c>4.0</c> used in <see cref="AllIsolatedClosedForm"/>.</summary>
    public bool Pi2DoublingConsistent() =>
        Math.Abs(DecayRateCoefficient - 4.0) < 1e-15
        && Math.Abs(OscillationFrequencyCoefficient - 4.0) < 1e-15;

    /// <summary>Probe-only closed form: S(0) = (N − 1) / N for ρ_cc on N qubits.
    /// Bond-set independent (depends only on the probe).</summary>
    public static double S0ClosedForm(int n)
    {
        if (n < 2)
            throw new ArgumentOutOfRangeException(nameof(n), n, "F89 requires N ≥ 2.");
        return (double)(n - 1) / n;
    }

    /// <summary>The canonical chain-topology-class label for a bond set on a chain:
    /// the sorted multiset of connected-path-lengths in the bond-graph induced by
    /// B ⊂ {NN-bonds}. Two bond sets are in the same S_N-orbit iff their topology
    /// classes are equal as sorted arrays.</summary>
    public static int[] ChainTopologyClass(int n, IReadOnlyList<int> bonds)
    {
        if (n < 2)
            throw new ArgumentOutOfRangeException(nameof(n), n, "F89 requires N ≥ 2.");
        if (bonds is null) throw new ArgumentNullException(nameof(bonds));
        if (bonds.Count == 0) return Array.Empty<int>();

        var seen = new HashSet<int>(bonds.Count);
        foreach (int b in bonds)
        {
            if (b < 0 || b > n - 2)
                throw new ArgumentOutOfRangeException(nameof(bonds), b,
                    $"Bond index {b} out of range [0, {n - 2}].");
            if (!seen.Add(b))
                throw new ArgumentException("Bond indices must be distinct.", nameof(bonds));
        }

        var sorted = bonds.OrderBy(b => b).ToList();
        var components = new List<int>();
        int currentLen = 1;
        for (int i = 1; i < sorted.Count; i++)
        {
            if (sorted[i] == sorted[i - 1] + 1)
                currentLen++;
            else
            {
                components.Add(currentLen);
                currentLen = 1;
            }
        }
        components.Add(currentLen);
        components.Sort();
        return components.ToArray();
    }

    /// <summary>True iff two bond sets share the same chain-topology class
    /// (S_N-orbit equivalent). This is the predicate that S(t; bondsA) =
    /// S(t; bondsB) per F89 for the ρ_cc + uniform-J multi-bond XY setup.</summary>
    public static bool AreInSameChainOrbit(int n,
                                            IReadOnlyList<int> bondsA,
                                            IReadOnlyList<int> bondsB)
    {
        var clsA = ChainTopologyClass(n, bondsA);
        var clsB = ChainTopologyClass(n, bondsB);
        return clsA.SequenceEqual(clsB);
    }

    /// <summary>Closed form for the all-isolated topology subclass (1)^m on N qubits
    /// with m disjoint NN-bonds, N − 2m bare sites, uniform XY coupling J, uniform
    /// Z-dephasing γ₀:
    ///
    /// <code>
    ///   S_(1)^m, N(t) = [(N − 1)/N + 4m(N − 2)(cos(4Jt) − 1)/(N²(N − 1))] · exp(−4γ₀ t)
    /// </code>
    ///
    /// <para>Asymptotic rate 4γ₀ universal across m (matches F73's vac-SE rate). The
    /// m-correction has frequency 4J and vanishes at cos(4Jt) = 1 (period π/(2J)),
    /// the in-phase moment where all (1)^m classes collapse to S(0)·exp(−4γ₀ t).</para>
    ///
    /// <para>Derivation: Lindbladian factorises over disjoint 2-qubit blocks plus
    /// bare sites; in each block only (vac, SE)_B and (SE, DE)_B-overlap coherences
    /// are populated (both rate 2γ₀ per coherence; 4γ₀ on |·|²). The 6γ₀ no-overlap
    /// channel requires a SE site outside the doubly-excited pair, impossible on a
    /// 2-qubit block alone. H_B-eigenstate phase tracking yields the cos(4Jt)
    /// interference. Verified bit-exact against bond-isolate CSVs N=7 m=1,2,3 in
    /// <c>simulations/_f89_all_isolated_closed_form_verify.py</c>.</para></summary>
    public static double AllIsolatedClosedForm(int n, int m, double j, double gammaZero, double t)
    {
        if (n < 2)
            throw new ArgumentOutOfRangeException(nameof(n), n, "F89 requires N ≥ 2.");
        if (m < 0)
            throw new ArgumentOutOfRangeException(nameof(m), m, "Bond count m must be ≥ 0.");
        if (2 * m > n)
            throw new ArgumentOutOfRangeException(nameof(m), m,
                $"All-isolated topology (1)^m needs 2m ≤ N sites; got m={m}, N={n}.");
        if (gammaZero < 0.0)
            throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        if (t < 0.0)
            throw new ArgumentOutOfRangeException(nameof(t), t, "t must be ≥ 0.");

        double s0 = (double)(n - 1) / n;
        double correction = m == 0
            ? 0.0
            : 4.0 * m * (n - 2) * (Math.Cos(4.0 * j * t) - 1.0) / ((double)n * n * (n - 1));
        return (s0 + correction) * Math.Exp(-4.0 * gammaZero * t);
    }

    public F89TopologyOrbitClosure(Pi2DyadicLadderClaim ladder)
        : base("F89 topology orbit closure: S(t) for ρ_cc + uniform-J multi-bond XY depends only on the S_N-orbit of the bond set; for chain-restricted B, orbit = topology class (sorted multiset of connected-path-lengths)",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F89 + " +
               "experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md + " +
               "simulations/_bond_isolate_compare_n7.py + " +
               "simulations/_bond_isolate_long_range_verify.py + " +
               "simulations/_bond_isolate_topology_classes_n7.py + " +
               "compute/RCPsiSquared.Propagate (bond-isolate --bonds i,j,...) + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs (a_{-1}=4 anchor for both decay and oscillation coefficients) + " +
               "compute/RCPsiSquared.Core/Symmetry/F73SpatialSumPurityClosurePi2Inheritance.cs (cited; same a_{-1} anchor) + " +
               "compute/RCPsiSquared.Core/Symmetry/F71MirrorSymmetryPi2Inheritance.cs (cited)")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
    }

    public override string DisplayName =>
        "F89 topology orbit closure for ρ_cc + uniform-J multi-bond XY";

    public override string Summary =>
        $"S(t) for ρ_cc depends only on S_N-orbit of bond set; for chain B, orbit = topology class (sorted path-length multiset); S(0) = (N−1)/N closed form ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("S(0) at N=7 (= 6/7)", S0ClosedForm(7));
            yield return InspectableNode.RealScalar("S(0) at N=4 (= 3/4)", S0ClosedForm(4));
            yield return InspectableNode.RealScalar("DecayRateCoefficient (= a_{-1} = 4, γ-axis)", DecayRateCoefficient);
            yield return InspectableNode.RealScalar("OscillationFrequencyCoefficient (= a_{-1} = 4, J-axis mirror)", OscillationFrequencyCoefficient);
            yield return new InspectableNode("Pi2 inheritance (a_{-1}=4 doubling on both energy axes)",
                summary: "γ-side: per-coherence Z-deph rate 2γ₀ = a_{0}·γ₀ doubles to S-decay rate 4γ₀ = a_{-1}·γ₀ on |·|² (identical to F73). J-side: H_B-eigenstate frequency 2J = a_{0}·J doubles to S-oscillation frequency 4J = a_{-1}·J on |·|² (NEW to F89). Same Pi2 ladder anchor a_{-1}=4 governs both.");
            yield return new InspectableNode("Verification anchors",
                summary: "N=7 multi-bond: 28 runs covering all 14 topology classes for k=1..6; 10 classes with ≥2 reps all show 0.00e+00 within-class diff. N=4 single-pair: 6 NN+LR pairs within 5.55e−17 via direct expm. N=7 single-NN-bond: 30 ordered pairs at 0.00e+00.");
            yield return new InspectableNode("All-isolated subclass closed form (Tier 1 derived)",
                summary: "S_(1)^m, N(t) = [(N−1)/N + 4m(N−2)(cos(4Jt)−1)/(N²(N−1))] · exp(−4γ₀t). Asymptotic rate 4γ₀ universal; m-correction has frequency 4J and vanishes at cos(4Jt)=1 (period π/(2J)). The empirical 'late-tail clustering at t≈20' is the in-phase moment t≈π/(2J)≈21 for J=0.075. See AllIsolatedClosedForm helper.");
            yield return new InspectableNode("Mixed-topology and pure-path classes (Tier 2 empirical)",
                summary: "Closed form open for non-all-isolated classes ((1,2), (2,2), (1,1,2), (3), (4), (5), (6) at N=7). Pure paths decay faster than 4γ₀ on visible time scales due to populated no-overlap-SE-DE coherences (rate 6γ₀).");
            yield return new InspectableNode("F86 contrast",
                summary: "F86's per-bond Q_peak fan (linear response ∂S/∂J_b at a chosen bond inside the full chain) breaks S_N differently than F89's uniform-J multi-bond setup. F89 does not predict the F86 fan; the two operate on different observables.");
        }
    }
}
