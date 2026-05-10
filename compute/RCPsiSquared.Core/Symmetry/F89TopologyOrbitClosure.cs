using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F89 topology orbit closure (Tier 1 derived, verified bit-identical at
/// N=7 across 12 topology classes and N=4 across all C(4,2)=6 site pairs).
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
/// <para>Cited scaffolding (constructor-injected as parents):</para>
/// <list type="bullet">
///   <item><see cref="F73SpatialSumPurityClosurePi2Inheritance"/>: closely
///         related closure for the (vac, SE) coherence block. F73 uses the
///         same orbit-style argument plus a uniform per-element 2γ₀ rate to
///         yield the closed exponential (1/2)·exp(−4γ₀t). The (S_1, S_2)
///         block has non-uniform per-element decay (rate 2γ₀ on overlap,
///         6γ₀ off overlap), so the F89 closure is orbit-only: it fixes the
///         bond-position dependence (constant in orbit) but not the time
///         dependence (no closed exponential).</item>
///   <item><see cref="F71MirrorSymmetryPi2Inheritance"/>: F71 is the spatial
///         Z₂ subgroup b ↔ N−2−b sitting inside the full S_N argument used
///         here. Any F71 mirror pair is in the same F89 orbit (live-checked
///         via <see cref="F71MirrorIsInSameOrbit"/>).</item>
/// </list>
///
/// <para>Tier1Derived: the orbit transitivity argument is elementary group
/// theory applied to the Lindblad equation; numerical verification at N=7
/// (multi-bond, 24 configurations across 12 topology classes, 8 with ≥ 2
/// representatives all showing 0.00e+00 within-class diff) and N=4
/// (single-pair, all 6 NN + long-range pairs identical within 5.55e−17 via
/// direct expm) corroborates the proof at machine precision.</para>
///
/// <para>The empirical late-tail clustering observation (at t ≫ 1, S(t)
/// classes cluster by isolated-edge count; all-isolated classes share a slow
/// tail close to F73's vac-SE rate 4γ₀) is Tier 2 empirical (open
/// derivation) and lives in the experiment writeup, not as a separate
/// Claim.</para>
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
    private readonly F73SpatialSumPurityClosurePi2Inheritance _f73;
    private readonly F71MirrorSymmetryPi2Inheritance _f71;

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

        foreach (int b in bonds)
            if (b < 0 || b > n - 2)
                throw new ArgumentOutOfRangeException(nameof(bonds), b,
                    $"Bond index {b} out of range [0, {n - 2}].");
        if (bonds.Distinct().Count() != bonds.Count)
            throw new ArgumentException("Bond indices must be distinct.", nameof(bonds));

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

    /// <summary>Live drift check: S(0) closed form at N=7 equals 6/7 to machine
    /// precision.</summary>
    public bool S0AtN7IsSixSevenths() => Math.Abs(S0ClosedForm(7) - 6.0 / 7.0) < 1e-15;

    /// <summary>Live drift check: any F71 mirror partner is in the same F89 orbit
    /// (F71 is the spatial-Z₂ subgroup of the full S_N argument used here).</summary>
    public bool F71MirrorIsInSameOrbit(int n, int b)
    {
        int partner = _f71.MirrorPair(n, b);
        return AreInSameChainOrbit(n, new[] { b }, new[] { partner });
    }

    /// <summary>Live drift check: F73's analogous closure at t=0 evaluates to 1/2 = a_2,
    /// confirming the cited-parent F73 still resolves correctly.</summary>
    public bool F73AnalogConsistent() => _f73.ClosureAtTZeroIsHalf();

    public F89TopologyOrbitClosure(
        F73SpatialSumPurityClosurePi2Inheritance f73,
        F71MirrorSymmetryPi2Inheritance f71)
        : base("F89 topology orbit closure: S(t) for ρ_cc + uniform-J multi-bond XY depends only on the S_N-orbit of the bond set; for chain-restricted B, orbit = topology class (sorted multiset of connected-path-lengths)",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F89 + " +
               "experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md + " +
               "simulations/_bond_isolate_compare_n7.py + " +
               "simulations/_bond_isolate_long_range_verify.py + " +
               "simulations/_bond_isolate_topology_classes_n7.py + " +
               "compute/RCPsiSquared.Propagate (bond-isolate --bonds i,j,...) + " +
               "compute/RCPsiSquared.Core/Symmetry/F73SpatialSumPurityClosurePi2Inheritance.cs (cited) + " +
               "compute/RCPsiSquared.Core/Symmetry/F71MirrorSymmetryPi2Inheritance.cs (cited)")
    {
        _f73 = f73 ?? throw new ArgumentNullException(nameof(f73));
        _f71 = f71 ?? throw new ArgumentNullException(nameof(f71));
    }

    public override string DisplayName =>
        "F89 topology orbit closure for ρ_cc + uniform-J multi-bond XY";

    public override string Summary =>
        $"S(t) for ρ_cc depends only on S_N-orbit of bond set; for chain B, orbit = topology class (sorted path-length multiset); S(0) = (N−1)/N closed form ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F89 topology orbit closure",
                summary: "S(t; H_B, ρ_cc) depends only on the S_N-orbit of B; for chain-restricted B, orbit = topology class (sorted multiset of connected-path-lengths)");
            yield return new InspectableNode("Cited scaffolding",
                summary: "F73 (closely related closed-form closure for vac-SE block, same orbit-style argument); F71 (spatial Z₂ subgroup of S_N, b ↔ N−2−b)");
            yield return InspectableNode.RealScalar("S(0) at N=7 (= 6/7)", S0ClosedForm(7));
            yield return InspectableNode.RealScalar("S(0) at N=4 (= 3/4)", S0ClosedForm(4));
            yield return new InspectableNode("Verification anchors",
                summary: "N=7 multi-bond: 24 runs covering 12 topology classes for k=1..6; 8 classes with ≥2 reps all show 0.00e+00 within-class diff. N=4 single-pair: 6 NN+LR pairs within 5.55e−17 via direct expm. N=7 single-NN-bond: 30 ordered pairs at 0.00e+00.");
            yield return new InspectableNode("Late-tail clustering (Tier 2 empirical)",
                summary: "Open derivation: at t ≫ 1, S(t) classes cluster by isolated-edge count in topology; all-isolated classes ((1), (1,1), (1,1,1) at N=7) share a slow tail close to F73 vac-SE rate 4γ₀ = 0.20 (empirical Γ ≈ 0.191).");
            yield return new InspectableNode("F71 mirror is F89 orbit subgroup",
                summary: "F71's b ↔ N−2−b mirror is one Z₂ inside the full S_N. Any F71 partner pair lives in the same F89 orbit (single-bond at any position is class (1)).");
            yield return new InspectableNode("F86 contrast",
                summary: "F86's per-bond Q_peak fan (linear response ∂S/∂J_b at a chosen bond inside the full chain) breaks S_N differently than F89's uniform-J multi-bond setup. F89 does not predict the F86 fan; the two operate on different observables.");
        }
    }
}
