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
/// <para>Late-tail clustering by isolated-edge count is documented as Tier 2
/// empirical (open derivation) in the experiment writeup.</para>
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

    public F89TopologyOrbitClosure()
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
    { }

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
            yield return new InspectableNode("Verification anchors",
                summary: "N=7 multi-bond: 24 runs covering 12 topology classes for k=1..6; 8 classes with ≥2 reps all show 0.00e+00 within-class diff. N=4 single-pair: 6 NN+LR pairs within 5.55e−17 via direct expm. N=7 single-NN-bond: 30 ordered pairs at 0.00e+00.");
            yield return new InspectableNode("Late-tail clustering (Tier 2 empirical)",
                summary: "Open derivation: at t ≫ 1, S(t) classes cluster by isolated-edge count in topology; all-isolated classes ((1), (1,1), (1,1,1) at N=7) share a slow tail close to F73 vac-SE rate 4γ₀ = 0.20 (empirical Γ ≈ 0.191).");
            yield return new InspectableNode("F86 contrast",
                summary: "F86's per-bond Q_peak fan (linear response ∂S/∂J_b at a chosen bond inside the full chain) breaks S_N differently than F89's uniform-J multi-bond setup. F89 does not predict the F86 fan; the two operate on different observables.");
        }
    }
}
