using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F89 mixed-topology additive identity (Tier 1 derived from
/// Lindbladian factorisation; 27/27 N=7 CSVs verified at 5.013·10⁻⁷):
///
/// <code>
///   S_T(t) = Σ_i S_(k_i)(t) − (m − 1)·N·S_bare(t; N)
///   S_bare(t; N) = (N − 1)/N² · exp(−4γ₀ t)
/// </code>
///
/// <para>Derivation, empirical anchor, and reduction-of-open-work argument:
/// see <c>experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md</c> § "Mixed-topology
/// additive identity" and <c>simulations/_f89_mixed_topology_additive.py</c>.</para></summary>
public sealed class F89AdditiveIdentityClaim : Claim
{
    private readonly F89TopologyOrbitClosure _f89;

    /// <summary>Per-bare-site closed form S_bare(t; N) = (N − 1)/N² · exp(−4γ₀ t).
    /// Returns the contribution of ONE bare site to S(t) at the given (N, γ₀, t).</summary>
    public static double BarePerSite(int n, double gammaZero, double t)
    {
        if (n < 2) throw new ArgumentOutOfRangeException(nameof(n), n, "N must be ≥ 2.");
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        if (t < 0) throw new ArgumentOutOfRangeException(nameof(t), t, "t must be ≥ 0.");
        return (double)(n - 1) / (n * n) * Math.Exp(-4.0 * gammaZero * t);
    }

    /// <summary>The "(m − 1)·N" overcounting-cancellation coefficient for a
    /// topology T = (k_1, ..., k_m) at N qubits. The S_bare(t; N) per-bare-site
    /// formula is multiplied by this coefficient and SUBTRACTED from
    /// Σ_i S_(k_i)(t) to recover S_T(t).</summary>
    public static int OvercountingCoefficient(int m, int n)
    {
        if (m < 1) throw new ArgumentOutOfRangeException(nameof(m), m, "Block count m must be ≥ 1.");
        if (n < 2) throw new ArgumentOutOfRangeException(nameof(n), n, "N must be ≥ 2.");
        return (m - 1) * n;
    }

    /// <summary>Bare-site count for a topology T = (k_1, ..., k_m) at N qubits:
    /// N − Σ_i (k_i + 1). Throws if the topology requires more sites than N.</summary>
    public static int BareSiteCount(IReadOnlyList<int> kValues, int n)
    {
        if (kValues is null) throw new ArgumentNullException(nameof(kValues));
        if (n < 2) throw new ArgumentOutOfRangeException(nameof(n), n, "N must be ≥ 2.");
        int blockSites = kValues.Sum(k => k + 1);
        int bare = n - blockSites;
        if (bare < 0) throw new ArgumentException(
            $"Topology ({string.Join(",", kValues)}) needs {blockSites} block sites; only {n} qubits available.",
            nameof(kValues));
        return bare;
    }

    /// <summary>Apply the additive identity: given the per-pure-path-k_i closed
    /// forms `S_paths[i](t)`, return S_T(t) for topology T = kValues at N qubits.
    ///
    /// <para>S_T(t) = Σ_i S_paths[i](t) − (m − 1)·N·S_bare(t; N)</para></summary>
    public static double Combine(
        IReadOnlyList<int> kValues,
        int n,
        double gammaZero,
        double t,
        Func<int, int, double, double, double> sPathK)
    {
        if (kValues is null) throw new ArgumentNullException(nameof(kValues));
        if (sPathK is null) throw new ArgumentNullException(nameof(sPathK));
        _ = BareSiteCount(kValues, n);  // validates topology fits in n qubits
        double total = 0;
        foreach (var k in kValues)
            total += sPathK(k, n, gammaZero, t);
        int m = kValues.Count;
        total -= OvercountingCoefficient(m, n) * BarePerSite(n, gammaZero, t);
        return total;
    }

    public F89AdditiveIdentityClaim(F89TopologyOrbitClosure f89)
        : base("F89 mixed-topology additive identity: S_T(t) = Σ_i S_(k_i)(t) − (m−1)·N·S_bare(t; N), with S_bare = (N−1)/N²·exp(−4γ₀t); reduces 14 per-class closed forms to 6 pure-path-k forms + 1 rule",
               Tier.Tier1Derived,
               "experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md + " +
               "simulations/_f89_mixed_topology_additive.py + " +
               "compute/RCPsiSquared.Core/Symmetry/F89TopologyOrbitClosure.cs")
    {
        _f89 = f89 ?? throw new ArgumentNullException(nameof(f89));
    }

    public override string DisplayName =>
        "F89 mixed-topology additive identity from Lindbladian factorisation";

    public override string Summary =>
        $"S_T(t) = Σ_i S_(k_i)(t) − (m−1)·N·S_bare(t; N); per-bare-site formula (N−1)/N²·exp(−4γ₀t); empirical anchor 27/27 N=7 CSVs at 5.013e-7 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("Sample S_bare at γ₀=0.05, t=0, N=7", BarePerSite(7, 0.05, 0));
            yield return InspectableNode.RealScalar("Sample S_bare at γ₀=0.05, t=10, N=7", BarePerSite(7, 0.05, 10));
            yield return InspectableNode.RealScalar("OvercountingCoefficient(m=2, N=7)", OvercountingCoefficient(2, 7));
        }
    }
}
