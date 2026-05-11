using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F89 universal AT-lock mechanism for path-k (SE, DE) sub-blocks
/// (Tier 1 derived; verified across path-3..6 numerically; commits 17c4187,
/// 0d97a14, etc.):
///
/// <para>For any path-k topology with N_block = k+1 ≥ 2 sites, the (SE, DE)
/// sector of L_super contains AT-locked eigenvectors satisfying:</para>
///
/// <list type="number">
///   <item>F_a modes (rate 2γ): eigenvector supported entirely on overlap
///   basis pairs (i ∈ {j, k} for SE site i and DE pair (j, k)). Frequencies
///   match SE-anti single-particle Bloch eigenvalues
///   E_n = 4J·cos(πn/(N_block+1)) for n in the S_2-anti orbit
///   {2, 4, ..., 2·floor(N_block/2)}.</item>
///
///   <item>F_b modes (rate 6γ): eigenvector supported entirely on no-overlap
///   basis pairs (i ∉ {j, k}). Universally INVISIBLE to S(t) per-site
///   reduction: w[l] picks elements requiring overlap, so sigs[F_b] ≈ 0
///   to machine precision (10⁻³⁰..10⁻⁶³ across paths verified).</item>
///
///   <item>F_a count = floor(N_block/2) (number of SE-anti single-particle
///   Bloch modes that contribute non-trivial amplitudes).</item>
/// </list>
///
/// <para>Anchors: <c>simulations/_f89_path3_at_lock_mechanism.py</c> +
/// <c>_f89_path4_path5_at_lock_scan.py</c> + <c>experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md</c>
/// § "AT-lock mechanism".</para></summary>
public sealed class F89PathKAtLockMechanismClaim : Claim
{
    private readonly F89TopologyOrbitClosure _f89;

    /// <summary>F_a mode count for path-k with block size N_block = k+1.
    /// Equals floor(N_block/2) = number of S_2-anti single-particle Bloch modes
    /// E_n = 4J·cos(πn/(N_block+1)) for n in the orbit {2, 4, ..., 2·floor(N_block/2)}.</summary>
    public static int FaCount(int nBlock)
    {
        if (nBlock < 2) throw new ArgumentOutOfRangeException(nameof(nBlock), nBlock, "N_block must be ≥ 2.");
        return nBlock / 2;
    }

    /// <summary>The S_2-antisymmetric Bloch index orbit {2, 4, ..., 2·floor(N_block/2)}.
    /// These are the n values for which E_n = 4J·cos(πn/(N_block+1)) gives F_a frequencies.</summary>
    public static IReadOnlyList<int> SeAntiBlochOrbit(int nBlock)
    {
        if (nBlock < 2) throw new ArgumentOutOfRangeException(nameof(nBlock), nBlock, "N_block must be ≥ 2.");
        var orbit = new List<int>();
        for (int n = 2; n <= nBlock; n += 2) orbit.Add(n);
        return orbit;
    }

    /// <summary>SE single-particle Bloch eigenvalue y_n = 4·cos(πn/(N_block+1))
    /// (in units where J=1). At J ≠ 1, multiply by J to get physical eigenvalue.</summary>
    public static double BlochEigenvalueY(int nBlock, int n)
    {
        if (nBlock < 2) throw new ArgumentOutOfRangeException(nameof(nBlock), nBlock, "N_block must be ≥ 2.");
        if (n < 1 || n > nBlock) throw new ArgumentOutOfRangeException(nameof(n), n, $"n must be in [1, {nBlock}].");
        return 4.0 * Math.Cos(Math.PI * n / (nBlock + 1));
    }

    public F89PathKAtLockMechanismClaim(F89TopologyOrbitClosure f89)
        : base("F89 universal AT-lock mechanism for path-k (SE, DE) sub-blocks: F_a (rate 2γ) eigvecs supported entirely on overlap basis pairs at SE-anti single-particle Bloch frequencies y_n = 4cos(πn/(N_block+1)) for n ∈ {2, 4, ..., 2·floor(N_block/2)}; F_b (rate 6γ) eigvecs supported entirely on no-overlap, universally invisible to S(t) per-site reduction; F_a count = floor(N_block/2)",
               Tier.Tier1Derived,
               "experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md + " +
               "simulations/_f89_path3_at_lock_mechanism.py + " +
               "simulations/_f89_path4_path5_at_lock_scan.py + " +
               "compute/RCPsiSquared.Core/Symmetry/F89TopologyOrbitClosure.cs")
    {
        _f89 = f89 ?? throw new ArgumentNullException(nameof(f89));
    }

    public override string DisplayName =>
        "F89 universal AT-lock mechanism: overlap-only (F_a) / no-overlap-only (F_b) eigvec support across path-k";

    public override string Summary =>
        $"F_a count = floor(N_block/2); F_a freqs = SE-anti Bloch y_n = 4cos(πn/(N_block+1)); F_b sigs ≈ 0 to S(t) (universal); verified path-3..6 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F_a count, path-3 (N_block=4)",
                summary: $"{FaCount(4)} (= floor(4/2))");
            yield return new InspectableNode("F_a count, path-6 (N_block=7)",
                summary: $"{FaCount(7)} (= floor(7/2))");
            yield return new InspectableNode("Bloch y at path-3 n=2",
                summary: $"4cos(2π/5) = {BlochEigenvalueY(4, 2):F6}");
        }
    }
}
