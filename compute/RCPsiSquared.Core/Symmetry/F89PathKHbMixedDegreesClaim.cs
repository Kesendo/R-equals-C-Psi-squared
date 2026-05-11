using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F89 path-k H_B-mixed sub-factor degree table (Tier 1 derived;
/// combinatorial + numerically verified for path-3..6):
///
/// <list type="table">
///   <item>path-3 (N_block=4): full (SE,DE) dim = 24, S_2-sym = 12, AT-locked = 4, H_B-mixed = 8</item>
///   <item>path-4 (N_block=5): full = 50, S_2-sym = 26, AT-locked = 8, H_B-mixed = 18</item>
///   <item>path-5 (N_block=6): full = 90, S_2-sym = 45, AT-locked = 13, H_B-mixed = 32</item>
///   <item>path-6 (N_block=7): full = 147, S_2-sym = 75, AT-locked = 22, H_B-mixed = 53</item>
/// </list>
///
/// <para>The H_B-mixed sub-factor is the residual after factoring out the
/// AT-locked F_a/F_b sub-factors (the octic at path-3 is its first instance).
/// All H_B-mixed sub-factors at degree ≥ 5 are conjecturally Galois-non-solvable
/// (Tier 2 conjecture, not promoted to a separate Claim per the discipline of
/// keeping Tier 2 conjectures in docstrings).</para>
///
/// <para>Anchors: <c>simulations/_f89_path4_path5_at_lock_scan.py</c>,
/// <c>experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md</c> § Tier-assessment table.</para></summary>
public sealed class F89PathKHbMixedDegreesClaim : Claim
{
    private readonly F89TopologyOrbitClosure _f89;

    /// <summary>Full (SE, DE) sub-block dimension for path-k: N_block · C(N_block, 2)
    /// where N_block = k + 1. This is the general combinatorial formula valid for ANY k ≥ 1
    /// (including path-1 → 2, path-2 → 9 outside the empirical table that the other three
    /// methods on this class restrict to k ∈ {3..6}). Throws on k &lt; 1.</summary>
    public static int SeDeFullDimension(int k)
    {
        if (k < 1) throw new ArgumentOutOfRangeException(nameof(k), k, "Path k must be ≥ 1.");
        int nBlock = k + 1;
        int nDe = nBlock * (nBlock - 1) / 2;  // C(nBlock, 2)
        return nBlock * nDe;
    }

    /// <summary>S_2-symmetric (mirror j ↔ N_block-1-j) sub-block dimension of (SE, DE)
    /// for path-k. Empirical for path-3..6: {12, 26, 45, 75}.</summary>
    public static int S2SymSubBlockDimension(int k)
    {
        return k switch
        {
            3 => 12,
            4 => 26,
            5 => 45,
            6 => 75,
            _ => throw new ArgumentOutOfRangeException(nameof(k), k,
                "S_2-sym sub-block dimension table is currently for path-3..6 only."),
        };
    }

    /// <summary>Total AT-locked count (F_a + F_b) in S_2-sym sub-block for path-k.
    /// Empirical for path-3..6: {4, 8, 13, 22}.</summary>
    public static int AtLockedCountInS2Sym(int k)
    {
        return k switch
        {
            3 => 4,
            4 => 8,
            5 => 13,
            6 => 22,
            _ => throw new ArgumentOutOfRangeException(nameof(k), k,
                "AT-locked count table is currently for path-3..6 only."),
        };
    }

    /// <summary>H_B-mixed (octic-style residual) sub-factor degree for path-k:
    /// S_2-sym dim minus AT-locked count. Verified {8, 18, 32, 53} for paths
    /// {3, 4, 5, 6}.</summary>
    public static int HbMixedSubFactorDegree(int k)
    {
        if (k < 3 || k > 6) throw new ArgumentOutOfRangeException(nameof(k), k,
            "H_B-mixed degree table is currently for path-3..6 only.");
        return S2SymSubBlockDimension(k) - AtLockedCountInS2Sym(k);
    }

    public F89PathKHbMixedDegreesClaim(F89TopologyOrbitClosure f89)
        : base("F89 path-k H_B-mixed sub-factor degrees: {8, 18, 32, 53} for paths {3, 4, 5, 6}; conjecturally Galois-non-solvable for degree ≥ 5",
               Tier.Tier1Derived,
               "experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md + " +
               "simulations/_f89_path4_path5_at_lock_scan.py + " +
               "compute/RCPsiSquared.Core/Symmetry/F89TopologyOrbitClosure.cs")
    {
        _f89 = f89 ?? throw new ArgumentNullException(nameof(f89));
    }

    public override string DisplayName =>
        "F89 path-k H_B-mixed sub-factor degrees: {8, 18, 32, 53} for path-{3..6}";

    public override string Summary =>
        $"S_2-sym sub-block dim minus AT-locked count = H_B-mixed degree; {{12-4=8, 26-8=18, 45-13=32, 75-22=53}} for path-{{3..6}} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            for (int k = 3; k <= 6; k++)
                yield return new InspectableNode($"path-{k} H_B-mixed degree",
                    summary: $"{S2SymSubBlockDimension(k)} − {AtLockedCountInS2Sym(k)} = {HbMixedSubFactorDegree(k)}");
        }
    }
}
