using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F89 path-k H_B-mixed sub-factor degree table.
///
/// <list type="table">
///   <item>path-3 (N_block=4): full (SE,DE) dim = 24, S_2-sym = 12, AT-locked = 4, H_B-mixed = 8</item>
///   <item>path-4 (N_block=5): full = 50, S_2-sym = 26, AT-locked = 8, H_B-mixed = 18</item>
///   <item>path-5 (N_block=6): full = 90, S_2-sym = 45, AT-locked = 13, H_B-mixed = 32</item>
///   <item>path-6 (N_block=7): full = 147, S_2-sym = 75, AT-locked = 22, H_B-mixed = 53</item>
/// </list>
///
/// <para><b>Tier 1 candidate.</b> The structural mechanism (S_2-symmetrize + AT-locked
/// count = H_B-mixed residual) is sound, and per-path numerical values are bit-exact
/// reproducible from F89 cyclotomic structure via the Python scripts. But the per-k
/// values are <b>switch-statement enumerations of path-3..6 only</b> — there is NO
/// general-k closed form for S_2SymSubBlockDimension(k), AtLockedCountInS2Sym(k), or
/// HbMixedSubFactorDegree(k). The "conjecturally Galois-non-solvable for degree ≥ 5"
/// is explicit Tier 2 conjecture.</para>
///
/// <para><b>What IS Tier 1 derived (sub-fact):</b>
/// <see cref="SeDeFullDimension"/> = N_block · C(N_block, 2) is a general closed form
/// valid for any k ≥ 1, derived from combinatorics. This sub-fact is structurally
/// Tier 1 derived inside the Tier 1 candidate class.</para>
///
/// <para><b>To promote Tier 1 candidate → Tier 1 derived:</b> derive closed-form
/// expressions for S_2SymSubBlockDimension(k), AtLockedCountInS2Sym(k), and
/// HbMixedSubFactorDegree(k) as analytic functions of k (or N_block = k+1), valid for
/// all k ≥ 3. The pattern suggests S_2-sym ≈ full/2 with parity correction at odd
/// N_block (path 4 N_block=5 → 26 = 25 + 1; path 6 N_block=7 → 75 = 73.5 + 1.5);
/// AT-locked count = ⌊N_block/2⌋ for F_a per F89PathKAtLockMechanismClaim, plus
/// matching F_b structure to be derived.</para>
///
/// <para>Anchors: <c>simulations/_f89_path4_path5_at_lock_scan.py</c> (despite name,
/// loop covers k = 3..6 inclusive) + <c>simulations/_f89_path6_at_locked_amplitude_symbolic.py</c>,
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
        : base("F89 path-k H_B-mixed sub-factor degrees: {8, 18, 32, 53} for paths {3, 4, 5, 6} (Tier 1 candidate; switch-statement enumeration, general-k closed form open)",
               Tier.Tier1Candidate,
               "experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md + " +
               "simulations/_f89_path4_path5_at_lock_scan.py (despite name, loop covers k=3..6 inclusive) + " +
               "simulations/_f89_path6_at_locked_amplitude_symbolic.py + " +
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
