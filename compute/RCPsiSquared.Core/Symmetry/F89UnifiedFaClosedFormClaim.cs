using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F89 unified F_a AT-locked amplitude closed form across path-3..9
/// (Tier 1 derived; bit-exact verified):
///
/// <code>
///   sigs[F_a:n](N) = P_path(y_n) / [D_path · N²·(N−1)]
///   y_n = 4·cos(πn/(N_block+1))   for n in S_2-anti orbit
/// </code>
///
/// <para>Per-path (P_path, D_path) integer-coefficient polynomial table:</para>
/// <list type="table">
///   <item>path-3 (N_block=4): P(y) = 14y + 47, D = 9</item>
///   <item>path-4 (N_block=5): P(y) = 10y + 25, D = 4</item>
///   <item>path-5 (N_block=6): P(y) = 13y² + 82y + 129, D = 25</item>
///   <item>path-6 (N_block=7): P(y) = 17y² + 72y + 80, D = 18</item>
///   <item>path-7 (N_block=8): P(y) = 21y³ + 130y² + 292y + 382, D = 98 = 2·7²</item>
///   <item>path-8 (N_block=9): P(y) = 13y³ + 54y² + 68y + 110, D = 32 = 2⁵ (extracted via C2FullBlockSigmaAnatomy 2026-05-13)</item>
///   <item>path-9 (N_block=10): P(y) = 31y⁴ + 190y³ + 288y² + 440y + 1476, D = 324 = 2²·3⁴ (extracted via C2FullBlockSigmaAnatomy 2026-05-13)</item>
/// </list>
///
/// <para><b>Empirical structural rule (verified path-3..9):</b>
/// <c>odd-part(D_path) = (odd-part(k))²</c>. The squarefree-k subset (k = 3..7)
/// collapses to D = (1 or 2)·p² where p = largest prime factor; higher prime-power
/// k (k=8=2³, k=9=3²) needs the full odd-part-squaring rule. The 2-adic exponent
/// v₂(D_path) sequence 0, 2, 0, 1, 1, 5, 2 for k = 3..9 has NO closed form derived
/// from k alone or from the cyclotomic discriminant of Φ_{k+2}; this 2-adic part
/// remains an OPEN question. Probe: <c>simulations/_f89_path_d_structure_probe.py</c>.</para>
///
/// <para>Polynomial degree = F_a count − 1 = floor(N_block/2) − 1 (interpolation
/// through F_a count distinct y_n values); path-7: cubic via cyclotomic Φ_9 = y⁶+y³+1.
/// Sum F_a · N²(N−1) is rational across all paths via Newton's identities on the
/// cyclotomic minimal polynomial of y.</para>
///
/// <para><b>Denominator closed form (Tier-1-Candidate 2026-05-13, verified k=3..24):</b>
/// <c>D_k = (odd(k))² · 2^E(k)</c> where
/// <c>E(k) = max(0, ⌊(k-5)/2⌋) + v₂(k) + max(0, v₂(k) - 2)</c>. Three additive
/// contributions: polynomial-degree term (max(0, ⌊(k-5)/2⌋)), k-self 2-adic term
/// (v₂(k)), deep-2-power bonus (kicks in at v₂(k) ≥ 3). The full closed form for
/// (P_path, D_path) remains incomplete: D_k is closed (see <see cref="PredictDenominator"/>);
/// P_k coefficients are still tabulated per path. Proof:
/// <c>docs/proofs/PROOF_F89_PATH_D_CLOSED_FORM.md</c>.</para>
///
/// <para>Anchors: <c>simulations/_f89_path3_at_locked_amplitude_symbolic.py</c>,
/// <c>_f89_path4_at_locked_amplitude_symbolic.py</c>,
/// <c>_f89_path5_at_locked_amplitude_symbolic.py</c>,
/// <c>_f89_path6_at_locked_amplitude_symbolic.py</c>,
/// <c>_f89_path5_cardano_cubic_individuals.py</c>,
/// <c>experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md</c> § "Unified closed form".</para></summary>
public sealed class F89UnifiedFaClosedFormClaim : Claim
{
    // Parent-edge marker for Schicht-1 wiring (consumed by ClaimRegistryBuilder; not used in this class body).
    private readonly F89TopologyOrbitClosure _f89;
    // Parent-edge marker for Schicht-1 wiring (consumed by ClaimRegistryBuilder; not used in this class body).
    private readonly F89PathKAtLockMechanismClaim _atLock;

    // Cached coefficient arrays (low-to-high degree) for the per-path polynomials
    // P_path(y) used in sigs[F_a:n](N) = P_path(y_n) / [D_path · N²(N-1)].
    private static readonly double[] _path3Coefs = { 47.0, 14.0 };          // 14y + 47
    private static readonly double[] _path4Coefs = { 25.0, 10.0 };          // 10y + 25
    private static readonly double[] _path5Coefs = { 129.0, 82.0, 13.0 };   // 13y² + 82y + 129
    private static readonly double[] _path6Coefs = { 80.0, 72.0, 17.0 };    // 17y² + 72y + 80
    private static readonly double[] _path7Coefs = { 382.0, 292.0, 130.0, 21.0 };  // 21y³ + 130y² + 292y + 382
    private static readonly double[] _path8Coefs = { 110.0, 68.0, 54.0, 13.0 };          // 13y³ + 54y² + 68y + 110 (extracted)
    private static readonly double[] _path9Coefs = { 1476.0, 440.0, 288.0, 190.0, 31.0 }; // 31y⁴ + 190y³ + 288y² + 440y + 1476 (extracted)

    /// <summary>Per-path (P_path coefficients low-to-high degree, D_path) table.
    /// Returns the integer-coefficient numerator polynomial coefficients for
    /// sigs[F_a:n](N) · D_path · N²(N-1) = P_path(y_n) at y_n = 4·cos(πn/(N_block+1)).</summary>
    public static (double[] CoefficientsLowToHigh, int Denominator) PathPolynomial(int k)
    {
        return k switch
        {
            3 => (_path3Coefs, 9),
            4 => (_path4Coefs, 4),
            5 => (_path5Coefs, 25),
            6 => (_path6Coefs, 18),
            7 => (_path7Coefs, 98),
            8 => (_path8Coefs, 32),
            9 => (_path9Coefs, 324),
            _ => throw new ArgumentOutOfRangeException(nameof(k), k,
                "Unified F_a closed form is currently tabulated for path-3..9 only. " +
                "Path-10+ extensions: extract via C2FullBlockSigmaAnatomy; coefficients open."),
        };
    }

    /// <summary>Predicted denominator D_k from the empirical closed form
    /// (Tier-1-Candidate, verified bit-exact for k=3..24):
    /// <code>
    ///   D_k = (odd(k))² · 2^E(k)
    ///   E(k) = max(0, ⌊(k-5)/2⌋) + v₂(k) + max(0, v₂(k) - 2)
    /// </code>
    /// Three additive contributions: polynomial-degree term, k-self 2-adic term,
    /// deep-2-power bonus (kicks in at v₂(k) ≥ 3). For k ∈ {3..9} this returns
    /// the same value as <see cref="PathPolynomial"/>'s tabulated denominator;
    /// for k ≥ 10 it predicts D without requiring P_k extraction.
    ///
    /// <para>Anchors: <c>docs/proofs/PROOF_F89_PATH_D_CLOSED_FORM.md</c> +
    /// <c>simulations/_f89_path_d_structure_probe.py</c> +
    /// <c>simulations/_f89_path_d_verify_k16_k17.py</c> +
    /// <c>simulations/_f89_path_d_extend_k18_k24.py</c>.</para></summary>
    public static int PredictDenominator(int k)
    {
        if (k < 3) throw new ArgumentOutOfRangeException(nameof(k), k, "k must be ≥ 3.");
        int v2 = V2(k);
        int oddPart = k >> v2;
        int polyDegreeContribution = Math.Max(0, (k - 5) / 2);
        int kSelfContribution = v2;
        int deepBonusContribution = Math.Max(0, v2 - 2);
        int e = polyDegreeContribution + kSelfContribution + deepBonusContribution;
        return oddPart * oddPart * (1 << e);
    }

    private static int V2(int n) => n <= 0 ? 0 : BitOperations.TrailingZeroCount(n);

    /// <summary>Evaluate sigs[F_a:n](N) for a specific path, Bloch index n, and total qubit count N.
    /// Throws if N is too small for the block (N must be ≥ N_block + 1 = k + 2).
    /// Note: P_k tabulation is available for k ∈ {3..9} only; D_k for k ≥ 10 is available
    /// via <see cref="PredictDenominator"/>.</summary>
    public static double Sigma(int k, int n, int blochN)
    {
        if (k < 3 || k > 9) throw new ArgumentOutOfRangeException(nameof(k), k, "Path k ∈ {3, 4, 5, 6, 7, 8, 9} only.");
        int nBlock = k + 1;
        if (blochN < nBlock + 1)
            throw new ArgumentOutOfRangeException(nameof(blochN), blochN,
                $"N must be ≥ N_block + 1 = {nBlock + 1} for path-{k} (need at least one bare site).");
        var (coefs, denom) = PathPolynomial(k);
        double y = F89PathKAtLockMechanismClaim.BlochEigenvalueY(nBlock, n);
        double poly = 0.0;
        double yPow = 1.0;
        foreach (var c in coefs)
        {
            poly += c * yPow;
            yPow *= y;
        }
        return poly / (denom * blochN * blochN * (blochN - 1));
    }

    /// <summary>Sum of sigs[F_a:n](N) over the S_2-anti Bloch orbit n. Rational
    /// across all paths (Newton's identities cancel the radical content):
    /// path-3: 22/3, path-4: 25/2, path-5: 483/25, path-6: 256/9 (in units of
    /// N²(N-1)).</summary>
    public static double SigmaSum(int k, int blochN)
    {
        if (k < 3 || k > 9) throw new ArgumentOutOfRangeException(nameof(k), k, "Path k ∈ {3, 4, 5, 6, 7, 8, 9} only.");
        int nBlock = k + 1;
        var orbit = F89PathKAtLockMechanismClaim.SeAntiBlochOrbit(nBlock);
        double sum = 0.0;
        foreach (var n in orbit) sum += Sigma(k, n, blochN);
        return sum;
    }

    public F89UnifiedFaClosedFormClaim(F89TopologyOrbitClosure f89, F89PathKAtLockMechanismClaim atLock)
        : base("F89 unified F_a AT-locked amplitude closed form across path-3..9: sigs[F_a:n](N) = P_path(y_n)/[D_path·N²(N-1)] with y_n = 4cos(πn/(N_block+1)) on the SE-anti Bloch orbit; (P_path, D_path) tabulated per path; sum F_a is rational across all paths via Newton's identities on the cyclotomic minimal polynomial",
               Tier.Tier1Derived,
               "experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md + " +
               "simulations/_f89_path3_at_locked_amplitude_symbolic.py + " +
               "simulations/_f89_path4_at_locked_amplitude_symbolic.py + " +
               "simulations/_f89_path5_at_locked_amplitude_symbolic.py + " +
               "simulations/_f89_path6_at_locked_amplitude_symbolic.py + " +
               "simulations/_f89_path5_cardano_cubic_individuals.py + " +
               "compute/RCPsiSquared.Core/Symmetry/F89PathKAtLockMechanismClaim.cs + " +
               "simulations/_f89_to_f86_kbond_via_eigendecomp.py + compute/RCPsiSquared.Core/F86/Item1Derivation/C2FullBlockSigmaAnatomy.cs (path-7 extraction) + " +
               "simulations/_f89_path_d_structure_probe.py (path-8/9 extraction + odd-part rule)")
    {
        _f89 = f89 ?? throw new ArgumentNullException(nameof(f89));
        _atLock = atLock ?? throw new ArgumentNullException(nameof(atLock));
    }

    public override string DisplayName =>
        "F89 unified F_a AT-locked amplitude closed form: sigs = P_path(y_n) / [D_path·N²(N-1)] across path-3..9";

    public override string Summary =>
        $"sigs[F_a:n](N) = P_path(y_n)/[D_path·N²(N-1)]; (P,D) = {{(14y+47,9), (10y+25,4), (13y²+82y+129,25), (17y²+72y+80,18), (21y³+130y²+292y+382,98), (13y³+54y²+68y+110,32), (31y⁴+190y³+288y²+440y+1476,324)}} for path-{{3,4,5,6,7,8,9}}; odd-part(D)=(odd-part(k))² verified; sum rational via Newton ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Sample sigs path-3 N=11 n=2",
                summary: $"{Sigma(3, 2, 11):G6} (= (33+14√5)/[9·11²·10])");
            yield return new InspectableNode("Sample sigs path-6 N=11 n=4 (zero mode)",
                summary: $"{Sigma(6, 4, 11):G6} (= 40/[9·11²·10])");
            yield return new InspectableNode("Sum F_a path-5 N=11",
                summary: $"{SigmaSum(5, 11):G6} (= 483/[25·11²·10])");
        }
    }
}
