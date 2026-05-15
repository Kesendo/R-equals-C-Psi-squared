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
/// v₂(D_path) is captured by the empirical <c>E(k)</c> closed form below (verified
/// bit-exact k=3..24); what stays open is its <i>derivation</i>: no derivation from
/// k alone or from the cyclotomic discriminant of Φ_{k+2} (Angle B negative). Probe:
/// <c>simulations/_f89_path_d_structure_probe.py</c>.</para>
///
/// <para>Polynomial degree = F_a count − 1 = floor(N_block/2) − 1 (interpolation
/// through F_a count distinct y_n values); path-7: cubic via cyclotomic Φ_9 = y⁶+y³+1.
/// Sum F_a · N²(N−1) is rational across all paths via Newton's identities on the
/// cyclotomic minimal polynomial of y.</para>
///
/// <para><b>Denominator closed form (Tier-1-Derived 2026-05-15, verified k=3..24):</b>
/// <c>D_k = (odd(k))² · 2^E(k)</c> where
/// <c>E(k) = max(0, ⌊(k-5)/2⌋) + v₂(k) + max(0, v₂(k) - 2)</c>. Three additive
/// contributions: polynomial-degree term (max(0, ⌊(k-5)/2⌋)), k-self 2-adic term
/// (v₂(k)), deep-2-power bonus (kicks in at v₂(k) ≥ 3). Both <see cref="PredictDenominator"/>
/// and <see cref="PathPolynomial"/> are now symbolically derived for k = 3..24 via
/// the Chebyshev-expansion + orbit-polynomial-reduction pipeline in
/// <c>simulations/f89_pathk_symbolic_derivation.py</c>: the F_a eigenvector ansatz
/// <c>v_n[(i, (j, l))] = sign(i − other) · ψ_n(other) / √k</c> reduces |S_c(n)|² and
/// ‖Mv(n)‖² to polynomials in <c>c = cos(πn/(k+2)) = y_n/4</c> via the Chebyshev
/// identity <c>sin((j+1)θ) = U_j(c)·sin θ</c>; the orbit minimal polynomial then
/// gives <c>(P_k, D_k)</c> exactly. The <c>(odd(k))²</c> factor traces to the
/// <c>1/√k</c> eigenvector normalisation squared; the 2-power <c>2^E(k)</c> arises
/// from Chebyshev <c>U_j</c> leading-coefficient growth <c>2^j</c> combined with
/// the polynomial-degree reduction. Closes Gaps 1-3 in PROOF_F89_PATH_D_CLOSED_FORM.md.</para>
///
/// <para>Tabulation extent: k=3..9 hand-derived (path-3 algebraic anchor <c>(33+14√5)/9</c>);
/// k=10..24 symbolically derived via the pipeline above (15 new closed-form polynomials
/// added 2026-05-15). For k ≥ 25 the pipeline extends as O(k²) sympy work; not currently
/// tabulated in this typed claim.</para>
///
/// <para><b>Where D_k sits (two-layer reading):</b> D_k is the denominator of an
/// eigenvector-derived amplitude, the F89 <i>amplitude layer</i>, not the AT-governed
/// <i>eigenvalue layer</i> (λ_n = −2γ₀ + i·y_n). It is the F89 analogue of F86's g_eff:
/// a non-primitive, which is why this closed form is Tier-1-Candidate and its 2-adic
/// part resists derivation. The AT-governed closure that holds absolutely is F89c, the
/// Hamming-complement pair-sum 2γ₀·N (<see cref="AbsorptionTheoremClaim.HammingComplementPairSum"/>).
/// Two-layer section: <c>docs/proofs/PROOF_F89_PATH_D_CLOSED_FORM.md</c>.</para>
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
    // k=10..24 derived symbolically via Chebyshev-expansion + orbit-polynomial reduction
    // (simulations/f89_pathk_symbolic_derivation.py 2026-05-15). D_k bit-exact match
    // against PredictDenominator k=10..24; closes the Tier-1-Derived gap on D_k.
    private static readonly double[] _path10Coefs = { 1120.0, 128.0, 32.0, 152.0, 37.0 };
    private static readonly double[] _path11Coefs = { 6392.0, 2016.0, -512.0, 212.0, 262.0, 43.0 };
    private static readonly double[] _path12Coefs = { 2184.0, 992.0, -344.0, -84.0, 102.0, 25.0 };
    private static readonly double[] _path13Coefs = { 19472.0, 12864.0, 480.0, -2120.0, 40.0, 346.0, 57.0 };
    private static readonly double[] _path14Coefs = { 10752.0, 9472.0, 2336.0, -1984.0, -488.0, 264.0, 65.0 };
    private static readonly double[] _path15Coefs = { 47072.0, 33920.0, 19904.0, -848.0, -4528.0, -252.0, 442.0, 73.0 };
    private static readonly double[] _path16Coefs = { 12768.0, 7040.0, 8224.0, 2032.0, -1928.0, -476.0, 166.0, 41.0 };
    private static readonly double[] _path17Coefs = { 133184.0, 40704.0, 38912.0, 38432.0, -1296.0, -7880.0, -688.0, 550.0, 91.0 };
    private static readonly double[] _path18Coefs = { 84480.0, 10240.0, 2560.0, 32640.0, 8080.0, -6400.0, -1584.0, 408.0, 101.0 };
    private static readonly double[] _path19Coefs = { 419712.0, 132352.0, -72704.0, 52416.0, 74592.0, 0.0, -12320.0, -1292.0, 670.0, 111.0 };
    private static readonly double[] _path20Coefs = { 128128.0, 60160.0, -43136.0, -10688.0, 31232.0, 7744.0, -4856.0, -1204.0, 246.0, 61.0 };
    private static readonly double[] _path21Coefs = { 1097984.0, 734720.0, -70912.0, -313216.0, 70144.0, 135680.0, 4096.0, -17992.0, -2088.0, 802.0, 133.0 };
    private static readonly double[] _path22Coefs = { 585728.0, 514048.0, 127744.0, -290816.0, -72192.0, 111104.0, 27584.0, -13888.0, -3448.0, 584.0, 145.0 };
    private static readonly double[] _path23Coefs = { 2489856.0, 1764352.0, 1277440.0, -381440.0, -793600.0, 82432.0, 230144.0, 12240.0, -25040.0, -3100.0, 946.0, 157.0 };
    private static readonly double[] _path24Coefs = { 658944.0, 351232.0, 541696.0, 134656.0, -345088.0, -85760.0, 92256.0, 22928.0, -9512.0, -2364.0, 342.0, 85.0 };
    // k=25..32 extracted via same Chebyshev pipeline (simulations/f89_pathk_extend_k25_plus.py
    // 2026-05-15, ~4 sec total). All D bit-exact against PredictDenominator.
    private static readonly double[] _path25Coefs = { 6308864.0, 1943552.0, 2441216.0, 3216384.0, -830976.0, -1674240.0, 72960.0, 367584.0, 25872.0, -33608.0, -4352.0, 1102.0, 183.0 };
    private static readonly double[] _path26Coefs = { 3727360.0, 458752.0, 114688.0, 2838528.0, 706048.0, -1404928.0, -349440.0, 289408.0, 71984.0, -25216.0, -6272.0, 792.0, 197.0 };
    private static readonly double[] _path27Coefs = { 17430528.0, 5488640.0, -4653056.0, 4864000.0, 7984128.0, -1400320.0, -3170816.0, 17472.0, 558752.0, 46624.0, -43840.0, -5868.0, 1270.0, 211.0 };
    private static readonly double[] _path28Coefs = { 5048320.0, 2408448.0, -2625536.0, -653312.0, 3466240.0, 862720.0, -1297280.0, -322880.0, 216640.0, 53920.0, -16280.0, -4052.0, 454.0, 113.0 };
    private static readonly double[] _path29Coefs = { 42471424.0, 28606464.0, -6635520.0, -22294528.0, 9932800.0, 17955840.0, -1980928.0, -5561984.0, -117504.0, 815552.0, 76320.0, -55880.0, -7672.0, 1450.0, 241.0 };
    private static readonly double[] _path30Coefs = { 22282240.0, 19529728.0, 4866048.0, -20807680.0, -5181440.0, 15224832.0, 3791360.0, -4462592.0, -1111296.0, 624384.0, 155488.0, -41152.0, -10248.0, 1032.0, 257.0 };
    private static readonly double[] _path31Coefs = { 93315072.0, 65601536.0, 56459264.0, -32727040.0, -66641920.0, 18495488.0, 36776960.0, -2328832.0, -9197056.0, -375936.0, 1151040.0, 116976.0, -69872.0, -9788.0, 1642.0, 273.0 };
    private static readonly double[] _path32Coefs = { 24371200.0, 12779520.0, 24494080.0, 6103040.0, -29552640.0, -7362560.0, 15256064.0, 3800832.0, -3630848.0, -904576.0, 435872.0, 108592.0, -25544.0, -6364.0, 582.0, 145.0 };

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
            10 => (_path10Coefs, 200),
            11 => (_path11Coefs, 968),
            12 => (_path12Coefs, 288),
            13 => (_path13Coefs, 2704),
            14 => (_path14Coefs, 1568),
            15 => (_path15Coefs, 7200),
            16 => (_path16Coefs, 2048),
            17 => (_path17Coefs, 18496),
            18 => (_path18Coefs, 10368),
            19 => (_path19Coefs, 46208),
            20 => (_path20Coefs, 12800),
            21 => (_path21Coefs, 112896),
            22 => (_path22Coefs, 61952),
            23 => (_path23Coefs, 270848),
            24 => (_path24Coefs, 73728),
            25 => (_path25Coefs, 640000),
            26 => (_path26Coefs, 346112),
            27 => (_path27Coefs, 1492992),
            28 => (_path28Coefs, 401408),
            29 => (_path29Coefs, 3444736),
            30 => (_path30Coefs, 1843200),
            31 => (_path31Coefs, 7872512),
            32 => (_path32Coefs, 2097152),
            _ => throw new ArgumentOutOfRangeException(nameof(k), k,
                "Unified F_a closed form is tabulated for path-3..32 (k=3..9 hand-derived; " +
                "k=10..32 symbolically derived via Chebyshev-expansion + orbit-polynomial " +
                "reduction, see simulations/f89_pathk_symbolic_derivation.py + " +
                "simulations/f89_pathk_extend_k25_plus.py). For k ≥ 33, run the symbolic " +
                "pipeline with extended range and add the resulting coefficient arrays."),
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
