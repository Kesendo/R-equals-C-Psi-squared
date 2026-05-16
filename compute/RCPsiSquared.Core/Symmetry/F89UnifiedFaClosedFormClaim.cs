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
/// <para><b>Structural rule:</b> <c>odd-part(D_path) = (odd-part(k))²</c>. The
/// squarefree-k subset (k = 3..7) collapses to D = (1 or 2)·p² where p = largest
/// prime factor; higher prime-power k (k=8=2³, k=9=3²) needs the full
/// odd-part-squaring rule. The 2-adic exponent v₂(D_path) is captured by the
/// <c>E(k)</c> closed form below; derived from the Chebyshev pipeline (closed
/// 2026-05-15, see <c>docs/proofs/PROOF_F89_PATH_D_CLOSED_FORM.md</c> § "Tier-1-Derived
/// closure"). Probe history: <c>simulations/_f89_path_d_structure_probe.py</c>.</para>
///
/// <para>Polynomial degree = F_a count − 1 = floor(N_block/2) − 1 (interpolation
/// through F_a count distinct y_n values); path-7: cubic via cyclotomic Φ_9 = y⁶+y³+1.
/// Sum F_a · N²(N−1) is rational across all paths via Newton's identities on the
/// cyclotomic minimal polynomial of y.</para>
///
/// <para><b>Denominator closed form (Tier-1-Derived 2026-05-15):</b>
/// <c>D_k = (odd(k))² · 2^E(k)</c> where
/// <c>E(k) = max(0, ⌊(k-5)/2⌋) + v₂(k) + max(0, v₂(k) - 2)</c>. Three additive
/// contributions: polynomial-degree term (max(0, ⌊(k-5)/2⌋)), k-self 2-adic term
/// (v₂(k)), deep-2-power bonus (kicks in at v₂(k) ≥ 3). Both <see cref="PredictDenominator"/>
/// and <see cref="PathPolynomial"/> are symbolically derived via the Chebyshev-expansion
/// + orbit-polynomial-reduction pipeline, available natively in C# as
/// <see cref="F89PathPolynomialPipeline"/> (<see cref="ComputePathPolynomialBig"/>);
/// originally prototyped in <c>simulations/f89_pathk_symbolic_derivation.py</c>. The F_a
/// eigenvector ansatz <c>v_n[(i, (j, l))] = sign(i − other) · ψ_n(other) / √k</c>
/// reduces |S_c(n)|² and ‖Mv(n)‖² to polynomials in <c>c = cos(πn/(k+2)) = y_n/4</c>
/// via the Chebyshev identity <c>sin((j+1)θ) = U_j(c)·sin θ</c>; the orbit minimal
/// polynomial then gives <c>(P_k, D_k)</c> exactly. The <c>(odd(k))²</c> factor traces
/// to the <c>1/√k</c> eigenvector normalisation squared; the 2-power <c>2^E(k)</c>
/// arises from Chebyshev <c>U_j</c> leading-coefficient growth <c>2^j</c> combined with
/// the polynomial-degree reduction. All three E(k) terms structurally accounted for
/// per <c>docs/proofs/PROOF_F89_PATH_D_CLOSED_FORM.md</c> § "Structural origin of D_k".</para>
///
/// <para>Tabulation extent: k=3..9 hand-derived (path-3 algebraic anchor <c>(33+14√5)/9</c>);
/// k=10..46 cached as int-typed (<see cref="PathPolynomial"/>), all bit-exact match against
/// the native C# pipeline. For k ≥ 47 the int-typed signature overflows (D_47 = 4,632,608,768
/// exceeds int.MaxValue); use <see cref="ComputePathPolynomialBig"/> for BigInteger output
/// at arbitrary k.</para>
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
    // k=33..46 extracted via Chebyshev pipeline (simulations/f89_pathk_extend_k25_plus.py
    // 2026-05-15, ~12 sec total for k=33..46). All D bit-exact against PredictDenominator.
    // k=46 is the last int-safe path: D_46 = 1,109,393,408 < int.MaxValue; D_47 = 4,632,608,768
    // exceeds int.MaxValue and would require a long-typed Denominator signature refactor.
    private static readonly double[] _path33Coefs = { 223363072.0, 69074944.0, 107479040.0, 177192960.0, -87449600.0, -166555648.0, 30879744.0, 69736960.0, -2010880.0, -14503680.0, -813568.0, 1579424.0, 170800.0, -85960.0, -12240.0, 1846.0, 307.0 };
    private static readonly double[] _path34Coefs = { 127008768.0, 15728640.0, 3932160.0, 160235520.0, 39936000.0, -143327232.0, -35721216.0, 56788992.0, 14153472.0, -11298816.0, -2816000.0, 1185152.0, 295376.0, -62464.0, -15568.0, 1304.0, 325.0 };
    private static readonly double[] _path35Coefs = { 574062592.0, 180551680.0, -207093760.0, 284590080.0, 533217280.0, -190676992.0, -371122176.0, 46173184.0, 124189184.0, -343296.0, -21995520.0, -1499200.0, 2116064.0, 240192.0, -104288.0, -15052.0, 2062.0, 343.0 };
    private static readonly double[] _path36Coefs = { 161251328.0, 77660160.0, -113704960.0, -28344320.0, 236912640.0, 59064320.0, -156059648.0, -38906880.0, 49769984.0, 12408064.0, -8470144.0, -2111680.0, 787584.0, 196352.0, -37688.0, -9396.0, 726.0, 181.0 };
    private static readonly double[] _path37Coefs = { 1342242816.0, 907673600.0, -333905920.0, -1110507520.0, 750387200.0, 1419837440.0, -365068288.0, -758757376.0, 61208576.0, 210016768.0, 3677440.0, -32279936.0, -2515968.0, 2777472.0, 327744.0, -125000.0, -18248.0, 2290.0, 381.0 };
    private static readonly double[] _path38Coefs = { 697303040.0, 610795520.0, 152371200.0, -1043333120.0, -260177920.0, 1232732160.0, 307412992.0, -626212864.0, -156162048.0, 166062080.0, 41411840.0, -24616960.0, -6138880.0, 2053120.0, 512000.0, -89920.0, -22424.0, 1608.0, 401.0 };
    private static readonly double[] _path39Coefs = { 2893938688.0, 2025062400.0, 2019164160.0, -1767178240.0, -3846963200.0, 1765539840.0, 3392405504.0, -629952512.0, -1448644608.0, 69211136.0, 340144640.0, 11451392.0, -46065664.0, -3962624.0, 3581312.0, 436240.0, -148240.0, -21852.0, 2530.0, 421.0 };
    private static readonly double[] _path40Coefs = { 749600768.0, 389283840.0, 891944960.0, 222494720.0, -1732771840.0, -432209920.0, 1442029568.0, 359690240.0, -588335104.0, -146750464.0, 132909056.0, 33152000.0, -17414656.0, -4343808.0, 1315552.0, 328144.0, -53096.0, -13244.0, 886.0, 221.0 };
    private static readonly double[] _path41Coefs = { 6690177024.0, 2073559040.0, 3848273920.0, 7621181440.0, -5592514560.0, -11136925696.0, 3716349952.0, 7421632512.0, -989974528.0, -2614765568.0, 58064896.0, 531097600.0, 24859648.0, -64170496.0, -5954816.0, 4546400.0, 568656.0, -174152.0, -25888.0, 2782.0, 463.0 };
    private static readonly double[] _path42Coefs = { 3714056192.0, 461373440.0, 115343360.0, 7007109120.0, 1748172800.0, -9769582592.0, -2437349376.0, 6195740672.0, 1545736192.0, -2094923776.0, -522649600.0, 410820608.0, 102493184.0, -48152576.0, -12013312.0, 3321984.0, 828784.0, -124288.0, -31008.0, 1944.0, 485.0 };
    private static readonly double[] _path43Coefs = { 16424370176.0, 5161091072.0, -7465861120.0, 12728401920.0, 26886144000.0, -14630977536.0, -28546301952.0, 7103692800.0, 15109447680.0, -1416724480.0, -4502781952.0, 8155136.0, 803604480.0, 46356480.0, -87528960.0, -8626368.0, 5692704.0, 728160.0, -202880.0, -30380.0, 3046.0, 507.0 };
    private static readonly double[] _path44Coefs = { 4521984000.0, 2191523840.0, -4029153280.0, -1005322240.0, 12143165440.0, 3030056960.0, -12257624064.0, -3058614272.0, 6208421888.0, 1549172736.0, -1782390784.0, -444755968.0, 308023296.0, 76860416.0, -32620928.0, -8139840.0, 2069696.0, 516448.0, -72152.0, -18004.0, 1062.0, 265.0 };
    private static readonly double[] _path45Coefs = { 37380685824.0, 25346179072.0, -12773752832.0, -44572344320.0, 40608727040.0, 82557927424.0, -33715978240.0, -66611675136.0, 12476219392.0, 28979167232.0, -1821343744.0, -7450025984.0, -110256128.0, 1183248384.0, 79069184.0, -117200000.0, -12130560.0, 7041344.0, 918112.0, -234568.0, -35352.0, 3322.0, 553.0 };
    private static readonly double[] _path46Coefs = { 19293798400.0, 16894656512.0, 4217372672.0, -42131783680.0, -10514595840.0, 72948908032.0, 18205573120.0, -56184799232.0, -14021820416.0, 23490985984.0, 5862555648.0, -5836570624.0, -1456611328.0, 899883008.0, 224580608.0, -86837248.0, -21671680.0, 5097728.0, 1272224.0, -166336.0, -41512.0, 2312.0, 577.0 };

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
            33 => (_path33Coefs, 17842176),
            34 => (_path34Coefs, 9469952),
            35 => (_path35Coefs, 40140800),
            36 => (_path36Coefs, 10616832),
            37 => (_path37Coefs, 89718784),
            38 => (_path38Coefs, 47316992),
            39 => (_path39Coefs, 199360512),
            40 => (_path40Coefs, 52428800),
            41 => (_path41Coefs, 440664064),
            42 => (_path42Coefs, 231211008),
            43 => (_path43Coefs, 969408512),
            44 => (_path44Coefs, 253755392),
            45 => (_path45Coefs, 2123366400),
            46 => (_path46Coefs, 1109393408),
            _ => throw new ArgumentOutOfRangeException(nameof(k), k,
                "PathPolynomial (int-typed Denominator) is tabulated for path-3..46. " +
                "k=46 is the last int-safe path: D_47 = 4,632,608,768 exceeds int.MaxValue. " +
                "For k ≥ 47, use ComputePathPolynomialBig(k) which returns BigInteger " +
                "coefficients via the native F89PathPolynomialPipeline (no tabulation cap)."),
        };
    }

    /// <summary>Native runtime computation of (P_k coefficients low-to-high, D_k)
    /// for arbitrary k ≥ 3 via the F89 Chebyshev pipeline
    /// (<see cref="F89PathPolynomialPipeline.Compute"/>). Returns BigInteger so it
    /// extends past the int.MaxValue boundary at k=47. Bit-exact match against the
    /// tabulated <see cref="PathPolynomial"/> for k=3..46; sole source for k ≥ 47.
    ///
    /// <para>For k ≤ 46 prefer <see cref="PathPolynomial"/> (cached, double-typed,
    /// no BigInteger overhead). Use this method when you need k ≥ 47, when you
    /// want to verify the pipeline against the tabulation, or when working
    /// natively in BigInteger arithmetic.</para></summary>
    public static (BigInteger[] CoefficientsLowToHigh, BigInteger Denominator) ComputePathPolynomialBig(int k)
        => F89PathPolynomialPipeline.Compute(k);

    /// <summary>Predicted denominator D_k from the closed-form pattern
    /// (Tier-1-Derived 2026-05-15 via the Chebyshev pipeline):
    /// <code>
    ///   D_k = (odd(k))² · 2^E(k)
    ///   E(k) = max(0, ⌊(k-5)/2⌋) + v₂(k) + max(0, v₂(k) - 2)
    /// </code>
    /// Three additive contributions: polynomial-degree term, k-self 2-adic term,
    /// deep-2-power bonus (kicks in at v₂(k) ≥ 3). For k ∈ {3..9} this returns
    /// the same value as <see cref="PathPolynomial"/>'s tabulated denominator;
    /// for k ≥ 10 it predicts D without requiring P_k extraction. Bit-exact
    /// match against <see cref="ComputePathPolynomialBig"/> at every tested k
    /// including v₂(k) = 3 (k=200) and v₂(k) = 2 with large k (k=100, 300).
    ///
    /// <para>Anchors: <c>docs/proofs/PROOF_F89_PATH_D_CLOSED_FORM.md</c> §
    /// "Tier-1-Derived closure achieved via Chebyshev pipeline" +
    /// <c>simulations/_f89_path_d_structure_probe.py</c> +
    /// <c>simulations/_f89_path_d_verify_k16_k17.py</c> +
    /// <c>simulations/_f89_path_d_extend_k18_k24.py</c>.</para></summary>
    public static int PredictDenominator(int k)
    {
        var (oddPart, e) = DenominatorExponent(k);
        return oddPart * oddPart * (1 << e);
    }

    /// <summary>Shared (odd(k), E(k)) factorisation used by both <see cref="PredictDenominator"/>
    /// (int-typed) and <see cref="PredictDenominatorBig"/> (BigInteger-typed). Throws if k &lt; 3.</summary>
    private static (int OddPart, int Exponent) DenominatorExponent(int k)
    {
        if (k < 3) throw new ArgumentOutOfRangeException(nameof(k), k, "k must be ≥ 3.");
        int v2 = V2(k);
        int oddPart = k >> v2;
        int e = Math.Max(0, (k - 5) / 2) + v2 + Math.Max(0, v2 - 2);
        return (oddPart, e);
    }

    private static int V2(int n) => n <= 0 ? 0 : BitOperations.TrailingZeroCount(n);

    /// <summary>Evaluate sigs[F_a:n](N) for a specific path, Bloch index n, and total
    /// qubit count N. Throws if N is too small for the block (N must be ≥ N_block + 1 = k + 2).
    ///
    /// <para>Path range: k=3..46 uses the cached <see cref="PathPolynomial"/> tabulation
    /// (double-typed denominator); k ≥ 47 uses the native Chebyshev pipeline
    /// (<see cref="ComputePathPolynomialBig"/>, BigInteger). Both routes are bit-exact
    /// equivalent for k=3..46; the BigInteger denominator is cast to double for the
    /// final σ_n evaluation, which is fine in practice for k ≲ 50 and degrades
    /// gracefully past that.</para></summary>
    public static double Sigma(int k, int n, int blochN)
    {
        if (k < 3) throw new ArgumentOutOfRangeException(nameof(k), k, "Path k must be ≥ 3.");
        int nBlock = k + 1;
        if (blochN < nBlock + 1)
            throw new ArgumentOutOfRangeException(nameof(blochN), blochN,
                $"N must be ≥ N_block + 1 = {nBlock + 1} for path-{k} (need at least one bare site).");
        double y = F89PathKAtLockMechanismClaim.BlochEigenvalueY(nBlock, n);
        double poly;
        double denomD;
        if (k <= 46)
        {
            var (coefs, denom) = PathPolynomial(k);
            poly = EvaluatePolynomialAtY(coefs, y);
            denomD = denom;
        }
        else
        {
            var (bigCoefs, bigDenom) = GetCachedPipelineResult(k);
            poly = EvaluateBigIntPolynomialAtY(bigCoefs, y);
            denomD = (double)bigDenom;
        }
        return poly / (denomD * blochN * blochN * (blochN - 1));
    }

    // Cache pipeline output for k ≥ 47. SigmaSum loops Sigma over the orbit (FA = ⌊(k+1)/2⌋
    // entries per k), so without the cache each orbit summation triggers FA full pipeline runs
    // for the same k. Memoising by k cuts that to a single run per k.
    private static readonly System.Collections.Concurrent.ConcurrentDictionary<int, (BigInteger[] Coefs, BigInteger Denominator)> _pipelineCache = new();
    private static (BigInteger[] Coefs, BigInteger Denominator) GetCachedPipelineResult(int k) =>
        _pipelineCache.GetOrAdd(k, ComputePathPolynomialBig);

    /// <summary>Sum of sigs[F_a:n](N) over the S_2-anti Bloch orbit n. Rational
    /// across all paths (Newton's identities cancel the radical content):
    /// path-3: 22/3, path-4: 25/2, path-5: 483/25, path-6: 256/9 (in units of
    /// N²(N-1)). Works for arbitrary k ≥ 3 via the same routing as <see cref="Sigma"/>.</summary>
    public static double SigmaSum(int k, int blochN)
    {
        if (k < 3) throw new ArgumentOutOfRangeException(nameof(k), k, "Path k must be ≥ 3.");
        int nBlock = k + 1;
        var orbit = F89PathKAtLockMechanismClaim.SeAntiBlochOrbit(nBlock);
        double sum = 0.0;
        foreach (var n in orbit) sum += Sigma(k, n, blochN);
        return sum;
    }

    private static double EvaluatePolynomialAtY(double[] coefsLowToHigh, double y)
    {
        double poly = 0.0;
        double yPow = 1.0;
        foreach (var c in coefsLowToHigh)
        {
            poly += c * yPow;
            yPow *= y;
        }
        return poly;
    }

    private static double EvaluateBigIntPolynomialAtY(BigInteger[] coefsLowToHigh, double y)
    {
        double poly = 0.0;
        double yPow = 1.0;
        foreach (var c in coefsLowToHigh)
        {
            poly += (double)c * yPow;
            yPow *= y;
        }
        return poly;
    }

    /// <summary>BigInteger-typed denominator closed form D_k = (odd(k))²·2^E(k) for
    /// arbitrary k ≥ 3. Use when k ≥ 47 (the int-typed <see cref="PredictDenominator"/>
    /// overflows at D_47 = 4,632,608,768 > int.MaxValue). For k ≤ 46 returns the
    /// same value as <see cref="PredictDenominator"/>.</summary>
    public static BigInteger PredictDenominatorBig(int k)
    {
        var (oddPart, e) = DenominatorExponent(k);
        BigInteger oddBig = new BigInteger(oddPart);
        return oddBig * oddBig * BigInteger.Pow(2, e);
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
