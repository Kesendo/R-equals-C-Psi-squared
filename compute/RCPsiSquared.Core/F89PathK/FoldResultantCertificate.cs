using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.F89PathK;

/// <summary>The ONE-WAY fold-resultant certificate for remainder R1 of the codim-1-by-additivity proof
/// (docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md §6, the sectorbraid arc): at every branch locus q* ≠ 0 of
/// the (1,2) residual factor F_res (F_18 R-even / the degree-17 R-odd factor at N=5), the holomorphic
/// fold μ = −λ_A − 2N is NOT an eigenvalue of the corner block (p_c+1, p_c+1). Formally, for
///   R(q) = Res_Λ(F_res(Λ, q), F2_corner(−Λ − 4N, q))   and   D(q) = disc_Λ(F_res)(q)
/// (×2-cleared variable Λ = 2λ, so the physical fold λ ↦ −λ − 2N reads Λ ↦ −Λ − 4N), the certificate
/// proves gcd(R, D) = c·q^e over Q(i): the ONLY common root is q = 0. D vanishes at every locus where
/// F_res has a repeated root, so this covers every branch locus at once, defective AND diabolic (the
/// radical trick, no discriminant split extracted); q = 0 itself is harmless for R1 (the block at q = 0
/// is the real diagonal A, semisimple, so no Jordan chain exists there to transport). The certificate is
/// ONE-WAY: a q-power gcd proves absence; anything else is refine-not-refute.
///
/// <para>THE EXACT BIVARIATE LAYER. F_full(Λ, q) is computed EXACTLY over Z[i][q] by the division-free
/// Berkowitz algorithm on the integer pencil A + q·C; the AT factor is exact and bivariate (its strands
/// are q-linear on the q-independent rate-confined invariant subspace,
/// <see cref="F89AtFactorReconstruction.ClearedAtSectors"/>), and F_res = F_full/AT by exact Λ-long
/// division (remainder asserted zero). The corner polynomial is the exact Berkowitz charpoly of the
/// ×2-cleared corner pencil, fold-composed. The bivariate objects are cross-checked at q0 = 2, 3 against
/// the independent per-point path (<see cref="ExactSample"/>: GaussianMatrixCharpoly + DivMod).</para>
///
/// <para>THE ALL-MOD-p MULTI-PRIME PROOF. Per split prime p ≡ 1 (mod 4) (i ↦ r, r² ≡ −1), R mod p and
/// D mod p are recovered exactly by sampling the Sylvester resultant/discriminant of the reduced
/// bivariate polynomials at q0 = 0, 1, 2, … and interpolating: their q-degrees are PROVEN bounded by
/// deg_q R ≤ deg(F_res)·deg(corner) − m_R and deg_q D ≤ deg(F_res)·(deg(F_res)−1) − 2·m_D, where m_R,
/// m_D count the multiplicity collisions of the exact q → ∞ leading forms (each branch obeys
/// Λ_i(q) = c_i·q + O(1) by Bauer-Fike, both hopping directions being diagonal-similar to i×(real
/// symmetric), verified exactly per run; coincident-c pairs contribute O(1) to the root-difference
/// products). These bounds size the interpolation; they need NOT be tight (the corner-fold bound is
/// tight, the identity-composition bound is loose: the identity leading forms share a common factor
/// beyond the isolated-root collisions m_R counts). THE DEGREE CERTIFICATE (guard G2, landed
/// 2026-07-03): deg(R mod p) &lt; deg_q R exactly when the CHOSEN prime ideal π_p = (p, i − r) divides
/// lc_q(R) ∈ Z[i] (the reduction map's kernel); distinct sampled primes choose ideals over distinct p,
/// pairwise coprime, so ∏ N(π_p) = ∏ p ≤ N(lc_q R) ≤ ‖R‖² (Hadamard), and with every p ≥ 2³⁰ at most
/// 2·log₂‖R‖/30 sampled primes can lose the degree. Sampling MORE distinct split primes than that
/// bound forces trueDeg = max_p deg(R mod p) = deg_q R EXACTLY, so a prime is USED only if
/// deg(R mod p) attains the certified trueDeg (⟹ π_p ∤ lc_q(R); G1 is automatic, everything is monic
/// in Λ). At every used prime the certificate demands gcd(R mod p, D mod p) = c·q^e. THE LIFT: the
/// primitive H = gcd(R, D) ∈ Z[i][q] has lc_q(H) | lc_q(R) (Gauss), so at every used prime (π_p ∤
/// lc_q(R) ⟹ π_p ∤ lc_q(H)) H mod p keeps its degree and divides the q-power gcd, i.e. EVERY
/// sub-leading coefficient of H is divisible by every used prime π. Since H | R, Mignotte/Landau bounds H's coefficient heights
/// by B_H = 2^{deg R}·‖R‖ with ‖R‖ ≤ the Hadamard row-norm product of the exact Sylvester matrix
/// (computed from the exact bivariate coefficients); once ∏ p exceeds the Gaussian-norm bound 4·B_H²,
/// the sub-leading coefficients vanish identically and H = lc·q^e EXACTLY. This is a COMPLETE proof, not
/// a heuristic: no common root of R and D exists off q = 0.</para>
///
/// <para>THE CERTIFIED DISC-LAYER READING (the β-exotic exclusion; landed 2026-07-09 after three empty
/// reviews; spec in ClaudeTasks/BETA_EXOTIC_DISC_MULTIPLICITY_CERTIFICATE.md). The report carries D's
/// squarefree layer degrees, the N=5 image of the N=4 disc structure disc = q²⁴·(diabolic)²·P_20
/// (experiments/F89_BRANCH_LOCUS_PALINDROME.md): the multiplicity-1 layer carries the √-branch (defective)
/// loci, the multiplicity-2 layer the squared diabolic factor, and R shares nothing with either (the layer
/// gcds are reported). Beyond diagnostics, MaxDiscMultiplicity is a PROOF once DiscLayersCertified holds:
/// a branch locus with Puiseux exponent 3/2 (the β-exotic, normal form [[0,s],[s²,0]], eigenvalues ±s^{3/2})
/// forces ord_{q*} disc_Λ(F_res) ≥ 3, because every other colliding pair contributes NON-NEGATIVE order; so
/// MaxDiscMultiplicity ≤ 2 excludes the β-exotic at every branch locus q ≠ 0 of that block and parity.
/// Deliberately weaker than squarefreeness: the diabolic loci are genuine DOUBLE roots and must survive.
/// This is why the discriminant-Galois route (dead: a full Galois group sees which sheets swap, never the
/// ½-vs-3/2 exponent) is not needed. THE LIFT is one-way: reduction mod π_p ∤ lc_q(D) is a degree-preserving
/// homomorphism, so a factor (q − α)^m survives and distinct roots can only MERGE, giving
/// max-mult(D mod p) ≥ max-mult(D). The layer prime is therefore certified at BOTH ends of the q-line: it
/// attains TrueDiscriminantDegree (no root escaped to q = ∞) and TrueQValuationD (no nonzero root collapsed
/// onto q = 0, where StripQ would have discarded it with the q-power), both certified by the same
/// Hadamard/lc-divisor device as trueDegR, against LcDivisorBoundD. DiscriminantDegree and QValuationD are
/// reported AT the layer prime. Companion fact: the AT subspace is q-independent and invariant for every q,
/// so D and K each preserve it and (both being Hermitian) also its orthogonal complement: AT ⊕ F_res is
/// genuinely block-DIAGONAL, no Jordan chain crosses the seam, and, the AT side being itself semisimple
/// at every q (AT UNIFORM γ, where per rate sector D is one scalar, so the restriction is r₀·I + i·q·K|_W
/// with K diagonal-similar to real symmetric, diagonalizable on every invariant subspace), disc_Λ(F_res)
/// misses no defect. THE RATE SCOPE, and this sentence is the ORIGINAL that two derivative sites quote
/// (<c>experiments/F89_PATH_K_DIABOLIC.md</c> and the arc <c>sideways_spin_ladder</c>, both fenced
/// 2026-08-22): "D is one scalar per rate sector" IS the uniform-γ hypothesis, unnamed. Under a per-site
/// profile the sector diagonal carries one rate-weighted subset sum per cell instead of one scalar, the
/// restriction is no longer r₀·I + i·q·K, and this leg of the certificate's soundness argument fails. The
/// certificates in this file are uniform-γ objects for a second and independent reason: they are EXACT
/// integer computations with γ = 1 built in (the target pencil's diagonal −4·n_diff, and the fold constant
/// λ₀ = −2N, which is −2σ evaluated at γ = 1), so a profile would require recomputing them rather than
/// reinterpreting them. What this
/// does NOT exclude: a cubic branch point (3×3 Jordan) has ord disc = 2 and hides in the multiplicity-2
/// layer; at a count-drop it is ruled out separately by parity (an even-order zero gives no sign change, so
/// it cannot flip two reals into a conjugate pair).</para>
///
/// <para>THE Re/Im GCD CERTIFICATE (even N, landed 2026-08-10): where the disc is genuinely COMPLEX
/// (the (1,2) block at even N ≥ 6 carries no within-block antiunitary; the self-fold is
/// block-scoped and sits on the (p, N/2) blocks there, experiments/F89_PATH_K_DIABOLIC.md's
/// mechanism paragraph), <see cref="CertifyDiscReImGcd"/> proves, when GcdDeg = 0 at a certified
/// prime and one-way as ever (anything else is refine-not-refute), the window-free statement that
/// NO real q ≠ 0 carries ANY repeated Λ-root of F_res, stronger than the scanned-strip verdict that
/// route replaced, by recovering Re D and Im D mod p separately (two embeddings of one split prime)
/// and certifying gcd(Re D, Im D) trivial: the sideways_spin_ladder arc's Sturm/gcd closing route
/// (OpenArcsRegistry, and the "Verdict, scoped" paragraph of the experiment doc). First certified
/// instance: N=6, both parities. Run: <c>dotnet test "compute/RCPsiSquared.Core.Tests" --filter
/// "Category=SLOW_EVEN_DISC|Category=EVEN_DISC_SCOUT"</c> (EvenNDiscLayerScoutTests; the scout
/// category carries the N=4 fail-closed control; the ~9 min N=8 layer scout sits apart as
/// Category=SLOW_EVEN_DISC8).</para>
///
/// <para>Run: <c>dotnet test "compute/RCPsiSquared.Diagnostics.Tests" --filter "Category=FOLDRESULTANT"
/// --logger "console;verbosity=detailed"</c> (FoldResultantCertificateTests).</para></summary>
public static class FoldResultantCertificate
{
    /// <summary>The multi-prime certificate's verdict and guard trail. Complete = true is the finished
    /// one-way proof (gcd(R, D) is exactly a q-power over Q(i)); SharedIsQPowerAtEveryPrime = false
    /// means some certified-good prime saw a shared factor beyond q^e: refine, not refute.
    /// ResultantDegree is the CERTIFIED true deg_q R (= max_p deg(R mod p), exact once PrimesSampled &gt;
    /// LcDivisorBound, the max count of sampled prime ideals that can divide lc_q(R)); it equals
    /// ResultantDegreeBound for the corner-fold and sits below it for the identity composition.</summary>
    public sealed record CompleteReport(
        int N, bool ROdd,
        int BlockDimension, int AtDegree, int ResidualDegree, int CornerDegree,
        int InfinityCollisionsR, int InfinityRepeatedD,
        int ResultantDegreeBound, int DiscriminantDegreeBound,
        int ResultantDegree, int DiscriminantDegree,
        int QValuationR, int QValuationD,
        int[] DiscLayerDegrees, int[] LayerGcdDegrees,
        int PrimesUsed, int PrimesSkipped, int PrimesSampled, int LcDivisorBound,
        long FirstPrime, long LastPrime,
        int ProofBoundDigits, int PrimeProductDigits,
        bool SharedIsQPowerAtEveryPrime, bool Complete,
        int TrueDiscriminantDegree, int TrueQValuationD, int LcDivisorBoundD,
        long LayerPrime, int MaxDiscMultiplicity, bool DiscLayersCertified);

    /// <summary>Runs the complete multi-prime fold-resultant proof for R1: the (1,2)-family residual of the
    /// N-site chain (R-even sym / R-odd 2-cycle sector) against the corner (p_c+1, p_c+1), fold-composed
    /// (μ = −λ_A − 2N ∉ spec(corner)). A thin corner-fold wrapper over <see cref="CertifyCore"/>.</summary>
    public static CompleteReport CertifyComplete(int n, bool rOdd, int maxPrimes = 5000, Action<string>? log = null)
        => CertifyCore(n, rOdd, (n + 1) / 2 + 1, (n + 1) / 2 + 1,
                       new GaussianInteger(-1, 0), new GaussianInteger(-4L * n, 0), maxPrimes, log);

    /// <summary>Runs the complete multi-prime exclusion proof for a general TARGET block (tWKet, tWBra) and
    /// Λ-composition (composeA·Λ + composeB): gcd(Res_Λ(F_res, F2_target(composeA·Λ + composeB)), disc_Λ(F_res))
    /// over Q(i) is a pure q-power, so no (1,2) residual root's composed image is a target eigenvalue at any
    /// branch locus q ≠ 0. Remainder R1 is the corner-fold case (w, w, −1, −4N); remainder R4's interior-core
    /// exclusion is the identity case, e.g. (1,1) with (1, 1, +1, 0): the braid eigenvalue λ_A itself (an
    /// F_res root) is not carried by the non-member (1,1) block.</summary>
    public static CompleteReport CertifyBlockExclusion(int n, bool rOdd, int tWKet, int tWBra,
        GaussianInteger composeA, GaussianInteger composeB, int maxPrimes = 5000, Action<string>? log = null)
        => CertifyCore(n, rOdd, tWKet, tWBra, composeA, composeB, maxPrimes, log);

    /// <summary>The certified disc-layer reading ALONE, with no corner block, no resultant and no
    /// Mignotte lift: the β-exotic exclusion at a chain length where <see cref="CertifyComplete"/> is
    /// unaffordable. Same statement as that method's MaxDiscMultiplicity, same one-way lift, same
    /// both-ends layer prime, but it never touches R.
    ///
    /// <para>Why a separate entry point (measured, gate Category=BETAN7PROBE): CertifyComplete also
    /// proves the remainder-R1 gcd, whose resultant runs against the corner block (p_c+1, p_c+1) of
    /// dimension C(N, (N+1)/2+1)² = 25 at N = 5 but 441 at N = 7. Its exact ℤ[i] charpoly and the
    /// degree-53 × 441 resultant dominate everything; the β-exclusion needs none of it. What it does
    /// need at n = 7 costs 357 ms (block charpoly) + 40 ms (divide by AT) per exact evaluation, and
    /// the discriminant is only ever taken mod p, at O(resDeg²) machine-int ops per node.</para>
    ///
    /// <para>The Mignotte/Hadamard PROOF bound is absent on purpose: it certifies the gcd lift, a
    /// statement this method does not make. What remains is the lc-divisor bound, which certifies
    /// deg_q D and v_q(D) separately; a prime attaining both is searched for and the report fails
    /// closed (<c>DiscLayersCertified = false</c>) if none of the sampled primes does.</para></summary>
    public static DiscMultiplicityReport CertifyDiscMultiplicity(
        int n, bool rOdd, int maxPrimes = 20000, Action<string>? log = null)
    {
        // n ≥ 4, EITHER parity of n: the D-only path never touches the corner block (p_c+1, p_c+1),
        // which is what CertifyCore's odd-N guard protects. Every structural assumption this path
        // does make is enforced by a throw on THIS path: the q=0 block diagonal and real and C pure
        // imaginary (ClearedAtSectors), C weighted-symmetric (AssertHoppingSymmetry, the Bauer-Fike
        // premise behind dBound), block linear in q (BlockPencil), AT divides exactly with the
        // residual degree checked, Gauss-lemma integrality of the sector charpolys. Even N first
        // exercised 2026-08-10 (the sideways_spin_ladder arc's Sturm route needs the N=6/8 disc
        // layers); N=4 doubles as the control, its disc factorization being known in closed form
        // (F89Path3OcticGaloisClaim, R-even; the R-odd twin follows by the even-N conjugation).
        if (n < 4)
            throw new ArgumentException("n ≥ 4 (the octic block).", nameof(n));

        var (aBlk, cBlk) = BlockPencil(n, rOdd);
        AssertHoppingSymmetry(cBlk, rOdd ? null : F89PathKSeDeBlock.SymOrbitSizes(n));

        var sectors = F89AtFactorReconstruction.ClearedAtSectors(n - 1, rOdd);
        int atDeg = 0;
        foreach (var s in sectors) atDeg += s.KCharpoly.Length - 1;
        int blockDim = aBlk.GetLength(0);
        int resDeg = blockDim - atDeg;

        // exact bivariate F_res over Z[i][q]: charpoly of the pencil, divided by the AT factor
        var chi = PencilCharpolyBivariate(aBlk, cBlk);
        var fRes = BivariateDivideExact(chi, BivariateAtFactor(sectors));
        if (fRes.Length - 1 != resDeg)
            throw new InvalidOperationException($"bivariate residual degree {fRes.Length - 1} ≠ {resDeg}.");

        // cross-check the bivariate layer against the independent D-only per-point path at q0 = 2, 3
        foreach (int q0 in new[] { 2, 3 })
            AssertPolyEqual(EvalBivariateAtQ(fRes, q0), ResidualPerPointAt(n, rOdd, q0),
                $"bivariate F_res at q0={q0}");

        // G2 for D alone: the q → ∞ leading form of F_res fixes m_D and the proven degree bound
        var (fResLead, lcD) = ResidualLeadingForm(n, rOdd);
        var (_, mD) = SquarefreeLayers(FromZi(fResLead));
        if ((mD == 0) != !lcD.Equals(GaussianInteger.Zero))
            throw new InvalidOperationException("m_D inconsistent with disc_x(f_res).");
        int dBound = resDeg * (resDeg - 1) - 2 * mD;

        // the lc-divisor bound for D: Hadamard on the Sylvester matrix of (F_res, F_res')
        BigInteger rowF = BigInteger.Zero;
        foreach (var coeff in fRes) rowF += HNormQ(coeff);
        BigInteger rowFp = resDeg * rowF;
        BigInteger normD = BigInteger.Pow(rowF, Math.Max(resDeg - 1, 0)) * BigInteger.Pow(rowFp, resDeg);
        int lcDivisorBoundD = (int)(2 * normD.GetBitLength() / 30) + 1;

        const int extra = 24;
        int samples = dBound + 1 + extra;
        var perPrime = new List<(int P, int DegD, int VD, int[] DStrip)>();
        int sampled = 0;

        var clock = System.Diagnostics.Stopwatch.StartNew();
        long candidate = (1L << 30) + 1;
        while (candidate % 4 != 1) candidate += 2;
        for (int tried = 0; tried < maxPrimes && sampled <= lcDivisorBoundD; candidate += 4)
        {
            if (!IsPrime(candidate)) continue;
            tried++;
            int p = checked((int)candidate);
            int? root = SqrtMinusOneFast(p);
            if (root is null) continue;

            var fResP = ReduceBivariate(fRes, p, root.Value);
            var dSamp = new int[samples];
            for (int q0 = 0; q0 < samples; q0++)
            {
                var fL = EvalBivariateModP(fResP, q0, p);
                if (DegP(fL) != resDeg || fL[resDeg] != 1)
                    throw new InvalidOperationException($"reduced residual not monic of degree {resDeg} at q0={q0}.");
                dSamp[q0] = DiscriminantModP(fL, p);
            }
            var dNodes = new int[dBound + 1];
            Array.Copy(dSamp, dNodes, dBound + 1);
            var dp = InterpolateAtIntegerNodes(dNodes, p);
            for (int q0 = dBound + 1; q0 < samples; q0++)
                if (EvalModP(dp, q0, p) != dSamp[q0])
                    throw new InvalidOperationException($"D interpolant fails verification at q0={q0} (p={p}).");

            // D ≡ 0 mod p (p divides the whole content) carries no information about either
            // invariant: skip it entirely, uncounted. Such a prime divides lc_q(D) too, so it is
            // one of the ≤ lcDivisorBoundD bad primes the sampling budget already tolerates;
            // not counting it only demands one more prime, never weakens the certificate.
            if (DegP(dp) < 0) continue;
            sampled++;
            int vDp = QValuation(dp);
            perPrime.Add((p, DegP(dp), vDp, StripQ(dp, vDp)));

            if (log is not null && sampled % 25 == 0)
                log($"sampled={sampled} (bound {lcDivisorBoundD}), p={p}, {clock.ElapsedMilliseconds} ms");
        }

        // trueDegD and eD are each certified once sampled > lcDivisorBoundD; a prime attaining BOTH is
        // then searched for (their union is NOT covered by one bound), and we fail closed if none does.
        int trueDegD = -1, eD = int.MaxValue;
        foreach (var s in perPrime)
        {
            if (s.DegD > trueDegD) trueDegD = s.DegD;
            if (s.VD < eD) eD = s.VD;
        }
        if (eD == int.MaxValue) eD = -1;

        long layerPrime = 0;
        int[] layerDegrees = Array.Empty<int>();
        foreach (var s in perPrime)
        {
            if (s.DegD != trueDegD || s.VD != eD) continue;
            layerPrime = s.P;
            var layers = FpSquarefreeLayers(s.DStrip, s.P);
            layerDegrees = new int[layers.Count];
            for (int k = 0; k < layers.Count; k++) layerDegrees[k] = DegP(layers[k]);
            break;
        }
        bool certified = layerPrime != 0 && sampled > lcDivisorBoundD && trueDegD >= 0;

        return new DiscMultiplicityReport(
            n, rOdd, blockDim, atDeg, resDeg, mD, dBound,
            trueDegD, eD, lcDivisorBoundD, layerPrime,
            layerDegrees, layerDegrees.Length, sampled, certified);
    }

    /// <summary>The verdict of <see cref="CertifyDiscMultiplicity"/>. <c>MaxDiscMultiplicity ≤ 2</c> with
    /// <c>DiscLayersCertified</c> is the β-exotic exclusion at that (n, parity): a Puiseux-3/2 point would
    /// force a disc root of multiplicity ≥ 3 off q = 0.</summary>
    public sealed record DiscMultiplicityReport(
        int N, bool ROdd, int BlockDimension, int AtDegree, int ResidualDegree,
        int InfinityRepeatedD, int DiscriminantDegreeBound,
        int TrueDiscriminantDegree, int TrueQValuationD, int LcDivisorBoundD,
        long LayerPrime, int[] DiscLayerDegrees, int MaxDiscMultiplicity,
        int PrimesSampled, bool DiscLayersCertified);

    /// <summary>Independent falsifier for <see cref="CertifyDiscMultiplicity"/>: read deg_q D, v_q(D) and
    /// the squarefree layer degrees at the <paramref name="nth"/> sampled split prime (0 = the first, the
    /// one the certificate uses). A degree, reduction or single-prime-fluke bug in the layer reading shows
    /// up as a disagreement between two primes; the certificate's soundness argument says every prime that
    /// attains trueDegD and the minimal v_q(D) must report the same layers.</summary>
    public static (int DegD, int VD, int[] LayerDegrees) DiscLayersAtNthPrime(int n, bool rOdd, int nth)
    {
        if (n < 4) throw new ArgumentException("n ≥ 4 (either parity; see CertifyDiscMultiplicity).", nameof(n));
        if (nth < 0) throw new ArgumentOutOfRangeException(nameof(nth));

        var (aBlk, cBlk) = BlockPencil(n, rOdd);
        // The falsifier must hold the same premise the certificate's dBound rests on (Bauer-Fike /
        // O(1)-branch needs C = i·K with K weighted-symmetric); CertifyDiscMultiplicity checks it,
        // and before 2026-08-10 this method silently inherited it.
        AssertHoppingSymmetry(cBlk, rOdd ? null : F89PathKSeDeBlock.SymOrbitSizes(n));
        var sectors = F89AtFactorReconstruction.ClearedAtSectors(n - 1, rOdd);
        int atDeg = 0;
        foreach (var s in sectors) atDeg += s.KCharpoly.Length - 1;
        int resDeg = aBlk.GetLength(0) - atDeg;

        var chi = PencilCharpolyBivariate(aBlk, cBlk);
        var fRes = BivariateDivideExact(chi, BivariateAtFactor(sectors));
        var (fResLead, _) = ResidualLeadingForm(n, rOdd);
        var (_, mD) = SquarefreeLayers(FromZi(fResLead));
        int dBound = resDeg * (resDeg - 1) - 2 * mD;
        int samples = dBound + 1 + 24;

        long candidate = (1L << 30) + 1;
        while (candidate % 4 != 1) candidate += 2;
        int seen = 0;
        for (; ; candidate += 4)
        {
            if (!IsPrime(candidate)) continue;
            int p = checked((int)candidate);
            int? root = SqrtMinusOneFast(p);
            if (root is null) continue;
            if (seen++ < nth) continue;

            var fResP = ReduceBivariate(fRes, p, root.Value);
            var dSamp = new int[samples];
            for (int q0 = 0; q0 < samples; q0++)
            {
                var fL = EvalBivariateModP(fResP, q0, p);
                if (DegP(fL) != resDeg || fL[resDeg] != 1)
                    throw new InvalidOperationException($"reduced residual not monic of degree {resDeg} at q0={q0}.");
                dSamp[q0] = DiscriminantModP(fL, p);
            }
            var dNodes = new int[dBound + 1];
            Array.Copy(dSamp, dNodes, dBound + 1);
            var dp = InterpolateAtIntegerNodes(dNodes, p);
            for (int q0 = dBound + 1; q0 < samples; q0++)
                if (EvalModP(dp, q0, p) != dSamp[q0])
                    throw new InvalidOperationException($"D interpolant fails verification at q0={q0} (p={p}).");

            int vDp = QValuation(dp);
            var layers = FpSquarefreeLayers(StripQ(dp, vDp), p);
            var degs = new int[layers.Count];
            for (int k = 0; k < layers.Count; k++) degs[k] = DegP(layers[k]);
            return (DegP(dp), vDp, degs);
        }
    }

    /// <summary>The ×2-cleared R-parity sector pencil A + q·C of the general (wKet, wBra) chain
    /// coherence block (<see cref="WeightCoherenceBlock.Build(int,int,int,Complex)"/>, γ = 1, J = q).
    /// SCOPE: blocks whose reflection action has NO fixed basis pair, i.e. every orbit a 2-cycle —
    /// guaranteed whenever one weight side has no reflection-symmetric configuration (any ODD weight
    /// at EVEN N: a reflection-fixed mask needs its bits in mirror pairs, impossible at odd popcount
    /// when no site is self-mirrored). Then both parities carry the uniform metric BᵀB = 2·I (sector
    /// columns e_t ± e_{Rt}, one per 2-cycle {t &lt; Rt} in increasing t), the ×2-cleared sector
    /// operator is BᵀLB = 2·UᵀLU for the orthonormal U = B/√2 with exact Gaussian-integer entries,
    /// and <see cref="AssertHoppingSymmetry"/> runs with null weights (no orbit-size skew). A fixed
    /// point throws: such a block needs the (1,2)-style orbit-sum metric instead. Both (1,2) and
    /// (1,3) at N = 6 qualify (ket weight 1 is odd), which is what lets the generalized path be
    /// validated against the committed (1,2) scout pins. Built for the sideways_spin_ladder arc's
    /// (1,3) certificate step; the (1,2) methods above keep their original
    /// <see cref="F89PathKSeDeBlock"/> route untouched.</summary>
    public static (GaussianInteger[,] A, GaussianInteger[,] C) WeightSectorPencil(
        int n, int wKet, int wBra, bool rOdd)
    {
        if (n < 2 || wKet < 0 || wKet > n || wBra < 0 || wBra > n)
            throw new ArgumentException($"need n ≥ 2 and weights in 0..n; got ({n}, {wKet}, {wBra}).");
        var b = new GaussianInteger[3][,];
        for (int q0 = 0; q0 <= 2; q0++) b[q0] = WeightSectorTwoTimes(n, wKet, wBra, rOdd, q0);
        int d = b[0].GetLength(0);
        var c = new GaussianInteger[d, d];
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
            {
                c[i, j] = b[1][i, j] - b[0][i, j];
                if (!(b[0][i, j] + c[i, j] + c[i, j]).Equals(b[2][i, j]))
                    throw new InvalidOperationException("the sector block is not linear in q0.");
            }
        return (b[0], c);
    }

    /// <summary>BᵀL(q0)B on the 2-cycle ± basis (see <see cref="WeightSectorPencil"/>): four
    /// full-block entries per sector entry, cross terms signed by the parity.</summary>
    private static GaussianInteger[,] WeightSectorTwoTimes(int n, int wKet, int wBra, bool rOdd, int q0)
    {
        var l = WeightCoherenceBlock.Build(n, wKet, wBra, new Complex(q0, 0));
        var perm = WeightCoherenceBlock.ReflectionPermutation(n, wKet, wBra);
        // R must commute with L ENTRY-WISE (each entry a single hop contribution of the
        // reflection-symmetric bond set), else BᵀLB is a projection, not a restriction, and every
        // downstream "certified" invariant would be about the wrong object (a review round named
        // this as the one uncovered spot the Sturm step must not inherit). Exact comparison: both
        // sides are the identical float construction.
        for (int t = 0; t < perm.Length; t++)
            for (int s = 0; s < perm.Length; s++)
                if (l[perm[t], perm[s]] != l[t, s])
                    throw new InvalidOperationException(
                        $"reflection does not commute entry-wise with L at ({t},{s}), q0={q0}.");
        var pairs = new List<(int T, int T2)>();
        for (int t = 0; t < perm.Length; t++)
        {
            if (perm[t] == t)
                throw new InvalidOperationException(
                    $"reflection fixed point at basis index {t}: WeightSectorPencil covers all-2-cycle " +
                    "blocks only (one odd weight side at even N guarantees that); a fixed point needs " +
                    "the orbit-sum metric of the (1,2) path.");
            if (perm[t] > t) pairs.Add((t, perm[t]));
        }
        int d = pairs.Count;
        double sgn = rOdd ? -1.0 : 1.0;
        var blk = new GaussianInteger[d, d];
        for (int r = 0; r < d; r++)
            for (int col = 0; col < d; col++)
            {
                Complex v = l[pairs[r].T, pairs[col].T]
                          + sgn * l[pairs[r].T, pairs[col].T2]
                          + sgn * l[pairs[r].T2, pairs[col].T]
                          + l[pairs[r].T2, pairs[col].T2];
                blk[r, col] = ExactGaussian(v);
            }
        return blk;
    }

    /// <summary>EXACT evenness of the bivariate about Λ₀: every odd Taylor coefficient
    /// a_m(q) = Σ_{k≥m} C(k,m)·c_k(q)·Λ₀^{k−m} of F_res at Λ₀ identically zero over ℤ[i].
    /// BigInteger binomials and powers, no reduction: a true verdict is a theorem about the
    /// exact residual, not a mod-p scout read.</summary>
    private static bool ExactOddTaylorAllZero(GaussianInteger[][] fRes, int lambda0)
    {
        int kMax = fRes.Length - 1;
        var binom = new System.Numerics.BigInteger[kMax + 1][];
        for (int k = 0; k <= kMax; k++)
        {
            binom[k] = new System.Numerics.BigInteger[k + 1];
            binom[k][0] = binom[k][k] = System.Numerics.BigInteger.One;
            for (int m = 1; m < k; m++)
                binom[k][m] = binom[k - 1][m - 1] + binom[k - 1][m];
        }
        var lamPow = new System.Numerics.BigInteger[kMax + 1];
        lamPow[0] = System.Numerics.BigInteger.One;
        for (int e = 1; e <= kMax; e++) lamPow[e] = lamPow[e - 1] * lambda0;
        int len = 0;
        foreach (var qs in fRes) len = Math.Max(len, qs.Length);
        for (int m = 1; m <= kMax; m += 2)
            for (int i = 0; i < len; i++)
            {
                var accRe = System.Numerics.BigInteger.Zero;
                var accIm = System.Numerics.BigInteger.Zero;
                for (int k = m; k <= kMax; k++)
                {
                    if (i >= fRes[k].Length) continue;
                    var c = fRes[k][i];
                    if (c.Re.IsZero && c.Im.IsZero) continue;
                    var w = binom[k][m] * lamPow[k - m];
                    accRe += c.Re * w;
                    accIm += c.Im * w;
                }
                if (!accRe.IsZero || !accIm.IsZero) return false;
            }
        return true;
    }

    /// <summary>Weight-general per-point cross-check of the bivariate residual, the guard the (1,2)
    /// siblings carry (<see cref="CertifyDiscMultiplicity"/> / <see cref="CertifyCore"/>) and the
    /// weight-general path must not drop, being exactly where a new block shape first appears: at
    /// q0 = 2, 3 the sector block is REBUILT from scratch (<see cref="WeightSectorTwoTimes"/>, a
    /// fresh WeightCoherenceBlock build; at q0 = 3 genuinely outside the pencil's construction
    /// points 0, 1, 2), its charpoly (the per-point Faddeev route, independent of the bivariate
    /// Berkowitz) divided by the AT factor at q0, the remainder asserted empty and the quotient
    /// compared to the bivariate F_res at q0. Also pins deg_Λ F_res = resDeg. One leg narrower
    /// than the (1,2) original, on purpose and recorded: the AT factor at q0 is the evaluated
    /// <see cref="BivariateAtFactor"/> rather than a second AT construction route (the (1,2)
    /// AtFactorAt is weight-hardwired); the AT is still non-trivially tested by the exact
    /// division of a fresh charpoly.</summary>
    private static void WeightCrossCheckResidual(
        GaussianInteger[][] fRes, int n, int wKet, int wBra, bool rOdd,
        IReadOnlyList<AtSector> sectors, int resDeg)
    {
        if (fRes.Length - 1 != resDeg)
            throw new InvalidOperationException($"bivariate residual degree {fRes.Length - 1} ≠ {resDeg}.");
        var at = BivariateAtFactor(sectors);
        foreach (int q0 in new[] { 2, 3 })
        {
            var blk = WeightSectorTwoTimes(n, wKet, wBra, rOdd, q0);
            var (res, rem) = GaussianPolynomial.DivMod(GaussianMatrixCharpoly.Characteristic(blk),
                                                       EvalBivariateAtQ(at, q0));
            if (rem.Length != 0)
                throw new InvalidOperationException($"AT factor does not divide the sector charpoly at q0={q0}.");
            AssertPolyEqual(EvalBivariateAtQ(fRes, q0), res, $"weight-general bivariate F_res at q0={q0}");
        }
    }

    /// <summary>The weight-general sibling of <see cref="DiscLayersAtNthPrime"/>: deg_q D, v_q(D) and
    /// the squarefree layer degrees of the AT-stripped residual discriminant of the (wKet, wBra)
    /// R-sector at the <paramref name="nth"/> split prime. The AT sectors are DERIVED from the q=0
    /// diagonal (<see cref="F89AtFactorReconstruction.ClearedAtSectorsFromPencil"/>), so the (1,2)
    /// rate-literal trap ({−4,−12} fed a {−8,−16} block silently clears nothing) cannot recur. The
    /// interpolation bound is the crude Bauer-Fike resDeg·(resDeg−1) with NO leading-form m_D
    /// subtraction (this scout skips the leading form on purpose, the extra sample nodes being
    /// cheap; the certified sibling <see cref="WeightCertifyDiscMultiplicity"/> carries the
    /// generalized leading form and the sharpened bound). Scout-grade like its sibling: one prime, layers certified
    /// only by cross-prime agreement (run nth = 0 and 1). For the (1,3)@N=6 census block the
    /// multiplicity-2 layer is the separating object: away from Λ = −2N a doubled defective locus is
    /// a DOUBLE root of the sector disc while a doubled diabolic crossing is a QUADRUPLE one (the
    /// Klein orbits of the census section in experiments/F89_PATH_K_DIABOLIC.md), so the layer
    /// profile reads the two species apart — a squarefree part or sign count would erase exactly
    /// that (docs/CAUGHT_ERRORS.md 2026-07-01 B1). The Λ = −2N stratum is NOT separated here; it
    /// gets its own read (<see cref="WeightLambdaStratumAtNthPrime"/>).</summary>
    public static (int P, int DegD, int VD, int[] LayerDegrees, int AtDeg, int ResDeg, bool DiscIsRealModP,
        int Mult2StratumGcdDeg)
        WeightDiscLayersAtNthPrime(int n, int wKet, int wBra, bool rOdd, int nth)
    {
        if (nth < 0) throw new ArgumentOutOfRangeException(nameof(nth));

        var (aBlk, cBlk) = WeightSectorPencil(n, wKet, wBra, rOdd);
        AssertHoppingSymmetry(cBlk, null);      // uniform metric: all orbits are 2-cycles
        var sectors = F89AtFactorReconstruction.ClearedAtSectorsFromPencil(aBlk, cBlk);
        int atDeg = 0;
        foreach (var s in sectors) atDeg += s.KCharpoly.Length - 1;
        int resDeg = aBlk.GetLength(0) - atDeg;

        var chi = PencilCharpolyBivariate(aBlk, cBlk);
        var fRes = BivariateDivideExact(chi, BivariateAtFactor(sectors));
        WeightCrossCheckResidual(fRes, n, wKet, wBra, rOdd, sectors, resDeg);
        int dBound = resDeg * (resDeg - 1);
        int samples = dBound + 1 + 24;

        long candidate = (1L << 30) + 1;
        while (candidate % 4 != 1) candidate += 2;
        int seen = 0;
        for (; ; candidate += 4)
        {
            if (!IsPrime(candidate)) continue;
            int p = checked((int)candidate);
            int? root = SqrtMinusOneFast(p);
            if (root is null) continue;
            if (seen++ < nth) continue;
            int r = root.Value;

            // both embeddings i ↦ r and i ↦ p−r: their D's agree iff Im D ≡ 0 mod p, the exact
            // form of the census's sector disc-reality (the float read was ≤ 5e-12); layers are
            // read from the first embedding.
            var dp = new int[2][];
            foreach (var (emb, k) in new[] { (r, 0), (p - r, 1) })
            {
                var fResP = ReduceBivariate(fRes, p, emb);
                var dSamp = new int[samples];
                for (int q0 = 0; q0 < samples; q0++)
                {
                    var fL = EvalBivariateModP(fResP, q0, p);
                    if (DegP(fL) != resDeg || fL[resDeg] != 1)
                        throw new InvalidOperationException($"reduced residual not monic of degree {resDeg} at q0={q0} (i↦{emb}).");
                    dSamp[q0] = DiscriminantModP(fL, p);
                }
                var dNodes = new int[dBound + 1];
                Array.Copy(dSamp, dNodes, dBound + 1);
                dp[k] = InterpolateAtIntegerNodes(dNodes, p);
                for (int q0 = dBound + 1; q0 < samples; q0++)
                    if (EvalModP(dp[k], q0, p) != dSamp[q0])
                        throw new InvalidOperationException($"D interpolant fails verification at q0={q0} (p={p}, i↦{emb}).");
            }
            bool discReal = dp[0].Length == dp[1].Length;
            if (discReal)
                for (int i = 0; i < dp[0].Length; i++)
                    if (dp[0][i] != dp[1][i]) { discReal = false; break; }

            int vDp = QValuation(dp[0]);
            var layers = FpSquarefreeLayers(StripQ(dp[0], vDp), p);
            var degs = new int[layers.Count];
            for (int k = 0; k < layers.Count; k++) degs[k] = DegP(layers[k]);

            // The gcd(u, disc_M G) = 1 hypothesis, gated as u | mult-2-layer (a review round found
            // the μ = 0 reading needs it: an f-root of order 2 could be one transversal crossing OR
            // two coincident defective EPs on the line; the 1+1 case would push a u-root into the
            // mult-4 layer, so deg gcd(u, L₂) = deg u excludes it and also fixes the mult-2 split
            // exactly). Computed only when f = F_res(−2n, q) is a perfect square (f-layers [0, k]);
            // −1 otherwise (e.g. the squarefree N=4 control, where the evenness analysis is void).
            int mult2Gcd = -1;
            {
                var fResP0 = ReduceBivariate(fRes, p, r);
                int lenF = 0;
                foreach (var qs in fResP0) lenF = Math.Max(lenF, qs.Length);
                var f = new int[lenF];
                long lam0 = (((-2L * n) % p) + p) % p;
                long pw = 1;
                for (int k = 0; k < fResP0.Length; k++)
                {
                    for (int i = 0; i < fResP0[k].Length; i++)
                        f[i] = (int)((f[i] + fResP0[k][i] * pw) % p);
                    pw = pw * lam0 % p;
                }
                if (DegP(f) >= 0)
                {
                    var fl = FpSquarefreeLayers(StripQ(f, QValuation(f)), p);
                    if (fl.Count == 2 && DegP(fl[0]) <= 0 && layers.Count >= 2)
                        mult2Gcd = DegP(GcdModP(fl[1], layers[1], p));
                }
            }
            return (p, DegP(dp[0]), vDp, degs, atDeg, resDeg, discReal, mult2Gcd);
        }
    }

    /// <summary>The verdict of <see cref="WeightSturmScout"/>: steps 0-2 of the Sturm design
    /// (the arc's real-root-count closing step, design note local under docs/superpowers/plans/).
    /// Sentinels when the evenness is void (the N=4 (1,2) control): MDeg = UDeg = −1,
    /// FIsPerfectSquare/CompositionIdentityAtNodes = false, MsPerNodeDiscM = 0.</summary>
    public sealed record WeightSturmScoutReport(
        int N, int WKet, int WBra, bool ROdd,
        int AtDeg, int ResDeg,
        bool ResidualIsReal,
        bool EvenInShiftedLambda,
        bool EvenInQ,
        int MDeg, int FDeg, int VQf,
        bool FIsPerfectSquare, int UDeg, bool UIsReal,
        bool CompositionIdentityAtNodes,
        double MsPerNodeDiscM,
        int FoldCheckerboardSign, int CheckerboardOddCellCount, int CheckerboardEvenCellCount,
        int FullCharpolyCheckerboardSign);

    /// <summary>The checkerboard support read of a Taylor-shifted bivariate: Sign = +1 when every
    /// coefficient with j + k odd vanishes and the survivors are real on even j / imaginary on odd
    /// j, -1 for the mirrored pattern (j + k even vanishing), 0 when neither holds. Cells are
    /// counted by j parity, which equals k parity only at Sign = +1.</summary>
    private static (int Sign, int OddCells, int EvenCells) CheckerboardRead(GaussianInteger[][] h)
    {
        bool plus = true, minus = true;
        for (int j = 0; j < h.Length && (plus || minus); j++)
            for (int k = 0; k < h[j].Length; k++)
            {
                var c = h[j][k];
                bool re = !c.Re.IsZero, im = !c.Im.IsZero;
                if (((j + k) & 1) != 0) { if (re || im) plus = false; }
                else if ((j & 1) == 0) { if (im) plus = false; }
                else if (re) plus = false;

                if (((j + k) & 1) == 0) { if (re || im) minus = false; }
                else if ((k & 1) == 0) { if (im) minus = false; }
                else if (re) minus = false;
            }
        int sign = plus ? 1 : (minus ? -1 : 0);
        int odd = 0, even = 0;
        if (sign != 0)
            for (int j = 0; j < h.Length; j++)
                for (int k = 0; k < h[j].Length; k++)
                {
                    var c = h[j][k];
                    if (c.Re.IsZero && c.Im.IsZero) continue;
                    if ((j & 1) != 0) odd++; else even++;
                }
        return (sign, odd, even);
    }

    /// <summary>Exactly the part of the checkerboard that licenses an even-in-q discriminant, read
    /// over ℤ[i] and nothing rounded: does every non-vanishing coefficient of the Λ₀-shifted
    /// bivariate share ONE parity of j + k (j the u-exponent, k the q-exponent)? Then
    /// F_res(Λ₀−u, −q) = (−1)^s·F_res(Λ₀+u, q), so the Λ-roots reflect about Λ₀ under q → −q, so
    /// disc_Λ is unchanged and D is EVEN in q. BOTH parities s license it, which is why this is a
    /// separate read from <see cref="CheckerboardRead"/>: that one also demands the realness
    /// pattern (real on even j, imaginary on odd j) and reports ±1 by which parity survives, and
    /// its Sign = 0 would decline blocks whose support parity is perfectly good. Monicity forces
    /// s ≡ deg_Λ (mod 2), so nothing here needs the fold-leg premise as an input: the support
    /// parity is the fact, fold-fixedness only the reason to expect it
    /// (<see cref="IsFoldFixedBlock"/>).</summary>
    private static bool QNegationSupportParityExact(GaussianInteger[][] fRes, int lambda0)
    {
        var h = TaylorShiftBivariate(fRes, lambda0);
        int parity = -1;
        for (int j = 0; j < h.Length; j++)
            for (int k = 0; k < h[j].Length; k++)
            {
                if (h[j][k].Re.IsZero && h[j][k].Im.IsZero) continue;
                int cell = (j + k) & 1;
                if (parity < 0) parity = cell;
                else if (parity != cell) return false;
            }
        return true;                                        // the zero polynomial cannot occur (monic)
    }

    /// <summary>The q-even sampling licence, in ONE place because the device and its cost probe
    /// must never disagree about which path a block takes: the probe exists to price the run the
    /// device will actually make, so a divergence here would price a run that then refuses to
    /// start. Returns whether to halve; throws where a FOLD-FIXED block fails the exact read,
    /// since the checkerboard's descent to F_res is measured per block and not proven, so that
    /// failure is a finding and not a nuisance. <paramref name="forceDense"/> suppresses the throw
    /// as well as the halving, which is what makes it a usable way past that finding.</summary>
    private static bool QEvenSamplingLicence(
        GaussianInteger[][] fRes, int n, int wKet, int wBra, bool forceDense)
    {
        if (forceDense) return false;
        bool qEven = QNegationSupportParityExact(fRes, -2 * n);
        if (!qEven && IsFoldFixedBlock(n, wKet, wBra))
            throw new InvalidOperationException(
                $"({wKet},{wBra})@N={n} is fold-fixed but its Λ₀-shifted F_res support is not of one " +
                "j+k parity: the checkerboard does not descend to F_res here. Re-run with " +
                "forceDenseSampling: true to certify anyway, and record the block.");
        return qEven;
    }

    /// <summary>The square abscissae m = q0² must stay pairwise distinct mod p across EVERY
    /// sampled q0, interpolation nodes and verification points alike, or a verification point
    /// would silently re-test a node it coincides with instead of checking anything; p divides
    /// (i−j)(i+j) otherwise, so 2·q0 &lt; p suffices. Hence <paramref name="lastQ0"/> is the last
    /// SAMPLED q0, not the last node. Checked per prime in both the device and the probe, so the
    /// failure names the prime rather than surfacing later as two coinciding nodes.</summary>
    private static void AssertSquareAbscissaeFit(int p, int lastQ0)
    {
        if (2L * lastQ0 >= p)
            throw new InvalidOperationException(
                $"prime {p} is too small for the square abscissae: q0² collides at or below q0={lastQ0}.");
    }

    /// <summary>Scout steps 0-2 of the Sturm design (the design note is local under
    /// docs/superpowers/plans/, NOT in the repo) on the (wKet, wBra) R-sector: §0 exact
    /// realness of F_res; §1 the exact G/f/ũ chain (G the even-Taylor reindex of F_res at
    /// Λ₀ = −2n, f = G(0, q) with v_q(f) read exactly, ũ the monic square root over ℚ(i) with
    /// the exact identity lc(f)·ũ² = f); §2 the composition identity
    /// disc_Λ(F_res) = (−4)^m·f·disc_M(G)² gated exactly at the fresh nodes q0 = 2, 3, plus the
    /// per-node exact disc_M(G) timing over <paramref name="timingNodes"/> further nodes (the
    /// Route-A cost input of the design's §5).</summary>
    public static WeightSturmScoutReport WeightSturmScout(
        int n, int wKet, int wBra, bool rOdd, int timingNodes = 10)
    {
        var (aBlk, cBlk) = WeightSectorPencil(n, wKet, wBra, rOdd);
        AssertHoppingSymmetry(cBlk, null);
        var sectors = F89AtFactorReconstruction.ClearedAtSectorsFromPencil(aBlk, cBlk);
        int atDeg = 0;
        foreach (var s in sectors) atDeg += s.KCharpoly.Length - 1;
        int resDeg = aBlk.GetLength(0) - atDeg;
        var chi = PencilCharpolyBivariate(aBlk, cBlk);
        var fRes = BivariateDivideExact(chi, BivariateAtFactor(sectors));
        WeightCrossCheckResidual(fRes, n, wKet, wBra, rOdd, sectors, resDeg);

        // §0: exact realness of F_res (a free assertion on the exact object; the mod-p read saw
        // both embeddings agree, this is the theorem-grade form).
        bool residualIsReal = true;
        foreach (var qs in fRes)
        {
            foreach (var c in qs)
                if (!c.Im.IsZero) { residualIsReal = false; break; }
            if (!residualIsReal) break;
        }

        // §4 of the design: q-evenness of F_res, gated exactly (every odd-q coefficient zero).
        // The derivation it backs: the bipartite gauge T = diag((−1)^(σ(a)+σ(b))) satisfies
        // T·L(q)·T = L(−q) entry-wise (hops flip the bipartite parity, so the q-linear entries
        // change sign; the dissipator diagonal is q-free), and it acts within each R-sector
        // whenever (p+q)(N−1) is even (the census's commutator character), so the sector
        // charpoly is even in q; the AT strands pair a+bq ↔ a−bq under T, hence their divisor
        // is even too, hence F_res = charpoly/AT is. Even F_res makes G even in its
        // q-coefficients and disc_M(G) even in q, which halves the Route-A evaluations
        // (the negative node half comes free). Where the character is odd (the N=4 (1,2)
        // control: (p+q)(N−1) = 9) the gauge swaps the sectors and nothing is derived; the
        // field then records the measurement.
        bool evenInQ = true;
        foreach (var qs in fRes)
        {
            for (int i = 1; i < qs.Length; i += 2)
                if (!qs[i].Re.IsZero || !qs[i].Im.IsZero) { evenInQ = false; break; }
            if (!evenInQ) break;
        }

        int lambda0 = -2 * n;
        bool even = ExactOddTaylorAllZero(fRes, lambda0);

        // §1: the exact Taylor shift about Λ₀ materialised; f = H[0] = F_res(Λ₀, q).
        var shifted = TaylorShiftBivariate(fRes, lambda0);
        var f = GaussianPolynomial.Trim(shifted[0]);
        int fDeg = GaussianPolynomial.Degree(f);
        int vQf = -1;
        for (int i = 0; i < f.Length; i++)
            if (!f[i].Re.IsZero || !f[i].Im.IsZero) { vQf = i; break; }

        // The Λ-side support parity of a FOLD-FIXED block (ket OR bra weight = N/2 at even N;
        // PROOF_CODIM1_BY_ADDITIVITY §7(b)), read off the SAME Taylor shift at no extra cost.
        // DERIVATION, and it is a determinant identity, not an antiunitary one: the BARE fold leg
        // is the LINEAR, q-REVERSING rearrangement P·S(q)·P = 2Λ₀ − S(−q) (the registry's
        // P·L(q)·P = −2N − L(−q) in the uncleared scale; P commutes with the site reflection, so
        // it stays inside the R-sector). Taking determinants,
        //     χ_q(Λ₀ + u) = (−1)^d · χ_{−q}(Λ₀ − u),   d = deg_Λ,
        // at EVERY q, real or complex, over ℤ[i]. Writing Σ c_{j,k}·u^j·q^k that forces
        // c_{j,k} = 0 unless j + k ≡ d (mod 2): a CHECKERBOARD support. The reality half is the
        // separate, unconditional q-side antiparity (c real for even k, imaginary for odd k), so
        // on the surviving cells even/even are real and odd/odd purely imaginary.
        // ε IS FORCED, not measured: ε = (−1)^d, and since the residual is monic the cell (d, 0)
        // decides the branch before any other coefficient is looked at. Both branches are still
        // evaluated, as a self-check that they are mutually exclusive, and the verdict is
        // cross-gated against (−1)^resDeg by the caller.
        // WHY BOTH χ AND F_res: on the full charpoly the identity above is a theorem with no
        // premise; on F_res it additionally needs the AT factor to carry the same parity (this
        // arc's ungated AT-swap premise). Reading both turns that premise into a measurement.
        var cbF = CheckerboardRead(TaylorShiftBivariate(fRes, lambda0));
        var cbFull = CheckerboardRead(TaylorShiftBivariate(chi, lambda0));
        int cbSign = cbF.Sign, cbOddCells = cbF.OddCells, cbEvenCells = cbF.EvenCells;
        int cbFullSign = cbFull.Sign;

        var (isSquare, uMonic) = MonicSquareRoot(f);
        int uDeg = isSquare ? QcDeg(uMonic) : -1;
        bool uIsReal = isSquare;
        if (isSquare)
            foreach (var c in uMonic)
                if (!c.Im.IsZero) { uIsReal = false; break; }

        int mDeg = -1;
        bool compId = false;
        double msPerNode = 0.0;
        if (even)
        {
            if ((resDeg & 1) != 0)
                throw new InvalidOperationException("even residual with odd Λ-degree.");
            mDeg = resDeg / 2;
            var g = new GaussianInteger[mDeg + 1][];
            for (int j = 0; j <= mDeg; j++) g[j] = GaussianPolynomial.Trim(shifted[2 * j]);
            if (GaussianPolynomial.Degree(g[mDeg]) != 0 || !g[mDeg][0].Equals(GaussianInteger.One))
                throw new InvalidOperationException("G is not monic in M (the reindex is wrong).");

            // §2a: the composition identity disc_Λ(F_res) = (−4)^m·f·disc_M(G)² EXACT at fresh
            // integer nodes, sign included ((−1)^m: m = 21 R-even is odd, a genuine minus).
            compId = true;
            foreach (int q0 in new[] { 2, 3 })
            {
                var d0 = GaussianPolynomial.Discriminant(EvalBivariateAtQ(fRes, q0));
                var dm = GaussianPolynomial.Discriminant(EvalBivariateAtQ(g, q0));
                var f0 = GaussianPolynomial.Evaluate(f, new GaussianInteger(q0, 0));
                var rhs = f0 * dm * dm;
                var minusFour = new GaussianInteger(-4, 0);
                for (int e = 0; e < mDeg; e++) rhs *= minusFour;
                if (!d0.Equals(rhs)) { compId = false; break; }
            }

            // §2b: the Route-A per-node cost, timed at REPRESENTATIVE nodes (q0 ≈ 700, the top
            // half of the ~769-node grid, where the evaluated coefficients carry their full bit
            // length; small nodes would flatter the extrapolation). Recorded, not gated to a
            // number: quiet-machine cost is re-measured before it is quoted anywhere.
            if (timingNodes > 0)
            {
                var sw = System.Diagnostics.Stopwatch.StartNew();
                for (int t = 0; t < timingNodes; t++)
                {
                    var dm = GaussianPolynomial.Discriminant(EvalBivariateAtQ(g, 700 + t));
                    if (dm.Re.IsZero && dm.Im.IsZero)
                        throw new InvalidOperationException(
                            $"disc_M(G) vanished at the integer node q0={700 + t}.");
                }
                sw.Stop();
                msPerNode = sw.Elapsed.TotalMilliseconds / timingNodes;
            }
        }

        return new WeightSturmScoutReport(n, wKet, wBra, rOdd, atDeg, resDeg, residualIsReal, even,
            evenInQ, mDeg, fDeg, vQf, isSquare, uDeg, uIsReal, compId, msPerNode,
            cbSign, cbOddCells, cbEvenCells, cbFullSign);
    }

    /// <summary>The verdict of <see cref="WeightDiscMExactLanding"/>: disc_M(G) landed EXACTLY
    /// over ℤ (Route A of the Sturm design, §2 + §5 step 3), with every §1 gate value. The
    /// coefficients (low-to-high in octic q) are carried for the downstream squarefree split and
    /// Sturm count; DiscVq is read off the landed object, and the layer degrees are the mod-p
    /// squarefree layers of the STRIPPED landed polynomial at the two committed primes (a floor
    /// read of the true split, which the exact squarefree step pins).</summary>
    public sealed record WeightDiscMLandingReport(
        int N, int WKet, int WBra, bool ROdd,
        int MDeg,
        int DiscDeg, int DiscVq,
        bool VerificationNodesPass,
        bool CompositionIdentityModBothPrimes,
        bool UGcdOne,
        int[] StrippedLayerDegreesFirstPrime,
        int[] StrippedLayerDegreesSecondPrime,
        BigInteger[] DiscCoefficients,
        BigInteger[] UCoefficients);

    /// <summary>Route A: land disc_M(G) exactly over ℤ for the (wKet, wBra) R-sector.
    /// <paramref name="certifiedDegQD"/> is the certified deg_q of D = disc_Λ(F_res) (the
    /// D-device's TrueDiscriminantDegree, 1572 / 2124 at (1,3)@N=6): the landing is certified
    /// RELATIVE to it via deg disc_M(G) = (deg_q D − deg f)/2, the design's first grade.</summary>
    public static WeightDiscMLandingReport WeightDiscMExactLanding(
        int n, int wKet, int wBra, bool rOdd, int certifiedDegQD)
    {
        var (aBlk, cBlk) = WeightSectorPencil(n, wKet, wBra, rOdd);
        AssertHoppingSymmetry(cBlk, null);
        var sectors = F89AtFactorReconstruction.ClearedAtSectorsFromPencil(aBlk, cBlk);
        var chi = PencilCharpolyBivariate(aBlk, cBlk);
        var fRes = BivariateDivideExact(chi, BivariateAtFactor(sectors));
        int resDeg = aBlk.GetLength(0);
        foreach (var s in sectors) resDeg -= s.KCharpoly.Length - 1;
        WeightCrossCheckResidual(fRes, n, wKet, wBra, rOdd, sectors, resDeg);
        return DiscMLandingCore(fRes, -2 * n, certifiedDegQD, n, wKet, wBra, rOdd);
    }

    /// <summary>The synthetic-control entry (design §6): the same landing core on a hand-built
    /// bivariate F(Λ, q) = G₀((Λ − Λ₀)², q), so a sympy-pinned reference can be compared
    /// coefficient for coefficient. The only control that can catch a G-reindex bug: the degree
    /// gates cannot see one.</summary>
    public static WeightDiscMLandingReport DiscMExactLandingFromBivariate(
        GaussianInteger[][] fRes, int lambda0, int certifiedDegQD)
        => DiscMLandingCore(fRes, lambda0, certifiedDegQD, 0, 0, 0, false);

    private static WeightDiscMLandingReport DiscMLandingCore(
        GaussianInteger[][] fRes, int lambda0, int certifiedDegQD,
        int n, int wKet, int wBra, bool rOdd)
    {
        // The three entry gates, all exact on the exact F_res: realness (§0), q-evenness (the
        // derived bipartite-gauge parity; it is what LICENSES the mirrored grid below, so it is
        // a hard gate, not a report field), and the (Λ−Λ₀)² evenness (the certified read).
        foreach (var qs in fRes)
        {
            for (int i = 0; i < qs.Length; i++)
            {
                if (!qs[i].Im.IsZero)
                    throw new InvalidOperationException("F_res is not real: the ℤ landing does not apply (§3 contingency).");
                if ((i & 1) == 1 && !qs[i].Re.IsZero)
                    throw new InvalidOperationException("F_res is not even in q: the mirrored grid is not licensed.");
            }
        }
        if (!ExactOddTaylorAllZero(fRes, lambda0))
            throw new InvalidOperationException("F_res is not even in (Λ−Λ₀): no G exists.");

        var shifted = TaylorShiftBivariate(fRes, lambda0);
        var f = GaussianPolynomial.Trim(shifted[0]);
        int fDeg = GaussianPolynomial.Degree(f);
        if (fDeg < 0 || f[0].Re.IsZero)
            throw new InvalidOperationException("v_q(f) ≠ 0: the deg/v chaining of the composition identity fails.");

        int resDeg = fRes.Length - 1;
        int mDeg = resDeg / 2;
        var g = new GaussianInteger[mDeg + 1][];
        for (int j = 0; j <= mDeg; j++) g[j] = GaussianPolynomial.Trim(shifted[2 * j]);
        if (GaussianPolynomial.Degree(g[mDeg]) != 0 || !g[mDeg][0].Equals(GaussianInteger.One))
            throw new InvalidOperationException("G is not monic in M (the reindex is wrong).");

        if ((certifiedDegQD - fDeg) % 2 != 0)
            throw new InvalidOperationException("certified deg_q D − deg f is odd: the composition identity cannot hold.");
        int discDeg = (certifiedDegQD - fDeg) / 2;
        if ((discDeg & 1) != 0)
            throw new InvalidOperationException("deg disc_M(G) is odd yet disc_M(G) must be even in q.");
        int k = discDeg / 2;

        // Node values at q0 = 0..k; the negative half is the SAME values mirrored (disc_M(G) is
        // even in q because G's coefficients are, gated above). Realness of each node value is
        // asserted: G is real, so a nonzero Im is a pipeline bug.
        int mm = 2 * k;
        var vals = new BigInteger[mm + 1];
        for (int q0 = 0; q0 <= k; q0++)
        {
            var d = GaussianPolynomial.Discriminant(EvalBivariateAtQ(g, q0));
            if (!d.Im.IsZero)
                throw new InvalidOperationException($"disc_M(G)({q0}) is not real.");
            vals[k + q0] = d.Re;
            vals[k - q0] = d.Re;
        }

        // Integer forward differences over the unit-spaced grid t = 0..2k (t = q + k), then the
        // Newton form P(t) = Σ_j Δ^j·C(t, j) accumulated in the (2k)!-scaled falling-factorial
        // basis: A_m = c_m, A_j = A_{j+1}·(t − j) + ((2k)!/j!)·c_j, so A_0 = (2k)!·P with every
        // intermediate an integer (disc_M(G) is integer-valued on ℤ, so the differences are
        // integers and the final division by (2k)! must be exact; a remainder is a bug).
        var c = (BigInteger[])vals.Clone();
        for (int j = 1; j <= mm; j++)
            for (int i = mm; i >= j; i--)
                c[i] -= c[i - 1];
        var acc = new BigInteger[mm + 1];
        int accLen = 1;
        acc[0] = c[mm];
        BigInteger scale = BigInteger.One;
        for (int j = mm - 1; j >= 0; j--)
        {
            scale *= j + 1;
            for (int i = accLen; i >= 1; i--)
                acc[i] = acc[i - 1] - j * acc[i];
            acc[0] = -j * acc[0] + scale * c[j];
            accLen++;
        }
        var pt = new BigInteger[mm + 1];
        for (int i = 0; i <= mm; i++)
        {
            pt[i] = BigInteger.DivRem(acc[i], scale, out var rem);
            if (!rem.IsZero)
                throw new InvalidOperationException($"interpolant coefficient t^{i} is not integral.");
        }

        // Taylor shift back to q (t = q + k): in-place Ruffini-Horner, then the landed gates.
        for (int i = 0; i < mm; i++)
            for (int j = mm - 1; j >= i; j--)
                pt[j] += k * pt[j + 1];

        int landedDeg = -1;
        for (int i = mm; i >= 0; i--)
            if (!pt[i].IsZero) { landedDeg = i; break; }
        if (landedDeg != discDeg)
            throw new InvalidOperationException($"landed deg disc_M(G) = {landedDeg}, certified-backed {discDeg}.");
        int vq = -1;
        for (int i = 0; i <= landedDeg; i++)
            if (!pt[i].IsZero) { vq = i; break; }
        for (int i = 1; i <= landedDeg; i += 2)
            if (!pt[i].IsZero)
                throw new InvalidOperationException($"landed disc_M(G) has a nonzero odd coefficient q^{i}.");
        var disc = new BigInteger[landedDeg + 1];
        Array.Copy(pt, disc, landedDeg + 1);

        // Verification nodes OFF the grid, one negative: the only end-to-end check of the
        // mirroring (the on-grid evenness is construction; a fresh q0 < 0 is data).
        bool verify = true;
        foreach (int q0 in new[] { k + 1, -(k + 1), k + 2 })
        {
            var d = GaussianPolynomial.Discriminant(EvalBivariateAtQ(g, q0));
            if (!d.Im.IsZero)
                throw new InvalidOperationException($"disc_M(G)({q0}) is not real.");
            BigInteger hv = BigInteger.Zero;
            for (int i = landedDeg; i >= 0; i--) hv = hv * q0 + disc[i];
            if (hv != d.Re) { verify = false; break; }
        }

        // The composition identity at the POLYNOMIAL level mod both committed primes, sign
        // included: D ≡ (−4)^m·f·disc_M(G)² with D sampled and interpolated mod p from the
        // exact F_res (the design's machinery map names this the gap the degree-level test
        // leaves open: right degrees cannot certify the right polynomial).
        bool compOk = true;
        {
            long candidate = (1L << 30) + 1;
            while (candidate % 4 != 1) candidate += 2;
            int seen = 0;
            for (; seen < 2; candidate += 4)
            {
                if (!IsPrime(candidate)) continue;
                int p = checked((int)candidate);
                int? root = SqrtMinusOneFast(p);
                if (root is null) continue;
                seen++;

                var fResP = ReduceBivariate(fRes, p, root.Value);
                int samples = certifiedDegQD + 1 + 8;
                var dSamp = new int[samples];
                for (int q0 = 0; q0 < samples; q0++)
                    dSamp[q0] = DiscriminantModP(EvalBivariateModP(fResP, q0, p), p);
                var dNodes = new int[certifiedDegQD + 1];
                Array.Copy(dSamp, dNodes, dNodes.Length);
                var dp = InterpolateAtIntegerNodes(dNodes, p);
                for (int q0 = dNodes.Length; q0 < samples; q0++)
                    if (EvalModP(dp, q0, p) != dSamp[q0])
                        throw new InvalidOperationException($"D interpolant fails verification at q0={q0} (p={p}).");

                var discP = new int[landedDeg + 1];
                for (int i = 0; i <= landedDeg; i++) discP[i] = BigIntMod(disc[i], p);
                var fP = new int[fDeg + 1];
                for (int i = 0; i <= fDeg; i++) fP[i] = BigIntMod(f[i].Re, p);
                var rhs = MulModP(MulModP(discP, discP, p), fP, p);
                long sign = PowModP(p - 4, mDeg, p);
                for (int i = 0; i < rhs.Length; i++) rhs[i] = (int)(rhs[i] * sign % p);
                int dDeg = DegP(dp), rDeg = DegP(rhs);
                if (dDeg != rDeg) { compOk = false; }
                else
                    for (int i = 0; i <= dDeg; i++)
                        if (dp[i] != rhs[i]) { compOk = false; break; }
            }
        }

        // gcd(u, disc_M G) = 1 EXACT over ℚ (the u-gate promoted from mod-p): Euclid on the
        // monic square root ũ of f against the landed polynomial.
        var (isSquare, uMonic) = MonicSquareRoot(f);
        bool uGcdOne = false;
        var uInt = Array.Empty<BigInteger>();
        if (isSquare)
        {
            var pQc = new Qc[landedDeg + 1];
            for (int i = 0; i <= landedDeg; i++)
                pQc[i] = new Qc(new BigRational(disc[i]), new BigRational(0));
            var a = uMonic;
            var b = QcRem(pQc, uMonic);
            while (QcDeg(b) >= 0)
            {
                var t = QcRem(a, b);
                a = b;
                b = t;
            }
            uGcdOne = QcDeg(a) == 0;

            // u as a primitive integer polynomial (the unique integer multiple of ũ up to sign),
            // for the downstream real-μ=0-diabolic count.
            BigInteger den = BigInteger.One;
            for (int i = 0; i <= QcDeg(uMonic); i++)
            {
                if (!uMonic[i].Im.IsZero)
                    throw new InvalidOperationException("ũ is not real in the real world.");
                den = BigRational.Lcm(den, uMonic[i].Re.Denominator);
            }
            uInt = new BigInteger[QcDeg(uMonic) + 1];
            for (int i = 0; i < uInt.Length; i++)
                uInt[i] = uMonic[i].Re.Numerator * (den / uMonic[i].Re.Denominator);
            uInt = PrimitiveZ(uInt);
        }

        // The mod-p squarefree layers of the STRIPPED landed polynomial at both committed
        // primes (a floor read of the true split; the exact squarefree step pins it).
        var layerDegs = new int[2][];
        {
            var stripped = new BigInteger[landedDeg - vq + 1];
            Array.Copy(disc, vq, stripped, 0, stripped.Length);
            long candidate = (1L << 30) + 1;
            while (candidate % 4 != 1) candidate += 2;
            int seen = 0;
            for (; seen < 2; candidate += 4)
            {
                if (!IsPrime(candidate)) continue;
                int p = checked((int)candidate);
                if (SqrtMinusOneFast(p) is null) continue;
                var sp = new int[stripped.Length];
                for (int i = 0; i < sp.Length; i++) sp[i] = BigIntMod(stripped[i], p);
                var layers = FpSquarefreeLayers(sp, p);
                layerDegs[seen] = new int[layers.Count];
                for (int i = 0; i < layers.Count; i++) layerDegs[seen][i] = DegP(layers[i]);
                seen++;
            }
        }

        return new WeightDiscMLandingReport(n, wKet, wBra, rOdd, mDeg, landedDeg, vq,
            verify, compOk, uGcdOne, layerDegs[0], layerDegs[1], disc, uInt);
    }

    /// <summary>The verdict of <see cref="SturmCountFromLanding"/>: the EXACT squarefree split
    /// stripped disc_M(G) = A·B² over ℤ (A the codim-1 defective species layer, B the doubled
    /// off-line content) with the split proven by the reconstruction identity P = A·B·B and the
    /// pairwise gcd/squarefree gates (so DegA/DegB are the TRUE layer degrees, replacing the
    /// mod-p floors; the merge parameter k of the repaired ceiling claims is read off), and the
    /// Sturm counts in BOTH census conventions (q &gt; 0 and all-real; conversion exactly 2 by
    /// the evenness of A, gated). PinCountsA holds one Sturm count per census pin interval of
    /// this parity, each expected exactly 1; CountU is the real-μ=0-diabolic corollary.</summary>
    public sealed record WeightSturmCountReport(
        int N, int WKet, int WBra, bool ROdd,
        int DegA, int DegB,
        bool SplitIdentityOk, bool AEvenInQ,
        int CountAPositive, int CountAAllReal,
        int CountBPositive, int CountBAllReal,
        int CountUPositive, int CountUAllReal,
        int[] PinCountsA,
        BigInteger[] ACoefficients, BigInteger[] BCoefficients);

    /// <summary>Design §3 + §4 on a landed disc_M(G): primitive-PRS squarefree split and the
    /// exact Sturm count. <paramref name="censusPins"/> are the committed census loci of THIS
    /// parity (octic q, floats from the experiment table); each is bracketed by the rational
    /// window ±1/4000 and Sturm-counted on A.</summary>
    public static WeightSturmCountReport SturmCountFromLanding(
        WeightDiscMLandingReport landing, double[] censusPins)
    {
        var stripped = new BigInteger[landing.DiscDeg - landing.DiscVq + 1];
        Array.Copy(landing.DiscCoefficients, landing.DiscVq, stripped, 0, stripped.Length);
        var (a, b, splitOk) = SquarefreeSplitAB2(stripped);

        bool aEven = true;
        for (int i = 1; i < a.Length; i += 2)
            if (!a[i].IsZero) { aEven = false; break; }

        var chainA = SturmChainZ(a);
        var chainB = SturmChainZ(b);
        int aPos = SturmCountPositive(chainA), aAll = SturmCountAllReal(chainA);
        int bPos = SturmCountPositive(chainB), bAll = SturmCountAllReal(chainB);
        int uPos = 0, uAll = 0;
        if (landing.UCoefficients.Length > 0)
        {
            var chainU = SturmChainZ(landing.UCoefficients);
            uPos = SturmCountPositive(chainU);
            uAll = SturmCountAllReal(chainU);
        }

        var pins = new int[censusPins.Length];
        for (int i = 0; i < censusPins.Length; i++)
        {
            var mid = new BigRational(new BigInteger(Math.Round(censusPins[i] * 1e7)), new BigInteger(10_000_000));
            // Half-width 1/4000: the R-odd close pair 0.4627729 / 0.4635508 sits 7.8·10⁻⁴
            // apart, so ±1/1000 brackets would overlap and each read 2.
            var half = new BigRational(BigInteger.One, new BigInteger(4000));
            pins[i] = SturmCountInterval(chainA, mid - half, mid + half);
        }

        return new WeightSturmCountReport(landing.N, landing.WKet, landing.WBra, landing.ROdd,
            DegZ(a), DegZ(b), splitOk, aEven, aPos, aAll, bPos, bAll, uPos, uAll, pins, a, b);
    }

    /// <summary>The N=4 (1,2) R-even control of the design's §6, the committed story end to end:
    /// D = disc_Λ(F_res) landed exactly over ℤ by full-grid integer Newton (NO evenness assumed;
    /// the a-priori node bound resDeg·(resDeg−1) makes the landing standalone), with the census's
    /// float sector-disc-reality upgraded to the theorem grade on the way (every node value has
    /// Im exactly 0, and the node count exceeds the degree bound, so Im D ≡ 0 identically). Then
    /// the same split + Sturm machinery: the simple layer must carry EXACTLY the four committed
    /// √-branch loci on q &gt; 0 and the double layer the diabolic 0.659 (which must NOT be
    /// counted in A). Closes the recorded gap "no exact Sturm run on P₁₀ exists".</summary>
    public sealed record WeightSturmN4ControlReport(
        int DegD, int VqD,
        int DegA, int DegB,
        bool SplitIdentityOk,
        int CountAPositive, int CountBPositive,
        int[] PinCountsA, int PinCountB659);

    public static WeightSturmN4ControlReport WeightSturmN4Control()
    {
        var (aBlk, cBlk) = WeightSectorPencil(4, 1, 2, rOdd: false);
        AssertHoppingSymmetry(cBlk, null);
        var sectors = F89AtFactorReconstruction.ClearedAtSectorsFromPencil(aBlk, cBlk);
        var chi = PencilCharpolyBivariate(aBlk, cBlk);
        var fRes = BivariateDivideExact(chi, BivariateAtFactor(sectors));
        int resDeg = fRes.Length - 1;
        int dBound = resDeg * (resDeg - 1);

        var vals = new BigInteger[dBound + 1];
        for (int q0 = 0; q0 <= dBound; q0++)
        {
            var d = GaussianPolynomial.Discriminant(EvalBivariateAtQ(fRes, q0));
            if (!d.Im.IsZero)
                throw new InvalidOperationException(
                    $"disc_Λ(F_res)({q0}) is not real at the N=4 control: the census reality read fails exactly.");
            vals[q0] = d.Re;
        }
        var dPoly = IntegerNewtonFromUnitGrid(vals, shiftBack: 0);
        int degD = DegZ(dPoly);
        int vq = 0;
        while (dPoly[vq].IsZero) vq++;
        foreach (int q0 in new[] { dBound + 1, dBound + 2 })
        {
            var d = GaussianPolynomial.Discriminant(EvalBivariateAtQ(fRes, q0));
            BigInteger hv = BigInteger.Zero;
            for (int i = degD; i >= 0; i--) hv = hv * q0 + dPoly[i];
            if (!d.Im.IsZero || hv != d.Re)
                throw new InvalidOperationException($"N=4 D interpolant fails verification at q0={q0}.");
        }

        var stripped = new BigInteger[degD - vq + 1];
        Array.Copy(dPoly, vq, stripped, 0, stripped.Length);
        var (a, b, splitOk) = SquarefreeSplitAB2(stripped);
        var chainA = SturmChainZ(a);
        var chainB = SturmChainZ(b);

        var committed = new[] { 0.460212, 0.854438, 0.857458, 1.738181 };
        var pins = new int[committed.Length];
        for (int i = 0; i < committed.Length; i++)
        {
            var mid = new BigRational(new BigInteger(Math.Round(committed[i] * 1e6)), new BigInteger(1_000_000));
            var half = new BigRational(BigInteger.One, new BigInteger(1000));
            pins[i] = SturmCountInterval(chainA, mid - half, mid + half);
        }
        var mid659 = new BigRational(new BigInteger(659), new BigInteger(1000));
        var half659 = new BigRational(BigInteger.One, new BigInteger(100));
        int pinB = SturmCountInterval(chainB, mid659 - half659, mid659 + half659);

        return new WeightSturmN4ControlReport(degD, vq, DegZ(a), DegZ(b), splitOk,
            SturmCountPositive(chainA), SturmCountPositive(chainB), pins, pinB);
    }

    /// <summary>The verdict of <see cref="WeightPscDevice"/>: the design's §7 leading option at
    /// the theorem grade, everything over ℤ. Where G has two double M-roots or a triple M-root,
    /// psc₀ AND psc₁ of (G, ∂G/∂M) vanish together; on disc roots psc₀ = 0 already, so the
    /// device is gcd(psc₁, stripped disc_M(G)), computed EXACTLY (psc₁ landed over ℤ by the
    /// same mirrored-grid integer Newton as the disc, the gcd against the landed stripped disc):
    /// deg 0 means NO q ≠ 0, real or complex, carries two double M-roots or a triple, which
    /// closes the mult-4 hideouts (i) codim-2 double-pair coincidences outright and (ii) 3×3
    /// Jordan cubic branch points off the μ = 0 line (a 3×3 Jordan is a triple Λ- hence triple
    /// M-collision; ON the line the u-gate closes it: a Λ = −2n multiplicity ≥ 3 needs M = 0
    /// double in G, i.e. a common root of u and disc_M(G), excluded by gcd(u, disc_M G) = 1;
    /// parity alone would still admit a J₃ ⊕ J₁ inside multiplicity 4).
    /// PinSpeciesSigns is the per-pin sign of the double M-root M* (−1 fold pair on the line,
    /// +1 conj pair on the axis), read EXACTLY by a Tarski query: the generalized Sturm chain of
    /// (A, A′·h) with h ≡ −s₁₀·psc₁ mod A gives Σ sign(h) over the roots of A in any interval,
    /// and at a root q* of A, sign(h) = sign(M*) (M* = −s₁₀/psc₁ there; psc₁(q*) ≠ 0 is the
    /// device's own coprimality, and s₁₀(q*) ≠ 0 because M* = 0 would make q* a u-root, which
    /// the landing's u-gate excludes). AxisTotal/FoldTotal is the same query over (0, ∞): the
    /// species split of the WHOLE simple layer, not just the pinned loci. On a RETURNED report
    /// DeviceCoprimeExact is always true, the non-coprime path throwing before the species read
    /// (whose validity needs it); the field is the record, the throw is the gate.</summary>
    public sealed record WeightPscDeviceReport(
        int N, int WKet, int WBra, bool ROdd,
        int MDeg,
        int PscDeg, int S10Deg,
        bool DeviceCoprimeExact,
        int[] PinSpeciesSigns,
        int AxisTotalPositive, int FoldTotalPositive);

    /// <summary>Design §7: the psc device in M (not Λ: the Λ-form is blind at (1,3), every locus
    /// being Klein-doubled; in M the doubling is divided out). Consumes the landed disc and the
    /// exact simple layer A (so it runs after <see cref="WeightDiscMExactLanding"/> and
    /// <see cref="SturmCountFromLanding"/> without repeating their cost) and rebuilds only the
    /// cheap exact G. psc₁ and s₁₀ are landed over ℤ from per-node Bareiss determinants of the
    /// (2m−3)×(2m−2) subresultant matrix on the mirrored even grid, verified at fresh nodes
    /// including a negative one.</summary>
    public static WeightPscDeviceReport WeightPscDevice(
        WeightDiscMLandingReport landing, WeightSturmCountReport count, double[] censusPins)
    {
        var (aBlk, cBlk) = WeightSectorPencil(landing.N, landing.WKet, landing.WBra, landing.ROdd);
        AssertHoppingSymmetry(cBlk, null);
        var sectors = F89AtFactorReconstruction.ClearedAtSectorsFromPencil(aBlk, cBlk);
        var chi = PencilCharpolyBivariate(aBlk, cBlk);
        var fRes = BivariateDivideExact(chi, BivariateAtFactor(sectors));
        foreach (var qs in fRes)
            foreach (var cf in qs)
                if (!cf.Im.IsZero)
                    throw new InvalidOperationException("F_res is not real: the psc device's ℤ form does not apply.");
        int lambda0 = -2 * landing.N;
        if (!ExactOddTaylorAllZero(fRes, lambda0))
            throw new InvalidOperationException("F_res is not even in (Λ−Λ₀): no G exists.");
        var shifted = TaylorShiftBivariate(fRes, lambda0);
        int resDeg = fRes.Length - 1;
        int mDeg = resDeg / 2;
        if (mDeg != landing.MDeg)
            throw new InvalidOperationException("G degree disagrees with the landing report.");
        var g = new GaussianInteger[mDeg + 1][];
        for (int j = 0; j <= mDeg; j++) g[j] = GaussianPolynomial.Trim(shifted[2 * j]);
        if (GaussianPolynomial.Degree(g[mDeg]) != 0 || !g[mDeg][0].Equals(GaussianInteger.One))
            throw new InvalidOperationException("G is not monic in M (the reindex is wrong).");

        // Exact psc₁ and s₁₀ over ℤ: both even in q (G's coefficients are), so the mirrored
        // even grid halves the Bareiss evaluations, exactly as the disc landing's grid did.
        int dMax = 0;
        for (int j = 0; j <= mDeg; j++) dMax = Math.Max(dMax, GaussianPolynomial.Degree(g[j]));
        int bound = (2 * mDeg - 3) * dMax;
        if ((bound & 1) != 0) bound++;
        int half = bound / 2;
        var valsP = new BigInteger[bound + 1];
        var valsS = new BigInteger[bound + 1];
        for (int q0 = 0; q0 <= half; q0++)
        {
            var (p1, s0) = Psc1PairExactAtNode(g, mDeg, q0);
            valsP[half + q0] = p1;
            valsP[half - q0] = p1;
            valsS[half + q0] = s0;
            valsS[half - q0] = s0;
        }
        var psc = TrimZ(IntegerNewtonFromUnitGrid(valsP, half));
        var s10 = TrimZ(IntegerNewtonFromUnitGrid(valsS, half));
        foreach (int q0 in new[] { half + 1, -(half + 1) })
        {
            var (p1, s0) = Psc1PairExactAtNode(g, mDeg, q0);
            if (EvalZ(psc, q0) != p1 || EvalZ(s10, q0) != s0)
                throw new InvalidOperationException($"psc₁/s₁₀ interpolant fails verification at q0={q0}.");
        }

        // The Tarski read is valid only on the full precondition stack; each is a hard gate, a
        // silent skip would return plausible species numbers with no meaning (review 2026-08-11).
        if (!landing.UGcdOne)
            throw new InvalidOperationException("u-gate failed: s₁₀ could vanish at a root of A.");
        if (!count.SplitIdentityOk)
            throw new InvalidOperationException("the squarefree split is unproven: Sturm–Tarski does not apply.");

        // The device, exact over ℤ: gcd(psc₁, stripped landed disc) must have degree 0. It is
        // ALSO the third precondition: psc₁ ≠ 0 at the roots of A.
        var stripped = new BigInteger[landing.DiscDeg - landing.DiscVq + 1];
        Array.Copy(landing.DiscCoefficients, landing.DiscVq, stripped, 0, stripped.Length);
        bool coprime = DegZ(GcdZ(psc, stripped)) == 0;
        if (!coprime)
            throw new InvalidOperationException(
                "gcd(psc₁, stripped disc) ≠ 1: a two-double/triple-M point exists and the species read is invalid.");

        // The Tarski chain: h ≡ −s₁₀·psc₁ mod A (a positive multiple; reductions and the
        // product each carry only positive scalars, and the final negation sets the sign).
        var a = count.ACoefficients;
        var s10m = RemPositiveMultipleZ(s10, a);
        var pscm = RemPositiveMultipleZ(psc, a);
        var h = RemPositiveMultipleZ(MulZ(s10m, pscm), a);
        if (DegZ(h) < 0)
            throw new InvalidOperationException(
                "A divides s₁₀·psc₁: the Tarski chain would collapse and count zeros silently.");
        for (int i = 0; i < h.Length; i++) h[i] = -h[i];
        var chain = GeneralizedSturmChainZ(a, PrimitiveZ(MulZ(DerivativeZ(a), h)));

        var pinSigns = new int[censusPins.Length];
        for (int i = 0; i < censusPins.Length; i++)
        {
            var mid = new BigRational(new BigInteger(Math.Round(censusPins[i] * 1e7)), new BigInteger(10_000_000));
            var hw = new BigRational(BigInteger.One, new BigInteger(4000));
            pinSigns[i] = SturmCountInterval(chain, mid - hw, mid + hw);
        }
        int tarski = SturmCountPositive(chain);      // Σ sign(M*) over ALL q > 0 roots of A
        int n = count.CountAPositive;
        if (((n + tarski) & 1) != 0)
            throw new InvalidOperationException("Tarski total and root count have different parity.");
        int axis = (n + tarski) / 2, fold = (n - tarski) / 2;

        return new WeightPscDeviceReport(landing.N, landing.WKet, landing.WBra, landing.ROdd,
            mDeg, DegZ(psc), DegZ(s10), coprime, pinSigns, axis, fold);
    }

    private static BigInteger EvalZ(BigInteger[] p, long x)
    {
        BigInteger acc = BigInteger.Zero;
        var bx = new BigInteger(x);
        for (int i = p.Length - 1; i >= 0; i--) acc = acc * bx + p[i];
        return acc;
    }

    /// <summary>psc₁ and s₁₀ of (G(·, q0), ∂G/∂M(·, q0)) as EXACT integers: fraction-free
    /// Bareiss on the (2m−3)×(2m−2) subresultant matrix (rows x^s·f for s = m−3..0 and x^t·f′
    /// for t = m−2..0 in the basis x^(2m−3)..x^0); after eliminating the first 2m−4 columns the
    /// last row holds the two bordered-minor determinants exactly (Bareiss invariant), with the
    /// row-swap sign tracked, so S₁ = psc₁·x + s₁₀ falls out canonically: the per-node values
    /// interpolate to the true polynomials (the mod-p Gaussian variant of an earlier draft
    /// returned node-dependent scalar multiples, which the interpolation cannot survive; caught
    /// by its own verification node).</summary>
    private static (BigInteger Psc1, BigInteger S10) Psc1PairExactAtNode(
        GaussianInteger[][] g, int mDeg, long q0)
    {
        var bq = new BigInteger(q0);
        var coefs = new BigInteger[mDeg + 1];
        for (int j = 0; j <= mDeg; j++)
        {
            BigInteger acc = BigInteger.Zero;
            int dj = GaussianPolynomial.Degree(g[j]);
            for (int i = dj; i >= 0; i--) acc = acc * bq + g[j][i].Re;
            coefs[j] = acc;
        }

        int m = mDeg;
        int rows = 2 * m - 3, cols = 2 * m - 2;
        var mat = new BigInteger[rows][];
        int r = 0;
        var df = new BigInteger[m];
        for (int i = 1; i <= m; i++) df[i - 1] = i * coefs[i];
        for (int s = m - 3; s >= 0; s--, r++)
        {
            mat[r] = new BigInteger[cols];
            for (int i = 0; i <= m; i++) mat[r][cols - 1 - (i + s)] = coefs[i];
        }
        for (int t = m - 2; t >= 0; t--, r++)
        {
            mat[r] = new BigInteger[cols];
            for (int i = 0; i <= m - 1; i++) mat[r][cols - 1 - (i + t)] = df[i];
        }

        int sign = 1;
        BigInteger prev = BigInteger.One;
        for (int c = 0; c < rows - 1; c++)
        {
            int piv = -1;
            for (int i = c; i < rows; i++)
                if (!mat[i][c].IsZero) { piv = i; break; }
            if (piv < 0) return (BigInteger.Zero, BigInteger.Zero);
            if (piv != c)
            {
                (mat[c], mat[piv]) = (mat[piv], mat[c]);
                sign = -sign;
            }
            for (int i = c + 1; i < rows; i++)
            {
                for (int j = c + 1; j < cols; j++)
                    mat[i][j] = (mat[c][c] * mat[i][j] - mat[i][c] * mat[c][j]) / prev;
                mat[i][c] = BigInteger.Zero;
            }
            prev = mat[c][c];
        }
        return (sign * mat[rows - 1][cols - 2], sign * mat[rows - 1][cols - 1]);
    }

    // ---- ℤ[x] arithmetic for the squarefree split and the Sturm chain (design §3 + §4) ----

    /// <summary>Exact integer Newton interpolation on the unit grid t = 0..m (values given in
    /// that order), returned in the variable q = t − shiftBack. Runs entirely over ℤ: the
    /// forward differences are integers, the Newton form is accumulated in the m!-scaled
    /// falling-factorial basis (A_m = c_m, A_j = A_{j+1}·(t − j) + (m!/j!)·c_j, so A_0 = m!·P),
    /// and the final division by m! is asserted exact.</summary>
    private static BigInteger[] IntegerNewtonFromUnitGrid(BigInteger[] values, int shiftBack)
    {
        int mm = values.Length - 1;
        var c = (BigInteger[])values.Clone();
        for (int j = 1; j <= mm; j++)
            for (int i = mm; i >= j; i--)
                c[i] -= c[i - 1];
        var acc = new BigInteger[mm + 1];
        int accLen = 1;
        acc[0] = c[mm];
        BigInteger scale = BigInteger.One;
        for (int j = mm - 1; j >= 0; j--)
        {
            scale *= j + 1;
            for (int i = accLen; i >= 1; i--)
                acc[i] = acc[i - 1] - j * acc[i];
            acc[0] = -j * acc[0] + scale * c[j];
            accLen++;
        }
        var pt = new BigInteger[mm + 1];
        for (int i = 0; i <= mm; i++)
        {
            pt[i] = BigInteger.DivRem(acc[i], scale, out var rem);
            if (!rem.IsZero)
                throw new InvalidOperationException($"interpolant coefficient t^{i} is not integral.");
        }
        if (shiftBack != 0)
            for (int i = 0; i < mm; i++)
                for (int j = mm - 1; j >= i; j--)
                    pt[j] += shiftBack * pt[j + 1];
        return pt;
    }

    private static int DegZ(BigInteger[] a)
    {
        for (int i = a.Length - 1; i >= 0; i--) if (!a[i].IsZero) return i;
        return -1;
    }

    private static BigInteger[] TrimZ(BigInteger[] a)
    {
        int d = DegZ(a);
        if (d < 0) return Array.Empty<BigInteger>();
        var r = new BigInteger[d + 1];
        Array.Copy(a, r, d + 1);
        return r;
    }

    private static BigInteger[] DerivativeZ(BigInteger[] a)
    {
        int d = DegZ(a);
        if (d <= 0) return Array.Empty<BigInteger>();
        var r = new BigInteger[d];
        for (int i = 1; i <= d; i++) r[i - 1] = i * a[i];
        return r;
    }

    /// <summary>Divide by the POSITIVE content: the sign pattern (which the Sturm chain lives
    /// on) is preserved exactly.</summary>
    private static BigInteger[] PrimitiveZ(BigInteger[] a)
    {
        int d = DegZ(a);
        if (d < 0) return Array.Empty<BigInteger>();
        BigInteger g = BigInteger.Zero;
        for (int i = 0; i <= d; i++) g = BigInteger.GreatestCommonDivisor(g, a[i]);
        var r = new BigInteger[d + 1];
        for (int i = 0; i <= d; i++) r[i] = a[i] / g;
        return r;
    }

    private static BigInteger[] MulZ(BigInteger[] a, BigInteger[] b)
    {
        int da = DegZ(a), db = DegZ(b);
        if (da < 0 || db < 0) return Array.Empty<BigInteger>();
        var r = new BigInteger[da + db + 1];
        for (int i = 0; i <= da; i++)
        {
            if (a[i].IsZero) continue;
            for (int j = 0; j <= db; j++)
                r[i + j] += a[i] * b[j];
        }
        return r;
    }

    /// <summary>Remainder of f by g up to a POSITIVE scalar multiple (repeated lc(g)-scaled
    /// elimination with the accumulated sign tracked and fixed at the end). A positive multiple
    /// is all a Sturm chain or a gcd needs, and keeping the multiplier's sign explicit is
    /// exactly where hand-rolled chains go wrong (design §4's warning).</summary>
    private static BigInteger[] RemPositiveMultipleZ(BigInteger[] f, BigInteger[] g)
    {
        int dg = DegZ(g);
        if (dg < 0) throw new DivideByZeroException("polynomial remainder by zero.");
        var lc = g[dg];
        var rem = TrimZ(f);
        bool negate = false;
        int dr = DegZ(rem);
        int steps = 0;
        while (dr >= dg)
        {
            var c = rem[dr];
            var next = new BigInteger[dr];
            for (int j = 0; j < dr; j++)
            {
                next[j] = rem[j] * lc;
                int gj = j - (dr - dg);
                if (gj >= 0 && gj < dg)
                    next[j] -= c * g[gj];
            }
            if (lc.Sign < 0) negate = !negate;
            rem = TrimZ(next);
            dr = DegZ(rem);
            // Periodic positive-content strip: a large degree gap (e.g. reducing the deg-2100
            // Tarski product mod the deg-108 A) multiplies by lc once per step, and without the
            // strip the coefficients grow multiplicatively; positive content preserves the
            // positive-multiple semantics exactly.
            if (++steps % 8 == 0 && dr >= dg) rem = PrimitiveZ(rem);
        }
        if (negate)
            for (int j = 0; j < rem.Length; j++) rem[j] = -rem[j];
        return rem;
    }

    /// <summary>Primitive-PRS gcd over ℤ (design §3: plain ℚ-Euclid is the blowup risk),
    /// normalised to a positive leading coefficient.</summary>
    private static BigInteger[] GcdZ(BigInteger[] f, BigInteger[] g)
    {
        var a = PrimitiveZ(TrimZ(f));
        var b = PrimitiveZ(TrimZ(g));
        if (DegZ(a) < DegZ(b)) (a, b) = (b, a);
        while (DegZ(b) > 0)
        {
            var r = RemPositiveMultipleZ(a, b);
            a = b;
            b = PrimitiveZ(r);
        }
        if (DegZ(b) == 0) a = new[] { BigInteger.One };
        if (DegZ(a) >= 0 && a[DegZ(a)].Sign < 0)
            for (int i = 0; i < a.Length; i++) a[i] = -a[i];
        return a;
    }

    /// <summary>Exact division in ℤ[x] (throws if not exact; when f = g·h with h ∈ ℤ[x], long
    /// division produces exactly h's integer coefficients).</summary>
    private static BigInteger[] ExactDivZ(BigInteger[] f, BigInteger[] g)
    {
        int df = DegZ(f), dg = DegZ(g);
        if (dg < 0) throw new DivideByZeroException("polynomial division by zero.");
        if (df < dg) throw new InvalidOperationException("exact ℤ[x] division with smaller dividend.");
        var rem = new BigInteger[df + 1];
        Array.Copy(f, rem, df + 1);
        var quo = new BigInteger[df - dg + 1];
        for (int i = df; i >= dg; i--)
        {
            var c = BigInteger.DivRem(rem[i], g[dg], out var rr);
            if (!rr.IsZero)
                throw new InvalidOperationException("expected-exact ℤ[x] division has a non-integer step.");
            quo[i - dg] = c;
            if (c.IsZero) continue;
            for (int j = 0; j <= dg; j++)
                rem[i - dg + j] -= c * g[j];
        }
        if (DegZ(rem) >= 0)
            throw new InvalidOperationException("expected-exact ℤ[x] division left a remainder.");
        return quo;
    }

    /// <summary>The exact squarefree split P = A·B² with A, B squarefree and coprime (the shape
    /// MaxDiscMultiplicity ≤ 4 forces on the stripped disc; a multiplicity-3 part throws in
    /// ExactDivZ rather than reaching the gates). Returns (A, B, ok) where ok bundles the reconstruction identity P = A·B·B checked
    /// coefficient for coefficient over ℤ (no-rounding case 1: an exact route exists) with the
    /// pairwise gcd and squarefree gates gcd(A,B) = gcd(A,A′) = gcd(B,B′) = 1.</summary>
    private static (BigInteger[] A, BigInteger[] B, bool Ok) SquarefreeSplitAB2(BigInteger[] p)
    {
        var pp = PrimitiveZ(TrimZ(p));
        var b = GcdZ(pp, DerivativeZ(pp));
        var ab = PrimitiveZ(ExactDivZ(pp, b));
        var a = PrimitiveZ(ExactDivZ(ab, b));

        // Reconstruction identity P = A·B·B over ℤ, sign included (A carries P's sign, B is
        // normalised positive by GcdZ; both quotients of primitives are primitive by Gauss, so
        // the PrimitiveZ wrappers are no-ops and the compare is exact, not up-to-content).
        var recon = MulZ(MulZ(a, b), b);
        bool ok = recon.Length == pp.Length;
        if (ok)
            for (int i = 0; i < recon.Length; i++)
                if (recon[i] != pp[i]) { ok = false; break; }
        ok = ok && DegZ(GcdZ(a, b)) == 0
                && (DegZ(a) <= 0 || DegZ(GcdZ(a, DerivativeZ(a))) == 0)
                && (DegZ(b) <= 0 || DegZ(GcdZ(b, DerivativeZ(b))) == 0);
        return (a, b, ok);
    }

    /// <summary>The Sturm chain of p with every element a POSITIVE-scalar multiple of the
    /// classical one (primitive parts of sign-fixed pseudo-remainders): the variation counts are
    /// those of the classical chain exactly.</summary>
    private static List<BigInteger[]> SturmChainZ(BigInteger[] p)
    {
        var p0 = PrimitiveZ(TrimZ(p));
        return GeneralizedSturmChainZ(p0, PrimitiveZ(DerivativeZ(p0)));
    }

    /// <summary>The Sylvester generalized chain (p0, p1, −rem, …): with p1 = p0′·h it computes
    /// Tarski queries, Σ sign(h) over the roots of p0 in an interval; with p1 = p0′ it is the
    /// classical Sturm chain.</summary>
    private static List<BigInteger[]> GeneralizedSturmChainZ(BigInteger[] p0, BigInteger[] p1)
    {
        var chain = new List<BigInteger[]> { p0 };
        if (DegZ(p1) < 0) return chain;
        chain.Add(p1);
        while (true)
        {
            var r = RemPositiveMultipleZ(chain[^2], chain[^1]);
            if (DegZ(r) < 0) break;
            var next = PrimitiveZ(r);
            for (int i = 0; i < next.Length; i++) next[i] = -next[i];
            chain.Add(next);
        }
        return chain;
    }

    private static int SignVariations(IEnumerable<int> signs)
    {
        int v = 0, prev = 0;
        foreach (var s in signs)
        {
            if (s == 0) continue;
            if (prev != 0 && s != prev) v++;
            prev = s;
        }
        return v;
    }

    private static int VariationsAtInfinity(List<BigInteger[]> chain, bool positive)
        => SignVariations(chain.Select(p =>
        {
            int d = DegZ(p);
            if (d < 0) return 0;
            int s = p[d].Sign;
            return positive || (d & 1) == 0 ? s : -s;
        }));

    private static int VariationsAtRational(List<BigInteger[]> chain, BigRational x)
        => SignVariations(chain.Select(p =>
        {
            int d = DegZ(p);
            if (d < 0) return 0;
            var acc = new BigRational(0);
            for (int i = d; i >= 0; i--) acc = acc * x + new BigRational(p[i]);
            return acc.Sign;
        }));

    /// <summary>Distinct real roots with q &gt; 0 (the census convention; 0 must not be a root
    /// of the chain's head, which the stripped inputs guarantee).</summary>
    private static int SturmCountPositive(List<BigInteger[]> chain)
        => VariationsAtRational(chain, new BigRational(0)) - VariationsAtInfinity(chain, positive: true);

    private static int SturmCountAllReal(List<BigInteger[]> chain)
        => VariationsAtInfinity(chain, positive: false) - VariationsAtInfinity(chain, positive: true);

    private static int SturmCountInterval(List<BigInteger[]> chain, BigRational lo, BigRational hi)
        => VariationsAtRational(chain, lo) - VariationsAtRational(chain, hi);

    /// <summary>The exact Taylor shift of the bivariate about Λ₀:
    /// H[m](q) = Σ_{k≥m} C(k,m)·Λ₀^{k−m}·f[k](q) over ℤ[i], the coefficient form
    /// <see cref="ExactOddTaylorAllZero"/> tests against, materialised (H[0] = f(Λ₀, q); the even
    /// rows H[2j] are the compression G with M = (Λ−Λ₀)² when the odd rows vanish).</summary>
    private static GaussianInteger[][] TaylorShiftBivariate(GaussianInteger[][] f, int lambda0)
    {
        int kMax = f.Length - 1;
        var binom = new BigInteger[kMax + 1][];
        for (int k = 0; k <= kMax; k++)
        {
            binom[k] = new BigInteger[k + 1];
            binom[k][0] = binom[k][k] = BigInteger.One;
            for (int m = 1; m < k; m++)
                binom[k][m] = binom[k - 1][m - 1] + binom[k - 1][m];
        }
        var lamPow = new BigInteger[kMax + 1];
        lamPow[0] = BigInteger.One;
        for (int e = 1; e <= kMax; e++) lamPow[e] = lamPow[e - 1] * lambda0;
        int len = 0;
        foreach (var qs in f) len = Math.Max(len, qs.Length);
        var h = new GaussianInteger[kMax + 1][];
        for (int m = 0; m <= kMax; m++)
        {
            var row = new GaussianInteger[len];
            for (int k = m; k <= kMax; k++)
            {
                var w = binom[k][m] * lamPow[k - m];
                if (w.IsZero) continue;
                for (int i = 0; i < f[k].Length; i++)
                    row[i] = row[i] + f[k][i] * w;
            }
            h[m] = row;
        }
        return h;
    }

    /// <summary>The exact monic square root: f = lc(f)·ũ² over ℚ(i) if it exists (the Sturm
    /// design's §1: normalise monic FIRST, because the recursion from the leading coefficient
    /// would need √lc, which need not exist in ℚ(i); c := lc(f) is then not an unknown). ũ is
    /// monic of degree deg f / 2, determined coefficient by coefficient from the top; the
    /// verification compares ALL coefficients of ũ² against f/lc exactly, so a squarefree f (the
    /// N=4 control) returns false rather than a wrong root. The load-bearing use: a true verdict
    /// upgrades the two-prime perfect-square read of the μ = 0 stratum to algebra.</summary>
    private static (bool IsSquare, Qc[] UMonic) MonicSquareRoot(GaussianInteger[] f)
    {
        int d = GaussianPolynomial.Degree(f);
        if (d < 0 || (d & 1) != 0) return (false, Array.Empty<Qc>());
        var fq = FromZi(f);
        var lc = fq[d];
        var ft = new Qc[d + 1];
        for (int i = 0; i <= d; i++) ft[i] = fq[i] / lc;
        int half = d / 2;
        var one = new Qc(new BigRational(1), new BigRational(0));
        var two = new Qc(new BigRational(2), new BigRational(0));
        var u = new Qc[half + 1];
        for (int i = 0; i < half; i++) u[i] = Qc.Zero;
        u[half] = one;
        for (int j = half - 1; j >= 0; j--)
        {
            // coefficient of x^{half+j} in ũ²: 2·u[j] (against the monic top) + the known middle.
            var acc = Qc.Zero;
            for (int a = j + 1; a <= half - 1; a++)
                acc += u[a] * u[half + j - a];
            u[j] = (ft[half + j] - acc) / two;
        }
        var sq = new Qc[d + 1];
        for (int i = 0; i <= d; i++) sq[i] = Qc.Zero;
        for (int a = 0; a <= half; a++)
            for (int b = 0; b <= half; b++)
                sq[a + b] += u[a] * u[b];
        for (int i = 0; i <= d; i++)
            if (!(sq[i] - ft[i]).IsZero) return (false, Array.Empty<Qc>());
        return (true, u);
    }

    /// <summary>The verdict of <see cref="WeightCertifyDiscMultiplicity"/>: the weight-general
    /// sibling of <see cref="DiscMultiplicityReport"/>, carrying the block weights.</summary>
    public sealed record WeightDiscMultiplicityReport(
        int N, int WKet, int WBra, bool ROdd, int SectorDimension, int AtDegree, int ResidualDegree,
        int InfinityRepeatedD, int DiscriminantDegreeBound,
        int TrueDiscriminantDegree, int TrueQValuationD, int LcDivisorBoundD,
        long LayerPrime, int[] DiscLayerDegrees, int MaxDiscMultiplicity,
        int PrimesSampled, bool DiscLayersCertified,
        int SamplesPerPrime, bool QEvenHalvingUsed, bool SharpeningSkipped);

    /// <summary>Is the (wKet, wBra) block fixed by a fold leg? Exactly when N is even and one weight
    /// side sits at N/2: f_P: (p, q̃) ↦ (p, N−q̃) fixes bra weight N/2, f_Q: (p, q̃) ↦ (N−p, q̃) fixes
    /// ket weight N/2 (PROOF_CODIM1_BY_ADDITIVITY §7, the fold lattice; EITHER leg suffices, the
    /// (2,1)@N=4 control having refuted the first draft's "bra only"). This is the premise of the
    /// checkerboard law, so it is the reason to EXPECT the q-even halving in
    /// <see cref="WeightCertifyDiscMultiplicity"/> to be licensed; it is not the licence itself,
    /// which is that method's own exact read (<see cref="QNegationSupportParityExact"/>), because
    /// the checkerboard's descent from the full charpoly to F_res is measured per block and not
    /// proven. A statement about the block's index pair alone, so it costs nothing to ask.</summary>
    public static bool IsFoldFixedBlock(int n, int wKet, int wBra)
        => n % 2 == 0 && (wKet == n / 2 || wBra == n / 2);

    /// <summary>The D-device's cost inputs WITHOUT the sampling loop (the E0 probe of the N=8
    /// design, 2026-08-12): everything <see cref="WeightCertifyDiscMultiplicity"/> computes
    /// BEFORE its prime loop — the degree bound dBound = resDeg(resDeg−1) − 2m_D, the lc-divisor
    /// prime cap, and the row-norm bit length the cap scales with. The device's total cost is
    /// (cap + 1) primes × one full <c>SamplesPerPrime</c>-node mod-p disc pass each, so this probe
    /// is what decides whether the device is affordable at a given (n, weights) before anything
    /// runs for days. <c>SamplesPerPrime</c> is dBound + 25 on the dense path and dBound/2 + 25
    /// where the q-even halving is licensed, and the probe takes the SAME exact licence read as
    /// the device (<c>QEvenHalving</c>), so a fold-fixed block is not priced at twice its cost.
    /// Deliberately duplicates the ~15 pre-loop lines of the device (same building blocks, same
    /// order) rather than refactoring the certified path; skips WeightCrossCheckResidual
    /// (validation of fRes, gated separately on the scout path).</summary>
    public static (int BlockDim, int AtDeg, int ResDeg, int MD, int DBound,
        int LcDivisorBoundD, long NormDBits, long RowFBits, double MsPerPrime,
        int SamplesPerPrime, bool QEvenHalving)
        WeightDiscDeviceBoundProbe(int n, int wKet, int wBra, bool rOdd, int timingPrimes = 0,
            bool skipLeadingFormSharpening = false)
    {
        var (aBlk, cBlk) = WeightSectorPencil(n, wKet, wBra, rOdd);
        AssertHoppingSymmetry(cBlk, null);
        var sectors = F89AtFactorReconstruction.ClearedAtSectorsFromPencil(aBlk, cBlk);
        int atDeg = 0;
        foreach (var s in sectors) atDeg += s.KCharpoly.Length - 1;
        int blockDim = aBlk.GetLength(0);
        int resDeg = blockDim - atDeg;

        var chi = PencilCharpolyBivariate(aBlk, cBlk);
        var fRes = BivariateDivideExact(chi, BivariateAtFactor(sectors));

        // mirrors the device's branch, and must: a probe that priced the sharpening the device
        // will skip would be the more expensive of the two at exactly the block it exists for
        var (fResLead, lcD) = ResidualLeadingFormFromPencil(cBlk, sectors);
        int mD = -1;
        int dBound = resDeg * (resDeg - 1);
        if (!skipLeadingFormSharpening)
        {
            (_, mD) = SquarefreeLayers(FromZi(fResLead));
            if ((mD == 0) != !lcD.Equals(GaussianInteger.Zero))
                throw new InvalidOperationException("m_D inconsistent with disc_x(f_res).");
            dBound -= 2 * mD;
        }

        BigInteger rowF = BigInteger.Zero;
        foreach (var coeff in fRes) rowF += HNormQ(coeff);
        BigInteger rowFp = resDeg * rowF;
        BigInteger normD = BigInteger.Pow(rowF, Math.Max(resDeg - 1, 0)) * BigInteger.Pow(rowFp, resDeg);
        int lcDivisorBoundD = (int)(2 * normD.GetBitLength() / 30) + 1;

        // the per-prime marginal (the timingNodes pattern of WeightSturmScout): run the
        // device's own per-prime sampling body for a few split primes, timed AFTER the
        // one-time pre-loop above, so the figure is the marginal and not the build
        // the probe must price the path the device will actually take, so it takes the same
        // licence read: on a block where D is even in q the device samples half as many points,
        // and a probe that priced the dense pass would overstate a fold-fixed block by ~2×
        bool qEven = QEvenSamplingLicence(fRes, n, wKet, wBra, forceDense: false);
        const int extra = 24;
        int samples = (qEven ? dBound / 2 : dBound) + 1 + extra;

        double msPerPrime = 0.0;
        if (timingPrimes > 0)
        {
            var clock = System.Diagnostics.Stopwatch.StartNew();
            int timed = 0;
            long candidate = (1L << 30) + 1;
            while (candidate % 4 != 1) candidate += 2;
            for (; timed < timingPrimes; candidate += 4)
            {
                if (!IsPrime(candidate)) continue;
                int p = checked((int)candidate);
                int? root = SqrtMinusOneFast(p);
                if (root is null) continue;
                var fResP = ReduceBivariate(fRes, p, root.Value);
                var dSamp = new int[samples];
                for (int q0 = 0; q0 < samples; q0++)
                {
                    var fL = EvalBivariateModP(fResP, q0, p);
                    if (DegP(fL) != resDeg || fL[resDeg] != 1)
                        throw new InvalidOperationException($"reduced residual not monic of degree {resDeg} at q0={q0}.");
                    dSamp[q0] = DiscriminantModP(fL, p);
                }
                int nodeCount = samples - extra;
                var dNodes = new int[nodeCount];
                Array.Copy(dSamp, dNodes, nodeCount);
                int[] dp;
                if (qEven)
                {
                    AssertSquareAbscissaeFit(p, samples - 1);
                    var abscissae = new int[nodeCount];
                    for (int i = 0; i < nodeCount; i++) abscissae[i] = (int)((long)i * i % p);
                    dp = SpreadEven(InterpolateAtNodes(abscissae, dNodes, p));
                }
                else dp = InterpolateAtIntegerNodes(dNodes, p);
                for (int q0 = nodeCount; q0 < samples; q0++)
                    if (EvalModP(dp, q0, p) != dSamp[q0])
                        throw new InvalidOperationException($"D interpolant fails verification at q0={q0} (p={p}).");
                timed++;
            }
            msPerPrime = clock.ElapsedMilliseconds / (double)timingPrimes;
        }

        return (blockDim, atDeg, resDeg, mD, dBound, lcDivisorBoundD,
                normD.GetBitLength(), rowF.GetBitLength(), msPerPrime, samples, qEven);
    }

    /// <summary>The certificate-grade D-device on the weight-general path: the exact sibling of
    /// <see cref="CertifyDiscMultiplicity"/> for any all-2-cycle (wKet, wBra) R-sector (see
    /// <see cref="WeightSectorPencil"/> for the scope). Same proof shape, nothing weakened: the
    /// q → ∞ leading form fixes m_D and the PROVEN degree bound resDeg·(resDeg−1) − 2·m_D (the
    /// scout methods use the crude bound, and skipLeadingFormSharpening takes it here too);
    /// the lc-divisor bound (Hadamard on the Sylvester matrix
    /// of (F_res, F_res′)) caps how many primes can lie about deg_q D or v_q(D); primes are
    /// sampled past that cap, deg_q D is certified as the max and v_q(D) as the min over the
    /// sampled primes, and the layer profile is read at a prime attaining BOTH (fail closed,
    /// <c>DiscLayersCertified = false</c>, if none does). Per-point cross-check at q0 = 2, 3
    /// (<see cref="WeightCrossCheckResidual"/>) guards the bivariate on every new block shape.
    /// For (1,3)@N=6 this upgrades the scout's two-prime layer agreement to certified
    /// TrueDiscriminantDegree / TrueQValuationD, the invariants the arc's Sturm closing step
    /// will lean on exactly as the (1,2) gcd certificate leaned on the (1,2) D-device's.
    /// <para>Since 2026-08-13 the sampling halves where D is even in q (the checkerboard
    /// corollary): see <see cref="QNegationSupportParityExact"/> for the licence, which is an
    /// exact ℤ[i] read and not the block's index pair. Nothing about the CERTIFICATE changes:
    /// same bound, same cap, same primes, same certified invariants, and the N=4 A/B gate holds
    /// the two paths to bit-identical reports; only the number of sample points per prime does.
    /// The interpolant is spread back to D before anything downstream reads it, so valuation,
    /// q-strip and layer profile stay on the object they always read. WHERE THE SOUNDNESS SITS,
    /// because it moved: on the halved path it rests wholly on the exact read plus the
    /// pre-existing proven dBound (D even and deg D ≤ dBound give deg E ≤ dBound/2, so dBound/2+1
    /// values determine E), and NOT on the 24 verification nodes. Those over-determine on the
    /// dense path and under-determine on the halved one, where a wrongly-assumed evenness leaves
    /// dBound/2 − 1 free roots against 24 test points; they are a per-prime smoke test for the
    /// interpolation code, never a falsifier of the premise. Two things also change on calls that
    /// do not end up halving: the pre-loop pays one Taylor shift of F_res for the licence read
    /// (skipped entirely when forceDenseSampling short-circuits it), and a fold-fixed block
    /// failing that read now throws where it previously certified.</para></summary>
    /// <param name="forceDenseSampling">Take the dense q0 = 0…dBound path even where the exact
    /// even-in-q read licenses the halved one. The A/B control of
    /// <c>Disc13CertificateTests.N4_QEvenHalving_ReproducesTheDenseCertificateExactly</c>, and the
    /// way out if a fold-fixed block ever fails the support-parity read (which otherwise throws,
    /// that failure being a finding about the checkerboard's descent to F_res and not a nuisance).
    /// That second purpose is UNEXERCISED: no block is known that reaches the throw, so the escape
    /// is reasoned, not gated.</param>
    /// <param name="skipLeadingFormSharpening">Use the crude Bauer-Fike bound resDeg·(resDeg−1)
    /// instead of paying m_D for the 2·m_D sharpening. Under 1% more sample points once resDeg
    /// reaches the hundreds, against a Qc SquarefreeLayers call measured at ~8 min for degree 80.
    /// The leading form itself is still built, so its AT-divisibility cross-check still runs; only
    /// the sharpening is skipped. The certified invariants are unmoved: DiscriminantDegreeBound
    /// rises and InfinityRepeatedD becomes −1.</param>
    public static WeightDiscMultiplicityReport WeightCertifyDiscMultiplicity(
        int n, int wKet, int wBra, bool rOdd, int maxPrimes = 20000, Action<string>? log = null,
        bool forceDenseSampling = false, bool skipLeadingFormSharpening = false)
    {
        var (aBlk, cBlk) = WeightSectorPencil(n, wKet, wBra, rOdd);
        AssertHoppingSymmetry(cBlk, null);      // uniform metric: all orbits are 2-cycles
        var sectors = F89AtFactorReconstruction.ClearedAtSectorsFromPencil(aBlk, cBlk);
        int atDeg = 0;
        foreach (var s in sectors) atDeg += s.KCharpoly.Length - 1;
        int blockDim = aBlk.GetLength(0);
        int resDeg = blockDim - atDeg;

        var chi = PencilCharpolyBivariate(aBlk, cBlk);
        var fRes = BivariateDivideExact(chi, BivariateAtFactor(sectors));
        WeightCrossCheckResidual(fRes, n, wKet, wBra, rOdd, sectors, resDeg);

        // m_D buys the SHARPENING 2·m_D off the crude Bauer-Fike bound resDeg(resDeg−1), and at a
        // large block it is priced in the wrong direction: the Qc SquarefreeLayers call behind it
        // was measured at ~8 min for degree 80 (a gcd chain over ℚ(i) with the coefficient blowup
        // that makes one reach for subresultants, plus a reducing GCD on every BigRational
        // operation; the Qc section is banner-scoped to "small degrees" and degree 80 is already
        // not small), while what it saves is 2·m_D sample points out of resDeg(resDeg−1), under
        // 1% once resDeg reaches the hundreds. skipLeadingFormSharpening takes the crude bound
        // instead, the same choice the scout methods of this file already make, and it weakens
        // the CERTIFICATE not at all: dBound is an UPPER bound either way, and deg_q D and v_q D
        // are read off the interpolant, never off the bound. What it must NOT skip is the leading
        // form itself, whose construction carries the q → ∞ cross-check that the AT leading form
        // divides charpoly(C_block); that guard is cheap beside the bivariate build and matters
        // most on a block shape nobody has run, so it stays outside the branch. m_D is then
        // reported as −1 rather than guessed.
        int mD = -1;
        int dBound = resDeg * (resDeg - 1);
        var (fResLead, lcD) = ResidualLeadingFormFromPencil(cBlk, sectors);
        if (!skipLeadingFormSharpening)
        {
            (_, mD) = SquarefreeLayers(FromZi(fResLead));
            if ((mD == 0) != !lcD.Equals(GaussianInteger.Zero))
                throw new InvalidOperationException("m_D inconsistent with disc_x(f_res).");
            dBound -= 2 * mD;
        }

        BigInteger rowF = BigInteger.Zero;
        foreach (var coeff in fRes) rowF += HNormQ(coeff);
        BigInteger rowFp = resDeg * rowF;
        BigInteger normD = BigInteger.Pow(rowF, Math.Max(resDeg - 1, 0)) * BigInteger.Pow(rowFp, resDeg);
        int lcDivisorBoundD = (int)(2 * normD.GetBitLength() / 30) + 1;

        // The q-even halving (the checkerboard corollary of 2026-08-12). Licensed by an EXACT
        // ℤ[i] read of F_res, never by the block's index pair: where the Λ₀-shifted support has
        // one j + k parity, D(q) = E(q²) with deg E = deg D / 2, so E is interpolated from the
        // values at q0 = 0…dBound/2 fed in at the abscissae m = q0², and the sample points per
        // prime HALVE. dBound is even on BOTH bound branches (resDeg(resDeg−1) is a product of
        // consecutive integers, and the sharpening subtracts the even 2·m_D), so dBound/2 bounds
        // deg E whether or not the sharpening ran.
        // Fold-fixedness is the reason to EXPECT the read (the checkerboard is a theorem there,
        // on the full charpoly, its descent to F_res measured per block), so a fold-fixed block
        // that fails the read contradicts that measurement and throws rather than quietly
        // costing twice; forceDenseSampling is the way past it.
        bool halve = QEvenSamplingLicence(fRes, n, wKet, wBra, forceDenseSampling);

        const int extra = 24;
        int nodeCount = (halve ? dBound / 2 : dBound) + 1;
        int samples = nodeCount + extra;
        var perPrime = new List<(int P, int DegD, int VD, int[] DStrip)>();
        int sampled = 0;

        // THE PRIME LOOP RUNS IN PARALLEL (2026-08-13), because it is embarrassingly so: each
        // prime's pass reads the shared immutable fRes and allocates only locals
        // (SampleDiscAtPrime). What is NOT left to the scheduler is the ANSWER. A batch holds
        // exactly the number of primes still wanted, its results are folded in PRIME order, and
        // counting stops at the same sampled = cap + 1 the sequential loop stopped at, so
        // perPrime, PrimesSampled and the chosen LayerPrime are identical to the sequential
        // device; the committed value pins at N=4 and N=6 are what hold that. Chunking keeps the
        // progress log alive and bounds the retained DStrip arrays, which at a large block are
        // the memory story (~76 KB per prime at resDeg 224).
        var clock = System.Diagnostics.Stopwatch.StartNew();
        long candidate = (1L << 30) + 1;
        while (candidate % 4 != 1) candidate += 2;
        int tried = 0;
        const int chunk = 256;

        while (sampled <= lcDivisorBoundD && tried < maxPrimes)
        {
            int want = Math.Min(lcDivisorBoundD + 1 - sampled, chunk);
            var batch = new List<(int P, int Root)>(want);
            for (; batch.Count < want && tried < maxPrimes; candidate += 4)
            {
                if (!IsPrime(candidate)) continue;
                tried++;
                int p = checked((int)candidate);
                int? root = SqrtMinusOneFast(p);
                if (root is null) continue;
                batch.Add((p, root.Value));
            }
            if (batch.Count == 0) break;

            var results = new (int P, int DegD, int VD, int[] DStrip)?[batch.Count];
            try
            {
                System.Threading.Tasks.Parallel.For(0, batch.Count, i =>
                    results[i] = SampleDiscAtPrime(
                        fRes, batch[i].P, batch[i].Root, resDeg, samples, nodeCount, halve));
            }
            catch (AggregateException ae) when (ae.InnerException is not null)
            {
                // the per-prime throws are diagnostics (non-monic reduction, failed verification,
                // an undersized prime); rethrow one as itself so the message survives. WHICH one
                // is the scheduler's, not the lowest prime, so a systemic failure names a
                // different p on each run; the siblings are dropped. That is the one place this
                // device is not deterministic, and it is in the diagnosis, never in the answer.
                System.Runtime.ExceptionServices.ExceptionDispatchInfo.Capture(ae.InnerException).Throw();
                throw;
            }

            foreach (var r in results)
            {
                if (r is null) continue;        // content prime: uncounted, see the (1,2) device
                sampled++;
                perPrime.Add(r.Value);
                if (log is not null && sampled % 25 == 0)
                    log($"sampled={sampled} (bound {lcDivisorBoundD}), p={r.Value.P}, {clock.ElapsedMilliseconds} ms");
            }
        }

        int trueDegD = -1, eD = int.MaxValue;
        foreach (var s in perPrime)
        {
            if (s.DegD > trueDegD) trueDegD = s.DegD;
            if (s.VD < eD) eD = s.VD;
        }
        if (eD == int.MaxValue) eD = -1;

        long layerPrime = 0;
        int[] layerDegrees = Array.Empty<int>();
        foreach (var s in perPrime)
        {
            if (s.DegD != trueDegD || s.VD != eD) continue;
            layerPrime = s.P;
            var layers = FpSquarefreeLayers(s.DStrip, s.P);
            layerDegrees = new int[layers.Count];
            for (int k = 0; k < layers.Count; k++) layerDegrees[k] = DegP(layers[k]);
            break;
        }
        bool certified = layerPrime != 0 && sampled > lcDivisorBoundD && trueDegD >= 0;

        return new WeightDiscMultiplicityReport(
            n, wKet, wBra, rOdd, blockDim, atDeg, resDeg, mD, dBound,
            trueDegD, eD, lcDivisorBoundD, layerPrime,
            layerDegrees, layerDegrees.Length, sampled, certified,
            samples, halve, skipLeadingFormSharpening);
    }

    /// <summary>The Λ-stratum read the layer profile cannot supply: at the ×2-cleared value
    /// Λ₀ = <paramref name="lambdaCleared"/> (the census's μ = 0 singleton line is physical
    /// λ = −N, i.e. Λ₀ = −2N), the q-locus where Λ₀ is a REPEATED Λ-root of the residual is
    /// exactly the common-root set of f(q) = F_res(Λ₀, q) and g(q) = (∂F_res/∂Λ)(Λ₀, q), read
    /// mod p as deg gcd(f, g). On that line an UNDOUBLED transversal crossing is already an
    /// order-2 disc zero (the Klein orbit degenerates to a singleton), so the multiplicity-2
    /// disc layer cannot tell it from an off-line doubled defective locus — this gcd isolates
    /// the stratum so the layer read stays clean. Returns the gcd's q-valuation beside its
    /// degree (roots at q = 0 are not physical couplings), the f-layer profile (a perfect square,
    /// f-layers [0, k], says every Λ₀-touch is transversal: NO defective locus sits on the line,
    /// the only source of an odd-order f-root being a √-branch exactly there), and
    /// <c>EvenInShiftedLambda</c>: whether EVERY odd Taylor coefficient of F_res at Λ₀ vanishes
    /// mod p, i.e. F_res ∈ GF(p)[(Λ−Λ₀)², q] — the full evenness read (GDeg = −1 alone certifies
    /// only the first odd coefficient a₁; a review round caught that gap and this member closes
    /// it). Scout-grade, one prime.</summary>
    public static (int P, int FDeg, int GDeg, int GcdDeg, int GcdQVal, int[] FLayerDegrees, bool EvenInShiftedLambda)
        WeightLambdaStratumAtNthPrime(
        int n, int wKet, int wBra, bool rOdd, int lambdaCleared, int nth)
    {
        if (nth < 0) throw new ArgumentOutOfRangeException(nameof(nth));

        var (aBlk, cBlk) = WeightSectorPencil(n, wKet, wBra, rOdd);
        AssertHoppingSymmetry(cBlk, null);
        var sectors = F89AtFactorReconstruction.ClearedAtSectorsFromPencil(aBlk, cBlk);
        int atDeg = 0;
        foreach (var s in sectors) atDeg += s.KCharpoly.Length - 1;
        var chi = PencilCharpolyBivariate(aBlk, cBlk);
        var fRes = BivariateDivideExact(chi, BivariateAtFactor(sectors));
        WeightCrossCheckResidual(fRes, n, wKet, wBra, rOdd, sectors, aBlk.GetLength(0) - atDeg);

        long candidate = (1L << 30) + 1;
        while (candidate % 4 != 1) candidate += 2;
        int seen = 0;
        for (; ; candidate += 4)
        {
            if (!IsPrime(candidate)) continue;
            int p = checked((int)candidate);
            int? root = SqrtMinusOneFast(p);
            if (root is null) continue;
            if (seen++ < nth) continue;

            var fResP = ReduceBivariate(fRes, p, root.Value);
            int len = 0;
            foreach (var qs in fResP) len = Math.Max(len, qs.Length);
            var f = new int[len];
            var g = new int[len];
            long lam0 = ((lambdaCleared % p) + p) % p;
            long pw = 1, pwPrev = 0;
            for (int k = 0; k < fResP.Length; k++)
            {
                for (int i = 0; i < fResP[k].Length; i++)
                {
                    f[i] = (int)((f[i] + fResP[k][i] * pw) % p);
                    if (k >= 1)
                        g[i] = (int)((g[i] + fResP[k][i] * pwPrev % p * k) % p);
                }
                pwPrev = pw;
                pw = pw * lam0 % p;
            }
            var gcd = GcdModP(f, g, p);
            int gcdDeg = DegP(gcd);
            int gcdVal = gcdDeg < 0 ? -1 : QValuation(gcd);
            // The FULL evenness read, EXACT over ℤ[i] (a review round caught that g ≡ 0 alone
            // certifies only the FIRST odd Taylor coefficient a₁ at Λ₀; a second round noted a
            // mod-p read stays mod-p evidence, while the bivariate is exact and in hand):
            // Taylor-shift F_res to Λ₀ via a_m(q) = Σ_{k≥m} C(k,m)·c_k(q)·Λ₀^{k−m} over ℤ[i]
            // and require EVERY odd a_m to be identically zero. True ⟺ F_res is a polynomial in
            // (Λ − Λ₀)², the within-sector fold as algebra, no prime involved.
            bool evenShifted = ExactOddTaylorAllZero(fRes, lambdaCleared);
            int fDeg = DegP(f);
            int[] fLayerDegs;
            if (fDeg < 0) fLayerDegs = Array.Empty<int>();
            else
            {
                int vf = QValuation(f);
                var fLayers = FpSquarefreeLayers(StripQ(f, vf), p);
                fLayerDegs = new int[fLayers.Count];
                for (int k = 0; k < fLayers.Count; k++) fLayerDegs[k] = DegP(fLayers[k]);
            }
            return (p, fDeg, DegP(g), gcdDeg, gcdVal, fLayerDegs, evenShifted);
        }
    }

    /// <summary>Scout for the even-N closing certificate (the sideways_spin_ladder arc's Sturm/gcd
    /// route): recover Re D and Im D mod p SEPARATELY and read their gcd in GF(p)[q].
    ///
    /// <para>The two-embedding trick: a split prime p ≡ 1 mod 4 has TWO square roots of −1, r and
    /// p−r. Reducing the ℤ[i] discriminant under i ↦ r gives (Re D + r·Im D) mod p; under i ↦ p−r
    /// it gives (Re D − r·Im D) mod p; the half-sum recovers Re D mod p and the half-difference,
    /// divided by 2r, recovers Im D mod p, coefficient-wise. A real q with D(q) = 0 forces Re D and
    /// Im D (real polynomials) to vanish together, i.e. a nontrivial factor of gcd(Re D, Im D) over
    /// ℚ; by Gauss, a nontrivial ℚ-gcd of integer polynomials survives reduction mod any prime
    /// preserving the degrees (the lemma <c>simulations/o2b_gcd_certificate.py</c> states verbatim
    /// for the UNSTRIPPED gcd; this method takes the gcd of q-power-STRIPPED reductions, so it
    /// additionally needs mod-p valuation preservation, or a common root can collapse onto q = 0
    /// and be stripped away; the full premise structure is stated on
    /// <see cref="CertifyDiscReImGcd"/>). And D vanishes at EVERY repeated Λ-root, defective and diabolic
    /// alike, so <c>GcdDeg == 0</c> at a good prime certifies: no real q₀ ≠ 0 with ANY coalescence,
    /// window-free. This method is the mod-p half of that statement, a SCOUT with TWO uncertified
    /// premises: (a) degree preservation (the lc-divisor-bound argument of the certificate methods),
    /// and (b) valuation preservation — <see cref="StripQ"/> removes q^v with the MOD-P valuation,
    /// so if p divided all the low coefficients of Re or Im, a true nonzero common root could
    /// collapse onto q = 0 and be stripped away (the same caveat CertifyDiscMultiplicity certifies
    /// as TrueQValuationD). A nonzero-looking gcd can only shrink under more primes; a zero-looking
    /// one is certified only once (a) and (b) hold; <see cref="CertifyDiscReImGcd"/> discharges
    /// both from D's own certified degree and valuation (promotion landed 2026-08-10).</para>
    ///
    /// <para>Controls: at N=4 the disc is REAL (the self-fold antiunitary makes the codimension 1,
    /// which is WHY four real defective loci and the diabolic quartic exist there — the foil that
    /// makes an empty N=6 axis non-trivial; committed anchors: the B3 entry in
    /// <c>docs/CAUGHT_ERRORS.md</c>, the mechanism paragraph of
    /// <c>experiments/F89_PATH_K_DIABOLIC.md</c>, the disc-reality gate of
    /// <c>simulations/even_n_literal_real_count.py</c>), so <c>ImIsZero</c> must come back true
    /// there; at N=6 the disc is genuinely complex and the question is live. Consistency gate, per
    /// prime and non-vacuous: the i ↦ r reduction IS D mod p, so max(deg Re, deg Im) must equal
    /// deg(D mod p) and min(v_Re, v_Im) must equal v_q(D mod p); a throw guards both.</para></summary>
    public static (int P, int Root, int DegRe, int VRe, int DegIm, int VIm, bool ImIsZero, int GcdDeg)
        DiscReImGcdAtNthPrime(int n, bool rOdd, int nth)
    {
        if (n < 4) throw new ArgumentException("n ≥ 4 (either parity; see CertifyDiscMultiplicity).", nameof(n));
        if (nth < 0) throw new ArgumentOutOfRangeException(nameof(nth));

        var (aBlk, cBlk) = BlockPencil(n, rOdd);
        AssertHoppingSymmetry(cBlk, rOdd ? null : F89PathKSeDeBlock.SymOrbitSizes(n));
        var sectors = F89AtFactorReconstruction.ClearedAtSectors(n - 1, rOdd);
        int atDeg = 0;
        foreach (var s in sectors) atDeg += s.KCharpoly.Length - 1;
        int resDeg = aBlk.GetLength(0) - atDeg;

        var chi = PencilCharpolyBivariate(aBlk, cBlk);
        var fRes = BivariateDivideExact(chi, BivariateAtFactor(sectors));
        var (fResLead, _) = ResidualLeadingForm(n, rOdd);
        var (_, mD) = SquarefreeLayers(FromZi(fResLead));
        int dBound = resDeg * (resDeg - 1) - 2 * mD;
        int samples = dBound + 1 + 24;

        long candidate = (1L << 30) + 1;
        while (candidate % 4 != 1) candidate += 2;
        int seen = 0;
        for (; ; candidate += 4)
        {
            if (!IsPrime(candidate)) continue;
            int p = checked((int)candidate);
            int? root = SqrtMinusOneFast(p);
            if (root is null) continue;
            if (seen++ < nth) continue;
            int r = root.Value;

            var dp = new int[2][];
            foreach (var (emb, k) in new[] { (r, 0), (p - r, 1) })
            {
                var fResP = ReduceBivariate(fRes, p, emb);
                var dSamp = new int[samples];
                for (int q0 = 0; q0 < samples; q0++)
                {
                    var fL = EvalBivariateModP(fResP, q0, p);
                    if (DegP(fL) != resDeg || fL[resDeg] != 1)
                        throw new InvalidOperationException($"reduced residual not monic of degree {resDeg} at q0={q0} (i↦{emb}).");
                    dSamp[q0] = DiscriminantModP(fL, p);
                }
                var dNodes = new int[dBound + 1];
                Array.Copy(dSamp, dNodes, dBound + 1);
                dp[k] = InterpolateAtIntegerNodes(dNodes, p);
                for (int q0 = dBound + 1; q0 < samples; q0++)
                    if (EvalModP(dp[k], q0, p) != dSamp[q0])
                        throw new InvalidOperationException($"D interpolant fails verification at q0={q0} (p={p}, i↦{emb}).");
            }

            // half-sum / half-difference, all products in long (the 2^30 overflow trap is live here)
            int len = Math.Max(dp[0].Length, dp[1].Length);
            long inv2 = InvModP(2, p);
            long inv2r = InvModP((int)(2L * r % p), p);
            var re = new int[len];
            var im = new int[len];
            for (int i = 0; i < len; i++)
            {
                long a = i < dp[0].Length ? dp[0][i] : 0;
                long b = i < dp[1].Length ? dp[1][i] : 0;
                re[i] = (int)((a + b) % p * inv2 % p);
                im[i] = (int)(((a - b + p) % p) * inv2r % p);
            }

            bool imZero = DegP(im) < 0;
            int vRe = DegP(re) < 0 ? -1 : QValuation(re);
            int vIm = imZero ? -1 : QValuation(im);

            // the non-vacuous consistency gate: dp[0] IS D mod p, and over ℤ the leading terms of
            // the real Re/Im cannot cancel in D = Re + i·Im, so degrees and q-valuations must
            // agree. Mod p a coincidental cancellation Re_k ≡ −r·Im_k has probability ~1/p; it
            // throws here, fail-loud and TERMINAL (the prime sequence is deterministic, so a hit
            // would recur every run and is handled then, not silently skipped).
            int degD = DegP(dp[0]);
            int vD = degD < 0 ? -1 : QValuation(dp[0]);
            if (Math.Max(DegP(re), DegP(im)) != degD)
                throw new InvalidOperationException(
                    $"max(deg Re, deg Im) = {Math.Max(DegP(re), DegP(im))} != deg(D mod p) = {degD} (p={p}).");
            int vMin = imZero ? vRe : (DegP(re) < 0 ? vIm : Math.Min(vRe, vIm));
            if (vMin != vD)
                throw new InvalidOperationException(
                    $"min(v_Re, v_Im) = {vMin} != v_q(D mod p) = {vD} (p={p}).");

            int gcdDeg = -1;
            if (!imZero && DegP(re) >= 0)
                gcdDeg = DegP(GcdModP(StripQ(re, vRe), StripQ(im, vIm), p));
            return (p, r, DegP(re), vRe, DegP(im), vIm, imZero, gcdDeg);
        }
    }

    /// <summary>The verdict of <see cref="CertifyDiscReImGcd"/>. <c>Certified</c> with
    /// <c>GcdDeg == 0</c> (equivalently <c>NoRealNonzeroCoalescence</c>) is the window-free even-N
    /// closing statement of the sideways_spin_ladder arc's Sturm/gcd route: NO real q ≠ 0 exists at
    /// which the (1,2) residual F_res has ANY repeated Λ-root, defective or diabolic, at that
    /// (n, parity). One-way as ever: <c>GcdDeg &gt; 0</c> or <c>Certified == false</c> is
    /// refine-not-refute. The two fail-closed modes are distinguishable:
    /// <c>DiscLayersCertified == false</c> means the D-side multi-prime device itself failed, while
    /// <c>DiscLayersCertified &amp;&amp; !Certified</c> means no scanned prime saw a genuinely
    /// complex disc attaining both invariants (the N=4 real-disc control's mode), or only ones with
    /// a degenerate Re reduction (Re D ≡ 0 mod p, where the gcd would be meaningless). Field fine
    /// print: <c>PrimesSampled</c>, the True… invariants and <c>LcDivisorBoundD</c> come from the
    /// D-side device (PrimesSampled is NOT the ≤ 8 Re/Im scout primes); on a declined report the
    /// True… fields are whatever the device returned and CertPrime/Root/DegRe/VRe/DegIm/VIm/GcdDeg
    /// are 0/−1 sentinels.</summary>
    public sealed record DiscReImGcdReport(
        int N, bool ROdd,
        int TrueDiscriminantDegree, int TrueQValuationD, int LcDivisorBoundD, int PrimesSampled,
        long CertPrime, int Root, int DegRe, int VRe, int DegIm, int VIm, int GcdDeg,
        bool DiscLayersCertified, bool Certified, bool NoRealNonzeroCoalescence);

    /// <summary>The PROMOTED Re/Im gcd certificate (scout: <see cref="DiscReImGcdAtNthPrime"/>; its
    /// two named premises are discharged here, which is the whole method): a real q₀ ≠ 0 with
    /// D(q₀) = 0 forces Re D(q₀) = Im D(q₀) = 0 (real polynomials at a real point), so
    /// gcd(Re D, Im D) over ℚ is nontrivial; <c>GcdDeg == 0</c> mod a GOOD prime refutes that, and
    /// this method certifies the prime is good.
    ///
    /// <para>Both preservation premises reduce to D's OWN certified invariants, no separate Re/Im
    /// bound needed. Two exact ℤ-level identities: deg_q D = max(deg Re, deg Im) and
    /// v_q(D) = min(v_Re, v_Im), because at the leading (resp. trailing) coefficient the real and
    /// imaginary parts sit in orthogonal components of ℤ[i] and cannot cancel. The scout's per-prime
    /// cross-gate enforces the same two identities mod p (throw on the ~1/p accident). At a prime
    /// attaining <see cref="DiscMultiplicityReport.TrueDiscriminantDegree"/> AND
    /// <see cref="DiscMultiplicityReport.TrueQValuationD"/> (both certified by
    /// <see cref="CertifyDiscMultiplicity"/>'s Hadamard/lc-divisor device), the chain
    /// max(deg Re mod p, deg Im mod p) = deg(D mod p) = deg_q D = max(deg Re, deg Im), with mod-p
    /// degrees only able to DROP, forces the max-attaining side to preserve its degree; the mirrored
    /// chain on valuations (only able to RISE) forces the min-attaining side to preserve its
    /// q-valuation. Possibly different sides, and that suffices.</para>
    ///
    /// <para>The Gauss step, with the two preserved handles: suppose a real q₀ ≠ 0 with D(q₀) = 0,
    /// and let m ∈ ℤ[q] be the primitive minimal polynomial of q₀ (deg m ≥ 1, m ≠ q so m(0) ≠ 0).
    /// Then m | Re and m | Im over ℤ (Gauss; both are nonzero over ℤ, witnessed by their nonzero
    /// reductions at the accepted prime), hence m | Re₁ and m | Im₁, the q-power-stripped parts.
    /// The valuation-preserved side X has q ∤ (X₁ mod p), and m mod p divides it, so q ∤ (m mod p).
    /// The degree-preserved side Y has p ∤ lc(Y) = lc(m)·lc(cofactor), so p ∤ lc(m) and
    /// deg(m mod p) = deg m ≥ 1. So m mod p is a nonconstant common divisor, coprime to q, of BOTH
    /// stripped reductions: GcdDeg ≥ 1. Contrapositive: GcdDeg == 0 at the certified prime proves no
    /// such q₀ exists, window-free. The certifying prime must also see Im ≢ 0 mod p (a nonzero
    /// reduction certifies D genuinely complex; the N=4 self-fold control has Im ≡ 0 and MUST fail
    /// closed here, its real disc carrying four honest real defective loci).</para></summary>
    public static DiscReImGcdReport CertifyDiscReImGcd(
        int n, bool rOdd, int maxPrimes = 20000, Action<string>? log = null)
    {
        var d = CertifyDiscMultiplicity(n, rOdd, maxPrimes, log);
        if (d.DiscLayersCertified)
        {
            // at most LcDivisorBoundD primes can drop the degree and at most as many can raise the
            // valuation; the first few split primes are overwhelmingly good, so a short scan finds
            // one attaining both. The cap of 8 is a deliberate fail-closed budget, not a bound
            // (LcDivisorBoundD would be the bound): eight straight misses end in the declined
            // report below rather than a long loop, and one-way-ness makes declining sound.
            for (int nth = 0; nth < 8; nth++)
            {
                var s = DiscReImGcdAtNthPrime(n, rOdd, nth);
                if (s.ImIsZero || s.DegRe < 0) continue;
                if (Math.Max(s.DegRe, s.DegIm) != d.TrueDiscriminantDegree) continue;
                if (Math.Min(s.VRe, s.VIm) != d.TrueQValuationD) continue;
                log?.Invoke($"certifying prime p={s.P} (nth={nth}): gcd deg = {s.GcdDeg}");
                return new DiscReImGcdReport(
                    n, rOdd, d.TrueDiscriminantDegree, d.TrueQValuationD, d.LcDivisorBoundD,
                    d.PrimesSampled, s.P, s.Root, s.DegRe, s.VRe, s.DegIm, s.VIm, s.GcdDeg,
                    DiscLayersCertified: true, Certified: true,
                    NoRealNonzeroCoalescence: s.GcdDeg == 0);
            }
        }
        return new DiscReImGcdReport(
            n, rOdd, d.TrueDiscriminantDegree, d.TrueQValuationD, d.LcDivisorBoundD, d.PrimesSampled,
            CertPrime: 0, Root: 0, DegRe: -1, VRe: -1, DegIm: -1, VIm: -1, GcdDeg: -1,
            DiscLayersCertified: d.DiscLayersCertified, Certified: false,
            NoRealNonzeroCoalescence: false);
    }

    /// <summary>The D-only exact per-point residual: block charpoly at q0, divided by the AT factor.
    /// Independent of the bivariate path, and free of the corner block.</summary>
    private static GaussianInteger[] ResidualPerPointAt(int n, bool rOdd, int q0)
    {
        var blk = rOdd
            ? F89PathKSeDeBlock.BuildTwoTimesROddBlock(q0, n)
            : F89PathKSeDeBlock.BuildTwoTimesSymBlock(q0, n);
        var (res, rem) = GaussianPolynomial.DivMod(GaussianMatrixCharpoly.Characteristic(blk),
                                                   AtFactorAt(n, rOdd, q0));
        if (rem.Length != 0)
            throw new InvalidOperationException($"AT factor does not divide the block charpoly at q0={q0}.");
        return res;
    }

    /// <summary>The q → ∞ leading form of F_res and the leading coefficient of its discriminant, with no
    /// target block built (the D-only half of <see cref="LeadingFormsCore"/>).</summary>
    private static (GaussianInteger[] ResidualLeadingForm, GaussianInteger DiscriminantLeadingCoeff)
        ResidualLeadingForm(int n, bool rOdd)
    {
        var (_, cBlk) = BlockPencil(n, rOdd);
        var sectors = F89AtFactorReconstruction.ClearedAtSectors(n - 1, rOdd);
        return ResidualLeadingFormFromPencil(cBlk, sectors);
    }

    /// <summary>The pencil-general body of <see cref="ResidualLeadingForm"/>: charpoly of the hopping
    /// direction C divided by the AT leading form ∏_sectors Σ_m p_m·i^{d−m}·x^m (each AT strand is
    /// Λ = r0 + i·q·κ, so at q → ∞ its leading unit is i^{d−m}).</summary>
    private static (GaussianInteger[] ResidualLeadingForm, GaussianInteger DiscriminantLeadingCoeff)
        ResidualLeadingFormFromPencil(GaussianInteger[,] cBlk, IReadOnlyList<AtSector> sectors)
    {
        var fFull = GaussianMatrixCharpoly.Characteristic(cBlk);
        var fAt = new[] { GaussianInteger.One };
        foreach (var s in sectors)
        {
            int d = s.KCharpoly.Length - 1;
            var polyX = new GaussianInteger[d + 1];
            for (int m = 0; m <= d; m++)
            {
                var unit = GaussianInteger.One;
                for (int e = 0; e < d - m; e++) unit *= GaussianInteger.I;
                polyX[m] = new GaussianInteger(s.KCharpoly[m], 0) * unit;
            }
            fAt = GaussianPolynomial.Multiply(fAt, polyX);
        }
        var (fRes, rem) = GaussianPolynomial.DivMod(fFull, fAt);
        if (rem.Length != 0)
            throw new InvalidOperationException("the AT leading form does not divide charpoly(C_block).");
        return (fRes, GaussianPolynomial.Discriminant(fRes));
    }

    /// <summary>The parameterised multi-prime proof engine. Accumulates split primes p ≥ 2³⁰ until their
    /// product exceeds the Mignotte/Hadamard proof bound; every used prime must see gcd(R mod p, D mod p) =
    /// c·q^e (else the report returns early with SharedIsQPowerAtEveryPrime = false: refine, not refute).
    /// Throws on structural violations (AT division remainder, interpolant mismatch, bivariate/per-point
    /// cross-check failure), because no prime can repair those.</summary>
    private static CompleteReport CertifyCore(int n, bool rOdd, int tWKet, int tWBra,
        GaussianInteger composeA, GaussianInteger composeB, int maxPrimes, Action<string>? log)
    {
        if (n < 5 || n % 2 == 0)
            throw new ArgumentException("odd n ≥ 5 (the corner (p_c+1, p_c+1) needs odd N).", nameof(n));

        // ---- the exact integer pencils and the O(1)-branch premise ----
        var (aBlk, cBlk) = BlockPencil(n, rOdd);
        var (aTgt, cTgt) = TargetPencil(n, tWKet, tWBra);
        AssertHoppingSymmetry(cBlk, rOdd ? null : F89PathKSeDeBlock.SymOrbitSizes(n));
        AssertHoppingSymmetry(cTgt, null);

        var sectors = F89AtFactorReconstruction.ClearedAtSectors(n - 1, rOdd);
        int atDeg = 0;
        foreach (var s in sectors) atDeg += s.KCharpoly.Length - 1;
        int blockDim = aBlk.GetLength(0);
        int resDeg = blockDim - atDeg;
        int corDeg = aTgt.GetLength(0);

        // ---- the exact bivariate layer: F_res(Λ, q) and G(Λ, q) = F2_target(composeA·Λ + composeB, q) ----
        var (fRes, g) = BuildBivariate(aBlk, cBlk, aTgt, cTgt, sectors, composeA, composeB);
        if (fRes.Length - 1 != resDeg)
            throw new InvalidOperationException($"bivariate residual degree {fRes.Length - 1} ≠ {resDeg}.");

        // cross-check the bivariate layer against the independent per-point path at q0 = 2, 3
        foreach (int q0 in new[] { 2, 3 })
        {
            var (resPt, gPt, _, _) = ExactSampleCore(n, rOdd, q0, tWKet, tWBra, composeA, composeB);
            AssertPolyEqual(EvalBivariateAtQ(fRes, q0), resPt, $"bivariate F_res at q0={q0}");
            AssertPolyEqual(EvalBivariateAtQ(g, q0), gPt, $"bivariate G at q0={q0}");
        }

        // ---- G2: the q → ∞ leading forms and the proven degree bounds ----
        var (fResLead, gLead, lcR, lcD) = LeadingFormsCore(n, rOdd, tWKet, tWBra, composeA);
        var (fLayers, fChainSum) = SquarefreeLayers(FromZi(fResLead));
        var (gLayers, _) = SquarefreeLayers(FromZi(gLead));
        int mD = fChainSum;                                 // Σ_c C(mult_{f_res}(c), 2)
        int mR = 0;                                         // Σ_c mult_{f_res}(c)·mult_g(c)
        foreach (var e in fLayers)
            foreach (var h in gLayers)
                mR += QcDeg(QcGcd(e, h));
        if ((mR == 0) != !lcR.Equals(GaussianInteger.Zero))
            throw new InvalidOperationException("m_R inconsistent with Res_x(f_res, g).");
        if ((mD == 0) != !lcD.Equals(GaussianInteger.Zero))
            throw new InvalidOperationException("m_D inconsistent with disc_x(f_res).");
        int rBound = resDeg * corDeg - mR;
        int dBound = resDeg * (resDeg - 1) - 2 * mD;

        // ---- the Mignotte/Hadamard proof bound on gcd(R, D)'s coefficient heights ----
        BigInteger rowF = BigInteger.Zero, rowG = BigInteger.Zero;
        foreach (var coeff in fRes) rowF += HNormQ(coeff);
        foreach (var coeff in g) rowG += HNormQ(coeff);
        BigInteger normR = BigInteger.Pow(rowF, corDeg) * BigInteger.Pow(rowG, resDeg);
        BigInteger heightH = BigInteger.Pow(2, rBound + 16) * normR;    // Mignotte 2^deg·M(R), with slack
        BigInteger proofBound = 4 * heightH * heightH;                  // Gaussian norm: N(c) ≤ 2·|c|²

        // ---- the same Hadamard device for D, so the DISC's degree and q-valuation are certified too ----
        // D = disc_Λ(F_res) = Res_Λ(F_res, ∂_Λ F_res) (F_res is monic in Λ). Its Sylvester matrix has
        // resDeg−1 rows of F_res coefficients and resDeg rows of ∂_Λ F_res coefficients, so Hadamard gives
        // ‖D‖ ≤ rowF^{resDeg−1} · rowFp^{resDeg} with rowFp ≤ resDeg·rowF. A sampled prime ideal π_p can
        // lose deg_q(D) only if π_p | lc_q(D), and can raise v_q(D) only if π_p divides the coefficient of
        // q^{v_q(D)}; both coefficients have 1-norm ≤ normD, the ideals lie over distinct rational primes
        // p ≥ 2³⁰, so at most lcDivisorBoundD sampled primes can lose the degree, and (separately) at most
        // lcDivisorBoundD can raise the valuation. Sampling past that bound therefore certifies EACH of
        // trueDegD = max_p deg(D mod p) = deg_q D and eD = min_p v_q(D mod p) = v_q(D): one bound per
        // statement, NOT one bound for their union. A prime good at BOTH ends is then not guaranteed a
        // priori by this count; it is SEARCHED for below and the report fails closed (layerPrime = 0,
        // DiscLayersCertified = false) if none of the sampled primes attains both. Only at such a prime
        // does "max multiplicity mod p" bound the true max multiplicity from above (no root escaped to ∞,
        // none collapsed onto q = 0 and got stripped with the q-power).
        BigInteger rowFp = resDeg * rowF;
        BigInteger normD = BigInteger.Pow(rowF, Math.Max(resDeg - 1, 0)) * BigInteger.Pow(rowFp, resDeg);
        int lcDivisorBoundD = (int)(2 * normD.GetBitLength() / 30) + 1;

        // ---- the multi-prime loop ----
        const int extra = 24;
        int samples = Math.Max(rBound, dBound) + 1 + extra;     // dBound can exceed rBound for small targets

        // THE DEGREE CERTIFICATE. deg(R mod p) < deg_q R exactly when the chosen prime ideal
        // π_p = (p, i − r) divides lc_q(R) ∈ Z[i] (the kernel of the reduction i ↦ r). The sampled
        // ideals lie over distinct rational primes, pairwise coprime, so ∏ N(π_p) = ∏ p ≤ N(lc_q R)
        // ≤ normR² (every coefficient of R has 1-norm ≤ the Hadamard product normR); with every
        // p ≥ 2³⁰, at most 2·log₂(normR)/30 sampled primes can lose the degree. Sampling MORE
        // distinct split primes than lcDivisorBound therefore forces trueDegR = max_p deg(R mod p)
        // = deg_q R EXACTLY: the empirical degree is itself certified, and rBound serves only as
        // the (possibly loose) interpolation length. Loose happens: the identity composition has
        // true 412 < 422 (R-even) / 384 < 394 (R-odd) at N=5, where the corner-fold bound is tight.
        int lcDivisorBound = (int)(2 * normR.GetBitLength() / 30) + 1;

        BigInteger primeProduct = BigInteger.One;
        int sampled = 0, skipped = 0;
        int trueDegR = -1;
        bool sharedIsQPower = true;
        var perPrime = new List<(int P, int DegR, int DegD, int VR, int VD, bool SharedOk, int[] RStrip, int[] DStrip)>();

        var clock = System.Diagnostics.Stopwatch.StartNew();
        long candidate = (1L << 30) + 1;
        while (candidate % 4 != 1) candidate += 2;
        for (int tried = 0; tried < maxPrimes
                            && !(primeProduct > proofBound && sampled > lcDivisorBound && sampled > lcDivisorBoundD);
             candidate += 4)
        {
            if (!IsPrime(candidate)) continue;
            tried++;
            int p = checked((int)candidate);
            int? root = SqrtMinusOneFast(p);
            if (root is null) { skipped++; continue; }
            int r = root.Value;

            var fResP = ReduceBivariate(fRes, p, r);
            var gP = ReduceBivariate(g, p, r);

            var rSamp = new int[samples];
            var dSamp = new int[samples];
            for (int q0 = 0; q0 < samples; q0++)
            {
                var fL = EvalBivariateModP(fResP, q0, p);
                var gL = EvalBivariateModP(gP, q0, p);
                if (DegP(fL) != resDeg || fL[resDeg] != 1)
                    throw new InvalidOperationException($"reduced residual not monic of degree {resDeg} at q0={q0}.");
                rSamp[q0] = ResultantModP(fL, gL, p);
                dSamp[q0] = DiscriminantModP(fL, p);
            }

            var rNodes = new int[rBound + 1];
            Array.Copy(rSamp, rNodes, rBound + 1);
            var rp = InterpolateAtIntegerNodes(rNodes, p);
            var dNodes = new int[dBound + 1];
            Array.Copy(dSamp, dNodes, dBound + 1);
            var dp = InterpolateAtIntegerNodes(dNodes, p);
            for (int q0 = 0; q0 < samples; q0++)
            {
                if (q0 > rBound && EvalModP(rp, q0, p) != rSamp[q0])
                    throw new InvalidOperationException($"R interpolant fails verification at q0={q0} (p={p}).");
                if (q0 > dBound && EvalModP(dp, q0, p) != dSamp[q0])
                    throw new InvalidOperationException($"D interpolant fails verification at q0={q0} (p={p}).");
            }
            // deg(rp) ≤ rBound and deg(dp) ≤ dBound hold BY CONSTRUCTION (bound+1 nodes); the
            // interpolant's correctness is the extras loop above. Whether THIS prime kept the true
            // degree is decided against trueDegR = max over all sampled primes, certified once
            // sampled > lcDivisorBound; a prime below the max divided lc_q(R) and is demoted.
            sampled++;
            int degRp = DegP(rp);
            int vRp = QValuation(rp), vDp = QValuation(dp);
            var rStrip = StripQ(rp, vRp);
            var dStrip = StripQ(dp, vDp);
            bool sharedOk = DegP(GcdModP(rStrip, dStrip, p)) == 0;
            perPrime.Add((p, degRp, DegP(dp), vRp, vDp, sharedOk, rStrip, dStrip));

            if (degRp > trueDegR)
            {
                trueDegR = degRp;                                   // demotes every earlier lower-degree prime
                primeProduct = BigInteger.One;
                foreach (var s in perPrime)
                    if (s.DegR == trueDegR && s.SharedOk) primeProduct *= s.P;
            }
            else if (degRp == trueDegR && sharedOk)
                primeProduct *= p;

            if (!sharedOk && degRp == trueDegR && sampled > lcDivisorBound)
                break;                                              // a certified-good prime saw a shared factor: refine, not refute

            if (log is not null && sampled % 50 == 0)
                log($"sampled={sampled} (bound {lcDivisorBound}), skipped={skipped}, p={p}, trueDegR={trueDegR}, product digits={primeProduct.ToString().Length}/{proofBound.ToString().Length}, {clock.ElapsedMilliseconds} ms");
        }

        // ---- the verdict over the certified-good primes (those attaining trueDegR) ----
        int used = 0;
        long firstPrime = 0, lastPrime = 0;
        int degRpFirst = -1, degDpFirst = -1, vR = -1, vD = -1;
        foreach (var s in perPrime)
        {
            if (s.DegR != trueDegR) { skipped++; continue; }        // divided lc_q(R): demoted
            if (firstPrime == 0)
            {
                firstPrime = s.P;
                degRpFirst = s.DegR; degDpFirst = s.DegD; vR = s.VR; vD = s.VD;
            }
            if (!s.SharedOk) { sharedIsQPower = false; lastPrime = s.P; continue; }
            used++;
            lastPrime = s.P;
        }

        // ---- the CERTIFIED disc-layer reading (the β-exotic exclusion; see the class docstring) ----
        // trueDegD and eD are certified over ALL sampled primes (an R-demoted prime is still a perfectly
        // good reduction for D). The layer prime must attain BOTH: deg D preserved (no root escaped to ∞)
        // and v_q(D) minimal (no nonzero root collapsed onto q = 0 and got stripped with the q-power).
        int trueDegD = -1, eD = int.MaxValue;
        foreach (var s in perPrime)
        {
            if (s.DegD > trueDegD) trueDegD = s.DegD;
            if (s.VD < eD) eD = s.VD;
        }
        if (eD == int.MaxValue) eD = -1;

        long layerPrime = 0;
        int[] layerDegrees = Array.Empty<int>(), layerGcdDegrees = Array.Empty<int>();
        foreach (var s in perPrime)
        {
            if (s.DegD != trueDegD || s.VD != eD) continue;
            layerPrime = s.P;
            var layers = FpSquarefreeLayers(s.DStrip, s.P);
            layerDegrees = new int[layers.Count];
            layerGcdDegrees = new int[layers.Count];
            for (int k = 0; k < layers.Count; k++)
            {
                layerDegrees[k] = DegP(layers[k]);
                layerGcdDegrees[k] = DegP(GcdModP(s.RStrip, layers[k], s.P));
            }
            break;
        }
        // The layer prime's own report fields must be the certified ones, so the test can compare.
        if (layerPrime != 0) { degDpFirst = trueDegD; vD = eD; }
        bool discLayersCertified = layerPrime != 0 && sampled > lcDivisorBoundD && trueDegD >= 0;
        int maxDiscMultiplicity = layerDegrees.Length;

        // trueDegR ≥ 0 is load-bearing (2026-07-03 adversarial review): R ≡ 0 (all rp zero, trueDegR = −1)
        // voids both the lc-divisor count (lc undefined) and the Mignotte lift (H ~ D, unbounded by R's
        // heights), so the degenerate path must never report Complete.
        bool complete = sharedIsQPower && primeProduct > proofBound && sampled > lcDivisorBound && trueDegR >= 0;
        return new CompleteReport(
            n, rOdd,
            blockDim, atDeg, resDeg, corDeg,
            mR, mD,
            rBound, dBound, degRpFirst, degDpFirst,
            vR, vD, layerDegrees, layerGcdDegrees,
            used, skipped, sampled, lcDivisorBound, firstPrime, lastPrime,
            proofBound.ToString().Length, primeProduct.ToString().Length,
            sharedIsQPower, complete,
            trueDegD, eD, lcDivisorBoundD, layerPrime, maxDiscMultiplicity, discLayersCertified);
    }

    /// <summary>Diagnostic: the true deg_q(R) vs the analytic upper bound rBound for a given target block and
    /// Λ-composition, at one split prime (interpolating extra nodes so the true degree shows). The corner-fold
    /// bound is TIGHT (true = bound); the (1,1) identity-composition bound is LOOSE (true 10 below bound), because
    /// the identity leading forms share a common factor beyond the isolated-root collisions the mR count sees.
    /// This diagnosis is why CertifyCore's guard is the CERTIFIED empirical trueDegR (max over more sampled
    /// primes than can divide lc_q(R); landed 2026-07-03) rather than a tightness assumption on rBound.
    /// Reports resDeg, targetDeg, mR, rBound, lcR, trueDegR.</summary>
    public static string DebugDegreeReport(int n, bool rOdd, int tWKet, int tWBra, GaussianInteger cA, GaussianInteger cB)
    {
        var (aBlk, cBlk) = BlockPencil(n, rOdd);
        var (aTgt, cTgt) = TargetPencil(n, tWKet, tWBra);
        var sectors = F89AtFactorReconstruction.ClearedAtSectors(n - 1, rOdd);
        var (fRes, g) = BuildBivariate(aBlk, cBlk, aTgt, cTgt, sectors, cA, cB);
        var (fResLead, gLead, lcR, lcD) = LeadingFormsCore(n, rOdd, tWKet, tWBra, cA);
        int resDeg = fRes.Length - 1, targetDeg = aTgt.GetLength(0);
        var (fLayers, _) = SquarefreeLayers(FromZi(fResLead));
        var (gLayers, _) = SquarefreeLayers(FromZi(gLead));
        int mR = 0;
        foreach (var e in fLayers) foreach (var h in gLayers) mR += QcDeg(QcGcd(e, h));
        int rBound = resDeg * targetDeg - mR;
        long cand = (1L << 30) + 1; while (cand % 4 != 1) cand += 2;
        while (!(IsPrime(cand) && SqrtMinusOneFast((int)cand) is not null)) cand += 4;
        int p = (int)cand, root = SqrtMinusOneFast(p)!.Value;
        var fResP = ReduceBivariate(fRes, p, root); var gP = ReduceBivariate(g, p, root);
        int samples = rBound + 40;
        var rSamp = new int[samples];
        for (int q0 = 0; q0 < samples; q0++)
            rSamp[q0] = ResultantModP(EvalBivariateModP(fResP, q0, p), EvalBivariateModP(gP, q0, p), p);
        var rp = InterpolateAtIntegerNodes(rSamp, p);
        return $"resDeg={resDeg} targetDeg={targetDeg} mR={mR} rBound={rBound} " +
               $"lcR={(lcR.Equals(GaussianInteger.Zero) ? "ZERO" : "nonzero")} trueDegR={DegP(rp)} (p={p})";
    }

    private static (GaussianInteger[][] FRes, GaussianInteger[][] G) BuildBivariate(
        GaussianInteger[,] aBlk, GaussianInteger[,] cBlk,
        GaussianInteger[,] aTarget, GaussianInteger[,] cTarget, IReadOnlyList<AtSector> sectors,
        GaussianInteger composeA, GaussianInteger composeB)
    {
        var fFull = PencilCharpolyBivariate(aBlk, cBlk);
        var at = BivariateAtFactor(sectors);
        var fRes = BivariateDivideExact(fFull, at);
        var fTarget = PencilCharpolyBivariate(aTarget, cTarget);
        var g = BivariateComposeLinearLambda(fTarget, composeA, composeB);   // fold (−Λ−4N) OR identity (Λ)
        return (fRes, g);
    }

    /// <summary>Debug window: the exact bivariate F_res and G evaluated at integer q0 over Z[i].</summary>
    public static (GaussianInteger[] Residual, GaussianInteger[] FoldedCorner) DebugBivariateAt(int n, bool rOdd, int q0)
    {
        var (aBlk, cBlk) = BlockPencil(n, rOdd);
        var (aCor, cCor) = CornerPencil(n);
        var sectors = F89AtFactorReconstruction.ClearedAtSectors(n - 1, rOdd);
        var (fRes, g) = BuildBivariate(aBlk, cBlk, aCor, cCor, sectors,
            new GaussianInteger(-1, 0), new GaussianInteger(-4L * n, 0));
        return (EvalBivariateAtQ(fRes, q0), EvalBivariateAtQ(g, q0));
    }

    /// <summary>Debug window: one mod-p pipeline sample (resultant + discriminant) against the exact
    /// Z[i] values reduced mod p, at a chosen (p, q0).</summary>
    public static (bool Ok, int ExpectedR, int GotR, int ExpectedD, int GotD) DebugSampleCheck(
        int n, bool rOdd, int p, int q0)
    {
        var (aBlk, cBlk) = BlockPencil(n, rOdd);
        var (aCor, cCor) = CornerPencil(n);
        var sectors = F89AtFactorReconstruction.ClearedAtSectors(n - 1, rOdd);
        var (fRes, g) = BuildBivariate(aBlk, cBlk, aCor, cCor, sectors,
            new GaussianInteger(-1, 0), new GaussianInteger(-4L * n, 0));
        int r = SqrtMinusOneFast(p) ?? throw new ArgumentException("prime does not split", nameof(p));

        var fL = EvalBivariateModP(ReduceBivariate(fRes, p, r), q0, p);
        var gL = EvalBivariateModP(ReduceBivariate(g, p, r), q0, p);
        int gotR = ResultantModP(fL, gL, p);
        int gotD = DiscriminantModP(fL, p);

        var resQ = EvalBivariateAtQ(fRes, q0);
        var gQ = EvalBivariateAtQ(g, q0);
        int expR = ReduceScalar(GaussianPolynomial.Resultant(resQ, gQ), p, r);
        int expD = ReduceScalar(GaussianPolynomial.Discriminant(resQ), p, r);
        return (expR == gotR && expD == gotD, expR, gotR, expD, gotD);
    }

    /// <summary>The fully exact per-point Z[i] pipeline at ONE integer q0: the residual factor
    /// F_res(Λ, q0) (the block charpoly divided by the AT factor), the fold-composed corner polynomial
    /// G(Λ) = F2_corner(−Λ − 4N), and the two certificate values Res_Λ(F_res, G) and disc_Λ(F_res).
    /// The independent arithmetic path the bivariate layer is cross-checked against, and the gate's
    /// negative-control substrate (Res(F_res, F_res) = 0).</summary>
    public static (GaussianInteger[] Residual, GaussianInteger[] FoldedCorner,
        GaussianInteger ResultantValue, GaussianInteger DiscriminantValue)
        ExactSample(int n, bool rOdd, int q0)
        => ExactSampleCore(n, rOdd, q0, (n + 1) / 2 + 1, (n + 1) / 2 + 1,
                           new GaussianInteger(-1, 0), new GaussianInteger(-4L * n, 0));

    /// <summary>The fully exact per-point Z[i] pipeline at ONE integer q0, generalised over the target block
    /// (tWKet, tWBra) and the Λ-composition (composeA·Λ + composeB): F_res(Λ, q0) = block charpoly / AT factor,
    /// G(Λ) = F2_target(composeA·Λ + composeB), and Res_Λ(F_res, G), disc_Λ(F_res). The corner-fold path is
    /// (w, w, −1, −4N); the interior-(1,1) identity path is (1, 1, +1, 0).</summary>
    private static (GaussianInteger[] Residual, GaussianInteger[] Composed,
        GaussianInteger ResultantValue, GaussianInteger DiscriminantValue)
        ExactSampleCore(int n, bool rOdd, int q0, int tWKet, int tWBra,
            GaussianInteger composeA, GaussianInteger composeB)
    {
        var blk = rOdd
            ? F89PathKSeDeBlock.BuildTwoTimesROddBlock(q0, n)
            : F89PathKSeDeBlock.BuildTwoTimesSymBlock(q0, n);
        var chFull = GaussianMatrixCharpoly.Characteristic(blk);
        var at = AtFactorAt(n, rOdd, q0);
        var (res, rem) = GaussianPolynomial.DivMod(chFull, at);
        if (rem.Length != 0)
            throw new InvalidOperationException($"AT factor does not divide the block charpoly at q0={q0}.");

        var (aTgt, cTgt) = TargetPencil(n, tWKet, tWBra);
        int dim = aTgt.GetLength(0);
        var tgt = new GaussianInteger[dim, dim];
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
                tgt[i, j] = aTgt[i, j] + new GaussianInteger(q0, 0) * cTgt[i, j];
        var f2tgt = GaussianMatrixCharpoly.Characteristic(tgt);
        var g = GaussianPolynomial.ComposeLinear(f2tgt, composeA, composeB);

        return (res, g, GaussianPolynomial.Resultant(res, g), GaussianPolynomial.Discriminant(res));
    }

    /// <summary>The exact bivariate AT factor evaluated at integer q0, over Z[i]: per rate sector,
    /// Σ_m p_m·(i·q0)^{d−m}·(Λ−r0)^m (the q-linear strands Λ = r0 + i·q·κ_j on the q-independent
    /// invariant subspace), multiplied over sectors. Monic of the AT degree; at q0=2 it reproduces
    /// <see cref="F89AtFactorReconstruction.ForPathK"/> exactly (pinned by the gate).</summary>
    public static GaussianInteger[] AtFactorAt(int n, bool rOdd, int q0)
    {
        var sectors = F89AtFactorReconstruction.ClearedAtSectors(n - 1, rOdd);
        var at = new[] { GaussianInteger.One };
        foreach (var s in sectors)
        {
            int d = s.KCharpoly.Length - 1;
            var polyX = new GaussianInteger[d + 1];
            var iq = new GaussianInteger(0, q0);
            for (int m = 0; m <= d; m++)
            {
                var scale = GaussianInteger.One;
                for (int e = 0; e < d - m; e++) scale *= iq;
                polyX[m] = new GaussianInteger(s.KCharpoly[m], 0) * scale;
            }
            var sector = GaussianPolynomial.ComposeLinear(
                polyX, GaussianInteger.One, new GaussianInteger(-s.Rate, 0));
            at = GaussianPolynomial.Multiply(at, sector);
        }
        return at;
    }

    /// <summary>The exact q → ∞ leading forms and their two resultant-type invariants: f_res(x) = the
    /// residual part of charpoly(C_block) (the AT strands' hopping part charpoly(C|W) = Σ_m p_m·i^{d−m}·x^m
    /// divided out), g(x) = charpoly(C2_corner)(−x), and Res_x(f_res, g), disc_x(f_res). Each branch obeys
    /// Λ_i(q) = c_i·q + O(1) with c_i a root of its leading form, so the multiplicity collisions of these
    /// forms give the proven q-degree bounds of R and D (guard G2).</summary>
    public static (GaussianInteger[] ResidualLeadingForm, GaussianInteger[] FoldedCornerLeadingForm,
        GaussianInteger ResultantLeadingCoeff, GaussianInteger DiscriminantLeadingCoeff)
        LeadingForms(int n, bool rOdd)
        => LeadingFormsCore(n, rOdd, (n + 1) / 2 + 1, (n + 1) / 2 + 1, new GaussianInteger(-1, 0));

    /// <summary>Leading forms generalised over the target block and the composeA sign (composeB drops in the
    /// q→∞ leading form): g(x) = charpoly(C_target)(composeA·x). Corner-fold composeA=−1; interior identity +1.</summary>
    private static (GaussianInteger[] ResidualLeadingForm, GaussianInteger[] ComposedLeadingForm,
        GaussianInteger ResultantLeadingCoeff, GaussianInteger DiscriminantLeadingCoeff)
        LeadingFormsCore(int n, bool rOdd, int tWKet, int tWBra, GaussianInteger composeA)
    {
        var (_, cBlk) = BlockPencil(n, rOdd);
        var fFull = GaussianMatrixCharpoly.Characteristic(cBlk);

        var sectors = F89AtFactorReconstruction.ClearedAtSectors(n - 1, rOdd);
        var fAt = new[] { GaussianInteger.One };
        foreach (var s in sectors)
        {
            int d = s.KCharpoly.Length - 1;
            var polyX = new GaussianInteger[d + 1];
            for (int m = 0; m <= d; m++)
            {
                var unit = GaussianInteger.One;
                for (int e = 0; e < d - m; e++) unit *= GaussianInteger.I;
                polyX[m] = new GaussianInteger(s.KCharpoly[m], 0) * unit;
            }
            fAt = GaussianPolynomial.Multiply(fAt, polyX);
        }
        var (fRes, rem) = GaussianPolynomial.DivMod(fFull, fAt);
        if (rem.Length != 0)
            throw new InvalidOperationException("the AT leading form does not divide charpoly(C_block).");

        var (_, cTgt) = TargetPencil(n, tWKet, tWBra);
        var fTgt = GaussianMatrixCharpoly.Characteristic(cTgt);
        var g = GaussianPolynomial.ComposeLinear(fTgt, composeA, GaussianInteger.Zero);

        return (fRes, g,
            GaussianPolynomial.Resultant(fRes, g),
            GaussianPolynomial.Discriminant(fRes));
    }

    // ---- exact integer pencils ----

    /// <summary>The ×2-cleared sector block as an exact integer pencil (A, C), block(q) = A + q·C, built
    /// from the q0 = 0, 1 builders; linearity is re-verified at q0 = 2 (throws on violation).</summary>
    private static (GaussianInteger[,] A, GaussianInteger[,] C) BlockPencil(int n, bool rOdd)
    {
        var b0 = rOdd ? F89PathKSeDeBlock.BuildTwoTimesROddBlock(0, n) : F89PathKSeDeBlock.BuildTwoTimesSymBlock(0, n);
        var b1 = rOdd ? F89PathKSeDeBlock.BuildTwoTimesROddBlock(1, n) : F89PathKSeDeBlock.BuildTwoTimesSymBlock(1, n);
        var b2 = rOdd ? F89PathKSeDeBlock.BuildTwoTimesROddBlock(2, n) : F89PathKSeDeBlock.BuildTwoTimesSymBlock(2, n);
        int dim = b0.GetLength(0);
        var c = new GaussianInteger[dim, dim];
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
            {
                c[i, j] = b1[i, j] - b0[i, j];
                if (!(b0[i, j] + c[i, j] + c[i, j]).Equals(b2[i, j]))
                    throw new InvalidOperationException("the sector block is not linear in q0.");
            }
        return (b0, c);
    }

    /// <summary>The ×2-cleared corner block (p_c+1, p_c+1) as an exact integer pencil (A, C).</summary>
    private static (GaussianInteger[,] A, GaussianInteger[,] C) CornerPencil(int n)
        => TargetPencil(n, (n + 1) / 2 + 1, (n + 1) / 2 + 1);

    /// <summary>The ×2-cleared (wKet, wBra) TARGET block as an exact integer pencil (A, C): twice the
    /// <see cref="WeightCoherenceBlock"/> at q = 0 and the q-slope, entries checked to be exact Gaussian
    /// integers (they are: diagonal −4·n_diff, hops ±4i). Serves the corner (p_c+1, p_c+1) for remainder R1's
    /// fold certificate AND the interior cores (e.g. (1,1)) for remainder R4's identity-composition exclusion
    /// certificate (the braid eigenvalue λ_A, an F_res root, is not carried by the non-member block).</summary>
    private static (GaussianInteger[,] A, GaussianInteger[,] C) TargetPencil(int n, int wKet, int wBra)
    {
        var l0 = WeightCoherenceBlock.Build(n, wKet, wBra, Complex.Zero);
        var l1 = WeightCoherenceBlock.Build(n, wKet, wBra, Complex.One);
        int dim = l0.GetLength(0);
        var a = new GaussianInteger[dim, dim];
        var c = new GaussianInteger[dim, dim];
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
            {
                a[i, j] = ExactGaussian(2.0 * l0[i, j]);
                c[i, j] = ExactGaussian(2.0 * (l1[i, j] - l0[i, j]));
            }
        return (a, c);
    }

    /// <summary>Checked Complex → Z[i] for the exactly-integer entries of the small coherence blocks.</summary>
    private static GaussianInteger ExactGaussian(Complex v)
    {
        double re = Math.Round(v.Real), im = Math.Round(v.Imaginary);
        if (Math.Abs(v.Real - re) > 1e-9 || Math.Abs(v.Imaginary - im) > 1e-9)
            throw new InvalidOperationException($"expected an exact Gaussian-integer entry, got {v}.");
        return new GaussianInteger(new BigInteger(re), new BigInteger(im));
    }

    // ---- the exact bivariate layer: F(Λ, q) as Λ-indexed q-polynomials over Z[i] ----

    private static GaussianInteger[] AddQ(GaussianInteger[] a, GaussianInteger[] b)
    {
        int len = Math.Max(a.Length, b.Length);
        var r = new GaussianInteger[len];
        for (int i = 0; i < len; i++)
        {
            var va = i < a.Length ? a[i] : GaussianInteger.Zero;
            var vb = i < b.Length ? b[i] : GaussianInteger.Zero;
            r[i] = va + vb;
        }
        return GaussianPolynomial.Trim(r);
    }

    private static GaussianInteger[] NegQ(GaussianInteger[] a)
    {
        var r = new GaussianInteger[a.Length];
        for (int i = 0; i < a.Length; i++) r[i] = -a[i];
        return r;
    }

    private static GaussianInteger[] SubQ(GaussianInteger[] a, GaussianInteger[] b) => AddQ(a, NegQ(b));

    private static GaussianInteger[] ScaleQ(GaussianInteger[] a, GaussianInteger s)
    {
        var r = new GaussianInteger[a.Length];
        for (int i = 0; i < a.Length; i++) r[i] = a[i] * s;
        return GaussianPolynomial.Trim(r);
    }

    /// <summary>det(ΛI − (A + q·C)) over Z[i][q] by the division-free Berkowitz algorithm with q-poly
    /// entries (the GaussianMatrixCharpoly algorithm, elements lifted from Z[i] to Z[i][q]). Returns the
    /// bivariate charpoly as Λ-indexed q-polynomials, monic in Λ.</summary>
    private static GaussianInteger[][] PencilCharpolyBivariate(GaussianInteger[,] aM, GaussianInteger[,] cM)
    {
        int n = aM.GetLength(0);
        GaussianInteger[] E(int i, int j) => GaussianPolynomial.Trim(new[] { aM[i, j], cM[i, j] });

        var one = new[] { GaussianInteger.One };
        var p = new[] { one };                              // highest-first in Λ during the build
        for (int i = n; i >= 1; i--)
        {
            int m = n - i + 1;
            int sub = m - 1;

            var q = new GaussianInteger[m + 1][];
            q[0] = one;
            q[1] = NegQ(E(i - 1, i - 1));
            if (sub > 0)
            {
                var vec = new GaussianInteger[sub][];
                for (int row = 0; row < sub; row++) vec[row] = E(i + row, i - 1);
                for (int k = 0; k < sub; k++)
                {
                    var term = Array.Empty<GaussianInteger>();
                    for (int j = 0; j < sub; j++)
                        term = AddQ(term, GaussianPolynomial.Multiply(E(i - 1, i + j), vec[j]));
                    q[k + 2] = NegQ(term);
                    if (k < sub - 1)
                    {
                        var next = new GaussianInteger[sub][];
                        for (int row = 0; row < sub; row++)
                        {
                            var s = Array.Empty<GaussianInteger>();
                            for (int col = 0; col < sub; col++)
                                s = AddQ(s, GaussianPolynomial.Multiply(E(i + row, i + col), vec[col]));
                            next[row] = s;
                        }
                        vec = next;
                    }
                }
            }
            var pNew = new GaussianInteger[m + 1][];
            for (int rr = 0; rr <= m; rr++)
            {
                var s = Array.Empty<GaussianInteger>();
                for (int c = 0; c < p.Length; c++)
                {
                    int d = rr - c;
                    if (d >= 0 && d <= m) s = AddQ(s, GaussianPolynomial.Multiply(q[d], p[c]));
                }
                pNew[rr] = s;
            }
            p = pNew;
        }
        Array.Reverse(p);                                   // highest-first → lowest-first in Λ
        return p;
    }

    /// <summary>The exact bivariate AT factor ∏_sectors Σ_m p_m·i^{d−m}·(Λ−r0)^m·q^{d−m} ∈ Z[i][Λ, q],
    /// monic in Λ (per sector by Horner in Λ over the q-monomial coefficients).</summary>
    private static GaussianInteger[][] BivariateAtFactor(IReadOnlyList<AtSector> sectors)
    {
        var at = new[] { new[] { GaussianInteger.One } };
        foreach (var s in sectors)
        {
            int d = s.KCharpoly.Length - 1;
            // coefficient of (Λ−r0)^m is p_m·i^{d−m}·q^{d−m}, a q-monomial
            GaussianInteger[] C(int m)
            {
                var unit = GaussianInteger.One;
                for (int e = 0; e < d - m; e++) unit *= GaussianInteger.I;
                var mono = new GaussianInteger[d - m + 1];
                for (int k = 0; k < mono.Length; k++) mono[k] = GaussianInteger.Zero;
                mono[d - m] = new GaussianInteger(s.KCharpoly[m], 0) * unit;
                return mono;
            }
            var negR0 = new GaussianInteger(-s.Rate, 0);
            var acc = new[] { C(d) };                       // Horner in x = Λ − r0
            for (int m = d - 1; m >= 0; m--)
            {
                var next = new GaussianInteger[acc.Length + 1][];
                next[0] = Array.Empty<GaussianInteger>();
                for (int l = 0; l < acc.Length; l++) next[l + 1] = acc[l];               // ·Λ
                for (int l = 0; l < acc.Length; l++) next[l] = AddQ(next[l], ScaleQ(acc[l], negR0)); // −r0·
                next[0] = AddQ(next[0], C(m));
                acc = next;
            }
            at = BivariateMultiply(at, acc);
        }
        return at;
    }

    private static GaussianInteger[][] BivariateMultiply(GaussianInteger[][] a, GaussianInteger[][] b)
    {
        var r = new GaussianInteger[a.Length + b.Length - 1][];
        for (int l = 0; l < r.Length; l++) r[l] = Array.Empty<GaussianInteger>();
        for (int i = 0; i < a.Length; i++)
            for (int j = 0; j < b.Length; j++)
                r[i + j] = AddQ(r[i + j], GaussianPolynomial.Multiply(a[i], b[j]));
        return r;
    }

    /// <summary>Λ-long division of bivariate polynomials by a divisor MONIC in Λ; the remainder must
    /// vanish identically (throws otherwise): the AT factor divides the block charpoly exactly.</summary>
    private static GaussianInteger[][] BivariateDivideExact(GaussianInteger[][] dividend, GaussianInteger[][] divisor)
    {
        int dd = divisor.Length - 1;
        if (GaussianPolynomial.Degree(divisor[dd]) != 0 || !divisor[dd][0].Equals(GaussianInteger.One))
            throw new ArgumentException("bivariate divisor must be monic in Λ.");
        var rem = new GaussianInteger[dividend.Length][];
        for (int l = 0; l < rem.Length; l++) rem[l] = dividend[l];
        int da = dividend.Length - 1;
        var quo = new GaussianInteger[da - dd + 1][];
        for (int i = da; i >= dd; i--)
        {
            var coef = rem[i];
            quo[i - dd] = coef;
            if (GaussianPolynomial.Degree(coef) >= 0)
                for (int j = 0; j <= dd; j++)
                    rem[i - dd + j] = SubQ(rem[i - dd + j], GaussianPolynomial.Multiply(coef, divisor[j]));
        }
        for (int l = 0; l < dd; l++)
            if (GaussianPolynomial.Degree(rem[l]) >= 0)
                throw new InvalidOperationException("bivariate AT division left a remainder.");
        return quo;
    }

    /// <summary>f(αΛ + β, q) by Horner in Λ (α, β Gaussian-integer constants): the fold composition.</summary>
    private static GaussianInteger[][] BivariateComposeLinearLambda(
        GaussianInteger[][] f, GaussianInteger alpha, GaussianInteger beta)
    {
        int d = f.Length - 1;
        var acc = new[] { f[d] };
        for (int k = d - 1; k >= 0; k--)
        {
            var next = new GaussianInteger[acc.Length + 1][];
            next[0] = Array.Empty<GaussianInteger>();
            for (int l = 0; l < acc.Length; l++) next[l + 1] = ScaleQ(acc[l], alpha);
            for (int l = 0; l < acc.Length; l++) next[l] = AddQ(next[l], ScaleQ(acc[l], beta));
            next[0] = AddQ(next[0], f[k]);
            acc = next;
        }
        return acc;
    }

    /// <summary>The bivariate polynomial at integer q0, as an exact Λ-polynomial over Z[i].</summary>
    private static GaussianInteger[] EvalBivariateAtQ(GaussianInteger[][] f, int q0)
    {
        var r = new GaussianInteger[f.Length];
        var x = new GaussianInteger(q0, 0);
        for (int l = 0; l < f.Length; l++) r[l] = GaussianPolynomial.Evaluate(f[l], x);
        return GaussianPolynomial.Trim(r);
    }

    private static void AssertPolyEqual(GaussianInteger[] a, GaussianInteger[] b, string what)
    {
        if (a.Length != b.Length)
            throw new InvalidOperationException($"{what}: degree mismatch {a.Length - 1} vs {b.Length - 1}.");
        for (int i = 0; i < a.Length; i++)
            if (!a[i].Equals(b[i]))
                throw new InvalidOperationException($"{what}: coefficient {i} mismatch.");
    }

    /// <summary>The 1-norm surrogate Σ(|Re| + |Im|) over a q-poly's coefficients (≥ the complex moduli;
    /// submultiplicative), feeding the Hadamard row-product bound on ‖R‖.</summary>
    private static BigInteger HNormQ(GaussianInteger[] poly)
    {
        BigInteger s = BigInteger.Zero;
        foreach (var c in poly) s += BigInteger.Abs(c.Re) + BigInteger.Abs(c.Im);
        return s;
    }

    // ---- mod-p reduction (i ↦ r, r² ≡ −1; the OcticGaloisCertificate idiom) ----

    private static int BigIntMod(BigInteger a, int p)
    {
        int m = (int)(a % p);
        return m < 0 ? m + p : m;
    }

    private static int ReduceScalar(GaussianInteger g, int p, int r) => BigIntMod(g.Re + g.Im * r, p);

    private static int[][] ReduceBivariate(GaussianInteger[][] f, int p, int r)
    {
        var res = new int[f.Length][];
        for (int l = 0; l < f.Length; l++)
        {
            res[l] = new int[Math.Max(f[l].Length, 1)];
            for (int k = 0; k < f[l].Length; k++) res[l][k] = ReduceScalar(f[l][k], p, r);
        }
        return res;
    }

    /// <summary>The reduced bivariate polynomial at sample q0, as a Λ-polynomial over F_p (per Λ-power a
    /// Horner evaluation of the q-poly).</summary>
    private static int[] EvalBivariateModP(int[][] fp, int q0, int p)
    {
        var r = new int[fp.Length];
        for (int l = 0; l < fp.Length; l++)
        {
            long acc = 0;
            for (int k = fp[l].Length - 1; k >= 0; k--) acc = (acc * q0 + fp[l][k]) % p;
            r[l] = (int)acc;
        }
        return TrimP(r);
    }

    // ---- F_p[x] arithmetic (int[] lowest-first; mirrors OcticGaloisCertificate's private idiom) ----

    private static int DegP(int[] a)
    {
        for (int i = a.Length - 1; i >= 0; i--) if (a[i] != 0) return i;
        return -1;
    }

    private static int[] TrimP(int[] a)
    {
        int d = DegP(a);
        if (d < 0) return new[] { 0 };
        var r = new int[d + 1];
        Array.Copy(a, r, d + 1);
        return r;
    }

    private static int InvModP(int a, int p)
    {
        long b = (long)a % p;                               // long BEFORE the shift: a % p + p can
        if (b < 0) b += p;                                  // exceed int.MaxValue for p ≈ 2³⁰
        long result = 1;
        long e = p - 2;
        while (e > 0) { if ((e & 1) == 1) result = result * b % p; b = b * b % p; e >>= 1; }
        return (int)result;
    }

    private static int PowModP(long a, long e, int p)
    {
        long result = 1, b = ((a % p) + p) % p;
        while (e > 0) { if ((e & 1) == 1) result = result * b % p; b = b * b % p; e >>= 1; }
        return (int)result;
    }

    private static int[] MulModP(int[] a, int[] b, int p)
    {
        int da = DegP(a), db = DegP(b);
        if (da < 0 || db < 0) return new[] { 0 };
        var r = new long[da + db + 1];
        for (int i = 0; i <= da; i++)
            if (a[i] != 0)
                for (int j = 0; j <= db; j++) r[i + j] = (r[i + j] + (long)a[i] * b[j]) % p;
        var ri = new int[r.Length];
        for (int k = 0; k < r.Length; k++) ri[k] = (int)r[k];
        return TrimP(ri);
    }

    /// <summary>Long division by a divisor with invertible leading coefficient.</summary>
    private static (int[] Quo, int[] Rem) DivModP(int[] dividend, int[] divisor, int p)
    {
        int df = DegP(divisor);
        if (df < 0) throw new DivideByZeroException("polynomial divisor is zero.");
        int da = DegP(dividend);
        if (da < df) return (new[] { 0 }, TrimP(dividend));
        var rem = new int[da + 1];
        Array.Copy(dividend, rem, Math.Min(dividend.Length, rem.Length));
        var quo = new int[da - df + 1];
        int invLead = InvModP(divisor[df], p);
        for (int i = da; i >= df; i--)
        {
            int coef = (int)((long)rem[i] * invLead % p);
            quo[i - df] = coef;
            if (coef != 0)
                for (int j = 0; j <= df; j++)
                    rem[i - df + j] = (int)((((long)rem[i - df + j] - (long)coef * divisor[j]) % p + p) % p);
        }
        return (TrimP(quo), TrimP(rem));
    }

    private static int[] MonicModP(int[] a, int p)
    {
        int d = DegP(a);
        if (d < 0) return new[] { 0 };
        int inv = InvModP(a[d], p);
        var r = new int[d + 1];
        for (int i = 0; i <= d; i++) r[i] = (int)((long)a[i] * inv % p);
        return r;
    }

    private static int[] GcdModP(int[] a, int[] b, int p)
    {
        var x = TrimP(a);
        var y = TrimP(b);
        while (DegP(y) >= 0)
        {
            var (_, rem) = DivModP(x, y, p);
            x = y;
            y = rem;
        }
        return MonicModP(x, p);
    }

    private static int[] DerivativeModP(int[] a, int p)
    {
        if (DegP(a) <= 0) return new[] { 0 };
        var r = new int[a.Length - 1];
        for (int i = 1; i < a.Length; i++) r[i - 1] = (int)((long)a[i] * (i % p) % p);
        return TrimP(r);
    }

    /// <summary>Res(a, b) over F_p by the Euclidean remainder sequence: Res(f, g) =
    /// (−1)^{deg f·deg g}·lc(g)^{deg f − deg r}·Res(g, r) with r = f mod g, terminating at
    /// Res(f, c) = c^{deg f}. Chosen over a Sylvester-determinant
    /// route for speed (O(n²)) and simplicity; validated against the exact Z[i] resultant in the gate.</summary>
    private static int ResultantModP(int[] a, int[] b, int p)
    {
        var f = TrimP(a);
        var g = TrimP(b);
        int df = DegP(f), dg = DegP(g);
        if (df < 0 || dg < 0) return 0;
        long res = 1;
        while (dg > 0)
        {
            var (_, r) = DivModP(f, g, p);
            int dr = DegP(r);
            if (dr >= dg)
                throw new InvalidOperationException(
                    $"Euclid resultant: remainder degree {dr} did not drop below {dg} (p={p}): non-field arithmetic?");
            if ((df & 1) == 1 && (dg & 1) == 1) res = (p - res) % p;        // (−1)^{df·dg}
            res = res * PowModP(g[dg], df - (dr < 0 ? 0 : dr), p) % p;
            if (dr < 0) return 0;                                           // g | f: shared root
            f = g;
            g = r;
            df = dg;
            dg = dr;
        }
        return (int)(res * PowModP(g[0], df, p) % p);
    }

    /// <summary>disc(f) = (−1)^{n(n−1)/2}·Res(f, f′) for monic f over F_p (matches
    /// <see cref="GaussianPolynomial.Discriminant"/>'s convention, so cross-checks are sign-exact).</summary>
    private static int DiscriminantModP(int[] f, int p)
    {
        int n = DegP(f);
        if (n < 1) return 0;
        int res = ResultantModP(f, DerivativeModP(f, p), p);
        return ((long)n * (n - 1) / 2 % 2 == 0) ? res : (p - res) % p;
    }

    private static int EvalModP(int[] poly, int x, int p)
    {
        int d = DegP(poly);
        if (d < 0) return 0;
        long acc = poly[d];
        for (int k = d - 1; k >= 0; k--) acc = (acc * x + poly[k]) % p;
        return (int)acc;
    }

    /// <summary>Newton interpolation at the unit-spaced nodes x = 0..k−1 over F_p (exact; the divided
    /// differences at unit spacing divide by the level only).</summary>
    private static int[] InterpolateAtIntegerNodes(int[] values, int p)
    {
        int k = values.Length;
        var dd = new long[k];
        for (int i = 0; i < k; i++) dd[i] = values[i];
        for (int lvl = 1; lvl < k; lvl++)
        {
            long invLvl = InvModP(lvl, p);
            for (int i = k - 1; i >= lvl; i--)
                dd[i] = ((dd[i] - dd[i - 1]) % p + p) % p * invLvl % p;
        }
        var poly = new long[k];
        poly[0] = dd[k - 1];
        int deg = 0;
        for (int i = k - 2; i >= 0; i--)
        {
            long negI = ((-(long)i) % p + p) % p;
            for (int c = deg + 1; c >= 1; c--) poly[c] = (poly[c - 1] + poly[c] * negI) % p;
            poly[0] = (poly[0] * negI + dd[i]) % p;
            deg++;
        }
        var res = new int[k];
        for (int c = 0; c < k; c++) res[c] = (int)poly[c];
        return TrimP(res);
    }

    /// <summary>Newton interpolation at arbitrary DISTINCT nodes over F_p; the unit-spaced
    /// <see cref="InterpolateAtIntegerNodes"/> is the special case that needs no inverses beyond
    /// the level. One BATCHED inversion per level (prefix products, one <see cref="InvModP"/>,
    /// back-substitution), so the cost stays O(k²) multiplications and O(k) inversions rather
    /// than the O(k²) inversions the naive divided differences would ask for.</summary>
    private static int[] InterpolateAtNodes(int[] nodes, int[] values, int p)
    {
        int k = values.Length;
        if (nodes.Length != k)
            throw new ArgumentException($"nodes ({nodes.Length}) and values ({k}) disagree in length.");
        var dd = new long[k];
        for (int i = 0; i < k; i++) dd[i] = values[i];

        var den = new long[k];
        var pre = new long[k];
        for (int lvl = 1; lvl < k; lvl++)
        {
            long acc = 1;
            for (int i = lvl; i < k; i++)
            {
                long d = ((long)nodes[i] - nodes[i - lvl]) % p;
                if (d < 0) d += p;
                if (d == 0)
                    throw new InvalidOperationException($"interpolation nodes {i - lvl} and {i} coincide mod {p}.");
                den[i] = d;
                pre[i] = acc;
                acc = acc * d % p;
            }
            long inv = InvModP((int)acc, p);
            for (int i = k - 1; i >= lvl; i--)
            {
                long invD = inv * pre[i] % p;               // = 1 / den[i], by the prefix product
                inv = inv * den[i] % p;
                dd[i] = ((dd[i] - dd[i - 1]) % p + p) % p * invD % p;
            }
        }

        var poly = new long[k];
        poly[0] = dd[k - 1];
        int deg = 0;
        for (int i = k - 2; i >= 0; i--)
        {
            long negX = ((-(long)nodes[i]) % p + p) % p;
            for (int c = deg + 1; c >= 1; c--) poly[c] = (poly[c - 1] + poly[c] * negX) % p;
            poly[0] = (poly[0] * negX + dd[i]) % p;
            deg++;
        }
        var res = new int[k];
        for (int c = 0; c < k; c++) res[c] = (int)poly[c];
        return TrimP(res);
    }

    /// <summary>One prime's pass of the D-device: reduce F_res mod p, sample the Λ-discriminant at
    /// every q0, interpolate D (halved through E(q²) where licensed), verify at the trailing nodes,
    /// and return the prime's (deg_q D, v_q D, q-stripped D). Null means a CONTENT prime, D ≡ 0
    /// mod p, which carries no information about either invariant and is left uncounted.
    /// <para>Pure and allocation-local by construction, which is what makes the caller's
    /// <c>Parallel.For</c> legitimate: the only shared state is <paramref name="fRes"/>, read and
    /// never written. Keep it that way. Anything added here that touches shared mutable state, a
    /// cache or a logger included, silently turns a certificate into a race.</para></summary>
    private static (int P, int DegD, int VD, int[] DStrip)? SampleDiscAtPrime(
        GaussianInteger[][] fRes, int p, int root, int resDeg, int samples, int nodeCount, bool halve)
    {
        if (halve) AssertSquareAbscissaeFit(p, samples - 1);
        var fResP = ReduceBivariate(fRes, p, root);
        var dSamp = new int[samples];
        for (int q0 = 0; q0 < samples; q0++)
        {
            var fL = EvalBivariateModP(fResP, q0, p);
            if (DegP(fL) != resDeg || fL[resDeg] != 1)
                throw new InvalidOperationException($"reduced residual not monic of degree {resDeg} at q0={q0}.");
            dSamp[q0] = DiscriminantModP(fL, p);
        }
        var dNodes = new int[nodeCount];
        Array.Copy(dSamp, dNodes, nodeCount);
        int[] dp;
        if (halve)
        {
            var abscissae = new int[nodeCount];
            for (int i = 0; i < nodeCount; i++) abscissae[i] = (int)((long)i * i % p);
            dp = SpreadEven(InterpolateAtNodes(abscissae, dNodes, p));
        }
        else dp = InterpolateAtIntegerNodes(dNodes, p);
        for (int q0 = nodeCount; q0 < samples; q0++)
            if (EvalModP(dp, q0, p) != dSamp[q0])
                throw new InvalidOperationException($"D interpolant fails verification at q0={q0} (p={p}).");

        if (DegP(dp) < 0) return null;
        int vDp = QValuation(dp);
        return (p, DegP(dp), vDp, StripQ(dp, vDp));
    }

    /// <summary>E(m) ↦ D(q) = E(q²): the even-support spread that turns the halved interpolant
    /// back into the object every downstream reader (valuation, q-strip, squarefree layers)
    /// already knows how to read, so the halving changes the sampling and NOTHING else.</summary>
    private static int[] SpreadEven(int[] e)
    {
        int d = DegP(e);
        if (d < 0) return new[] { 0 };
        var r = new int[2 * d + 1];
        for (int j = 0; j <= d; j++) r[2 * j] = e[j];
        return r;
    }

    private static int QValuation(int[] poly)
    {
        if (DegP(poly) < 0) return 0;
        int v = 0;
        while (v < poly.Length && poly[v] == 0) v++;
        return v;
    }

    private static int[] StripQ(int[] poly, int v)
    {
        if (v == 0) return TrimP(poly);
        int d = DegP(poly);
        if (d < v) return new[] { 0 };
        var r = new int[d - v + 1];
        Array.Copy(poly, v, r, 0, d - v + 1);
        return r;
    }

    /// <summary>The squarefree layers u_1, u_2, … of f over F_p (u_k = roots of multiplicity exactly k,
    /// each once; f = lc·∏u_k^k), by the gcd chain. Valid because every multiplicity here is far below p.</summary>
    private static List<int[]> FpSquarefreeLayers(int[] f, int p)
    {
        var chain = new List<int[]>();                      // f_0 = f, f_i = gcd(f_{i−1}, f_{i−1}')
        var cur = MonicModP(f, p);
        chain.Add(cur);
        while (DegP(cur) > 0)
        {
            cur = GcdModP(cur, DerivativeModP(cur, p), p);
            chain.Add(cur);
        }
        var rad = new List<int[]>();                        // rad_i = f_{i−1}/f_i = roots of mult ≥ i
        for (int i = 1; i < chain.Count; i++)
            rad.Add(DivModP(chain[i - 1], chain[i], p).Quo);
        var layers = new List<int[]>();                     // u_k = rad_k / rad_{k+1}
        for (int k = 0; k < rad.Count; k++)
            layers.Add(k + 1 < rad.Count ? DivModP(rad[k], rad[k + 1], p).Quo : rad[k]);
        return layers;
    }

    // ---- prime generation (p ≡ 1 mod 4, p ≥ 2³⁰; int-ranged so the long-product F_p ops stay exact) ----

    private static bool IsPrime(long n)
    {
        if (n < 2) return false;
        foreach (long sp in new long[] { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 })
        {
            if (n == sp) return true;
            if (n % sp == 0) return false;
        }
        long d = n - 1;
        int s = 0;
        while ((d & 1) == 0) { d >>= 1; s++; }
        foreach (long a in new long[] { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 })
        {
            long x = PowModLong(a, d, n);
            if (x == 1 || x == n - 1) continue;
            bool witness = true;
            for (int i = 1; i < s; i++)
            {
                x = MulModLong(x, x, n);
                if (x == n - 1) { witness = false; break; }
            }
            if (witness) return false;
        }
        return true;
    }

    private static long MulModLong(long a, long b, long m) => (long)((UInt128)a * (UInt128)b % (UInt128)m);

    private static long PowModLong(long a, long e, long m)
    {
        long result = 1, b = a % m;
        while (e > 0) { if ((e & 1) == 1) result = MulModLong(result, b, m); b = MulModLong(b, b, m); e >>= 1; }
        return result;
    }

    /// <summary>r with r² ≡ −1 (mod p) for a split prime p ≡ 1 (mod 4): a^{(p−1)/4} for the first
    /// quadratic non-residue a (checked; a residue yields ±1 and is skipped).</summary>
    private static int? SqrtMinusOneFast(int p)
    {
        for (int a = 2; a < 1000; a++)
        {
            long r = PowModP(a, (p - 1) / 4, p);
            if (r * r % p == p - 1) return (int)r;
        }
        return null;
    }

    // ---- the O(1)-branch premise check ----

    /// <summary>Asserts the hopping direction is C = i·K with K REAL and weighted-symmetric
    /// (K[r,c]·d_r = K[c,r]·d_c, positive metric weights d; null = uniform). Then C is anti-self-adjoint
    /// in the positive-definite metric diag(d), hence NORMAL there: diagonalizable with self-adjoint
    /// compressions, so every eigenvalue branch of A + qC is Λ = c·q + O(1) with NO Puiseux growth
    /// (Bauer-Fike), the premise of the multiplicity-adjusted degree bounds. Both parts are checked:
    /// mere weighted COMPLEX symmetry would not imply diagonalizability (a complex-symmetric matrix can
    /// be defective), the purely-imaginary part is load-bearing (2026-07-03 adversarial review).</summary>
    private static void AssertHoppingSymmetry(GaussianInteger[,] c, int[]? weights)
    {
        int d = c.GetLength(0);
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
            {
                if (!c[i, j].Re.IsZero)
                    throw new InvalidOperationException(
                        "the hopping direction is not purely imaginary: the normality premise fails.");
                var lhs = weights is null ? c[i, j] : c[i, j] * new GaussianInteger(weights[i], 0);
                var rhs = weights is null ? c[j, i] : c[j, i] * new GaussianInteger(weights[j], 0);
                if (!lhs.Equals(rhs))
                    throw new InvalidOperationException(
                        "the hopping direction is not (weighted-)symmetric: the O(1) branch-bound premise fails.");
            }
    }

    // ---- exact Q(i)[x] arithmetic for the leading forms' multiplicity structure (small degrees) ----

    private readonly struct Qc
    {
        public readonly BigRational Re;
        public readonly BigRational Im;
        public Qc(BigRational re, BigRational im) { Re = re; Im = im; }
        public static Qc Zero => new(new BigRational(0), new BigRational(0));
        public bool IsZero => Re.IsZero && Im.IsZero;
        public static Qc operator +(Qc a, Qc b) => new(a.Re + b.Re, a.Im + b.Im);
        public static Qc operator -(Qc a, Qc b) => new(a.Re - b.Re, a.Im - b.Im);
        public static Qc operator *(Qc a, Qc b)
            => new(a.Re * b.Re - a.Im * b.Im, a.Re * b.Im + a.Im * b.Re);
        public static Qc operator /(Qc a, Qc b)
        {
            var norm = b.Re * b.Re + b.Im * b.Im;
            return new((a.Re * b.Re + a.Im * b.Im) / norm, (a.Im * b.Re - a.Re * b.Im) / norm);
        }
    }

    private static Qc[] FromZi(GaussianInteger[] p)
    {
        var r = new Qc[Math.Max(p.Length, 1)];
        for (int k = 0; k < r.Length; k++) r[k] = Qc.Zero;
        for (int k = 0; k < p.Length; k++) r[k] = new Qc(new BigRational(p[k].Re), new BigRational(p[k].Im));
        return r;
    }

    private static int QcDeg(Qc[] a)
    {
        for (int i = a.Length - 1; i >= 0; i--) if (!a[i].IsZero) return i;
        return -1;
    }

    private static Qc[] QcMonic(Qc[] a)
    {
        int d = QcDeg(a);
        if (d < 0) return new[] { Qc.Zero };
        var r = new Qc[d + 1];
        for (int i = 0; i <= d; i++) r[i] = a[i] / a[d];
        return r;
    }

    private static Qc[] QcDerivative(Qc[] a)
    {
        int d = QcDeg(a);
        if (d <= 0) return new[] { Qc.Zero };
        var r = new Qc[d];
        for (int i = 1; i <= d; i++)
            r[i - 1] = a[i] * new Qc(new BigRational(i), new BigRational(0));
        return r;
    }

    private static Qc[] QcRem(Qc[] dividend, Qc[] divisor)
    {
        int df = QcDeg(divisor);
        if (df < 0) throw new DivideByZeroException("polynomial divisor is zero.");
        var rem = new Qc[Math.Max(QcDeg(dividend) + 1, 1)];
        for (int i = 0; i < rem.Length; i++) rem[i] = i < dividend.Length ? dividend[i] : Qc.Zero;
        for (int i = QcDeg(rem); i >= df; i--)
        {
            var coef = rem[i] / divisor[df];
            if (coef.IsZero) continue;
            for (int j = 0; j <= df; j++)
                rem[i - df + j] -= coef * divisor[j];
        }
        return rem;
    }

    private static Qc[] QcDivExact(Qc[] dividend, Qc[] divisor)
    {
        int df = QcDeg(divisor), da = QcDeg(dividend);
        if (da < df) throw new InvalidOperationException("exact division with smaller dividend.");
        var rem = new Qc[da + 1];
        for (int i = 0; i <= da; i++) rem[i] = i < dividend.Length ? dividend[i] : Qc.Zero;
        var quo = new Qc[da - df + 1];
        for (int i = da; i >= df; i--)
        {
            var coef = rem[i] / divisor[df];
            quo[i - df] = coef;
            if (coef.IsZero) continue;
            for (int j = 0; j <= df; j++)
                rem[i - df + j] -= coef * divisor[j];
        }
        if (QcDeg(rem) >= 0)
            throw new InvalidOperationException("division expected to be exact left a remainder.");
        return quo;
    }

    private static Qc[] QcGcd(Qc[] a, Qc[] b)
    {
        var x = a;
        var y = b;
        while (QcDeg(y) >= 0)
        {
            var rem = QcRem(x, y);
            x = y;
            y = rem;
        }
        return QcMonic(x);
    }

    /// <summary>The squarefree layers e_i (roots of multiplicity ≥ i, each once) of f, by the gcd chain
    /// f_0 = f, f_i = gcd(f_{i−1}, f_{i−1}′), e_i = f_{i−1}/f_i. Also returns Σ_{i≥1} deg f_i, which for
    /// f = ∏ u_k^k equals Σ_k C(k,2)·deg u_k: the number of same-root pairs (the m_D count).</summary>
    private static (List<Qc[]> Layers, int ChainDegreeSum) SquarefreeLayers(Qc[] f)
    {
        var layers = new List<Qc[]>();
        int chainSum = 0;
        var cur = QcMonic(f);
        while (QcDeg(cur) > 0)
        {
            var g = QcGcd(cur, QcDerivative(cur));
            chainSum += QcDeg(g);
            layers.Add(QcDivExact(cur, g));
            cur = g;
        }
        return (layers, chainSum);
    }
}
