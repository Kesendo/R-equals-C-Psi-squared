using System;
using System.Collections.Generic;
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
/// <para>DIAGNOSTICS kept in the report (first prime): the discriminant's squarefree layers, the N=5
/// image of the N=4 disc structure disc = q²⁴·(diabolic)²·P_20 (experiments/F89_BRANCH_LOCUS_PALINDROME.md):
/// the multiplicity-1 layer carries the √-branch (defective) loci, the multiplicity-2 layer the squared
/// diabolic factor, and R shares nothing with either (the layer gcds are reported).</para>
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
        bool SharedIsQPowerAtEveryPrime, bool Complete);

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
        for (int tried = 0; tried < maxPrimes && !(primeProduct > proofBound && sampled > lcDivisorBound); candidate += 4)
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
        int[] layerDegrees = Array.Empty<int>(), layerGcdDegrees = Array.Empty<int>();
        foreach (var s in perPrime)
        {
            if (s.DegR != trueDegR) { skipped++; continue; }        // divided lc_q(R): demoted
            if (firstPrime == 0)
            {
                firstPrime = s.P;
                degRpFirst = s.DegR; degDpFirst = s.DegD; vR = s.VR; vD = s.VD;
                var layers = FpSquarefreeLayers(s.DStrip, s.P);     // diagnostics: the disc's layer structure
                layerDegrees = new int[layers.Count];
                layerGcdDegrees = new int[layers.Count];
                for (int k = 0; k < layers.Count; k++)
                {
                    layerDegrees[k] = DegP(layers[k]);
                    layerGcdDegrees[k] = DegP(GcdModP(s.RStrip, layers[k], s.P));
                }
            }
            if (!s.SharedOk) { sharedIsQPower = false; lastPrime = s.P; continue; }
            used++;
            lastPrime = s.P;
        }

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
            sharedIsQPower, complete);
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
