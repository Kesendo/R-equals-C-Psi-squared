using System;
using System.Collections.Generic;
using System.Numerics;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>A Liouvillian-free, one-sided SOFT-certifier (PROOF_F103 §7.12). The true soft criterion is
/// the bipartiteness of the BASIS-STATE hopping graph (a 2^N, letter-dependent object); this class tries
/// scalable SUFFICIENT conditions and certifies soft if any holds. It never claims hard: NotCertified
/// means no scalable strategy applies (the chain-scope hard proxy stays in
/// <see cref="PalindromeMaskClassifier"/>). A certificate is correct for any N and any topology.
///
/// <para>The certifier carries THREE soft mechanisms. (1) The DIAGONAL chiral K: the structured
/// 2-colourings of the basis-state graph, certifying the diagonal −N-mode soft cases. These are linear
/// (the chiral K, <see cref="CertifyByLinearSiteColoring"/>), pure-pairing (⌊n/2⌋ mod 2,
/// <see cref="CertifyByExcitationPairing"/>), excitation-parity (n mod 2,
/// <see cref="CertifyByExcitationParity"/>), and the bit_b-MIXED site-swap reflection
/// (<see cref="CertifyBySiteSwapSymmetry"/>). (2) The NON-DIAGONAL 2-body hidden-Q routing (Stufe A,
/// <see cref="CertifyByRouting"/>): a per-site product Q from the P1/P4 families that palindromizes a sum
/// of 2-body bilinears sharing a uniform Q-family, reaching soft cases whose basis-state graph is
/// NON-bipartite (so no colouring exists at any degree). (3) The DERIVED k-body per-term routing (Stufe B,
/// <see cref="CertifyByRoutingKBody"/>, detailed below): the per-term k-site anticommutator condition that
/// reaches the routable k-body cases the 2-body table misses. The colourings are SOUND but not complete, in two
/// layers. (a) A scalability gap: some soft Hamiltonians are bipartite only through a non-structured
/// colouring no scalable colouring reaches (XY+YX+XZ+ZX on a triangle is soft, its basis-state graph
/// bipartite at ANF-degree 2, but neither linear nor an excitation grading). (b) A structural ceiling,
/// deeper: some soft Hamiltonians have a NON-bipartite basis-state graph (XX+XZ on the chain is soft at
/// N=3..6 by the spectral authority <see cref="PauliPairTrichotomy"/>, yet <see cref="BipartiteChirality"/>
/// reports its basis-state graph non-bipartite). That class is NOT open and NO LONGER beyond the certifier:
/// the hidden-Q routing reaches it (XX+XZ routes to the uniform family {P4}), so it is now certified, not
/// merely named. So NotCertified does not imply not-soft.</para>
///
/// <para>The three colouring strategies all gate on bit_b-homogeneity (one Klein cell); a bit_b-MIXED set
/// they simply decline, even when it is soft. That gate is no longer blunt: the fourth strategy,
/// <see cref="CertifyBySiteSwapSymmetry"/>, recovers the soft bit_b-MIXED cases whose geometry is a spatial
/// bond reflection. Its certificate is that the term-set is a sum of CONTIGUOUS (adjacent) 2-body
/// bilinears, mask-bipartite, and reversal-symmetric (the multiset of terms is invariant under reversing
/// each label, <see cref="IsReversalSymmetric"/>); then the site-swap that reverses the chain is a soft
/// symmetry the chiral colourings miss (XX+XY+YX is soft this way, bit_b-MIXED yet certified). This branch
/// is EMPIRICALLY VERIFIED, not derived: zero false-positives over all ADJACENT 2-body bilinear sums
/// (k = 2..9 terms, N = 3, 4, 5). It is hard-gated to adjacent 2-body bilinears: a 3-body set (XXX+XXY+YXX)
/// or a NON-adjacent 2-body label (XIIX) can be reversal-symmetric, bit_b-MIXED, mask-bipartite yet
/// spectrally HARD, so both are rejected.</para>
///
/// <para>(3) The DERIVED k-body per-term routing (Stufe B, <see cref="CertifyByRoutingKBody"/>): a periodic
/// per-site product Q palindromizes a k-body (k ≥ 3) chain Hamiltonian iff each per-site map class-swaps
/// {I,Z} ↔ {X,Y} (automatic, the dissipator leg) AND the per-term anticommutator {Q_k, [T,·]_k} = 0 at
/// every window-parity (the Hamiltonian leg). This is checked on 4^k (the term's span), is Liouvillian-free,
/// constructive-sound (it exhibits Q), and N-independent by additivity (<see cref="KBodyPalindromeRouting"/>).
/// It reaches the k-body routed-soft cases the 2-body family table misses (XIX+XXY+YXX routes via the P4
/// pattern, IYI+XZY+YZX via a period-2 alternating Q), span-bounded by <see cref="KBodyPalindromeRouting.MaxBody"/>.</para>
///
/// <para>§7.12 ceiling (CLOSED at zero): with both routing mechanisms added, the non-bipartite-soft
/// 2-body class (XX+XZ, Stufe A) and the routable k-body cases (Stufe B) are certified, no longer the
/// ceiling. And the 2 Z-middle cases XZX+XZY+YZX, YZY+XZY+YZX, formerly read as "palindromized only by a
/// non-local Π", route after all: the period-4 GOLDEN per-site router palindromizes each under the
/// WINDOW-SUMMED anticommutator condition (Stufe B′, <see cref="CertifyByRoutingWindowSummed"/>; the
/// cancellation is cross-template inside one window, which the per-term Stufe B cannot see), so they are
/// now certified, see docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md + F116. The formerly-counted cases were
/// already LOCAL: XIX+XIY+YIX, YIY+XIY+YIX route via a continuous-uniform per-site Q (continuous-sum, the
/// 6 to 4 step), and IXI+IIY+YII, IYI+IIX+XII route via a site-varying product of single-site crossover
/// maps, certified by <see cref="CertifyBySingleSiteField"/> (the 4 to 2 step). The arc completes
/// 6 → 4 → 2 → 0: within the k=3 sliding-window soft family no case needs a non-local mirror, see
/// experiments/CEILING_FOUR_NONLOCAL_CASES.md. What remains is coverage/scalability (soft cases whose
/// routers sit outside the scalable strategies), not locality; NotCertified still does not imply
/// not-soft.</para>
///
/// <para>Two-sided front door: <see cref="Certify"/> stays one-sided soft, while the new <see cref="Decide"/>
/// adds the N-free HARD verdict for the diagonal cell (F115, <see cref="WindowedObstructionScan.IsHardPair"/>),
/// the symmetric twin of the soft strategies above. It is gated to two-term Klein-(0,1) Mixed pairs and defers
/// everything else to the spectral authority <see cref="PauliPairTrichotomy"/>. Anchor: F115 /
/// <c>WindowedHardnessClaim</c>.</para></summary>
public static class PalindromeSoftCertifier
{
    /// <summary>Which scalable soft strategy certified the Hamiltonian (None = not certified).</summary>
    public enum SoftStrategy { None, LinearSiteColoring, ExcitationPairing, ExcitationParity, SiteSwapSymmetry, Routing, RoutingKBody, SingleSiteField, RoutingWindowSummed }

    /// <summary>Result of <see cref="Certify"/>: whether soft is certified, and by which strategy.</summary>
    public readonly record struct SoftCertificate(bool Certified, SoftStrategy Strategy);

    /// <summary>The two-sided verdict of <see cref="Decide"/>: Soft (a scalable soft pattern applies),
    /// Hard (the F115 diagonal-cell valuation proves it), or Undetermined (defer to the spectral authority).</summary>
    public enum Decision { Soft, Hard, Undetermined }

    /// <summary>Which scalable HARD strategy decided (None = not hard-certified). Symmetric to
    /// <see cref="SoftStrategy"/>.</summary>
    public enum HardStrategy { None, DiagonalCellValuation }

    /// <summary>Result of <see cref="Decide"/>: the verdict, the deciding strategy (soft or hard), and a
    /// human-readable reason (the soft strategy name, or the exhibited (1+x)-valuation obstruction).</summary>
    public readonly record struct PalindromeDecision(
        Decision Verdict, SoftStrategy SoftStrategy, HardStrategy HardStrategy, string Reason);

    /// <summary>True iff the summed Hamiltonian is a pure pairing (every basis-edge Δn = ±2), detected
    /// by a σ± decomposition: the mixed (hopping) pieces must cancel. N-independent.</summary>
    public static bool IsPurePairing(IReadOnlyList<PauliTerm> terms)
    {
        // The σ± coefficients are the input coefficients scaled by ±1 and ±i, so a true zero is exact;
        // this tolerance only absorbs float round-off in the ±i accumulation.
        const double CoefficientTolerance = 1e-12;
        // Accumulate σ± coefficients keyed by (X/Y mask, Z mask, sign pattern ε). ε bit set = σ_- there.
        var coeffs = new Dictionary<(ulong Xy, ulong Z, ulong Eps), Complex>();
        foreach (var t in terms)
        {
            ulong xyMask = 0, zMask = 0;
            var xyPositions = new List<int>();
            for (int i = 0; i < t.Letters.Count; i++)
            {
                var letter = t.Letters[i];
                if (letter == PauliLetter.X || letter == PauliLetter.Y) { xyMask |= 1UL << i; xyPositions.Add(i); }
                else if (letter == PauliLetter.Z) zMask |= 1UL << i;
            }
            if (xyPositions.Count != 2) return false;   // the Δn=±2 colouring needs exactly 2 X/Y flips per term (this also rejects pure-diagonal terms)
            for (ulong bits = 0; bits < 4; bits++)       // the 4 sign patterns over the 2 X/Y positions pinned above
            {
                ulong eps = 0;
                Complex coeff = t.Coefficient;
                for (int p = 0; p < 2; p++)
                {
                    int pos = xyPositions[p];
                    bool minus = ((bits >> p) & 1UL) != 0;        // this position takes σ_-
                    if (minus) eps |= 1UL << pos;
                    // X: coeff 1 for both signs. Y = -i σ_+ + i σ_-: -i for σ_+, +i for σ_-.
                    if (t.Letters[pos] == PauliLetter.Y)
                        coeff *= minus ? Complex.ImaginaryOne : -Complex.ImaginaryOne;
                }
                var key = (xyMask, zMask, eps);
                coeffs[key] = coeffs.GetValueOrDefault(key) + coeff;
            }
        }
        bool anyPure = false;
        foreach (var kv in coeffs)
        {
            bool allPlus = kv.Key.Eps == 0;
            bool allMinus = kv.Key.Eps == kv.Key.Xy;
            if (allPlus || allMinus)
            {
                if (kv.Value.Magnitude > CoefficientTolerance) anyPure = true;
            }
            else if (kv.Value.Magnitude > CoefficientTolerance)
            {
                return false;                                     // a surviving mixed (hopping) piece
            }
        }
        return anyPure;                                           // pure pairing iff only ± pieces survive
    }

    /// <summary>The excitation-number strategy: certify soft iff the Hamiltonian is a pure pairing
    /// (then ⌊n/2⌋ mod 2 two-colours the basis-state graph, soft on any topology).</summary>
    public static bool CertifyByExcitationPairing(IReadOnlyList<PauliTerm> terms) => IsPurePairing(terms);

    /// <summary>True iff every term flips an ODD number of sites (odd k_xy = #X/Y per term), i.e. the
    /// Hamiltonian sits in the bit_a = 1 Klein-cell row. Then every basis-edge has odd Δn. N-independent.</summary>
    public static bool IsAllOddFlip(IReadOnlyList<PauliTerm> terms)
    {
        if (terms.Count == 0) return false;
        foreach (var t in terms)
        {
            int kxy = 0;
            foreach (var letter in t.Letters)
                if (letter == PauliLetter.X || letter == PauliLetter.Y) kxy++;
            if (kxy % 2 == 0) return false;   // an even (incl. zero) X/Y count gives an even-Δn edge
        }
        return true;
    }

    /// <summary>The excitation-parity strategy: certify soft iff every term has odd k_xy AND the terms
    /// share one Klein cell (bit_b-homogeneous). Then every basis-edge has odd Δn, so n mod 2 two-colours
    /// the basis-state graph; but a bipartite basis-state graph only certifies soft WITHIN a single Klein
    /// cell (the §7.6 diagonal-cell scope). A bit_b-MIXED all-odd-flip set can be hard despite the n mod 2
    /// colouring (e.g. XZ + ZXZ, bit_b = {1, 0}, is spectrally hard at N = 3..4), so the bit_b gate is a
    /// soundness requirement, not a refinement. The odd sibling of the pure-pairing strategy. N-independent.</summary>
    public static bool CertifyByExcitationParity(IReadOnlyList<PauliTerm> terms) =>
        IsAllOddFlip(terms) && PalindromeMaskClassifier.IsBitBHomogeneous(terms);

    /// <summary>The linear site-colouring strategy: certify soft iff the chain flip-mask set is bipartite
    /// (the chiral K). Reuses <see cref="PalindromeMaskClassifier"/>.</summary>
    public static bool CertifyByLinearSiteColoring(IReadOnlyList<PauliTerm> terms, int n)
    {
        // The mask-bipartite test is a valid soft certificate only WITHIN a single Klein cell; a mixed
        // cell (terms of different bit_b = #(Y+Z) parity) can be hard while still mask-bipartite. Gate on
        // bit_b-homogeneity, matching PalindromeMaskClassifier.Classify, to avoid a false positive.
        if (!PalindromeMaskClassifier.IsBitBHomogeneous(terms)) return false;
        // A pure-diagonal term (no X/Y) lifts the diagonal, which the diagonal chiral K cannot negate,
        // so the bipartite flip-graph would not actually certify soft. Reject to avoid a false positive.
        foreach (var t in terms)
        {
            bool hasFlip = false;
            foreach (var letter in t.Letters)
                if (letter == PauliLetter.X || letter == PauliLetter.Y) { hasFlip = true; break; }
            if (!hasFlip) return false;
        }
        var masks = PalindromeMaskClassifier.FlipMasks(terms, n);
        return masks.Count > 0 && PalindromeMaskClassifier.MaskSetIsBipartite(masks);
    }

    /// <summary>True iff the term-set is invariant under reversing each term's label (the site-swap that
    /// reverses the chain): the multiset of (reversed-label, coefficient) equals the multiset of
    /// (label, coefficient). A term's reversal reverses its <see cref="PauliLetter"/> sequence
    /// ("XY" → "YX", "XYZ" → "ZYX"). The spatial-reflection witness used by the site-swap strategy.
    /// N-independent.</summary>
    public static bool IsReversalSymmetric(IReadOnlyList<PauliTerm> terms)
    {
        // Multiset over (label, coefficient): a reversal symmetry is the two multisets coinciding.
        var original = new Dictionary<(string Label, Complex Coefficient), int>();
        var reversed = new Dictionary<(string Label, Complex Coefficient), int>();
        foreach (var t in terms)
        {
            int count = t.Letters.Count;
            var rev = new PauliLetter[count];
            for (int i = 0; i < count; i++) rev[i] = t.Letters[count - 1 - i];   // reverse the letter sequence
            var oKey = (t.Label, t.Coefficient);
            var rKey = (PauliLabel.Format(rev), t.Coefficient);
            original[oKey] = original.GetValueOrDefault(oKey) + 1;
            reversed[rKey] = reversed.GetValueOrDefault(rKey) + 1;
        }
        if (original.Count != reversed.Count) return false;
        foreach (var kv in original)
            if (reversed.GetValueOrDefault(kv.Key) != kv.Value) return false;
        return true;
    }

    /// <summary>The site-swap-symmetry strategy: certify soft for a bit_b-MIXED set whose geometry is a
    /// spatial bond reflection. Certifies iff ALL hold: every term is a CONTIGUOUS bilinear (exactly two
    /// non-identity letters on adjacent sites), the set is bit_b-MIXED (the gate the three colourings
    /// decline), no term is pure-diagonal (some X/Y), the chain flip-mask set is bipartite, and the set is
    /// reversal-symmetric (<see cref="IsReversalSymmetric"/>). The reversal symmetry is load-bearing:
    /// without it 2-body mask-bipartite mixed sets can be hard. EMPIRICALLY VERIFIED (zero false-positives
    /// over all ADJACENT 2-body bilinear sums, k = 2..9 terms, N = 3, 4, 5), not derived. The contiguity
    /// and 2-body gates are soundness requirements: at 3-body the rule false-positives (XXX+XXY+YXX is
    /// reversal-symmetric, bit_b-MIXED, mask-bipartite, yet spectrally HARD at N = 4, 5), and a NON-adjacent
    /// 2-body label false-positives too ({XIIX, XY, YX} is hard at N = 5), so both are rejected here.</summary>
    public static bool CertifyBySiteSwapSymmetry(IReadOnlyList<PauliTerm> terms, int n)
    {
        // 2-body gate: every term is a CONTIGUOUS bilinear, exactly two non-identity letters on ADJACENT
        // sites. This is the verified scope (adjacent 2-body bilinear sums). It rejects the 3-body killer
        // XXX+XXY+YXX, AND non-adjacent 2-body labels like XIIX that a bare non-identity count would admit:
        // those were never in the sweep and can be spectrally hard ({XIIX, XY, YX} is hard at N=5).
        foreach (var t in terms)
        {
            int first = -1, last = -1, nonId = 0;
            for (int i = 0; i < t.Letters.Count; i++)
                if (t.Letters[i] != PauliLetter.I)
                {
                    if (first < 0) first = i;
                    last = i;
                    nonId++;
                }
            if (nonId != 2 || last - first != 1) return false;   // exactly two non-identity letters, adjacent
        }
        // This strategy exists for the cells the three colourings decline: require bit_b-MIXED.
        if (PalindromeMaskClassifier.IsBitBHomogeneous(terms)) return false;
        // A pure-diagonal term (no X/Y) lifts the diagonal the reflection cannot negate; reject it,
        // matching the guard in CertifyByLinearSiteColoring.
        foreach (var t in terms)
        {
            bool hasFlip = false;
            foreach (var letter in t.Letters)
                if (letter == PauliLetter.X || letter == PauliLetter.Y) { hasFlip = true; break; }
            if (!hasFlip) return false;
        }
        var masks = PalindromeMaskClassifier.FlipMasks(terms, n);
        if (masks.Count == 0 || !PalindromeMaskClassifier.MaskSetIsBipartite(masks)) return false;
        return IsReversalSymmetric(terms);
    }

    /// <summary>True iff the label is a Mother bilinear {XX, YY, ZZ} (Klein cell (0,0)).</summary>
    private static bool IsMotherLabel(string l) => l is "XX" or "YY" or "ZZ";

    /// <summary>The hidden-Q routing strategy: the NON-DIAGONAL soft mechanism the colourings lack. A sum
    /// of 2-body bilinears is certified soft when one uniform per-site product Q (a P1/P4 family member,
    /// see <see cref="TwoTermPalindromeRouting"/>) palindromizes EVERY term, hence the whole sum.
    ///
    /// <para>The rule reads LABELS only (coefficient-independent: one uniform Q palindromizes any real
    /// linear combination, Q(Σ cᵢ tᵢ)Q⁻¹ = Σ cᵢ(−tᵢ) = −H). For each term it looks up the per-term family
    /// mask via <see cref="TwoTermPalindromeRouting.TryGetUniformFamilyMask"/>; if ANY term is not a
    /// recognized 2-body bilinear (a k-body or padded label), the set is out of scope and declined (Stufe B,
    /// the k-body routed-soft frontier is deferred). Otherwise it certifies in two cases:</para>
    /// <list type="number">
    ///   <item>The uniform-routing certificate: the family masks share a member (their bitwise-AND is
    ///     non-zero) AND the set is not all-Mother ({XX, YY, ZZ}, which is truly, not soft, the canonical
    ///     Π already pairs the spectrum). That shared Q palindromizes every term ⟹ the sum.</item>
    ///   <item>The two-term alternating/continuous escapes (XY+ZZ, XZ+YZ, ...): at exactly two terms, the
    ///     fate is read directly from <see cref="TwoTermPalindromeRouting.Classify"/> (these are
    ///     TWO-TERM-SPECIFIC, not generalized to 3+ terms).</item>
    /// </list>
    /// <para>A single recognized non-Mother bilinear is in scope and certified via its own uniform family
    /// (e.g. XZ alone routes to {P4}), while a single Mother bilinear (XX, YY, or ZZ) is excluded by the
    /// all-Mother gate: it is truly, not soft, the canonical Π already pairs its spectrum.</para>
    ///
    /// <para>SOUND by the additivity argument (the Q is exhibited as a uniform per-site product, not an
    /// N-dependent graph property, so the certificate is N-stable) and EMPIRICALLY VERIFIED: the soundness
    /// sweep found zero Hard false-positives over all multi-term 2-body bilinear sums (k = 2..5, N = 4, 5),
    /// N-stable through N = 6.</para></summary>
    public static bool CertifyByRouting(IReadOnlyList<PauliTerm> terms)
    {
        if (terms.Count == 0) return false;

        int intersection = 0;
        bool allMother = true;
        bool first = true;
        foreach (var t in terms)
        {
            if (!TwoTermPalindromeRouting.TryGetUniformFamilyMask(t.Label, out int m))
                return false;   // a non-recognized / k-body label: out of the 2-body scope (Stufe B)
            intersection = first ? m : intersection & m;
            first = false;
            if (!IsMotherLabel(t.Label)) allMother = false;
        }

        // The uniform-routing certificate: one shared product Q palindromizes every term (the sum), unless
        // the set is all-Mother (truly, the canonical Π already pairs it).
        if (intersection != 0 && !allMother) return true;

        // The two-term alternating/continuous escapes (read directly, two-term-specific).
        if (terms.Count == 2)
        {
            var r = TwoTermPalindromeRouting.Classify(terms[0].Label, terms[1].Label);
            if (r.Fate == TrichotomyClass.Soft) return true;
        }

        return false;
    }

    /// <summary>The DERIVED k-body hidden-Q routing strategy (Stufe B): the per-term, k-site-local,
    /// Liouvillian-free palindrome condition that reaches the k-body (k ≥ 3) routed-soft cases the 2-body
    /// family table (<see cref="CertifyByRouting"/>) cannot. The palindrome condition Q L Q⁻¹ = −L − 2σ
    /// decomposes EXACTLY into (a) the automatic per-site dissipator swap {I,Z} ↔ {X,Y} (every candidate
    /// map is class-swapping, n_XY → N − n_XY gives the −2σ shift) AND (b) the per-term anticommutator
    /// {Q_k, [T,·]_k} = 0 for every placed term at every window-parity. Condition (b) is checked on 4^k
    /// (the term's span, NOT the 2^N Liouvillian) and additivity over windows gives the palindrome at EVERY
    /// N, so the certificate is CONSTRUCTIVE (it exhibits the periodic per-site product Q), sound by
    /// derivation, and N-independent. See <see cref="KBodyPalindromeRouting"/>.
    ///
    /// <para>Span-bounded: a term of span k yields a 4^k × 4^k check, so the strategy declines a set with
    /// any term outside [2, <see cref="KBodyPalindromeRouting.MaxBody"/>]. The 2 Z-middle ceiling cases
    /// (XZX+XZY+YZX, YZY+XZY+YZX) fail the PER-TERM condition (each template's anticommutator alone is
    /// nonzero), so <see cref="KBodyPalindromeRouting.Routes"/> correctly returns false for them; they ARE
    /// per-site routable, by the golden period-4 router under the window-summed condition, and are
    /// certified by <see cref="CertifyByRoutingWindowSummed"/> (Stufe B′, F116). The two I-heavy cases
    /// IXI+IIY+YII, IYI+IIX+XII also return false here (the per-term router does not see their
    /// single-site-field router), but they ARE local, certified by
    /// <see cref="CertifyBySingleSiteField"/>.</para></summary>
    public static bool CertifyByRoutingKBody(IReadOnlyList<PauliTerm> terms, int n)
    {
        if (terms.Count == 0) return false;

        // Span gate: every term's span k = label length must be in [2, MaxBody]. A span-1 term is not a
        // routed bilinear; a span > MaxBody term makes the 4^k anticommutator check large (out of scope).
        foreach (var t in terms)
        {
            int k = t.Letters.Count;
            if (k < 2 || k > KBodyPalindromeRouting.MaxBody) return false;
        }

        return KBodyPalindromeRouting.Routes(terms, n);
    }

    /// <summary>The single-site-field strategy: certify soft iff every term is weight-1 with letter X or Y
    /// (a transverse single-site field). Then the chain Hamiltonian is a sum of single-site transverse fields
    /// H = Σ_i (a_i X_i + b_i Y_i), so L = Σ_i L_i over COMMUTING single-site Liouvillians, and the per-site
    /// product Q = ⊗_i M_i palindromizes the whole chain (each M_i the per-site crossover map Ad_{R_z(θ_i)},
    /// θ_i = atan2(b_i, a_i)). Constructive, N-independent, sound by derivation.
    ///
    /// <para>SOUND because the detection is gated to TRANSVERSE: a single-site X/Y field is soft (its
    /// Liouvillian spectrum {0, −2γ, −γ ± 2i} is palindromic about −γ), and a sum of commuting soft single-site
    /// Liouvillians is soft. Z is EXCLUDED: a single-site Z (longitudinal) field has spectrum {0, 0, −2γ ± 2i},
    /// whose 0 eigenvalue has no partner −2γ about −γ, so it is HARD; certifying it would break the one-sided
    /// soundness. This reaches the two I-heavy cases (IXI+IIY+YII, IYI+IIX+XII) the other strategies decline
    /// (bit_b-MIXED weight-1), correcting the §7.12 ceiling from 4 to 2.</para></summary>
    public static bool CertifyBySingleSiteField(IReadOnlyList<PauliTerm> terms)
    {
        if (terms.Count == 0) return false;
        foreach (var t in terms)
        {
            int nonId = 0;
            PauliLetter only = PauliLetter.I;
            foreach (var letter in t.Letters)
                if (letter != PauliLetter.I) { nonId++; only = letter; }
            if (nonId != 1) return false;                                      // not weight-1
            if (only != PauliLetter.X && only != PauliLetter.Y) return false;  // Z is longitudinal => hard, exclude
        }
        return true;
    }

    /// <summary>The DERIVED window-summed k-body routing strategy (Stufe B′, F116): the golden period-4
    /// per-site router for the term-sets whose palindromizer cancels CROSS-TEMPLATE inside one window,
    /// which the per-term Stufe B (<see cref="CertifyByRoutingKBody"/>) cannot see. The palindrome
    /// condition Q L Q⁻¹ = −L − 2σ decomposes exactly as in Stufe B into (a) the automatic per-site
    /// dissipator class-swap {I,Z} ↔ {X,Y} AND (b′) the TEMPLATE-SUMMED anticommutator
    /// {Q_k, Σ_T [T,·]_k} = 0 at every window offset 0..3, the sharp window-level condition (the per-term
    /// (b) is sufficient, not necessary). Checked on 4^k; window additivity gives the palindrome at EVERY
    /// N ≥ k, so the certificate is CONSTRUCTIVE (the candidate is exhibited by
    /// <see cref="KBodyPalindromeRouting.RoutesWindowSummed"/>), sound by derivation, and N-independent.
    /// This certifies the two formerly "non-local" Z-middle ceiling cases, XZX+XZY+YZX (the golden
    /// [a, a, b, b] pattern) and YZY+XZY+YZX (its X↔Y conjugate), closing the §7.12 ceiling at zero. See
    /// docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md + F116.</summary>
    public static bool CertifyByRoutingWindowSummed(IReadOnlyList<PauliTerm> terms, int n)
    {
        if (terms.Count == 0) return false;
        return KBodyPalindromeRouting.RoutesWindowSummed(terms, n) is not null;
    }

    /// <summary>The HARD-side strategy (PROOF_F103 §7.7 / F115): a two-term Klein-(0,1) Mixed pair is
    /// hard iff its two X/Y flip-masks have different (1+x)-adic valuations
    /// (<see cref="WindowedObstructionScan.IsHardPair"/>). N-free, O(k); the symmetric twin of the soft
    /// strategies. Gated to the proven scope (exactly two diagonal-cell Mixed, bit_b-homogeneous,
    /// y_par-homogeneous templates under Z-dephasing); anything else returns false so the caller defers
    /// to the spectral authority. Soundness gated by PalindromeHardSweepTests.</summary>
    public static bool CertifyHardByDiagonalCellValuation(IReadOnlyList<PauliTerm> terms)
    {
        if (terms.Count != 2) return false;
        if (!PalindromeMaskClassifier.IsBitBHomogeneous(terms)) return false;
        if (!TryDiagonalCellMixedMask(terms[0], out ulong m0, out int yPar0)) return false;
        if (!TryDiagonalCellMixedMask(terms[1], out ulong m1, out int yPar1)) return false;
        if (yPar0 != yPar1) return false;                       // y_par-homogeneous
        return WindowedObstructionScan.IsHardPair(m0, m1);
    }

    /// <summary>If <paramref name="term"/> is a Klein-(0,1) diagonal-cell Mixed string (X/Y count even
    /// and >= 2, #(Y/Z) odd), output its X/Y flip-mask and #Y parity and return true; else false. Same
    /// gate as <see cref="WindowedObstructionScan.CellTerms"/>.</summary>
    private static bool TryDiagonalCellMixedMask(PauliTerm term, out ulong xyMask, out int yParity)
    {
        xyMask = 0;
        int na = 0, nb = 0, ny = 0;
        for (int i = 0; i < term.Letters.Count; i++)
        {
            var l = term.Letters[i];
            if (l == PauliLetter.X || l == PauliLetter.Y) { na++; xyMask |= 1UL << i; }
            if (l == PauliLetter.Y || l == PauliLetter.Z) nb++;
            if (l == PauliLetter.Y) ny++;
        }
        yParity = ny & 1;
        return (na & 1) == 0 && na >= 2 && (nb & 1) == 1;
    }

    /// <summary>Try the stronger, topology-independent excitation strategies first (pairing, then
    /// parity; a term-set is at most one of them), then the chain-only linear one, then the bit_b-MIXED
    /// site-swap-symmetry one, then the 2-body hidden-Q routing (Stufe A), then the derived k-body per-term
    /// routing residual (Stufe B), then the single-site-field products, and LAST the window-summed golden
    /// routing (Stufe B′); return the certificate. A certified set is at most one strategy (the earlier
    /// strategies take precedence on overlap, e.g. XY+YX stays ExcitationPairing, XZ+ZX stays
    /// ExcitationParity, and a pure-2-body routed set stays Routing; k ≥ 3 and mixed-span sets fall to
    /// RoutingKBody; weight-1 transverse single-site sums fall to SingleSiteField; the two Z-middle golden
    /// cases fall to RoutingWindowSummed, and its last place keeps every pre-existing certificate on its
    /// original strategy label).</summary>
    public static SoftCertificate Certify(IReadOnlyList<PauliTerm> terms, int n)
    {
        if (CertifyByExcitationPairing(terms)) return new SoftCertificate(true, SoftStrategy.ExcitationPairing);
        if (CertifyByExcitationParity(terms)) return new SoftCertificate(true, SoftStrategy.ExcitationParity);
        if (CertifyByLinearSiteColoring(terms, n)) return new SoftCertificate(true, SoftStrategy.LinearSiteColoring);
        if (CertifyBySiteSwapSymmetry(terms, n)) return new SoftCertificate(true, SoftStrategy.SiteSwapSymmetry);
        if (CertifyByRouting(terms)) return new SoftCertificate(true, SoftStrategy.Routing);
        if (CertifyByRoutingKBody(terms, n)) return new SoftCertificate(true, SoftStrategy.RoutingKBody);
        if (CertifyBySingleSiteField(terms)) return new SoftCertificate(true, SoftStrategy.SingleSiteField);
        if (CertifyByRoutingWindowSummed(terms, n)) return new SoftCertificate(true, SoftStrategy.RoutingWindowSummed);
        return new SoftCertificate(false, SoftStrategy.None);
    }

    /// <summary>Two-sided decider: certify Soft (the existing <see cref="Certify"/> strategies), else
    /// certify Hard (the F115 diagonal-cell valuation, <see cref="CertifyHardByDiagonalCellValuation"/>),
    /// else Undetermined (the caller defers to <see cref="PauliPairTrichotomy"/>). <see cref="Certify"/>
    /// stays soft-only; this is the additive, more-powerful front door (N-free on both sides where it
    /// decides).</summary>
    public static PalindromeDecision Decide(IReadOnlyList<PauliTerm> terms, int n)
    {
        var soft = Certify(terms, n);
        if (soft.Certified)
            return new PalindromeDecision(Decision.Soft, soft.Strategy, HardStrategy.None, $"soft: {soft.Strategy}");

        if (CertifyHardByDiagonalCellValuation(terms))
        {
            TryDiagonalCellMixedMask(terms[0], out ulong m0, out _);
            TryDiagonalCellMixedMask(terms[1], out ulong m1, out _);
            int v0 = WindowedObstructionScan.ValuationAtOnePlusX(m0);
            int v1 = WindowedObstructionScan.ValuationAtOnePlusX(m1);
            return new PalindromeDecision(Decision.Hard, SoftStrategy.None, HardStrategy.DiagonalCellValuation,
                $"hard: diagonal-cell (1+x)-valuation {v0} != {v1} (F115 odd-cycle obstruction)");
        }

        return new PalindromeDecision(Decision.Undetermined, SoftStrategy.None, HardStrategy.None,
            "undetermined: no scalable soft pattern, out of F115 hard scope; defer to PauliPairTrichotomy");
    }
}
