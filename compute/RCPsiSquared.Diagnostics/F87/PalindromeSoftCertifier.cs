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
/// <para>The certifier carries BOTH soft mechanisms. (1) The DIAGONAL chiral K: the structured
/// 2-colourings of the basis-state graph, certifying the diagonal −N-mode soft cases. These are linear
/// (the chiral K, <see cref="CertifyByLinearSiteColoring"/>), pure-pairing (⌊n/2⌋ mod 2,
/// <see cref="CertifyByExcitationPairing"/>), excitation-parity (n mod 2,
/// <see cref="CertifyByExcitationParity"/>), and the bit_b-MIXED site-swap reflection
/// (<see cref="CertifyBySiteSwapSymmetry"/>). (2) The NON-DIAGONAL hidden-Q routing
/// (<see cref="CertifyByRouting"/>): a per-site product Q from the P1/P4 families that palindromizes a sum
/// of 2-body bilinears sharing a uniform Q-family, reaching soft cases whose basis-state graph is
/// NON-bipartite (so no colouring exists at any degree). The colourings are SOUND but not complete, in two
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
/// <para>§7.12 ceiling (the remaining frontier): with the hidden-Q routing added, the non-bipartite-soft
/// 2-body class (XX+XZ) is no longer the ceiling, it is certified. What stays beyond the certifier is the
/// k-body routed-soft frontier (Stufe B): the routing family table is 2-body, so a 3-body routed-soft case
/// like XZX+XZY+YZX (soft by the spectral authority at N=4,5,6, NotCertified) is the current ceiling.</para></summary>
public static class PalindromeSoftCertifier
{
    /// <summary>Which scalable soft strategy certified the Hamiltonian (None = not certified).</summary>
    public enum SoftStrategy { None, LinearSiteColoring, ExcitationPairing, ExcitationParity, SiteSwapSymmetry, Routing }

    /// <summary>Result of <see cref="Certify"/>: whether soft is certified, and by which strategy.</summary>
    public readonly record struct SoftCertificate(bool Certified, SoftStrategy Strategy);

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

    /// <summary>Try the stronger, topology-independent excitation strategies first (pairing, then
    /// parity; a term-set is at most one of them), then the chain-only linear one, then the bit_b-MIXED
    /// site-swap-symmetry one, then the non-diagonal hidden-Q routing residual; return the certificate.
    /// A certified set is at most one strategy (the existing strategies take precedence on overlap, e.g.
    /// XY+YX stays ExcitationPairing, XZ+ZX stays ExcitationParity).</summary>
    public static SoftCertificate Certify(IReadOnlyList<PauliTerm> terms, int n)
    {
        if (CertifyByExcitationPairing(terms)) return new SoftCertificate(true, SoftStrategy.ExcitationPairing);
        if (CertifyByExcitationParity(terms)) return new SoftCertificate(true, SoftStrategy.ExcitationParity);
        if (CertifyByLinearSiteColoring(terms, n)) return new SoftCertificate(true, SoftStrategy.LinearSiteColoring);
        if (CertifyBySiteSwapSymmetry(terms, n)) return new SoftCertificate(true, SoftStrategy.SiteSwapSymmetry);
        if (CertifyByRouting(terms)) return new SoftCertificate(true, SoftStrategy.Routing);
        return new SoftCertificate(false, SoftStrategy.None);
    }
}
