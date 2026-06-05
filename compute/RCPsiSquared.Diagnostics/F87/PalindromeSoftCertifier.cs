using System;
using System.Collections.Generic;
using System.Numerics;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>A Liouvillian-free, one-sided SOFT-certifier (PROOF_F103 §7.12). The true soft criterion is
/// the bipartiteness of the BASIS-STATE hopping graph (a 2^N, letter-dependent object); this class tries
/// scalable SUFFICIENT conditions ("soft colourings") and certifies soft if any holds. It never claims
/// hard: NotCertified means no scalable strategy applies (the chain-scope hard proxy stays in
/// <see cref="PalindromeMaskClassifier"/>). A certificate is correct for any N and any topology.
///
/// <para>The strategies are the structured 2-colourings of the basis-state graph: linear (the chiral K,
/// <see cref="CertifyByLinearSiteColoring"/>), pure-pairing (⌊n/2⌋ mod 2, <see cref="CertifyByExcitationPairing"/>),
/// and excitation-parity (n mod 2, <see cref="CertifyByExcitationParity"/>). They are SOUND but not
/// complete: some soft Hamiltonians have a non-structured basis-graph 2-colouring that no scalable
/// strategy reaches (XY+YX+XZ+ZX placed on a triangle is soft, yet bipartite only via a colouring that is
/// neither linear nor an excitation grading), so NotCertified does not imply not-soft. That residual is
/// the price of staying Liouvillian-free; the full criterion is the 2^N basis-state graph itself.</para></summary>
public static class PalindromeSoftCertifier
{
    /// <summary>Which scalable soft-colouring certified the Hamiltonian (None = not certified).</summary>
    public enum SoftStrategy { None, LinearSiteColoring, ExcitationPairing, ExcitationParity }

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

    /// <summary>The excitation-parity strategy: certify soft iff every term has odd k_xy (so every
    /// basis-edge has odd Δn and the excitation parity n mod 2 two-colours the basis-state graph, soft on
    /// any topology). The odd sibling of the pure-pairing strategy. N-independent.</summary>
    public static bool CertifyByExcitationParity(IReadOnlyList<PauliTerm> terms) => IsAllOddFlip(terms);

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

    /// <summary>Try the stronger, topology-independent excitation strategies first (pairing, then
    /// parity; a term-set is at most one of them), then the chain-only linear one; return the certificate.</summary>
    public static SoftCertificate Certify(IReadOnlyList<PauliTerm> terms, int n)
    {
        if (CertifyByExcitationPairing(terms)) return new SoftCertificate(true, SoftStrategy.ExcitationPairing);
        if (CertifyByExcitationParity(terms)) return new SoftCertificate(true, SoftStrategy.ExcitationParity);
        if (CertifyByLinearSiteColoring(terms, n)) return new SoftCertificate(true, SoftStrategy.LinearSiteColoring);
        return new SoftCertificate(false, SoftStrategy.None);
    }
}
