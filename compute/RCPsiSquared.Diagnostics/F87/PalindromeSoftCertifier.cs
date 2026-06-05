using System;
using System.Collections.Generic;
using System.Numerics;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>A Liouvillian-free, one-sided SOFT-certifier (PROOF_F103 §7.12). The true soft criterion is
/// the bipartiteness of the BASIS-STATE hopping graph (a 2^N, letter-dependent object); this class tries
/// scalable SUFFICIENT conditions ("soft colourings") and certifies soft if any holds. It never claims
/// hard: NotCertified means no scalable strategy applies (the chain-scope hard proxy stays in
/// <see cref="PalindromeMaskClassifier"/>). A certificate is correct for any N and any topology.</summary>
public static class PalindromeSoftCertifier
{
    /// <summary>Which scalable soft-colouring certified the Hamiltonian (None = not certified).</summary>
    public enum SoftStrategy { None, LinearSiteColoring, ExcitationPairing }

    /// <summary>Result of <see cref="Certify"/>: whether soft is certified, and by which strategy.</summary>
    public readonly record struct SoftCertificate(bool Certified, SoftStrategy Strategy);

    /// <summary>True iff the summed Hamiltonian is a pure pairing (every basis-edge Δn = ±2), detected
    /// by a σ± decomposition: the mixed (hopping) pieces must cancel. N-independent.</summary>
    public static bool IsPurePairing(IReadOnlyList<PauliTerm> terms)
    {
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
            if (xyPositions.Count == 0) return false;   // a pure-diagonal term: not a pairing
            int k = xyPositions.Count;
            for (ulong bits = 0; bits < (1UL << k); bits++)
            {
                ulong eps = 0;
                Complex coeff = t.Coefficient;
                for (int p = 0; p < k; p++)
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
                if (kv.Value.Magnitude > 1e-12) anyPure = true;
            }
            else if (kv.Value.Magnitude > 1e-12)
            {
                return false;                                     // a surviving mixed (hopping) piece
            }
        }
        return anyPure;                                           // pure pairing iff only ± pieces survive
    }

    /// <summary>The excitation-number strategy: certify soft iff the Hamiltonian is a pure pairing
    /// (then ⌊n/2⌋ mod 2 two-colours the basis-state graph, soft on any topology).</summary>
    public static bool CertifyByExcitationPairing(IReadOnlyList<PauliTerm> terms) => IsPurePairing(terms);

    /// <summary>The linear site-colouring strategy: certify soft iff the chain flip-mask set is bipartite
    /// (the chiral K). Reuses <see cref="PalindromeMaskClassifier"/>.</summary>
    public static bool CertifyByLinearSiteColoring(IReadOnlyList<PauliTerm> terms, int n)
    {
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

    /// <summary>Try the stronger excitation strategy first, then the linear one; return the certificate.</summary>
    public static SoftCertificate Certify(IReadOnlyList<PauliTerm> terms, int n) => throw new NotImplementedException();
}
