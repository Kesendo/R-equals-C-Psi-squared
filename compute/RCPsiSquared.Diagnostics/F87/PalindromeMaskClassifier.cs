using System.Collections.Generic;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>A Liouvillian-free palindrome classifier for bit_b-homogeneous hopping Hamiltonians.
///
/// <para>For a single-Klein-cell hopping Hamiltonian the dephased spectrum is palindromic (soft) iff
/// the set of flip-masks the terms produce carries no odd 𝔽₂-relation, equivalently iff a chiral
/// K = diag(±1) with K H K = −H exists (the §7.1 test). That is a pure GF(2) question on N-bit masks,
/// so it classifies H at any N without building the 2^N Liouvillian; only the flip-mask set is needed.</para>
///
/// <para><b>Scope (verified sharp at N = 4 against the trichotomy):</b> the mask test gives the correct
/// soft/hard verdict exactly when the terms are bit_b-homogeneous (a single Klein cell) AND have hopping
/// (some X/Y). Mixed-cell Hamiltonians (terms of different #(Y+Z) parity) and pure-diagonal lifts fall
/// outside, returning <see cref="Verdict.OutOfScope"/>: their soft/hard is not fixed by the flip-mask set
/// alone (XX+YY is truly but XX+YY+XY is hard, and they share the same flip-mask set). Anchor: PROOF_F103
/// §7.1, F87; the cancellation caveat is §7.10.</para></summary>
public static class PalindromeMaskClassifier
{
    /// <summary>Soft (palindrome restorable), Hard (palindrome broken), or OutOfScope (the flip-mask
    /// set does not determine the verdict: mixed Klein cell, or no hopping).</summary>
    public enum Verdict { Soft, Hard, OutOfScope }

    /// <summary>The distinct nonzero flip-masks of the k-body term templates placed across an N-site
    /// chain by the sliding-window builder: the X/Y positions of each term at each window 0 … N−k.</summary>
    public static HashSet<ulong> FlipMasks(IReadOnlyList<PauliTerm> terms, int n)
    {
        var masks = new HashSet<ulong>();
        foreach (var t in terms)
        {
            int k = t.Letters.Count;
            for (int w = 0; w + k <= n; w++)
            {
                ulong m = 0;
                for (int p = 0; p < k; p++)
                    if (t.Letters[p] == PauliLetter.X || t.Letters[p] == PauliLetter.Y)
                        m |= 1UL << (w + p);
                if (m != 0) masks.Add(m);
            }
        }
        return masks;
    }

    /// <summary>True iff all terms share the same bit_b = #(Y+Z) mod 2 parity (a single Klein cell).</summary>
    public static bool IsBitBHomogeneous(IReadOnlyList<PauliTerm> terms)
    {
        if (terms.Count == 0) return true;
        int b = terms[0].Pi2Parity;
        for (int i = 1; i < terms.Count; i++)
            if (terms[i].Pi2Parity != b) return false;
        return true;
    }

    /// <summary>True iff a linear φ over GF(2) with φ(m) = 1 for every mask exists (the chiral-K
    /// existence = no odd 𝔽₂-relation among the masks). Gaussian elimination over GF(2) on the rows
    /// [m | 1] (target bit appended in bit 0); a reduction to the bare target (value 1) is 0 = 1,
    /// i.e. an odd relation, i.e. non-bipartite.</summary>
    public static bool MaskSetIsBipartite(IReadOnlyCollection<ulong> masks)
    {
        var basis = new List<ulong>();
        foreach (ulong mask in masks)
        {
            ulong r = (mask << 1) | 1UL;
            foreach (ulong b in basis)
            {
                ulong x = r ^ b;
                if (x < r) r = x;
            }
            if (r == 1UL) return false;     // 0 = 1: an odd relation
            if (r != 0)
            {
                basis.Add(r);
                basis.Sort((a, c) => c.CompareTo(a));
            }
        }
        return true;
    }

    /// <summary>Classify the dephased palindrome of a hopping Hamiltonian without building the 2^N
    /// Liouvillian. Soft / Hard for a bit_b-homogeneous hopping H; <see cref="Verdict.OutOfScope"/> for
    /// a mixed Klein cell or a Hamiltonian with no hopping.</summary>
    public static Verdict Classify(IReadOnlyList<PauliTerm> terms, int n)
    {
        var masks = FlipMasks(terms, n);
        if (masks.Count == 0 || !IsBitBHomogeneous(terms)) return Verdict.OutOfScope;
        return MaskSetIsBipartite(masks) ? Verdict.Soft : Verdict.Hard;
    }
}
