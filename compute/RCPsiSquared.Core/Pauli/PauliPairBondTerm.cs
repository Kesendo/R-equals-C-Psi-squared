using System.Numerics;

namespace RCPsiSquared.Core.Pauli;

/// <summary>A single Pauli-pair bond term: σ_a on site i, σ_b on site i+1 (or whichever bond).
///
/// The canonical 2-letter representation used across the F-diagnostics (F49, F87, F80, F81,
/// F83). For raw N-letter terms see <see cref="PauliTerm"/>.
/// </summary>
public sealed record PauliPairBondTerm(PauliLetter LetterA, PauliLetter LetterB)
{
    /// <summary>Truly criterion at 2-body: #Y even AND #Z even on the (a, b) pair.
    /// Equivalent to "a == b OR {a, b} ⊆ {I, X}".</summary>
    public bool IsTruly
    {
        get
        {
            int nY = (LetterA == PauliLetter.Y ? 1 : 0) + (LetterB == PauliLetter.Y ? 1 : 0);
            int nZ = (LetterA == PauliLetter.Z ? 1 : 0) + (LetterB == PauliLetter.Z ? 1 : 0);
            return (nY & 1) == 0 && (nZ & 1) == 0;
        }
    }

    /// <summary>Π²-parity (Σ bit_b mod 2) of the bond term. 0 = Π²-even, 1 = Π²-odd.</summary>
    public int Pi2Parity => (LetterA.BitB() + LetterB.BitB()) & 1;
}

/// <summary>Algebraic 4-way Π²-class of a 2-body Pauli-pair Hamiltonian (the term-list
/// aggregate of <see cref="PauliPairBondTerm.IsTruly"/> and <see cref="PauliPairBondTerm.Pi2Parity"/>).
///
/// <para>This is the **algebraic** cut over a Hamiltonian, computed from the per-term
/// bit_b parities and #Y/#Z counts. It is complementary to the **spectral** 3-way
/// <c>RCPsiSquared.Diagnostics.F87.TrichotomyClass</c> (which classifies via the F1 residual
/// + spectrum-pairing test on the actual Liouvillian). The two cuts agree on all canonical
/// cases verified to date but use independent definitions; the algebraic one needs no L-build.</para>
/// </summary>
public enum Pi2Class
{
    /// <summary>Every term is truly (#Y even AND #Z even per term). M = 0 by Master Lemma.</summary>
    Truly,
    /// <summary>Every term has Π²-parity = 1 (Π²-odd). For 2-body bilinears these are the
    /// 4 cases (X,Y), (X,Z), (Y,X), (Z,X). F80 sign-walk applies.</summary>
    Pi2OddPure,
    /// <summary>Every term has Π²-parity = 0 (Π²-even) and at least one is non-truly.
    /// For 2-body bilinears the non-truly Π²-even cases are (Y,Z) and (Z,Y).</summary>
    Pi2EvenNonTruly,
    /// <summary>Term list contains both Π²-odd and Π²-even terms. Includes the
    /// truly + Π²-odd combinations (e.g. XX+XY).</summary>
    Mixed,
}

public static class PauliPairBondTermExtensions
{
    /// <summary>Convert a list of bond terms into the bilinear spec consumed by
    /// <see cref="PauliHamiltonian.Bilinear"/>: each (la, lb, coupling) entry.</summary>
    public static IReadOnlyList<(PauliLetter, PauliLetter, Complex)> ToBilinearSpec(
        this IEnumerable<PauliPairBondTerm> terms, double coupling) =>
        terms.Select(t => (t.LetterA, t.LetterB, (Complex)coupling)).ToList();

    /// <summary>Algebraic 4-way Π²-class of a term list. Uses only per-term bit_b parities
    /// and #Y/#Z counts; no L-build, no spectrum.</summary>
    public static Pi2Class Pi2ClassOf(this IReadOnlyList<PauliPairBondTerm> terms)
    {
        if (terms.Count == 0) return Pi2Class.Truly;

        bool hasOdd = false, hasEvenNonTruly = false, hasEvenTruly = false;
        foreach (var t in terms)
        {
            if (t.Pi2Parity == 1) hasOdd = true;
            else if (t.IsTruly) hasEvenTruly = true;
            else hasEvenNonTruly = true;
        }

        if (!hasOdd && !hasEvenNonTruly) return Pi2Class.Truly;
        if (hasOdd && (hasEvenNonTruly || hasEvenTruly)) return Pi2Class.Mixed;
        if (hasOdd) return Pi2Class.Pi2OddPure;
        return Pi2Class.Pi2EvenNonTruly;
    }
}
