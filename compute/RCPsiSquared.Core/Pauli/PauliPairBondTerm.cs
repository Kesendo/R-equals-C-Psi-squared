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

public static class PauliPairBondTermExtensions
{
    /// <summary>Convert a list of bond terms into the bilinear spec consumed by
    /// <see cref="PauliHamiltonian.Bilinear"/>: each (la, lb, coupling) entry.</summary>
    public static IReadOnlyList<(PauliLetter, PauliLetter, Complex)> ToBilinearSpec(
        this IEnumerable<PauliPairBondTerm> terms, double coupling) =>
        terms.Select(t => (t.LetterA, t.LetterB, (Complex)coupling)).ToList();
}
