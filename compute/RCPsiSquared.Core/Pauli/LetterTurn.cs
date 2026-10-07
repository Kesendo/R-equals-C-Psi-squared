using System.Numerics;

namespace RCPsiSquared.Core.Pauli;

/// <summary>The cube turn that carries a dephasing letter to Z, applied exactly to Pauli letters.
///
/// <para>For a jump letter J there is a single-qubit Clifford U with U·J·U† = Z; conjugation by
/// U^⊗N permutes the letters with signs, so a Hamiltonian given as Pauli terms is turned exactly,
/// with no dense U·H·U†: the letters are relabelled and each coefficient picks up the product of the
/// signs. A Lindbladian with J-dephasing on every site is unitarily equivalent
/// to the Z-dephasing Lindbladian of the turned Hamiltonian, with the same rates and the same
/// spectrum.</para>
///
/// <para>The turns: for J = Z the identity; for J = X the Hadamard (X → Z, Z → X, Y → −Y); for
/// J = Y the quarter turn about X (X → X, Y → Z, Z → −Y), the inverse of the turn F108 Part 3 uses
/// to carry Z-dephasing to Y-dephasing. Both are proper rotations (each preserves XY = iZ). MirrorWorld's
/// <c>EndCount.TurnToZ</c> holds an independent copy with a cyclic Y turn; any proper turn gives
/// the same spectrum.</para></summary>
public static class LetterTurn
{
    /// <summary>The image of <paramref name="letter"/> under the turn that carries
    /// <paramref name="jump"/> to Z, with its sign.</summary>
    public static (PauliLetter Letter, int Sign) TurnToZ(PauliLetter letter, PauliLetter jump)
    {
        if (letter == PauliLetter.I) return (PauliLetter.I, 1);
        return jump switch
        {
            PauliLetter.Z => (letter, 1),
            PauliLetter.X => letter switch
            {
                PauliLetter.X => (PauliLetter.Z, 1),
                PauliLetter.Z => (PauliLetter.X, 1),
                _ => (PauliLetter.Y, -1),
            },
            PauliLetter.Y => letter switch
            {
                PauliLetter.X => (PauliLetter.X, 1),
                PauliLetter.Y => (PauliLetter.Z, 1),
                _ => (PauliLetter.Y, -1),
            },
            _ => throw new ArgumentException("the identity is not a dephasing letter", nameof(jump)),
        };
    }

    /// <summary>The term with every letter turned; the coefficient carries the product of the signs.</summary>
    public static PauliTerm Turn(PauliTerm term, PauliLetter jump)
    {
        if (term is null) throw new ArgumentNullException(nameof(term));
        var letters = new PauliLetter[term.N];
        int sign = 1;
        for (int i = 0; i < term.N; i++)
        {
            var (l, s) = TurnToZ(term.Letters[i], jump);
            letters[i] = l;
            sign *= s;
        }
        return new PauliTerm(letters, sign * term.Coefficient);
    }

    /// <summary>The Hamiltonian with every term turned, in the same term order.</summary>
    public static PauliHamiltonian Turn(PauliHamiltonian H, PauliLetter jump)
    {
        if (H is null) throw new ArgumentNullException(nameof(H));
        return new PauliHamiltonian(H.N, H.Terms.Select(t => Turn(t, jump)).ToList());
    }
}
