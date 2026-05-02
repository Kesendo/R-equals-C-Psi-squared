namespace RCPsiSquared.Core.Pauli;

/// <summary>The four single-qubit Pauli operators, indexed by (bit_a, bit_b) ∈ {0,1}²
/// and packed as <c>a + 2·b</c>. The naming follows the F77/F81 framework convention:
///
///   I = (0, 0) = 0    bit_a = 0 (no XY-content)         bit_b = 0 (Π²-even)
///   X = (1, 0) = 1    bit_a = 1 (XY content)            bit_b = 0 (Π²-even)
///   Z = (0, 1) = 2    bit_a = 0 (no XY)                 bit_b = 1 (Π²-odd)
///   Y = (1, 1) = 3    bit_a = 1 (XY content)            bit_b = 1 (Π²-odd)
///
/// bit_a counts X+Y (the n_XY weight); bit_b counts Y+Z (the Π²-parity / n_YZ).
/// Y-parity (#Y mod 2) is independent only at k ≥ 3-body terms.
/// </summary>
public enum PauliLetter
{
    I = 0,
    X = 1,
    Z = 2,
    Y = 3,
}

public static class PauliLetterExtensions
{
    /// <summary>bit_a (XY-content): 1 for X and Y, 0 for I and Z.</summary>
    public static int BitA(this PauliLetter letter) => (int)letter & 1;

    /// <summary>bit_b (Π²-parity / YZ-content): 1 for Y and Z, 0 for I and X.</summary>
    public static int BitB(this PauliLetter letter) => ((int)letter >> 1) & 1;

    /// <summary>Single-character label: I/X/Y/Z.</summary>
    public static char Symbol(this PauliLetter letter) => letter switch
    {
        PauliLetter.I => 'I',
        PauliLetter.X => 'X',
        PauliLetter.Z => 'Z',
        PauliLetter.Y => 'Y',
        _ => '?',
    };

    public static PauliLetter FromSymbol(char c) => c switch
    {
        'I' or 'i' => PauliLetter.I,
        'X' or 'x' => PauliLetter.X,
        'Y' or 'y' => PauliLetter.Y,
        'Z' or 'z' => PauliLetter.Z,
        _ => throw new ArgumentException($"unknown Pauli letter: '{c}'"),
    };
}
