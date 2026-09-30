using System.Numerics;

namespace RCPsiSquared.Core.Pauli;

/// <summary>A Pauli string as two bit masks, with its exact algebra: the product with its phase and
/// whether two strings commute, both read off popcounts, no matrix and no floating point.
///
/// Site l carries <see cref="PauliLetter"/>'s own two bits, the X part (bit_a) at bit l of
/// <see cref="X"/> and the Z part (bit_b) at bit l of <see cref="Z"/>, so I = (0,0), X = (1,0),
/// Z = (0,1), Y = (1,1). The string is P(x, z) = i^|x AND z| · X^x · Z^z: every string is
/// Hermitian, squares to the identity with phase +1, and Y = i·X·Z site by site, the matrices of
/// <see cref="PauliString.Build"/>. Moving the Z part of the first string past the X part of the
/// second costs (−1)^|z_a AND x_b|, which gives
///
/// <code>
///     P_a · P_b = i^k · P_c,   c = (x_a XOR x_b, z_a XOR z_b),
///     k = |x_a AND z_a| + |x_b AND z_b| + 2·|z_a AND x_b| − |x_c AND z_c|   (mod 4).
/// </code>
///
/// Two strings commute exactly when |x_a AND z_b| + |z_a AND x_b| is even (then k is even) and
/// anticommute exactly when it is odd (then k is odd). The tests judge both against
/// <see cref="PauliString.Build"/>, a dense route that never sees a mask.
///
/// <para>The same algebra is written once more in compute/MirrorWorld/PauliString.cs, the
/// sober base, which references nothing in RCPsiSquared.*; the two are independent by design.</para>
/// </summary>
public readonly record struct PauliMask(ulong X, ulong Z)
{
    /// <summary>The bound of the masks.</summary>
    public const int MaxSites = 64;

    public static readonly PauliMask Identity = new(0, 0);

    public bool IsIdentity => X == 0 && Z == 0;

    public static PauliMask FromLetters(IReadOnlyList<PauliLetter> letters)
    {
        if (letters.Count > MaxSites)
            throw new ArgumentOutOfRangeException(nameof(letters), letters.Count, $"at most {MaxSites} sites");
        ulong x = 0, z = 0;
        for (int l = 0; l < letters.Count; l++)
        {
            if (letters[l].BitA() == 1) x |= 1UL << l;
            if (letters[l].BitB() == 1) z |= 1UL << l;
        }
        return new PauliMask(x, z);
    }

    /// <summary>The string read letter by letter from I/X/Y/Z, site 0 first.</summary>
    public static PauliMask Parse(string letters)
        => FromLetters(letters.Select(PauliLetterExtensions.FromSymbol).ToArray());

    public PauliLetter Letter(int site)
    {
        if (site < 0 || site >= MaxSites)
            throw new ArgumentOutOfRangeException(nameof(site), site, $"sites run from 0 to {MaxSites - 1}");
        return (PauliLetter)((int)((X >> site) & 1) + 2 * (int)((Z >> site) & 1));
    }

    public PauliLetter[] ToLetters(int n)
    {
        if (n < 1 || n > MaxSites)
            throw new ArgumentOutOfRangeException(nameof(n), n, $"1 to {MaxSites} sites");
        if (n < MaxSites && ((X >> n) != 0 || (Z >> n) != 0))
            throw new ArgumentOutOfRangeException(nameof(n), n, $"this string has letters beyond site {n - 1}");
        return Enumerable.Range(0, n).Select(Letter).ToArray();
    }

    public string ToString(int n) => new(ToLetters(n).Select(l => l.Symbol()).ToArray());

    public static bool Commute(PauliMask a, PauliMask b)
        => ((BitOperations.PopCount(a.X & b.Z) + BitOperations.PopCount(a.Z & b.X)) & 1) == 0;

    /// <summary>a · b = i^PhasePower · Product, PhasePower in 0..3.</summary>
    public static (PauliMask Product, int PhasePower) Multiply(PauliMask a, PauliMask b)
    {
        ulong x = a.X ^ b.X, z = a.Z ^ b.Z;
        int k = BitOperations.PopCount(a.X & a.Z) + BitOperations.PopCount(b.X & b.Z)
              + 2 * BitOperations.PopCount(a.Z & b.X) - BitOperations.PopCount(x & z);
        return (new PauliMask(x, z), ((k % 4) + 4) % 4);
    }
}
