using System.Numerics;

namespace MirrorWorld;

/// <summary>A Pauli string as two bit masks, and its exact algebra: the product with its phase, and
/// whether two strings commute. Machinery, like ModP, not a GameObject: it states no physics and
/// appears in no Own or Inherited bucket.
///
/// Site l carries the letter (x_l, z_l): I = (0,0), X = (1,0), Z = (0,1), Y = (1,1), and the string
/// is P(x, z) = i^|x AND z| · X^x · Z^z, so every string is Hermitian, squares to the identity with
/// phase +1, and Y = i·X·Z site by site. The product follows from moving the Z part of the first
/// string past the X part of the second, which costs (-1)^|z_a AND x_b|:
///
///     P_a · P_b = i^k · P_(x_a XOR x_b, z_a XOR z_b),
///     k = |x_a AND z_a| + |x_b AND z_b| + 2·|z_a AND x_b| - |x_c AND z_c|   (mod 4),
///
/// with c the product string. Two strings commute exactly when |x_a AND z_b| + |z_a AND x_b| is even,
/// and then k is even; they anticommute exactly when it is odd, and then k is odd (the product of two
/// anticommuting Hermitian strings is anti-Hermitian). PauliStringTests judges both against dense
/// 2^N matrices built letter by letter, a route that never sees the masks.
///
/// Up to 32 sites, the bound EndCount's packing of a string into one 2N-bit word sets (the masks
/// alone would hold 64). The letter order of Parse and ToString is the site order,
/// site 0 first, the order PauliMode's Letters use.</summary>
public readonly struct PauliString : IEquatable<PauliString>
{
    public const int MaxSites = 32;

    public ulong X { get; }
    public ulong Z { get; }

    public PauliString(ulong x, ulong z)
    {
        if ((x >> MaxSites) != 0)
            throw new ArgumentOutOfRangeException(nameof(x), $"a string carries at most {MaxSites} sites");
        if ((z >> MaxSites) != 0)
            throw new ArgumentOutOfRangeException(nameof(z), $"a string carries at most {MaxSites} sites");
        X = x;
        Z = z;
    }

    public static readonly PauliString Identity = new(0, 0);

    public bool IsIdentity => X == 0 && Z == 0;

    /// <summary>The string read letter by letter, site 0 first; letters I, X, Y, Z only.</summary>
    public static PauliString Parse(string letters)
    {
        if (letters.Length == 0 || letters.Length > MaxSites)
            throw new ArgumentException($"a string has 1 to {MaxSites} letters; got {letters.Length}", nameof(letters));
        ulong x = 0, z = 0;
        for (int l = 0; l < letters.Length; l++)
        {
            switch (letters[l])
            {
                case 'I': break;
                case 'X': x |= 1UL << l; break;
                case 'Z': z |= 1UL << l; break;
                case 'Y': x |= 1UL << l; z |= 1UL << l; break;
                default: throw new ArgumentException($"'{letters[l]}' at site {l} is not a Pauli letter", nameof(letters));
            }
        }
        return new PauliString(x, z);
    }

    public char Letter(int site)
    {
        if (site < 0 || site >= MaxSites)
            throw new ArgumentOutOfRangeException(nameof(site), site, $"sites run from 0 to {MaxSites - 1}");
        bool xb = ((X >> site) & 1) == 1, zb = ((Z >> site) & 1) == 1;
        return xb ? (zb ? 'Y' : 'X') : (zb ? 'Z' : 'I');
    }

    public string ToString(int n)
    {
        if (n < 1 || n > MaxSites)
            throw new ArgumentOutOfRangeException(nameof(n), n, $"1 to {MaxSites} sites");
        if ((X >> n) != 0 || (Z >> n) != 0)
            throw new ArgumentOutOfRangeException(nameof(n), $"this string has letters beyond site {n - 1}");
        var c = new char[n];
        for (int l = 0; l < n; l++) c[l] = Letter(l);
        return new string(c);
    }

    public static bool Commute(PauliString a, PauliString b)
        => ((BitOperations.PopCount(a.X & b.Z) + BitOperations.PopCount(a.Z & b.X)) & 1) == 0;

    /// <summary>a · b = i^PhasePower · Product, PhasePower in 0..3.</summary>
    public static (PauliString Product, int PhasePower) Multiply(PauliString a, PauliString b)
    {
        ulong x = a.X ^ b.X, z = a.Z ^ b.Z;
        int k = BitOperations.PopCount(a.X & a.Z) + BitOperations.PopCount(b.X & b.Z)
              + 2 * BitOperations.PopCount(a.Z & b.X) - BitOperations.PopCount(x & z);
        return (new PauliString(x, z), ((k % 4) + 4) % 4);
    }

    public bool Equals(PauliString other) => X == other.X && Z == other.Z;
    public override bool Equals(object? obj) => obj is PauliString o && Equals(o);
    public override int GetHashCode() => HashCode.Combine(X, Z);
    public static bool operator ==(PauliString a, PauliString b) => a.Equals(b);
    public static bool operator !=(PauliString a, PauliString b) => !a.Equals(b);
}
