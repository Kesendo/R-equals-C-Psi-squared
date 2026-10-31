namespace RCPsiSquared.Core.Numerics;

/// <summary>A monomial x_0^e_0 ⋯ x_31^e_31 in at most <see cref="MaxVariables"/> variables, every exponent at most
/// <see cref="MaxExponent"/>, packed eight exponents to a 64-bit word (variable v in word v/8, byte v mod 8). The
/// packing keeps each byte below 128, so divisibility, quotient and lcm are word operations with no borrow
/// across bytes. Ordered by grevlex (<see cref="CompareGrevlex"/>), the order the Gröbner engine
/// <see cref="GroebnerTraceModP"/> works in.</summary>
public readonly struct PackedMonomial : IEquatable<PackedMonomial>
{
    public const int MaxVariables = 32;
    public const int MaxExponent = 127;
    private const ulong High = 0x8080808080808080UL;

    private readonly ulong _w0, _w1, _w2, _w3;

    /// <summary>The total degree Σ e_v.</summary>
    public int Degree { get; }

    private PackedMonomial(ulong w0, ulong w1, ulong w2, ulong w3, int degree)
    {
        _w0 = w0; _w1 = w1; _w2 = w2; _w3 = w3; Degree = degree;
    }

    /// <summary>The constant monomial 1.</summary>
    public static readonly PackedMonomial One = default;

    public static PackedMonomial Variable(int v, int power = 1)
    {
        if (v < 0 || v >= MaxVariables) throw new ArgumentOutOfRangeException(nameof(v), $"variable {v} outside 0..{MaxVariables - 1}");
        if (power < 0 || power > MaxExponent) throw new ArgumentOutOfRangeException(nameof(power), $"exponent {power} outside 0..{MaxExponent}");
        ulong w = (ulong)power << ((v & 7) * 8);
        return (v >> 3) switch
        {
            0 => new PackedMonomial(w, 0, 0, 0, power),
            1 => new PackedMonomial(0, w, 0, 0, power),
            2 => new PackedMonomial(0, 0, w, 0, power),
            _ => new PackedMonomial(0, 0, 0, w, power),
        };
    }

    private ulong Word(int k) => k switch { 0 => _w0, 1 => _w1, 2 => _w2, _ => _w3 };

    /// <summary>The exponent of variable v.</summary>
    public int this[int v] =>
        v is >= 0 and < MaxVariables
            ? (int)((Word(v >> 3) >> ((v & 7) * 8)) & 0xFF)
            : throw new ArgumentOutOfRangeException(nameof(v), $"variable {v} outside 0..{MaxVariables - 1}");

    public static PackedMonomial operator *(PackedMonomial a, PackedMonomial b)
    {
        ulong w0 = a._w0 + b._w0, w1 = a._w1 + b._w1, w2 = a._w2 + b._w2, w3 = a._w3 + b._w3;
        // both bytes below 128, so a sum never carries; a sum at or past 128 sets the byte's high bit
        if (((w0 | w1 | w2 | w3) & High) != 0)
            throw new OverflowException($"a monomial exponent passed {MaxExponent}");
        return new PackedMonomial(w0, w1, w2, w3, a.Degree + b.Degree);
    }

    // byte x ≥ byte y, for every byte, when all bytes are below 128: (x | 0x80) − y keeps the high bit exactly then
    private static bool BytesAtLeast(ulong x, ulong y) => (((x | High) - y) & High) == High;

    /// <summary>Whether this monomial divides m.</summary>
    public bool Divides(PackedMonomial m) =>
        Degree <= m.Degree && BytesAtLeast(m._w0, _w0) && BytesAtLeast(m._w1, _w1)
        && BytesAtLeast(m._w2, _w2) && BytesAtLeast(m._w3, _w3);

    /// <summary>m / this; the caller has checked <see cref="Divides"/>.</summary>
    public PackedMonomial DivideInto(PackedMonomial m) =>
        new(m._w0 - _w0, m._w1 - _w1, m._w2 - _w2, m._w3 - _w3, m.Degree - Degree);

    private static ulong ByteMax(ulong x, ulong y)
    {
        ulong mask = ((((x | High) - y) & High) >> 7) * 0xFF;   // 0xFF in the bytes where x ≥ y
        return (x & mask) | (y & ~mask);
    }

    private static int ByteSum(ulong w)
    {
        int s = 0;
        for (; w != 0; w >>= 8) s += (int)(w & 0xFF);
        return s;
    }

    public static PackedMonomial Lcm(PackedMonomial a, PackedMonomial b)
    {
        ulong w0 = ByteMax(a._w0, b._w0), w1 = ByteMax(a._w1, b._w1), w2 = ByteMax(a._w2, b._w2), w3 = ByteMax(a._w3, b._w3);
        return new PackedMonomial(w0, w1, w2, w3, ByteSum(w0) + ByteSum(w1) + ByteSum(w2) + ByteSum(w3));
    }

    /// <summary>Whether a and b share no variable (Buchberger's first criterion).</summary>
    public static bool AreCoprime(PackedMonomial a, PackedMonomial b) => Lcm(a, b).Degree == a.Degree + b.Degree;

    /// <summary>Graded reverse lexicographic order with x_0 &gt; x_1 &gt; ⋯: the higher degree is larger; at equal
    /// degree the monomial with the smaller exponent at the last variable where they differ is larger.</summary>
    public static int CompareGrevlex(PackedMonomial a, PackedMonomial b)
    {
        if (a.Degree != b.Degree) return a.Degree.CompareTo(b.Degree);
        for (int k = 3; k >= 0; k--)
        {
            ulong x = a.Word(k), y = b.Word(k);
            if (x == y) continue;
            int shift = (63 - System.Numerics.BitOperations.LeadingZeroCount(x ^ y)) & ~7;
            ulong ex = (x >> shift) & 0xFF, ey = (y >> shift) & 0xFF;
            return ex < ey ? 1 : -1;
        }
        return 0;
    }

    /// <summary>Grevlex, largest first: the order polynomial terms are stored in.</summary>
    public static readonly IComparer<PackedMonomial> Descending =
        Comparer<PackedMonomial>.Create((a, b) => CompareGrevlex(b, a));

    public bool Equals(PackedMonomial other) =>
        _w0 == other._w0 && _w1 == other._w1 && _w2 == other._w2 && _w3 == other._w3;

    public override bool Equals(object? obj) => obj is PackedMonomial m && Equals(m);
    public override int GetHashCode() => HashCode.Combine(_w0, _w1, _w2, _w3);
    public static bool operator ==(PackedMonomial a, PackedMonomial b) => a.Equals(b);
    public static bool operator !=(PackedMonomial a, PackedMonomial b) => !a.Equals(b);

    /// <summary>The monomial in the given variable names, "1" for the constant.</summary>
    public string ToString(IReadOnlyList<string> names)
    {
        var parts = new List<string>();
        for (int v = 0; v < MaxVariables; v++)
        {
            int e = this[v];
            if (e == 0) continue;
            string name = v < names.Count ? names[v] : $"x{v}";
            parts.Add(e == 1 ? name : $"{name}^{e}");
        }
        return parts.Count == 0 ? "1" : string.Join("·", parts);
    }

    public override string ToString() => ToString(Array.Empty<string>());
}
