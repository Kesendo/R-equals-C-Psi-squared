using System.Numerics;

namespace RCPsiSquared.Core.Numerics;

/// <summary>Arbitrary-precision rational using <see cref="BigInteger"/>. Stored in
/// reduced form (gcd(Numerator, |Denominator|) = 1, Denominator > 0). Used by the
/// F89 Chebyshev pipeline (<see cref="Symmetry.F89PathPolynomialPipeline"/>) for
/// runtime computation of (P_k, D_k) at arbitrary k beyond the int.MaxValue
/// tabulation boundary.</summary>
public readonly struct BigRational : IEquatable<BigRational>
{
    public BigInteger Numerator { get; }
    public BigInteger Denominator { get; }

    public static readonly BigRational Zero = new(BigInteger.Zero, BigInteger.One, skipReduce: true);
    public static readonly BigRational One = new(BigInteger.One, BigInteger.One, skipReduce: true);

    public BigRational(BigInteger numerator, BigInteger denominator)
    {
        if (denominator.IsZero) throw new DivideByZeroException("BigRational denominator must be non-zero.");
        if (denominator.Sign < 0) { numerator = -numerator; denominator = -denominator; }
        if (numerator.IsZero) { Numerator = BigInteger.Zero; Denominator = BigInteger.One; return; }
        var gcd = BigInteger.GreatestCommonDivisor(BigInteger.Abs(numerator), denominator);
        Numerator = numerator / gcd;
        Denominator = denominator / gcd;
    }

    private BigRational(BigInteger numerator, BigInteger denominator, bool skipReduce)
    {
        Numerator = numerator;
        Denominator = denominator;
    }

    public BigRational(long value) : this(new BigInteger(value), BigInteger.One, skipReduce: true) { }
    public BigRational(BigInteger value) : this(value, BigInteger.One, skipReduce: true) { }

    public static implicit operator BigRational(long value) => new(value);
    public static implicit operator BigRational(int value) => new(value);
    public static implicit operator BigRational(BigInteger value) => new(value);

    public bool IsZero => Numerator.IsZero;
    public bool IsInteger => Denominator.IsOne;
    public int Sign => Numerator.Sign;

    public static BigRational operator +(BigRational a, BigRational b)
        => new(a.Numerator * b.Denominator + b.Numerator * a.Denominator, a.Denominator * b.Denominator);

    public static BigRational operator -(BigRational a, BigRational b)
        => new(a.Numerator * b.Denominator - b.Numerator * a.Denominator, a.Denominator * b.Denominator);

    public static BigRational operator -(BigRational a)
        => new(-a.Numerator, a.Denominator, skipReduce: true);

    public static BigRational operator *(BigRational a, BigRational b)
        => new(a.Numerator * b.Numerator, a.Denominator * b.Denominator);

    public static BigRational operator /(BigRational a, BigRational b)
    {
        if (b.IsZero) throw new DivideByZeroException();
        return new BigRational(a.Numerator * b.Denominator, a.Denominator * b.Numerator);
    }

    public static bool operator ==(BigRational a, BigRational b) => a.Equals(b);
    public static bool operator !=(BigRational a, BigRational b) => !a.Equals(b);

    public bool Equals(BigRational other) => Numerator == other.Numerator && Denominator == other.Denominator;
    public override bool Equals(object? obj) => obj is BigRational b && Equals(b);
    public override int GetHashCode() => HashCode.Combine(Numerator, Denominator);

    public override string ToString() =>
        Denominator.IsOne ? Numerator.ToString() : $"{Numerator}/{Denominator}";

    /// <summary>Returns LCM of two non-zero positive BigIntegers.</summary>
    public static BigInteger Lcm(BigInteger a, BigInteger b)
    {
        if (a.IsZero || b.IsZero) return BigInteger.Zero;
        return BigInteger.Abs(a / BigInteger.GreatestCommonDivisor(a, b) * b);
    }
}
