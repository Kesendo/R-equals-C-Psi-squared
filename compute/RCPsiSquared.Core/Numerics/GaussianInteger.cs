using System.Numerics;

namespace RCPsiSquared.Core.Numerics;

/// <summary>A Gaussian integer a + b·i (a, b ∈ ℤ, arbitrary precision over BigInteger) — the
/// exact ring Z[i] the F89 path-k characteristic polynomial and F_d isolation live in. A value
/// type with exact +, −, ×; no division (the algorithms stay division-free, Berkowitz +
/// exact polynomial division). Implicit lifts from int/BigInteger make matrix literals clean.</summary>
public readonly record struct GaussianInteger(BigInteger Re, BigInteger Im)
{
    public static readonly GaussianInteger Zero = new(BigInteger.Zero, BigInteger.Zero);
    public static readonly GaussianInteger One = new(BigInteger.One, BigInteger.Zero);
    public static readonly GaussianInteger I = new(BigInteger.Zero, BigInteger.One);

    public static GaussianInteger operator +(GaussianInteger a, GaussianInteger b)
        => new(a.Re + b.Re, a.Im + b.Im);

    public static GaussianInteger operator -(GaussianInteger a, GaussianInteger b)
        => new(a.Re - b.Re, a.Im - b.Im);

    public static GaussianInteger operator -(GaussianInteger a)
        => new(-a.Re, -a.Im);

    public static GaussianInteger operator *(GaussianInteger a, GaussianInteger b)
        => new(a.Re * b.Re - a.Im * b.Im, a.Re * b.Im + a.Im * b.Re);

    public static implicit operator GaussianInteger(BigInteger n) => new(n, BigInteger.Zero);

    public static implicit operator GaussianInteger(int n) => new(n, BigInteger.Zero);
}
