namespace RCPsiSquared.Core.Numerics;

/// <summary>A Gaussian rational a + b·i (a, b ∈ ℚ, arbitrary precision over <see cref="BigRational"/>) —
/// the exact field ℚ(i). It is <see cref="GaussianInteger"/> with division: the ring Z[i] stays
/// division-free on purpose (Berkowitz plus exact polynomial division), but a Taylor recurrence divides
/// by the step index k+1 at every order, so it needs the field rather than the ring.
/// <para>Division is by the norm, b⁻¹ = conj(b)/|b|², which is exact because <see cref="BigRational"/>
/// is exact; the only failure mode is division by zero, which throws rather than returning a nonsense
/// value. Note that ℚ(i) has no order, so there is no comparison and no absolute value here: the
/// magnitude of a Gaussian rational is in general irrational, and <see cref="NormSquared"/> (the exact
/// rational |z|²) is what a caller compares instead.</para></summary>
public readonly struct GaussianRational : IEquatable<GaussianRational>
{
    public BigRational Re { get; }
    public BigRational Im { get; }

    public GaussianRational(BigRational re, BigRational im) { Re = re; Im = im; }

    public static readonly GaussianRational Zero = new(BigRational.Zero, BigRational.Zero);
    public static readonly GaussianRational One = new(BigRational.One, BigRational.Zero);
    public static readonly GaussianRational I = new(BigRational.Zero, BigRational.One);

    public bool IsZero => Re.IsZero && Im.IsZero;

    /// <summary>The exact rational |z|² = Re² + Im². The magnitude itself is in general irrational,
    /// so this is the quantity that stays inside the field.</summary>
    public BigRational NormSquared => Re * Re + Im * Im;

    public GaussianRational Conjugate => new(Re, -Im);

    public static GaussianRational operator +(GaussianRational a, GaussianRational b)
        => new(a.Re + b.Re, a.Im + b.Im);

    public static GaussianRational operator -(GaussianRational a, GaussianRational b)
        => new(a.Re - b.Re, a.Im - b.Im);

    public static GaussianRational operator -(GaussianRational a)
        => new(-a.Re, -a.Im);

    public static GaussianRational operator *(GaussianRational a, GaussianRational b)
        => new(a.Re * b.Re - a.Im * b.Im, a.Re * b.Im + a.Im * b.Re);

    public static GaussianRational operator /(GaussianRational a, GaussianRational b)
    {
        var norm = b.NormSquared;
        if (norm.IsZero) throw new DivideByZeroException("Division by zero in Q(i).");
        return new((a.Re * b.Re + a.Im * b.Im) / norm, (a.Im * b.Re - a.Re * b.Im) / norm);
    }

    public static implicit operator GaussianRational(BigRational q) => new(q, BigRational.Zero);
    public static implicit operator GaussianRational(int n) => new(n, BigRational.Zero);
    public static implicit operator GaussianRational(GaussianInteger z)
        => new(new BigRational(z.Re), new BigRational(z.Im));

    public static bool operator ==(GaussianRational a, GaussianRational b) => a.Equals(b);
    public static bool operator !=(GaussianRational a, GaussianRational b) => !a.Equals(b);

    public bool Equals(GaussianRational other) => Re == other.Re && Im == other.Im;
    public override bool Equals(object? obj) => obj is GaussianRational z && Equals(z);
    public override int GetHashCode() => HashCode.Combine(Re, Im);

    public override string ToString()
        => Im.IsZero ? Re.ToString()
         : Re.IsZero ? $"{Im}i"
         : Im.Sign < 0 ? $"{Re} - {-Im}i"
         : $"{Re} + {Im}i";
}
