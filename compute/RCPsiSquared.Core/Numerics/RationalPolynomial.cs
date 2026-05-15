using System.Numerics;

namespace RCPsiSquared.Core.Numerics;

/// <summary>Univariate polynomial with <see cref="BigRational"/> coefficients,
/// stored low-to-high degree. Trailing zeros trimmed; zero polynomial = empty.
/// Used by the F89 Chebyshev pipeline.</summary>
public sealed class RationalPolynomial : IEquatable<RationalPolynomial>
{
    private readonly BigRational[] _coefs;

    public RationalPolynomial(IEnumerable<BigRational> coefsLowToHigh)
    {
        var list = coefsLowToHigh.ToList();
        while (list.Count > 0 && list[^1].IsZero) list.RemoveAt(list.Count - 1);
        _coefs = list.ToArray();
    }

    public RationalPolynomial(params BigRational[] coefsLowToHigh)
        : this((IEnumerable<BigRational>)coefsLowToHigh) { }

    /// <summary>Degree; -1 for the zero polynomial.</summary>
    public int Degree => _coefs.Length - 1;
    public bool IsZero => _coefs.Length == 0;

    public BigRational this[int degree] =>
        degree >= 0 && degree < _coefs.Length ? _coefs[degree] : BigRational.Zero;

    public IReadOnlyList<BigRational> CoefficientsLowToHigh => _coefs;

    public static readonly RationalPolynomial Zero = new(Array.Empty<BigRational>());
    public static readonly RationalPolynomial One = new(new[] { BigRational.One });

    /// <summary>Monomial coef·y^degree.</summary>
    public static RationalPolynomial Monomial(int degree, BigRational coef)
    {
        if (degree < 0) throw new ArgumentOutOfRangeException(nameof(degree));
        if (coef.IsZero) return Zero;
        var arr = new BigRational[degree + 1];
        for (int i = 0; i < degree; i++) arr[i] = BigRational.Zero;
        arr[degree] = coef;
        return new RationalPolynomial(arr);
    }

    /// <summary>Build (y - root) for a rational root.</summary>
    public static RationalPolynomial LinearFactor(BigRational root) =>
        new(new[] { -root, BigRational.One });

    public static RationalPolynomial operator +(RationalPolynomial a, RationalPolynomial b)
    {
        int n = Math.Max(a._coefs.Length, b._coefs.Length);
        var result = new BigRational[n];
        for (int i = 0; i < n; i++) result[i] = a[i] + b[i];
        return new RationalPolynomial(result);
    }

    public static RationalPolynomial operator -(RationalPolynomial a, RationalPolynomial b)
    {
        int n = Math.Max(a._coefs.Length, b._coefs.Length);
        var result = new BigRational[n];
        for (int i = 0; i < n; i++) result[i] = a[i] - b[i];
        return new RationalPolynomial(result);
    }

    public static RationalPolynomial operator -(RationalPolynomial a)
        => new(a._coefs.Select(c => -c));

    public static RationalPolynomial operator *(RationalPolynomial a, RationalPolynomial b)
    {
        if (a.IsZero || b.IsZero) return Zero;
        var result = new BigRational[a._coefs.Length + b._coefs.Length - 1];
        for (int i = 0; i < result.Length; i++) result[i] = BigRational.Zero;
        for (int i = 0; i < a._coefs.Length; i++)
            for (int j = 0; j < b._coefs.Length; j++)
                result[i + j] = result[i + j] + a._coefs[i] * b._coefs[j];
        return new RationalPolynomial(result);
    }

    public static RationalPolynomial operator *(BigRational scalar, RationalPolynomial p)
    {
        if (scalar.IsZero || p.IsZero) return Zero;
        return new RationalPolynomial(p._coefs.Select(c => scalar * c));
    }

    /// <summary>Euclidean division: returns (Quotient, Remainder) with
    /// Degree(Remainder) &lt; Degree(divisor).</summary>
    public (RationalPolynomial Quotient, RationalPolynomial Remainder) DivMod(RationalPolynomial divisor)
    {
        if (divisor.IsZero) throw new DivideByZeroException("Polynomial division by zero.");
        if (this.Degree < divisor.Degree) return (Zero, this);

        var remainder = _coefs.ToList();
        var quotientArr = new BigRational[Degree - divisor.Degree + 1];
        for (int i = 0; i < quotientArr.Length; i++) quotientArr[i] = BigRational.Zero;

        var divisorLead = divisor._coefs[^1];
        for (int i = remainder.Count - 1; i >= divisor.Degree; i--)
        {
            if (remainder[i].IsZero) continue;
            var coef = remainder[i] / divisorLead;
            int qDeg = i - divisor.Degree;
            quotientArr[qDeg] = coef;
            for (int j = 0; j <= divisor.Degree; j++)
                remainder[i - divisor.Degree + j] = remainder[i - divisor.Degree + j] - coef * divisor._coefs[j];
        }

        return (new RationalPolynomial(quotientArr), new RationalPolynomial(remainder));
    }

    public RationalPolynomial Mod(RationalPolynomial divisor) => DivMod(divisor).Remainder;

    public static implicit operator RationalPolynomial(BigRational coef) => Monomial(0, coef);

    /// <summary>LCM of denominators across all coefficients. Returns 1 for the
    /// zero polynomial.</summary>
    public BigInteger DenominatorLcm()
    {
        if (IsZero) return BigInteger.One;
        var lcm = _coefs[0].Denominator;
        for (int i = 1; i < _coefs.Length; i++)
            lcm = BigRational.Lcm(lcm, _coefs[i].Denominator);
        return lcm;
    }

    public bool Equals(RationalPolynomial? other)
    {
        if (other is null) return false;
        if (ReferenceEquals(this, other)) return true;
        if (_coefs.Length != other._coefs.Length) return false;
        for (int i = 0; i < _coefs.Length; i++)
            if (!_coefs[i].Equals(other._coefs[i])) return false;
        return true;
    }

    public override bool Equals(object? obj) => obj is RationalPolynomial p && Equals(p);
    public override int GetHashCode()
    {
        if (_coefs.Length == 0) return 0;
        return HashCode.Combine(_coefs.Length, _coefs[^1], _coefs[0]);
    }

    public override string ToString()
    {
        if (IsZero) return "0";
        var sb = new System.Text.StringBuilder();
        for (int i = _coefs.Length - 1; i >= 0; i--)
        {
            if (_coefs[i].IsZero) continue;
            if (sb.Length > 0) sb.Append(_coefs[i].Sign > 0 ? " + " : " - ");
            else if (_coefs[i].Sign < 0) sb.Append('-');
            var abs = _coefs[i].Sign < 0 ? -_coefs[i] : _coefs[i];
            if (!abs.Equals(BigRational.One) || i == 0) sb.Append(abs);
            if (i >= 1) sb.Append("·y");
            if (i >= 2) sb.Append('^').Append(i);
        }
        return sb.ToString();
    }
}
