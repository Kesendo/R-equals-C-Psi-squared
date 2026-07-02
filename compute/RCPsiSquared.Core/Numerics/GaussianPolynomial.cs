using System;

namespace RCPsiSquared.Core.Numerics;

/// <summary>Exact polynomial arithmetic over Z[i] (coefficients lowest-first, index = power,
/// trimmed of trailing zeros; the zero polynomial is the empty array). Division is only by a
/// MONIC divisor (leading coefficient 1), which is exact over the ring — the regime the F89
/// isolation needs: charpoly ÷ AT-factor, both monic.</summary>
public static class GaussianPolynomial
{
    /// <summary>Degree of p (−1 for the zero polynomial), ignoring trailing zeros.</summary>
    public static int Degree(GaussianInteger[] p)
    {
        for (int i = p.Length - 1; i >= 0; i--)
            if (!p[i].Equals(GaussianInteger.Zero)) return i;
        return -1;
    }

    /// <summary>p with trailing zeros removed; the zero polynomial becomes the empty array.</summary>
    public static GaussianInteger[] Trim(GaussianInteger[] p)
    {
        int d = Degree(p);
        if (d < 0) return Array.Empty<GaussianInteger>();
        if (d == p.Length - 1) return p;
        var r = new GaussianInteger[d + 1];
        Array.Copy(p, r, d + 1);
        return r;
    }

    /// <summary>Polynomial long division by a MONIC divisor. Returns (quotient, remainder),
    /// both trimmed; remainder has degree &lt; deg(divisor). Throws if the divisor is zero or not
    /// monic (dividing by a non-unit leading coefficient would leave the ring Z[i]).</summary>
    public static (GaussianInteger[] Quotient, GaussianInteger[] Remainder) DivMod(
        GaussianInteger[] dividend, GaussianInteger[] divisor)
    {
        int dd = Degree(divisor);
        if (dd < 0) throw new DivideByZeroException("divisor is the zero polynomial.");
        if (!divisor[dd].Equals(GaussianInteger.One))
            throw new ArgumentException("divisor must be monic for exact division over Z[i].");

        var rem = (GaussianInteger[])dividend.Clone();
        int dr = Degree(rem);
        if (dr < dd) return (Array.Empty<GaussianInteger>(), Trim(rem));

        var quo = new GaussianInteger[dr - dd + 1];
        for (int i = dr; i >= dd; i--)
        {
            var coef = rem[i];                       // monic divisor ⟹ quotient coeff = leading rem
            quo[i - dd] = coef;
            if (!coef.Equals(GaussianInteger.Zero))
                for (int j = 0; j <= dd; j++)
                    rem[i - dd + j] = rem[i - dd + j] - coef * divisor[j];
        }
        return (Trim(quo), Trim(rem));
    }

    /// <summary>Product a·b over Z[i] (trimmed).</summary>
    public static GaussianInteger[] Multiply(GaussianInteger[] a, GaussianInteger[] b)
    {
        int da = Degree(a), db = Degree(b);
        if (da < 0 || db < 0) return Array.Empty<GaussianInteger>();
        var r = new GaussianInteger[da + db + 1];           // default-initialised to Zero
        for (int i = 0; i <= da; i++)
            for (int j = 0; j <= db; j++)
                r[i + j] += a[i] * b[j];
        return Trim(r);
    }

    /// <summary>True iff gcd(a, b) = 1 over Q(i), tested by the resultant (det of the Sylvester
    /// matrix over Z[i]): the resultant is nonzero iff a and b share no root. This is the
    /// exact, Gaussian-rational-free coprimality leg of the F_d isolation triple.</summary>
    public static bool AreCoprime(GaussianInteger[] a, GaussianInteger[] b)
    {
        int m = Degree(a), n = Degree(b);
        if (m < 0 || n < 0) return false;                   // the zero polynomial shares every root
        if (m == 0 || n == 0) return true;                  // a nonzero constant is coprime to all
        return !Determinant(BuildSylvester(a, b, m, n)).Equals(GaussianInteger.Zero);
    }

    /// <summary>The (m+n)×(m+n) Sylvester matrix of a (deg m) and b (deg n); its determinant is
    /// the resultant. Coefficients laid out highest-first along the shifted rows.</summary>
    private static GaussianInteger[,] BuildSylvester(GaussianInteger[] a, GaussianInteger[] b, int m, int n)
    {
        var s = new GaussianInteger[m + n, m + n];          // default-initialised to Zero
        for (int r = 0; r < n; r++)                         // n rows of a, shifted right by r
            for (int c = 0; c <= m; c++)
                s[r, r + c] = a[m - c];
        for (int r = 0; r < m; r++)                         // m rows of b, shifted right by r
            for (int c = 0; c <= n; c++)
                s[n + r, r + c] = b[n - c];
        return s;
    }

    /// <summary>det(M) over Z[i], read off the Berkowitz characteristic polynomial:
    /// charpoly(M)[0] = det(−M) = (−1)^n det(M).</summary>
    private static GaussianInteger Determinant(GaussianInteger[,] matrix)
    {
        int n = matrix.GetLength(0);
        var charpoly = GaussianMatrixCharpoly.Characteristic(matrix);
        var c0 = charpoly.Length > 0 ? charpoly[0] : GaussianInteger.Zero;
        return (n % 2 == 0) ? c0 : -c0;
    }

    /// <summary>The formal derivative dp/dx = Σ k·p_k x^{k−1} (trimmed; empty for a constant).</summary>
    public static GaussianInteger[] Derivative(GaussianInteger[] p)
    {
        int d = Degree(p);
        if (d <= 0) return Array.Empty<GaussianInteger>();
        var r = new GaussianInteger[d];
        for (int k = 1; k <= d; k++) r[k - 1] = (GaussianInteger)k * p[k];
        return Trim(r);
    }

    /// <summary>p evaluated at x (Horner), exact over Z[i].</summary>
    public static GaussianInteger Evaluate(GaussianInteger[] p, GaussianInteger x)
    {
        int d = Degree(p);
        if (d < 0) return GaussianInteger.Zero;
        var acc = p[d];
        for (int k = d - 1; k >= 0; k--) acc = acc * x + p[k];
        return acc;
    }

    /// <summary>The resultant Res(a, b) over Z[i] = det of the Sylvester matrix: zero iff a and b
    /// share a root. Zero if either is the zero polynomial; a nonzero constant c against b of degree
    /// n gives c^n. <see cref="AreCoprime"/> is exactly Res ≠ 0; this exposes the value, which the
    /// fold-resultant absence check (Res_λ(F_res, F_corner(−λ−2N)) as a function of q) samples.</summary>
    public static GaussianInteger Resultant(GaussianInteger[] a, GaussianInteger[] b)
    {
        int m = Degree(a), n = Degree(b);
        if (m < 0 || n < 0) return GaussianInteger.Zero;    // the zero polynomial shares every root
        if (m == 0) return Power(a[0], n);                  // Res(c, b) = c^deg(b)
        if (n == 0) return Power(b[0], m);
        return Determinant(BuildSylvester(a, b, m, n));
    }

    private static GaussianInteger Power(GaussianInteger baseVal, int exp)
    {
        var acc = GaussianInteger.One;
        for (int i = 0; i < exp; i++) acc *= baseVal;
        return acc;
    }

    /// <summary>The discriminant of a MONIC polynomial p: (−1)^{n(n−1)/2}·Res(p, p′) (the leading-
    /// coefficient division of the general formula is trivial for monic p, the regime of the char
    /// polys this serves). Zero iff p has a repeated root; the branch loci of F_res are its zeros
    /// in q. Throws if p is not monic.</summary>
    public static GaussianInteger Discriminant(GaussianInteger[] p)
    {
        int n = Degree(p);
        if (n < 1) return GaussianInteger.Zero;
        if (!p[n].Equals(GaussianInteger.One))
            throw new ArgumentException("Discriminant is implemented for monic polynomials only.");
        var res = Resultant(p, Derivative(p));
        return ((n * (n - 1) / 2) % 2 == 0) ? res : -res;
    }

    /// <summary>Composition p(α·x + β) over Z[i] (Horner in the linear argument). The fold-resultant
    /// needs F_corner(−λ − 2N) = ComposeLinear(F_corner, −1, −2N): its roots are −(corner spectrum)
    /// − 2N, so it shares a root with F_res (the (1,2) residual) exactly when some residual root's fold
    /// μ = −λ_A − 2N is a corner eigenvalue.</summary>
    public static GaussianInteger[] ComposeLinear(GaussianInteger[] p, GaussianInteger alpha, GaussianInteger beta)
    {
        int d = Degree(p);
        if (d < 0) return Array.Empty<GaussianInteger>();
        var arg = new[] { beta, alpha };                    // α·x + β (lowest-first)
        var acc = new[] { p[d] };                           // Horner from the leading coefficient
        for (int k = d - 1; k >= 0; k--)
        {
            acc = Multiply(acc, arg);
            if (acc.Length == 0) acc = new[] { GaussianInteger.Zero };
            acc[0] = acc[0] + p[k];
        }
        return Trim(acc);
    }
}
