using System;
using System.Collections.Generic;
using System.Linq;

namespace RCPsiSquared.Core.Numerics;

/// <summary>Frobenius cycle types of a Gaussian-integer polynomial reduced modulo a split
/// prime p ≡ 1 (mod 4), via distinct-degree factorisation over F_p. The prime splits in
/// Z[i] (residue field F_p), i ↦ r with r²≡−1; by Dedekind, when the reduction is squarefree
/// the multiset of irreducible-factor degrees is the cycle type of a Frobenius element of the
/// Galois group over Q(i). This is the from-below engine that certifies the F89 octic group
/// (Gal = S_8); reusable across the H_B-mixed factors (path-3 octic, path-4/5/6 degrees
/// 18/32/53). Mirrors simulations/f89_path3_octic_galois.py.
///
/// <para>Polynomials are int[] coefficient arrays, index = power (lowest-first), reduced
/// to [0, p). Cost is trivial at these sizes (degree ≤ 53, p &lt; a few hundred).</para></summary>
public static class OcticGaloisCertificate
{
    /// <summary>The cycle type (irreducible-factor degrees, sorted descending) of the
    /// Gaussian-integer polynomial Σ_k (re[k] + i·im[k])·x^k reduced modulo a prime 𝔭 | p.
    /// Returns null when p is not usable: p ≢ 1 (mod 4) (does not split in Z[i]), the
    /// reduction drops degree, or it is not squarefree (p | disc ⟹ Dedekind n/a).</summary>
    public static int[]? CycleType(long[] re, long[] im, int p)
    {
        if (re.Length != im.Length) throw new ArgumentException("re and im must have equal length.");
        if (p < 5 || p % 4 != 1) return null;          // need a prime that splits in Z[i]
        int? root = SqrtMinusOne(p);
        if (root is null) return null;
        long r = root.Value;

        int topPower = re.Length - 1;
        var f = new int[re.Length];
        for (int k = 0; k < re.Length; k++)            // a + b·i  ↦  a + b·r  (mod p)
            f[k] = (int)ScalarMod(re[k] + im[k] * r, p);
        if (Degree(f) != topPower) return null;        // leading coefficient vanished: degree drop

        var fp = Monic(f, p);
        if (Degree(Gcd(fp, Derivative(fp, p), p)) > 0) return null;  // not squarefree ⟹ p | disc

        return DistinctDegreeFactorisation(fp, p).OrderByDescending(d => d).ToArray();
    }

    /// <summary>Multiset of irreducible-factor degrees of a squarefree monic f over F_p.</summary>
    private static List<int> DistinctDegreeFactorisation(int[] f, int p)
    {
        var degrees = new List<int>();
        var work = (int[])f.Clone();
        var xqi = new[] { 0, 1 };                       // x^(p^0) = x
        for (int d = 1; Degree(work) >= 2 * d; d++)
        {
            xqi = PowMod(xqi, p, work, p);              // x^(p^d) mod work (Frobenius iterate)
            var g = Gcd(work, Sub(xqi, new[] { 0, 1 }, p), p);   // gcd(work, x^(p^d) − x)
            int gd = Degree(g);
            if (gd > 0)
            {
                for (int c = 0; c < gd / d; c++) degrees.Add(d);
                work = Div(work, g, p);                 // exact: g | work
            }
        }
        if (Degree(work) > 0) degrees.Add(Degree(work));   // remainder is a single irreducible
        return degrees;
    }

    // ---- F_p polynomial arithmetic (index = power, lowest-first) ----

    private static long ScalarMod(long a, int p) { long m = a % p; return m < 0 ? m + p : m; }

    private static int Degree(int[] a)
    {
        for (int i = a.Length - 1; i >= 0; i--) if (a[i] != 0) return i;
        return -1;                                       // zero polynomial
    }

    private static int[] Trim(int[] a)
    {
        int d = Degree(a);
        if (d < 0) return new[] { 0 };
        var r = new int[d + 1];
        Array.Copy(a, r, d + 1);
        return r;
    }

    private static int Inv(int a, int p)                 // a^(p−2) mod p (Fermat)
    {
        long result = 1, b = ScalarMod(a, p), e = p - 2;
        while (e > 0) { if ((e & 1) == 1) result = result * b % p; b = b * b % p; e >>= 1; }
        return (int)result;
    }

    private static int[] Monic(int[] a, int p)
    {
        int d = Degree(a);
        if (d < 0) return new[] { 0 };
        int inv = Inv(a[d], p);
        var r = new int[d + 1];
        for (int i = 0; i <= d; i++) r[i] = (int)ScalarMod((long)a[i] * inv, p);
        return r;
    }

    private static int[] Sub(int[] a, int[] b, int p)
    {
        int n = Math.Max(a.Length, b.Length);
        var r = new int[n];
        for (int i = 0; i < n; i++)
        {
            long va = i < a.Length ? a[i] : 0, vb = i < b.Length ? b[i] : 0;
            r[i] = (int)ScalarMod(va - vb, p);
        }
        return Trim(r);
    }

    private static int[] Mul(int[] a, int[] b, int p)
    {
        int da = Degree(a), db = Degree(b);
        if (da < 0 || db < 0) return new[] { 0 };
        var r = new long[da + db + 1];
        for (int i = 0; i <= da; i++)
            if (a[i] != 0)
                for (int j = 0; j <= db; j++) r[i + j] = (r[i + j] + (long)a[i] * b[j]) % p;
        var ri = new int[r.Length];
        for (int k = 0; k < r.Length; k++) ri[k] = (int)r[k];
        return Trim(ri);
    }

    /// <summary>Remainder a mod f.</summary>
    private static int[] Mod(int[] a, int[] f, int p)
    {
        int df = Degree(f);
        if (df < 0) throw new DivideByZeroException("polynomial modulus is zero.");
        var rem = new int[Math.Max(Degree(a) + 1, 1)];
        Array.Copy(a, rem, Math.Min(a.Length, rem.Length));
        int invLead = Inv(f[df], p);
        for (int i = Degree(rem); i >= df; i = Degree(rem))
        {
            int factor = (int)ScalarMod((long)rem[i] * invLead, p);
            for (int j = 0; j <= df; j++)
                rem[i - df + j] = (int)ScalarMod(rem[i - df + j] - (long)factor * f[j], p);
            if (Degree(rem) >= i) break;                 // safety: degree must drop each step
        }
        return Trim(rem);
    }

    /// <summary>Exact quotient a / f (caller guarantees f | a).</summary>
    private static int[] Div(int[] a, int[] f, int p)
    {
        int df = Degree(f), da = Degree(a);
        if (da < df) return new[] { 0 };
        var rem = new int[da + 1];
        Array.Copy(a, rem, Math.Min(a.Length, rem.Length));
        var quo = new int[da - df + 1];
        int invLead = Inv(f[df], p);
        for (int i = da; i >= df; i--)
        {
            int coef = (int)ScalarMod((long)rem[i] * invLead, p);
            quo[i - df] = coef;
            if (coef != 0)
                for (int j = 0; j <= df; j++)
                    rem[i - df + j] = (int)ScalarMod(rem[i - df + j] - (long)coef * f[j], p);
        }
        return Trim(quo);
    }

    private static int[] Gcd(int[] a, int[] b, int p)
    {
        var x = Trim(a);
        var y = Trim(b);
        while (Degree(y) >= 0)                            // y != 0
        {
            var rem = Mod(x, y, p);
            x = y;
            y = rem;
        }
        return Monic(x, p);
    }

    private static int[] Derivative(int[] a, int p)
    {
        if (Degree(a) <= 0) return new[] { 0 };
        var r = new int[a.Length - 1];
        for (int i = 1; i < a.Length; i++) r[i - 1] = (int)ScalarMod((long)a[i] * i, p);
        return Trim(r);
    }

    /// <summary>base^exp mod modulus over F_p, by binary exponentiation.</summary>
    private static int[] PowMod(int[] @base, long exp, int[] modulus, int p)
    {
        var result = new[] { 1 };
        var b = Mod(@base, modulus, p);
        while (exp > 0)
        {
            if ((exp & 1) == 1) result = Mod(Mul(result, b, p), modulus, p);
            b = Mod(Mul(b, b, p), modulus, p);
            exp >>= 1;
        }
        return result;
    }

    private static int? SqrtMinusOne(int p)              // r with r² ≡ −1 (mod p)
    {
        for (int r = 1; r < p; r++) if ((long)r * r % p == p - 1) return r;
        return null;
    }
}
