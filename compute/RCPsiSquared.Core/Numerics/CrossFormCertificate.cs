namespace RCPsiSquared.Core.Numerics;

/// <summary>The GF(p) cross-form certificate slice: an independent C# transcription of the
/// F127 six-angle cross form 𝔉 evaluated over F_p (p ≡ 1 mod 4 so i = √−1 exists), with the
/// Conway-Jones double-constraint variety entered through the quadratic z₃² + S·z₃ + 1 = 0
/// (both roots, product 1, solved by Tonelli-Shanks). Mirrors the committed Python
/// <c>simulations/cross_triple_orthogonality.py</c> (<c>_ModPField</c>, <c>_cross_form_generic</c>,
/// <c>certify_cross_form_gfp</c>) so the live witness recomputes a genuinely discriminating slice
/// of F127 at inspect time: on-variety points must evaluate to 0, and the off-variety CONTROL
/// must be nonzero (without the control a zero function would pass).
///
/// <para>Determinism: sampling uses a fixed xorshift64* stream seeded from the prime, never
/// <c>System.Random</c> and never the clock, so an inspect run is reproducible bit-for-bit.
/// This is a SLICE, not the wall: the F127 proof object is the 527/527 full-grid + CRT sweep
/// (<c>simulations/grid_proof_sweep.py --assert</c>); this class is the independent second
/// implementation the code-trust caveat asks for, at sample scale.</para></summary>
public static class CrossFormCertificate
{
    /// <summary>Deterministic xorshift64* stream (never System.Random: its algorithm is not
    /// guaranteed stable across .NET versions, and a witness must reproduce bit-for-bit).</summary>
    private sealed class Stream(ulong seed)
    {
        private ulong _s = seed == 0 ? 0x9E3779B97F4A7C15UL : seed;

        public long Next(long loInclusive, long hiExclusive)
        {
            _s ^= _s >> 12; _s ^= _s << 25; _s ^= _s >> 27;
            ulong r = _s * 0x2545F4914F6CDD1DUL;
            return loInclusive + (long)(r % (ulong)(hiExclusive - loInclusive));
        }
    }

    public static long PowMod(long a, long e, long p)
    {
        a %= p; if (a < 0) a += p;
        long r = 1;
        while (e > 0)
        {
            if ((e & 1) == 1) r = r * a % p;
            a = a * a % p;
            e >>= 1;
        }
        return r;
    }

    public static long Inv(long a, long p) => PowMod(a, p - 2, p);

    /// <summary>Tonelli-Shanks square root mod an odd prime p; −1 if a is a non-residue.</summary>
    public static long Tonelli(long a, long p)
    {
        a %= p; if (a < 0) a += p;
        if (a == 0) return 0;
        if (PowMod(a, (p - 1) / 2, p) != 1) return -1;
        if (p % 4 == 3) return PowMod(a, (p + 1) / 4, p);
        long q = p - 1, s = 0;
        while (q % 2 == 0) { q /= 2; s++; }
        long z = 2;
        while (PowMod(z, (p - 1) / 2, p) != p - 1) z++;
        long m = s, c = PowMod(z, q, p), t = PowMod(a, q, p), r = PowMod(a, (q + 1) / 2, p);
        while (t != 1)
        {
            long i = 0, t2 = t;
            while (t2 != 1) { t2 = t2 * t2 % p; i++; }
            long b = PowMod(c, 1L << (int)(m - i - 1), p);
            m = i; c = b * b % p;
            t = t * c % p;
            r = r * b % p;
        }
        return r;
    }

    /// <summary>i = √−1 mod p (requires p ≡ 1 mod 4).</summary>
    public static long SqrtMinusOne(long p) => Tonelli(p - 1, p);

    /// <summary>The two roots of z² + S·z + 1 = 0 mod p for S = z₁+1/z₁+z₂+1/z₂ (the
    /// Conway-Jones constraint Σcos = 0 in z-coordinates); null if the discriminant is a
    /// non-residue or a root is 0. The roots multiply to 1 (checked).</summary>
    public static (long R1, long R2)? BothRoots(long z1, long z2, long p)
    {
        long s = (z1 + Inv(z1, p) + z2 + Inv(z2, p)) % p;
        long d = Tonelli(((s * s - 4) % p + p) % p, p);
        if (d < 0) return null;
        long half = Inv(2, p);
        long r1 = ((p - s + d) % p) * half % p;
        long r2 = (((p - s - d) % p + p) % p) * half % p;
        if (r1 == 0 || r2 == 0) return null;
        if (r1 * r2 % p != 1) throw new InvalidOperationException("root product != 1");
        return (r1, r2);
    }

    /// <summary>𝔉 evaluated over F_p at z-coordinates zs = (z₁,z₂,z₃,w₁,w₂,w₃), the direct port
    /// of the Python <c>_cross_form_generic</c> (one transcription, used on-variety and for the
    /// control). Throws <see cref="DivideByZeroException"/> on a pole (the caller skips).</summary>
    public static long Evaluate(long p, long iRoot, long[] zs)
    {
        long Mul(long a, long b) => a * b % p;
        long Add(long a, long b) => (a + b) % p;
        long Sub(long a, long b) => ((a - b) % p + p) % p;
        long Div(long a, long b)
        {
            if (b % p == 0) throw new DivideByZeroException();
            return Mul(a, Inv(b, p));
        }

        long IPow(long x, int k)
        {
            if (k < 0) return IPow(Div(1, x), -k);
            long r = 1;
            for (int t = 0; t < k; t++) r = Mul(r, x);
            return r;
        }

        long SinV(int j) => Div(Sub(zs[j], Div(1, zs[j])), Mul(2, iRoot));

        long Mono(int[] d)
        {
            long r = 1;
            for (int j = 0; j < 6; j++)
                if (d[j] != 0) r = Mul(r, IPow(zs[j], d[j]));
            return r;
        }

        long CotHalf(int[] d)
        {
            long z = Mono(d);
            return Mul(iRoot, Div(Add(z, 1), Sub(z, 1)));
        }

        int[] AddD(params int[][] ds)
        {
            var r = new int[6];
            foreach (var d in ds)
                for (int j = 0; j < 6; j++) r[j] += d[j];
            return r;
        }
        int[] Neg(int[] d)
        {
            var r = new int[6];
            for (int j = 0; j < 6; j++) r[j] = -d[j];
            return r;
        }
        int[] E(int j)
        {
            var r = new int[6];
            r[j] = 1;
            return r;
        }

        (int[] Psum, int[] Pdif, long Alpha, long Beta, long Lap)[] Pieces(int off)
        {
            var res = new (int[], int[], long, long, long)[3];
            for (int i = 0; i < 3; i++)
            {
                int j = -1, l = -1;
                for (int t = 0; t < 3; t++)
                    if (t != i) { if (j < 0) j = t; else l = t; }
                res[i] = (AddD(E(off + j), E(off + l)), AddD(E(off + l), Neg(E(off + j))),
                          Sub(SinV(off + l), SinV(off + j)),
                          Sub(0, Add(SinV(off + l), SinV(off + j))),
                          i % 2 == 0 ? 1 : p - 1);
            }
            return res;
        }

        var pa = Pieces(0);
        var pb = Pieces(3);
        long quarter = Inv(4, p);

        long Xh(int[] mu, int[] xi, int[] up) =>
            Mul(quarter, Sub(Add(CotHalf(AddD(mu, xi, Neg(up))), CotHalf(AddD(mu, Neg(xi), up))),
                             Add(CotHalf(AddD(mu, xi, up)), CotHalf(AddD(mu, Neg(xi), Neg(up))))));

        long tot = 0;
        for (int i = 0; i < 3; i++)
        {
            var pi = pa[i];
            for (int j = 0; j < 3; j++)
            {
                var pj = pb[j];
                long Xs(int[] mu)
                {
                    long r = Mul(Mul(pi.Alpha, pj.Alpha), Xh(mu, pi.Psum, pj.Psum));
                    r = Add(r, Mul(Mul(pi.Alpha, pj.Beta), Xh(mu, pi.Psum, pj.Pdif)));
                    r = Add(r, Mul(Mul(pi.Beta, pj.Alpha), Xh(mu, pi.Pdif, pj.Psum)));
                    return Add(r, Mul(Mul(pi.Beta, pj.Beta), Xh(mu, pi.Pdif, pj.Pdif)));
                }
                var muP = AddD(E(i), E(3 + j));
                var muM = AddD(E(i), Neg(E(3 + j)));
                long term = Sub(Mul(CotHalf(muP), Xs(muP)), Mul(CotHalf(muM), Xs(muM)));
                tot = Add(tot, Mul(pi.Lap, Mul(pj.Lap, term)));
            }
        }
        return Div(tot, 2);
    }

    /// <summary>One prime's certificate slice: <paramref name="samples"/> deterministic variety
    /// points (all four (z₃, w₃) root combinations each) plus one off-variety control per point.
    /// Returns (onVarietyEvaluations, nonZeroOnVariety, nonZeroControls).</summary>
    public static (int OnVariety, int BadOnVariety, int NonZeroControls, int Samples)
        CertifySlice(long p, int samples)
    {
        long iRoot = SqrtMinusOne(p);
        if (iRoot < 0 || iRoot * iRoot % p != p - 1)
            throw new InvalidOperationException($"p = {p} is not split (need p = 1 mod 4)");
        var stream = new Stream((ulong)p * 0x9E3779B97F4A7C15UL);
        int good = 0, bad = 0, ctrl = 0, evals = 0;
        while (good < samples)
        {
            long z1 = stream.Next(2, p), z2 = stream.Next(2, p);
            long w1 = stream.Next(2, p), w2 = stream.Next(2, p);
            long wCtrl = stream.Next(2, p);
            var rz = BothRoots(z1, z2, p);
            var rw = BothRoots(w1, w2, p);
            if (rz is null || rw is null) continue;
            try
            {
                foreach (long z3 in new[] { rz.Value.R1, rz.Value.R2 })
                    foreach (long w3 in new[] { rw.Value.R1, rw.Value.R2 })
                    {
                        evals++;
                        if (Evaluate(p, iRoot, new[] { z1, z2, z3, w1, w2, w3 }) != 0) bad++;
                    }
                if (Evaluate(p, iRoot, new[] { z1, z2, rz.Value.R1, w1, w2, wCtrl }) != 0) ctrl++;
            }
            catch (DivideByZeroException)
            {
                continue;   // a pole in the sample: skip, like the Python certificate
            }
            good++;
        }
        return (evals, bad, ctrl, good);
    }
}
