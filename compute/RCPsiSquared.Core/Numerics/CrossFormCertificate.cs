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

    /// <summary>The §3 core function T = Σ_{i,j} (−1)^{i+j} α_i(z) α_j(w)·cot((a_i+b_j)/2) over F_p,
    /// with α_i(z) = sin(a_l) − sin(a_j) on the ascending complement {j &lt; l}, sin x = (x − 1/x)/(2i),
    /// and cot((a_i+b_j)/2) = i(Z+1)/(Z−1), Z = z_i·w_j. This is an INDEPENDENT transcription (not a
    /// slice of <see cref="Evaluate"/>, which is the full 𝔉): T is the residue every sheet's nine
    /// events share (<see cref="SheetLattice"/>). Throws <see cref="DivideByZeroException"/> on a pole
    /// Z = 1 (the caller skips), which also covers the §3.3 singular cotangent sublocus z₃w₃ = 1.</summary>
    public static long EvaluateCoreT(long p, long iRoot, long[] zs)
    {
        long Mul(long a, long b) => a * b % p;
        long Add(long a, long b) => (a + b) % p;
        long Sub(long a, long b) => ((a - b) % p + p) % p;
        long Div(long a, long b)
        {
            if (b % p == 0) throw new DivideByZeroException();
            return Mul(a, Inv(b, p));
        }

        long Sin(long z) => Div(Sub(z, Inv(z, p)), Mul(2, iRoot));

        long Alpha(int off, int i)
        {
            int j = -1, l = -1;
            for (int t = 0; t < 3; t++)
                if (t != i) { if (j < 0) j = t; else l = t; }
            return Sub(Sin(zs[off + l]), Sin(zs[off + j]));
        }

        long CotHalf(long z) => Mul(iRoot, Div(Add(z, 1), Sub(z, 1)));

        long tot = 0;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                long term = Mul(Mul(Alpha(0, i), Alpha(3, j)), CotHalf(Mul(zs[i], zs[3 + j])));
                tot = Add(tot, (i + j) % 2 == 0 ? term : Sub(0, term));
            }
        return tot;
    }

    /// <summary>One prime's CORE-IDENTITY slice (PROOF_F127_RESIDUE_COLLAPSE §3): deterministic
    /// points on the three-constraint variety {Qz(z₃)=0, Qw(w₃)=0, z₁z₂z₃w₁w₂w₃=1}, sampled by
    /// choosing z₁, z₃, w₁, solving z₂ from Ca (z₂+1/z₂ = −(z₁+1/z₁)−(z₃+1/z₃)), setting the sheet
    /// w₃ = K/w₂ with K = 1/(z₁z₂z₃w₁), and solving w₂ from Cb, which reduces to
    /// (K+1)w₂² + K(w₁+1/w₁)w₂ + K(K+1) = 0. All three constraints are re-checked exactly per point
    /// (a violation is a construction bug and throws). T must vanish on-variety; the control
    /// re-randomizes w₁ (breaking Cb and the sheet), where T must be nonzero. Returns
    /// (onVarietyPoints, nonZeroOnVariety, nonZeroControls, controlsEvaluated).</summary>
    public static (int OnVarietyEvals, int BadOnVariety, int ControlsNonzero, int ControlsEvaluated)
        CertifyCoreIdentitySlice(long p, int samples)
    {
        long iRoot = SqrtMinusOne(p);
        if (iRoot < 0 || iRoot * iRoot % p != p - 1)
            throw new InvalidOperationException($"p = {p} is not split (need p = 1 mod 4)");
        var stream = new Stream((ulong)p * 0xD1B54A32D192ED03UL);   // independent of the 𝔉 slice's stream
        long inv2 = Inv(2, p);
        long CosSum(long a, long b, long c) =>
            ((a + Inv(a, p)) % p + (b + Inv(b, p)) % p + (c + Inv(c, p)) % p) % p;

        int good = 0, bad = 0, ctrlNonzero = 0, ctrlEval = 0;
        while (good < samples)
        {
            long z1 = stream.Next(2, p), z3 = stream.Next(2, p), w1 = stream.Next(2, p);
            long vz = ((p - (z1 + Inv(z1, p)) % p) + (p - (z3 + Inv(z3, p)) % p)) % p;
            long dz = Tonelli(((vz * vz - 4) % p + p) % p, p);
            if (dz < 0) continue;
            long z2 = (vz + dz) % p * inv2 % p;
            if (z2 == 0) continue;

            long k = Inv(z1 * z2 % p * z3 % p * w1 % p, p);
            long w1s = (w1 + Inv(w1, p)) % p;
            long aq = (k + 1) % p;
            if (aq == 0) continue;
            long bq = k * w1s % p;
            long cq = k * aq % p;
            long disc = ((bq * bq - 4 * aq % p * cq) % p + p) % p;
            long dw = Tonelli(disc, p);
            if (dw < 0) continue;
            long w2 = (p - bq + dw) % p * Inv(2 * aq % p, p) % p;
            if (w2 == 0) continue;
            long w3 = k * Inv(w2, p) % p;
            if (w3 == 0) continue;

            if (CosSum(z1, z2, z3) != 0 || CosSum(w1, w2, w3) != 0 ||
                z1 * z2 % p * z3 % p * w1 % p * w2 % p * w3 % p != 1)
                throw new InvalidOperationException("constructed point violates a variety constraint");

            long tval;
            try { tval = EvaluateCoreT(p, iRoot, new[] { z1, z2, z3, w1, w2, w3 }); }
            catch (DivideByZeroException) { continue; }   // a pole z_i w_j = 1: skip
            good++;
            if (tval != 0) bad++;

            for (int attempt = 0; attempt < 4; attempt++)
            {
                long w1c = stream.Next(2, p);
                long cval;
                try { cval = EvaluateCoreT(p, iRoot, new[] { z1, z2, z3, w1c, w2, w3 }); }
                catch (DivideByZeroException) { continue; }
                ctrlEval++;
                if (cval != 0) ctrlNonzero++;
                break;
            }
        }
        return (good, bad, ctrlNonzero, ctrlEval);
    }

    /// <summary>The §3 CLOSED FORM (2026-07-14, <c>simulations/f127_closed_form.py</c> gate G1)
    /// checked live at GENERIC points: T·P = ⅛[2cos s((e₁−f₁)² − 2sin²s) + sin s(Σsin2a + Σsin2b)]
    /// ·V_a·V_b is an UNCONDITIONAL identity, so no variety construction is needed; the strongest
    /// check is at free half-angle points ζᵢ, ωⱼ (zᵢ = ζᵢ², wⱼ = ωⱼ², e^{is} = Πζ·Πω). P, V_a, V_b
    /// are the Cauchy sine products; sin x = (m − 1/m)/(2i) for e^{ix} = m. The corruption control
    /// sign-flips the sin s term; it must BREAK at nearly every point, else the comparison is
    /// vacuous. Returns (points, mismatches, controlMismatches, controlEvals).</summary>
    public static (int Points, int Mismatches, int ControlMismatches, int ControlEvals)
        CertifyClosedFormSlice(long p, int samples)
    {
        long iRoot = SqrtMinusOne(p);
        if (iRoot < 0 || iRoot * iRoot % p != p - 1)
            throw new InvalidOperationException($"p = {p} is not split (need p = 1 mod 4)");
        var stream = new Stream((ulong)p * 0xA24BAED4963EE407UL);   // own stream, independent of the other slices
        long Mul(long a, long b) => a * b % p;
        long Add(long a, long b) => (a + b) % p;
        long Sub(long a, long b) => ((a - b) % p + p) % p;
        long inv2i = Inv(Mul(2, iRoot), p);
        long inv2 = Inv(2, p);
        long inv8 = Inv(8, p);
        long SinOf(long m) => Mul(Sub(m, Inv(m, p)), inv2i);
        long CosOf(long m) => Mul(Add(m, Inv(m, p)), inv2);

        int points = 0, mismatches = 0, ctrlMismatch = 0, ctrlEval = 0;
        while (points < samples)
        {
            var zeta = new long[3];
            var omega = new long[3];
            for (int t = 0; t < 3; t++) { zeta[t] = stream.Next(2, p); omega[t] = stream.Next(2, p); }
            var zs = new long[6];
            for (int t = 0; t < 3; t++) { zs[t] = Mul(zeta[t], zeta[t]); zs[3 + t] = Mul(omega[t], omega[t]); }

            long tVal;
            try { tVal = EvaluateCoreT(p, iRoot, zs); }
            catch (DivideByZeroException) { continue; }   // a cot pole zᵢwⱼ = 1 (⟺ P = 0): skip

            long pProd = 1, va = 1, vb = 1;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    pProd = Mul(pProd, SinOf(Mul(zeta[i], omega[j])));
            for (int i = 0; i < 3; i++)
                for (int j = i + 1; j < 3; j++)
                {
                    va = Mul(va, SinOf(Mul(zeta[i], Inv(zeta[j], p))));
                    vb = Mul(vb, SinOf(Mul(omega[i], Inv(omega[j], p))));
                }
            // degenerate sample: skip (pProd = 0 ⟺ a cot pole, already thrown above; the live
            // legs are va/vb = 0, the z_i = z_j coincidences)
            if (pProd == 0 || va == 0 || vb == 0) continue;

            long mu = 1;
            for (int t = 0; t < 3; t++) mu = Mul(Mul(mu, zeta[t]), omega[t]);
            long cosS = CosOf(mu), sinS = SinOf(mu);
            long e1 = 0, f1 = 0, s2a = 0, s2b = 0;
            for (int t = 0; t < 3; t++)
            {
                e1 = Add(e1, CosOf(zs[t]));
                f1 = Add(f1, CosOf(zs[3 + t]));
                s2a = Add(s2a, SinOf(Mul(zs[t], zs[t])));
                s2b = Add(s2b, SinOf(Mul(zs[3 + t], zs[3 + t])));
            }
            long dSq = Mul(Sub(e1, f1), Sub(e1, f1));
            long bracket = Add(Mul(Mul(2, cosS), Sub(dSq, Mul(2, Mul(sinS, sinS)))),
                               Mul(sinS, Add(s2a, s2b)));
            long rhs = Mul(Mul(Mul(inv8, bracket), va), vb);
            long lhs = Mul(tVal, pProd);
            points++;
            if (lhs != rhs) mismatches++;

            // corruption control: sign-flip the sin s term; must disagree with the LHS
            long bracketBad = Add(Mul(Mul(2, cosS), Sub(dSq, Mul(2, Mul(sinS, sinS)))),
                                  Sub(0, Mul(sinS, Add(s2a, s2b))));
            long rhsBad = Mul(Mul(Mul(inv8, bracketBad), va), vb);
            ctrlEval++;
            if (lhs != rhsBad) ctrlMismatch++;
        }
        return (points, mismatches, ctrlMismatch, ctrlEval);
    }

    /// <summary>Corollary 2 of the closed form (2026-07-14, NEW surface): T = 0 already on the
    /// SHARPER locus {Σcos a = Σcos b, sheet Πzw = 1}, the equality of the two cosine sums
    /// sufficing. Construction: sample z₁,z₂,z₃,w₁ free with S_a := Σ(zᵢ+1/zᵢ) ≠ 0 (nonzero sums
    /// make the slice discriminating against the F127 variety, where both sums are 0); the sheet
    /// fixes w₂w₃ = K = 1/(Πz·w₁) and the equality fixes w₂+w₃ = (S_a − w₁ − 1/w₁)·K/(K+1), so
    /// (w₂,w₃) are the roots of x² − ux + K (Tonelli). Both constraints are re-checked exactly
    /// per point (a violation throws). The control re-randomizes w₁ (breaking sheet and equality);
    /// T must be nonzero there. Returns (points, nonZeroOnLocus, controlsNonzero, controlEvals).</summary>
    public static (int Points, int BadOnLocus, int ControlsNonzero, int ControlEvals)
        CertifySharperLocusSlice(long p, int samples)
    {
        long iRoot = SqrtMinusOne(p);
        if (iRoot < 0 || iRoot * iRoot % p != p - 1)
            throw new InvalidOperationException($"p = {p} is not split (need p = 1 mod 4)");
        var stream = new Stream((ulong)p * 0xE7037ED1A0B428DBUL);   // own stream
        long Mul(long a, long b) => a * b % p;
        long Add(long a, long b) => (a + b) % p;
        long Sub(long a, long b) => ((a - b) % p + p) % p;
        long inv2 = Inv(2, p);

        int points = 0, bad = 0, ctrlNonzero = 0, ctrlEval = 0;
        while (points < samples)
        {
            long z1 = stream.Next(2, p), z2 = stream.Next(2, p), z3 = stream.Next(2, p);
            long w1 = stream.Next(2, p);
            long sa = Add(Add(Add(z1, Inv(z1, p)), Add(z2, Inv(z2, p))), Add(z3, Inv(z3, p)));
            if (sa == 0) continue;                       // want Σcos a ≠ 0: the discriminating half
            long k = Inv(Mul(Mul(Mul(z1, z2), z3), w1), p);
            if (Add(k, 1) == 0) continue;
            long v = Sub(sa, Add(w1, Inv(w1, p)));
            long u = Mul(Mul(v, k), Inv(Add(k, 1), p));
            long disc = Sub(Mul(u, u), Mul(4, k));
            long d = Tonelli(disc, p);
            if (d < 0) continue;
            long w2 = Mul(Add(u, d), inv2);
            if (w2 == 0) continue;
            long w3 = Mul(k, Inv(w2, p));
            if (w3 == 0) continue;

            long sb = Add(Add(Add(w1, Inv(w1, p)), Add(w2, Inv(w2, p))), Add(w3, Inv(w3, p)));
            long sheet = Mul(Mul(Mul(Mul(Mul(z1, z2), z3), w1), w2), w3);
            if (sb != sa || sheet != 1)
                throw new InvalidOperationException("constructed point violates a sharper-locus constraint");

            long tVal;
            try { tVal = EvaluateCoreT(p, iRoot, new[] { z1, z2, z3, w1, w2, w3 }); }
            catch (DivideByZeroException) { continue; }
            points++;
            if (tVal != 0) bad++;

            for (int attempt = 0; attempt < 4; attempt++)
            {
                long w1c = stream.Next(2, p);
                long cval;
                try { cval = EvaluateCoreT(p, iRoot, new[] { z1, z2, z3, w1c, w2, w3 }); }
                catch (DivideByZeroException) { continue; }
                ctrlEval++;
                if (cval != 0) ctrlNonzero++;
                break;
            }
        }
        return (points, bad, ctrlNonzero, ctrlEval);
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
