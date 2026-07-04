using System.Numerics;

namespace MirrorWorld;

// The golden ceiling router (adopted 2026-07-04 from F116, docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md
// sections 1-3 + 8): the two Z-middle ceiling cases -- the open-chain sliding-window Hamiltonians of
// XZX + XZY + YZX and its X<->Y sibling under Z-dephasing -- are palindromized by a LOCAL, period-4,
// per-site router:
//
//   W = (x)_l q_{l mod 4},   W L W^-1 = -L - 2 sigma,   sigma = sum_l gamma_l  (site rates arbitrary)
//
// The frame walks [a, a, b, b] with a = phi X + Y and b = X - phi Y, the two projective roots of the
// golden locus alpha^2 - alpha beta - beta^2 = 0; h_l = q_l(Z) = (-1)^(l+1) i R(g_l) (R the quarter
// turn in the X-Y plane); every q_l is class-swapping ({I,Z} <-> {X,Y}) with q_l^2 = -(2+phi) I, so
// W is a scalar times a unitary (condition number 1). The mechanism is a WINDOW LEMMA: the
// window-summed anticommutator {q(x)q(x)q, [XZX+XZY+YZX, .]_3} vanishes identically at all four
// offsets (cross-template cancellation inside one window -- per-term it fails, which is why the
// per-term certifier never saw it), and window additivity lifts it to every N >= 3. The golden case
// is c = 1 of the METALLIC family: for t1 XZX + t2 XZY + t3 YZX the soft set is exactly the line
// t2 = t3, and with c = t1/t2 the same [a,a,b,b] router works at a = (r,1), b = (1,-r),
// r(c) = (c + sqrt(c^2+4))/2 the metallic mean (F116_MetallicMean), q_l^2 = -(1+r^2) I.
//
// Two faces verified from below, no eigensolver: the window lemma on the 64-dim window space
// (coefficient algebra), and the two-sided dense end-to-end of the proof's section 2,
// W(rho) = (2+phi)^(N/2) P rho Q with P = (x)(aZ | Z) and Q = (x)(Z | Za) alternating over sites
// (a-hat the normalized frame axis), P H P^-1 = Q H Q^-1 = -H -- the conjugation checked on the
// full Pauli basis with site-dependent rates.
public sealed class Router : GameObject
{
    public Router(World world) : base(world) { }

    // left: what the router itself produces.
    public override IReadOnlyList<string> Own => new[] { "frame", "window", "conjugation" };

    // ---- the four per-site maps (coefficient space, rows/cols = I,X,Y,Z), the [a,a,b,b] frame:
    // g = (r,1),(r,1),(1,-r),(1,-r); h = (-1)^(l+1) i R(g). Golden = r(1) = phi. ----
    public static Complex[][,] SiteMaps(double r)
    {
        var im = Complex.ImaginaryOne;
        return new[]
        {
            BuildSiteMap(r, 1.0, im, -im * r),
            BuildSiteMap(r, 1.0, -im, im * r),
            BuildSiteMap(1.0, -r, -im * r, -im),
            BuildSiteMap(1.0, -r, im * r, im),
        };
    }

    static Complex[,] BuildSiteMap(Complex gx, Complex gy, Complex hx, Complex hy) => new[,]
    {
        { Complex.Zero, -gx, -gy, Complex.Zero },
        { gx, Complex.Zero, Complex.Zero, hx },
        { gy, Complex.Zero, Complex.Zero, hy },
        { Complex.Zero, hx, hy, Complex.Zero },
    };

    // the sibling's maps: conjugation by s (I -> I, X <-> Y, Z -> -Z), per site.
    static Complex[,] MirrorMap(Complex[,] q)
    {
        var s = new Complex[,]
        {
            { 1, 0, 0, 0 },
            { 0, 0, 1, 0 },
            { 0, 1, 0, 0 },
            { 0, 0, 0, -1 },
        };
        return Dense.Mul(s, Dense.Mul(q, s));
    }

    // ---- the window lemma: worst-offset Frobenius norm of {Q_k, S} on the 64-dim window space,
    // S the summed commutator superoperator of the weighted templates. Zero exactly on the
    // metallic locus r = r(c); O(1) off it. ----
    public static double WindowAnticommutatorNorm(double c, double r, bool sibling)
    {
        var maps = SiteMaps(r);
        if (sibling) maps = maps.Select(MirrorMap).ToArray();
        var templates = sibling
            ? new[] { ("YZY", c), ("YZX", 1.0), ("XZY", 1.0) }
            : new[] { ("XZX", c), ("XZY", 1.0), ("YZX", 1.0) };
        var s = SummedCommutator(templates);
        double worst = 0;
        for (int offset = 0; offset < 4; offset++)
        {
            var qk = Dense.Kron(maps[offset % 4], Dense.Kron(maps[(offset + 1) % 4], maps[(offset + 2) % 4]));
            var qs = Dense.Mul(qk, s);
            var sq = Dense.Mul(s, qk);
            double norm = 0;
            for (int i = 0; i < 64; i++)
                for (int j = 0; j < 64; j++)
                {
                    var v = qs[i, j] + sq[i, j];
                    norm += v.Real * v.Real + v.Imaginary * v.Imaginary;
                }
            worst = Math.Max(worst, Math.Sqrt(norm));
        }
        return worst;
    }

    // S = sum_T w_T [T, .]_3 on the 4^3 window coefficient space (site 0 = the most significant
    // base-4 digit, matching the Kronecker order): column s carries 2 w phase at row T*s exactly
    // when T and s anticommute (an odd number of anticommuting sites).
    static Complex[,] SummedCommutator((string Letters, double Weight)[] templates)
    {
        var s = new Complex[64, 64];
        foreach (var (letters, weight) in templates)
            for (int col = 0; col < 64; col++)
            {
                Complex phase = Complex.One;
                int row = 0, anti = 0;
                for (int site = 0; site < 3; site++)
                {
                    int a = Idx(letters[site]);
                    int b = (col >> (2 * (2 - site))) & 3;
                    var (product, ph) = Product(a, b);
                    row |= product << (2 * (2 - site));
                    phase *= ph;
                    if (a != 0 && b != 0 && a != b) anti++;
                }
                if (anti % 2 == 1) s[row, col] += 2.0 * weight * phase;
            }
        return s;
    }

    static int Idx(char letter) => letter switch { 'I' => 0, 'X' => 1, 'Y' => 2, 'Z' => 3, _ => throw new ArgumentException($"{letter}") };

    // the single-site Pauli product a*b = phase * letter.
    static (int Letter, Complex Phase) Product(int a, int b)
    {
        if (a == 0) return (b, Complex.One);
        if (b == 0) return (a, Complex.One);
        if (a == b) return (0, Complex.One);
        int product = 6 - a - b;                                     // the third letter
        bool cyclic = (a, b) is (1, 2) or (2, 3) or (3, 1);          // XY=iZ, YZ=iX, ZX=iY
        return (product, cyclic ? Complex.ImaginaryOne : -Complex.ImaginaryOne);
    }

    // ---- the two-sided dense face (the proof's section 2, previously verified only in Python):
    // P = (x)_l (a-hat Z at even l, Z at odd l), Q = (x)_l (Z at even l, Z a-hat at odd l), with
    // a-hat the normalized [a,a,b,b] frame axis. Returns the worst entries of the chiral pair
    // P H P^-1 + H, Q H Q^-1 + H, and of W L W^-1 + L + 2 sigma over the full Pauli basis. ----
    public (double ChiralP, double ChiralQ, double Conjugation) DenseResiduals(int n, double c, double[] gammas)
    {
        double r = Formulas.F116_MetallicMean(c);
        int d = 1 << n;

        var h = new Complex[d, d];
        foreach (var (letters, weight) in new[] { ("XZX", c), ("XZY", 1.0), ("YZX", 1.0) })
            for (int w = 0; w + 3 <= n; w++)
            {
                var full = Enumerable.Repeat('I', n).ToArray();
                for (int site = 0; site < 3; site++) full[w + site] = letters[site];
                var term = Dense.PauliString(full);
                for (int row = 0; row < d; row++)
                    for (int col = 0; col < d; col++)
                        h[row, col] += weight * term[row, col];
            }

        // per-site factors: the frame axis a-hat, then aZ / Z on the P side, Z / Za on the Q side.
        var x = new Complex[,] { { 0, 1 }, { 1, 0 } };
        var y = new Complex[,] { { 0, -Complex.ImaginaryOne }, { Complex.ImaginaryOne, 0 } };
        var z = new Complex[,] { { 1, 0 }, { 0, -1 } };
        double norm = Math.Sqrt(1.0 + r * r);
        Complex[,] Axis(int l)
        {
            var (ax, ay) = l % 4 <= 1 ? (r, 1.0) : (1.0, -r);
            var m = new Complex[2, 2];
            for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                    m[i, j] = (ax * x[i, j] + ay * y[i, j]) / norm;
            return m;
        }

        Complex[,] p = null!, pInv = null!, q2 = null!, qInv = null!;
        for (int l = 0; l < n; l++)
        {
            var axis = Axis(l);
            var (pl, plInv) = l % 2 == 0 ? (Dense.Mul(axis, z), Dense.Mul(z, axis)) : (z, z);
            var (ql, qlInv) = l % 2 == 0 ? (z, z) : (Dense.Mul(z, axis), Dense.Mul(axis, z));
            p = l == 0 ? pl : Dense.Kron(p, pl);
            pInv = l == 0 ? plInv : Dense.Kron(pInv, plInv);
            q2 = l == 0 ? ql : Dense.Kron(q2, ql);
            qInv = l == 0 ? qlInv : Dense.Kron(qInv, qlInv);
        }

        double chiralP = Dense.MaxAbsSum(Dense.Mul(p, Dense.Mul(h, pInv)), h);
        double chiralQ = Dense.MaxAbsSum(Dense.Mul(q2, Dense.Mul(h, qInv)), h);

        var zs = Enumerable.Range(0, n).Select(l => Dense.SiteZ(n, l)).ToArray();
        double sigma = gammas.Sum();
        double conjugation = 0;
        foreach (var letters in Dense.AllStrings(n))
        {
            var rho = Dense.PauliString(letters);
            var inner = Dense.Mul(pInv, Dense.Mul(rho, qInv));
            var lInner = Add(Dense.CommutatorFlow(h, inner), Dense.Dephasing(zs, gammas, inner));
            var lhs = Dense.Mul(p, Dense.Mul(lInner, q2));
            var lRho = Add(Dense.CommutatorFlow(h, rho), Dense.Dephasing(zs, gammas, rho));
            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                    conjugation = Math.Max(conjugation,
                        (lhs[i, j] + lRho[i, j] + 2.0 * sigma * rho[i, j]).Magnitude);
        }
        return (chiralP, chiralQ, conjugation);
    }

    static Complex[,] Add(Complex[,] a, Complex[,] b)
    {
        int d = a.GetLength(0);
        var m = new Complex[d, d];
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                m[i, j] = a[i, j] + b[i, j];
        return m;
    }
}
