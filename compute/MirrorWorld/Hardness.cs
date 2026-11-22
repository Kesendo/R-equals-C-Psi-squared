using System.Numerics;

namespace MirrorWorld;

// The hardness of the palindrome (adopted 2026-07-04 from the F87-hardness bloc: F102/F103 with
// F105/F106 stability, F107/F109 purity, F110/F111 cell rules, F115 the GF(2)[x] valuation theory
// of PROOF_F103_F87_Z2_CUBED_REFINEMENT.md section 7.7, F117 the Pascal-Gram converse of
// PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md section 5; the window floor, F173's Theorem B of
// PROOF_SIGNED_MASK_BALANCE.md, adopted 2026-11-22).
//
// F87 classifies a windowed Pauli pair by how the palindrome closes: truly (the operator identity
// holds, M = 0), soft (only the spectrum still pairs about -sigma), hard (even the spectrum
// breaks). The DEFINITION is spectral -- and the bloc's whole achievement is that it never has to
// be solved again. MirrorWorld adopts the EIGENSOLVER-FREE certificates the bloc built:
//
//   the cube    -- y_par = #Y mod 2 is the third Z2 beside the Klein axes (F102), and the purity
//                  rules are pure letter parity: truly forces y_par = 0 under every dephase
//                  letter against its canonical mirror (F107; other mirrors of the same
//                  dissipator put part of truly at y_par = 1, F103 section 8); the mother
//                  sector's non-truly side is all-odd, y_par = 1 (F109);
//                  hard Klein-homogeneous pairs live only in the dephase letter's own diagonal Klein cell, Y-inverted
//                  (F110), at k=N=4 exactly on the pure-D templates (F111).
//   the valuation -- a diagonal-cell pair's X/Y window masks, read as GF(2)[x] polynomials, decide
//                  hardness: hard iff the (1+x)-adic valuations differ (F115), on every chain from the
//                  pair's window floor N* on and soft below it (N* exact by F173's Theorem B, adopted
//                  here as WindowFloor; N >= 2k-2 always suffices), so on chains of N >= 2k-2 sites in ONE SUBTRACTION.
//                  Classes of size 2^(k-1-v), the A203241 hard-pair count there, the obstruction ceiling
//                  min(2W-1, 2k-3). Pure bit arithmetic, checked in the main repo to k = 20.
//   the traces  -- the all-gamma converse (F117) reads hardness from odd power-sums of the
//                  recentered M = A + gamma Q: soft means every odd Tr(M^m) vanishes; hard fires
//                  at m* = the FIRST odd m with Tr(M^m) != 0, a Pascal-Gram sum of squares, so
//                  positive at EVERY gamma. The deg-1 girth face fires at 2*girth+1 only when
//                  its moment t_girth != 0; when it is silent (silence is not softness) a higher
//                  class fires, m* = 2*girth + deg with deg odd -- the K3 pair fires at
//                  m* = 2*3 + 3 = 9. Traces, never eigenvalues -- Newton's identities stand in
//                  for spectra.
//
// The one thing deliberately NOT adopted: the spectral classifier itself (F87's definition, the
// F104 engine). It is the thing the certificates replaced; it stays in the main repo.
public sealed class Hardness : GameObject
{
    public Hardness(World world) : base(world) { }

    // left: what the hardness itself produces.
    public override IReadOnlyList<string> Own => new[] { "valuation", "cube", "traces" };

    // ---- the valuation face (F115): GF(2)[x], the mask's (1+x)-adic valuation ----

    // v_{1+x}(mask): how often (1+x) divides the mask read as a polynomial (bit j = x^j).
    // An odd-popcount mask is not divisible at all (v = 0): bit_a is the lowest valuation bit.
    public static int Valuation(ulong mask)
    {
        int v = 0;
        while (mask != 0 && BitOperations.PopCount(mask) % 2 == 0)
        {
            mask = DivideByOnePlusX(mask);
            v++;
        }
        return v;
    }

    // hard iff the two valuations differ, on a chain from the pair's window floor on (N >= 2k-2 always suffices);
    // the method takes no N and reads that verdict (IsHardPairAt below takes one), on chains of N >= 2k-2 sites
    // the whole section-7 verdict in one subtraction.
    public static bool IsHardPair(ulong p1, ulong p2) => Valuation(p1) != Valuation(p2);

    // the window floor (F173 Theorem B, adopted 2026-11-22): two Mixed diagonal-cell templates of widths w1, w2 with
    // X/Y masks m1, m2 on the open chain. Every relation among their placements is a multiple r*(q2, q1) of the reduced
    // pair, q_i = m_i / gcd(m1, m2), and fits a chain of N sites exactly when N >= max(w1 + deg q2, w2 + deg q1) + deg r.
    // In the diagonal cell both templates carry an odd #(Y+Z), so a relation is odd exactly when r(1) = 1 and
    // q1(1) + q2(1) = 1, that is when the valuations differ, and r = 1 fits first: the pair is soft below the floor
    // and hard from it on. The method sees only the masks, so the odd #(Y+Z) is the caller's promise: XX + XIX,
    // both even, never breaks, though its masks read a floor of 3. Equal valuations give no floor, soft at every
    // length; translates, outside Theorem B's hypothesis, are among them (every relation even, PROOF_F103 section
    // 7.7). (1+x) divides the gcd, so the floor is at most w1 + w2 - 2, which is 2k - 2 for two templates of width k.
    public static int? WindowFloor(ulong m1, int w1, ulong m2, int w2)
    {
        DiagonalCellMask(m1, w1);
        DiagonalCellMask(m2, w2);
        if (!IsHardPair(m1, m2)) return null;
        ulong g = Gcd2(m1, m2), q1 = Div2(m1, g), q2 = Div2(m2, g);
        return Math.Max(w1 + Deg(q2), w2 + Deg(q1));
    }

    // the verdict at one chain length: hard exactly from the window floor on.
    public static bool IsHardPairAt(ulong m1, int w1, ulong m2, int w2, int n) =>
        WindowFloor(m1, w1, m2, w2) is int floor && n >= floor;

    // the floor's scope: a Mixed diagonal-cell X/Y mask is nonzero with an even popcount and lies inside its template
    static void DiagonalCellMask(ulong m, int w)
    {
        if (m == 0 || BitOperations.PopCount(m) % 2 != 0)
            throw new ArgumentException("a Mixed diagonal-cell mask is nonzero with an even popcount");
        if (w < 1 || w > 63)
            throw new ArgumentException("a template's width lies between 1 and 63");
        if (Deg(m) >= w)
            throw new ArgumentException("the mask must lie inside its template's width");
    }

    static int Deg(ulong p) => 63 - BitOperations.LeadingZeroCount(p);

    // the gcd of two polynomials over GF(2), bit i the coefficient of x^i (the stitch's readers use it too)
    internal static ulong Gcd2(ulong a, ulong b)
    {
        while (b != 0)
        {
            int db = Deg(b);
            while (a != 0 && Deg(a) >= db) a ^= b << (Deg(a) - db);
            (a, b) = (b, a);
        }
        return a;
    }

    // the exact quotient a / b over GF(2); b divides a
    static ulong Div2(ulong a, ulong b)
    {
        ulong q = 0;
        int db = Deg(b);
        while (a != 0 && Deg(a) >= db)
        {
            int s = Deg(a) - db;
            q |= 1UL << s;
            a ^= b << s;
        }
        if (a != 0) throw new ArgumentException("the divisor does not divide");
        return q;
    }

    // the cross-class pair count over even-popcount nonzero k-bit masks: (4^(k-1) - 3 2^(k-1) + 2)/3
    // (OEIS A203241), the hard mask pairs on a chain of N >= 2k-2 sites, where every pair has reached its
    // window floor, and the Klein/y-parity dressing 2^(2k-3).
    public static long HardMaskPairCount(int k) => ((1L << (2 * (k - 1))) - 3L * (1L << (k - 1)) + 2) / 3;
    public static long DressedHardCount(int k) => HardMaskPairCount(k) << (2 * k - 3);

    // the obstruction-size law (the shared-factor-free face): the largest minimal odd cycle over
    // hard pairs is min(2W - 1, 2k - 3); at k = 3 every obstruction is the triangle.
    public static int ObstructionCeiling(int k, int w) => Math.Min(2 * w - 1, 2 * k - 3);

    static ulong DivideByOnePlusX(ulong p)
    {
        ulong q = 0;
        while (p > 1)
        {
            int d = 63 - BitOperations.LeadingZeroCount(p);
            q |= 1UL << (d - 1);
            p ^= 0b11UL << (d - 1);
        }
        return q;
    }

    // ---- the cube face (F102 + the purity rules F107/F109 + the cell rules F110/F111) ----

    // (bit_a, bit_b, y_par) = ((#X+#Y)%2, (#Y+#Z)%2, #Y%2): the two Klein axes and the third.
    public static (int A, int B, int Y) Cube(char[] letters)
    {
        int nx = letters.Count(c => c == 'X'), ny = letters.Count(c => c == 'Y'), nz = letters.Count(c => c == 'Z');
        return ((nx + ny) % 2, (ny + nz) % 2, ny % 2);
    }

    // the per-dephase truly criteria (F107): the two letters other than the canonical mirror's flip
    // letter (X for Z and Y, Z for X) occur an even number of times, and all three contain
    // "#Y even" -- truly forces y_par = 0.
    public static bool Truly(char[] letters, char dephase)
    {
        int nx = letters.Count(c => c == 'X'), ny = letters.Count(c => c == 'Y'), nz = letters.Count(c => c == 'Z');
        return dephase switch
        {
            'Z' => ny % 2 == 0 && nz % 2 == 0,
            'X' => nx % 2 == 0 && ny % 2 == 0,
            'Y' => ny % 2 == 0 && nz % 2 == 0,
            _ => throw new ArgumentException($"not a dephase letter: {dephase}"),
        };
    }

    // hard Klein-homogeneous pairs live only in the dephase letter's OWN Klein cell (F110): Z -> (0,1), X -> (1,0),
    // Y -> (1,1) -- the cell of the letter itself.
    public static (int A, int B) DiagonalCell(char dephase)
    {
        var (a, b, _) = Cube(new[] { dephase });
        return (a, b);
    }

    // a pure-D template (F111): only the dephase letter and I -- at k = N = 4 these are exactly
    // the hard carriers (a pair is hard iff it touches one).
    public static bool IsPureTemplate(char[] letters, char dephase)
        => letters.All(c => c == 'I' || c == dephase) && letters.Any(c => c == dephase);

    // the adopted census constants: the k=3 diagonal hard split (N-stable from N=4 up,
    // F103 = F105 bit-exact; at N=3 the same cell reads 34:0, the split without the
    // adjacency half of F103 §6's rule, which needs the term placed at two windows)
    // and the k=N=4 template decomposition 36 + 192 + 300 = 528 with 228:0 purity (F106/F111).
    public static (int Dominant, int Recessive) HardSplitK3 => (42, 8);
    public static (int PurePure, int PureMixed, int MixedMixed) TemplateDecompositionK4 => (36, 192, 300);

    // ---- the trace face (F117): odd power-sums of the recentered M = A + gamma Q ----
    // A = -i[H,.] as a d^2 x d^2 matrix, Q = the dephasing frequency diagonal N - 2 popcount(i^j);
    // spec(M) symmetric about 0 (soft) iff every odd Tr(M^m) vanishes; a hard pair fires at the
    // FIRST odd m with Tr(M^m) != 0 (m* = 2*girth + deg, deg odd; the deg-1 face only when its
    // girth moment is nonzero), positive at every gamma (the Pascal-Gram sum of squares). Traces only.
    public double[] OddPowerSums(char[][] templates, int n, double j, double gamma, int upToOdd)
    {
        int d = 1 << n;
        var h = new Complex[d, d];
        foreach (var template in templates)
        {
            int k = template.Length;
            for (int w = 0; w + k <= n; w++)
            {
                var letters = Enumerable.Repeat('I', n).ToArray();
                for (int s = 0; s < k; s++) letters[w + s] = template[s];
                var term = Dense.PauliString(letters);
                for (int p = 0; p < d; p++)
                    for (int q = 0; q < d; q++)
                        h[p, q] += j * term[p, q];
            }
        }

        int dim = d * d;
        var m = new Complex[dim, dim];
        for (int i = 0; i < d; i++)                                          // A = -i (H x I - I x H^T)
            for (int jdx = 0; jdx < d; jdx++)
            {
                int row = i * d + jdx;
                for (int i2 = 0; i2 < d; i2++)
                    m[row, i2 * d + jdx] += -Complex.ImaginaryOne * h[i, i2];
                for (int j2 = 0; j2 < d; j2++)
                    m[row, i * d + j2] += Complex.ImaginaryOne * h[j2, jdx];
                m[row, row] += gamma * (n - 2 * BitOperations.PopCount((uint)(i ^ jdx)));
            }

        var sums = new List<double>();
        var power = m;
        for (int exponent = 1; exponent <= upToOdd; exponent++)
        {
            if (exponent > 1) power = Dense.Mul(power, m);
            if (exponent % 2 == 1) sums.Add(Dense.Trace(power).Real);
        }
        return sums.ToArray();
    }
}
