using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using QPoly = System.Collections.Generic.Dictionary<RCPsiSquared.Core.Numerics.PackedMonomial, RCPsiSquared.Core.Numerics.GaussianRational>;
using ZPoly = System.Collections.Generic.Dictionary<RCPsiSquared.Core.Numerics.PackedMonomial, RCPsiSquared.Core.Numerics.GaussianInteger>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The least Schmidt rank of a thermal palindrome element on two sites, the exact half
/// (<c>experiments/THE_CUBE_SITE_BY_SITE.md</c>, finding 1; claim <see cref="Core.Symmetry.MirrorLocalityBoundClaim"/>,
/// its node "the least rank at two sites"). Two sites, hopping (J/2)(XX + YY), optionally (Δ/2)·ZZ, cooling |0⟩⟨1|
/// at γ↓_l and heating |1⟩⟨0| at γ↑_l on each site, c = Σ(γ↓ + γ↑). An element is an invertible superoperator E with
/// E·L·E⁻¹ = −L − c; the Z corner bounds its Schmidt rank across the cut below by the steady state's, 4 here.
///
/// <para>What is decided live, exactly: that no element of Schmidt rank 4 keeps the coherence order
/// d = n_ket − n_bra (without ZZ and with it), and that none reverses it (without ZZ; with ZZ none reverses it at
/// all). The route, all in the computational basis over ℚ(i): the element space as an exact null space (the maps E
/// with E·L + (L + c)·E = 0 that keep d; for reversal the maps U = E·S_XY⁻¹ with L·U = U·L† that keep d, F137's S_XY
/// carrying L onto −L† − c); the realigned matrix of a map that keeps d is block diagonal by the order shift δ of
/// the first site's factor (blocks 1, 4, 6, 4, 1), so rank 4 with the steady state's split 1 + 2 + 1 forces
/// R_{±2} = 0, rank R_{±1} ≤ 1 and rank R_0 ≤ 2: linear equations and 2 × 2 and 3 × 3 minors on the coordinates t;
/// invertibility enters as one Rabinowitsch variable per order sector, z_m·det(E_m) = 1; the two scalings that keep
/// both rank and the equations (the global scale, and ω^d on the order d) are fixed by charts, the first nonzero
/// coordinate of the order-0 and of the order-1 sector set to 1.</para>
///
/// <para>Each chart is decided by <see cref="GroebnerTraceModP"/> over GF(p), p ≡ 1 (mod 4) with i ↦ √−1, and when
/// that run reaches a constant, by <see cref="GroebnerExactReplay"/> over ℚ(i): the steps that led to the constant,
/// repeated on the exact equations, give a nonzero constant as an ideal combination of them, so 1 ∈ I over
/// ℚ(i) and the chart has no point over ℂ. That is a proof; a prime whose replay mismatches is set aside and the
/// next one tried. A chart whose modular run closes WITHOUT a constant is reported as such and is only evidence of
/// a point: that is the anti-vacuity reading (the sectors d ≥ 0 alone are not empty), not a proof of it.</para>
///
/// <para>Also live: the dimensions of the element spaces (exact null spaces), the steady state's split, the
/// reason with ZZ no element reverses d (the characteristic polynomial of L on the order-1 sector is coprime to its
/// conjugate, a resultant over ℤ[i], so the Sylvester equation of the maps U has only the zero solution there), and
/// one control on the equations at one temperature: an exact product element satisfies every one of them with every
/// sector invertible. The full system is never decided on a case with points (stage A of the survivor chart has
/// points by the second route's exact basis; here its modular runs close without a constant at two primes, which is
/// evidence, not a proof).</para>
///
/// <para>Second route: <c>simulations/kms_element_least_rank.py</c> (sympy's Gröbner over ℚ(i), the same
/// coordinates up to one common factor), whose counts (32 charts, 320, 320 and 352 equations, the one stage-A survivor t1 = t3 = 1) the
/// tests pin against this witness. The rank-5 point and the reading of block D stay there, numerical.</para></summary>
public sealed class ThermalElementLeastRankWitness : IInspectable
{
    // ----- the two-site index tables (vec index k = i + 4j for |i⟩⟨j|, site 0 the high bit) -----

    private static readonly (int I, int J)[] Idx = Enumerable.Range(0, 16).Select(k => (k % 4, k / 4)).ToArray();
    private static int Pop(int x) => BitOperations.PopCount((uint)x);

    /// <summary>The coherence order d = n_ket − n_bra of each vec index (THE_ONE_SQUARE §10).</summary>
    public static readonly int[] Charge = Idx.Select(ij => Pop(ij.I) - Pop(ij.J)).ToArray();

    private static readonly (int A, int B)[] Unknowns =
        (from a in Enumerable.Range(0, 16) from b in Enumerable.Range(0, 16) where Charge[a] == Charge[b] select (a, b)).ToArray();

    private static int Bit(int x, int s) => (x >> (1 - s)) & 1;
    private static (int, int) Unit(int k, int s) => (Bit(Idx[k].I, s), Bit(Idx[k].J, s));
    private static readonly (int, int)[] Units = { (0, 0), (0, 1), (1, 0), (1, 1) };
    private static readonly ((int, int) Out, (int, int) In)[] Ops = (from o in Units from i in Units select (o, i)).ToArray();
    private static int Loc((int, int) u) => u.Item1 - u.Item2;
    private static int Shift(((int, int) Out, (int, int) In) op) => Loc(op.Out) - Loc(op.In);
    private static int[] Sector(int m) => Enumerable.Range(0, 16).Where(k => Charge[k] == m).ToArray();

    // ----- the decision for one case -----

    /// <summary>One decided case: the frame (E keeping d, or U keeping d, i.e. elements reversing it), whether its
    /// premises hold on its generator (Z⊗Z an eigenvector at −c and a unique steady state, so E(Z⊗Z) ∝ ρ_ss, split
    /// (1, 2, 1), which forces the rank-4 pattern, and for the frame U the reflection by S_XY), the equations, and per chart whether it was proven empty (a replayed unit
    /// over ℚ(i)). Where a premise fails the case is not decided at all.</summary>
    public sealed record Case(
        string Name, bool PremisesHold, int Dimension, int EquationCount, int ChartCount, int ChartsSkipped, int ChartsProvenEmpty,
        IReadOnlyList<string> StageASurvivors, IReadOnlyList<long> PrimesSetAside, int ReplayedPolynomials, long MaxCoefficientBits)
    {
        /// <summary>The premises hold and every chart is proven empty: no element of Schmidt rank 4 in this frame.</summary>
        public bool RankFourExcluded => PremisesHold && ChartCount > 0 && ChartsProvenEmpty == ChartCount;
    }

    public IReadOnlyList<long> Down { get; }
    public IReadOnlyList<long> Up { get; }
    public long J { get; }

    /// <summary>The ZZ strength Δ of the case "with ZZ".</summary>
    public long Delta { get; }

    public long C { get; }

    /// <summary>All elements (no condition on d): dimension without ZZ, with ZZ.</summary>
    public (int WithoutZz, int WithZz) FullElementDimension { get; }

    /// <summary>The maps U keeping d (the elements reversing d): dimension without ZZ, with ZZ.</summary>
    public (int WithoutZz, int WithZz) ReversingDimension { get; }

    /// <summary>The steady state's ranks on the blocks of the first site's order −1, 0, 1, without ZZ; null where
    /// the steady state is not unique.</summary>
    public IReadOnlyList<int>? SteadyStateSplit { get; }

    /// <summary>The same with ZZ at Δ.</summary>
    public IReadOnlyList<int>? SteadyStateSplitWithZz { get; }

    /// <summary>Whether Z⊗Z is an eigenvector of L at −c (the Z corner): (without ZZ, with ZZ), exactly.</summary>
    public (bool WithoutZz, bool WithZz) ZCorner { get; }

    /// <summary>Whether S_XY, ρ ↦ X⊗X·ρ·Y⊗Y, carries L onto −L† − c (F137): (without ZZ, with ZZ), checked exactly as
    /// S·L + (L† + c)·S = 0.</summary>
    public (bool WithoutZz, bool WithZz) SxyReflects { get; }

    /// <summary>Whether the characteristic polynomial of L on the order-1 sector is coprime to its conjugate:
    /// (without ZZ, with ZZ). Coprime means L there shares no eigenvalue with its conjugate.</summary>
    public (bool WithoutZz, bool WithZz) OrderOneCoprimeToConjugate { get; }

    /// <summary>Whether some map U keeping d is nonzero on the orders ±1: (without ZZ, with ZZ).</summary>
    public (bool WithoutZz, bool WithZz) ReversingTouchesOrderOne { get; }

    /// <summary>The control at one temperature: the product element is one, is invertible, satisfies every
    /// rank-4 equation with every sector determinant nonzero.</summary>
    public bool ControlHolds { get; }

    public int ControlEquationCount { get; }

    /// <summary>Elements keeping d without ZZ; elements keeping d with ZZ; maps U keeping d without ZZ.</summary>
    public IReadOnlyList<Case> Cases { get; }

    /// <summary>The primes the charts are decided at, p ≡ 1 (mod 4) below 2^31, found by trial division.</summary>
    public static readonly long[] Primes = GaussianModP.PrimesOneModFour(6);

    public ThermalElementLeastRankWitness(IReadOnlyList<long>? down = null, IReadOnlyList<long>? up = null, long j = 1, long delta = 1)
    {
        Down = (down ?? new long[] { 3, 2 }).ToArray();
        Up = (up ?? new long[] { 1, 1 }).ToArray();
        if (Down.Count != 2 || Up.Count != 2) throw new ArgumentOutOfRangeException(nameof(down), "two sites: one cooling and one heating rate each");
        J = j; Delta = delta;
        C = Down.Sum() + Up.Sum();
        var l0 = Generator(Down, Up, j, 0);
        var l1 = Generator(Down, Up, j, delta);
        FullElementDimension = (FullDimension(l0, C), FullDimension(l1, C));
        var keep0 = ElementSpace(l0, C, 'E');
        var keep1 = ElementSpace(l1, C, 'E');
        var rev0 = ElementSpace(l0, C, 'U');
        var rev1 = ElementSpace(l1, C, 'U');
        ReversingDimension = (rev0.Count, rev1.Count);
        SteadyStateSplit = SplitOfSteadyState(l0);
        SteadyStateSplitWithZz = SplitOfSteadyState(l1);
        SxyReflects = (Reflects(l0, C), Reflects(l1, C));
        ZCorner = (IsZCorner(l0, C), IsZCorner(l1, C));
        bool forced0 = ZCorner.WithoutZz && SteadyStateSplit is { } s0 && s0.SequenceEqual(new[] { 1, 2, 1 });
        bool forced1 = ZCorner.WithZz && SteadyStateSplitWithZz is { } s1 && s1.SequenceEqual(new[] { 1, 2, 1 });
        OrderOneCoprimeToConjugate = (OrderOneCoprime(l0), OrderOneCoprime(l1));
        ReversingTouchesOrderOne = (TouchesOrderOne(rev0), TouchesOrderOne(rev1));

        (ControlHolds, ControlEquationCount) = Control(j);

        Cases = new[]
        {
            Decide("elements keeping d, without ZZ", keep0, forced0),
            Decide($"elements keeping d, with ZZ at Δ = {delta}", keep1, forced1),
            Decide("maps U keeping d, i.e. elements U·S_XY reversing d, without ZZ", rev0, forced0 && SxyReflects.WithoutZz),
        };
    }

    // ----- the generator, exact over ℚ(i) -----

    private static GaussianRational G(long re, long im = 0) => new(new BigRational(re), new BigRational(im));

    private static GaussianRational[,] Mat4(Func<int, int, GaussianRational> f)
    {
        var m = new GaussianRational[4, 4];
        for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) m[a, b] = f(a, b);
        return m;
    }

    private static GaussianRational[,] MatMul(GaussianRational[,] x, GaussianRational[,] y) =>
        Mat4((a, b) => { var s = GaussianRational.Zero; for (int k = 0; k < 4; k++) s += x[a, k] * y[k, b]; return s; });

    private static GaussianRational[,] Dagger(GaussianRational[,] x) => Mat4((a, b) => x[b, a].Conjugate);

    /// <summary>A one-site 2 × 2 matrix on site l of the two (site 0 the high bit).</summary>
    private static GaussianRational[,] On(int l, GaussianRational[,] a2) =>
        Mat4((r, c) =>
        {
            int r0 = r >> 1, r1 = r & 1, c0 = c >> 1, c1 = c & 1;
            if (l == 0) return r1 == c1 ? a2[r0, c0] : GaussianRational.Zero;
            return r0 == c0 ? a2[r1, c1] : GaussianRational.Zero;
        });

    private static GaussianRational[,] Two(long a, long b, long c, long d) => new[,] { { G(a), G(b) }, { G(c), G(d) } };

    /// <summary>L on vec(ρ) (column-major, index i + 4j for |i⟩⟨j|), column k + 4l the image of |k⟩⟨l|.</summary>
    public static GaussianRational[,] Generator(IReadOnlyList<long> down, IReadOnlyList<long> up, long j, long delta)
    {
        var half = new GaussianRational(new BigRational(1, 2), BigRational.Zero);
        var x = Two(0, 1, 1, 0);
        var y = new[,] { { G(0), G(0, -1) }, { G(0, 1), G(0) } };
        var z = Two(1, 0, 0, -1);
        var xx = MatMul(On(0, x), On(1, x)); var yy = MatMul(On(0, y), On(1, y)); var zz = MatMul(On(0, z), On(1, z));
        var h = Mat4((a, b) => half * (G(j) * (xx[a, b] + yy[a, b]) + G(delta) * zz[a, b]));
        var jumps = new List<(GaussianRational[,] A, long Rate)>();
        for (int l = 0; l < 2; l++)
        {
            jumps.Add((On(l, Two(0, 1, 0, 0)), down[l]));   // |0⟩⟨1|, cooling
            jumps.Add((On(l, Two(0, 0, 1, 0)), up[l]));     // |1⟩⟨0|, heating
        }
        var L = new GaussianRational[16, 16];
        for (int col = 0; col < 16; col++)
        {
            var rho = Mat4((a, b) => a == Idx[col].I && b == Idx[col].J ? GaussianRational.One : GaussianRational.Zero);
            var hr = MatMul(h, rho); var rh = MatMul(rho, h);
            var outM = Mat4((a, b) => GaussianRational.I * (rh[a, b] - hr[a, b]));   // −i[H, ρ]
            foreach (var (A, rate) in jumps)
            {
                if (rate == 0) continue;
                var ad = Dagger(A);
                var ada = MatMul(ad, A);
                var sandwich = MatMul(MatMul(A, rho), ad);
                var anti = MatMul(ada, rho); var anti2 = MatMul(rho, ada);
                var r = G(rate);
                outM = Mat4((a, b) => outM[a, b] + r * (sandwich[a, b] - half * (anti[a, b] + anti2[a, b])));
            }
            for (int row = 0; row < 16; row++) L[row, col] = outM[Idx[row].I, Idx[row].J];
        }
        return L;
    }

    // ----- exact linear algebra over ℚ(i) -----

    /// <summary>Reduced row echelon form in place; returns the pivot column of each pivot row.</summary>
    private static List<int> Rref(List<GaussianRational[]> rows, int cols)
    {
        var pivots = new List<int>();
        int r = 0;
        for (int c = 0; c < cols && r < rows.Count; c++)
        {
            int pr = -1;
            for (int i = r; i < rows.Count; i++) if (!rows[i][c].IsZero) { pr = i; break; }
            if (pr < 0) continue;
            (rows[r], rows[pr]) = (rows[pr], rows[r]);
            var inv = GaussianRational.One / rows[r][c];
            for (int k = c; k < cols; k++) if (!rows[r][k].IsZero) rows[r][k] = rows[r][k] * inv;
            for (int i = 0; i < rows.Count; i++)
            {
                if (i == r || rows[i][c].IsZero) continue;
                var f = rows[i][c];
                for (int k = c; k < cols; k++) if (!rows[r][k].IsZero) rows[i][k] = rows[i][k] - f * rows[r][k];
            }
            pivots.Add(c);
            r++;
        }
        return pivots;
    }

    /// <summary>The canonical null-space basis (one vector per free column, in column order, the free entry 1): the
    /// second route's basis up to one common factor (sympy's rref_den keeps a denominator), which the global-scale
    /// chart absorbs, so the charts, the equations and the survivors are the second route's.</summary>
    private static List<GaussianRational[]> Nullspace(List<GaussianRational[]> rows, int cols)
    {
        var pivots = Rref(rows, cols);
        var free = Enumerable.Range(0, cols).Where(c => !pivots.Contains(c)).ToList();
        var basis = new List<GaussianRational[]>();
        foreach (int f in free)
        {
            var v = Enumerable.Repeat(GaussianRational.Zero, cols).ToArray();
            v[f] = GaussianRational.One;
            for (int i = 0; i < pivots.Count; i++) v[pivots[i]] = -rows[i][f];
            basis.Add(v);
        }
        return basis;
    }

    /// <summary>The maps keeping d as coefficient vectors on the unknowns: frame 'E', E·L + (L + c)·E = 0; frame 'U',
    /// L·U − U·L† = 0.</summary>
    private static List<GaussianRational[]> ElementSpace(GaussianRational[,] L, long c, char frame)
    {
        var pos = new Dictionary<(int, int), int>();
        for (int q = 0; q < Unknowns.Length; q++) pos[Unknowns[q]] = q;
        var rows = new List<GaussianRational[]>();
        for (int a = 0; a < 16; a++)
            for (int b = 0; b < 16; b++)
            {
                var row = Enumerable.Repeat(GaussianRational.Zero, Unknowns.Length).ToArray();
                bool any = false;
                for (int k = 0; k < 16; k++)
                {
                    if (frame == 'E')
                    {
                        if (pos.TryGetValue((a, k), out int p1) && !L[k, b].IsZero) { row[p1] += L[k, b]; any = true; }
                        if (pos.TryGetValue((k, b), out int p2) && !L[a, k].IsZero) { row[p2] += L[a, k]; any = true; }
                    }
                    else
                    {
                        if (pos.TryGetValue((k, b), out int p2) && !L[a, k].IsZero) { row[p2] += L[a, k]; any = true; }
                        var ld = L[b, k].Conjugate;                    // L†[k, b]
                        if (pos.TryGetValue((a, k), out int p1) && !ld.IsZero) { row[p1] -= ld; any = true; }
                    }
                }
                if (frame == 'E' && pos.TryGetValue((a, b), out int p0)) { row[p0] += G(c); any = true; }
                if (any && row.Any(v => !v.IsZero)) rows.Add(row);
            }
        return Nullspace(rows, Unknowns.Length);
    }

    /// <summary>dim of all solutions of E·L + (L + c)·E = 0, no condition on d: L keeps d, so the equations split
    /// into one Sylvester system per pair of orders (m', m), the block of E from order m to order m'.</summary>
    private static int FullDimension(GaussianRational[,] L, long c)
    {
        int total = 0;
        for (int mo = -2; mo <= 2; mo++)
            for (int mi = -2; mi <= 2; mi++)
            {
                var so = Sector(mo); var si = Sector(mi);
                var cols = (from a in so from b in si select (a, b)).ToArray();
                var pos = cols.Select((ab, q) => (ab, q)).ToDictionary(x => x.ab, x => x.q);
                var rows = new List<GaussianRational[]>();
                foreach (var (a, b) in cols)
                {
                    var row = Enumerable.Repeat(GaussianRational.Zero, cols.Length).ToArray();
                    foreach (int k in si) if (!L[k, b].IsZero) row[pos[(a, k)]] += L[k, b];
                    foreach (int k in so) if (!L[a, k].IsZero) row[pos[(k, b)]] += L[a, k];
                    row[pos[(a, b)]] += G(c);
                    rows.Add(row);
                }
                total += cols.Length - Rref(rows, cols.Length).Count;
            }
        return total;
    }

    private static List<int>? SplitOfSteadyState(GaussianRational[,] L)
    {
        var rows = Enumerable.Range(0, 16).Select(r => Enumerable.Range(0, 16).Select(k => L[r, k]).ToArray()).ToList();
        var ns = Nullspace(rows, 16);
        if (ns.Count != 1) return null;                         // not unique: the Z corner's bound does not apply
        var v = ns[0];
        GaussianRational Rho(int i, int jj) => v[i + 4 * jj];
        var split = new List<int>();
        foreach (int d in new[] { -1, 0, 1 })
        {
            var r1 = (from i1 in new[] { 0, 1 } from j1 in new[] { 0, 1 } where i1 - j1 == d select (i1, j1)).ToArray();
            var c2 = (from i2 in new[] { 0, 1 } from j2 in new[] { 0, 1 } where i2 - j2 == -d select (i2, j2)).ToArray();
            var m = r1.Select(x => c2.Select(y => Rho(2 * x.i1 + y.i2, 2 * x.j1 + y.j2)).ToArray()).ToList();
            split.Add(Rref(m, c2.Length).Count);
        }
        return split;
    }

    /// <summary>L·vec(Z⊗Z) = −c·vec(Z⊗Z), entry for entry.</summary>
    private static bool IsZCorner(GaussianRational[,] L, long c)
    {
        var z = Enumerable.Range(0, 16).Select(k => Idx[k].I != Idx[k].J ? 0 : (Pop(Idx[k].I) % 2 == 0 ? 1 : -1)).ToArray();
        for (int a = 0; a < 16; a++)
        {
            var v = G(c * z[a]);
            for (int k = 0; k < 16; k++) if (z[k] != 0) v += L[a, k] * G(z[k]);
            if (!v.IsZero) return false;
        }
        return true;
    }

    /// <summary>S·L + (L† + c)·S = 0 for S the superoperator of ρ ↦ X⊗X·ρ·Y⊗Y (S² = 1, so S·L·S⁻¹ = −L† − c).</summary>
    private static bool Reflects(GaussianRational[,] L, long c)
    {
        var x = Two(0, 1, 1, 0);
        var y = new[,] { { G(0), G(0, -1) }, { G(0, 1), G(0) } };
        var xx = MatMul(On(0, x), On(1, x)); var yy = MatMul(On(0, y), On(1, y));
        var S = new GaussianRational[16, 16];
        for (int col = 0; col < 16; col++)
        {
            var rho = Mat4((a, b) => a == Idx[col].I && b == Idx[col].J ? GaussianRational.One : GaussianRational.Zero);
            var img = MatMul(MatMul(xx, rho), yy);
            for (int row = 0; row < 16; row++) S[row, col] = img[Idx[row].I, Idx[row].J];
        }
        for (int a = 0; a < 16; a++)
            for (int b = 0; b < 16; b++)
            {
                var v = G(c) * S[a, b];
                for (int k = 0; k < 16; k++) v += S[a, k] * L[k, b] + L[k, a].Conjugate * S[k, b];
                if (!v.IsZero) return false;
            }
        return true;
    }

    private static bool OrderOneCoprime(GaussianRational[,] L)
    {
        // 2L has entries in ℤ[i] (the only halves are those of (Δ/2)·ZZ and of the anticommutator); scaling does not move coprimality
        var s1 = Sector(1);
        var m = new GaussianInteger[s1.Length, s1.Length];
        var mc = new GaussianInteger[s1.Length, s1.Length];
        for (int a = 0; a < s1.Length; a++)
            for (int b = 0; b < s1.Length; b++)
            {
                var e = G(2) * L[s1[a], s1[b]];
                if (!e.Re.IsInteger || !e.Im.IsInteger) throw new InvalidOperationException("2L is not integral on the order-1 sector");
                m[a, b] = new GaussianInteger(e.Re.Numerator, e.Im.Numerator);
                mc[a, b] = new GaussianInteger(e.Re.Numerator, -e.Im.Numerator);
            }
        return GaussianPolynomial.AreCoprime(GaussianMatrixCharpoly.Characteristic(m), GaussianMatrixCharpoly.Characteristic(mc));
    }

    private static bool TouchesOrderOne(List<GaussianRational[]> basis) =>
        basis.Any(v => Enumerable.Range(0, Unknowns.Length).Any(q => !v[q].IsZero && Math.Abs(Charge[Unknowns[q].A]) == 1));

    // ----- polynomials over ℚ(i) on the coordinates -----

    private static QPoly Lin(GaussianRational[] form)
    {
        var p = new QPoly();
        for (int q = 0; q < form.Length; q++) if (!form[q].IsZero) p[PackedMonomial.Variable(q)] = form[q];
        return p;
    }

    private static void AddTo(QPoly acc, QPoly b, GaussianRational s)
    {
        foreach (var (m, c) in b)
        {
            var v = (acc.TryGetValue(m, out var old) ? old : GaussianRational.Zero) + s * c;
            if (v.IsZero) acc.Remove(m); else acc[m] = v;
        }
    }

    private static QPoly Mul(QPoly a, QPoly b)
    {
        var r = new QPoly();
        foreach (var (ma, ca) in a)
            foreach (var (mb, cb) in b)
            {
                var m = ma * mb;
                var v = (r.TryGetValue(m, out var old) ? old : GaussianRational.Zero) + ca * cb;
                if (v.IsZero) r.Remove(m); else r[m] = v;
            }
        return r;
    }

    private static QPoly Det(IReadOnlyList<IReadOnlyList<QPoly>> m)
    {
        if (m.Count == 1) return m[0][0];
        var total = new QPoly();
        for (int jj = 0; jj < m.Count; jj++)
        {
            if (m[0][jj].Count == 0) continue;
            var minor = m.Skip(1).Select(row => (IReadOnlyList<QPoly>)row.Where((_, k) => k != jj).ToList()).ToList();
            AddTo(total, Mul(m[0][jj], Det(minor)), jj % 2 == 0 ? GaussianRational.One : -GaussianRational.One);
        }
        return total;
    }

    private static string Key(QPoly f) =>
        string.Join(";", f.OrderBy(kv => kv.Key, PackedMonomial.Descending).Select(kv => $"{kv.Key}:{kv.Value}"));

    private static QPoly Substitute(QPoly f, IReadOnlyDictionary<int, int> values)
    {
        var r = new QPoly();
        foreach (var (m, c) in f)
        {
            var mm = m;
            bool vanish = false;
            foreach (var (v, val) in values)
            {
                int e = m[v];
                if (e == 0) continue;
                if (val == 0) { vanish = true; break; }
                mm = PackedMonomial.Variable(v, e).DivideInto(mm);   // value 1
            }
            if (vanish) continue;
            var nv = (r.TryGetValue(mm, out var old) ? old : GaussianRational.Zero) + c;
            if (nv.IsZero) r.Remove(mm); else r[mm] = nv;
        }
        return r;
    }

    private static GaussianRational Evaluate(QPoly f, IReadOnlyList<GaussianRational> x)
    {
        var s = GaussianRational.Zero;
        foreach (var (m, c) in f)
        {
            var t = c;
            for (int v = 0; v < x.Count; v++) for (int e = m[v]; e > 0; e--) t *= x[v];
            s += t;
        }
        return s;
    }

    /// <summary>f times the lcm of its denominators: the same zero set, coefficients in ℤ[i].</summary>
    private static ZPoly ToGaussianIntegers(QPoly f)
    {
        BigInteger lcm = BigInteger.One;
        foreach (var c in f.Values) { lcm = BigRational.Lcm(lcm, c.Re.Denominator); lcm = BigRational.Lcm(lcm, c.Im.Denominator); }
        var r = new ZPoly();
        foreach (var (m, c) in f)
            r[m] = new GaussianInteger(c.Re.Numerator * (lcm / c.Re.Denominator), c.Im.Numerator * (lcm / c.Im.Denominator));
        return r;
    }

    /// <summary>f with variable v renamed to position[v].</summary>
    private static QPoly Rename(QPoly f, int[] position)
    {
        var r = new QPoly();
        foreach (var (m, c) in f)
        {
            var mm = PackedMonomial.One;
            for (int v = 0; v < position.Length; v++) if (m[v] > 0) mm *= PackedMonomial.Variable(position[v], m[v]);
            r[mm] = c;
        }
        return r;
    }

    private static HashSet<int> Variables(QPoly f)
    {
        var s = new HashSet<int>();
        foreach (var m in f.Keys) for (int v = 0; v < PackedMonomial.MaxVariables; v++) if (m[v] > 0) s.Add(v);
        return s;
    }

    // ----- the rank-4 system -----

    private sealed record RankFourSystem(int N, List<QPoly> Equations, Dictionary<int, List<int>> SectorVariables, Dictionary<int, QPoly> SectorDeterminants);

    private static RankFourSystem BuildSystem(List<GaussianRational[]> basis)
    {
        int n = basis.Count;
        var forms = new Dictionary<(int, int), GaussianRational[]>();
        for (int p = 0; p < Unknowns.Length; p++) forms[Unknowns[p]] = basis.Select(v => v[p]).ToArray();
        var zero = Enumerable.Repeat(GaussianRational.Zero, n).ToArray();
        GaussianRational[] Form(int a, int b) => forms.TryGetValue((a, b), out var f) ? f : zero;

        // the realigned matrix: entry (r, cl) with r the site-0 op and cl the site-1 op of E[a, b]
        var rm = new Dictionary<(((int, int), (int, int)), ((int, int), (int, int))), (int, int)>();
        for (int a = 0; a < 16; a++) for (int b = 0; b < 16; b++) rm[((Unit(a, 0), Unit(b, 0)), (Unit(a, 1), Unit(b, 1)))] = (a, b);
        List<List<QPoly>> Block(int d)
        {
            var rows = Ops.Where(op => Shift(op) == d).ToArray();
            var cols = Ops.Where(op => Shift(op) == -d).ToArray();
            return rows.Select(r => cols.Select(cl => { var (a, b) = rm[(r, cl)]; return Lin(Form(a, b)); }).ToList()).ToList();
        }

        var eqs = new List<QPoly>();
        foreach (int d in new[] { 2, -2 })
            foreach (var row in Block(d)) foreach (var e in row) if (e.Count > 0) eqs.Add(e);
        foreach (var (d, r) in new[] { (1, 2), (-1, 2), (0, 3) })
        {
            var m = Block(d);
            foreach (var ri in Combinations(m.Count, r))
                foreach (var ci in Combinations(m[0].Count, r))
                {
                    var minor = ri.Select(i => (IReadOnlyList<QPoly>)ci.Select(jj => m[i][jj]).ToList()).ToList();
                    var p = Det(minor);
                    if (p.Count > 0) eqs.Add(p);
                }
        }
        var seen = new HashSet<string>();
        eqs = eqs.Where(e => seen.Add(Key(e))).ToList();

        var secvars = new Dictionary<int, List<int>>();
        var dets = new Dictionary<int, QPoly>();
        for (int mm = -2; mm <= 2; mm++)
        {
            var s = Sector(mm);
            secvars[mm] = Enumerable.Range(0, n).Where(q => s.Any(a => s.Any(b => !Form(a, b)[q].IsZero))).ToList();
            dets[mm] = Det(s.Select(a => (IReadOnlyList<QPoly>)s.Select(b => Lin(Form(a, b))).ToList()).ToList());
        }
        return new RankFourSystem(n, eqs, secvars, dets);
    }

    private static IEnumerable<int[]> Combinations(int n, int k)
    {
        var idx = Enumerable.Range(0, k).ToArray();
        while (true)
        {
            yield return idx.ToArray();
            int i = k - 1;
            while (i >= 0 && idx[i] == n - k + i) i--;
            if (i < 0) yield break;
            idx[i]++;
            for (int t = i + 1; t < k; t++) idx[t] = idx[t - 1] + 1;
        }
    }

    private static IEnumerable<Dictionary<int, int>> Charts(List<int> group)
    {
        for (int jj = 0; jj < group.Count; jj++)
        {
            var c = new Dictionary<int, int>();
            for (int i = 0; i < jj; i++) c[group[i]] = 0;
            c[group[jj]] = 1;
            yield return c;
        }
    }

    /// <summary>Decide one input set: true when some prime's trace reaches a constant and its replay over ℚ(i) proves
    /// it (the counts are that replay's); false when two primes' modular runs close without a constant (evidence of a
    /// point, no proof; one prime alone could be a bad one); primes whose replay mismatched are set aside.</summary>
    private static bool DecideIdeal(List<ZPoly> exact, int nvars, List<long> setAside, ref int replayed, ref long maxBits)
    {
        int open = 0;
        foreach (long p in Primes)
        {
            var field = new GaussianModP(p);
            var modular = exact.Select(f => ModPPolynomial.FromTerms(f.Select(kv => (kv.Key, field.Map(kv.Value))), p)).ToList();
            var trace = GroebnerTraceModP.Compute(modular, p, nvars);
            if (!trace.ContainsOneModP)
            {
                if (++open == 2) return false;
                continue;
            }
            var replay = GroebnerExactReplay.ReplayUnit(trace, exact, field);
            if (replay.ProvesUnit)
            {
                replayed += replay.Replayed;
                maxBits = Math.Max(maxBits, replay.MaxBits);
                return true;
            }
            setAside.Add(p);
        }
        throw new InvalidOperationException($"the {Primes.Length} primes neither proved a unit nor closed twice without one");
    }

    private static string ChartLabel(Dictionary<int, int> sub) =>
        string.Join(", ", sub.OrderBy(kv => kv.Key).Select(kv => $"t{kv.Key} = {kv.Value}"));

    private static Case Decide(string name, List<GaussianRational[]> basis, bool premises)
    {
        if (!premises) return new Case(name, false, basis.Count, 0, 0, 0, 0, Array.Empty<string>(), Array.Empty<long>(), 0, 0);
        var sys = BuildSystem(basis);
        int n = sys.N;
        int nvars = n + 5;                                       // t_0..t_{n−1}, then z_m for m = −2..2
        // the variable order of the second route: sympy's generators sorted as strings (t0 < t1 < t10 < … < t9 < z_0 < … < z_m2)
        var names = Enumerable.Range(0, n).Select(q => $"t{q}").Concat(new[] { "z_m2", "z_m1", "z_0", "z_1", "z_2" }).ToArray();
        var order = names.Select((nm, v) => (nm, v)).OrderBy(x => x.nm, StringComparer.Ordinal).Select(x => x.v).ToArray();
        var position = new int[nvars];
        for (int r = 0; r < nvars; r++) position[order[r]] = r;
        ZPoly Exact(QPoly f) => ToGaussianIntegers(Rename(f, position));
        var inA = new HashSet<int>(sys.SectorVariables[0].Concat(sys.SectorVariables[1]).Concat(sys.SectorVariables[2]));
        var eqVars = sys.Equations.Select(Variables).ToList();
        int charts = 0, skipped = 0, proven = 0, replayed = 0;
        long maxBits = 0;
        var survivors = new List<string>();
        var setAside = new List<long>();
        foreach (var c0 in Charts(sys.SectorVariables[0]))
            foreach (var c1 in Charts(sys.SectorVariables[1]))
            {
                var sub = c0.Concat(c1).ToDictionary(kv => kv.Key, kv => kv.Value);
                var ds = Enumerable.Range(-2, 5).ToDictionary(m => m, m => Substitute(sys.SectorDeterminants[m], sub));
                if (ds.Values.Any(d => d.Count == 0)) { skipped++; continue; }
                charts++;
                QPoly Rabinowitsch(int m)
                {
                    var r = Mul(new QPoly { [PackedMonomial.Variable(n + m + 2)] = GaussianRational.One }, ds[m]);
                    AddTo(r, new QPoly { [PackedMonomial.One] = GaussianRational.One }, -GaussianRational.One);
                    return r;
                }
                var stageA = new List<QPoly>();
                var stageB = new List<QPoly>();
                for (int q = 0; q < sys.Equations.Count; q++)
                {
                    var e = Substitute(sys.Equations[q], sub);
                    if (e.Count == 0) continue;
                    (eqVars[q].IsSubsetOf(inA) ? stageA : stageB).Add(e);
                }
                for (int m = 0; m <= 2; m++) stageA.Add(Rabinowitsch(m));
                for (int m = -2; m < 0; m++) stageB.Add(Rabinowitsch(m));
                var exactA = stageA.Select(Exact).ToList();
                bool emptyA = DecideIdeal(exactA, nvars, setAside, ref replayed, ref maxBits);
                if (!emptyA) survivors.Add(ChartLabel(sub));
                bool empty = emptyA || DecideIdeal(exactA.Concat(stageB.Select(Exact)).ToList(), nvars, setAside, ref replayed, ref maxBits);
                if (empty) proven++;
            }
        return new Case(name, true, n, sys.Equations.Count, charts, skipped, proven, survivors, setAside, replayed, maxBits);
    }

    // ----- the control: one temperature, a product element -----

    private (bool Holds, int Equations) Control(long j)
    {
        long[] downB = { 3, 3 }, upB = { 1, 1 };
        long cb = downB.Sum() + upB.Sum();
        var lb = Generator(downB, upB, j, 0);
        var basis = ElementSpace(lb, cb, 'E');
        var sys = BuildSystem(basis);
        // on each site n0 ↦ n1, n1 ↦ −3·n0, p ↦ p, q ↦ −3·q (units (out, in)), the product over the sites
        var site = new Dictionary<((int, int), (int, int)), long>
        {
            [((0, 0), (1, 1))] = -3, [((1, 1), (0, 0))] = 1, [((1, 0), (1, 0))] = 1, [((0, 1), (0, 1))] = -3,
        };
        long Site(((int, int), (int, int)) key) => site.TryGetValue(key, out long v) ? v : 0;
        var eb = new GaussianRational[16, 16];
        for (int a = 0; a < 16; a++) for (int b = 0; b < 16; b++) eb[a, b] = G(Site((Unit(a, 0), Unit(b, 0))) * Site((Unit(a, 1), Unit(b, 1))));
        bool isElement = true;
        for (int a = 0; a < 16 && isElement; a++)
            for (int b = 0; b < 16; b++)
            {
                var s = G(cb) * eb[a, b];
                for (int k = 0; k < 16; k++) s += eb[a, k] * lb[k, b] + lb[a, k] * eb[k, b];
                if (!s.IsZero) { isElement = false; break; }
            }
        var full = Enumerable.Range(0, 16).Select(a => Enumerable.Range(0, 16).Select(b => eb[a, b]).ToArray()).ToList();
        bool invertible = Rref(full, 16).Count == 16;
        // coordinates on the basis: solve Σ_q x_q basis[q] = E_b on the unknowns, exactly
        int nq = basis.Count;
        var aug = Enumerable.Range(0, Unknowns.Length)
            .Select(p => basis.Select(v => v[p]).Append(eb[Unknowns[p].A, Unknowns[p].B]).ToArray()).ToList();
        var piv = Rref(aug, nq + 1);
        if (piv.Contains(nq) || piv.Count != nq) return (false, sys.Equations.Count);
        var x = new GaussianRational[nq + 5];
        for (int i = 0; i < nq; i++) x[piv[i]] = aug[i][nq];
        for (int i = nq; i < nq + 5; i++) x[i] = GaussianRational.Zero;
        bool equations = sys.Equations.All(e => Evaluate(e, x).IsZero);
        bool sectors = sys.SectorDeterminants.Values.All(d => !Evaluate(d, x).IsZero);
        return (isElement && invertible && equations && sectors, sys.Equations.Count);
    }

    // ----- reading -----

    private string Rates => $"γ↓ = ({string.Join(", ", Down)}), γ↑ = ({string.Join(", ", Up)}), J = {J}";

    public string DisplayName => $"the least Schmidt rank of a thermal element on two sites ({Rates}, c = {C})";

    public string Summary =>
        $"rank 4 excluded exactly: {string.Join("; ", Cases.Select(cs => $"{cs.Name} {(cs.RankFourExcluded ? "EXCLUDED" : "NOT DECIDED")} ({cs.ChartsProvenEmpty} of {cs.ChartCount} charts)"))}; " +
        $"with ZZ no element reverses d ({(OrderOneCoprimeToConjugate.WithZz && !ReversingTouchesOrderOne.WithZz && SxyReflects.WithZz ? "exact" : "NOT reproduced")}); " +
        $"control {(ControlHolds ? "holds" : "FAILS")}";

    private static string Split(IReadOnlyList<int>? split) =>
        split is null ? "not unique" : $"unique, split ({string.Join(", ", split)})";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode(
                "1. the element spaces and the premises",
                $"all elements: {FullElementDimension.WithoutZz} without ZZ, {FullElementDimension.WithZz} with ZZ at Δ = {Delta}; " +
                $"keeping d: {Cases[0].Dimension} and {Cases[1].Dimension}; reversing d (maps U): {ReversingDimension.WithoutZz} and {ReversingDimension.WithZz} " +
                $"(exact null spaces over ℚ(i)). Z⊗Z an eigenvector at −c: without ZZ {(ZCorner.WithoutZz ? "yes" : "NO")}, with ZZ " +
                $"{(ZCorner.WithZz ? "yes" : "NO")}. The steady state by the first site's order −1, 0, 1: {Split(SteadyStateSplit)} without ZZ, " +
                $"{Split(SteadyStateSplitWithZz)} with ZZ; a unique one split (1, 2, 1) makes rank 4 force " +
                "R_±2 = 0, rank R_±1 ≤ 1, rank R_0 ≤ 2. S_XY carries L onto −L† − c (exact): without ZZ " +
                $"{(SxyReflects.WithoutZz ? "yes" : "NO")}, with ZZ {(SxyReflects.WithZz ? "yes" : "NO")}.",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                "2. with ZZ no element reverses d",
                $"the characteristic polynomial of L on the order-1 sector is coprime to its conjugate: without ZZ " +
                $"{(OrderOneCoprimeToConjugate.WithoutZz ? "yes" : "no")}, with ZZ {(OrderOneCoprimeToConjugate.WithZz ? "yes" : "no")} " +
                "(a resultant over ℤ[i]); coprime means the maps U keeping d vanish on the orders ±1, and then every U·S_XY is singular " +
                $"(the maps U touch them: without ZZ {(ReversingTouchesOrderOne.WithoutZz ? "yes" : "no")}, with ZZ {(ReversingTouchesOrderOne.WithZz ? "yes" : "no")}). " +
                (OrderOneCoprimeToConjugate.WithZz && !ReversingTouchesOrderOne.WithZz && SxyReflects.WithZz
                    ? "With ZZ no element reverses d." : "With ZZ this is NOT reproduced here."),
                provenance: NodeProvenance.Live);
            int idx = 3;
            foreach (var cs in Cases)
                yield return new InspectableNode(
                    $"{idx++}. {cs.Name}",
                    !cs.PremisesHold
                        ? "NOT DECIDED: a premise fails on this generator (the steady state is not unique or does not split (1, 2, 1), " +
                          "or S_XY does not reflect), so rank 4 does not reduce to the equations here."
                        : $"{cs.Dimension} coordinates, {cs.EquationCount} rank-4 equations, {cs.ChartCount} charts decided" +
                          (cs.ChartsSkipped > 0 ? $" and {cs.ChartsSkipped} more skipped (a sector determinant identically zero, empty without a computation)" : "") +
                          $"; proven empty over ℚ(i): {cs.ChartsProvenEmpty} of {cs.ChartCount} " +
                          $"(modular trace, then the exact replay over ℚ(i) of {cs.ReplayedPolynomials} polynomials, monic, coefficients up to {cs.MaxCoefficientBits} bits" +
                          (cs.PrimesSetAside.Count > 0 ? $"; primes set aside after a mismatched replay: {string.Join(", ", cs.PrimesSetAside)}" : "") + "). " +
                          (cs.RankFourExcluded ? "No element of Schmidt rank 4 here. " : "NOT every chart is proven empty. ") +
                          $"Anti-vacuity, the sectors d ≥ 0 alone: the modular runs close without a constant at two primes in {cs.StageASurvivors.Count} chart(s)" +
                          (cs.StageASurvivors.Count > 0 ? $" ({string.Join("; ", cs.StageASurvivors)}), evidence of a point, not a proof" : "") + ".",
                    provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                $"{idx}. control: one temperature",
                $"(3, 3 | 1, 1): that the product map n0 ↦ n1, n1 ↦ −3·n0, p ↦ p, q ↦ −3·q on each site is an invertible element, " +
                $"satisfies all {ControlEquationCount} rank-4 equations and keeps every sector invertible (exact evaluation): " +
                (ControlHolds ? "holds." : "FAILS.") + " A product element has Schmidt rank 1, so this checks the construction of the " +
                "equations, not the decision.",
                provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
