using System.Globalization;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The block from one excitation to N − 1, live (claim <see cref="OneDephasedSeatClaim"/>, the bullet of that
/// name in the paragraph "The sheet and blindness, at N = 3 and N = 4" of the section on blind dephased seats of
/// <c>docs/proofs/PROOF_ONE_DEPHASED_SEAT.md</c>). Heisenberg bonds on a connected graph of N sites (signed couplings),
/// the jump Z on seat d = 0 alone, fields perpendicular to z on the other sites and none on d.
///
/// <para><b>The lemma.</b> A far element W ([H, W] = 0, {Z_d, W} = 0) homogeneous under Z^⊗N commutes with the bond
/// part and the field part separately. Its block 𝒯 from the one-excitation sector to the (N − 1)-excitation sector then
/// solves B_(N−1)·𝒯 = 𝒯·B_1 and Z_d·𝒯 = −𝒯·Z_d, and Ξ = X^⊗N carries the solutions one to one onto the operators
/// on the one-excitation sector that commute with the signed Laplacian 𝓛_J and the projector P_d onto e_d, whose
/// dimension is F157's Theorem A count, 1 + dim commutant(𝓛_J on the Krylov complement of e_d). So the solution space
/// is one-dimensional, spanned by Ξ's block, exactly where d is not blind.</para>
///
/// <para><b>A road of its own.</b> The proof's gate, <c>simulations/blind_seat_block_lemma_gate.py</c>, restricts the
/// bond part from the full 2^N matrix of Pauli products and solves in sympy, and counts far elements by a rank modulo a
/// prime. This witness builds the two sectors' bond blocks from the bit action (ZZ a sign per bond, XX + YY the flip of
/// an unequal pair with weight 2J), solves for the block over ℚ, and counts the other side of the bijection, the
/// commutant of 𝓛_J and P_d, separately from the couplings, without the identity B_1 = constant − 2𝓛_J; the blind
/// dimension is N minus the Krylov rank of 𝓛_J. The far elements homogeneous under Z^⊗N are counted exactly over ℚ
/// as the kernel of W ↦ [H, W] on the Pauli strings with X or Y at d and a fixed parity of height (the number of letters
/// X or Y).</para>
///
/// <para><b>What is read.</b> At every row: the block's solution dimension against the commutant count, against
/// blindness, and Ξ's block as a solution. At N = 3 with d not blind, the proof's two conclusions: no homogeneous far
/// element of even height, and one of odd height exactly when the fields lie on one line through the origin (class A's
/// colouring then supplies it). At N = 4 with d not blind, the bullet "At N = 4, a seat that is not blind": none of
/// either height with the fields off one line, and with them on one line, not all zero, exactly one, of even height
/// (class A's colouring), but two of even height on the antiparallel locus (a field zero, the other two opposite, and the
/// field-free site coupled to the other two exactly as site 0 is, oppositely to each). Beside it the descent's two conditions at the row, computed exactly: (I) the dimension
/// of the cyclic space alg(H_b, Z_d, F⁺)|vac⟩ over ℚ(i), realified, and (II) per sector the vectors with d excited in
/// every term at the vacuum's energy; at N = 5 the counts and conditions are reported, (II) predicted only at the
/// sectors 1, 2 and N − 1. Fixed grids at N = 3, 4 and 5 for the block and at N = 3 and 4 for the far counts.</para>
///
/// <para>Args: <c>--N</c> (3 to 5; default 3), <c>--bonds</c> (i-j:J entries, comma-separated, J an integer or p/q;
/// default 0-1:1,0-2:2,1-2:3), <c>--fields</c> (j:x,y entries separated by semicolons, x and y integers or p/q, no
/// entry for site 0; default 1:3,4;2:5,0, quoted in a shell for the semicolons).</para></summary>
public sealed class DephasedSeatFlipBlockWitness : IInspectable
{
    public const int MinN = 3, MaxN = 5;

    public int N { get; }
    public IReadOnlyDictionary<(int, int), BigRational> Bonds { get; }
    public IReadOnlyDictionary<int, (BigRational X, BigRational Y)> Fields { get; }

    public DephasedSeatFlipBlockWitness(int n = 3, string? bonds = null, string? fields = null)
    {
        if (n < MinN || n > MaxN)
            throw new ArgumentOutOfRangeException(nameof(n), n, $"--N runs from {MinN} to {MaxN}.");
        N = n;
        Bonds = ParseBonds(bonds ?? "0-1:1,0-2:2,1-2:3", n);
        if (!Connected(n, Bonds))
            throw new ArgumentException("--bonds leave the graph disconnected; the lemma is stated on a connected graph.");
        Fields = ParseFields(fields ?? (n == 3 ? "1:3,4;2:5,0" : ""), n);
    }

    // ---- parsing ----------------------------------------------------------------------------------------

    private static BigRational Rat(string s, string arg)
    {
        var q = s.Trim().Split('/');
        if (q.Length > 2 || !long.TryParse(q[0], NumberStyles.Integer, CultureInfo.InvariantCulture, out var p)
            || (q.Length == 2 && !long.TryParse(q[1], NumberStyles.Integer, CultureInfo.InvariantCulture, out _)))
            throw new ArgumentException($"--{arg}: cannot read \"{s}\" as an integer or p/q.");
        long den = q.Length == 2 ? long.Parse(q[1], CultureInfo.InvariantCulture) : 1;
        if (den == 0) throw new ArgumentException($"--{arg}: zero denominator in \"{s}\".");
        return new BigRational(p, den);
    }

    public static Dictionary<(int, int), BigRational> ParseBonds(string spec, int n)
    {
        var d = new Dictionary<(int, int), BigRational>();
        foreach (var entry in spec.Split(',', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries))
        {
            var parts = entry.Split(':');
            var ends = parts[0].Split('-');
            if (parts.Length != 2 || ends.Length != 2 || !int.TryParse(ends[0], out var a) || !int.TryParse(ends[1], out var b)
                || a < 0 || b < 0 || a >= n || b >= n || a == b)
                throw new ArgumentException($"--bonds: cannot read \"{entry}\" as i-j:J with 0 <= i != j < {n}.");
            var key = (Math.Min(a, b), Math.Max(a, b));
            if (d.ContainsKey(key)) throw new ArgumentException($"--bonds: the bond {key} is given twice.");
            var j = Rat(parts[1], "bonds");
            if (!j.IsZero) d[key] = j;
        }
        return d;
    }

    public static Dictionary<int, (BigRational, BigRational)> ParseFields(string spec, int n)
    {
        var d = new Dictionary<int, (BigRational, BigRational)>();
        foreach (var entry in spec.Split(';', StringSplitOptions.RemoveEmptyEntries | StringSplitOptions.TrimEntries))
        {
            var parts = entry.Split(':');
            var xy = parts.Length == 2 ? parts[1].Split(',') : Array.Empty<string>();
            if (parts.Length != 2 || !int.TryParse(parts[0], out var site) || site < 0 || site >= n || xy.Length != 2)
                throw new ArgumentException($"--fields: cannot read \"{entry}\" as j:x,y with 0 < j < {n}.");
            if (site == 0) throw new ArgumentException("--fields: the lemma's setting has no field on the dephased seat 0.");
            d[site] = (Rat(xy[0], "fields"), Rat(xy[1], "fields"));
        }
        return d;
    }

    private static bool Connected(int n, IReadOnlyDictionary<(int, int), BigRational> bonds)
    {
        var seen = new HashSet<int> { 0 };
        var stack = new Stack<int>(new[] { 0 });
        while (stack.Count > 0)
        {
            int x = stack.Pop();
            foreach (var ((a, b), _) in bonds)
            {
                int y = a == x ? b : b == x ? a : -1;
                if (y >= 0 && seen.Add(y)) stack.Push(y);
            }
        }
        return seen.Count == n;
    }

    // ---- the block and the commutant ------------------------------------------------------------------

    private static List<int> Sector(int n, int k) =>
        Enumerable.Range(0, 1 << n).Where(x => System.Numerics.BitOperations.PopCount((uint)x) == k).ToList();

    /// <summary>The bond part on the sector of k excitations, from the bit action: ZZ is +J on an equal pair and −J on
    /// an unequal one, XX + YY sends an unequal pair to its flip with weight 2J. Bit j of a basis index is site j.</summary>
    public static BigRational[,] BondSector(int n, IReadOnlyDictionary<(int, int), BigRational> bonds, int k, out List<int> basis)
    {
        basis = Sector(n, k);
        var idx = basis.Select((x, i) => (x, i)).ToDictionary(p => p.x, p => p.i);
        var m = Zero(basis.Count, basis.Count);
        foreach (var x in basis)
            foreach (var ((i, j), J) in bonds)
            {
                bool bi = ((x >> i) & 1) == 1, bj = ((x >> j) & 1) == 1;
                m[idx[x], idx[x]] += bi == bj ? J : -J;
                if (bi != bj) m[idx[x ^ (1 << i) ^ (1 << j)], idx[x]] += 2 * J;
            }
        return m;
    }

    private static BigRational[,] Zero(int r, int c)
    {
        var m = new BigRational[r, c];
        for (int a = 0; a < r; a++) for (int b = 0; b < c; b++) m[a, b] = BigRational.Zero;
        return m;
    }

    /// <summary>Unknowns X (rows × cols, flattened row-major) solving Σ equations A·X − X·B = 0 or A·X + X·B = 0; each
    /// pair (A, B, sign) adds rows·cols equations. Returns the dimension of the solution space and a basis.</summary>
    private static List<BigRational[]> Solve(int rows, int cols, params (BigRational[,] A, BigRational[,] B, int Sign)[] eqs)
    {
        int nEq = eqs.Length * rows * cols, nUn = rows * cols;
        var m = Zero(nEq, nUn);
        int e = 0;
        foreach (var (A, B, sign) in eqs)
            for (int r = 0; r < rows; r++)
                for (int c = 0; c < cols; c++, e++)
                {
                    for (int k = 0; k < rows; k++) m[e, k * cols + c] += A[r, k];         // (A X)_rc = Σ_k A_rk X_kc
                    for (int k = 0; k < cols; k++) m[e, r * cols + k] += sign * B[k, c];  // (X B)_rc = Σ_k X_rk B_kc
                }
        return BigRationalLinearAlgebra.Nullspace(m);
    }

    public static BigRational[,] Laplacian(int n, IReadOnlyDictionary<(int, int), BigRational> bonds)
    {
        var l = Zero(n, n);
        foreach (var ((a, b), J) in bonds) { l[a, a] += J; l[b, b] += J; l[a, b] -= J; l[b, a] -= J; }
        return l;
    }

    /// <summary>F157's count at seat 0: N minus the rank of [e_0, 𝓛e_0, …, 𝓛^(N−1)e_0].</summary>
    public static int BlindDimension(int n, IReadOnlyDictionary<(int, int), BigRational> bonds)
    {
        var l = Laplacian(n, bonds);
        var k = Zero(n, n);
        var v = new BigRational[n];
        for (int a = 0; a < n; a++) v[a] = a == 0 ? BigRational.One : BigRational.Zero;
        for (int c = 0; c < n; c++)
        {
            for (int a = 0; a < n; a++) k[a, c] = v[a];
            var w = new BigRational[n];
            for (int a = 0; a < n; a++) { w[a] = BigRational.Zero; for (int b = 0; b < n; b++) w[a] += l[a, b] * v[b]; }
            v = w;
        }
        return n - OperatorSchmidtRank.RankOverQ(k);
    }

    public sealed record BlockReading(int Blind, int BlockDimension, int CommutantDimension, bool XiSolves)
    {
        /// <summary>The lemma at this row: the two sides of the bijection agree, Ξ's block solves, and the space is
        /// one-dimensional exactly where d is not blind.</summary>
        public bool Agrees => BlockDimension == CommutantDimension && XiSolves && (BlockDimension == 1) == (Blind == 0);
    }

    /// <summary>The block reading. <paramref name="withoutJump"/> drops the jump's equation, a planted mutation for
    /// the tests (the solution space then grows where d is not blind).</summary>
    public static BlockReading ReadBlock(int n, IReadOnlyDictionary<(int, int), BigRational> bonds, bool withoutJump = false)
    {
        var b1 = BondSector(n, bonds, 1, out var s1);
        var bn = BondSector(n, bonds, n - 1, out var sn);
        var z1 = Zero(n, n); var zn = Zero(n, n);
        for (int i = 0; i < n; i++)
        {
            z1[i, i] = (s1[i] & 1) == 1 ? -BigRational.One : BigRational.One;
            zn[i, i] = (sn[i] & 1) == 1 ? -BigRational.One : BigRational.One;
        }
        var eqs = withoutJump ? new[] { (bn, b1, -1) } : new[] { (bn, b1, -1), (zn, z1, +1) };
        int block = Solve(n, n, eqs).Count;
        // Ξ's block: |x> -> |x xor 1..1>, from sector 1 onto sector N-1
        int full = (1 << n) - 1;
        var xi = Zero(n, n);
        for (int c = 0; c < n; c++) xi[sn.IndexOf(s1[c] ^ full), c] = BigRational.One;
        bool solves = IsZero(Sub(Mul(bn, xi), Mul(xi, b1))) && IsZero(Add(Mul(zn, xi), Mul(xi, z1)));
        // the other side: operators on one excitation commuting with 𝓛 and P_0 (sites as the basis, not the bit sector)
        var l = Laplacian(n, bonds);
        var p = Zero(n, n); p[0, 0] = BigRational.One;
        int comm = Solve(n, n, (l, l, -1), (p, p, -1)).Count;
        return new BlockReading(BlindDimension(n, bonds), block, comm, solves);
    }

    private static BigRational[,] Mul(BigRational[,] a, BigRational[,] b)
    {
        int r = a.GetLength(0), k = a.GetLength(1), c = b.GetLength(1);
        var m = Zero(r, c);
        for (int i = 0; i < r; i++) for (int j = 0; j < c; j++) for (int t = 0; t < k; t++) m[i, j] += a[i, t] * b[t, j];
        return m;
    }

    private static BigRational[,] Sub(BigRational[,] a, BigRational[,] b) => Combine(a, b, -1);
    private static BigRational[,] Add(BigRational[,] a, BigRational[,] b) => Combine(a, b, +1);

    private static BigRational[,] Combine(BigRational[,] a, BigRational[,] b, int s)
    {
        var m = Zero(a.GetLength(0), a.GetLength(1));
        for (int i = 0; i < a.GetLength(0); i++) for (int j = 0; j < a.GetLength(1); j++) m[i, j] = a[i, j] + s * b[i, j];
        return m;
    }

    private static bool IsZero(BigRational[,] a)
    {
        foreach (var x in a) if (!x.IsZero) return false;
        return true;
    }

    // ---- far elements homogeneous under Z^⊗N ---------------------------------------------------------

    private static (int Phase, int Letter) Mul1(int a, int b)
    {
        if (a == 0) return (0, b);
        if (b == 0) return (0, a);
        if (a == b) return (0, 0);
        return ((b - a + 3) % 3 == 1 ? 1 : 3, 6 - a - b);      // X·Y = iZ, Y·Z = iX, Z·X = iY
    }

    /// <summary>The far elements of the given height parity (0 even, 1 odd): the kernel of W ↦ [H, W] on the Pauli
    /// strings with X or Y at site 0 and that parity of letters X or Y, exactly over ℚ (divided by 2i every entry of the
    /// commutator is ± one coefficient).</summary>
    public static int FarCount(int n, IReadOnlyDictionary<(int, int), BigRational> bonds,
        IReadOnlyDictionary<int, (BigRational X, BigRational Y)> fields, int parity)
    {
        var terms = new List<(BigRational C, int[] S)>();
        foreach (var ((a, b), J) in bonds)
            for (int letter = 1; letter <= 3; letter++)
            {
                var s = new int[n]; s[a] = letter; s[b] = letter; terms.Add((J, s));
            }
        foreach (var (site, (x, y)) in fields)
        {
            if (!x.IsZero) { var s = new int[n]; s[site] = 1; terms.Add((x, s)); }
            if (!y.IsZero) { var s = new int[n]; s[site] = 2; terms.Add((y, s)); }
        }
        var cols = new List<int[]>();
        int total = 1;
        for (int i = 0; i < n; i++) total *= 4;
        for (int code = 0; code < total; code++)
        {
            var s = new int[n];
            int c = code, lit = 0;
            for (int i = n - 1; i >= 0; i--) { s[i] = c % 4; c /= 4; if (s[i] is 1 or 2) lit++; }
            if (s[0] is 1 or 2 && lit % 2 == parity) cols.Add(s);
        }
        var rows = new Dictionary<int, Dictionary<int, BigRational>>();
        for (int col = 0; col < cols.Count; col++)
            foreach (var (coef, t) in terms)
            {
                var s = cols[col];
                int clash = 0;
                for (int i = 0; i < n; i++) if (t[i] != 0 && s[i] != 0 && t[i] != s[i]) clash++;
                if (clash % 2 == 0) continue;
                int phase = 0, key = 0;
                for (int i = 0; i < n; i++) { var (p, l) = Mul1(t[i], s[i]); phase += p; key = key * 4 + l; }
                if (!rows.TryGetValue(key, out var row)) rows[key] = row = new Dictionary<int, BigRational>();
                row[col] = (row.TryGetValue(col, out var v) ? v : BigRational.Zero) + ((phase & 3) == 1 ? coef : -coef);
            }
        if (rows.Count == 0) return cols.Count;
        var m = Zero(rows.Count, cols.Count);
        int r = 0;
        foreach (var row in rows.Values) { foreach (var (c, v) in row) m[r, c] = v; r++; }
        return cols.Count - OperatorSchmidtRank.RankOverQ(m);
    }

    /// <summary>Whether every nonzero field lies on one line through the origin (true with fewer than two).</summary>
    public static bool OnOneLine(IReadOnlyDictionary<int, (BigRational X, BigRational Y)> fields)
    {
        var vs = fields.Values.Where(f => !f.X.IsZero || !f.Y.IsZero).ToList();
        return vs.All(a => vs.All(b => (a.X * b.Y - a.Y * b.X).IsZero));
    }

    public sealed record FarReading(int N, int Blind, bool OnLine, bool AnyField, bool OnLocus, int Odd, int Even)
    {
        /// <summary>Where the proof predicts: d not blind, and N = 3, or N = 4 with some field nonzero.</summary>
        public bool Predicted => Blind == 0 && (N == 3 || (N == 4 && AnyField));

        /// <summary>At N = 3: even 0, odd nonzero exactly on one line. At N = 4: (even, odd) = (0, 0) off one line, (1, 0)
        /// on it, and (2, 0) on the antiparallel locus. Null where nothing is predicted.</summary>
        public bool? Agrees => !Predicted ? null
            : N == 3 ? Even == 0 && (Odd > 0) == OnLine
            : OnLine ? Even == (OnLocus ? 2 : 1) && Odd == 0 : Even == 0 && Odd == 0;
    }

    /// <summary>The antiparallel locus at N = 4: one site k of 1, 2, 3 without a field, the fields on the other two, i
    /// and j, opposite and nonzero, and J_0i = J_ki = −J_0j = −J_kj ≠ 0 (site k coupled to i and j exactly as site 0 is).
    /// There the far elements homogeneous under Z^⊗4 at a seat that is not blind are class A's colouring and
    /// (X_0X_k + Y_0Y_k) ⊗ |𝓉⟩⟨𝓉|, 𝓉 the triplet state of i and j with one excitation.</summary>
    public static bool OnAntiparallelLocus(int n, IReadOnlyDictionary<(int, int), BigRational> bonds,
        IReadOnlyDictionary<int, (BigRational X, BigRational Y)> fields)
    {
        if (n != 4) return false;
        var nz = Enumerable.Range(1, 3).Where(j => fields.TryGetValue(j, out var f) && (!f.X.IsZero || !f.Y.IsZero)).ToList();
        if (nz.Count != 2) return false;
        int i = nz[0], j = nz[1], k = 6 - i - j;
        var (fi, fj) = (fields[i], fields[j]);
        if (!(fi.X + fj.X).IsZero || !(fi.Y + fj.Y).IsZero) return false;
        BigRational J(int a, int b) => bonds.TryGetValue((Math.Min(a, b), Math.Max(a, b)), out var v) ? v : BigRational.Zero;
        var c = J(0, i);
        return !c.IsZero && J(k, i) == c && J(0, j) == -c && J(k, j) == -c;
    }

    public static FarReading ReadFar(int n, IReadOnlyDictionary<(int, int), BigRational> bonds,
        IReadOnlyDictionary<int, (BigRational X, BigRational Y)> fields)
    {
        int blind = BlindDimension(n, bonds);
        bool line = OnOneLine(fields);
        bool any = fields.Values.Any(f => !f.X.IsZero || !f.Y.IsZero);
        int odd = FarCount(n, bonds, fields, 1), even = FarCount(n, bonds, fields, 0);
        return new FarReading(n, blind, line, any, OnAntiparallelLocus(n, bonds, fields), odd, even);
    }

    // ---- the descent's two conditions -----------------------------------------------------------------

    /// <summary>Condition (II): for each sector k = 1 .. N − 1 the dimension of the vectors with site 0 excited in every
    /// term at the vacuum's energy E0 (the sum of the couplings), exactly over ℚ.</summary>
    public static int[] FenceDimensions(int n, IReadOnlyDictionary<(int, int), BigRational> bonds)
    {
        var e0 = bonds.Values.Aggregate(BigRational.Zero, (a, b) => a + b);
        var dims = new int[n - 1];
        for (int k = 1; k < n; k++)
        {
            var bk = BondSector(n, bonds, k, out var basis);
            var cols = Enumerable.Range(0, basis.Count).Where(i => (basis[i] & 1) == 1).ToList();
            var m = Zero(basis.Count, cols.Count);
            for (int r = 0; r < basis.Count; r++)
                for (int c = 0; c < cols.Count; c++)
                    m[r, c] = bk[r, cols[c]] - (r == cols[c] ? e0 : BigRational.Zero);
            dims[k - 1] = BigRationalLinearAlgebra.Nullspace(m).Count;
        }
        return dims;
    }

    /// <summary>Condition (I): the complex dimension of the cyclic space alg(H_b, Z_0, F⁺)|vac⟩, F⁺ = Σ_j ω_j·σ⁻_j with
    /// ω_j = x_j + i·y_j, exactly: the ℚ-span of the realified orbit closed under H_b, Z_0, F⁺ and multiplication by i
    /// is the realification of the complex span, so its ℚ-dimension is twice the complex one.</summary>
    public static int CyclicDimension(int n, IReadOnlyDictionary<(int, int), BigRational> bonds,
        IReadOnlyDictionary<int, (BigRational X, BigRational Y)> fields)
    {
        int dim = 1 << n;
        var h = Zero(dim, dim);
        for (int x = 0; x < dim; x++)
            foreach (var ((i, j), J) in bonds)
            {
                bool bi = ((x >> i) & 1) == 1, bj = ((x >> j) & 1) == 1;
                h[x, x] += bi == bj ? J : -J;
                if (bi != bj) h[x ^ (1 << i) ^ (1 << j), x] += 2 * J;
            }
        BigRational[] Apply(BigRational[] v, int op)
        {
            var o = new BigRational[2 * dim];
            for (int a = 0; a < 2 * dim; a++) o[a] = BigRational.Zero;
            switch (op)
            {
                case 0:                                            // H_b, real
                    for (int a = 0; a < dim; a++)
                        for (int b = 0; b < dim; b++)
                            if (!h[a, b].IsZero) { o[a] += h[a, b] * v[b]; o[dim + a] += h[a, b] * v[dim + b]; }
                    break;
                case 1:                                            // Z_0
                    for (int a = 0; a < dim; a++) { int s = (a & 1) == 1 ? -1 : 1; o[a] = s * v[a]; o[dim + a] = s * v[dim + a]; }
                    break;
                case 2:                                            // F+ = A + iB: (A re − B im, A im + B re)
                    foreach (var (j, (fx, fy)) in fields)
                        for (int x = 0; x < dim; x++)
                            if (((x >> j) & 1) == 0)
                            {
                                int y = x | (1 << j);
                                o[y] += fx * v[x] - fy * v[dim + x];
                                o[dim + y] += fx * v[dim + x] + fy * v[x];
                            }
                    break;
                default:                                           // multiplication by i
                    for (int a = 0; a < dim; a++) { o[a] = -v[dim + a]; o[dim + a] = v[a]; }
                    break;
            }
            return o;
        }
        var basis = new List<(BigRational[] V, int Pivot)>();
        var queue = new Stack<BigRational[]>();
        var vac = new BigRational[2 * dim];
        for (int a = 0; a < 2 * dim; a++) vac[a] = a == 0 ? BigRational.One : BigRational.Zero;
        queue.Push(vac);
        while (queue.Count > 0)
        {
            var v = queue.Pop();
            var r = (BigRational[])v.Clone();
            foreach (var (bv, piv) in basis)
                if (!r[piv].IsZero) { var c = r[piv]; for (int a = 0; a < r.Length; a++) r[a] -= c * bv[a]; }
            int p = Array.FindIndex(r, z => !z.IsZero);
            if (p < 0) continue;
            var inv = BigRational.One / r[p];
            for (int a = 0; a < r.Length; a++) r[a] *= inv;
            foreach (var (bv, _) in basis)
                if (!bv[p].IsZero) { var c = bv[p]; for (int a = 0; a < r.Length; a++) bv[a] -= c * r[a]; }
            basis.Add((r, p));
            for (int op = 0; op < 4; op++) queue.Push(Apply(v, op));
        }
        return basis.Count / 2;                                   // the realified span is closed under i, so its dimension is even
    }

    public sealed record DescentReading(int Blind, int Cyclic, int Full, int[] Fence)
    {
        /// <summary>The sectors at which (II) is proven where d is not blind: 1, 2 and N − 1.</summary>
        public bool FenceAgrees => Blind > 0 || new[] { 1, 2, Fence.Length }.Distinct().All(k => Fence[k - 1] == 0);
    }

    public static DescentReading ReadDescent(int n, IReadOnlyDictionary<(int, int), BigRational> bonds,
        IReadOnlyDictionary<int, (BigRational X, BigRational Y)> fields) =>
        new(BlindDimension(n, bonds), CyclicDimension(n, bonds, fields), 1 << n, FenceDimensions(n, bonds));

    // ---- fixed grids -----------------------------------------------------------------------------------

    /// <summary>A deterministic generator (a 64-bit LCG), so the grids are the same on every machine.</summary>
    private sealed class Lcg(ulong seed)
    {
        private ulong _s = seed;
        public int Next(int bound) { _s = _s * 6364136223846793005UL + 1442695040888963407UL; return (int)((_s >> 33) % (ulong)bound); }
    }

    private static readonly int[] Couplings = { -3, -2, -1, 1, 2, 3 };
    private static readonly (int X, int Y)[] FieldList =
        { (1, 0), (0, 1), (-1, 0), (3, 4), (4, 3), (-3, 4), (5, 0), (1, 1), (-1, 1), (2, 0), (0, 2), (0, 0) };

    public static IEnumerable<Dictionary<(int, int), BigRational>> BondGrid(int n, int count, ulong seed)
    {
        var g = new Lcg(seed);
        int made = 0;
        while (made < count)
        {
            var b = new Dictionary<(int, int), BigRational>();
            for (int i = 0; i < n; i++)
                for (int j = i + 1; j < n; j++)
                    if (g.Next(10) < 6) b[(i, j)] = new BigRational(Couplings[g.Next(Couplings.Length)], 1);
            if (made % 3 == 2)                                    // every third row symmetric under the swap 1 <-> 2: d blind
            {
                var sym = new Dictionary<(int, int), BigRational>();
                foreach (var ((i, j), J) in b)
                {
                    int Map(int s) => s == 1 ? 2 : s == 2 ? 1 : s;
                    sym[(i, j)] = J;
                    sym[(Math.Min(Map(i), Map(j)), Math.Max(Map(i), Map(j)))] = J;
                }
                b = sym;
            }
            if (!Connected(n, b)) continue;
            made++;
            yield return b;
        }
    }

    public sealed record BlockGridReading(int Rows, int NotBlind, int Blind, int Disagreements);

    public static BlockGridReading ReadBlockGrid(int n, int count = 24, bool withoutJump = false)
    {
        int rows = 0, nb = 0, bl = 0, bad = 0;
        foreach (var b in BondGrid(n, count, 20261118UL + (ulong)n))
        {
            var r = ReadBlock(n, b, withoutJump);
            rows++;
            if (r.Blind == 0) nb++; else bl++;
            if (!r.Agrees) bad++;
        }
        return new BlockGridReading(rows, nb, bl, bad);
    }

    public sealed record FarGridReading(int Rows, int NotBlindOffLine, int NotBlindOnLine, int OffLineWithElement,
        int OnLineWithElement, int Disagreements);

    public static FarGridReading ReadFarGrid(int n, int count)
    {
        var g = new Lcg(777UL + (ulong)n);
        int rows = 0, off = 0, on = 0, offHit = 0, onHit = 0, bad = 0;
        foreach (var b in BondGrid(n, count, 4242UL + (ulong)n))
        {
            var f = new Dictionary<int, (BigRational, BigRational)>();
            for (int j = 1; j < n; j++)
            {
                var (x, y) = FieldList[g.Next(FieldList.Length)];
                f[j] = (new BigRational(x, 1), new BigRational(y, 1));
            }
            rows++;
            if (BlindDimension(n, b) > 0) continue;
            var r = ReadFar(n, b, f);
            bool hit = r.Odd + r.Even > 0;
            if (r.OnLine) { on++; if (hit) onHit++; } else { off++; if (hit) offHit++; }
            if (r.Agrees == false) bad++;
        }
        return new FarGridReading(rows, off, on, offHit, onHit, bad);
    }

    // ---- IInspectable -----------------------------------------------------------------------------------

    private BlockReading? _block;
    private FarReading? _far;
    private DescentReading? _descent;
    public BlockReading Block => _block ??= ReadBlock(N, Bonds);
    public FarReading Far => _far ??= ReadFar(N, Bonds, Fields);
    public DescentReading Descent => _descent ??= ReadDescent(N, Bonds, Fields);

    private string RowText => string.Join(", ", Bonds.Select(p => $"{p.Key.Item1}-{p.Key.Item2}:{p.Value}")) +
        (Fields.Count == 0 ? "; no fields" : "; fields " + string.Join(", ", Fields.Select(p => $"{p.Key}:({p.Value.X}, {p.Value.Y})")));

    public string DisplayName => $"the block from one excitation to N − 1 at N = {N} ({RowText})";

    public string Summary
    {
        get
        {
            var b = Block;
            return string.Format(CultureInfo.InvariantCulture,
                "seat 0 blind of dimension {0}; blocks solving the bond and jump conditions: {1}, the commutant of 𝓛 and P_0: {2}, " +
                "X^N's block {3}; the lemma {4}",
                b.Blind, b.BlockDimension, b.CommutantDimension, b.XiSolves ? "solves" : "does NOT solve", b.Agrees ? "holds" : "FAILS");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var b = Block;
            yield return new InspectableNode("the block at this row",
                summary: string.Format(CultureInfo.InvariantCulture,
                    "Solutions 𝒯 of B_(N−1)·𝒯 = 𝒯·B_1 and Z_0·𝒯 = −𝒯·Z_0, the bond blocks built from the bit action: {0}. " +
                    "Operators on one excitation commuting with the signed Laplacian and P_0 (F157's Theorem A count): {1}. " +
                    "X^N's block {2}. Seat 0 blind of dimension {3}, so the lemma predicts {4}; it {5}.",
                    b.BlockDimension, b.CommutantDimension, b.XiSolves ? "solves both" : "does NOT solve both", b.Blind,
                    b.Blind == 0 ? "one dimension, X^N's block" : "more than one, the commutant count",
                    b.Agrees ? "holds" : "FAILS"),
                provenance: NodeProvenance.Live);
            var f = Far;
            yield return new InspectableNode("far elements homogeneous under Z^N at this row",
                summary: string.Format(CultureInfo.InvariantCulture,
                    "Of odd height {0}, of even height {1} (exact over Q); the fields {2} on one line. {3}",
                    f.Odd, f.Even, f.OnLine ? "lie" : "do not lie",
                    !f.Predicted
                        ? (f.Blind > 0 ? "Seat 0 is blind, where the bullets predict nothing."
                           : N == 4 ? "No field is nonzero; nothing is predicted."
                           : "At N = 5 nothing is predicted; the proof measures that none exists with seat 0 not blind and the fields not on one line.")
                        : N == 3
                            ? $"At N = 3 with seat 0 not blind the bullet says: none of even height, and of odd height exactly with the fields on one line; it {(f.Agrees == true ? "holds" : "FAILS")}."
                            : $"At N = 4 with seat 0 not blind the bullet says: none of either height off one line, exactly one of even height (class A's colouring) on it, and two of even height on the antiparallel locus, class A's colouring and (X_0X_k + Y_0Y_k) (x) |t><t|; the row is {(f.OnLocus ? "on" : "off")} the locus, and the bullet {(f.Agrees == true ? "holds" : "FAILS")}."),
                provenance: NodeProvenance.Live);
            var dsc = Descent;
            yield return new InspectableNode("the descent's two conditions at this row",
                summary: string.Format(CultureInfo.InvariantCulture,
                    "(I) the cyclic space alg(H_b, Z_0, F+)|vac> has dimension {0} of {1} (exact over Q(i)); (II) vectors with " +
                    "seat 0 excited in every term at the vacuum's energy, by sector 1..{2}: [{3}] (exact). {4}",
                    dsc.Cyclic, dsc.Full, N - 1, string.Join(", ", dsc.Fence),
                    dsc.Blind > 0 ? "Seat 0 is blind; the bullet predicts nothing."
                    : $"Seat 0 not blind: (II) is proven empty at the sectors 1, 2 and N - 1, here {(dsc.FenceAgrees ? "empty" : "NOT empty")}" +
                      (N == 4 && f.Predicted
                          ? (f.OnLocus
                              ? $", and on the antiparallel locus (I) fails by the case analysis, the cyclic space 14 of 16, here {dsc.Cyclic}, so it {(dsc.Cyclic == 14 ? "holds" : "FAILS")}."
                              : $", and at N = 4 off the antiparallel locus (I) holds by the case analysis, here {(dsc.Cyclic == dsc.Full ? "full" : "NOT full")}.")
                          : ".")),
                provenance: NodeProvenance.Live);
            var grids = new[] { 3, 4, 5 }.Select(n => (n, ReadBlockGrid(n))).ToList();
            yield return new InspectableNode("the block on fixed grids at N = 3, 4, 5",
                summary: string.Join("; ", grids.Select(p => string.Format(CultureInfo.InvariantCulture,
                    "N = {0}: {1} rows ({2} with seat 0 not blind, {3} blind), {4} disagreeing",
                    p.n, p.Item2.Rows, p.Item2.NotBlind, p.Item2.Blind, p.Item2.Disagreements))),
                provenance: NodeProvenance.Live);
            var f3 = ReadFarGrid(3, 30);
            var f4 = ReadFarGrid(4, 20);
            yield return new InspectableNode("far elements on fixed grids at N = 3 and 4",
                summary: string.Format(CultureInfo.InvariantCulture,
                    "Rows with seat 0 not blind. N = 3: {0} with the fields off one line ({1} carrying a homogeneous far " +
                    "element), {2} on one line ({3} carrying one), {4} disagreeing with the bullet. N = 4: {5} off one line " +
                    "({6} carrying one), {7} on one line ({8} carrying one or more), {9} disagreeing with the bullet " +
                    "on the rows it predicts.",
                    f3.NotBlindOffLine, f3.OffLineWithElement, f3.NotBlindOnLine, f3.OnLineWithElement, f3.Disagreements,
                    f4.NotBlindOffLine, f4.OffLineWithElement, f4.NotBlindOnLine, f4.OnLineWithElement, f4.Disagreements),
                provenance: NodeProvenance.Live);
            yield return new InspectableNode("where the proof lives",
                summary: "The lemma and its N = 3 conclusions are argued by hand in the bullet 'The block from one excitation " +
                         "to N − 1, at every N' of docs/proofs/PROOF_ONE_DEPHASED_SEAT.md, gate " +
                         "simulations/blind_seat_block_lemma_gate.py, and N = 4 with the descent's two conditions in the bullet " +
                         "'At N = 4, a seat that is not blind', gate simulations/blind_seat_descent_gate.py; this witness reads " +
                         "them at points on a road of its own. " +
                         "E pairs at a seat that is not blind through a far element of no definite height, which neither reaches.");
        }
    }
}
