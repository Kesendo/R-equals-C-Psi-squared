using System.Numerics;
using MirrorWorld;
using static MirrorWorldTests.LetterCubeTests;

namespace MirrorWorldTests;

// From-below guard for the sandwich sum and the sextant's intertwining (built 2026-11-07): the two ladders Phi and S+ of F125 and F142 held exactly on the letters as sums of sandwiches,
// and Sextant.Intertwines reading M·L = Λ·M between two settings, c read off the identity. Two routes meet twice: the
// map on the letters against the same sandwiches built densely from Jordan-Wigner strings at N = 2 and 3, and the
// verdict against the identity computed densely in the string basis. The readings are the registry's: F125's
// Theorem B at any bond and rate profile with its two traps as the controls, F142's Lemma 2.3 with its breaks, F125
// §5 on the ring by height parity, F166's bra flip joining the two ladders, the spin ladder's loss under ZZ (F171's scope)
// and Φ's (F125's Δ-tightness), and the
// composition R∘Φ = S⁺∘R, a reflection with a kernel, read as one map on the full generator.
public class SandwichSumTests
{
    static readonly LetterCube.Gauss[,] Z = LetterCube.LetterMatrix('Z');

    // ---- settings: the XY chain or ring in the sextant's string basis ----

    // a term of H as a string: the named letters on their sites, I elsewhere (not SandwichSum's Term, a sandwich of the sum)
    static PauliString HTerm(int n, params (int Site, char Letter)[] letters)
    {
        var c = Enumerable.Repeat('I', n).ToArray();
        foreach (var (site, letter) in letters) c[site] = letter;
        return PauliString.Parse(new string(c));
    }

    // H = Σ_b J_b (XX + YY) on the chain's bonds (hop 2J_b), a wrap bond J_w (X_0 X_{n-1} + Y_0 Y_{n-1}) when given,
    // Z fields h_l, a ZZ term on every chain bond, and the YY coefficients replaced by bondsY when given (anisotropy);
    // Z dephasing at the rates, a zero rate removing its jump
    static Sextant.Setting XY(int n, long[] bonds, long[] rates, long? wrap = null, long[]? fields = null, long zz = 0,
                              long[]? bondsY = null, IEnumerable<(PauliString, long)>? more = null)
    {
        var h = new List<(PauliString, long)>();
        for (int b = 0; b < n - 1; b++)
        {
            h.Add((HTerm(n,(b, 'X'), (b + 1, 'X')), bonds[b]));
            h.Add((HTerm(n,(b, 'Y'), (b + 1, 'Y')), (bondsY ?? bonds)[b]));
            if (zz != 0) h.Add((HTerm(n,(b, 'Z'), (b + 1, 'Z')), zz));
        }
        if (wrap is { } w)
        {
            h.Add((HTerm(n,(0, 'X'), (n - 1, 'X')), w));
            h.Add((HTerm(n,(0, 'Y'), (n - 1, 'Y')), w));
        }
        if (fields is not null)
            for (int l = 0; l < n; l++) if (fields[l] != 0) h.Add((HTerm(n,(l, 'Z')), fields[l]));
        if (more is not null) h.AddRange(more);
        var jumps = Enumerable.Range(0, n).Where(l => rates[l] != 0).Select(l => (l, Z, rates[l]));
        return new Sextant.Setting(n, h.Where(t => t.Item2 != 0), jumps);
    }

    static long[] Draw(Random r, int count, int lo, int hi) => Enumerable.Range(0, count).Select(_ => (long)r.Next(lo, hi + 1)).ToArray();

    static long[] Uniform(int count, long v) => Enumerable.Repeat(v, count).ToArray();

    static LetterCube.Move R(int n) => LetterCube.TableMoves(n).Single(m => m.Name.StartsWith("R,"));

    static Sextant.Intertwining Read(Sextant.Setting from, Sextant.Setting to, SandwichSum m, int? parity = null)
        => Sextant.Intertwines(from, to, m, parity);

    static void AssertCopy(Sextant.Intertwining r)
    {
        Assert.Equal((Sextant.OnHamiltonian.Keeps, Sextant.OnDissipator.Keeps, Sextant.Verdict.Copy, ""), (r.OnH, r.OnD, r.Verdict, r.C));
        Assert.True(r.Rank > 0, $"{r.Move}: an empty identity");
    }

    // ---- the dense route: the sandwiches from Jordan-Wigner strings, the generator from H and the jumps ----

    static Complex[,] Kron(Complex[,] a, Complex[,] b)
    {
        int ra = a.GetLength(0), rb = b.GetLength(0);
        var k = new Complex[ra * rb, ra * rb];
        for (int i = 0; i < ra; i++) for (int j = 0; j < ra; j++)
            for (int p = 0; p < rb; p++) for (int q = 0; q < rb; q++)
                k[i * rb + p, j * rb + q] = a[i, j] * b[p, q];
        return k;
    }

    static Complex[,] Product(IEnumerable<LetterCube.Gauss[,]> factors)
    {
        Complex[,]? m = null;
        foreach (var f in factors) m = m is null ? ToComplex(f) : Kron(m, ToComplex(f));
        return m!;
    }

    static Complex[,] Transpose(Complex[,] a)
    {
        int d = a.GetLength(0);
        var r = new Complex[d, d];
        for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) r[i, j] = a[j, i];
        return r;
    }

    static Complex[,] Dagger(Complex[,] a)
    {
        int d = a.GetLength(0);
        var r = new Complex[d, d];
        for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) r[i, j] = Complex.Conjugate(a[j, i]);
        return r;
    }

    static Complex[,] Lin(Complex[,] a, Complex ka, Complex[,] b, Complex kb)
    {
        int d = a.GetLength(0);
        var r = new Complex[d, d];
        for (int i = 0; i < d; i++) for (int j = 0; j < d; j++) r[i, j] = ka * a[i, j] + kb * b[i, j];
        return r;
    }

    static Complex Trace(Complex[,] m) { Complex t = 0; for (int i = 0; i < m.GetLength(0); i++) t += m[i, i]; return t; }

    // the left and right factor of site m in A rho B with A = Z_{<l} a_l, B = Z_{<l'} b_{l'}: the Jordan-Wigner string in
    // the convention of the sources and of SpectatorIntertwiner (the string on the sites before l, the stagger +1 at
    // site 0); the dense route is independent of the piece algebra, not of the convention
    static (LetterCube.Gauss[,], LetterCube.Gauss[,]) Jw(int m, int l, LetterCube.Gauss[,] a, int lp, LetterCube.Gauss[,] b)
        => (m == l ? a : m < l ? Z : LetterCube.LetterMatrix('I'), m == lp ? b : m < lp ? Z : LetterCube.LetterMatrix('I'));

    // a dense sum of sandwiches: terms (weight, A, B)
    sealed record DenseSum(int N, List<(Complex W, Complex[,] A, Complex[,] B)> Terms)
    {
        public Complex[,] Apply(Complex[,] rho)
        {
            int d = 1 << N;
            var r = new Complex[d, d];
            foreach (var (w, a, b) in Terms) r = Lin(r, 1, Mul(a, Mul(rho, b)), w);
            return r;
        }
    }

    static DenseSum DensePhi(int n, bool strings = true)
        => new(n, Enumerable.Range(0, n).Select(l => (Complex.One,
               Product(Enumerable.Range(0, n).Select(m => strings ? Jw(m, l, LetterCube.Raising, l, LetterCube.Lowering).Item1 : m == l ? LetterCube.Raising : LetterCube.LetterMatrix('I'))),
               Product(Enumerable.Range(0, n).Select(m => strings ? Jw(m, l, LetterCube.Raising, l, LetterCube.Lowering).Item2 : m == l ? LetterCube.Lowering : LetterCube.LetterMatrix('I'))))).ToList());

    static DenseSum DenseSPlus(int n, bool staggered = true)
        => new(n, Enumerable.Range(0, n).Select(l => ((Complex)(staggered && l % 2 == 1 ? -1 : 1),
               Product(Enumerable.Range(0, n).Select(m => Jw(m, l, LetterCube.Raising, l, LetterCube.Raising).Item1)),
               Product(Enumerable.Range(0, n).Select(m => Jw(m, l, LetterCube.Raising, l, LetterCube.Raising).Item2)))).ToList());

    static DenseSum DenseSingleMode(int n, long[] phi)
    {
        var terms = new List<(Complex, Complex[,], Complex[,])>();
        for (int l = 0; l < n; l++)
            for (int lp = 0; lp < n; lp++)
                if (phi[l] * phi[lp] != 0)
                    terms.Add((phi[l] * phi[lp],
                               Product(Enumerable.Range(0, n).Select(m => Jw(m, l, LetterCube.Raising, lp, LetterCube.Lowering).Item1)),
                               Product(Enumerable.Range(0, n).Select(m => Jw(m, l, LetterCube.Raising, lp, LetterCube.Lowering).Item2))));
        return new DenseSum(n, terms);
    }

    // the generator and its adjoint, densely: L(rho) = -i[H, rho] + Σ r (J rho J† - {J†J, rho}/2)
    static Func<Complex[,], bool, Complex[,]> DenseGenerator(Sextant.Setting s)
    {
        int d = 1 << s.N;
        var h = new Complex[d, d];
        foreach (var (term, c) in s.H) h = Lin(h, 1, DenseString(term.ToString(s.N)), c);
        var jumps = s.Jumps.Select(j => (Product(Enumerable.Range(0, s.N).Select(m => m == j.Site ? j.Jump : LetterCube.LetterMatrix('I'))), (Complex)j.Rate)).ToList();
        return (rho, adjoint) =>
        {
            var r = Lin(Mul(h, rho), adjoint ? Complex.ImaginaryOne : -Complex.ImaginaryOne, Mul(rho, h), adjoint ? -Complex.ImaginaryOne : Complex.ImaginaryOne);
            foreach (var (j, rate) in jumps)
            {
                var jd = Dagger(j);
                var jdj = Mul(jd, j);
                var sandwich = adjoint ? Mul(jd, Mul(rho, j)) : Mul(j, Mul(rho, jd));
                r = Lin(r, 1, sandwich, rate);
                r = Lin(r, 1, Lin(Mul(jdj, rho), 1, Mul(rho, jdj), 1), -rate / 2);
            }
            return r;
        };
    }

    // the Pauli expansion of a dense matrix, times 2^N: coefficient Tr(S X)
    static Dictionary<PauliString, Complex> Expand(Complex[,] x, int n)
    {
        var r = new Dictionary<PauliString, Complex>();
        foreach (var s in LetterCube.Strings(n))
        {
            var c = Trace(Mul(DenseString(s.ToString(n)), x));
            if (c != Complex.Zero) r[s] = c;
        }
        return r;
    }

    static void AssertSameMap(IReadOnlyDictionary<PauliString, LetterCube.Gauss> a, IReadOnlyDictionary<PauliString, LetterCube.Gauss> b, string what)
    {
        Assert.True(a.Count == b.Count, $"{what}: {a.Count} strings against {b.Count}");
        foreach (var (s, g) in a)
            Assert.True(b.TryGetValue(s, out var h) && h == g, $"{what}: {s} carries {g} against {(b.TryGetValue(s, out var hh) ? hh.ToString() : "nothing")}");
    }

    static void AssertSameExpansion(IReadOnlyDictionary<PauliString, LetterCube.Gauss> letters, Dictionary<PauliString, Complex> dense, string what)
    {
        Assert.True(letters.Count == dense.Count, $"{what}: {letters.Count} strings on the letters, {dense.Count} densely");
        foreach (var (s, g) in letters)
            Assert.True(dense.TryGetValue(s, out var c) && c == new Complex(g.Re, g.Im), $"{what}: {s} carries {g} on the letters, {(dense.TryGetValue(s, out var cc) ? cc : 0)} densely");
    }

    // The sandwich sum on the letters against the same sandwiches built densely from the strings Z_{<l} and sigma^± at
    // N = 2 and 3: Phi, S+, both adjoints, the eta ladder without its strings, the spin ladder without its stagger,
    // the single mode and Phi composed with R, on every string, exactly (the dense entries are small integers).
    [Fact]
    public void The_Sandwich_Sum_Agrees_With_The_Dense_Sandwiches()
    {
        foreach (int n in new[] { 2, 3 })
        {
            long scale = 1L << n;
            var phi = DensePhi(n);
            var r = R(n);
            var cases = new List<(SandwichSum Sum, Func<Complex[,], Complex[,]> Dense)>
            {
                (SandwichSum.EtaLadder(n), rho => phi.Apply(rho)),
                (SandwichSum.SpinLadder(n), rho => DenseSPlus(n).Apply(rho)),
                (SandwichSum.EtaLadder(n, strings: false), rho => DensePhi(n, strings: false).Apply(rho)),
                (SandwichSum.SpinLadder(n, staggered: false), rho => DenseSPlus(n, staggered: false).Apply(rho)),
                (SandwichSum.SingleMode(n, n == 2 ? new long[] { 1, -1 } : new long[] { 1, 0, -1 }), rho => DenseSingleMode(n, n == 2 ? new long[] { 1, -1 } : new long[] { 1, 0, -1 }).Apply(rho)),
                // the adjoint: sigma -> Σ A† sigma B†
                (SandwichSum.EtaLadder(n).Adjoint, rho => new DenseSum(n, phi.Terms.Select(t => (Complex.Conjugate(t.W), Dagger(t.A), Dagger(t.B))).ToList()).Apply(rho)),
                (SandwichSum.SpinLadder(n).Adjoint, rho => new DenseSum(n, DenseSPlus(n).Terms.Select(t => (Complex.Conjugate(t.W), Dagger(t.A), Dagger(t.B))).ToList()).Apply(rho)),
                // Phi then R: rho -> Phi(rho)·X^N
                (SandwichSum.EtaLadder(n).Then(r), rho => Mul(phi.Apply(rho), DenseString(new string('X', n)))),
                // a transposed sandwich, ⊗σ⁺ · ρ^T · X^N, and its adjoint σ -> X^N · σ^T · ⊗σ⁺ (the factors are real)
                (new SandwichSum("a transposed sandwich", n, new[] { SandwichSum.Sandwich(LetterCube.Gauss.One,
                     Enumerable.Range(0, n).Select(_ => (LetterCube.Raising, LetterCube.LetterMatrix('X'))).ToList(), transpose: true) }),
                 rho => Mul(Product(Enumerable.Repeat(LetterCube.Raising, n)), Mul(Transpose(rho), DenseString(new string('X', n))))),
                (new SandwichSum("a transposed sandwich", n, new[] { SandwichSum.Sandwich(LetterCube.Gauss.One,
                     Enumerable.Range(0, n).Select(_ => (LetterCube.Raising, LetterCube.LetterMatrix('X'))).ToList(), transpose: true) }).Adjoint,
                 rho => Mul(DenseString(new string('X', n)), Mul(Transpose(rho), Product(Enumerable.Repeat(LetterCube.Raising, n))))),
            };
            foreach (var (sum, dense) in cases)
                foreach (var s in LetterCube.Strings(n))
                {
                    var image = dense(DenseString(s.ToString(n)));
                    AssertSameExpansion(sum.Apply(s), Expand(image, n), $"{sum.Name} on {s.ToString(n)} at N = {n}");
                }
            Assert.Equal(scale, SandwichSum.EtaLadder(n).Scale);
        }
    }

    // F166 on the letters: the shift R, rho -> rho X^N, turns the right-hand sigma^- into sigma^+ and every Z of the
    // right-hand string into -Z, one sign per site below l, the stagger; so R Phi R = S+ exactly, with the sites
    // counted from 0 (from 1 it is -S+), and R S+ R = Phi back. Read as maps, on every string, at N = 2 to 4.
    [Fact]
    public void The_Bra_Flip_Joins_The_Two_Ladders()
    {
        foreach (int n in new[] { 2, 3, 4 })
        {
            var r = SandwichSum.FromMove(R(n));
            var viaR = r.Then(SandwichSum.EtaLadder(n)).Then(r);
            var back = r.Then(SandwichSum.SpinLadder(n)).Then(r);
            foreach (var s in LetterCube.Strings(n))
            {
                AssertSameMap(SandwichSum.SpinLadder(n).Apply(s), viaR.Apply(s), $"R Phi R on {s.ToString(n)}");
                AssertSameMap(SandwichSum.EtaLadder(n).Apply(s), back.Apply(s), $"R S+ R on {s.ToString(n)}");
            }
        }
    }

    // F125 Theorem B read exactly: Phi intertwines the dephased XY chain at any bond profile, any rate profile and
    // any Z fields (on-site potentials), both parts, and so does its adjoint W†; the ZZ term (XXZ, Heisenberg) loses
    // the Hamiltonian part and keeps the dissipator (F125's delta-tightness; at N = 2 the ZZ energy z_0 z_1 is a function
    // of the particle number alone, so Φ, adding a particle to both sides, turns its sign: ad_ZZ flipped exactly, which
    // beside the kept hop breaks the whole), and a pairing term, the
    // anisotropy XX - YY, loses it too (the H class is the particle-conserving quadratic one). N = 2 to 5, drawn
    // profiles with zero bonds and zero rates among them.
    [Fact]
    public void Phi_Intertwines_The_XY_Chain_At_Any_Profile()
    {
        var rnd = new Random(2026_11_07);
        foreach (int n in new[] { 2, 3, 4, 5 })
            for (int draw = 0; draw < 3; draw++)
            {
                var bonds = Draw(rnd, n - 1, -3, 3);
                var rates = Draw(rnd, n, 0, 3);
                var fields = Draw(rnd, n, -2, 2);
                if (bonds.All(b => b == 0)) bonds[0] = 1;
                if (rates.All(g => g == 0)) rates[0] = 1;
                var chain = XY(n, bonds, rates, fields: fields);
                AssertCopy(Read(chain, chain, SandwichSum.EtaLadder(n)));
                AssertCopy(Read(chain, chain, SandwichSum.EtaLadder(n).Adjoint));
                var heisenberg = XY(n, bonds, rates, zz: 1);
                var lost = Read(heisenberg, heisenberg, SandwichSum.EtaLadder(n));
                Assert.Equal((Sextant.OnHamiltonian.Breaks, Sextant.OnDissipator.Keeps, Sextant.Verdict.Neither), (lost.OnH, lost.OnD, lost.Verdict));
                if (n == 2)
                {
                    var zzOnly = XY(2, new long[] { 0 }, rates, zz: 1);
                    Assert.Equal((Sextant.OnHamiltonian.Flips, Sextant.OnDissipator.Keeps), (Read(zzOnly, zzOnly, SandwichSum.EtaLadder(2)).OnH, Read(zzOnly, zzOnly, SandwichSum.EtaLadder(2)).OnD));
                }
                var paired = XY(n, bonds, rates, bondsY: bonds.Select(b => -b).ToArray());      // XX - YY = 2(σ⁺σ⁺ + σ⁻σ⁻)
                var pairing = Read(paired, paired, SandwichSum.EtaLadder(n));
                Assert.Equal((Sextant.OnHamiltonian.Breaks, Sextant.OnDissipator.Keeps, Sextant.Verdict.Neither), (pairing.OnH, pairing.OnD, pairing.Verdict));
            }
    }

    // F125's two traps as controls, on the cube: the eta ladder without its strings keeps the height, so the light
    // commutes, and loses the Hamiltonian part, where the strings place the Z a hop leaves behind; the single mode
    // c_k† rho c_k (integer modes: (1, 1) and (1, -1) at N = 2, (1, 0, -1) at N = 3 on the uniform chain, where
    // [H, c_k†] = E_k c_k†) keeps the Hamiltonian part and loses the dissipator, its cross terms one-sided at two
    // sites, moving the disagreement set. At N = 2 the string-free ladder does more than lose: its one string is Z_0,
    // a sublattice gauge in the sense of PROOF_K_PARTNERSHIP Lemma 1 (KHK = −H; its K counts the sites from 0 and is
    // diag(1, −1) on one excitation, so Z_0 is the other sublattice's gauge, −K there; either reverses the Hamiltonian
    // part and keeps the light, Lemma 4's first display from Lemmas 1 and 5, and ChiralKClaim), and the term at site 0
    // is even under Z_0 (σ⁺ and σ⁻ are both odd), so the string-free map is Φ followed by the gauge, Ad_{Z_0}∘Φ, and
    // FLIPS the Hamiltonian part exactly, a copy between the chain and the chain with its bond reversed. From N = 3
    // on no product of one-site diagonal turns does it: the l-th string-free term, followed by such a turn, is the
    // l-th term of Φ exactly when the turn is Z on every site below l, the identity on every site above it and
    // anything diagonal at l (σ⁺ρσ⁻ is invariant under any diagonal turn at its site); at N = 2 the two conditions
    // meet in Z_0, from N = 3 on the site 1 must be the identity for l = 0 and Z for l = 2.
    [Fact]
    public void The_Two_Traps_Of_F125_Fail_Where_The_Cube_Says()
    {
        foreach (int n in new[] { 2, 3, 4, 5 })
        {
            var chain = XY(n, Uniform(n - 1, 1), Uniform(n, 1));
            var bare = Read(chain, chain, SandwichSum.EtaLadder(n, strings: false));
            Assert.Equal((n == 2 ? Sextant.OnHamiltonian.Flips : Sextant.OnHamiltonian.Breaks, Sextant.OnDissipator.Keeps, Sextant.Verdict.Neither), (bare.OnH, bare.OnD, bare.Verdict));
            var turned = Read(chain, XY(n, Uniform(n - 1, -1), Uniform(n, 1)), SandwichSum.EtaLadder(n, strings: false));
            if (n == 2) AssertCopy(turned); else Assert.Equal(Sextant.Verdict.Neither, turned.Verdict);
        }
        foreach (var (n, phi) in new[] { (2, new long[] { 1, 1 }), (2, new long[] { 1, -1 }), (3, new long[] { 1, 0, -1 }) })
        {
            var chain = XY(n, Uniform(n - 1, 1), Uniform(n, 1));
            var mode = Read(chain, chain, SandwichSum.SingleMode(n, phi));
            Assert.Equal((Sextant.OnHamiltonian.Keeps, Sextant.OnDissipator.Breaks, Sextant.Verdict.Neither), (mode.OnH, mode.OnD, mode.Verdict));
            // the mode's H-part is exact only for an eigenmode: (1, 1, 1) at N = 3 is none, and both parts break
            if (n == 3)
            {
                var notAMode = Read(chain, chain, SandwichSum.SingleMode(3, new long[] { 1, 1, 1 }));
                Assert.Equal((Sextant.OnHamiltonian.Breaks, Sextant.OnDissipator.Breaks), (notAMode.OnH, notAMode.OnD));
            }
        }
    }

    // F142's Lemma 2.3 read exactly: S+ and its adjoint S- intertwine the dephased XY chain at any bond profile (bond
    // disorder survives) and any rate profile, and lose the Hamiltonian part, the dissipator kept, with a Z field on
    // one site, with a hop between sites of one parity (XZX + YZY, the second-neighbour hop c_0† c_2 + h.c.), with
    // the stagger dropped (F125's §6 map V), and with the ZZ term (F171's scope: the ZZ term loses the spin ladder; Φ's
    // loss is F125's Δ-tightness).
    [Fact]
    public void S_Plus_Intertwines_Exactly_When_The_Hops_Join_Opposite_Parities()
    {
        var rnd = new Random(142);
        foreach (int n in new[] { 2, 3, 4, 5 })
            for (int draw = 0; draw < 3; draw++)
            {
                var bonds = Draw(rnd, n - 1, -3, 3);
                var rates = Draw(rnd, n, 0, 3);
                if (bonds.All(b => b == 0)) bonds[0] = 1;
                if (rates.All(g => g == 0)) rates[0] = 1;
                var chain = XY(n, bonds, rates);
                AssertCopy(Read(chain, chain, SandwichSum.SpinLadder(n)));
                AssertCopy(Read(chain, chain, SandwichSum.SpinLadder(n).Adjoint));
                var unstaggered = Read(chain, chain, SandwichSum.SpinLadder(n, staggered: false));
                Assert.Equal((Sextant.OnHamiltonian.Breaks, Sextant.OnDissipator.Keeps, Sextant.Verdict.Neither), (unstaggered.OnH, unstaggered.OnD, unstaggered.Verdict));
                var fields = new long[n]; fields[rnd.Next(n)] = 1;
                var potential = Read(XY(n, bonds, rates, fields: fields), XY(n, bonds, rates, fields: fields), SandwichSum.SpinLadder(n));
                Assert.Equal((Sextant.OnHamiltonian.Breaks, Sextant.OnDissipator.Keeps, Sextant.Verdict.Neither), (potential.OnH, potential.OnD, potential.Verdict));
                if (n >= 3)
                {
                    var second = new[] { (HTerm(n,(0, 'X'), (1, 'Z'), (2, 'X')), 1L), (HTerm(n,(0, 'Y'), (1, 'Z'), (2, 'Y')), 1L) };
                    var sameParity = XY(n, bonds, rates, more: second);
                    var hop = Read(sameParity, sameParity, SandwichSum.SpinLadder(n));
                    Assert.Equal((Sextant.OnHamiltonian.Breaks, Sextant.OnDissipator.Keeps, Sextant.Verdict.Neither), (hop.OnH, hop.OnD, hop.Verdict));
                    AssertCopy(Read(sameParity, sameParity, SandwichSum.EtaLadder(n)));                // any quadratic conserving H serves Phi
                }
                var heisenberg = XY(n, bonds, rates, zz: 1);
                var lost = Read(heisenberg, heisenberg, SandwichSum.SpinLadder(n));
                Assert.Equal((Sextant.OnHamiltonian.Breaks, Sextant.OnDissipator.Keeps, Sextant.Verdict.Neither), (lost.OnH, lost.OnD, lost.Verdict));
            }
    }

    // F125 §5 on the ring, by height parity. A string's X and Y letters are where ket and bra disagree, so its height
    // k_Z has the parity of their particle numbers' difference (F61's n_XY parity, the weak symmetry Z^⊗N·ρ·Z^⊗N;
    // PROOF_CODIM1_BY_ADDITIVITY §5 reads it on the cube), and the blocks whose two sides see one Jordan-Wigner flux
    // are the strings of even height. Phi needs one hopping matrix on both sides, h_ket = h_bra, and S+ needs
    // Σ h_ket Σ = −h_bra (PROOF_GREY_GLASS_LADDER_IMAGES §7); each ladder steps the flux on both sides, so both land
    // on the ring with its wrap bond reversed, whose image block sees the fluxes in the source's order. Phi therefore
    // carries the ring onto the reversed ring on the even heights and Phi² back onto the ring, at every N, the odd
    // heights lost as a class (F125 keeps one edge pair at even N inside it, which a reading of the whole class does
    // not see). For S+ the parity of k that serves it is the ring's: on the even ring the wrap bond joins sites of
    // opposite parity and Σ h_f Σ = −h_f, so S+ reaches the reversed ring on the even heights; on the odd ring the
    // wrap bond joins two even sites and Σ h_f Σ = −h_{−f}, so S+ reaches the reversed ring on the ODD heights, where
    // the two sides see opposite fluxes, and neither ring on the even ones (F171 §7's "each ladder works on its own
    // parity of k", read here as a class and not from a corner). With the wrap bond zero the ring is the chain and
    // everything holds on every height. N = 3, 4 and 5, uniform light, the wrap bond at drawn couplings.
    [Fact]
    public void On_The_Ring_Each_Ladder_Reaches_The_Ring_With_Its_Wrap_Bond_Reversed_On_Its_Own_Parity_Of_The_Height()
    {
        var rnd = new Random(125);
        foreach (int n in new[] { 3, 4, 5 })
            for (int draw = 0; draw < 2; draw++)
            {
                var bonds = Draw(rnd, n - 1, 1, 3);
                long wrap = rnd.Next(1, 4) * (draw == 0 ? 1 : -1);
                var rates = Uniform(n, 1);
                var ring = XY(n, bonds, rates, wrap: wrap);
                var reversed = XY(n, bonds, rates, wrap: -wrap);
                var open = XY(n, bonds, rates, wrap: 0);
                var phi = SandwichSum.EtaLadder(n);
                AssertCopy(Read(ring, reversed, phi, parity: 0));
                AssertCopy(Read(reversed, ring, phi, parity: 0));
                AssertCopy(Read(ring, ring, phi.Then(phi, "Phi squared"), parity: 0));
                Assert.Equal(Sextant.Verdict.Neither, Read(ring, ring, phi, parity: 0).Verdict);
                Assert.Equal(Sextant.Verdict.Neither, Read(ring, reversed, phi, parity: 1).Verdict);
                Assert.Equal(Sextant.Verdict.Neither, Read(ring, ring, phi, parity: 1).Verdict);
                Assert.Equal(Sextant.Verdict.Neither, Read(ring, reversed, phi).Verdict);
                Assert.Equal(Sextant.OnDissipator.Keeps, Read(ring, reversed, phi, parity: 1).OnD);           // the light is blind to the wrap bond
                var sPlus = SandwichSum.SpinLadder(n);
                if (n % 2 == 0)
                {
                    AssertCopy(Read(ring, reversed, sPlus, parity: 0));
                    AssertCopy(Read(ring, reversed, sPlus.Adjoint, parity: 0));
                    Assert.Equal(Sextant.Verdict.Neither, Read(ring, ring, sPlus, parity: 0).Verdict);
                    Assert.Equal(Sextant.Verdict.Neither, Read(ring, reversed, sPlus, parity: 1).Verdict);
                    Assert.Equal(Sextant.Verdict.Neither, Read(ring, ring, sPlus, parity: 1).Verdict);
                }
                else
                {
                    Assert.Equal(Sextant.Verdict.Neither, Read(ring, reversed, sPlus, parity: 0).Verdict);
                    Assert.Equal(Sextant.Verdict.Neither, Read(ring, ring, sPlus, parity: 0).Verdict);
                    Assert.Equal(Sextant.OnDissipator.Keeps, Read(ring, ring, sPlus, parity: 0).OnD);
                    // the odd heights: the two sides see opposite fluxes, Σ h_f Σ = −h_{−f} is S+'s condition, and the
                    // image block of the reversed ring sees them in the source's order: the reversed ring is reached
                    AssertCopy(Read(ring, reversed, sPlus, parity: 1));
                    AssertCopy(Read(ring, reversed, sPlus.Adjoint, parity: 1));
                    Assert.Equal(Sextant.Verdict.Neither, Read(ring, ring, sPlus, parity: 1).Verdict);
                    Assert.Equal(Sextant.Verdict.Neither, Read(ring, reversed, sPlus).Verdict);
                }
                AssertCopy(Read(open, open, phi));
                AssertCopy(Read(open, open, sPlus));
            }
    }

    // The composition with R, the shift rho -> rho X^N, which reflects the dephased XY chain (§9): a copy followed by
    // a reflection is a reflection, so M = R∘Phi carries L onto -L† - c with c = 2Σγ_l, read off the identity and not
    // the trace, since M has a kernel; by F166 it is S⁺∘R as well, the reflection followed by the spin ladder (F125's
    // containment corollary and the shared-field page hold it on blocks and eigenvectors, W composed with a fold at
    // the mirrored eigenvalue; here it is one map on the full generator). With the transpose first, Π_Z∘Phi carries L
    // onto -L - c. The rank is Phi's. The spectrum pairs fully through R; the kernel limits what M transports.
    [Fact]
    public void Phi_Then_R_Reflects_With_A_Kernel()
    {
        var rnd = new Random(121);
        foreach (int n in new[] { 2, 3, 4 })
        {
            var bonds = Draw(rnd, n - 1, 1, 3);
            var rates = Draw(rnd, n, 1, 3);
            var chain = XY(n, bonds, rates);
            var phi = SandwichSum.EtaLadder(n);
            var viaR = Read(chain, chain, phi.Then(R(n)));
            Assert.Equal((Sextant.OnHamiltonian.Keeps, Sextant.OnDissipator.Reflects, Sextant.Verdict.Reflection), (viaR.OnH, viaR.OnD, viaR.Verdict));
            Assert.Equal($"{2 * rates.Sum()}", viaR.C);
            Assert.Equal(Read(chain, chain, phi).Rank, viaR.Rank);
            var moves1 = LetterCube.TableMoves(1);
            var piZ = LetterCube.Move.Uniform(moves1.Single(m => m.Name == "the transpose").Sites[0].Then(moves1.Single(m => m.Name.StartsWith("R,")).Sites[0], "Pi_Z"), n);
            var viaPi = Read(chain, chain, phi.Then(piZ));
            Assert.Equal((Sextant.OnHamiltonian.Flips, Sextant.OnDissipator.Reflects, Sextant.Verdict.PiReflection), (viaPi.OnH, viaPi.OnD, viaPi.Verdict));
            Assert.Equal($"{2 * rates.Sum()}", viaPi.C);
        }
    }

    static string Plain(double v) => (v == 0 ? 0 : v).ToString(System.Globalization.CultureInfo.InvariantCulture);   // no "-0"

    // the dense verdict: M L_from against Λ M on the sources, Λ = L_to, -L_to† - c, -L_to - c, c read off the first
    // nonzero entry of M and checked everywhere (the entries are small Gaussian integers, exact in doubles)
    static (Sextant.Verdict Verdict, string C) DenseVerdict(Sextant.Setting from, Sextant.Setting to, DenseSum m, int? parity)
    {
        int n = from.N;
        var lFrom = DenseGenerator(from);
        var lTo = DenseGenerator(to);
        var sources = LetterCube.Strings(n).Where(s => parity is null || (BitOperations.PopCount(s.X) & 1) == parity).Select(s => DenseString(s.ToString(n))).ToList();
        bool copy = true, refl = true, pi = true;
        Complex? kRefl = null, kPi = null;
        foreach (var src in sources)
        {
            var ms = m.Apply(src);
            var a = m.Apply(lFrom(src, false));
            var b = lTo(ms, false);
            var bd = lTo(ms, true);
            int d = 1 << n;
            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                {
                    if (a[i, j] != b[i, j]) copy = false;
                    var rRefl = a[i, j] + bd[i, j];               // = -c · M
                    var rPi = a[i, j] + b[i, j];
                    if (ms[i, j] == Complex.Zero)
                    {
                        if (rRefl != Complex.Zero) refl = false;
                        if (rPi != Complex.Zero) pi = false;
                        continue;
                    }
                    kRefl ??= rRefl / ms[i, j];
                    kPi ??= rPi / ms[i, j];
                    if (rRefl != kRefl.Value * ms[i, j]) refl = false;
                    if (rPi != kPi.Value * ms[i, j]) pi = false;
                }
        }
        if (kRefl is null || kPi is null) return (Sextant.Verdict.Neither, "");
        if (copy) return (Sextant.Verdict.Copy, "");
        if (refl && kRefl.Value.Imaginary == 0) return (Sextant.Verdict.Reflection, Plain(-kRefl.Value.Real));
        if (pi && kPi.Value.Imaginary == 0) return (Sextant.Verdict.PiReflection, Plain(-kPi.Value.Real));
        return (Sextant.Verdict.Neither, "");
    }

    // The relation against the identity computed densely at N = 2 and 3: the same sources, the same three forms, c
    // read the same way, on the chain (Phi, S+, the two traps, Phi then R), on the ring by parity (Phi onto the
    // reversed ring, Phi onto the ring, S+ on the odd ring) and on Heisenberg.
    [Fact]
    public void The_Intertwining_Agrees_With_The_Dense_Generator()
    {
        var rnd = new Random(3);
        int seen = 0;
        foreach (int n in new[] { 2, 3 })
        {
            var bonds = Draw(rnd, n - 1, 1, 3);
            var rates = Draw(rnd, n, 1, 2);
            var chain = XY(n, bonds, rates);
            var heisenberg = XY(n, bonds, rates, zz: 1);
            var r = R(n);
            var rDense = DenseString(new string('X', n));
            var cases = new List<(Sextant.Setting From, Sextant.Setting To, SandwichSum Sum, DenseSum Dense, int? Parity)>
            {
                (chain, chain, SandwichSum.EtaLadder(n), DensePhi(n), null),
                (chain, chain, SandwichSum.SpinLadder(n), DenseSPlus(n), null),
                (chain, chain, SandwichSum.EtaLadder(n, strings: false), DensePhi(n, strings: false), null),
                (chain, chain, SandwichSum.SingleMode(n, n == 2 ? new long[] { 1, 1 } : new long[] { 1, 0, -1 }), DenseSingleMode(n, n == 2 ? new long[] { 1, 1 } : new long[] { 1, 0, -1 }), null),
                (chain, chain, SandwichSum.EtaLadder(n).Then(r), new DenseSum(n, DensePhi(n).Terms.Select(t => (t.W, t.A, Mul(t.B, rDense))).ToList()), null),
                (heisenberg, heisenberg, SandwichSum.EtaLadder(n), DensePhi(n), null),
                (heisenberg, heisenberg, SandwichSum.SpinLadder(n), DenseSPlus(n), 0),
            };
            if (n == 3)
            {
                var ring = XY(n, bonds, rates, wrap: 2);
                var reversed = XY(n, bonds, rates, wrap: -2);
                cases.Add((ring, reversed, SandwichSum.EtaLadder(n), DensePhi(n), 0));
                cases.Add((ring, ring, SandwichSum.EtaLadder(n), DensePhi(n), 0));
                cases.Add((ring, reversed, SandwichSum.EtaLadder(n), DensePhi(n), 1));
                cases.Add((ring, reversed, SandwichSum.SpinLadder(n), DenseSPlus(n), 0));
                cases.Add((ring, ring, SandwichSum.SpinLadder(n), DenseSPlus(n), 0));
                cases.Add((ring, reversed, SandwichSum.SpinLadder(n), DenseSPlus(n), 1));
            }
            foreach (var (from, to, sum, dense, parity) in cases)
            {
                var reading = Read(from, to, sum, parity);
                var (verdict, c) = DenseVerdict(from, to, dense, parity);
                Assert.True(reading.Verdict == verdict, $"{sum.Name}, parity {parity}, N = {n}: the sextant reads {reading.Verdict}, the dense route {verdict}");
                if (verdict != Sextant.Verdict.Neither)                             // the dense route reads c for a form that holds
                    Assert.True(reading.C == c, $"{sum.Name}, parity {parity}, N = {n}: c = {reading.C} on the letters, {c} densely");
                if (verdict != Sextant.Verdict.Neither) seen++;
            }
        }
        Assert.True(seen >= 6);                                           // the positive branches are checked too
    }

    // Composition and the scale: Phi then Phi on the letters is the dense Phi applied twice, the scale 2^N kept, and a
    // move's pieces are its units times two, so FromMove(R) applied to a string is R's image times 2^N.
    [Fact]
    public void Composition_Keeps_The_Scale()
    {
        const int n = 3;
        var phi = SandwichSum.EtaLadder(n);
        var twice = phi.Then(phi);
        var dense = DensePhi(n);
        foreach (var s in LetterCube.Strings(n))
        {
            AssertSameExpansion(twice.Apply(s), Expand(dense.Apply(dense.Apply(DenseString(s.ToString(n)))), n), $"Phi squared on {s.ToString(n)}");
            var (image, phase) = R(n).Apply(s);
            AssertSameMap(new Dictionary<PauliString, LetterCube.Gauss> { [image] = (1L << n) * phase }, SandwichSum.FromMove(R(n)).Apply(s), $"R on {s.ToString(n)}");
        }
        Assert.Throws<ArgumentException>(() => phi.Then(SandwichSum.EtaLadder(2)));
        Assert.Throws<ArgumentOutOfRangeException>(() => Sextant.Intertwines(XY(2, new long[] { 1 }, new long[] { 1, 1 }), XY(2, new long[] { 1 }, new long[] { 1, 1 }), SandwichSum.EtaLadder(2), 2));
        Assert.Throws<ArgumentException>(() => Sextant.Intertwines(XY(2, new long[] { 1 }, new long[] { 1, 1 }), XY(2, new long[] { 1 }, new long[] { 1, 1 }), SandwichSum.EtaLadder(3)));
        // the class: a field along X carries a source out of its height parity, and the reading says so
        var fielded = XY(3, new long[] { 1, 1 }, new long[] { 1, 1, 1 }, more: new[] { (HTerm(3, (0, 'X')), 1L) });
        Assert.False(Sextant.Intertwines(fielded, fielded, SandwichSum.EtaLadder(3), 0).ClassKept);
        var ring3 = XY(3, new long[] { 1, 1 }, new long[] { 1, 1, 1 }, wrap: 1);
        Assert.True(Sextant.Intertwines(ring3, XY(3, new long[] { 1, 1 }, new long[] { 1, 1, 1 }, wrap: -1), SandwichSum.EtaLadder(3), 0).ClassKept);
    }
}
