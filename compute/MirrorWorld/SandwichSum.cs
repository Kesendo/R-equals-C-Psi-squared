namespace MirrorWorld;

using static MirrorWorld.LetterCube;

// A sandwich sum (built 2026-11-07 beside the sextant): a map on the operator space that is a SUM of sandwiches, rho -> sum_t w_t · A_t rho^(T) B_t, every A_t and B_t a product of
// one-site 2x2 Gaussian factors, held exactly on the letters. It is the move kind the cube's moves and the theodolite's
// sandwich do not reach:
// the cube's moves (LetterCube.SiteMove, Move) are monomial, one letter to a unit times one letter, the theodolite's
// sandwich is one product of turns, and neither takes a sum over the sites or a factor of rank one. The two ladders
// of the doubled picture are exactly that: the eta ladder Phi(rho) = sum_l c_l^dagger rho c_l of F125 (the spectator
// map W; F142 writes it Phi, docs/proofs/PROOF_FROZEN_BAND_SO4.md) and the spin ladder S+(rho) = sum_l (-1)^l
// c_l^dagger rho c_l^dagger of F142, with c_l = Z_0 ... Z_{l-1} sigma^-_l the Jordan-Wigner fermion of site l, its
// string on the sites before it (the glossary's Fermions section; sites counted from 0 here, as in
// RCPsiSquared.Diagnostics' SpectatorIntertwiner, which builds the same two maps on cells, so the stagger gives site 0
// the sign +1 where the proofs, counting from 1, give it -1: S+ differs between the two by a global sign, Phi does
// not). What the sextant asks of a sandwich sum is Sextant.Intertwines; this file holds the map alone.
//
// On the letters a term is a product over the sites of one-site pieces, each piece the 4x4 matrix of rho -> A P B
// in the letters I, X, Y, Z: piece[c, b] = Tr(P_c A P_b B) / 2, a Gaussian half-integer, so every piece is held
// TIMES TWO and a term over N sites, and the whole sum, times 2^N (the sextant's generator is held times four, and
// the identity M L = L' M is homogeneous in M, so the scale never enters a verdict). Composition divides the scale
// back out, exactly: by the completeness of the letters, sum_d P_d Tr(P_d X) = 2X, the product of two pieces times two
// is twice the composite's piece times two.
//
// Read on the cube, the reason the ladders commute with the light, which
// docs/THE_ONE_SQUARE.md §10 already places on the doubled site (the two doublets, raised by sigma^+ rho sigma^- and
// by sigma^+ rho sigma^+, I and Z against X and Y): the piece of Phi at its site, sigma^+ rho sigma^-, sends I and Z
// to (I - Z)/2 and X and Y to zero, so it keeps the site's quiet edge and empties its far one; the piece of S+,
// sigma^+ rho sigma^+, sends X to (X - iY)/2 and Y to -i(X - iY)/2 and empties the quiet edge; and the string Z_{<l}
// on both sides is a sign on the far letters below l. No piece moves a letter between a site's two edges, so every
// term keeps every site's edge, the disagreement set of a cell and not only its size, the height k_Z; Z dephasing
// charges a string -2·sum of gamma_l over its far letters, a function of that set, so it commutes with each term and
// with the sum at every rate profile (PROOF_FROZEN_BAND_SO4's Lemma 2.1, "neither ladder moves the disagreement set";
// F125's Lemma 2, "the sign squares away"), and at one rate the charge is one number per height plane. The
// single-mode trap of F125 fails for the same reason read backwards: c_k^dagger rho c_k for a mode c_k = sum_l phi_l
// c_l holds cross terms c_l^dagger rho c_m with l != m, one-sided at two sites, and a one-sided letter multiplication
// steps between the edges (sigma^+ · I = sigma^+ lies on the far edge, sigma^+ · X = (I - Z)/2 on the quiet one), so
// a cross term flips the edge at both sites, l and m, and the height moves by 0 or ±2: at one rate the ±2 images are
// what the light tells apart, at a profile the moved set as well. The Hamiltonian part is where
// the SUM works and not the term: a hop at a bond, acting on the one-sided factor of the term at one of its sites,
// leaves a Z behind at the site it empties (the commutator of sigma^+ with sigma^-, -Z here), and that Z is the string's letter
// of the neighbouring fermion, so the two cross terms a bond makes, from the left factor of one term and the right
// factor of the other, are one operator c_m^dagger rho c_l with opposite signs and cancel; without the strings the Z
// sits on the left in one and on the right in the other, and F125's first trap fails the Hamiltonian part. The stagger is the spin ladder's version of the same cancellation:
// its cross terms carry (-1)^l h_{ml} + (-1)^m h_{ml}, which vanish exactly when h joins only sites of opposite
// parity (F142's Lemma 2.3, Sigma h Sigma = -h, which for nearest-neighbour hopping PROOF_K_PARTNERSHIP's Lemma 1
// proves as KHK = -H, its K counted from 0). And
// the two ladders are one map seen through the bra's flip, F166's identity Pi_Z W Pi_Z^-1 = S+
// (experiments/THE_MIRROR_UNDER_A_SHARED_FIELD.md: the transpose leaves W fixed, and R, a particle-hole
// transformation of the bra species, carries W onto S+; CorrelatedZDephasingPriceClaim calls it a staggered Shiba
// map): on the letters the shift R,
// rho -> rho X^N, turns the right-hand sigma^- into sigma^+ and every Z of the right-hand string into -Z, X Z X = -Z,
// one sign per site below l, which is the stagger; so R Phi R = S+ with the sites counted from 0, read exactly in
// the tests. R here is the cube's shift by X^N on the right (LetterCube's and MirrorGroup's R), not F125's R, the
// reflection of the sites, which grades W by (-1)^(p+q).
//
// Words, fenced at the door: a PIECE is the one-site 4x4 matrix of a factor pair, not GammaFold's PieceCount (the
// pieces of a descended map), CollisionGap's pieces (Corollary G's) or LetterCube's Step and Arrow; a TERM is one
// sandwich of the sum, not a term of H (Setting.H holds those); the SITE SUM of F125 ("the site sum matters") is
// Phi's sum over one site on both sides, a member here, where the type is named for the wider thing, any sum of
// sandwiches, the single mode and a lone move among them; the STRING of a fermion is its Jordan-Wigner factor
// Z_{<l}, a Pauli string in the letters, and PauliString is the type (not "the Klein string": Klein already names
// PauliMode's parity cell, Mirror's leg, ParameterKlein, Lattice's four watchings and the antilinear triangle); the
// LADDERS are F142's two, as in GreyRoom, and not Divisor's ladder (how long a frozen mode holds) or CollisionGap's
// (the comb rung); the STAGGER is the sublattice sign (-1)^l with the sites counted from 0, the Sigma of F170 (F171's
// proof counts from 1, so its Sigma takes the other sublattice). A sandwich sum is machinery like PauliString: it
// states no physics on its own, and the readings live in the sextant and its tests.
public sealed class SandwichSum
{
    const string Letters = "IXYZ";

    // one sandwich of the sum: its weight and one piece per site, each times two
    public sealed record Term(Gauss Weight, Gauss[][,] PiecesTimesTwo);

    public string Name { get; }
    public int N { get; }
    public IReadOnlyList<Term> Terms { get; }

    public SandwichSum(string name, int n, IEnumerable<Term> terms)
    {
        if (n < 1 || n > PauliString.MaxSites)
            throw new ArgumentOutOfRangeException(nameof(n), n, $"1 to {PauliString.MaxSites} sites");
        Name = name;
        N = n;
        Terms = terms.ToList();
        foreach (var t in Terms)
            if (t.PiecesTimesTwo.Length != n) throw new ArgumentException("a term carries one piece per site", nameof(terms));
    }

    /// <summary>The scale Apply returns at: 2^N times the map.</summary>
    public long Scale => 1L << N;

    static Gauss[,] Mul2(Gauss[,] a, Gauss[,] b)
    {
        var c = new Gauss[2, 2];
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                c[i, j] = a[i, 0] * b[0, j] + a[i, 1] * b[1, j];
        return c;
    }

    // The piece of rho -> A rho^(T) B on one site, times two: piece[c, b] = Tr(P_c A P_b^(T) B), letters I, X, Y, Z.
    public static Gauss[,] PieceTimesTwo(Gauss[,] a, Gauss[,] b, bool transpose = false)
    {
        var piece = new Gauss[4, 4];
        for (int k = 0; k < 4; k++)
        {
            var p = LetterMatrix(Letters[k]);
            if (transpose) p = new[,] { { p[0, 0], p[1, 0] }, { p[0, 1], p[1, 1] } };
            var m = Mul2(a, Mul2(p, b));
            for (int c = 0; c < 4; c++)
            {
                var prod = Mul2(LetterMatrix(Letters[c]), m);
                piece[c, k] = prod[0, 0] + prod[1, 1];
            }
        }
        return piece;
    }

    public static Gauss[,] IdentityPieceTimesTwo()
    {
        var piece = new Gauss[4, 4];
        for (int b = 0; b < 4; b++) piece[b, b] = new Gauss(2, 0);
        return piece;
    }

    // the piece of a monomial one-site map: 2 · its unit at [image, letter]
    public static Gauss[,] PieceTimesTwo(SiteMove m)
    {
        var piece = new Gauss[4, 4];
        for (int b = 0; b < 4; b++) piece[m.Image[b], b] = 2 * m.Phase[b];
        return piece;
    }

    /// <summary>One sandwich w · A rho^(T) B with A = ⊗ factors[l].A and B = ⊗ factors[l].B.</summary>
    public static Term Sandwich(Gauss weight, IReadOnlyList<(Gauss[,] A, Gauss[,] B)> factors, bool transpose = false)
        => new(weight, factors.Select(f => PieceTimesTwo(f.A, f.B, transpose)).ToArray());

    public static SandwichSum FromMove(Move m)
        => new(m.Name, m.Sites.Length, new[] { new Term(Gauss.One, m.Sites.Select(PieceTimesTwo).ToArray()) });

    static Gauss[,] I2 => LetterMatrix('I');
    static Gauss[,] Z2 => LetterMatrix('Z');

    // the left and right factor of site m in the sandwich A rho B with A = Z_{<l} a_l and B = Z_{<l'} b_{l'}
    static (Gauss[,], Gauss[,]) Factor(int m, int l, Gauss[,] a, int lPrime, Gauss[,] b, bool strings)
        => (m == l ? a : m < l && strings ? Z2 : I2, m == lPrime ? b : m < lPrime && strings ? Z2 : I2);

    /// <summary>Phi, the eta ladder: sum_l c_l^dagger rho c_l (F125's W), the strings included; without them F125's
    /// first trap, sum_l sigma^+_l rho sigma^-_l, a different map.</summary>
    public static SandwichSum EtaLadder(int n, bool strings = true)
        => new(strings ? "Phi, the eta ladder" : "the eta ladder without its strings", n,
               Enumerable.Range(0, n).Select(l => Sandwich(Gauss.One,
                   Enumerable.Range(0, n).Select(m => Factor(m, l, Raising, l, Lowering, strings)).ToList())));

    /// <summary>S+, the spin ladder: sum_l (-1)^l c_l^dagger rho c_l^dagger (F142), or unstaggered, every weight one.</summary>
    public static SandwichSum SpinLadder(int n, bool staggered = true)
        => new(staggered ? "S+, the spin ladder" : "the spin ladder without its stagger", n,
               Enumerable.Range(0, n).Select(l => Sandwich(staggered && l % 2 == 1 ? -Gauss.One : Gauss.One,
                   Enumerable.Range(0, n).Select(m => Factor(m, l, Raising, l, Raising, true)).ToList())));

    /// <summary>F125's second trap, the single-mode spectator c^dagger rho c with c^dagger = sum_l phi_l c_l^dagger
    /// for integer amplitudes phi: the cross terms c_l^dagger rho c_m with l != m are what F125's site sum, Phi, lacks.</summary>
    public static SandwichSum SingleMode(int n, IReadOnlyList<long> amplitude)
    {
        if (amplitude.Count != n) throw new ArgumentException("one amplitude per site", nameof(amplitude));
        var terms = new List<Term>();
        for (int l = 0; l < n; l++)
            for (int lp = 0; lp < n; lp++)
            {
                long w = amplitude[l] * amplitude[lp];
                if (w == 0) continue;
                terms.Add(Sandwich(new Gauss(w, 0), Enumerable.Range(0, n).Select(m => Factor(m, l, Raising, lp, Lowering, true)).ToList()));
            }
        return new SandwichSum($"the single mode ({string.Join(", ", amplitude)})", n, terms);
    }

    /// <summary>The adjoint map under the trace inner product: every piece conjugate-transposed and every weight
    /// conjugated (the letters are orthogonal under the trace, all of one norm), which for a term A rho B is
    /// sigma -> conj(w) · A^dagger sigma B^dagger and for a transposed term A rho^T B is
    /// sigma -> conj(w) · conj(B) sigma^T conj(A);
    /// so Phi's adjoint is W^dagger = sum_l c_l sigma c_l^dagger and S+'s is S-.</summary>
    public SandwichSum Adjoint
    {
        get
        {
            var terms = Terms.Select(t => new Term(t.Weight.Conj, t.PiecesTimesTwo.Select(p =>
            {
                var q = new Gauss[4, 4];
                for (int c = 0; c < 4; c++) for (int b = 0; b < 4; b++) q[c, b] = p[b, c].Conj;
                return q;
            }).ToArray()));
            return new SandwichSum($"the adjoint of {Name}", N, terms);
        }
    }

    /// <summary>The composite next ∘ this, one term per pair of terms, the pieces multiplied and the scale divided out
    /// exactly (the completeness of the letters makes the product of two pieces times two twice a piece times two).</summary>
    public SandwichSum Then(SandwichSum next, string? name = null)
    {
        if (next.N != N) throw new ArgumentException("the two sums act on the same sites", nameof(next));
        var terms = new List<Term>();
        foreach (var t in Terms)
            foreach (var u in next.Terms)
            {
                var pieces = new Gauss[N][,];
                for (int l = 0; l < N; l++)
                {
                    var p = new Gauss[4, 4];
                    for (int c = 0; c < 4; c++)
                        for (int b = 0; b < 4; b++)
                        {
                            var v = Gauss.Zero;
                            for (int d = 0; d < 4; d++) v = v + u.PiecesTimesTwo[l][c, d] * t.PiecesTimesTwo[l][d, b];
                            if (v.Re % 2 != 0 || v.Im % 2 != 0)
                                throw new InvalidOperationException(
                                    "the product of two pieces times two is twice a piece times two; it is not, so a piece was not built from factors");
                            p[c, b] = new Gauss(v.Re / 2, v.Im / 2);
                        }
                    pieces[l] = p;
                }
                terms.Add(new Term(t.Weight * u.Weight, pieces));
            }
        return new SandwichSum(name ?? $"{Name}, then {next.Name}", N, terms);
    }

    public SandwichSum Then(Move m, string? name = null) => Then(FromMove(m), name);

    /// <summary>The image of a string, times 2^N: a sum of strings with Gaussian-integer coefficients, zeros dropped.</summary>
    public IReadOnlyDictionary<PauliString, Gauss> Apply(PauliString s)
    {
        var result = new Dictionary<PauliString, Gauss>();
        foreach (var t in Terms)
        {
            var partial = new Dictionary<PauliString, Gauss> { [s] = t.Weight };
            for (int l = 0; l < N && partial.Count > 0; l++)
            {
                int b = LetterIndex(s, l);
                var piece = t.PiecesTimesTwo[l];
                var next = new Dictionary<PauliString, Gauss>();
                for (int c = 0; c < 4; c++)
                {
                    var v = piece[c, b];
                    if (v.IsZero) continue;
                    foreach (var (str, coeff) in partial)
                    {
                        var image = WithLetter(str, l, c);
                        next[image] = next.TryGetValue(image, out var old) ? old + coeff * v : coeff * v;
                    }
                }
                partial = next;
            }
            foreach (var (str, coeff) in partial)
                result[str] = result.TryGetValue(str, out var old) ? old + coeff : coeff;
        }
        return result.Where(kv => !kv.Value.IsZero).ToDictionary(kv => kv.Key, kv => kv.Value);
    }
}
