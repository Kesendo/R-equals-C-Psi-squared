using System.Numerics;

namespace MirrorWorld;

// The mirror group (adopted 2026-07-04 from F118, docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md):
// the canonical palindromizer is not elementary -- it factors. With F = X^N, the ket reflection
// R(rho) = rho*F and the transpose D(rho) = rho^T generate, on the Pauli basis, a group of EIGHT
// signed permutations, the dihedral D4:
//
//   rotations        {1, Pi_Z = R o D, F = Pi_Z^2, Pi_Y = Pi_Z^3}   (F = the charge conjugation, the center)
//   diagonal mirrors {D = diag((-1)^nY), FD = diag((-1)^nZ)}        (the Klein cell signs PauliMode carries)
//   edge mirrors     {R: rho -> rho*F,  FR: rho -> F*rho}           (the one-sided reflections)
//
// Every member acts site-locally -- one letter rule, the same at every site, lifted by product --
// and every phase lives in {1, -1, i, -i}, so the whole group compares EXACTLY, no eigensolver,
// no tolerance. Pi_Z walks the April rule I -> X, X -> I, Y -> iZ, Z -> iY with no extra phase:
// the hard-won i falls out of Y^T = -Y meeting YX = -iZ.
//
// The palindrome splits along the generators: D flips the Hamiltonian commutator (D L_H D = -L_H)
// and fixes the dissipator; R fixes L_H and reflects the dissipator, carrying the entire constant
// (R L_diss R = -L_diss - 2*sigma, sigma = sum_l gamma_l -- the same constant Mirror pays as its
// price: this group is Mirror's fold lattice read at the operator level instead of the block level).
// And the polarity cube's three axes are characters: bit_a of Ad_{Z^N}, bit_b of Ad_{X^N} = F,
// y_par of the transpose D -- the truly cell is the joint-fixed cell of the two diagonal mirrors.
//
// Deliberately outside (named open in F118): the letter group S3 (adjoining it would assemble
// S3 x| D4), K1, the golden router W (F116), and F71's bond mirror.
public sealed class MirrorGroup : GameObject
{
    public int N { get; }

    public MirrorGroup(World world, int n) : base(world) => N = n;

    // left: what the group itself produces.
    public override IReadOnlyList<string> Own => new[] { "members", "split", "cube" };

    // ---- a member: a site-local signed permutation of the Pauli letters (I X Y Z), possibly
    // antilinear (conjugates coefficients; the F119 antilinear double needs it). ----
    public sealed class Member : IEquatable<Member>
    {
        internal readonly int[] To;        // letter index (I=0, X=1, Y=2, Z=3) -> letter index
        internal readonly Complex[] Ph;    // letter index -> phase in {1, -1, i, -i}
        public bool Antilinear { get; }
        public string Name { get; }

        internal Member(string name, int[] to, Complex[] ph, bool antilinear)
        {
            Name = name; To = to; Ph = ph; Antilinear = antilinear;
        }

        public (char Letter, Complex Phase) ApplySite(char letter)
        {
            int i = Idx(letter);
            return (Alphabet[To[i]], Ph[i]);
        }

        // lift to a Pauli string: per-site images, phase = the product of the per-site phases.
        public (char[] Letters, Complex Phase) Apply(char[] letters)
        {
            var image = new char[letters.Length];
            Complex phase = Complex.One;
            for (int s = 0; s < letters.Length; s++)
            {
                int i = Idx(letters[s]);
                image[s] = Alphabet[To[i]];
                phase *= Ph[i];
            }
            return (image, phase);
        }

        public bool Equals(Member? other)
        {
            if (other is null || Antilinear != other.Antilinear) return false;
            for (int i = 0; i < 4; i++)
                if (To[i] != other.To[i] || Ph[i] != other.Ph[i]) return false;
            return true;
        }

        public override bool Equals(object? obj) => Equals(obj as Member);

        public override int GetHashCode()
        {
            var h = new HashCode();
            h.Add(Antilinear);
            for (int i = 0; i < 4; i++) { h.Add(To[i]); h.Add(Ph[i].Real); h.Add(Ph[i].Imaginary); }
            return h.ToHashCode();
        }
    }

    static readonly char[] Alphabet = { 'I', 'X', 'Y', 'Z' };
    static int Idx(char c) => c switch
    {
        'I' => 0, 'X' => 1, 'Y' => 2, 'Z' => 3,
        _ => throw new ArgumentException($"not a Pauli letter: {c}")
    };

    static readonly Complex One = Complex.One, MinusOne = -Complex.One;
    static readonly Complex Im = Complex.ImaginaryOne, MinusIm = -Complex.ImaginaryOne;
    static readonly int[] Keep = { 0, 1, 2, 3 };     // letters untouched (diagonal members)
    static readonly int[] Swap = { 1, 0, 3, 2 };     // I <-> X, Y <-> Z (the F-multiplied members)

    // the eight (per-site tables, each verified against its operator form in EightFormsWorstResidual):
    public static readonly Member Identity = new("1", Keep, new[] { One, One, One, One }, false);
    public static readonly Member PiZ = new("Pi_Z", Swap, new[] { One, One, Im, Im }, false);         // rho -> rho^T * F
    public static readonly Member F = new("F", Keep, new[] { One, One, MinusOne, MinusOne }, false);  // rho -> F rho F
    public static readonly Member PiY = new("Pi_Y", Swap, new[] { One, One, MinusIm, MinusIm }, false); // rho -> F * rho^T
    public static readonly Member D = new("D", Keep, new[] { One, One, MinusOne, One }, false);       // rho -> rho^T
    public static readonly Member FD = new("FD", Keep, new[] { One, One, One, MinusOne }, false);     // rho -> F rho^T F
    public static readonly Member R = new("R", Swap, new[] { One, One, MinusIm, Im }, false);         // rho -> rho * F
    public static readonly Member FR = new("FR", Swap, new[] { One, One, Im, MinusIm }, false);       // rho -> F * rho

    // outside the eight but needed by the cube (bit_a's conjugation) and the double (the antilinear unit):
    public static readonly Member AdZ = new("AdZ", Keep, new[] { One, MinusOne, MinusOne, One }, false); // rho -> Z^N rho Z^N
    public static readonly Member K = new("K", Keep, new[] { One, One, One, One }, true);                // rho -> rho^dagger

    // a after b (b applied first). If a is antilinear it conjugates b's phase on the way through.
    public static Member Compose(Member a, Member b)
    {
        var to = new int[4];
        var ph = new Complex[4];
        for (int i = 0; i < 4; i++)
        {
            to[i] = a.To[b.To[i]];
            var carried = a.Antilinear ? Complex.Conjugate(b.Ph[i]) : b.Ph[i];
            ph[i] = carried * a.Ph[b.To[i]];
        }
        return new Member($"{a.Name} o {b.Name}", to, ph, a.Antilinear ^ b.Antilinear);
    }

    // brute closure from the generators (small groups; exact equality, so a plain list suffices).
    public static IReadOnlyList<Member> Closure(params Member[] generators)
    {
        var members = new List<Member> { Identity };
        bool grew = true;
        while (grew)
        {
            grew = false;
            foreach (var g in generators)
                foreach (var m in members.ToArray())
                {
                    var c = Compose(g, m);
                    if (!members.Any(x => x.Equals(c))) { members.Add(c); grew = true; }
                }
        }
        return members;
    }

    // ---- from-below verification: the eight Pauli-basis forms against their operator definitions,
    // on every Pauli string at N. A signed permutation compares exactly; worst |diff| is returned. ----
    public double EightFormsWorstResidual()
    {
        var f = Dense.PauliString(Enumerable.Repeat('X', N).ToArray());
        var forms = new (Member Member, Func<Complex[,], Complex[,]> Operator)[]
        {
            (Identity, rho => rho),
            (PiZ, rho => Dense.Mul(Dense.Transpose(rho), f)),
            (F,   rho => Dense.Mul(f, Dense.Mul(rho, f))),
            (PiY, rho => Dense.Mul(f, Dense.Transpose(rho))),
            (D,   rho => Dense.Transpose(rho)),
            (FD,  rho => Dense.Mul(f, Dense.Mul(Dense.Transpose(rho), f))),
            (R,   rho => Dense.Mul(rho, f)),
            (FR,  rho => Dense.Mul(f, rho)),
        };
        double worst = 0;
        foreach (var letters in Dense.AllStrings(N))
        {
            var sigma = Dense.PauliString(letters);
            foreach (var (member, op) in forms)
            {
                var (image, phase) = member.Apply(letters);
                var expected = Dense.Scale(phase, Dense.PauliString(image));
                worst = Math.Max(worst, Dense.MaxAbsDiff(op(sigma), expected));
            }
        }
        return worst;
    }

    // ---- the action on a full deterministic rho: Pi_Z applied in coefficient space (expand in the
    // Pauli basis, move every string by the signed permutation, resum) against rho^T * F -- and
    // against the wrong-sided F * rho^T, which must be rejected at O(1). ----
    public (double Right, double Wrong) PiZOnRho()
    {
        int d = 1 << N;
        var rho = Dense.TestMatrix(d, 3, 7);
        var viaGroup = new Complex[d, d];
        foreach (var letters in Dense.AllStrings(N))
        {
            var sigma = Dense.PauliString(letters);
            Complex c = Dense.Trace(Dense.Mul(sigma, rho)) / d;      // Tr(sigma rho)/2^N
            var (image, phase) = PiZ.Apply(letters);
            var img = Dense.PauliString(image);
            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                    viaGroup[i, j] += c * phase * img[i, j];
        }
        var f = Dense.PauliString(Enumerable.Repeat('X', N).ToArray());
        var rhoT = Dense.Transpose(rho);
        return (Dense.MaxAbsDiff(viaGroup, Dense.Mul(rhoT, f)),
                Dense.MaxAbsDiff(viaGroup, Dense.Mul(f, rhoT)));
    }

    // ---- the palindrome split, verified on the FULL Pauli basis (a superoperator identity checked
    // on a basis is checked, period): D L_H D = -L_H and D L_diss D = L_diss; R L_H R = L_H and
    // R L_diss R = -L_diss - 2*sigma. XXZ chain with anisotropy delta, site-dependent gamma. ----
    public (double DFlipsH, double DFixesDiss, double RFixesH, double RReflectsDiss)
        PalindromeSplitResiduals(double j, double delta, double[] gammas)
    {
        var h = Dense.XxzChain(N, j, delta);
        var f = Dense.PauliString(Enumerable.Repeat('X', N).ToArray());
        var zs = Enumerable.Range(0, N).Select(l => Dense.SiteZ(N, l)).ToArray();
        double sigmaTot = gammas.Sum();
        double dH = 0, dDiss = 0, rH = 0, rDiss = 0;
        foreach (var letters in Dense.AllStrings(N))
        {
            var sigma = Dense.PauliString(letters);
            var sigmaT = Dense.Transpose(sigma);
            var sigmaF = Dense.Mul(sigma, f);

            // D row: transpose conjugation flips L_H, fixes L_diss.
            var dOnH = Dense.Transpose(Dense.CommutatorFlow(h, sigmaT));
            var lh = Dense.CommutatorFlow(h, sigma);
            dH = Math.Max(dH, Dense.MaxAbsSum(dOnH, lh));                        // dOnH + lh = 0

            var dOnDiss = Dense.Transpose(Dense.Dephasing(zs, gammas, sigmaT));
            var ldiss = Dense.Dephasing(zs, gammas, sigma);
            dDiss = Math.Max(dDiss, Dense.MaxAbsDiff(dOnDiss, ldiss));           // dOnDiss - ldiss = 0

            // R row: the one-sided reflection fixes L_H, reflects L_diss and carries the constant.
            var rOnH = Dense.Mul(Dense.CommutatorFlow(h, sigmaF), f);
            rH = Math.Max(rH, Dense.MaxAbsDiff(rOnH, lh));                       // rOnH - lh = 0

            var rOnDiss = Dense.Mul(Dense.Dephasing(zs, gammas, sigmaF), f);
            double worstCell = 0;
            int d = 1 << N;
            for (int a = 0; a < d; a++)
                for (int b = 0; b < d; b++)
                    worstCell = Math.Max(worstCell,
                        (rOnDiss[a, b] + ldiss[a, b] + 2.0 * sigmaTot * sigma[a, b]).Magnitude);
            rDiss = Math.Max(rDiss, worstCell);                                  // rOnDiss + ldiss + 2 sigma = 0
        }
        return (dH, dDiss, rH, rDiss);
    }
}

// small dense pieces, shared by the group and the triangle (nothing external; exact integer /
// half-integer entries wherever the identity being pinned is exact).
internal static class Dense
{
    static readonly Complex[][,] Sigma =
    {
        new Complex[,] { { 1, 0 }, { 0, 1 } },                                                    // I
        new Complex[,] { { 0, 1 }, { 1, 0 } },                                                    // X
        new Complex[,] { { 0, -Complex.ImaginaryOne }, { Complex.ImaginaryOne, 0 } },             // Y
        new Complex[,] { { 1, 0 }, { 0, -1 } },                                                   // Z
    };

    static int Idx(char c) => c switch { 'I' => 0, 'X' => 1, 'Y' => 2, 'Z' => 3, _ => throw new ArgumentException($"{c}") };

    public static IEnumerable<char[]> AllStrings(int n)
    {
        char[] alphabet = { 'I', 'X', 'Y', 'Z' };
        int total = 1 << (2 * n);
        for (int idx = 0; idx < total; idx++)
        {
            var l = new char[n];
            int x = idx;
            for (int s = 0; s < n; s++) { l[s] = alphabet[x & 3]; x >>= 2; }
            yield return l;
        }
    }

    public static Complex[,] PauliString(char[] letters)
    {
        var m = Sigma[Idx(letters[0])];
        for (int s = 1; s < letters.Length; s++) m = Kron(m, Sigma[Idx(letters[s])]);
        return m;
    }

    public static Complex[,] SiteZ(int n, int l)
    {
        var letters = Enumerable.Repeat('I', n).ToArray();
        letters[l] = 'Z';
        return PauliString(letters);
    }

    public static Complex[,] XxzChain(int n, double j, double delta)
    {
        int d = 1 << n;
        var h = new Complex[d, d];
        foreach (var (a, b) in Topology.Chain(n))
            foreach (var (letter, weight) in new[] { ('X', j), ('Y', j), ('Z', j * delta) })
            {
                var letters = Enumerable.Repeat('I', n).ToArray();
                letters[a] = letter; letters[b] = letter;
                var term = PauliString(letters);
                for (int p = 0; p < d; p++)
                    for (int q = 0; q < d; q++)
                        h[p, q] += weight * term[p, q];
            }
        return h;
    }

    public static Complex[,] Kron(Complex[,] a, Complex[,] b)
    {
        int ra = a.GetLength(0), ca = a.GetLength(1), rb = b.GetLength(0), cb = b.GetLength(1);
        var m = new Complex[ra * rb, ca * cb];
        for (int i = 0; i < ra; i++)
            for (int j = 0; j < ca; j++)
                for (int k = 0; k < rb; k++)
                    for (int l = 0; l < cb; l++)
                        m[i * rb + k, j * cb + l] = a[i, j] * b[k, l];
        return m;
    }

    public static Complex[,] Mul(Complex[,] a, Complex[,] b)
    {
        int d = a.GetLength(0);
        var m = new Complex[d, d];
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
            {
                Complex s = Complex.Zero;
                for (int k = 0; k < d; k++) s += a[i, k] * b[k, j];
                m[i, j] = s;
            }
        return m;
    }

    public static Complex[,] Transpose(Complex[,] a)
    {
        int d = a.GetLength(0);
        var m = new Complex[d, d];
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                m[i, j] = a[j, i];
        return m;
    }

    public static Complex[,] Conjugate(Complex[,] a)
    {
        int d = a.GetLength(0);
        var m = new Complex[d, d];
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                m[i, j] = Complex.Conjugate(a[i, j]);
        return m;
    }

    public static Complex[,] Dagger(Complex[,] a) => Conjugate(Transpose(a));

    public static Complex[,] Scale(Complex s, Complex[,] a)
    {
        int d = a.GetLength(0);
        var m = new Complex[d, d];
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                m[i, j] = s * a[i, j];
        return m;
    }

    public static Complex Trace(Complex[,] a)
    {
        Complex t = Complex.Zero;
        for (int i = 0; i < a.GetLength(0); i++) t += a[i, i];
        return t;
    }

    // L_H = -i[H, .]
    public static Complex[,] CommutatorFlow(Complex[,] h, Complex[,] rho)
    {
        int d = h.GetLength(0);
        var hr = Mul(h, rho);
        var rh = Mul(rho, h);
        var m = new Complex[d, d];
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                m[i, j] = -Complex.ImaginaryOne * (hr[i, j] - rh[i, j]);
        return m;
    }

    // the Z-dephasing dissipator: sum_l gamma_l (Z_l rho Z_l - rho)
    public static Complex[,] Dephasing(Complex[][,] zs, double[] gammas, Complex[,] rho)
    {
        int d = rho.GetLength(0);
        var m = new Complex[d, d];
        for (int l = 0; l < gammas.Length; l++)
        {
            var zrz = Mul(zs[l], Mul(rho, zs[l]));
            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                    m[i, j] += gammas[l] * (zrz[i, j] - rho[i, j]);
        }
        return m;
    }

    public static double MaxAbsDiff(Complex[,] a, Complex[,] b)
    {
        double worst = 0;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++)
                worst = Math.Max(worst, (a[i, j] - b[i, j]).Magnitude);
        return worst;
    }

    public static double MaxAbsSum(Complex[,] a, Complex[,] b)
    {
        double worst = 0;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = 0; j < a.GetLength(1); j++)
                worst = Math.Max(worst, (a[i, j] + b[i, j]).Magnitude);
        return worst;
    }

    // a deterministic dense matrix (no RNG in this world): integer real parts in [-2, 2], integer
    // imaginary parts in [-1, 1], generically non-symmetric and non-Hermitian.
    public static Complex[,] TestMatrix(int d, int seedA, int seedB)
    {
        var m = new Complex[d, d];
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                m[i, j] = new Complex((i * seedA + j * seedB) % 5 - 2, (i * seedB + 2 * j + 1) % 3 - 1);
        return m;
    }
}
