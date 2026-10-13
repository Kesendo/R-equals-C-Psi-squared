using System.Numerics;

namespace MirrorWorld;

// The antilinear triangle (adopted 2026-07-04 from F119, docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md):
// the transpose theta(A) = A^T, the entrywise conjugation conj(A) = A-bar, and the adjoint
// dagger(A) = A-dagger satisfy dagger = theta o conj and form, with the identity, a Klein
// four-group graded by two +-1 characters: ell (linearity) and m (multiplicativity) --
// theta = (linear, reversing), conj = (antilinear, preserving), dagger = (antilinear, reversing).
// On a Pauli string: theta(sigma) = conj(sigma) = (-1)^nY sigma, dagger(sigma) = sigma.
//
// The engine is the TRANSPORT LAW
//
//   mu o L_H o mu = ell(mu)*m(mu) * L_{mu(H)},   L_H = -i[H, .],  any H,
//
// one sign from the -i (antilinearity), one from the commutator's order (reversal); the dephasing
// dissipator is fixed by all three. The fixed-point collapse ties the vertices: each vertex's
// fixed-point set is where the other two agree (H = H-dagger iff H^T = H-bar).
//
// In the Pauli basis the triangle docks onto the mirror group: theta = D, dagger = the antilinear
// unit K, conj = D o K, and the closure <R, D, K> is the ANTILINEAR DOUBLE, D4 x Z2 -- order 16,
// exactly eight antiunitary members. (conj commutes with R; dagger does not -- the side flips.)
//
// One level up, theta transposes the generator's matrix: for a unital L the transpose is the
// generator of the transposed system, and read through the trace pairing that is reciprocity
// (preparation and readout trade places; see the section below FixedPointCollapse).
public sealed class AntilinearTriangle : GameObject
{
    public int N { get; }

    public AntilinearTriangle(World world, int n) : base(world)
    {
        if (n < 1 || n > 15)
            throw new ArgumentOutOfRangeException(nameof(n), n, "dense Pauli enumeration needs 1 <= N <= 15");
        N = n;
    }

    // left: what the triangle itself produces.
    public override IReadOnlyList<string> Own => new[] { "vertices", "transport", "double", "reciprocity" };

    // the three vertices, as members of the (doubled) mirror group:
    public static MirrorGroup.Member Theta => MirrorGroup.D;                       // linear, reversing
    public static readonly MirrorGroup.Member Conj =
        MirrorGroup.Compose(MirrorGroup.D, MirrorGroup.K);                         // antilinear, preserving
    public static MirrorGroup.Member Dagger => MirrorGroup.K;                      // antilinear, reversing

    // ell: +1 linear, -1 antilinear.
    public static int Ell(MirrorGroup.Member vertex) => vertex.Antilinear ? -1 : +1;

    // m: +1 multiplicative (preserving), -1 antimultiplicative (reversing). The four V4 vertices only.
    public static int Em(MirrorGroup.Member vertex)
    {
        if (vertex.Equals(MirrorGroup.Identity) || vertex.Equals(Conj)) return +1;
        if (vertex.Equals(Theta) || vertex.Equals(Dagger)) return -1;
        throw new ArgumentException("m is graded on the triangle's four vertices only.");
    }

    // ---- the transport law on a deterministic NON-Hermitian H: for each vertex mu,
    // mu(L_H(mu(rho))) = ell*m * L_{mu(H)}(rho), machine-exact (integer entries). ----
    public double TransportWorstResidual()
    {
        int d = 1 << N;
        var h = Dense.TestMatrix(d, 3, 5);
        var rho = Dense.TestMatrix(d, 7, 2);
        var vertices = new (Func<Complex[,], Complex[,]> Mu, Complex[,] MuH, int Sign)[]
        {
            (Dense.Transpose, Dense.Transpose(h), Ell(Theta) * Em(Theta)),         // theta -> -L_{H^T}
            (Dense.Conjugate, Dense.Conjugate(h), Ell(Conj) * Em(Conj)),           // conj  -> -L_{H-bar}
            (Dense.Dagger, Dense.Dagger(h), Ell(Dagger) * Em(Dagger)),             // dagger -> +L_{H-dagger}
        };
        double worst = 0;
        foreach (var (mu, muH, sign) in vertices)
        {
            var lhs = mu(Dense.CommutatorFlow(h, mu(rho)));
            var rhs = Dense.Scale(sign, Dense.CommutatorFlow(muH, rho));
            worst = Math.Max(worst, Dense.MaxAbsDiff(lhs, rhs));
        }
        return worst;
    }

    // ---- the dephasing dissipator is fixed by all three vertices, verified on the full Pauli basis. ----
    public double DissipatorFixedWorstResidual(double[] gammas)
    {
        var zs = Enumerable.Range(0, N).Select(l => Dense.SiteZ(N, l)).ToArray();
        var vertices = new Func<Complex[,], Complex[,]>[] { Dense.Transpose, Dense.Conjugate, Dense.Dagger };
        double worst = 0;
        foreach (var letters in Dense.AllStrings(N))
        {
            var sigma = Dense.PauliString(letters);
            var ldiss = Dense.Dephasing(zs, gammas, sigma);
            foreach (var mu in vertices)
                worst = Math.Max(worst, Dense.MaxAbsDiff(mu(Dense.Dephasing(zs, gammas, mu(sigma))), ldiss));
        }
        return worst;
    }

    // ---- the fixed-point collapse: H^T = H-bar iff H = H-dagger. A Hermitian H agrees exactly;
    // the deterministic non-Hermitian H splits at O(1). ----
    public (double Hermitian, double NonHermitian) FixedPointCollapse()
    {
        int d = 1 << N;
        var m = Dense.TestMatrix(d, 3, 5);
        var hermitian = new Complex[d, d];
        var mDagger = Dense.Dagger(m);
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                hermitian[i, j] = m[i, j] + mDagger[i, j];
        return (Dense.MaxAbsDiff(Dense.Transpose(hermitian), Dense.Conjugate(hermitian)),
                Dense.MaxAbsDiff(Dense.Transpose(m), Dense.Conjugate(m)));
    }

    // ---- one level up: the generator's transpose (docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md §7). The matrix
    // of a superoperator S in the matrix-unit basis carries S(E_ij)[k,l] at row (k,l), column (i,j). Its
    // transpose is again a Lindbladian, the one of the TRANSPOSED system (H -> H^T, every jump J -> J^T),
    // exactly when L is unital (sum_k r_k [J_k, J_k^dagger] = 0); normal jumps are the per-jump case and
    // every Pauli string is one. Read through the trace pairing this trades preparation and readout,
    // Tr(B e^{tL}(A)) = Tr(A e^{tL'}(B)) for symmetric A and B, and a real H makes L' = L: preparing A and
    // reading B gives what preparing B and reading A gives. With a flux (H Hermitian, not real) the way
    // back is the system with the flux reversed, H^T = H-bar. ----

    // L(rho) = -i[H, rho] + sum_k r_k (J_k rho J_k^dagger - 1/2 {J_k^dagger J_k, rho}), any jumps.
    public static Complex[,] Lindblad(Complex[,] h, Complex[][,] jumps, double[] rates, Complex[,] rho)
    {
        int d = rho.GetLength(0);
        var m = Dense.CommutatorFlow(h, rho);
        for (int k = 0; k < jumps.Length; k++)
        {
            var jd = Dense.Dagger(jumps[k]);
            var jdj = Dense.Mul(jd, jumps[k]);
            var sandwich = Dense.Mul(jumps[k], Dense.Mul(rho, jd));
            var left = Dense.Mul(jdj, rho);
            var right = Dense.Mul(rho, jdj);
            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                    m[i, j] += rates[k] * (sandwich[i, j] - 0.5 * (left[i, j] + right[i, j]));
        }
        return m;
    }

    // max over all matrix units of |forward(E_ij)[k,l] - back(E_kl)[i,j]|: zero iff the matrix of back is
    // the transpose of the matrix of forward.
    public static double TransposeResidual(Func<Complex[,], Complex[,]> forward, Func<Complex[,], Complex[,]> back, int d)
    {
        var imageForward = new Complex[d * d][,];
        var imageBack = new Complex[d * d][,];
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
            {
                var unit = new Complex[d, d];
                unit[i, j] = Complex.One;
                imageForward[i * d + j] = forward(unit);
                imageBack[i * d + j] = back(unit);
            }
        double worst = 0;
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                for (int k = 0; k < d; k++)
                    for (int l = 0; l < d; l++)
                        worst = Math.Max(worst, (imageForward[i * d + j][k, l] - imageBack[k * d + l][i, j]).Magnitude);
        return worst;
    }

    // order by order in t: max over p = 0..maxOrder of |Tr(B forward^p(A)) - Tr(A back^p(B))|, the
    // Taylor coefficients of "prepared as A, read as B" against "prepared as B, read as A".
    public static double ReciprocityResidual(Func<Complex[,], Complex[,]> forward, Func<Complex[,], Complex[,]> back,
        Complex[,] a, Complex[,] b, int maxOrder)
    {
        var x = a;
        var y = b;
        double worst = 0;
        for (int p = 0; p <= maxOrder; p++)
        {
            worst = Math.Max(worst, (Dense.Trace(Dense.Mul(b, x)) - Dense.Trace(Dense.Mul(a, y))).Magnitude);
            x = forward(x);
            y = back(y);
        }
        return worst;
    }

    // what the swap is and is not. On a generic real H (no mirror symmetry) with unequal Z rates, take
    // "prepare an excitation on site 0 over a quiet rest, read site N-1 with the rest ignored". Its POSITIONAL
    // mirror (prepare site N-1, read site 0, same background habits) is a different experiment and differs
    // already in this linear readout; its OPERATOR swap (prepare what was read, read what was prepared: the
    // rest now starts maximally mixed and is read quiet) agrees exactly, order by order in t.
    public (double PositionalMirror, double OperatorSwap) BackgroundSwapResiduals(int maxOrder)
    {
        if (N < 2)
            throw new InvalidOperationException("two ends need N >= 2");
        int d = 1 << N;
        var h = RealSymmetricTestMatrix(d, 2, 3);
        var zs = SiteZs(N);
        var rates = Rates(N);
        Complex[,] L(Complex[,] r) => Lindblad(h, zs, rates, r);
        int first = 1 << (N - 1), last = 1;                       // the bit of site 0 and of site N-1
        var prepareFirst = new Complex[d, d];
        prepareFirst[first, first] = 1;                           // |1 0 ... 0>, the rest quiet
        var prepareLast = new Complex[d, d];
        prepareLast[last, last] = 1;                              // |0 ... 0 1>
        var readFirst = new Complex[d, d];
        var readLast = new Complex[d, d];
        for (int s = 0; s < d; s++)
        {
            if ((s & first) != 0) readFirst[s, s] = 1;            // site 0 excited, the rest ignored
            if ((s & last) != 0) readLast[s, s] = 1;              // site N-1 excited, the rest ignored
        }
        double positional = 0, swap = 0;
        Complex[,] x = prepareFirst, y = prepareLast, z = readLast;
        for (int p = 0; p <= maxOrder; p++)
        {
            var forward = Dense.Trace(Dense.Mul(readLast, x));
            positional = Math.Max(positional, (forward - Dense.Trace(Dense.Mul(readFirst, y))).Magnitude);
            swap = Math.Max(swap, (forward - Dense.Trace(Dense.Mul(prepareFirst, z))).Magnitude);
            x = L(x);
            y = L(y);
            z = L(z);
        }
        return (positional, swap);
    }

    // when positions DO trade: on a number-conserving real H (uneven XXZ bonds of weights 1, 2, 3, 1, ...,
    // anisotropy delta, along a chain or out of a star's hub at site N/2) with unequal Z rates, Z dephasing
    // keeps every (ket, bra) excitation block apart, and inside the block one excitation reaches, "read site
    // b, the rest ignored" IS the projector on the excitation at b. So sending one excitation from end site
    // 0 to end site N-1 over the rest in |0>, or a |+> (whose signal lives in the (0,1) coherence block), is
    // its own operator swap: the positional reversal agrees exactly. Read at an interior site instead, the
    // positional mirror is no longer the swap and differs. A Bell pair mixes the 0 and 2 excitation blocks;
    // the far pair's |11> population and |00><11| coherence are transfers the corollary fixes, but its |00>
    // population collects, in the two-excitation block, every state with the far pair empty, starting with
    // the preparation's own state (a return probability), which the corollary does not tie to the mirrored
    // experiment. So the Bell readout agrees only by some further symmetry: the spin flip X^N (which carries
    // |1100> to |0011> at N = 4, and which an uneven longitudinal field sum_l (l+1) Z_l breaks), or, at
    // delta = 0 on the chain, the free-fermion closure of the one-body correlations under Jordan-Wigner (Z
    // dephasing is quadratic in the fermions, so the far pair's occupation is a one-excitation signal the
    // corollary fixes). On longer uneven chains at delta = 1, 2 and on the XX star (whose Jordan-Wigner
    // strings break the closure) it differs.
    public (double SingleExcitation, double Interior, double Plus, double BellPair) NumberConservingPositionalResiduals(
        int delta, int maxOrder, bool star = false, bool field = false)
    {
        if (N < 4)
            throw new InvalidOperationException("two disjoint end pairs need N >= 4");
        int d = 1 << N;
        var h = NumberConservingHamiltonian(delta, star, field);
        var zs = SiteZs(N);
        var rates = Rates(N);
        Complex[,] L(Complex[,] r) => Lindblad(h, zs, rates, r);

        int Bit(int site) => 1 << (N - 1 - site);
        Complex[,] Population(int state)
        {
            var m = new Complex[d, d];
            m[state, state] = 1;
            return m;
        }
        Complex[,] Excited(int site)
        {
            var m = new Complex[d, d];
            for (int s = 0; s < d; s++)
                if ((s & Bit(site)) != 0) m[s, s] = 1;
            return m;
        }
        Complex[,] PlusOverQuiet(int site)
        {
            var m = new Complex[d, d];
            foreach (int r in new[] { 0, Bit(site) })
                foreach (int c in new[] { 0, Bit(site) })
                    m[r, c] = 0.5;
            return m;
        }
        Complex[,] SiteX(int site) => Dense.PauliString(Enumerable.Range(0, N).Select(l => l == site ? 'X' : 'I').ToArray());
        Complex[,] BellOverQuiet(int first)
        {
            var m = new Complex[d, d];
            int both = Bit(first) | Bit(first + 1);
            foreach (int r in new[] { 0, both })
                foreach (int c in new[] { 0, both })
                    m[r, c] = 0.5;
            return m;
        }
        Complex[,] BellRead(int first)
        {
            var m = new Complex[d, d];
            int both = Bit(first) | Bit(first + 1);
            for (int s = 0; s < d; s++)
            {
                int pair = s & both;
                if (pair != 0 && pair != both) continue;
                int partner = s ^ both;                           // |00> <-> |11> on the pair, the rest kept
                m[s, s] += 0.5;
                m[s, partner] += 0.5;
            }
            return m;
        }

        double Residual(Complex[,] prepareHere, Complex[,] readThere, Complex[,] prepareThere, Complex[,] readHere)
        {
            double worst = 0;
            Complex[,] x = prepareHere, y = prepareThere;
            for (int p = 0; p <= maxOrder; p++)
            {
                worst = Math.Max(worst, (Dense.Trace(Dense.Mul(readThere, x)) - Dense.Trace(Dense.Mul(readHere, y))).Magnitude);
                x = L(x);
                y = L(y);
            }
            return worst;
        }

        int last = N - 1;
        return (Residual(Population(Bit(0)), Excited(last), Population(Bit(last)), Excited(0)),
                Residual(Population(Bit(0)), Excited(last - 1), Population(Bit(last)), Excited(1)),
                Residual(PlusOverQuiet(0), SiteX(last), PlusOverQuiet(last), SiteX(0)),
                Residual(BellOverQuiet(0), BellRead(N - 2), BellOverQuiet(N - 2), BellRead(0)));
    }

    // the silence reading: on a number-conserving real H (any graph, unequal Z rates) send a Bell pair from
    // the pair at sites 0, 1 over the rest in |0>, and read the pair at sites N-2, N-1 WITH THE SILENCE (the
    // rest projected on |0...0>, the pair kept as it is). That pair lives in span{|00>, |11>}, built from the
    // stationary vacuum population, the |11> transfer and the vacuum-to-|11> coherence, all three fixed by
    // the corollary, so the same experiment run from the other pair gives the same pair state entry for
    // entry, and with it every function of it (CPsi, R = CPsi^2, a mutual information). Compared order by
    // order in t, every entry of the unnormalized silence-read pair; the habit reading (the rest traced out)
    // as the control.
    public (double Silence, double Habit) SilenceReadingResiduals(int delta, int maxOrder, bool star = false)
    {
        if (N < 4)
            throw new InvalidOperationException("two disjoint pairs need N >= 4");
        int d = 1 << N;
        var h = NumberConservingHamiltonian(delta, star, field: false);
        var zs = SiteZs(N);
        var rates = Rates(N);
        Complex[,] L(Complex[,] r) => Lindblad(h, zs, rates, r);
        int Bit(int site) => 1 << (N - 1 - site);
        Complex[,] BellOverQuiet(int first)
        {
            var m = new Complex[d, d];
            int both = Bit(first) | Bit(first + 1);
            foreach (int r in new[] { 0, both })
                foreach (int c in new[] { 0, both })
                    m[r, c] = 0.5;
            return m;
        }
        // the pair's basis state ab (a on the pair's first site) as a full basis index, the rest in |0>
        int Index(int first, int ab) => ((ab >> 1) & 1) * Bit(first) + (ab & 1) * Bit(first + 1);
        Complex[,] Habit(Complex[,] rho, int first)
        {
            var m = new Complex[4, 4];
            int pairMask = Bit(first) | Bit(first + 1);
            for (int rest = 0; rest < d; rest++)
            {
                if ((rest & pairMask) != 0) continue;             // each configuration of the rest once
                for (int ab = 0; ab < 4; ab++)
                    for (int cd = 0; cd < 4; cd++)
                        m[ab, cd] += rho[rest | Index(first, ab), rest | Index(first, cd)];
            }
            return m;
        }
        int near = 0, far = N - 2;
        Complex[,] x = BellOverQuiet(near), y = BellOverQuiet(far);
        double silence = 0, habit = 0;
        for (int p = 0; p <= maxOrder; p++)
        {
            for (int ab = 0; ab < 4; ab++)
                for (int cd = 0; cd < 4; cd++)
                    silence = Math.Max(silence,
                        (x[Index(far, ab), Index(far, cd)] - y[Index(near, ab), Index(near, cd)]).Magnitude);
            habit = Math.Max(habit, Dense.MaxAbsDiff(Habit(x, far), Habit(y, near)));
            x = L(x);
            y = L(y);
        }
        return (silence, habit);
    }

    // uneven XXZ bonds of weights 1, 2, 3, 1, ... along a chain or out of a star's hub at site N/2, with an
    // optional uneven longitudinal field sum_l (l+1) Z_l: real and number-conserving, integer entries.
    Complex[,] NumberConservingHamiltonian(int delta, bool star, bool field)
    {
        int d = 1 << N;
        var h = new Complex[d, d];
        int hub = N / 2;
        var bonds = star
            ? Enumerable.Range(0, N).Where(l => l != hub).Select(l => (Math.Min(l, hub), Math.Max(l, hub))).ToArray()
            : Enumerable.Range(0, N - 1).Select(a => (a, a + 1)).ToArray();
        for (int k = 0; k < bonds.Length; k++)
        {
            var (a, b) = bonds[k];
            int weight = k % 3 + 1;
            foreach (var (letter, scale) in new[] { ('X', 1), ('Y', 1), ('Z', delta) })
            {
                var letters = Enumerable.Repeat('I', N).ToArray();
                letters[a] = letter;
                letters[b] = letter;
                var term = Dense.PauliString(letters);
                for (int p = 0; p < d; p++)
                    for (int q = 0; q < d; q++)
                        h[p, q] += weight * scale * term[p, q];
            }
        }
        if (field)
        {
            var zs = SiteZs(N);
            for (int l = 0; l < N; l++)
                for (int p = 0; p < d; p++)
                    h[p, p] += (l + 1) * zs[l][p, p];
        }
        return h;
    }

    // the transpose law on a deterministic NON-Hermitian H at integer rates, with Pauli-string jumps (Z on
    // every site, the all-X string, a string carrying Y) and, from N = 2 on, the phased cycle: a normal jump
    // that is neither symmetric nor unitary, the one input here on which J^T, J, J-bar and J-dagger differ.
    // Integer entries, so exact.
    public double GeneratorTransposeWorstResidual() =>
        TransposedSystemResidual(Dense.TestMatrix(1 << N, 3, 5), LawJumps(), Dense.Transpose);

    // the jumps carried back unchanged, conjugated, or daggered instead of transposed: each breaks the law
    // on the phased cycle (N >= 2).
    public (double Kept, double Conjugated, double Daggered) WrongJumpMapResiduals()
    {
        if (N < 2)
            throw new InvalidOperationException("the phased cycle needs N >= 2");
        var h = Dense.TestMatrix(1 << N, 3, 5);
        return (TransposedSystemResidual(h, LawJumps(), j => j),
                TransposedSystemResidual(h, LawJumps(), Dense.Conjugate),
                TransposedSystemResidual(h, LawJumps(), Dense.Dagger));
    }

    // the same with an amplitude-damping jump (sigma-minus on site 0) added: L is no longer unital and the
    // law breaks.
    public double NonUnitalJumpTransposeResidual() =>
        TransposedSystemResidual(Dense.TestMatrix(1 << N, 3, 5), SiteZs(N).Append(SiteLowering()).ToArray(), Dense.Transpose);

    // unitality, not normality, is the condition: sigma-minus and sigma-plus on site 0 (neither normal) at
    // EQUAL rates keep L unital and the law exact; at unequal rates (a finite-temperature bath) it breaks.
    public (double Equal, double Unequal) UnitalPairTransposeResidual()
    {
        var h = Dense.TestMatrix(1 << N, 3, 5);
        var lowering = SiteLowering();
        var jumps = new[] { lowering, Dense.Dagger(lowering) };
        return (TransposedSystemResidual(h, jumps, Dense.Transpose, new[] { 1.0, 1.0 }),
                TransposedSystemResidual(h, jumps, Dense.Transpose, new[] { 1.0, 2.0 }));
    }

    // at a finite temperature the transpose is not L but L seen through the thermal weight. Cooling (sigma-minus)
    // at rate g_l and heating (sigma-plus) at q * g_l on every site l, a single ratio q = e^{-beta} for all
    // sites and the rates free per site, Z dephasing at any rates, a real number-conserving H: then
    // L^T = Gamma^-1 L Gamma, Gamma(X) = s X s, with s = q^{n/2} on a basis state of n excitations (pi_beta^{1/2}
    // up to its norm, pi_beta ~ e^{-beta n}). The weight commutes with H and with every Z, and turns sigma-minus
    // into q^{1/2} sigma-minus and sigma-plus into q^{-1/2} sigma-plus, which trades the two jumps' rates. At
    // q = 1 it is the unital law above. q = 1/4 makes s = 2^{-n} dyadic, so every entry is exact. Controls: a
    // different ratio on site 0, an X field on site 0 (H no longer number-conserving), X dephasing on site 0;
    // each breaks the uniform weight.
    public (double Thermal, double UnequalRatio, double XField, double XDephasing) ThermalTransposeResiduals()
    {
        int d = 1 << N;
        var h = NumberConservingHamiltonian(2, star: false, field: true);
        var cool = Enumerable.Range(0, N).Select(l => 4.0 * (l % 3 + 1)).ToArray();
        var xField = Dense.PauliString(Enumerable.Range(0, N).Select(l => l == 0 ? 'X' : 'I').ToArray());
        var hx = new Complex[d, d];
        for (int p = 0; p < d; p++)
            for (int q = 0; q < d; q++)
                hx[p, q] = h[p, q] + xField[p, q];
        double Residual(Complex[,] hh, Complex[][,] jumps, double[] rates)
        {
            Complex[,] L(Complex[,] r) => Lindblad(hh, jumps, rates, r);
            return TransposeResidual(L, r => ThermalWeigh(L(ThermalWeigh(r, inverse: false)), inverse: true), d);
        }
        var (jumps, rates) = ThermalJumps(cool, cool.Select(g => g / 4).ToArray());
        var heatOff = cool.Select(g => g / 4).ToArray();
        heatOff[0] = cool[0] / 16;
        var (jumpsOff, ratesOff) = ThermalJumps(cool, heatOff);
        return (Residual(h, jumps, rates),
                Residual(h, jumpsOff, ratesOff),
                Residual(hx, jumps, rates),
                Residual(h, jumps.Append(xField).ToArray(), rates.Append(1.0).ToArray()));
    }

    // what the thermal transpose says about the two ends. For all A and B,
    // Tr(B L^p(A)) = Tr(Gamma^-1(A^T) L^p(Gamma(B^T))); on basis dyads the weights are numbers: the transfer of
    // |x><y| to |x'><y'| is q^{-(n_x + n_y - n_x' - n_y')/2} times the transfer of |x'><y'| to |x><y|. So the
    // |11> transfer from the pair at sites 0, 1 to the pair at N-2, N-1 and the vacuum-to-|11> coherence carried
    // across (as a complex number) are the same both ways at a single ratio, and one excitation decaying to the
    // vacuum stands to its reverse as 1/q = e^{beta} = 4: detailed balance. The controls that must fail: the
    // decay compared without the weight, and the transfer and the coherence with a different ratio on site 0 (the
    // coherence only where the ratio reaches it within maxOrder, at N = 4 first at t^6). Order by order in t at q = 1/4,
    // exactly.
    public (double Transfer, double Coherence, double Decay, double DecayUnweighted, double TransferUnequalRatio,
        double CoherenceUnequalRatio) ThermalReciprocityResiduals(int maxOrder)
    {
        if (N < 4)
            throw new InvalidOperationException("two disjoint pairs need N >= 4");
        int d = 1 << N;
        var h = NumberConservingHamiltonian(2, star: false, field: true);
        var cool = Enumerable.Range(0, N).Select(l => 4.0 * (l % 3 + 1)).ToArray();
        var (jumps, rates) = ThermalJumps(cool, cool.Select(g => g / 4).ToArray());
        var heatOff = cool.Select(g => g / 4).ToArray();
        heatOff[0] = cool[0] / 16;
        var (jumpsOff, ratesOff) = ThermalJumps(cool, heatOff);
        int Bit(int site) => 1 << (N - 1 - site);
        Complex[,] Dyad(int r, int c)
        {
            var m = new Complex[d, d];
            m[r, c] = Complex.One;
            return m;
        }
        int near = Bit(0) | Bit(1), far = Bit(N - 2) | Bit(N - 1);
        // max over p of |Tr(b L^p(a)) - weight * Tr(a^T L^p(b^T))|
        double Residual(Complex[,] a, Complex[,] b, double weight, bool offRatio = false)
        {
            Complex[,] L(Complex[,] r) => offRatio ? Lindblad(h, jumpsOff, ratesOff, r) : Lindblad(h, jumps, rates, r);
            Complex[,] x = a, y = Dense.Transpose(b);
            var aT = Dense.Transpose(a);
            double worst = 0;
            for (int p = 0; p <= maxOrder; p++)
            {
                worst = Math.Max(worst, (Dense.Trace(Dense.Mul(b, x)) - weight * Dense.Trace(Dense.Mul(aT, y))).Magnitude);
                x = L(x);
                y = L(y);
            }
            return worst;
        }
        return (Residual(Dyad(near, near), Dyad(far, far), 1.0),
                Residual(Dyad(0, near), Dyad(far, 0), 1.0),
                Residual(Dyad(Bit(0), Bit(0)), Dyad(0, 0), 4.0),
                Residual(Dyad(Bit(0), Bit(0)), Dyad(0, 0), 1.0),
                Residual(Dyad(near, near), Dyad(far, far), 1.0, offRatio: true),
                Residual(Dyad(0, near), Dyad(far, 0), 1.0, offRatio: true));
    }

    // Z on every site at rates 1..N, sigma-minus at cool[l] and sigma-plus at heat[l] on every site l.
    (Complex[][,] Jumps, double[] Rates) ThermalJumps(double[] cool, double[] heat)
    {
        var lower = Enumerable.Range(0, N).Select(Lowering).ToArray();
        var jumps = SiteZs(N).Concat(lower).Concat(lower.Select(Dense.Dagger)).ToArray();
        var rates = Rates(N).Concat(cool).Concat(heat).ToArray();
        return (jumps, rates);
    }

    // sigma-minus on one site: clears that site's bit (site 0 the most significant).
    Complex[,] Lowering(int site)
    {
        int d = 1 << N, bit = 1 << (N - 1 - site);
        var m = new Complex[d, d];
        for (int p = 0; p < d; p++)
            if ((p & bit) != 0)
                m[p & ~bit, p] = Complex.One;
        return m;
    }

    // Gamma(X) = s X s with s = 2^{-n} (q = 1/4, s = pi_beta^{1/2} up to its norm), or its inverse.
    Complex[,] ThermalWeigh(Complex[,] x, bool inverse)
    {
        int d = x.GetLength(0);
        var y = new Complex[d, d];
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
            {
                double w = Math.Pow(2.0, BitOperations.PopCount((uint)i) + BitOperations.PopCount((uint)j));
                y[i, j] = inverse ? x[i, j] * w : x[i, j] / w;
            }
        return y;
    }

    // Lindblad is a Lindbladian: it preserves the trace exactly with the non-normal sigma-minus and the
    // phased cycle among its jumps, and with H = 0 and Z jumps it is the world's dephasing dissipator.
    public (double Trace, double Dephasing) LindbladFormResidual()
    {
        int d = 1 << N;
        var rho = Dense.TestMatrix(d, 7, 2);
        var jumps = LawJumps().Append(SiteLowering()).ToArray();
        var image = Lindblad(Dense.TestMatrix(d, 3, 5), jumps, Rates(jumps.Length), rho);
        var zs = SiteZs(N);
        var gammas = Rates(N);
        var dephasing = Lindblad(new Complex[d, d], zs, gammas, rho);
        return (Dense.Trace(image).Magnitude, Dense.MaxAbsDiff(dephasing, Dense.Dephasing(zs, gammas, rho)));
    }

    Complex[][,] LawJumps()
    {
        var jumps = SiteZs(N).Concat(new[]
        {
            Dense.PauliString(Enumerable.Repeat('X', N).ToArray()),
            Dense.PauliString(Enumerable.Range(0, N).Select(l => l == 0 ? 'Y' : 'X').ToArray()),
        });
        return (N >= 2 ? jumps.Append(PhasedCycle(1 << N)) : jumps).ToArray();
    }

    // |0> -> |1> -> i|2>, |2> -> -|0>, |3> -> 2|3>: J-dagger J = J J-dagger = diag(1, 1, 1, 4, 0, ...), so
    // normal; J^T is not +-J, and J is not unitary.
    static Complex[,] PhasedCycle(int d)
    {
        var j = new Complex[d, d];
        j[1, 0] = Complex.One;
        j[2, 1] = Complex.ImaginaryOne;
        j[0, 2] = -Complex.One;
        j[3, 3] = 2;
        return j;
    }

    Complex[,] SiteLowering()
    {
        var lower = new Complex[2, 2];
        lower[0, 1] = Complex.One;
        return N == 1 ? lower : Dense.Kron(lower, Dense.PauliString(Enumerable.Repeat('I', N - 1).ToArray()));
    }

    static double[] Rates(int count) => Enumerable.Range(1, count).Select(r => (double)r).ToArray();

    static double TransposedSystemResidual(Complex[,] h, Complex[][,] jumps, Func<Complex[,], Complex[,]> jumpBack,
        double[]? rates = null)
    {
        rates ??= Rates(jumps.Length);
        var hT = Dense.Transpose(h);
        var jumpsBack = jumps.Select(jumpBack).ToArray();
        return TransposeResidual(r => Lindblad(h, jumps, rates, r), r => Lindblad(hT, jumpsBack, rates, r), h.GetLength(0));
    }

    // a ring of n >= 3 sites, XX+YY+ZZ on every bond, the wrap bond's hopping turned by a quarter flux of
    // sign s: s * (X_{n-1} Y_0 - Y_{n-1} X_0) in place of its XX+YY. Hermitian, not real for s != 0; its
    // transpose is the ring with the flux reversed (s = 0 leaves the wrap bond with ZZ only).
    public static Complex[,] FluxRing(int n, int flux)
    {
        if (n < 3)
            throw new ArgumentOutOfRangeException(nameof(n), n, "a ring needs n >= 3");
        if (flux is < -1 or > 1)
            throw new ArgumentOutOfRangeException(nameof(flux), flux, "the flux sign is -1, 0 or +1");
        int d = 1 << n;
        var h = new Complex[d, d];
        void Add(int a, char la, int b, char lb, double weight)
        {
            var letters = Enumerable.Repeat('I', n).ToArray();
            letters[a] = la;
            letters[b] = lb;
            var term = Dense.PauliString(letters);
            for (int p = 0; p < d; p++)
                for (int q = 0; q < d; q++)
                    h[p, q] += weight * term[p, q];
        }
        for (int a = 0; a < n; a++)
        {
            int b = (a + 1) % n;
            Add(a, 'Z', b, 'Z', 1.0);
            if (b != 0)
            {
                Add(a, 'X', b, 'X', 1.0);
                Add(a, 'Y', b, 'Y', 1.0);
            }
            else
            {
                Add(a, 'X', b, 'Y', flux);
                Add(a, 'Y', b, 'X', -flux);
            }
        }
        return h;
    }

    public static Complex[][,] SiteZs(int n) => Enumerable.Range(0, n).Select(l => Dense.SiteZ(n, l)).ToArray();

    // a generic real symmetric integer matrix, entries in [-3, 3]. Dense.TestMatrix is not used here: the
    // real part of its seed pair (a, b) is the transpose of the real part of (b, a), so the two symmetrize
    // to the same matrix, and a swap test fed (1,4) and (4,1) would compare a matrix with itself.
    public static Complex[,] RealSymmetricTestMatrix(int d, int seedA, int seedB)
    {
        var m = new Complex[d, d];
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                m[i, j] = SymmetricEntry(i, j, seedA, seedB);
        return m;
    }

    // a generic Hermitian integer matrix that is not real: the symmetric real part above plus an
    // antisymmetric imaginary part.
    public static Complex[,] HermitianTestMatrix(int d, int seedA, int seedB)
    {
        int Skew(int i, int j) => (i * i * seedA + j * seedB) % 5;
        var m = new Complex[d, d];
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                m[i, j] = new Complex(SymmetricEntry(i, j, seedA, seedB), Skew(i, j) - Skew(j, i));
        return m;
    }

    static int SymmetricEntry(int i, int j, int seedA, int seedB) =>
        ((i + 1) * (j + 1) * seedA + (i + j) * seedB + i * i + j * j) % 7 - 3;

    public static double TransposeEqualsResidual(Complex[,] x, Complex[,] y) => Dense.MaxAbsDiff(Dense.Transpose(x), y);
}
