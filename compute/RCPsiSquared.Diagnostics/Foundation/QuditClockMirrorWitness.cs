using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for F172 (<c>docs/proofs/PROOF_QUDIT_CLOCK_MIRROR.md</c>): N sites of d levels, the
/// clock Z|j⟩ = ω^j|j⟩ (ω = e^{2πi/d}, the Weyl clock operator, not MirrorWorld's <c>Clock</c>), the shift
/// X|j⟩ = |j + 1⟩, U = (X^{d/2})^⊗N, the level reflection Q = P^⊗N, P|j⟩ = |−j⟩, and Y = (Z^{d/2})^⊗N; the
/// clock pair {Z_l, Z_l†} or the lone jump Z_l at one integer rate γ on every site.
///
/// <para>Every identity is read exactly in ℤ[ζ_m], m = lcm(d, 4), through <see cref="CyclotomicRing"/>: no primes, no
/// norm bound, no floats. A mirror M is a monomial map on the matrix units |j⟩⟨k| (a permutation with roots of unity,
/// after a transpose or not), and M·L·M⁻¹ = −L′ − c is checked as M(L(E)) + L′(M(E)) + c·M(E) = 0 on every matrix unit
/// E, L applied to a unit directly from H and the jumps. That is a different road from the page's gate
/// (<c>simulations/qudit_clock_mirror_gate.py</c>), which reduces whole superoperators modulo two primes and closes them
/// with a norm bound; the Hamiltonians drawn here come from their own seeds.</para>
///
/// <para>Children, the hidden qubit first: X^{d/2} Z X^{−d/2} = −Z on the operators and Z = Z_b ⊗ diag(ω^{j′}), beside
/// the reading that at odd d no ω^{−s} is −1; the pair's one-sided element ρ ↦ Uρ (−L† − 4Nγ) at a Hermitian H that is
/// not real, with a clock field as the control; its transpose mirror ρ ↦ ρᵀU (−L − 4Nγ) at a real symmetric H, with that
/// H that is not real as the control; the far end widened to UΔ (an H commuting with UΔ and not with U); the lone jump's
/// four routes, each later route's Hamiltonian read against the other routes as controls; and the reduction
/// ℒ_U L(H, Z) ℒ_U⁻¹ = −L(−H, Z) − 2Nγ. Every Hamiltonian is checked Hermitian (conjugation the Galois map ζ ↦ ζ⁻¹),
/// since the targets −L† are built as the generator with −H and the adjoint jumps. The counts of the two ends and the
/// moments are not recomputed here; the gate's G14 to G17 hold them.</para></summary>
public sealed class QuditClockMirrorWitness : IInspectable
{
    /// <summary>The largest Hilbert-space dimension d^N the witness builds (every one of the d^{2N} matrix units is
    /// read).</summary>
    public const int MaxDimension = 64;

    public int D { get; }
    public int N { get; }
    public long Gamma { get; }
    public int Dimension { get; }
    public int M { get; }

    /// <summary>One identity read on the matrix units, in order, until the first unit where it fails.</summary>
    public sealed record IdentityReading(string Name, string Map, string Target, bool Holds, int UnitsRead, (int J, int K)? FirstFailure);

    /// <summary>X^{d/2} Z X^{−d/2} = −Z, read on the d × d operators: the clock conjugated by the half-shift.</summary>
    public bool HalfShiftAnticommutes { get; }
    /// <summary>In the labels j = (d/2)·b + j′: ω^j = (−1)^b ω^{j′} at every level, and the half-shift flips b and keeps j′.</summary>
    public bool HiddenQubitFactorizes { get; }
    /// <summary>A reading of the roots of unity that cannot fail: at d = 3, 5, 7 no ω^{−s} equals −1.</summary>
    public IReadOnlyList<(int OddD, bool NoShiftAnticommutes)> OddDimensionReadings { get; }
    /// <summary>Every Hamiltonian built here satisfies H[a, b] = conj(H[b, a]), conjugation the Galois map ζ ↦ ζ⁻¹.</summary>
    public bool HamiltoniansHermitian { get; }
    /// <summary>Each route's own Hamiltonian meets that route's condition and commutes with U, read exactly.</summary>
    public bool RoutesConditionsOnOwnH { get; }

    public IdentityReading PairOneSided { get; }
    public IdentityReading PairOneSidedClockFieldControl { get; }
    public IdentityReading PairTranspose { get; }
    public IdentityReading PairTransposeControl { get; }
    public IdentityReading FarElement { get; }
    public bool FarElementHamiltonianMissesU { get; }

    /// <summary>The draws of (Δ, H₀) allowed for an H that commutes with UΔ and not with U, and the draws used.</summary>
    public const int FarDraws = 20;
    public int FarDrawsUsed { get; }

    public IdentityReading LoneParityRoute { get; }
    public IdentityReading LoneTransposeRoute { get; }
    /// <summary>A route read at another route's Hamiltonian, beside that route's condition on H, computed exactly; the
    /// identity must hold exactly where the condition does. On a short chain an H can meet two routes (at d = 4, N = 1
    /// every imaginary H commuting with U is also reversed by Y), and there the reading holds and is reported so.</summary>
    public sealed record ControlReading(IdentityReading Reading, bool ConditionHolds)
    {
        public bool Agrees => Reading.Holds == ConditionHolds;
    }

    public IReadOnlyList<ControlReading> LoneTransposeRouteControls { get; }
    public IdentityReading LoneHalfPowerRoute { get; }
    public IReadOnlyList<ControlReading> LoneHalfPowerControls { get; }
    public IdentityReading LoneConjugationRoute { get; }
    public IReadOnlyList<ControlReading> LoneConjugationControls { get; }
    public IdentityReading LoneReduction { get; }

    public QuditClockMirrorWitness(int d = 4, int N = 2, long gamma = 1)
    {
        if (d < 4 || d % 2 != 0)
            throw new ArgumentOutOfRangeException(nameof(d), d,
                "even d ≥ 4: at odd d no shift anticommutes with the clock, and at d = 2 the level reflection is the identity " +
                "and the lone jump Hermitian, so its routes do not separate");
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "N ≥ 1");
        if (gamma < 1 || gamma > 1000) throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "1 ≤ γ ≤ 1000, an integer rate");
        long dim = 1;
        for (int l = 0; l < N; l++) dim *= d;
        if (dim > MaxDimension)
            throw new ArgumentOutOfRangeException(nameof(N), N, $"d^N = {dim} exceeds {MaxDimension}");

        D = d;
        this.N = N;
        Gamma = gamma;
        Dimension = (int)dim;
        M = Lcm(d, 4);
        var r = new Ring(M);
        int dimI = Dimension;

        // ----- the site operators as monomials, tensored over the chain -----
        var shiftHalf = Mono.Site(d, j => (j + d / 2) % d, _ => 0);
        var reflect = Mono.Site(d, j => (d - j) % d, _ => 0);
        var halfPower = Mono.Site(d, j => j, j => (M / 2) * j);          // Z^{d/2} = diag((−1)^j)
        var U = Mono.Tensor(shiftHalf, N);
        var Q = Mono.Tensor(reflect, N);
        var Y = Mono.Tensor(halfPower, N);
        var I = Mono.Identity(dimI);

        // the hidden qubit, on the operators: X^{d/2} Z X^{-d/2} = -Z, and Z = Z_b (x) diag(w^{j'})
        var zSite = Ham.Zero(r, d);
        for (int j = 0; j < d; j++) zSite[j * d + j] = r.Root((M / d) * j);
        HalfShiftAnticommutes = Ham.Conjugate(r, zSite, shiftHalf).Zip(Ham.Neg(r, zSite), r.Equal).All(x => x);
        HiddenQubitFactorizes = Enumerable.Range(0, d).All(j =>
        {
            int b = j / (d / 2), jp = j % (d / 2);
            var rhs = r.Root((M / d) * jp);
            if (b == 1) rhs = r.Neg(rhs);
            int flipped = shiftHalf.Perm[j];
            return r.Equal(r.Root((M / d) * j), rhs) && flipped / (d / 2) == 1 - b && flipped % (d / 2) == jp;
        });
        OddDimensionReadings = new[] { 3, 5, 7 }.Select(od =>
        {
            int mo = Lcm(od, 4);
            var ro = new Ring(mo);
            bool none = Enumerable.Range(1, od - 1).All(s => !ro.Equal(ro.Root(-(mo / od) * s), ro.Neg(ro.One)));
            return (od, none);
        }).ToList();

        // ----- Hamiltonians -----
        var rng = new Random(20261107);
        var swap = Ham.SwapChain(r, d, N, 1);
        var realShift = Ham.ShiftField(r, d, N, 2, imaginary: false);
        var imagShift = Ham.ShiftField(r, d, N, 3, imaginary: true);
        var clockField = Ham.ClockField(r, d, N, 1);

        var pairJumps = Enumerable.Range(0, N).SelectMany(l => new[] { (l, 1), (l, -1) }).ToArray();
        var loneJumps = Enumerable.Range(0, N).Select(l => (l, 1)).ToArray();
        var loneDagger = Enumerable.Range(0, N).Select(l => (l, -1)).ToArray();
        long pairC = QuditClockMirrorClaim.PairTwiceCentre(N, gamma), loneC = QuditClockMirrorClaim.LoneTwiceCentre(N, gamma);

        // the pair, one-sided: rho -> U rho at a Hermitian H that is not real
        var hNotReal = Ham.Add(r, swap, imagShift);
        PairOneSided = Read(r, "the pair's one-sided element", "ρ ↦ Uρ", "−L† − 4Nγ",
            new Lv(r, d, N, hNotReal, pairJumps, gamma), Map.Of(U, I, false),
            new Lv(r, d, N, Ham.Neg(r, hNotReal), pairJumps, gamma), pairC);
        var hClock = Ham.Add(r, hNotReal, clockField);
        PairOneSidedClockFieldControl = Read(r, "control: a clock field h(Z + Z†) added", "ρ ↦ Uρ", "−L† − 4Nγ",
            new Lv(r, d, N, hClock, pairJumps, gamma), Map.Of(U, I, false),
            new Lv(r, d, N, Ham.Neg(r, hClock), pairJumps, gamma), pairC);

        // the pair's transpose mirror: rho -> rho^T U at a real symmetric H
        var hReal = Ham.Add(r, swap, realShift);
        PairTranspose = Read(r, "the pair's transpose mirror", "ρ ↦ ρᵀU", "−L − 4Nγ",
            new Lv(r, d, N, hReal, pairJumps, gamma), Map.Of(I, U, true),
            new Lv(r, d, N, hReal, pairJumps, gamma), pairC);
        PairTransposeControl = Read(r, "control: the transpose mirror at the H that is not real", "ρ ↦ ρᵀU", "−L − 4Nγ",
            new Lv(r, d, N, hNotReal, pairJumps, gamma), Map.Of(I, U, true),
            new Lv(r, d, N, hNotReal, pairJumps, gamma), pairC);

        // the far end widened: W = U Delta, H = sum over the cyclic group of W of W^k H0 W^-k
        // a small chain can draw a Delta whose group forces [H, U] = 0 as well; draw again, at most FarDraws times
        Mono W = U;
        long[][] hFar = Ham.Zero(r, dimI);
        for (FarDrawsUsed = 1; FarDrawsUsed <= FarDraws; FarDrawsUsed++)
        {
            var delta = Mono.Diagonal(dimI, _ => (M / d) * rng.Next(d));
            W = Mono.Compose(U, delta);
            hFar = Ham.Symmetrize(r, Ham.RandomReal(r, dimI, rng), W);
            if (!Ham.CommutesWith(r, hFar, U)) break;
        }
        FarElementHamiltonianMissesU = !Ham.CommutesWith(r, hFar, U);
        if (FarDrawsUsed > FarDraws) FarDrawsUsed = FarDraws;
        FarElement = Read(r, "the far end widened to UΔ", "ρ ↦ UΔρ", "−L† − 4Nγ",
            new Lv(r, d, N, hFar, pairJumps, gamma), Map.Of(W, I, false),
            new Lv(r, d, N, Ham.Neg(r, hFar), pairJumps, gamma), pairC);

        // ----- the lone jump -----
        // route 1: [Q, H] = 0, M = L_U Ad_Q, rho -> (UQ) rho Q, target -L^dag - 2N gamma
        var hQ = Ham.Add(r, swap, realShift);
        LoneParityRoute = ReadLoneParity(r, d, N, "route 1: [Q, H] = 0", hQ, U, Q, loneJumps, loneDagger, gamma, loneC);

        // route 2: Q H-bar Q = H, M(rho) = U Q rho^T Q, target -L - 2N gamma
        LoneTransposeRoute = ReadLoneTranspose(r, d, N, "route 2: Q H̄ Q = H", hNotReal, U, Q, loneJumps, gamma, loneC);

        // route 3: Y H Y = -H, M = L_U Ad_Y, target -L - 2N gamma
        var h1 = Ham.Symmetrize(r, Ham.RandomReal(r, dimI, rng), U);
        var hY = Ham.Add(r, h1, Ham.Neg(r, Ham.Conjugate(r, h1, Y)));
        LoneHalfPowerRoute = ReadLoneHalfPower(r, d, N, "route 3: Y H Y⁻¹ = −H, Y = (Z^{d/2})^⊗N", hY, U, Y, loneJumps, gamma, loneC);

        // route 4: H-bar = -H, M(rho) = U rho^T, target -L^dag - 2N gamma
        var a1 = Ham.Symmetrize(r, Ham.RandomAntisymmetric(r, dimI, rng), U);
        var hImag = Ham.Scale(r, a1, r.I);
        LoneConjugationRoute = ReadLoneConjugation(r, d, N, "route 4: H̄ = −H", hImag, U, I, loneJumps, loneDagger, gamma, loneC);

        // each later route's Hamiltonian read against the other routes, beside their conditions on H
        bool C1(long[][] h) => Ham.CommutesWith(r, h, Q);
        bool C2(long[][] h) => Ham.Equal(r, Ham.Conjugate(r, Ham.Conj(r, h), Q), h);
        bool C3(long[][] h) => Ham.Equal(r, Ham.Conjugate(r, h, Y), Ham.Neg(r, h));
        bool C4(long[][] h) => Ham.Equal(r, Ham.Conj(r, h), Ham.Neg(r, h));
        LoneTransposeRouteControls = new[]
        {
            new ControlReading(ReadLoneParity(r, d, N, "route 1 at route 2's H", hNotReal, U, Q, loneJumps, loneDagger, gamma, loneC), C1(hNotReal)),
            new ControlReading(ReadLoneHalfPower(r, d, N, "route 3 at route 2's H", hNotReal, U, Y, loneJumps, gamma, loneC), C3(hNotReal)),
            new ControlReading(ReadLoneConjugation(r, d, N, "route 4 at route 2's H", hNotReal, U, I, loneJumps, loneDagger, gamma, loneC), C4(hNotReal)),
        };
        LoneHalfPowerControls = new[]
        {
            new ControlReading(ReadLoneParity(r, d, N, "route 1 at route 3's H (route 2's condition there, H being real)", hY, U, Q, loneJumps, loneDagger, gamma, loneC), C1(hY)),
            new ControlReading(ReadLoneConjugation(r, d, N, "route 4 at route 3's H", hY, U, I, loneJumps, loneDagger, gamma, loneC), C4(hY)),
        };
        LoneConjugationControls = new[]
        {
            new ControlReading(ReadLoneParity(r, d, N, "route 1 at route 4's H", hImag, U, Q, loneJumps, loneDagger, gamma, loneC), C1(hImag)),
            new ControlReading(ReadLoneTranspose(r, d, N, "route 2 at route 4's H", hImag, U, Q, loneJumps, gamma, loneC), C2(hImag)),
            new ControlReading(ReadLoneHalfPower(r, d, N, "route 3 at route 4's H", hImag, U, Y, loneJumps, gamma, loneC), C3(hImag)),
        };

        // the reduction: L_U L(H, Z) L_U^-1 = -L(-H, Z) - 2N gamma at an H commuting with U only
        var hRed = Ham.Symmetrize(r, Ham.RandomReal(r, dimI, rng), U);
        LoneReduction = Read(r, "the reduction, kicks forward against kicks back", "ρ ↦ Uρ", "−L(−H, Z) − 2Nγ",
            new Lv(r, d, N, hRed, loneJumps, gamma), Map.Of(U, I, false),
            new Lv(r, d, N, Ham.Neg(r, hRed), loneJumps, gamma), loneC);

        RoutesConditionsOnOwnH = C1(hQ) && C2(hNotReal) && C3(hY) && C4(hImag) &&
            new[] { hQ, hNotReal, hY, hImag, hRed }.All(h => Ham.CommutesWith(r, h, U));
        HamiltoniansHermitian = new[] { hNotReal, hClock, hReal, hFar, hQ, hY, hImag, hRed }.All(h => Ham.IsHermitian(r, h, dimI));
    }

    private static IdentityReading ReadLoneHalfPower(Ring r, int d, int N, string name, long[][] h, Mono U, Mono Y,
        (int, int)[] jumps, long gamma, long c) =>
        Read(r, name, "ρ ↦ UYρY", "−L − 2Nγ",
            new Lv(r, d, N, h, jumps, gamma), Map.Of(Mono.Compose(U, Y), Y, false),
            new Lv(r, d, N, h, jumps, gamma), c);

    private static IdentityReading ReadLoneConjugation(Ring r, int d, int N, string name, long[][] h, Mono U, Mono I,
        (int, int)[] jumps, (int, int)[] dagger, long gamma, long c) =>
        Read(r, name, "ρ ↦ Uρᵀ", "−L† − 2Nγ",
            new Lv(r, d, N, h, jumps, gamma), Map.Of(U, I, true),
            new Lv(r, d, N, Ham.Neg(r, h), dagger, gamma), c);

    private static IdentityReading ReadLoneParity(Ring r, int d, int N, string name, long[][] h, Mono U, Mono Q,
        (int, int)[] jumps, (int, int)[] dagger, long gamma, long c) =>
        Read(r, name, "ρ ↦ (UQ)ρQ", "−L† − 2Nγ",
            new Lv(r, d, N, h, jumps, gamma), Map.Of(Mono.Compose(U, Q), Q, false),
            new Lv(r, d, N, Ham.Neg(r, h), dagger, gamma), c);

    private static IdentityReading ReadLoneTranspose(Ring r, int d, int N, string name, long[][] h, Mono U, Mono Q,
        (int, int)[] jumps, long gamma, long c) =>
        Read(r, name, "ρ ↦ UQρᵀQ", "−L − 2Nγ",
            new Lv(r, d, N, h, jumps, gamma), Map.Of(Mono.Compose(U, Q), Q, true),
            new Lv(r, d, N, h, jumps, gamma), c);

    /// <summary>M(L(E)) + L′(M(E)) + c·M(E) = 0 on every matrix unit E.</summary>
    private static IdentityReading Read(Ring r, string name, string map, string target, Lv L, Map mirror, Lv Lp, long c)
    {
        int dim = L.Dim;
        bool ok = true;
        int read = 0;
        (int, int)? fail = null;
        for (int j = 0; j < dim && ok; j++)
            for (int k = 0; k < dim && ok; k++)
            {
                var acc = new Dictionary<int, long[]>();
                foreach (var (idx, v) in L.ApplyUnit(j, k, r.One))
                {
                    var (mj, mk, ph) = mirror.ApplyUnit(idx / dim, idx % dim);
                    Accumulate(r, acc, mj * dim + mk, r.Mul(r.Root(ph), v));
                }
                var (j2, k2, ph2) = mirror.ApplyUnit(j, k);
                var s = r.Root(ph2);
                foreach (var (idx, v) in Lp.ApplyUnit(j2, k2, s)) Accumulate(r, acc, idx, v);
                Accumulate(r, acc, j2 * dim + k2, r.Scale(s, c));
                ok = acc.Values.All(r.IsZero);
                read++;
                if (!ok) fail = (j, k);
            }
        return new IdentityReading(name, map, target, ok, read, fail);
    }

    private static void Accumulate(Ring r, Dictionary<int, long[]> acc, int idx, long[] v)
    {
        acc[idx] = acc.TryGetValue(idx, out var cur) ? r.Add(cur, v) : v;
    }

    private static int Lcm(int a, int b)
    {
        int g = a, h = b;
        while (h != 0) (g, h) = (h, g % h);
        return a / g * b;
    }

    // ===== exact pieces =====

    /// <summary>ℤ[ζ_m] through <see cref="CyclotomicRing"/>.</summary>
    private sealed class Ring
    {
        public int Modulus { get; }
        public long[] One { get; }
        public long[] I { get; }
        public Ring(int m)
        {
            Modulus = m;
            One = Root(0);
            I = Root(m / 4);
        }
        public long[] Zero() => CyclotomicRing.Zero(Modulus);
        public long[] Root(long e)
        {
            var v = Zero();
            CyclotomicRing.AddRootPower(v, Modulus, e, 1);
            return v;
        }
        public long[] Int(long n) { var v = Zero(); v[0] = n; return v; }
        public long[] Add(long[] a, long[] b) => CyclotomicRing.Add(a, b);
        public long[] Neg(long[] a) => CyclotomicRing.Negate(a);
        public long[] Mul(long[] a, long[] b) => CyclotomicRing.Multiply(a, b, Modulus);
        public long[] Scale(long[] a, long n) { var v = new long[a.Length]; checked { for (int i = 0; i < a.Length; i++) v[i] = a[i] * n; } return v; }
        public bool IsZero(long[] a) => CyclotomicRing.IsZero(a);
        /// <summary>Complex conjugation, the Galois map ζ ↦ ζ⁻¹, on the coefficient vector.</summary>
        public long[] Conj(long[] a)
        {
            var v = Zero();
            for (int i = 0; i < a.Length; i++)
                if (a[i] != 0) CyclotomicRing.AddRootPower(v, Modulus, -i, a[i]);
            return v;
        }
        public bool Equal(long[] a, long[] b) => CyclotomicRing.AreEqual(a, b);
    }

    /// <summary>A monomial operator: |j⟩ ↦ ζ^{Phase[j]} |Perm[j]⟩.</summary>
    private sealed record Mono(int[] Perm, int[] Phase)
    {
        public int Dim => Perm.Length;
        public static Mono Identity(int dim) => new(Enumerable.Range(0, dim).ToArray(), new int[dim]);
        public static Mono Site(int d, Func<int, int> perm, Func<int, int> phase) =>
            new(Enumerable.Range(0, d).Select(perm).ToArray(), Enumerable.Range(0, d).Select(phase).ToArray());
        public static Mono Diagonal(int dim, Func<int, int> phase) =>
            new(Enumerable.Range(0, dim).ToArray(), Enumerable.Range(0, dim).Select(phase).ToArray());
        public static Mono Tensor(Mono site, int n)
        {
            var cur = new Mono(new[] { 0 }, new[] { 0 });
            for (int l = 0; l < n; l++)
            {
                int dc = cur.Dim, ds = site.Dim;
                var perm = new int[dc * ds];
                var phase = new int[dc * ds];
                for (int a = 0; a < dc; a++)
                    for (int b = 0; b < ds; b++)
                    {
                        perm[a * ds + b] = cur.Perm[a] * ds + site.Perm[b];
                        phase[a * ds + b] = cur.Phase[a] + site.Phase[b];
                    }
                cur = new Mono(perm, phase);
            }
            return cur;
        }
        /// <summary>(A·B)|j⟩ = A(B|j⟩).</summary>
        public static Mono Compose(Mono a, Mono b) =>
            new(Enumerable.Range(0, b.Dim).Select(j => a.Perm[b.Perm[j]]).ToArray(),
                Enumerable.Range(0, b.Dim).Select(j => b.Phase[j] + a.Phase[b.Perm[j]]).ToArray());
        public int[] InversePerm()
        {
            var inv = new int[Dim];
            for (int j = 0; j < Dim; j++) inv[Perm[j]] = j;
            return inv;
        }
    }

    /// <summary>ρ ↦ A·ρ·B or ρ ↦ A·ρᵀ·B on the matrix units, A and B monomial.</summary>
    private sealed class Map
    {
        private readonly Mono _a, _b;
        private readonly int[] _bInv;
        private readonly bool _transpose;
        private Map(Mono a, Mono b, bool transpose) { _a = a; _b = b; _bInv = b.InversePerm(); _transpose = transpose; }
        public static Map Of(Mono a, Mono b, bool transpose) => new(a, b, transpose);
        /// <summary>|j⟩⟨k| ↦ ζ^phase |j′⟩⟨k′|.</summary>
        public (int J, int K, long Phase) ApplyUnit(int j, int k)
        {
            if (_transpose) (j, k) = (k, j);
            // A|j><k|B: A|j> = z^{a_j}|pi(j)>; <k|B = z^{beta_b}<b| with B.Perm[b] = k
            int b = _bInv[k];
            return (_a.Perm[j], b, (long)_a.Phase[j] + _b.Phase[b]);
        }
    }

    /// <summary>L(ρ) = −i[H, ρ] + γ Σ (Z_l^s ρ Z_l^{s†} − ρ), applied to one matrix unit.</summary>
    private sealed class Lv
    {
        private readonly Ring _r;
        private readonly int _d, _n;
        private readonly (int Site, int Sign)[] _jumps;
        private readonly long _gamma;
        private readonly List<(int Row, long[] V)>[] _cols;   // column j of H: nonzero (a, H[a, j])
        private readonly List<(int Col, long[] V)>[] _rows;   // row k of H: nonzero (b, H[k, b])
        public int Dim { get; }
        public Lv(Ring r, int d, int n, long[][] h, (int, int)[] jumps, long gamma)
        {
            _r = r; _d = d; _n = n; _jumps = jumps; _gamma = gamma;
            Dim = 1;
            for (int l = 0; l < n; l++) Dim *= d;
            if (h.Length != Dim * Dim) throw new ArgumentException("H is not d^N × d^N", nameof(h));
            _cols = new List<(int, long[])>[Dim];
            _rows = new List<(int, long[])>[Dim];
            for (int x = 0; x < Dim; x++) { _cols[x] = new(); _rows[x] = new(); }
            for (int a = 0; a < Dim; a++)
                for (int b = 0; b < Dim; b++)
                {
                    var v = h[a * Dim + b];
                    if (r.IsZero(v)) continue;
                    _cols[b].Add((a, v));
                    _rows[a].Add((b, v));
                }
        }
        private int Level(int idx, int site)
        {
            for (int l = _n - 1; l > site; l--) idx /= _d;
            return idx % _d;
        }
        public IEnumerable<(int Idx, long[] V)> ApplyUnit(int j, int k, long[] s)
        {
            var minusI = _r.Neg(_r.I);
            var minusIs = _r.Mul(minusI, s);
            var plusIs = _r.Mul(_r.I, s);
            foreach (var (a, v) in _cols[j]) yield return (a * Dim + k, _r.Mul(minusIs, v));   // −i H|j><k|
            foreach (var (b, v) in _rows[k]) yield return (j * Dim + b, _r.Mul(plusIs, v));    // +i |j><k| H
            var coef = _r.Zero();
            foreach (var (site, sign) in _jumps)
            {
                int diff = Level(j, site) - Level(k, site);
                coef = _r.Add(coef, _r.Root((long)sign * (_r.Modulus / _d) * diff));
                coef = _r.Add(coef, _r.Neg(_r.One));
            }
            yield return (j * Dim + k, _r.Mul(_r.Scale(coef, _gamma), s));
        }
    }

    /// <summary>Hamiltonians as dense D × D arrays of ring elements, row-major.</summary>
    private static class Ham
    {
        public static long[][] Zero(Ring r, int dim) => Enumerable.Range(0, dim * dim).Select(_ => r.Zero()).ToArray();

        private static int Pow(int d, int n) { int p = 1; for (int i = 0; i < n; i++) p *= d; return p; }

        public static long[][] SwapChain(Ring r, int d, int n, long j)
        {
            int dim = Pow(d, n);
            var h = Zero(r, dim);
            for (int l = 0; l + 1 < n; l++)
            {
                int stride = Pow(d, n - 2 - l);
                for (int x = 0; x < dim; x++)
                {
                    int a = x / stride / d % d, b = x / stride % d;     // levels of sites l and l + 1
                    int y = x + (b - a) * stride * d + (a - b) * stride;
                    h[y * dim + x] = r.Add(h[y * dim + x], r.Int(j));
                }
            }
            return h;
        }

        /// <summary>h Σ_l (X_l + X_l†), or h Σ_l i(X_l − X_l†).</summary>
        public static long[][] ShiftField(Ring r, int d, int n, long hcoef, bool imaginary)
        {
            int dim = Pow(d, n);
            var h = Zero(r, dim);
            for (int l = 0; l < n; l++)
            {
                int stride = Pow(d, n - 1 - l);
                for (int x = 0; x < dim; x++)
                {
                    int a = x / stride % d;
                    int up = x + (((a + 1) % d) - a) * stride;       // X|a> = |a + 1>
                    int down = x + (((a + d - 1) % d) - a) * stride; // X†|a> = |a − 1>
                    var cu = imaginary ? r.Scale(r.I, hcoef) : r.Int(hcoef);
                    var cd = imaginary ? r.Neg(r.Scale(r.I, hcoef)) : r.Int(hcoef);
                    h[up * dim + x] = r.Add(h[up * dim + x], cu);
                    h[down * dim + x] = r.Add(h[down * dim + x], cd);
                }
            }
            return h;
        }

        /// <summary>h Σ_l (Z_l + Z_l†), diagonal.</summary>
        public static long[][] ClockField(Ring r, int d, int n, long hcoef)
        {
            int dim = Pow(d, n);
            var h = Zero(r, dim);
            for (int l = 0; l < n; l++)
            {
                int stride = Pow(d, n - 1 - l);
                for (int x = 0; x < dim; x++)
                {
                    int a = x / stride % d;
                    var v = r.Add(r.Root((r.Modulus / d) * a), r.Root(-(r.Modulus / d) * a));
                    h[x * dim + x] = r.Add(h[x * dim + x], r.Scale(v, hcoef));
                }
            }
            return h;
        }

        public static long[][] RandomReal(Ring r, int dim, Random rng)
        {
            var h = Zero(r, dim);
            for (int a = 0; a < dim; a++)
                for (int b = a; b < dim; b++)
                {
                    long v = rng.Next(-1, 2);
                    h[a * dim + b] = r.Int(v);
                    h[b * dim + a] = r.Int(v);
                }
            return h;
        }

        public static long[][] RandomAntisymmetric(Ring r, int dim, Random rng)
        {
            var h = Zero(r, dim);
            for (int a = 0; a < dim; a++)
                for (int b = a + 1; b < dim; b++)
                {
                    long v = rng.Next(-1, 2);
                    h[a * dim + b] = r.Int(v);
                    h[b * dim + a] = r.Int(-v);
                }
            return h;
        }

        public static long[][] Add(Ring r, long[][] a, long[][] b) => a.Zip(b, r.Add).ToArray();
        public static long[][] Neg(Ring r, long[][] a) => a.Select(r.Neg).ToArray();
        public static long[][] Scale(Ring r, long[][] a, long[] s) => a.Select(v => r.Mul(v, s)).ToArray();

        /// <summary>W H W⁻¹ for a monomial W: (W H W⁻¹)[π(a), π(b)] = ζ^{w_a − w_b} H[a, b].</summary>
        public static long[][] Conjugate(Ring r, long[][] h, Mono w)
        {
            int dim = w.Dim;
            var o = Zero(r, dim);
            for (int a = 0; a < dim; a++)
                for (int b = 0; b < dim; b++)
                {
                    var v = h[a * dim + b];
                    if (r.IsZero(v)) continue;
                    o[w.Perm[a] * dim + w.Perm[b]] = r.Mul(r.Root(w.Phase[a] - w.Phase[b]), v);
                }
            return o;
        }

        /// <summary>Σ_k W^k H W^{−k} over the cyclic group of W.</summary>
        public static long[][] Symmetrize(Ring r, long[][] h, Mono w)
        {
            var acc = h;
            var cur = w;
            for (int k = 1; k <= 4 * r.Modulus; k++)
            {
                if (IsIdentity(r, cur)) return acc;
                acc = Add(r, acc, Conjugate(r, h, cur));
                cur = Mono.Compose(w, cur);
            }
            throw new InvalidOperationException("the monomial's order exceeds its bound");
        }

        private static bool IsIdentity(Ring r, Mono w)
        {
            for (int j = 0; j < w.Dim; j++)
                if (w.Perm[j] != j || !r.Equal(r.Root(w.Phase[j]), r.One)) return false;
            return true;
        }

        public static bool IsHermitian(Ring r, long[][] h, int dim)
        {
            for (int a = 0; a < dim; a++)
                for (int b = 0; b < dim; b++)
                    if (!r.Equal(h[a * dim + b], r.Conj(h[b * dim + a]))) return false;
            return true;
        }

        public static long[][] Conj(Ring r, long[][] h) => h.Select(r.Conj).ToArray();
        public static bool Equal(Ring r, long[][] a, long[][] b) => a.Zip(b, r.Equal).All(x => x);

        public static bool CommutesWith(Ring r, long[][] h, Mono w)
        {
            var c = Conjugate(r, h, w);
            return h.Zip(c, r.Equal).All(x => x);
        }
    }

    // ===== IInspectable =====

    public IEnumerable<IdentityReading> MustHold => new[]
    {
        PairOneSided, PairTranspose, FarElement, LoneParityRoute, LoneTransposeRoute, LoneHalfPowerRoute,
        LoneConjugationRoute, LoneReduction,
    };

    public IEnumerable<ControlReading> LoneControls =>
        LoneTransposeRouteControls.Concat(LoneHalfPowerControls).Concat(LoneConjugationControls);

    /// <summary>The controls whose condition fails: the pair's two, and every lone-jump route read where its condition
    /// does not hold.</summary>
    public IEnumerable<IdentityReading> MustBreak =>
        new[] { PairOneSidedClockFieldControl, PairTransposeControl }
            .Concat(LoneControls.Where(c => !c.ConditionHolds).Select(c => c.Reading));

    /// <summary>Lone-jump routes whose condition holds at another route's Hamiltonian as well; their identity must hold.</summary>
    public IEnumerable<IdentityReading> AlsoReached => LoneControls.Where(c => c.ConditionHolds).Select(c => c.Reading);

    public bool AllHold =>
        HalfShiftAnticommutes && HiddenQubitFactorizes && HamiltoniansHermitian && RoutesConditionsOnOwnH &&
        LoneControls.All(c => c.Agrees) &&
        OddDimensionReadings.All(o => o.NoShiftAnticommutes) && FarElementHamiltonianMissesU &&
        MustHold.All(x => x.Holds) && MustBreak.All(x => !x.Holds);

    public string DisplayName => $"QuditClockMirrorWitness (d={D}, N={N}, γ={Gamma})";

    private static string Verdict(bool ok, string yes, string no) => ok ? yes : no;

    public string Summary =>
        $"exact over ℤ[ζ_{M}] on the {Dimension * Dimension} matrix units: the hidden qubit " +
        $"{Verdict(HalfShiftAnticommutes && HiddenQubitFactorizes, "holds on the operators (X^{d/2} Z X^{−d/2} = −Z, Z = Z_b ⊗ diag(ω^{j′}))", "FAILS")}; " +
        $"every Hamiltonian {Verdict(HamiltoniansHermitian, "Hermitian", "NOT Hermitian")}; " +
        $"{MustHold.Count(x => x.Holds)} of {MustHold.Count()} mirror identities hold on every unit (the pair's one-sided element, its transpose mirror, " +
        $"the far end UΔ, the lone jump's four routes, the reduction), {MustBreak.Count(x => !x.Holds)} of {MustBreak.Count()} controls " +
        $"break where their condition fails{(AlsoReached.Any() ? $", and {AlsoReached.Count(x => x.Holds)} of {AlsoReached.Count()} routes whose condition holds at another route's H as well hold there" : "")}" +
        $"{Verdict(RoutesConditionsOnOwnH, "", "; a route's own H MISSES its condition")}{Verdict(FarElementHamiltonianMissesU, "", $"; no far-end H missing U in {FarDraws} draws, so that row proves nothing")}: " +
        $"{Verdict(AllHold, "F172 read live", "F172 FAILS somewhere below")}.";

    private static IInspectable Node(IdentityReading x, bool mustHold) => new InspectableNode(
        displayName: x.Name,
        summary: $"{x.Map}: M·L·M⁻¹ = {x.Target} " +
                 (x.Holds ? $"holds exactly on all {x.UnitsRead} matrix units"
                          : $"fails at the unit |{x.FirstFailure!.Value.J}⟩⟨{x.FirstFailure!.Value.K}|, the {x.UnitsRead}th read") +
                 (x.Holds == mustHold ? "" : mustHold ? " (EXPECTED TO HOLD)" : " (EXPECTED TO BREAK)"),
        provenance: NodeProvenance.Live);

    private static IInspectable ControlNode(ControlReading c) => Node(
        c.Reading with { Name = $"{(c.ConditionHolds ? "also reached" : "control")}: {c.Reading.Name}, its condition {(c.ConditionHolds ? "holds" : "fails")} there" },
        c.ConditionHolds);

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode(
                displayName: "first: the hidden qubit, on the operators",
                summary: $"X^{{d/2}} Z X^{{−d/2}} {Verdict(HalfShiftAnticommutes, "= −Z", "≠ −Z")} at d = {D}, the clock conjugated by the half-shift, " +
                         $"so U anticommutes with every clock jump; in the labels j = (d/2)·b + j′, ω^j = (−1)^b ω^{{j′}} and the half-shift flips b " +
                         $"and keeps j′: {Verdict(HiddenQubitFactorizes, "Z = Z_b ⊗ diag(ω^{j′}), X^{d/2} = X_b ⊗ 1", "the factorisation FAILS")}. A reading that " +
                         "cannot fail beside it: " +
                         string.Join("; ", OddDimensionReadings.Select(o => $"d = {o.OddD}, {(o.NoShiftAnticommutes ? "no ω^{−s} is −1" : "an ω^{−s} IS −1")}")) +
                         " (the d-th roots of unity are closed under negation exactly at even d).",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "the clock pair {Z, Z†}",
                summary: "one-sided element, transpose mirror and the far end UΔ" +
                         Verdict(FarElementHamiltonianMissesU, $" (its H commutes with UΔ and not with U, found at draw {FarDrawsUsed} of {FarDraws})",
                             $" (no H commuting with UΔ and not with U in {FarDraws} draws)"),
                children: new[] { Node(PairOneSided, true), Node(PairOneSidedClockFieldControl, false), Node(PairTranspose, true),
                                  Node(PairTransposeControl, false), Node(FarElement, true) },
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                displayName: "the lone jump Z: four routes, one of each kind of symmetry",
                summary: "Q (keeps H, Z → Z†), Q after conjugation (keeps H and Z), Y (reverses H, keeps Z), conjugation (reverses H, Z → Z†); " +
                         $"each route's own H {Verdict(RoutesConditionsOnOwnH, "meets its condition and commutes with U", "MISSES its condition")}; " +
                         "the Hamiltonians of routes 2, 3 and 4 are read against the other routes, each beside that route's condition on H, " +
                         "and the identity must hold exactly where the condition does (at a real H route 2's condition is route 1's)",
                children: new[] { Node(LoneParityRoute, true), Node(LoneTransposeRoute, true) }
                    .Concat(LoneTransposeRouteControls.Select(ControlNode))
                    .Concat(new[] { Node(LoneHalfPowerRoute, true) })
                    .Concat(LoneHalfPowerControls.Select(ControlNode))
                    .Concat(new[] { Node(LoneConjugationRoute, true) })
                    .Concat(LoneConjugationControls.Select(ControlNode)).ToList(),
                provenance: NodeProvenance.Live);
            yield return Node(LoneReduction, true);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
