using System.Globalization;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live "one diagonal": the whole functioning of the disagreement-count diagonal,
/// recomputed at inspect time (<c>inspect --root diagonal</c>). Where the typed claims record
/// (<see cref="RCPsiSquared.Core.Symmetry.ThreeDephasingDiagonalsOrbitClaim"/>,
/// <see cref="RCPsiSquared.Core.Symmetry.AbsorptionTheoremClaim"/>) and the throwaway verifiers
/// (<c>simulations/mirror_inventory_d4.py</c>, <c>one_diagonal_mirror_group.py</c>,
/// <c>rung_dynamics.py</c>) compute once, this witness ports both verifiers into the live layer and
/// makes the structure queryable. The C#-witness-first discipline: a Python verifier that outlives its
/// session is a witness waiting to be ported.
///
/// <para>The one diagonal is the disagreement count k = popcount(i⊕j) of a coherence |i⟩⟨j|; under
/// Z-dephasing L_D = γ·(Q − N·I), the integer levels k are the rungs (Re λ = −2γk, AbsorptionTheorem
/// §4.7). Four live nodes expose its full functioning:</para>
/// <list type="number">
///   <item><b>the rungs</b>: k = 0..N, rate −2γk, multiplicity 2^N·C(N,k); the band-edge floor at k=1 is 2γ.</item>
///   <item><b>the three readings</b> (the SYMMETRY, the mirror group D₄ within one diagonal): D fixes Q
///         (D·Q·D = +Q, the rate reading), R reflects it (R·Q·R = −Q, the mirror/palindrome reading,
///         carrying −2Σγ), the {D, 𝓕D} joint-fixed cell (n_Y even ∧ n_Z even) is truly (the judge reading).</item>
///   <item><b>the orbit</b>: {Q_X, Q_Y, Q_Z} (Q_P = Σ_l kron(P_l, P_lᵀ)) is one orbit of the single-qubit
///         Clifford basis-change S₃ ⟨h_zx, h_yz⟩, hence same spectrum; the structure is S₃ ⋉ D₄.</item>
///   <item><b>the dynamics</b> (the COMPLEMENT): L_H = −i[H,·] connects rung k ↔ {k, k±2} — EVEN steps
///         (a hop flips two bits), never k±1 — so the disagreement-count PARITY is conserved (U(1) feature);
///         the band edge (k=1, odd) and the {0,2} survivor (even) are parity-split → a level crossing; and
///         the mirror R maps k → N−k (the palindrome = the count read from the other end).</item>
/// </list>
///
/// <para>Guard: the operator dimension is 4^N (the d²×d² coherence space, d = 2^N); built only for
/// 4^N ≤ <see cref="MaxDim"/> (= 1024), i.e. N ≤ 5. Convention: row-stacking vec, kron(A,B): ρ ↦ AρBᵀ,
/// matching the verifiers + the typed battery.</para></summary>
public sealed class DiagonalWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private const double Tol = 1e-9;

    /// <summary>Largest coherence-space dimension 4^N the live build will materialise.</summary>
    public const int MaxDim = 1024;

    public int N { get; }
    public double Gamma { get; }
    public double Delta { get; }
    public double J { get; }
    public int D { get; }       // 2^N
    public int Dim => D * D;    // 4^N coherence space

    // built once in the ctor (cheap at the default N=3; on-demand for larger N).
    private readonly ComplexMatrix _qX, _qY, _qZ, _R, _Dsuper, _hZX, _hYZ, _LH;
    private readonly int[] _rung;   // popcount(i^j) per coherence index m = i*d + j

    public DiagonalWitness(int n = 3, double gamma = 0.05, double delta = 0.7, double j = 1.0)
    {
        if (n < 2) throw new ArgumentOutOfRangeException(nameof(n), $"N must be ≥ 2; got {n}");
        if (gamma <= 0) throw new ArgumentOutOfRangeException(nameof(gamma), $"gamma must be positive; got {gamma}");
        long dim = 1L << (2 * n);   // 4^N
        if (dim > MaxDim)
            throw new ArgumentOutOfRangeException(nameof(n),
                $"4^N = {dim} exceeds the live-build guard {MaxDim} for N={n}; pick N ≤ 5.");

        N = n; Gamma = gamma; Delta = delta; J = j; D = 1 << n;

        _qX = DephasingDiagonal(N, PauliLetter.X);
        _qY = DephasingDiagonal(N, PauliLetter.Y);
        _qZ = DephasingDiagonal(N, PauliLetter.Z);
        _Dsuper = TransposeSuper(N);
        _R = KetReflection(N);
        _hZX = AdUnitary(Hadamard1(), N);          // Z↔X basis move
        _hYZ = AdUnitary(RxHalfPi1(), N);          // Z↔Y basis move
        _LH = HamiltonianSuper(BuildChainH(N, J, Delta), N);
        _rung = RungLabels(N);
    }

    public string DisplayName => $"the one diagonal (N={N}, γ={Gamma.ToString("0.###", Inv)})";

    public string Summary
    {
        get
        {
            double dRate = MaxAbsDiff(_Dsuper * _qZ * _Dsuper, _qZ);
            double dMirror = MaxAbsDiff(_R * _qZ * _R, _qZ.Multiply(-Complex.One));
            int orbit = OrbitSize(_qZ, new[] { _hZX, _hYZ });
            bool parity = ParityConserved();
            return $"N={N}: {N + 1} rungs (rate −2γk); readings hold (D-fix {dRate:0.0e+00}, R-anti {dMirror:0.0e+00}); "
                 + $"orbit(Q_Z)={orbit} (=3 ⟹ basis-S₃); L_H even-step, k-parity {(parity ? "CONSERVED" : "BROKEN")} (S₃⋉D₄)";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return TheRungs();
            yield return TheThreeReadings();
            yield return TheOrbit();
            yield return TheDynamics();
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    // ---- public live readouts (for tests + callers; the same computations the nodes surface) ----
    /// <summary>L_H connects rungs by these |Δk| (the even-step ladder: {0,2}, never an odd step).</summary>
    public IReadOnlyList<int> RungStepsOfLH() => RungSteps(_LH);
    /// <summary>The disagreement-count parity is conserved (L = L_H + L_D block-diagonal in k mod 2).</summary>
    public bool DisagreementParityConserved() => ParityConserved();
    /// <summary>The orbit of Q_Z under the basis-change S₃ ⟨h_zx, h_yz⟩ (= 3 ⟹ the three diagonals are one orbit).</summary>
    public int OrbitSizeOfQZ() => OrbitSize(_qZ, new[] { _hZX, _hYZ });
    /// <summary>spec(Q_X) = spec(Q_Y) = spec(Q_Z) (conjugate ⟹ co-spectral).</summary>
    public bool ThreeDiagonalsSameSpectrum() => SpectraEqual(_qX, _qY) && SpectraEqual(_qX, _qZ);
    /// <summary>The §3 Hamiltonian column of the palindrome split (PROOF_PI_FACTORS_AS_R_TIMES_D): the
    /// devs of D·L_H·D = −L_H (D flips L_H) and R·L_H·R = +L_H (R fixes L_H). The mirror group's action
    /// on L_H — the bridge between the symmetry (D₄) and the even-step rung dynamics (same L_H, two views).</summary>
    public (double DFlipLH, double RFixLH) MirrorActionOnLH() =>
        (MaxAbsDiff(_Dsuper * _LH * _Dsuper, _LH.Multiply(-Complex.One)),
         MaxAbsDiff(_R * _LH * _R, _LH));

    // ---------------------------------------------------------------- node 1: the rungs
    private IInspectable TheRungs()
    {
        var nodes = new List<IInspectable>();
        for (int k = 0; k <= N; k++)
        {
            long mult = (1L << N) * Binom(N, k);   // 2^N * C(N,k) coherences at rung k
            string note = k == 0 ? " (populations, dark)" : k == 1 ? " (the band-edge floor = 2γ)" : "";
            nodes.Add(new InspectableNode($"rung k={k}",
                summary: $"rate Re λ = −2γk = {(-2 * Gamma * k).ToString("0.####", Inv)}; "
                       + $"multiplicity 2^N·C(N,k) = {mult}{note}"));
        }
        return new InspectableNode("the rungs (the one diagonal)",
            summary: $"k = popcount(i⊕j) ∈ {{0..{N}}}; L_D = γ·(Q − N·I), Re λ = −2γk (AbsorptionTheorem §4.7); "
                   + "the light sets the rungs (topology-free).",
            children: nodes);
    }

    // ---------------------------------------------------------------- node 2: the three readings (D₄)
    private IInspectable TheThreeReadings()
    {
        double dRate = MaxAbsDiff(_Dsuper * _qZ * _Dsuper, _qZ);                       // D·Q·D = +Q
        double dMirror = MaxAbsDiff(_R * _qZ * _R, _qZ.Multiply(-Complex.One));        // R·Q·R = −Q
        int truly = TrulyCellCount(N);                                                 // n_Y even ∧ n_Z even
        long total = 1L << (2 * N);                                                    // 4^N Pauli strings
        return new InspectableNode("the three readings (the mirror group D₄, within one diagonal)",
            summary: "one diagonal, moved three ways by D₄ (PROOF_PI_FACTORS_AS_R_TIMES_D).",
            children: new IInspectable[]
            {
                new InspectableNode("rate: D fixes Q (D·Q·D = +Q)",
                    summary: $"the transpose fixes the diagonal — the rate is what the diagonal says; dev = {dRate:0.0e+00}"),
                new InspectableNode("mirror: R reflects Q (R·Q·R = −Q)",
                    summary: $"the ket-flip anti-fixes Q, carrying the −2Σγ palindrome shift; dev = {dMirror:0.0e+00}"),
                new InspectableNode("judge: the {D, 𝓕D} joint-fixed cell = truly",
                    summary: $"truly (n_Y even ∧ n_Z even) = {truly}/{total} Pauli strings at N={N} "
                           + "(the F87 trichotomy cell; 63/63 N=3 in mirror_inventory_d4.py block D)"),
            });
    }

    // ---------------------------------------------------------------- node 3: the orbit (basis-S₃)
    private IInspectable TheOrbit()
    {
        double dZX = MaxAbsDiff(_hZX * _qZ * _hZX.ConjugateTranspose(), _qX);
        double dYZ = MaxAbsDiff(_hYZ * _qZ * _hYZ.ConjugateTranspose(), _qY);
        int orbit = OrbitSize(_qZ, new[] { _hZX, _hYZ });
        bool sameSpec = SpectraEqual(_qX, _qY) && SpectraEqual(_qX, _qZ);
        return new InspectableNode("the orbit ({Q_X, Q_Y, Q_Z} one orbit, the three diagonals)",
            summary: $"orbit(Q_Z) under the basis-change S₃ ⟨h_zx, h_yz⟩ = {orbit} (=3 ⟹ the three diagonals are one orbit); "
                   + $"same spectrum = {sameSpec}; the structure is S₃ ⋉ D₄.",
            children: new IInspectable[]
            {
                new InspectableNode("h_zx : Q_Z → Q_X (the Z↔X Hadamard move)", summary: $"dev = {dZX:0.0e+00}"),
                new InspectableNode("h_yz : Q_Z → Q_Y (the Z↔Y R_x(π/2) move)", summary: $"dev = {dYZ:0.0e+00}"),
                new InspectableNode("same spectrum (conjugate ⟹ co-spectral)",
                    summary: $"spec(Q_X) = spec(Q_Y) = spec(Q_Z): {sameSpec}"),
            });
    }

    // ---------------------------------------------------------------- node 4: the dynamics (L_H even-step)
    private IInspectable TheDynamics()
    {
        var steps = RungSteps(_LH);                       // the set of |Δk| L_H connects
        bool parity = ParityConserved();
        int kBandEdge = _rung[0 * D + 1];                 // |vac⟩⟨1-magnon|: i=0, j=1
        int kBackwards = MirrorBackwardsWorst();          // max | rung(R m) − (N − rung(m)) |
        var (dFlipLH, rFixLH) = MirrorActionOnLH();       // §3: D·L_H·D = −L_H, R·L_H·R = +L_H
        return new InspectableNode("the dynamics (L_H is an even-step rung ladder — the complement)",
            summary: $"L_H connects rungs by Δk ∈ {{{string.Join(",", steps)}}} (even only); disagreement-count "
                   + $"parity {(parity ? "CONSERVED" : "BROKEN")}; band edge k={kBandEdge} (odd) ⟂ {{0,2}} (even) ⟹ level crossing.",
            children: new IInspectable[]
            {
                new InspectableNode("even-step ladder: L_H : k ↔ {k, k±2}",
                    summary: $"a hop flips two bits ⟹ Δk even; rung steps = {{{string.Join(",", steps)}}} (no k±1)"),
                new InspectableNode("disagreement-count parity conserved",
                    summary: $"L = L_H + L_D is block-diagonal in k mod 2: even↔odd block max = "
                           + $"{EvenOddBlockMax():0.0e+00} (U(1) feature; a transverse field would break it)"),
                new InspectableNode("the sector split (why the handover is a level crossing)",
                    summary: $"band edge |vac⟩⟨magnon| at k={kBandEdge} (odd); the {{0,2}} survivor even — "
                           + "parity-protected from hybridizing ⟹ a level crossing, not an EP coalescence"),
                new InspectableNode("the mirror reads the rungs backwards: R : k → N−k",
                    summary: $"max | rung(R·m) − (N − rung(m)) | = {kBackwards} (0 ⟹ the palindrome IS the "
                           + "disagreement count read from the other end)"),
                new InspectableNode("the mirror group acts on L_H (§3 of PROOF_PI_FACTORS — the bridge to the symmetry)",
                    summary: $"D flips it (D·L_H·D = −L_H, dev {dFlipLH:0.0e+00}); R fixes it (R·L_H·R = +L_H, dev "
                           + $"{rFixLH:0.0e+00}). With the rate/mirror readings (the L_diss column) this completes the FULL "
                           + "§3 2×2 palindrome split, live: the SAME L_H seen by the mirror group (D-flip/R-fix) and by the "
                           + "rung ladder (even-step) — the symmetry and the dynamics, welded."),
            });
    }

    // ============================ builders (reused conventions) ============================
    private static ComplexMatrix DephasingDiagonal(int n, PauliLetter letter)
    {
        int d = 1 << n, d2 = d * d;
        var acc = Matrix<Complex>.Build.Dense(d2, d2);
        for (int l = 0; l < n; l++)
        {
            var p = PauliString.SiteOp(n, l, letter);
            acc += p.KroneckerProduct(p.Transpose());   // Q_P = Σ_l kron(P_l, P_lᵀ)
        }
        return acc;
    }

    private static ComplexMatrix TransposeSuper(int n)
    {
        int d = 1 << n, d2 = d * d;
        var t = Matrix<Complex>.Build.Dense(d2, d2);
        for (int i = 0; i < d; i++)
            for (int jj = 0; jj < d; jj++)
                t[jj * d + i, i * d + jj] = Complex.One;   // vec(ρᵀ) = D·vec(ρ)
        return t;
    }

    private static ComplexMatrix KetReflection(int n)
    {
        int d = 1 << n;
        var idH = Matrix<Complex>.Build.DenseIdentity(d);
        var letters = new PauliLetter[n];
        for (int s = 0; s < n; s++) letters[s] = PauliLetter.X;
        var f = PauliString.Build(letters);                 // X^⊗N
        return idH.KroneckerProduct(f);                     // R: ρ ↦ ρ·F
    }

    private static ComplexMatrix AdUnitary(ComplexMatrix u1, int n)
    {
        var u = u1;
        for (int s = 1; s < n; s++) u = u.KroneckerProduct(u1);   // U = u1^⊗N
        return u.KroneckerProduct(u.Conjugate());                  // Ad_U: vec(UρU†) = (U⊗U*)vec(ρ)
    }

    private static ComplexMatrix Hadamard1()
    {
        double s = 1.0 / Math.Sqrt(2.0);
        return Matrix<Complex>.Build.DenseOfArray(new Complex[,] { { s, s }, { s, -s } });
    }

    private static ComplexMatrix RxHalfPi1()
    {
        double c = Math.Cos(Math.PI / 4.0), sn = Math.Sin(Math.PI / 4.0);
        return Matrix<Complex>.Build.DenseOfArray(new Complex[,]
            { { c, new Complex(0, -sn) }, { new Complex(0, -sn), c } });
    }

    private static ComplexMatrix BuildChainH(int n, double j, double delta)
    {
        int d = 1 << n;
        var h = Matrix<Complex>.Build.Dense(d, d);
        for (int b = 0; b < n - 1; b++)
        {
            var xx = PauliString.SiteOp(n, b, PauliLetter.X) * PauliString.SiteOp(n, b + 1, PauliLetter.X);
            var yy = PauliString.SiteOp(n, b, PauliLetter.Y) * PauliString.SiteOp(n, b + 1, PauliLetter.Y);
            var zz = PauliString.SiteOp(n, b, PauliLetter.Z) * PauliString.SiteOp(n, b + 1, PauliLetter.Z);
            h += (xx + yy).Multiply(new Complex(j, 0)) + zz.Multiply(new Complex(delta, 0));
        }
        return h;
    }

    private static ComplexMatrix HamiltonianSuper(ComplexMatrix h, int n)
    {
        int d = 1 << n;
        var idH = Matrix<Complex>.Build.DenseIdentity(d);
        var lh = h.KroneckerProduct(idH) - idH.KroneckerProduct(h.Transpose());   // [H,·]
        return lh.Multiply(new Complex(0, -1));                                     // −i[H,·]
    }

    private static int[] RungLabels(int n)
    {
        int d = 1 << n;
        var r = new int[d * d];
        for (int i = 0; i < d; i++)
            for (int jj = 0; jj < d; jj++)
                r[i * d + jj] = System.Numerics.BitOperations.PopCount((uint)(i ^ jj));
        return r;
    }

    // ============================ live readouts ============================
    private List<int> RungSteps(ComplexMatrix lh)
    {
        var set = new SortedSet<int>();
        for (int m = 0; m < lh.RowCount; m++)
            for (int nn = 0; nn < lh.ColumnCount; nn++)
                if (lh[m, nn].Magnitude > 1e-12)
                    set.Add(Math.Abs(_rung[m] - _rung[nn]));
        return set.ToList();
    }

    private bool ParityConserved() => EvenOddBlockMax() < 1e-12;

    private double EvenOddBlockMax()
    {
        var ld = Matrix<Complex>.Build.Dense(Dim, Dim);
        for (int m = 0; m < Dim; m++) ld[m, m] = new Complex(-2.0 * Gamma * _rung[m], 0);   // L_D diagonal
        var l = _LH + ld;
        double worst = 0.0;
        for (int m = 0; m < Dim; m++)
            for (int nn = 0; nn < Dim; nn++)
                if (((_rung[m] - _rung[nn]) & 1) == 1)          // even↔odd entry
                    worst = Math.Max(worst, l[m, nn].Magnitude);
        return worst;
    }

    private int MirrorBackwardsWorst()
    {
        int worst = 0;
        for (int m = 0; m < Dim; m++)
        {
            int mp = 0; double best = -1;
            for (int row = 0; row < Dim; row++)
                if (_R[row, m].Magnitude > best) { best = _R[row, m].Magnitude; mp = row; }
            worst = Math.Max(worst, Math.Abs(_rung[mp] - (N - _rung[m])));
        }
        return worst;
    }

    private int OrbitSize(ComplexMatrix q, ComplexMatrix[] gens)
    {
        var orbit = new List<ComplexMatrix> { q };
        bool changed = true;
        while (changed)
        {
            changed = false;
            foreach (var g in gens)
                foreach (var o in orbit.ToList())
                {
                    var qg = g * o * g.ConjugateTranspose();
                    if (!orbit.Any(x => MaxAbsDiff(qg, x) <= Tol)) { orbit.Add(qg); changed = true; }
                }
        }
        return orbit.Count;
    }

    private static int TrulyCellCount(int n)
    {
        int count = 0;
        long total = 1L << (2 * n);
        for (long s = 0; s < total; s++)
        {
            int nY = 0, nZ = 0;
            for (int site = 0; site < n; site++)
            {
                int letter = (int)((s >> (2 * site)) & 3);   // 0=I,1=X,2=Y,3=Z
                if (letter == 2) nY++;
                if (letter == 3) nZ++;
            }
            if ((nY & 1) == 0 && (nZ & 1) == 0) count++;
        }
        return count;
    }

    private static long Binom(int n, int k)
    {
        if (k < 0 || k > n) return 0;
        long r = 1;
        for (int i = 0; i < k; i++) r = r * (n - i) / (i + 1);
        return r;
    }

    private static bool SpectraEqual(ComplexMatrix a, ComplexMatrix b)
    {
        var ea = a.Evd().EigenValues.Enumerate().Select(z => z.Real).OrderBy(x => x).ToArray();
        var eb = b.Evd().EigenValues.Enumerate().Select(z => z.Real).OrderBy(x => x).ToArray();
        return ea.Length == eb.Length && ea.Zip(eb, (x, y) => Math.Abs(x - y)).All(v => v <= 1e-9);
    }

    private static double MaxAbsDiff(ComplexMatrix a, ComplexMatrix b)
    {
        double m = 0.0;
        for (int i = 0; i < a.RowCount; i++)
            for (int jj = 0; jj < a.ColumnCount; jj++)
            {
                double v = (a[i, jj] - b[i, jj]).Magnitude;
                if (v > m) m = v;
            }
        return m;
    }
}
