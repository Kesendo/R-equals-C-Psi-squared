using System.Globalization;
using MirrorWorld;

CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;   // dots, not commas

const double gamma = 0.5;
var world = new World();
Console.WriteLine($"empty world (Z-dephasing, no Hamiltonian), gamma={gamma}");
Console.WriteLine($"World.Own (left, the frame): {string.Join(",", world.Own)}");
Console.WriteLine();

foreach (int N in new[] { 2, 3, 4 })
{
    var modes = EnumeratePauli(N).Select(l => new PauliMode(world, l, gamma)).ToList();
    Console.WriteLine($"==== N={N}  ({modes.Count} Pauli-string modes = {1 << N}*sum C(N,k); the superposition basis) ====");

    for (int k = 0; k <= N; k++)
    {
        int count = modes.Count(m => m.K == k);
        string tag = k == 0 ? "  [populations / kernel]"
                   : k == N ? "  [drain]"
                   : k == N - k ? "  [SELF-MIRROR center]" : "";
        Console.WriteLine($"  k={k}: {count,3} modes   rate {-2.0 * gamma * k,5:0.0}   mirror k={N - k}{tag}");
    }

    Console.WriteLine(N % 2 == 0
        ? $"  fold k<->N-k about k=N/2={N / 2}: even N, the self-mirror sector EXISTS (k={N / 2})"
        : $"  fold k<->N-k about k=N/2={(N / 2.0):0.0}: odd N, axis BETWEEN rungs, NO self-mirror sector");

    int kc = (N % 2 == 0) ? N / 2 : (N - 1) / 2;   // central (even) or inner (odd)
    var cells = modes.Where(m => m.K == kc)
                     .GroupBy(m => m.Klein)
                     .OrderBy(g => g.Key.Y).ThenBy(g => g.Key.Z)
                     .Select(g => $"(nY%2={g.Key.Y},nZ%2={g.Key.Z}):{g.Count()}");
    Console.WriteLine($"  superposition / Klein quartering of sector k={kc}: {string.Join("  ", cells)}   [truly = (0,0)]");

    // Grading B (T1): the joint-popcount blocks a number-conserving Hamiltonian respects.
    Console.WriteLine($"  Grading B: {(N + 1) * (N + 1)} joint-popcount blocks (p,q), size C(N,p)*C(N,q):");
    long total = 0;
    for (int p = 0; p <= N; p++)
    {
        var sizes = new List<string>();
        for (int q = 0; q <= N; q++)
        {
            var b = new Block(world, N, p, q);
            total += b.Size;
            sizes.Add($"{b.Size,3}");
        }
        Console.WriteLine($"    p={p}: {string.Join(" ", sizes)}");
    }
    Console.WriteLine($"    sum = {total} = 4^{N}; the diagonal (p,p) carries k=0, the {N + 1} populations (the kernel)");
    Console.WriteLine();
}

// The special cases as T1 closed forms , adopted, not re-derived (all already Tier-1 in the repo).
Console.WriteLine("==== special cases (T1 closed forms, adopted) ====");
for (int n = 2; n <= 6; n++)
    Console.WriteLine($"  N={n}: V-effect 1+cos(pi/N) = {1.0 + Math.Cos(Math.PI / n):0.0000}   ceiling g2(K_N)=4/N = {4.0 / n:0.0000}");
Console.WriteLine($"  N=3: minimal EP, Q_EP = 2/g_eff (the whole rate ladder is HD=1 and HD=3)");
Console.WriteLine($"  N=4: self-fold diabolic q_EP = sqrt((-1+sqrt13)/6) = {Math.Sqrt((-1.0 + Math.Sqrt(13.0)) / 6.0):0.0000}  (real q, N=4 only); (2,2) ceiling K_4 = 2-2/sqrt3 = {2.0 - 2.0 / Math.Sqrt(3.0):0.0000}");
Console.WriteLine($"  N=5: diabolics leave the real axis into complex-conjugate pairs (only loud defective EPs remain real)");

// H on: the superposition leaves the integer grid (the on-grid d_total folds are T1, adopted).
Console.WriteLine();
Console.WriteLine("==== H on: the superposition leaves the grid (T1 d_total folds, adopted) ====");
foreach (int n in new[] { 2, 3, 4 })
{
    var bare = Redistribution.Bare(n);
    var grid = Redistribution.OnGrid(n)!;
    var diff = Enumerable.Range(0, n + 1).Select(k => grid[k] - bare[k]).ToArray();
    int off = (1 << (2 * n)) - grid.Sum();
    Console.WriteLine($"  N={n}: bare [{string.Join(",", bare)}]  ->  on-grid [{string.Join(",", grid)}]   diff [{string.Join(",", diff)}]   off-grid {off}");
}
Console.WriteLine("  edges stay N+1 (kernel + drain); even N spikes the center k=N/2; the rest bleeds off-grid (fractional <n_XY>)");

// The clock: theta = arctan(Q), Q = J/gamma (T1, F95). The two hands of one winding mode.
Console.WriteLine();
Console.WriteLine("==== the clock: theta = arctan(Q), Q = J/gamma (T1) ====");
foreach (double q in new[] { 0.0, 0.5, 1.0, 2.0, 5.0, 100.0 })
{
    var clk = new Clock(world, j: q, gamma: 1.0);   // J=q, gamma=1 => Q=q
    string mark = q == 0.0 ? "  [J=0: pure decay]" : Math.Abs(q - 1.0) < 1e-9 ? "  [theta=45deg <-> 1/4]" : q >= 100.0 ? "  [carbon, deep quantum]" : "";
    Console.WriteLine($"  Q={clk.Q,5:0.0}: theta = {clk.ThetaDeg,5:0.0} deg{mark}");
}
Console.WriteLine("  gamma=0 -> theta=90deg (pure circle, turning forever); 90deg <-> 1/2; the two hands are gamma (radial) and J (angular)");

// The survivor: the slowest non-stationary mode (T1), regime-dependent.
Console.WriteLine();
Console.WriteLine("==== the survivor: the slowest non-stationary mode (T1) ====");
foreach (int n in new[] { 2, 3, 4, 5 })
{
    var s = new Survivor(world, n);
    string lo = s.HasHalfFillingSurvivor
        ? $"low-Q: half-filling k={n / 2}, R-odd/X-odd, dark"
        : "low-Q: no (N/2,N/2) sector (odd N), no half-filling survivor";
    Console.WriteLine($"  N={n}: {lo};  hands over at Q*={s.Qstar:0.00} to the (0,1) band edge (<n_XY>=1, rate -2g)");
}

// R-parity and mod-4 (T1): the even/odd EP behaviour and the N=3-mod-4 amplification.
Console.WriteLine();
Console.WriteLine("==== R-parity & mod-4 (T1) ====");
foreach (int n in new[] { 3, 4, 5, 6 })
{
    bool epSplits = n % 2 == 1;     // odd N: EP splits by R-parity; even N: isospectral
    bool amplifies = n % 4 == 3;    // self-Pi SE mode R-odd at N=3 mod 4 -> amplifies
    Console.WriteLine($"  N={n}: {(epSplits ? "EP SPLITS by R-parity (sigma+ != sigma-)" : "EP isospectral (sigma+ = sigma-, no split)")};  " +
                      $"mod-4: {(amplifies ? "self-Pi SE mode R-odd, amplifies (~2.2x)" : "non-amplifying")}");
}
Console.WriteLine("  even N: PH partners (k,N+1-k) carry opposite R-parity -> isospectral; odd N: same R-parity + a lone R-fixed zero mode -> split");

// F-registry closed forms, adopted verbatim (J=1).
Console.WriteLine();
Console.WriteLine("==== F-registry closed forms (adopted verbatim, J=1) ====");
foreach (int n in new[] { 3, 4, 5 })
{
    Console.WriteLine($"  N={n}: F2  omega_k=4J(1-cos(pi k/N))   = [{string.Join(", ", Formulas.F2_Dispersion(n, 1.0).Select(x => x.ToString("0.00")))}]");
    Console.WriteLine($"        F2b E_k=2J cos(pi k/(N+1))       = [{string.Join(", ", Formulas.F2b_SingleExcitation(n, 1.0).Select(x => x.ToString("0.00")))}]");
}
Console.WriteLine($"  F1 residuals (N=4, gamma=0.5: sg=2, sg2=1): ||M(T1)||^2 = {Formulas.F1_T1Residual(4, 2.0, 1.0):0.0}, ||M(depol)||^2 = {Formulas.F1_DepolResidual(4, 2.0, 1.0):0.0}");
Console.WriteLine($"  F2b-corollary coherence hand omega_mem=2J cos(pi/(N+1)): N=3,4,5 = {Formulas.OmegaMem(3, 1, 0):0.000}, {Formulas.OmegaMem(4, 1, 0):0.000}, {Formulas.OmegaMem(5, 1, 0):0.000}  (sqrt2, phi, sqrt3)");
Console.WriteLine($"  coherence horizon Q*(N) exact: N=2..5 = {Formulas.Qstar(2):0.000}, {Formulas.Qstar(3):0.000}, {Formulas.Qstar(4):0.000}, {Formulas.Qstar(5):0.000}  (asymptotic 2N/pi)");
foreach (int n in new[] { 3, 4, 5 })
{
    var (mn, mx, bw) = Formulas.F3_RateBounds(n, 0.5);
    Console.WriteLine($"  N={n}: F3 min={mn:0.0}/max={mx:0.0}/bw={bw:0.0}  F4 kernel={Formulas.F4_KernelDim(n)}  F5 depol-err={Formulas.F5_DepolError(n, 0.5):0.000}  F23 XOR-frac={Formulas.F23_XorFraction(n):0.0000}  F50 w1-deg={Formulas.F50_Weight1Degeneracy(n)}");
}
Console.WriteLine($"  F33 N=3 exact rates (g=0.5) = [{string.Join(", ", Formulas.F33_N3Rates(0.5).Select(x => x.ToString("0.000")))}]  (<n_XY>=1, 4/3, 5/3)");
foreach (int n in new[] { 3, 4, 5 })
{
    var (qmax, qmin, qmean, qspread) = Formulas.F7_QSpectrum(n, 1.0, 0.5);
    Console.WriteLine($"  N={n}: F7 Q max/min/mean/spread = {qmax:0.00}/{qmin:0.00}/{qmean:0.0}/{qspread:0.00}");
}
Console.WriteLine($"  F8 2x law (N=4,g=0.5): unpaired={Formulas.F8_DecayLaw(4, 0.5).Unpaired:0.0}, paired-mean={Formulas.F8_DecayLaw(4, 0.5).PairedMean:0.0}, ratio=2");
Console.WriteLine($"  F12 single-qubit crossing t*/T2 = {Formulas.F12_CrossingFraction} (root of x^3+x=1/2)");
Console.WriteLine($"  F16 fold R=C(Psi+R)^2, boundary CPsi={Formulas.F16_FoldBoundary} (Mandelbrot u->u^2+c)");
Console.WriteLine($"  F25 CPsi Bell+ Z-deph: crossing f*={Formulas.F25_CrossingF}, K_Z={Formulas.F25_K};  F27 K_X=K_Y=ln2/8={Formulas.F27_KX:0.0000}, K_depol={Formulas.F27_KDepol}");
Console.WriteLine($"  F15 theta compass at CPsi=0.5: {Formulas.F15_ThetaDeg(0.5):0.0} deg (0 at the 1/4 crossing)");

// all 4^N Pauli strings, site 0 the least-significant base-4 digit
static IEnumerable<char[]> EnumeratePauli(int N)
{
    char[] alphabet = { 'I', 'X', 'Y', 'Z' };
    int total = 1 << (2 * N);
    for (int idx = 0; idx < total; idx++)
    {
        var l = new char[N];
        int x = idx;
        for (int s = 0; s < N; s++) { l[s] = alphabet[x & 3]; x >>= 2; }
        yield return l;
    }
}
