using System.Globalization;
using MirrorWorld;

CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;   // dots, not commas

// ---- run mode "grow": the diagonal protocol, step 1 (the world splits) ----
// ClaudeTasks/DIAGONAL_PROTOCOL_GAME.md, rules 1-3 running in time. No Hamiltonian yet: a field of
// possibilities under the one question, structure (the diagonal) staying while novelty (off-diagonal)
// is culled by its disagreement. The mirror is visible in the rate ladder (k decays, N-k its twin).
if (args.Length > 0 && args[0] == "grow")
{
    const int gn = 2;
    const double gg = 0.5, dt = 0.05;
    const int ticks = 20;
    var gworld = new World();
    var field = new Field(gworld, gn, gg);
    field.SeedUniform();

    Console.WriteLine("the world splits (diagonal protocol, step 1 -- rules 1-3 in time)");
    Console.WriteLine($"  source ClaudeTasks/DIAGONAL_PROTOCOL_GAME.md; N={gn}, gamma={gg}, dt={dt}");
    Console.WriteLine("  rule 2, the one question: how much do two possibilities disagree? more disagreement, faster fade.");
    Console.WriteLine("  structure = the diagonal (k=0) that stays; novelty = the off-diagonal that is culled.");
    Console.WriteLine($"  the knower's cut: {field.AliveCount} alive (upper triangle, stepped), {field.ImmortalCount} immortal (k=0, held), {field.AliveCount} mirrored (rho=rho-dagger) -- only {field.AliveCount} of {field.Dim * field.Dim} cells run.");
    Console.WriteLine($"  {"t",4} {"structure",9} {"novelty",9}   {string.Join(" ", Enumerable.Range(0, gn + 1).Select(k => $"k={k}".PadLeft(7)))}");
    for (int tick = 0; tick <= ticks; tick++)
    {
        var byK = field.WeightByDisagreement();
        Console.WriteLine($"  {field.T,4:0.00} {field.Structure,9:0.000} {field.Novelty,9:0.000}   {string.Join(" ", byK.Select(v => v.ToString("0.000").PadLeft(7)))}");
        field.Step(dt);
    }
    Console.WriteLine("  structure held flat (k=0 costs nothing); novelty fell, k=2 twice as fast as k=1 -- the mirror in the rates.");
    Console.WriteLine("  (step 2 -- run `-- live` -- adds the inner restlessness: novelty BORN from structure.)");
    return;
}

// ---- run mode "live": the diagonal protocol, step 2 (the world lives) ----
// The empty world's watching plus the inner restlessness: H (a flip-flop on the bond) makes coherence out
// of population -- novelty BORN from structure (rule 4) -- while the watching culls it. Novelty no longer
// only fades; it is born and culled, and the populations relax to a living balance.
if (args.Length > 0 && args[0] == "live")
{
    int gn = args.Length > 1 ? int.Parse(args[1]) : 2;
    string gtopo = args.Length > 2 ? args[2] : "chain";
    const double gj = 1.0, gg = 0.5, dt = 0.05;
    const int ticks = 40;
    var lworld = new World();
    var rest = new Restless(lworld, gn, gj, gg, Topology.Named(gtopo, gn));
    rest.Seed(1);                                   // |0...01>, one excitation at one end (block (1,1))

    Console.WriteLine("the world lives (diagonal protocol, step 2 -- the inner restlessness)");
    Console.WriteLine($"  source ClaudeTasks/DIAGONAL_PROTOCOL_GAME.md; N={gn}, topology={gtopo}, J={gj}, gamma={gg}, dt={dt}");
    Console.WriteLine("  rule 4, the inner motion: H (the flip-flop bond) makes coherence out of population --");
    Console.WriteLine("  novelty BORN from structure, then culled by the watching. seed = |0...01>, one excitation.");
    Console.WriteLine($"  cut (c): the loop stays in joint-popcount block (1,1) -- {rest.AliveCount} of {rest.Dim * rest.Dim} cells run, {rest.ForbiddenCount} forbidden (F63). the k=1 coherences are among them, never computed.");
    Console.WriteLine($"  {"t",4} {"structure",9} {"novelty",9}   {string.Join(" ", Enumerable.Range(0, gn + 1).Select(k => $"k={k}".PadLeft(7)))}");
    for (int tick = 0; tick <= ticks; tick++)
    {
        var byK = rest.WeightByDisagreement();
        Console.WriteLine($"  {rest.T,4:0.00} {rest.Structure,9:0.000} {rest.Novelty,9:0.000}   {string.Join(" ", byK.Select(v => v.ToString("0.000").PadLeft(7)))}");
        rest.Step(dt);
    }
    Console.WriteLine("  novelty rose from 0 (born by H from the population), churned, then damped -- the living");
    Console.WriteLine("  balance: the excitation spreads over the chain toward 1/N per site (structure persists, redistributed).");
    return;
}

// ---- run mode "scale": the complexity wall and the cut ----
// Why the full-spectrum side hit the wall and a state's dynamics did not. full = 4^N (the whole Liouvillian
// an eigendecomposition tackles); single-exc = the (1,1) block a one-excitation state lives in (N^2,
// polynomial); half-fill = the largest block C(N,N/2)^2 (the spectrum's irreducible core, still exponential).
if (args.Length > 0 && args[0] == "scale")
{
    var sworld = new World();
    Console.WriteLine("the complexity wall and the cut (why the full spectrum hit it, a state's dynamics did not)");
    Console.WriteLine($"  {"N",3} {"full 4^N",13} {"single-exc",11} {"forbidden",10} {"half-fill",11} {"hf forbidden",13}");
    for (int n = 2; n <= 12; n++)
    {
        long full = 1L << (2 * n);
        long se = new Block(sworld, n, 1, 1).Size;          // N^2, the one-excitation block
        long hf = new Block(sworld, n, n / 2, n / 2).Size;  // C(N,N/2)^2, the largest block
        Console.WriteLine($"  {n,3} {full,13} {se,11} {100.0 * (full - se) / full,9:0.00}% {hf,11} {100.0 * (full - hf) / full,11:0.000}%");
    }
    Console.WriteLine("  single-exc: N^2 -- a one-excitation state runs in a polynomial corner, tractable at any N.");
    Console.WriteLine("  half-fill: ~4^N/sqrt(N) -- the full spectrum needs it, so the spectrum stays exponential.");
    Console.WriteLine("  the wall was global (all 4^N, the spectrum); a single state's dynamics is block-local.");
    return;
}

// ---- run mode "topo": play with the geometry (who can reach whom) ----
// Same handshake and seed, different bonds. All excitation-conserving, so all share the block cut (F63);
// only the dynamics -- how fast novelty is born and spreads -- differs.
if (args.Length > 0 && args[0] == "topo")
{
    int gn = args.Length > 1 ? int.Parse(args[1]) : 4;
    const double gj = 1.0, gg = 0.5, dt = 0.05;
    var tworld = new World();
    Console.WriteLine($"playing with the geometry (who can reach whom): N={gn}, seed |0...01>, J={gj}, gamma={gg}");
    Console.WriteLine("  same handshake, different bonds. all excitation-conserving -> the same block cut (F63);");
    Console.WriteLine("  only the dynamics differs. (seed site 0 is the chain's endpoint, but the star's hub.)");
    Console.WriteLine($"  {"topology",9} {"bonds",6} {"alive",6} {"forbidden",10}   novelty at t = 0.25 / 0.50 / 1.00");
    foreach (var name in new[] { "chain", "ring", "star", "complete" })
    {
        var rest = new Restless(tworld, gn, gj, gg, Topology.Named(name, gn));
        rest.Seed(1);
        var snap = new double[3];
        int[] marks = { 5, 10, 20 };
        int mi = 0;
        for (int tick = 0; tick <= 20; tick++)
        {
            if (mi < marks.Length && tick == marks[mi]) snap[mi++] = rest.Novelty;
            rest.Step(dt);
        }
        Console.WriteLine($"  {name,9} {Topology.Named(name, gn).Length,6} {rest.AliveCount,6} {rest.ForbiddenCount,10}   {snap[0],6:0.000} / {snap[1],6:0.000} / {snap[2],6:0.000}");
    }
    Console.WriteLine("  alive is identical across all four -- the cut is topology-invariant. the geometry only");
    Console.WriteLine("  steers the dynamics: more bonds from the seed (star hub, complete) -> novelty born faster.");
    return;
}

// ---- run mode "regime": how alive at each Q = J/gamma (the clock) ----
// The same world at different ratios of restlessness to watching. Low Q: the watching wins (overdamped,
// novelty barely born and quickly drained -- near death). High Q: the restlessness wins (underdamped,
// novelty churns on and on). The clock angle theta = arctan(Q) reads the regime (reused from Clock).
if (args.Length > 0 && args[0] == "regime")
{
    int gn = args.Length > 1 ? int.Parse(args[1]) : 4;
    const double gg = 0.5, dt = 0.05;
    const int ticks = 80;                            // to t = 4
    var rworld = new World();
    Console.WriteLine($"how alive at each regime Q = J/gamma (N={gn} chain, seed |0...01>, gamma={gg})");
    Console.WriteLine("  the clock: theta = arctan(Q). low Q = the watching wins (overdamped, ~death);");
    Console.WriteLine("  high Q = the restlessness wins (underdamped, novelty churns on and on).");
    Console.WriteLine($"  {"Q",5} {"theta",7} {"peak nov",9} {"t@peak",7} {"late churn (t3-4)",17}");
    foreach (double q in new[] { 0.2, 0.5, 1.0, 2.0, 5.0, 20.0 })
    {
        double j = q * gg;
        var clk = new Clock(rworld, j, gg);          // theta = arctan(Q), Q = J/gamma (the adopted clock)
        var rest = new Restless(rworld, gn, j, gg, Topology.Chain(gn));
        rest.Seed(1);
        double peak = 0, tpeak = 0, lateMin = double.MaxValue, lateMax = 0;
        for (int tick = 0; tick <= ticks; tick++)
        {
            double nov = rest.Novelty;
            if (nov > peak) { peak = nov; tpeak = rest.T; }
            if (tick >= ticks - 20) { lateMin = Math.Min(lateMin, nov); lateMax = Math.Max(lateMax, nov); }
            rest.Step(dt);
        }
        Console.WriteLine($"  {clk.Q,5:0.0} {clk.ThetaDeg,6:0.0} {peak,9:0.000} {tpeak,7:0.00} {lateMax - lateMin,17:0.000}");
    }
    Console.WriteLine("  peak rises and t@peak falls, both monotonic -- more restlessness births more novelty, sooner.");
    Console.WriteLine("  late churn is ~0 below Q=1 (overdamped: a slow bump that just decays) and switches on above it");
    Console.WriteLine("  (underdamped: novelty oscillates). It is not monotonic at high Q -- the chain is multi-mode, so");
    Console.WriteLine("  the underdamped churn beats; a single scalar cannot hold it. The clean turn is Q=1, theta=45deg.");
    return;
}

const double gamma = 0.5;
var world = new World();
Console.WriteLine($"empty world (Z-dephasing, no Hamiltonian), gamma={gamma}");
Console.WriteLine($"World.Own (left, the frame): {string.Join(",", world.Own)}");
Console.WriteLine();

foreach (int N in new[] { 2, 3, 4 })
{
    var modes = PauliMode.Enumerate(world, N, gamma).ToList();
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
Console.WriteLine($"  F34 qubit necessity: d^2-2d=0 -> d in {{{string.Join(", ", Formulas.F34_QubitNecessity())}}} (d=0 nothing, d=2 the qubit, the polarity root)");
Console.WriteLine($"  D1 BW=8J cos(pi/N) (N=4,J=1) = {Formulas.D1_Bandwidth(4, 1.0):0.000} -> 8J at large N");
Console.WriteLine($"  D4 crossing rhs (d-1)/2: d=2 -> {Formulas.D4_CrossingRhs(2):0.0}, d=4 -> {Formulas.D4_CrossingRhs(4):0.0}");
Console.WriteLine($"  D6 gap=2g={Formulas.D6_Gap(0.5):0.0}, mixing time (N=4) <= {Formulas.D6_MixingTime(4, 0.5):0.00};  F38 Pi^2=(-1)^(nY+nZ)=X^N (order 4)");
Console.WriteLine($"  F18 fold threshold Sg_crit/J: Bell={Formulas.F18_FoldThresholdBell}, product={Formulas.F18_FoldThresholdProduct}");
Console.WriteLine($"  F36/F37 neural (Wilson-Cowan): Q J Q + J + 2S = 0, pairing mu+mu'=-(1/tE+1/tI); C.elegans 0.013 vs 0.108 random");
Console.WriteLine($"  F61 bit_a parity Pi^2_X=Z^N=(-1)^(nX+nY)=(-1)^k (the k-parity H conserves); F63 [L,Pi^2]=0 (conserved, all N)");
Console.WriteLine($"  F39 det(Pi)=(-1)^(N 4^(N-1)): N=1->{Formulas.F39_DetPi(1)}, N>=2->{Formulas.F39_DetPi(4)};  F41 t_Pi=pi/(4J sin^2(pi/2N)) (N=4) = {Formulas.F41_PalindromicTime(4, 1.0):0.000}");
Console.WriteLine($"  F49 cross-term R(N)=sqrt((N-2)/(N 4^(N-1))): N=2,3,4 = {Formulas.F49_CrossTerm(2):0.000}, {Formulas.F49_CrossTerm(3):0.0000}, {Formulas.F49_CrossTerm(4):0.0000} (0, 1/sqrt48, 1/sqrt128);  F49b ||L_Dc||^2=g^2 4^N N (N=4,g=0.5) = {Formulas.F49b_CenteredDissipatorNormSq(4, 0.5):0.0}");
Console.WriteLine($"  F49c shadow-crossing R(N)=sqrt((N-1)/(N 4^(N-1))): N=3,4 = {Formulas.F49c_CrossTermCrossing(3):0.0000}, {Formulas.F49c_CrossTermCrossing(4):0.0000}");
Console.WriteLine($"  F55 K_death=ln10={Formulas.F55_KDeath:0.000} (99% absorption); immortal modes=N+1={Formulas.F55_ImmortalModes(4)} (N=4);  F56 cusp K(eps=1e-4,tol=1e-12) = {Formulas.F56_CriticalSlowing(1e-4, 1e-12):0.000}");
Console.WriteLine($"  F57 Bell+ dwell prefactor K_dwell/delta = {Formulas.F57_DwellPrefactorBell} (gamma-independent)");
Console.WriteLine($"  F60 GHZ_N CPsi(0)=1/(2^N-1): N=2,3,4,5 = {Formulas.F60_GhzCPsi0(2):0.000}, {Formulas.F60_GhzCPsi0(3):0.000}, {Formulas.F60_GhzCPsi0(4):0.000}, {Formulas.F60_GhzCPsi0(5):0.000} (<1/4 for N>=3)");
Console.WriteLine($"  F62 W_N CPsi(0)=2(N^2-4N+8)/(3N^3): N=2,3,4 = {Formulas.F62_WstateCPsi0(2):0.000}, {Formulas.F62_WstateCPsi0(3):0.000}, {Formulas.F62_WstateCPsi0(4):0.000};  F59 prefactor(W0=1/2,k=2,W2=0.371) = {Formulas.F59_DwellPrefactor(2, 0.5, 0.3709):0.000}");
Console.WriteLine($"  F63 [L,Pi^2]=0, 4 blocks dim 4^(N-1)={Formulas.F63_BlockDim(4)} (N=4); conserved/sector (even,odd) N=4 = {Formulas.F63_ConservedPerSector(4)}");
Console.WriteLine($"  F65 SE spectrum (4/(N+1))sin^2(k pi/(N+1)): N=3 = [{string.Join(", ", Formulas.F65_SingleExcitationRates(3).Select(x => x.ToString("0.000")))}], N=5 = [{string.Join(", ", Formulas.F65_SingleExcitationRates(5).Select(x => x.ToString("0.000")))}]");
Console.WriteLine($"  F65 Niven: rational iff N in {{0,1,2,3,5}} (N=4 golden-irrational); F66 pole multiplicity N+1 = {Formulas.F66_PoleMultiplicity(4)} (N=4)");
Console.WriteLine($"  F68 partner rate alpha_p=2g0-alpha_b (g0=0.05, alpha_b=0.0138) = {Formulas.F68_PartnerRate(0.0138, 0.05):0.0000};  F69 GHZ3+W3 optimum pair-CPsi(0) = {Formulas.F69_N3Optimum:0.000} (ratio {Formulas.F69_RatioToQuarter} to 1/4, irreducible sextic)");
Console.WriteLine($"  F70 k-local sees |dN|<=k (single-site<=1, pair<=2); F71 c1 components=floor(N/2)={Formulas.F71_C1IndependentComponents(5)} (N=5), R|psi_k>=(-1)^(k+1) [k=1->{Formulas.F71_ReflectionParity(1)}, k=2->{Formulas.F71_ReflectionParity(2)}]");
Console.WriteLine($"  F98 Dicke asymptote (N+2)/(4(N+1)): N=4,6,8 = {Formulas.F98_DickeAsymptote(4):0.000}, {Formulas.F98_DickeAsymptote(6):0.000}, {Formulas.F98_DickeAsymptote(8):0.000} (-> 1/4)");
Console.WriteLine($"  F121 qudit palindrome paired(d,N)/d^2N: d=2,N=2 = {Formulas.F121_PairedCeiling(2, 2)}/16 (full); d=3,N=2 = {Formulas.F121_PairedCeiling(3, 2)}/81 (partial, d^2-2d=0 re-seen)");
Console.WriteLine($"  F122 ceiling g2: K_5=4/N={Formulas.F122_CompleteCeiling(5):0.000}, star_6=4/(N-1)={Formulas.F122_StarCeiling(6):0.000}, K_4=2-2/sqrt3={Formulas.F122_K4Ceiling():0.000}, ring(1,1) N=5={Formulas.F122_RingCommutant(5):0.000}");

// (the 4^N Pauli-string enumeration now lives in PauliMode.Enumerate, shared with the tests)
