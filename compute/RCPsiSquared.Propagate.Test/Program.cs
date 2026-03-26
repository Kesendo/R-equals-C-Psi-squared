using System.Diagnostics;
using System.Globalization;
using System.Numerics;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex;
using RCPsiSquared.Propagate;

// ============================================================
// PROFILE MODE: lightweight single-evaluation, stdout only
// ============================================================
if (args.Length >= 1 && args[0] == "profile")
{
    try { Control.UseNativeMKL(); } catch { }
    Control.MaxDegreeOfParallelism = Environment.ProcessorCount;
    RunProfileEvaluation(args);
    return;
}

// ============================================================
// STAGED MODE: two-phase gamma profile switching
// ============================================================
if (args.Length >= 1 && args[0] == "staged")
{
    try { Control.UseNativeMKL(); } catch { }
    Control.MaxDegreeOfParallelism = Environment.ProcessorCount;
    RunStagedEvaluation(args);
    return;
}

// ============================================================
// SPECTRAL MODE: eigendecomposition + midpoint analysis
// ============================================================
if (args.Length >= 1 && args[0] == "spectral")
{
    try { Control.UseNativeMKL(); } catch { }
    Control.MaxDegreeOfParallelism = Environment.ProcessorCount;
    RunSpectralAnalysis(args);
    return;
}

// ============================================================
// WAVE MODE: per-pair MI to see energy flow through chain
// ============================================================
if (args.Length >= 1 && args[0] == "wave")
{
    try { Control.UseNativeMKL(); } catch { }
    Control.MaxDegreeOfParallelism = Environment.ProcessorCount;
    RunWaveAnalysis(args);
    return;
}

// ============================================================
// RESONANCE MODE: CΨ oscillation around ¼ with structured bath
// ============================================================
if (args.Length >= 1 && args[0] == "resonance")
{
    try { Control.UseNativeMKL(); } catch { }
    Control.MaxDegreeOfParallelism = Environment.ProcessorCount;
    RunResonanceTest(args);
    return;
}

// ============================================================
// SWEEP MODE: multi-phase gamma profile sweep
// ============================================================
if (args.Length >= 1 && args[0] == "sweep")
{
    try { Control.UseNativeMKL(); } catch { }
    Control.MaxDegreeOfParallelism = Environment.ProcessorCount;
    RunSweepEvaluation(args);
    return;
}

// ============================================================
// OUTPUT
// ============================================================
var outFileName = args.Contains("pull") ? "pull_principle.txt" : "mediator_bridge_scale.txt";
var outPath = Path.Combine(
    Path.GetDirectoryName(AppContext.BaseDirectory)!,
    "..", "..", "..", "..", "..", "simulations", "results", outFileName);
outPath = Path.GetFullPath(outPath);

if (!Directory.Exists(Path.GetDirectoryName(outPath)))
{
    outPath = Path.Combine(AppContext.BaseDirectory, outFileName);
}

using var outFile = new StreamWriter(outPath, false, System.Text.Encoding.UTF8);

void Log(string msg = "")
{
    Console.WriteLine(msg);
    outFile.WriteLine(msg);
    outFile.Flush();
}

// ============================================================
// MKL INIT
// ============================================================
try { Control.UseNativeMKL(); Log("MKL: active"); }
catch { Log("MKL: not available, using managed provider"); }
Control.MaxDegreeOfParallelism = Environment.ProcessorCount;
Log($"Threads: {Environment.ProcessorCount}");

var totalSw = Stopwatch.StartNew();
Log("The Mediator Bridge Scales (Level 3, N=11, Time Propagation)");
Log($"Started: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
Log();

// ============================================================
// HELPER: Build initial states
// ============================================================
Complex[] BellPair(int q0, int q1, int nQubits)
{
    int d = 1 << nQubits;
    var psi = new Complex[d];
    psi[0] = 1.0 / Math.Sqrt(2);
    psi[(1 << (nQubits - 1 - q0)) | (1 << (nQubits - 1 - q1))] = 1.0 / Math.Sqrt(2);
    return psi;
}

// ============================================================
// TEST 0: Validation (5-qubit Level 2, compare with Python results)
// ============================================================
void RunTest0()
{
    Log(new string('=', 70));
    Log("TEST 0: VALIDATION (Level 2, N=5)");
    Log(new string('=', 70));
    Log();

    int n = 5;
    int d = 1 << n;
    double gamma = 0.05;
    var gammas = Enumerable.Repeat(gamma, n).ToArray();

    var sw = Stopwatch.StartNew();
    Log("Building Hamiltonian...");
    var bonds = Topology.MediatorBridge(2);
    var H = Topology.BuildHamiltonian(n, bonds);
    Log($"H built: {sw.Elapsed.TotalSeconds:F1}s");

    var prop = new LindbladPropagator(H, gammas, n);

    // |psi> = (|00> + |11>)/sqrt(2) ⊗ |0> ⊗ |+> ⊗ |+>
    var psi = new Complex[d];
    for (int q3 = 0; q3 < 2; q3++)
        for (int q4 = 0; q4 < 2; q4++)
        {
            int idx00 = (0 << 4) | (0 << 3) | (0 << 2) | (q3 << 1) | q4;
            int idx11 = (1 << 4) | (1 << 3) | (0 << 2) | (q3 << 1) | q4;
            psi[idx00] += 1.0 / Math.Sqrt(2) * 0.5;  // 1/sqrt(2) * 1/2 for each +
            psi[idx11] += 1.0 / Math.Sqrt(2) * 0.5;
        }

    var rho0 = DensityMatrixTools.PureState(psi);

    // Propagate and measure at key times
    var tMeasure = new[] { 0.0, 1.0, 3.0, 5.0, 10.0, 20.0 };

    Log();
    Log($"  {"t",6}  {"MI(A:B)",10}  {"CPsi_24",10}  {"Purity",10}");
    Log($"  {new string('-', 42)}");

    prop.Propagate(rho0, 20.0, 0.05, tMeasure, (t, rho) =>
    {
        double mi = DensityMatrixTools.MutualInformation(rho, n,
            new[] { 0, 1 }, new[] { 3, 4 });
        var rho24 = DensityMatrixTools.PartialTrace(rho, n, new[] { 1, 3 });
        var (cpsi, _, _) = DensityMatrixTools.ComputeCPsi(rho24);
        double pur = DensityMatrixTools.Purity(rho);
        Log($"  {t,6:F1}  {mi,10:F6}  {cpsi,10:F6}  {pur,10:F6}");
    });

    Log();
    Log($"Expected from Python: MI(A:B) max ~ 0.86, CPsi_24 max ~ 0.24");
    Log($"[Test 0 completed in {sw.Elapsed.TotalSeconds:F1}s]");
    Log();
}

// ============================================================
// TEST 1: Cross-Bridge Information Flow (Level 3, N=11)
// ============================================================
void RunTest1()
{
    Log(new string('=', 70));
    Log("TEST 1: CROSS-BRIDGE INFORMATION FLOW (Level 3, N=11)");
    Log(new string('=', 70));
    Log();

    int n = 11;
    int d = 1 << n; // 2048
    double gamma = 0.05;
    var gammas = Enumerable.Repeat(gamma, n).ToArray();

    Log($"N={n}, d={d}, d^2={d * d}");
    Log("Topology: (0-1)-m1-(3-4) -- M(5) -- (6-7)-m2-(9-10)");
    Log();

    var sw = Stopwatch.StartNew();
    Log("Building Hamiltonian...");
    var bonds = Topology.MediatorBridge(3);
    Log($"Bonds: {bonds.Length}");
    foreach (var b in bonds)
        Log($"  {b.QubitA} -- {b.QubitB}  J={b.J}");
    var H = Topology.BuildHamiltonian(n, bonds);
    Log($"H built: {sw.Elapsed.TotalSeconds:F1}s");
    Log();

    var prop = new LindbladPropagator(H, gammas, n);

    // Initial state: Bell(0,1) + neutral everywhere else
    // Bell(0,1) ⊗ |0>^9
    var psi = BellPair(0, 1, n);

    var rho0 = DensityMatrixTools.PureState(psi);

    // Measure points
    var tMeasure = new double[41];
    for (int i = 0; i <= 40; i++) tMeasure[i] = i * 0.5;

    // Bridge A = {0,1,2,3,4}, Bridge B = {6,7,8,9,10}, Meta-M = {5}
    int[] bridgeA = { 0, 1, 2, 3, 4 };
    int[] bridgeB = { 6, 7, 8, 9, 10 };
    int[] pairA = { 0, 1 };
    int[] pairD = { 9, 10 };

    Log("--- Bell(0,1) initial state ---");
    Log();
    Log($"  {"t",6}  {"MI(BA:BB)",10}  {"MI(A:D)",10}  {"CPsi_AD",10}  {"Conc_AD",10}");
    Log($"  {new string('-', 52)}");

    double maxMI_AB = 0, maxMI_AD = 0;

    prop.Propagate(rho0, 20.0, 0.05, tMeasure, (t, rho) =>
    {
        double miBridges = DensityMatrixTools.MutualInformation(rho, n, bridgeA, bridgeB);
        double miAD = DensityMatrixTools.MutualInformation(rho, n, pairA, pairD);
        var rhoAD = DensityMatrixTools.PartialTrace(rho, n, new[] { 0, 1, 9, 10 });
        // For CPsi, extract 2-qubit state of qubits 0 and 10
        var rho_0_10 = DensityMatrixTools.PartialTrace(rho, n, new[] { 0, 10 });
        var (cpsi, _, _) = DensityMatrixTools.ComputeCPsi(rho_0_10);
        double conc = 0;
        try { conc = DensityMatrixTools.Concurrence2Q(rho_0_10); } catch { }

        if (miBridges > maxMI_AB) maxMI_AB = miBridges;
        if (miAD > maxMI_AD) maxMI_AD = miAD;

        // Print every 4th point to keep output manageable
        if (t < 0.01 || Math.Abs(t % 2.0) < 0.01 || Math.Abs(t - 20) < 0.01)
        {
            Log($"  {t,6:F1}  {miBridges,10:F6}  {miAD,10:F6}  {cpsi,10:F6}  {conc,10:F6}");
        }
    });

    Log();
    Log($"Peak MI(BridgeA : BridgeB): {maxMI_AB:F6}");
    Log($"Peak MI(PairA : PairD): {maxMI_AD:F6}");

    if (maxMI_AB > 0.001)
        Log("*** CROSS-BRIDGE INFORMATION FLOWS ***");
    else
        Log("No significant cross-bridge information flow detected.");

    Log();
    Log($"[Test 1 completed in {sw.Elapsed.TotalSeconds:F1}s]");
    Log();
}

// ============================================================
// TEST 2: Three Conditions at Level 3
// ============================================================
void RunTest2()
{
    Log(new string('=', 70));
    Log("TEST 2: THREE CONDITIONS AT LEVEL 3");
    Log(new string('=', 70));
    Log();

    int n = 11;
    int d = 1 << n;
    double gamma = 0.05;

    // (a) Meta-mediator coupling sweep
    Log("--- (a) J_meta sweep ---");
    Log();
    Log($"  {"J_meta",8}  {"MI(BA:BB)",10}  {"MI(A:D)",10}");
    Log($"  {new string('-', 32)}");

    // Initial state: Bell(0,1) ⊗ |0>^9
    var psi = BellPair(0, 1, n);
    var rho0 = DensityMatrixTools.PureState(psi);

    int[] bridgeA = { 0, 1, 2, 3, 4 };
    int[] bridgeB = { 6, 7, 8, 9, 10 };
    int[] pairA = { 0, 1 };
    int[] pairD = { 9, 10 };

    foreach (double jMeta in new[] { 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0 })
    {
        var bonds = Topology.MediatorBridge(3, jMeta: jMeta);
        var H = Topology.BuildHamiltonian(n, bonds);
        var gammas = Enumerable.Repeat(gamma, n).ToArray();
        var prop = new LindbladPropagator(H, gammas, n);

        double maxMI_AB = 0, maxMI_AD = 0;
        var tMeas = new[] { 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0 };

        prop.Propagate(rho0, 20.0, 0.05, tMeas, (t, rho) =>
        {
            double mi1 = DensityMatrixTools.MutualInformation(rho, n, bridgeA, bridgeB);
            double mi2 = DensityMatrixTools.MutualInformation(rho, n, pairA, pairD);
            if (mi1 > maxMI_AB) maxMI_AB = mi1;
            if (mi2 > maxMI_AD) maxMI_AD = mi2;
        });

        Log($"  {jMeta,8:F2}  {maxMI_AB,10:F6}  {maxMI_AD,10:F6}");
    }
    Log();

    // (b) Meta-mediator noise sweep
    Log("--- (b) gamma_M (qubit 5) sweep ---");
    Log();
    Log($"  {"gM",8}  {"MI(BA:BB)",10}  {"MI(A:D)",10}");
    Log($"  {new string('-', 32)}");

    foreach (double gM in new[] { 0.0, 0.01, 0.05, 0.1, 0.2, 0.5 })
    {
        var bonds = Topology.MediatorBridge(3);
        var H = Topology.BuildHamiltonian(n, bonds);
        var gammas = Enumerable.Repeat(gamma, n).ToArray();
        gammas[5] = gM;  // meta-mediator noise
        var prop = new LindbladPropagator(H, gammas, n);

        double maxMI_AB = 0, maxMI_AD = 0;
        var tMeas = new[] { 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0 };

        prop.Propagate(rho0, 20.0, 0.05, tMeas, (t, rho) =>
        {
            double mi1 = DensityMatrixTools.MutualInformation(rho, n, bridgeA, bridgeB);
            double mi2 = DensityMatrixTools.MutualInformation(rho, n, pairA, pairD);
            if (mi1 > maxMI_AB) maxMI_AB = mi1;
            if (mi2 > maxMI_AD) maxMI_AD = mi2;
        });

        Log($"  {gM,8:F3}  {maxMI_AB,10:F6}  {maxMI_AD,10:F6}");
    }
    Log();
}

// ============================================================
// TEST 3: Standing Wave (correlators across the bridge)
// ============================================================
void RunTest3()
{
    Log(new string('=', 70));
    Log("TEST 3: STANDING WAVE ACROSS THE BRIDGE");
    Log(new string('=', 70));
    Log();

    int n = 11;
    int d = 1 << n;
    double gamma = 0.05;
    var gammas = Enumerable.Repeat(gamma, n).ToArray();

    var bonds = Topology.MediatorBridge(3);
    var H = Topology.BuildHamiltonian(n, bonds);
    var prop = new LindbladPropagator(H, gammas, n);

    // Bell(0,1) initial
    var psi = BellPair(0, 1, n);
    var rho0 = DensityMatrixTools.PureState(psi);

    // Build correlator operators (only a few key ones)
    Log("Building correlator operators...");
    var XX_01 = PauliOps.At(PauliOps.X, 0, n) * PauliOps.At(PauliOps.X, 1, n);
    var XX_4_6 = PauliOps.At(PauliOps.X, 4, n) * PauliOps.At(PauliOps.X, 6, n);
    var XX_0_10 = PauliOps.At(PauliOps.X, 0, n) * PauliOps.At(PauliOps.X, 10, n);
    var ZZ_01 = PauliOps.At(PauliOps.Z, 0, n) * PauliOps.At(PauliOps.Z, 1, n);
    var ZZ_4_6 = PauliOps.At(PauliOps.Z, 4, n) * PauliOps.At(PauliOps.Z, 6, n);
    var ZZ_0_10 = PauliOps.At(PauliOps.Z, 0, n) * PauliOps.At(PauliOps.Z, 10, n);
    Log("Operators built.");
    Log();

    Log($"  {"t",6}  {"XX_01",10}  {"ZZ_01",10}  {"XX_4_6",10}  {"ZZ_4_6",10}  {"XX_0_10",10}  {"ZZ_0_10",10}");
    Log($"  {new string('-', 70)}");

    var tMeas = new double[11];
    for (int i = 0; i <= 10; i++) tMeas[i] = i * 2.0;

    prop.Propagate(rho0, 20.0, 0.05, tMeas, (t, rho) =>
    {
        var xx01 = DensityMatrixTools.ExpectationValue(rho, XX_01).Real;
        var zz01 = DensityMatrixTools.ExpectationValue(rho, ZZ_01).Real;
        var xx46 = DensityMatrixTools.ExpectationValue(rho, XX_4_6).Real;
        var zz46 = DensityMatrixTools.ExpectationValue(rho, ZZ_4_6).Real;
        var xx010 = DensityMatrixTools.ExpectationValue(rho, XX_0_10).Real;
        var zz010 = DensityMatrixTools.ExpectationValue(rho, ZZ_0_10).Real;
        Log($"  {t,6:F1}  {xx01,10:F6}  {zz01,10:F6}  {xx46,10:F6}  {zz46,10:F6}  {xx010,10:F6}  {zz010,10:F6}");
    });

    Log();
}

// ============================================================
// PULL PRINCIPLE TESTS
// ============================================================

/// <summary>
/// Part 1: Scaling curve - MI vs N for uniform chain and hierarchical mediator.
/// </summary>
void RunPull1()
{
    Log(new string('=', 70));
    Log("PULL PART 1: SCALING CURVE");
    Log(new string('=', 70));
    Log();

    double gamma = 0.05;

    Log($"  {"N",4}  {"Type",12}  {"MI_peak",10}  {"t_peak",8}  {"MI_ss",10}");
    Log($"  {new string('-', 50)}");

    foreach (int n in new[] { 3, 5, 7, 9, 11 })
    {
        int d = 1 << n;

        var rho0 = DensityMatrixTools.PureState(BellPair(0, 1, n));

        int[] pairFirst = { 0, 1 };
        int[] pairLast = { n - 2, n - 1 };

        // Series A: uniform chain
        {
            var couplings = Enumerable.Repeat(1.0, n - 1).ToArray();
            var bonds = Topology.Chain(n, couplings);
            var H = Topology.BuildHamiltonian(n, bonds);
            var gammas = Enumerable.Repeat(gamma, n).ToArray();
            var prop = new LindbladPropagator(H, gammas, n);

            double maxMI = 0, tPeak = 0, miSS = 0;
            var tMeas = new double[21];
            for (int i = 0; i <= 20; i++) tMeas[i] = i;

            prop.Propagate(rho0, 20.0, 0.05, tMeas, (t, rho) =>
            {
                double mi = DensityMatrixTools.MutualInformation(rho, n, pairFirst, pairLast);
                if (mi > maxMI) { maxMI = mi; tPeak = t; }
                if (Math.Abs(t - 20) < 0.01) miSS = mi;
            });

            Log($"  {n,4}  {"uniform",12}  {maxMI,10:F6}  {tPeak,8:F1}  {miSS,10:F6}");
        }

        // Series B: hierarchical mediator bridge
        if (n >= 5)
        {
            int level = n == 5 ? 2 : (n == 11 ? 3 : 0);
            Bond[] bonds;
            if (level > 0)
            {
                bonds = Topology.MediatorBridge(level);
            }
            else
            {
                // N=7, 9: extended chain, not perfect hierarchy
                var couplings = Enumerable.Repeat(1.0, n - 1).ToArray();
                bonds = Topology.Chain(n, couplings);
            }
            var H = Topology.BuildHamiltonian(n, bonds);
            var gammas = Enumerable.Repeat(gamma, n).ToArray();
            var prop = new LindbladPropagator(H, gammas, n);

            double maxMI = 0, tPeak = 0, miSS = 0;
            var tMeas = new double[21];
            for (int i = 0; i <= 20; i++) tMeas[i] = i;

            prop.Propagate(rho0, 20.0, 0.05, tMeas, (t, rho) =>
            {
                double mi = DensityMatrixTools.MutualInformation(rho, n, pairFirst, pairLast);
                if (mi > maxMI) { maxMI = mi; tPeak = t; }
                if (Math.Abs(t - 20) < 0.01) miSS = mi;
            });

            Log($"  {n,4}  {"hierarchical",12}  {maxMI,10:F6}  {tPeak,8:F1}  {miSS,10:F6}");
        }
    }
    Log();
}

/// <summary>
/// Part 2: Pull optimization - 2:1 coupling at various levels.
/// </summary>
void RunPull2()
{
    Log(new string('=', 70));
    Log("PULL PART 2: THE PULL OPTIMIZATION (2:1)");
    Log(new string('=', 70));
    Log();

    int n = 11;
    int d = 1 << n;
    double gamma = 0.05;
    var gammas = Enumerable.Repeat(gamma, n).ToArray();

    // Bell(0,1) initial
    var psi = BellPair(0, 1, n);
    var rho0 = DensityMatrixTools.PureState(psi);

    int[] bridgeA = { 0, 1, 2, 3, 4 };
    int[] bridgeB = { 6, 7, 8, 9, 10 };
    int[] pairA = { 0, 1 };
    int[] pairD = { 9, 10 };

    // Five configurations
    var configs = new (string Name, double[] Couplings)[]
    {
        // Bonds: 0-1, 1-2, 2-3, 3-4, 4-5, 5-6, 6-7, 7-8, 8-9, 9-10
        ("symmetric",    new[] { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 }),
        ("2:1 bridges",  new[] { 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0 }),
        ("2:1 meta",     new[] { 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0 }),
        ("2:1 all pull", new[] { 1.0, 1.0, 2.0, 1.0, 1.0, 2.0, 1.0, 1.0, 2.0, 1.0 }),
        ("2:1 reversed", new[] { 1.0, 2.0, 1.0, 1.0, 2.0, 1.0, 1.0, 2.0, 1.0, 1.0 }),
    };

    Log($"  {"Config",15}  {"MI(BA:BB)",10}  {"MI(A:D)",10}  {"t_peak",8}");
    Log($"  {new string('-', 48)}");

    foreach (var (name, couplings) in configs)
    {
        // Build chain with specified couplings
        var bonds = new Bond[10];
        for (int i = 0; i < 10; i++)
            bonds[i] = new Bond(i, i + 1, couplings[i], new[] { "X", "Y", "Z" });

        var H = Topology.BuildHamiltonian(n, bonds);
        var prop = new LindbladPropagator(H, gammas, n);

        double maxMI_AB = 0, maxMI_AD = 0, tPeak = 0;
        var tMeas = new double[21];
        for (int i = 0; i <= 20; i++) tMeas[i] = i;

        prop.Propagate(rho0, 20.0, 0.05, tMeas, (t, rho) =>
        {
            double mi1 = DensityMatrixTools.MutualInformation(rho, n, bridgeA, bridgeB);
            double mi2 = DensityMatrixTools.MutualInformation(rho, n, pairA, pairD);
            if (mi1 > maxMI_AB) { maxMI_AB = mi1; tPeak = t; }
            if (mi2 > maxMI_AD) maxMI_AD = mi2;
        });

        Log($"  {name,15}  {maxMI_AB,10:F6}  {maxMI_AD,10:F6}  {tPeak,8:F1}");
    }
    Log();
}

/// <summary>
/// Part 3: Relay protocol - time-dependent gamma.
/// </summary>
void RunPull3()
{
    Log(new string('=', 70));
    Log("PULL PART 3: RELAY PROTOCOL");
    Log(new string('=', 70));
    Log();

    int n = 11;
    int d = 1 << n;
    double gamma = 0.05;
    double gammaQuiet = 0.005;
    double stageTime = 0.78;    // 0.039/gamma = 0.039/0.05

    // Bell(0,1) initial
    var psi = BellPair(0, 1, n);
    var rho0 = DensityMatrixTools.PureState(psi);

    int[] pairA = { 0, 1 };
    int[] pairD = { 9, 10 };
    int[] bridgeA = { 0, 1, 2, 3, 4 };
    int[] bridgeB = { 6, 7, 8, 9, 10 };

    var bonds = Topology.MediatorBridge(3);
    var H = Topology.BuildHamiltonian(n, bonds);

    // (a) Passive baseline
    {
        var gammas = Enumerable.Repeat(gamma, n).ToArray();
        var prop = new LindbladPropagator(H, gammas, n);

        double maxMI_AD = 0, maxMI_AB = 0;
        var tMeas = new double[21];
        for (int i = 0; i <= 20; i++) tMeas[i] = i;

        prop.Propagate(rho0, 20.0, 0.05, tMeas, (t, rho) =>
        {
            double mi1 = DensityMatrixTools.MutualInformation(rho, n, bridgeA, bridgeB);
            double mi2 = DensityMatrixTools.MutualInformation(rho, n, pairA, pairD);
            if (mi1 > maxMI_AB) maxMI_AB = mi1;
            if (mi2 > maxMI_AD) maxMI_AD = mi2;
        });

        Log($"  Passive:  MI(BA:BB)={maxMI_AB:F6}  MI(A:D)={maxMI_AD:F6}");
    }

    // (b) Relay protocol (staged gamma)
    // 6 stages, each ~0.8 time units
    // Stage k: qubit receiving[k] has low gamma, rest normal
    int[][] relayReceivers = {
        new[]{ 2 },       // m1 receives from Pair A
        new[]{ 3, 4 },    // Pair B receives from m1
        new[]{ 5 },       // Meta-M receives from Pair B
        new[]{ 6, 7 },    // Pair C receives from Meta-M
        new[]{ 8 },       // m2 receives from Pair C
        new[]{ 9, 10 },   // Pair D receives from m2
    };

    {
        // Run staged: for each stage, build propagator with correct gamma
        var rho = rho0.Clone();
        double t = 0;

        for (int stage = 0; stage < relayReceivers.Length; stage++)
        {
            var gammas = Enumerable.Repeat(gamma, n).ToArray();
            foreach (int q in relayReceivers[stage])
                gammas[q] = gammaQuiet;

            var prop = new LindbladPropagator(H, gammas, n);

            // Propagate this stage
            double tEnd = (stage + 1) * stageTime;
            int nSteps = (int)(stageTime / 0.05);
            for (int step = 0; step < nSteps; step++)
            {
                var k1 = prop.EvalRHS(rho);
                var k2 = prop.EvalRHS(rho + 0.05 / 2 * k1);
                var k3 = prop.EvalRHS(rho + 0.05 / 2 * k2);
                var k4 = prop.EvalRHS(rho + 0.05 * k3);
                rho = rho + (0.05 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
                rho = (rho + rho.ConjugateTranspose()) / 2.0;
            }
            t = tEnd;
        }

        double miAB = DensityMatrixTools.MutualInformation(rho, n, bridgeA, bridgeB);
        double miAD = DensityMatrixTools.MutualInformation(rho, n, pairA, pairD);
        Log($"  Relay:    MI(BA:BB)={miAB:F6}  MI(A:D)={miAD:F6}  (t={t:F1})");
    }

    // (c) Relay + 2:1 pull
    {
        var bonds21 = new Bond[10];
        var couplings = new[] { 1.0, 1.0, 2.0, 1.0, 1.0, 2.0, 1.0, 1.0, 2.0, 1.0 };
        for (int i = 0; i < 10; i++)
            bonds21[i] = new Bond(i, i + 1, couplings[i], new[] { "X", "Y", "Z" });
        var H21 = Topology.BuildHamiltonian(n, bonds21);

        var rho = rho0.Clone();
        double t = 0;

        for (int stage = 0; stage < relayReceivers.Length; stage++)
        {
            var gammas = Enumerable.Repeat(gamma, n).ToArray();
            foreach (int q in relayReceivers[stage])
                gammas[q] = gammaQuiet;

            var prop = new LindbladPropagator(H21, gammas, n);

            int nSteps = (int)(stageTime / 0.05);
            for (int step = 0; step < nSteps; step++)
            {
                var k1 = prop.EvalRHS(rho);
                var k2 = prop.EvalRHS(rho + 0.05 / 2 * k1);
                var k3 = prop.EvalRHS(rho + 0.05 / 2 * k2);
                var k4 = prop.EvalRHS(rho + 0.05 * k3);
                rho = rho + (0.05 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
                rho = (rho + rho.ConjugateTranspose()) / 2.0;
            }
            t = (stage + 1) * stageTime;
        }

        double miAB = DensityMatrixTools.MutualInformation(rho, n, bridgeA, bridgeB);
        double miAD = DensityMatrixTools.MutualInformation(rho, n, pairA, pairD);
        Log($"  Relay+2:1 MI(BA:BB)={miAB:F6}  MI(A:D)={miAD:F6}  (t={t:F1})");
    }

    Log();
}

// ============================================================
// PROFILE EVALUATION MODE
// ============================================================
void RunProfileEvaluation(string[] pArgs)
{
    // Usage: profile <N> <g1,g2,...,gN> [--tmax 20] [--dt 0.05]
    if (pArgs.Length < 3)
    {
        Console.Error.WriteLine("Usage: profile <N> <g1,g2,...,gN> [--tmax 20] [--dt 0.05]");
        Environment.Exit(1);
        return;
    }

    int n = int.Parse(pArgs[1]);
    var gammas = pArgs[2].Split(',')
        .Select(s => double.Parse(s, CultureInfo.InvariantCulture))
        .ToArray();

    if (gammas.Length != n)
    {
        Console.Error.WriteLine($"ERROR: Expected {n} gamma values, got {gammas.Length}");
        Environment.Exit(1);
        return;
    }

    double tMax = 20.0;
    double dt = 0.05;
    for (int i = 3; i < pArgs.Length - 1; i++)
    {
        if (pArgs[i] == "--tmax")
            tMax = double.Parse(pArgs[i + 1], CultureInfo.InvariantCulture);
        if (pArgs[i] == "--dt")
            dt = double.Parse(pArgs[i + 1], CultureInfo.InvariantCulture);
    }

    int d = 1 << n;

    // Heisenberg chain Hamiltonian (J=1.0, N-1 bonds)
    var couplings = Enumerable.Repeat(1.0, n - 1).ToArray();
    var bonds = Topology.Chain(n, couplings);
    var H = Topology.BuildHamiltonian(n, bonds);

    // Initial state: |+>^N  (matches resonant_return.py Test 1 setup)
    var psi = new Complex[d];
    double psiNorm = 1.0 / Math.Sqrt(d);
    for (int idx = 0; idx < d; idx++) psi[idx] = psiNorm;
    var rho0 = DensityMatrixTools.PureState(psi);

    // Propagator with user-specified gamma profile
    var prop = new LindbladPropagator(H, gammas, n);

    // Measure every 0.5 time units
    double mInterval = 0.5;
    int nMeas = (int)(tMax / mInterval) + 1;
    var tMeas = new double[nMeas];
    for (int i = 0; i < nMeas; i++) tMeas[i] = i * mInterval;

    double bestSumMI = -1, bestT = 0, bestPeakMI = -1;
    double cpsi01Best = 0, purBest = 0;
    double sumMI5 = 0;

    prop.Propagate(rho0, tMax, dt, tMeas, (t, rho) =>
    {
        // Sum-MI: sum of MI(k, k+1) for all adjacent pairs
        double sumMI = 0;
        for (int k = 0; k < n - 1; k++)
            sumMI += DensityMatrixTools.MutualInformation(rho, n,
                new[] { k }, new[] { k + 1 });

        // MI(0, N-1)
        double mi0N = DensityMatrixTools.MutualInformation(rho, n,
            new[] { 0 }, new[] { n - 1 });

        if (sumMI > bestSumMI)
        {
            bestSumMI = sumMI;
            bestT = t;
            var rho01 = DensityMatrixTools.PartialTrace(rho, n, new[] { 0, 1 });
            var (cpsiVal, _, _) = DensityMatrixTools.ComputeCPsi(rho01);
            cpsi01Best = cpsiVal;
            purBest = DensityMatrixTools.Purity(rho);
        }

        if (mi0N > bestPeakMI) bestPeakMI = mi0N;
        if (Math.Abs(t - 5.0) < 0.01) sumMI5 = sumMI;
    });

    // Output: single machine-parseable line
    Console.WriteLine(FormattableString.Invariant(
        $"RESULT SumMI={bestSumMI:F6} PeakMI={bestPeakMI:F6} PeakT={bestT:F2} CPsi01={cpsi01Best:F6} Purity={purBest:F6} SumMI5={sumMI5:F6}"));
}

// ============================================================
// STAGED EVALUATION MODE
// ============================================================
void RunStagedEvaluation(string[] pArgs)
{
    // Usage: staged <N> <t_switch> <g1,...,gN> <g1,...,gN> [--tmax 20] [--dt 0.05]
    if (pArgs.Length < 5)
    {
        Console.Error.WriteLine("Usage: staged <N> <t_switch> <g1,...,gN> <g1,...,gN> [--tmax 20] [--dt 0.05]");
        Environment.Exit(1);
        return;
    }

    int n = int.Parse(pArgs[1]);
    double tSwitch = double.Parse(pArgs[2], CultureInfo.InvariantCulture);
    var gammas1 = pArgs[3].Split(',')
        .Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray();
    var gammas2 = pArgs[4].Split(',')
        .Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray();

    if (gammas1.Length != n || gammas2.Length != n)
    {
        Console.Error.WriteLine($"ERROR: Expected {n} gamma values per phase, got {gammas1.Length} and {gammas2.Length}");
        Environment.Exit(1);
        return;
    }

    double tMax = 20.0;
    double dt = 0.05;
    for (int i = 5; i < pArgs.Length - 1; i++)
    {
        if (pArgs[i] == "--tmax")
            tMax = double.Parse(pArgs[i + 1], CultureInfo.InvariantCulture);
        if (pArgs[i] == "--dt")
            dt = double.Parse(pArgs[i + 1], CultureInfo.InvariantCulture);
    }

    // Heisenberg chain (same as profile mode)
    var couplings = Enumerable.Repeat(1.0, n - 1).ToArray();
    var bonds = Topology.Chain(n, couplings);
    var H = Topology.BuildHamiltonian(n, bonds);

    // Initial state |+>^N
    int d = 1 << n;
    var psi = new Complex[d];
    double psiNorm = 1.0 / Math.Sqrt(d);
    for (int idx = 0; idx < d; idx++) psi[idx] = psiNorm;
    var rho = DensityMatrixTools.PureState(psi);

    // Tracking
    double bestSumMI = -1, bestPeakMI = -1, bestT = 0;
    double mInterval = 0.5;
    double nextMeas = mInterval;
    int totalSteps = (int)(tMax / dt);
    int switchStep = (int)(tSwitch / dt);

    // Build initial propagator
    var prop = new LindbladPropagator(H, tSwitch > 0 ? gammas1 : gammas2, n);

    for (int step = 0; step < totalSteps; step++)
    {
        // Switch propagator at the boundary
        if (step == switchStep && tSwitch > 0)
            prop = new LindbladPropagator(H, gammas2, n);

        // RK4 step (same pattern as relay protocol)
        var k1 = prop.EvalRHS(rho);
        var k2 = prop.EvalRHS(rho + dt / 2 * k1);
        var k3 = prop.EvalRHS(rho + dt / 2 * k2);
        var k4 = prop.EvalRHS(rho + dt * k3);
        rho = rho + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
        rho = (rho + rho.ConjugateTranspose()) / 2.0;

        double t = (step + 1) * dt;

        // Measure at mInterval intervals
        if (t >= nextMeas - dt / 2)
        {
            int phase = step < switchStep ? 1 : 2;

            double sumMI = 0;
            for (int k = 0; k < n - 1; k++)
                sumMI += DensityMatrixTools.MutualInformation(rho, n,
                    new[] { k }, new[] { k + 1 });
            double mi0N = DensityMatrixTools.MutualInformation(rho, n,
                new[] { 0 }, new[] { n - 1 });

            Console.WriteLine(FormattableString.Invariant(
                $"T={nextMeas:F2} SumMI={sumMI:F6} PeakMI={mi0N:F6} Phase={phase}"));

            if (sumMI > bestSumMI) { bestSumMI = sumMI; bestT = nextMeas; }
            if (mi0N > bestPeakMI) bestPeakMI = mi0N;

            nextMeas += mInterval;
        }
    }

    // Final RESULT line (machine-parseable)
    Console.WriteLine(FormattableString.Invariant(
        $"RESULT SumMI={bestSumMI:F6} PeakMI={bestPeakMI:F6} PeakT={bestT:F2} SwitchT={tSwitch:F2}"));
}

// ============================================================
// SWEEP EVALUATION MODE: multi-phase gamma profile sweep
// ============================================================
void RunSweepEvaluation(string[] pArgs)
{
    // Usage: sweep <N> <stage_time> <profile1> <profile2> ... [--tmax 20] [--dt 0.05]
    if (pArgs.Length < 4)
    {
        Console.Error.WriteLine("Usage: sweep <N> <stage_time> <profile1> <profile2> ... [--tmax 20] [--dt 0.05]");
        Environment.Exit(1);
        return;
    }

    int n = int.Parse(pArgs[1]);
    double stageTime = double.Parse(pArgs[2], CultureInfo.InvariantCulture);

    // Collect profiles until we hit a flag
    var profiles = new List<double[]>();
    double tMax = 20.0;
    double dt = 0.05;

    for (int i = 3; i < pArgs.Length; i++)
    {
        if (pArgs[i] == "--tmax") { tMax = double.Parse(pArgs[++i], CultureInfo.InvariantCulture); continue; }
        if (pArgs[i] == "--dt") { dt = double.Parse(pArgs[++i], CultureInfo.InvariantCulture); continue; }
        if (pArgs[i] == "--diag") continue;

        var g = pArgs[i].Split(',')
            .Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray();
        if (g.Length != n)
        {
            Console.Error.WriteLine($"ERROR: Profile {profiles.Count + 1}: expected {n} gammas, got {g.Length}");
            Environment.Exit(1);
            return;
        }
        profiles.Add(g);
    }

    if (profiles.Count == 0) { Console.Error.WriteLine("ERROR: No profiles provided"); Environment.Exit(1); return; }

    // Heisenberg chain (same as profile mode)
    var couplings = Enumerable.Repeat(1.0, n - 1).ToArray();
    var bonds = Topology.Chain(n, couplings);
    var H = Topology.BuildHamiltonian(n, bonds);

    // Initial state |+>^N
    int d = 1 << n;
    var psi = new Complex[d];
    double psiNorm = 1.0 / Math.Sqrt(d);
    for (int idx = 0; idx < d; idx++) psi[idx] = psiNorm;
    var rho = DensityMatrixTools.PureState(psi);

    // Diagnostic mode: output CΨ for all pairs
    bool diagMode = pArgs.Any(a => a == "--diag");

    // Tracking
    double bestSumMI = -1, bestPeakMI = -1, bestT = 0;
    double mInterval = 0.5;
    double nextMeas = mInterval;
    int totalSteps = (int)(tMax / dt);
    int stepsPerStage = (int)(stageTime / dt);

    // Build first propagator
    int currentProfile = 0;
    var prop = new LindbladPropagator(H, profiles[0], n);

    for (int step = 0; step < totalSteps; step++)
    {
        // Switch to next profile when stage boundary is crossed
        int profileIdx = Math.Min(step / stepsPerStage, profiles.Count - 1);
        if (profileIdx != currentProfile)
        {
            currentProfile = profileIdx;
            prop = new LindbladPropagator(H, profiles[currentProfile], n);
        }

        // RK4 step
        var k1 = prop.EvalRHS(rho);
        var k2 = prop.EvalRHS(rho + dt / 2 * k1);
        var k3 = prop.EvalRHS(rho + dt / 2 * k2);
        var k4 = prop.EvalRHS(rho + dt * k3);
        rho = rho + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
        rho = (rho + rho.ConjugateTranspose()) / 2.0;

        double t = (step + 1) * dt;

        // Measure at mInterval intervals
        if (t >= nextMeas - dt / 2)
        {
            double sumMI = 0;
            var miPairs = new double[n - 1];
            for (int k = 0; k < n - 1; k++)
            {
                miPairs[k] = DensityMatrixTools.MutualInformation(rho, n,
                    new[] { k }, new[] { k + 1 });
                sumMI += miPairs[k];
            }
            double mi0N = DensityMatrixTools.MutualInformation(rho, n,
                new[] { 0 }, new[] { n - 1 });

            if (diagMode)
            {
                // CΨ for all adjacent pairs and endpoints
                var cpsiValues = new double[n - 1];
                for (int k = 0; k < n - 1; k++)
                {
                    var rhoK = DensityMatrixTools.PartialTrace(rho, n, new[] { k, k + 1 });
                    var (cpsi, _, _) = DensityMatrixTools.ComputeCPsi(rhoK);
                    cpsiValues[k] = cpsi;
                }
                var rho0N = DensityMatrixTools.PartialTrace(rho, n, new[] { 0, n - 1 });
                var (cpsi0N, _, _) = DensityMatrixTools.ComputeCPsi(rho0N);
                double pur = DensityMatrixTools.Purity(rho);

                var cpsiStr = string.Join(" ",
                    cpsiValues.Select((c, i) => FormattableString.Invariant($"CPsi{i}{i+1}={c:F4}")));
                Console.WriteLine(FormattableString.Invariant(
                    $"T={nextMeas:F2} SumMI={sumMI:F6} PeakMI={mi0N:F6} Stage={currentProfile + 1} {cpsiStr} CPsi0{n-1}={cpsi0N:F4} Pur={pur:F6}"));
            }
            else
            {
                Console.WriteLine(FormattableString.Invariant(
                    $"T={nextMeas:F2} SumMI={sumMI:F6} PeakMI={mi0N:F6} Stage={currentProfile + 1}"));
            }

            if (sumMI > bestSumMI) { bestSumMI = sumMI; bestT = nextMeas; }
            if (mi0N > bestPeakMI) bestPeakMI = mi0N;

            nextMeas += mInterval;
        }
    }

    // Final RESULT line
    Console.WriteLine(FormattableString.Invariant(
        $"RESULT SumMI={bestSumMI:F6} PeakMI={bestPeakMI:F6} PeakT={bestT:F2} Stages={profiles.Count} StageTime={stageTime:F2}"));
}

// ============================================================
// WAVE ANALYSIS: per-pair MI to see energy propagation
// ============================================================
void RunWaveAnalysis(string[] pArgs)
{
    // Usage: wave <N> <g1,...,gN> [--tmax 20]
    if (pArgs.Length < 3)
    {
        Console.Error.WriteLine("Usage: wave <N> <g1,...,gN> [--tmax 20]");
        Environment.Exit(1);
        return;
    }

    int n = int.Parse(pArgs[1]);
    var gammas = pArgs[2].Split(',')
        .Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray();

    double tMax = 20.0;
    double dt = 0.05;
    for (int i = 3; i < pArgs.Length - 1; i++)
    {
        if (pArgs[i] == "--tmax") tMax = double.Parse(pArgs[i + 1], CultureInfo.InvariantCulture);
    }

    int d = 1 << n;

    // Heisenberg chain
    var couplings = Enumerable.Repeat(1.0, n - 1).ToArray();
    var bonds = Topology.Chain(n, couplings);
    var H = Topology.BuildHamiltonian(n, bonds);

    // Initial state |+>^N
    var psi = new Complex[d];
    double psiNorm = 1.0 / Math.Sqrt(d);
    for (int idx = 0; idx < d; idx++) psi[idx] = psiNorm;
    var rho = DensityMatrixTools.PureState(psi);

    var prop = new LindbladPropagator(H, gammas, n);

    // Header
    Console.Write(FormattableString.Invariant($"{"T",5}"));
    for (int k = 0; k < n - 1; k++)
        Console.Write(FormattableString.Invariant($" {"MI" + k + "" + (k+1),8}"));
    Console.Write(FormattableString.Invariant($" {"SumMI",8} {"Purity",8}"));
    Console.WriteLine();

    double mInterval = 0.5;
    double nextMeas = 0.0;
    int totalSteps = (int)(tMax / dt);

    for (int step = 0; step <= totalSteps; step++)
    {
        if (step > 0)
        {
            var k1 = prop.EvalRHS(rho);
            var k2 = prop.EvalRHS(rho + dt / 2 * k1);
            var k3 = prop.EvalRHS(rho + dt / 2 * k2);
            var k4 = prop.EvalRHS(rho + dt * k3);
            rho = rho + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
            rho = (rho + rho.ConjugateTranspose()) / 2.0;
        }

        double t = step * dt;

        if (t >= nextMeas - dt / 2)
        {
            Console.Write(FormattableString.Invariant($"{t,5:F1}"));

            double sumMI = 0;
            for (int k = 0; k < n - 1; k++)
            {
                double mi = DensityMatrixTools.MutualInformation(rho, n,
                    new[] { k }, new[] { k + 1 });
                sumMI += mi;
                Console.Write(FormattableString.Invariant($" {mi,8:F4}"));
            }

            double pur = DensityMatrixTools.Purity(rho);
            Console.Write(FormattableString.Invariant($" {sumMI,8:F4} {pur,8:F4}"));
            Console.WriteLine();

            nextMeas += mInterval;
        }
    }
}

// ============================================================
// RESONANCE TEST: CΨ oscillation around ¼ with structured bath
// ============================================================
void RunResonanceTest(string[] pArgs)
{
    // Usage: resonance <N> <g1,...,gN> [--tmax 20] [--j J_coupling]
    int n = pArgs.Length >= 2 ? int.Parse(pArgs[1]) : 3;
    double[] gammas;
    if (pArgs.Length >= 3 && pArgs[2].Contains(','))
        gammas = pArgs[2].Split(',')
            .Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray();
    else
        gammas = new[] { 0.05, 0.05, 0.001 }; // default: 2 system + 1 quiet bath

    double tMax = 20.0;
    double jCoupling = 1.0;
    for (int i = 3; i < pArgs.Length - 1; i++)
    {
        if (pArgs[i] == "--tmax") tMax = double.Parse(pArgs[i + 1], CultureInfo.InvariantCulture);
        if (pArgs[i] == "--j") jCoupling = double.Parse(pArgs[i + 1], CultureInfo.InvariantCulture);
    }

    int d = 1 << n;
    double dt = 0.01; // fine resolution for oscillation detection
    double mInterval = 0.05; // measure every 0.05

    Console.WriteLine($"=== RESONANCE TEST: CΨ OSCILLATION ===");
    Console.WriteLine(FormattableString.Invariant(
        $"N={n}, J={jCoupling:F2}, γ=[{string.Join(",", gammas.Select(g => FormattableString.Invariant($"{g:F4}")))}]"));
    Console.WriteLine(FormattableString.Invariant($"Σγ={gammas.Sum():F4}, dt={dt}, tmax={tMax}"));
    Console.WriteLine();

    // Build Hamiltonian
    var couplings = Enumerable.Repeat(jCoupling, n - 1).ToArray();
    var bonds = Topology.Chain(n, couplings);
    var H = Topology.BuildHamiltonian(n, bonds);

    // Initial state: Bell(0,1)⊗|+>_rest if --bell, else |+>^N
    bool bellInit = pArgs.Any(a => a == "--bell");
    Complex[] psi;
    if (bellInit && n >= 2)
    {
        // Bell(0,1): (|00...0> + |11...0>) / √2, with remaining qubits in |+>
        // First build Bell pair on qubits 0,1
        var bellPsi = new Complex[d];
        bellPsi[0] = 1.0 / Math.Sqrt(2);
        bellPsi[(1 << (n - 1)) | (1 << (n - 2))] = 1.0 / Math.Sqrt(2);
        // Tensor with |+> on bath qubits: apply H gate to each bath qubit
        // For simplicity: build Bell(0,1) ⊗ |+>^(N-2) directly
        int dSys = 4; // 2 system qubits
        int dBath = 1 << (n - 2);
        psi = new Complex[d];
        double bathNorm = 1.0 / Math.Sqrt(dBath);
        for (int bIdx = 0; bIdx < dBath; bIdx++)
        {
            // |00>_sys ⊗ |bIdx>_bath
            psi[(0 << (n - 2)) | bIdx] += (1.0 / Math.Sqrt(2)) * bathNorm;
            // |11>_sys ⊗ |bIdx>_bath
            psi[(3 << (n - 2)) | bIdx] += (1.0 / Math.Sqrt(2)) * bathNorm;
        }
        Console.WriteLine("Initial state: Bell(0,1) ⊗ |+>^(N-2)");
    }
    else
    {
        psi = new Complex[d];
        double psiNorm = 1.0 / Math.Sqrt(d);
        for (int idx = 0; idx < d; idx++) psi[idx] = psiNorm;
        Console.WriteLine("Initial state: |+>^N");
    }
    var rho = DensityMatrixTools.PureState(psi);

    var prop = new LindbladPropagator(H, gammas, n);

    // Track crossings
    double prevCpsi = 1.0;
    int crossDown = 0, crossUp = 0;
    double nextMeas = mInterval;
    int totalSteps = (int)(tMax / dt);

    Console.WriteLine(FormattableString.Invariant(
        $"{"T",6} | {"CΨ(01)",8} | {"CΨ(0B)",8} | {"MI(01)",8} | {"MI(0B)",8} | {"Pur",7} | event"));

    for (int step = 0; step < totalSteps; step++)
    {
        // RK4 step
        var k1 = prop.EvalRHS(rho);
        var k2 = prop.EvalRHS(rho + dt / 2 * k1);
        var k3 = prop.EvalRHS(rho + dt / 2 * k2);
        var k4 = prop.EvalRHS(rho + dt * k3);
        rho = rho + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
        rho = (rho + rho.ConjugateTranspose()) / 2.0;

        double t = (step + 1) * dt;

        if (t >= nextMeas - dt / 2)
        {
            // System CΨ (qubits 0-1, tracing out bath)
            var rho01 = DensityMatrixTools.PartialTrace(rho, n, new[] { 0, 1 });
            var (cpsi01, _, _) = DensityMatrixTools.ComputeCPsi(rho01);

            // System-bath CΨ
            var rho0B = DensityMatrixTools.PartialTrace(rho, n, new[] { 0, n - 1 });
            var (cpsi0B, _, _) = DensityMatrixTools.ComputeCPsi(rho0B);

            // MI at crossings
            double mi01 = DensityMatrixTools.MutualInformation(rho, n, new[] { 0 }, new[] { 1 });
            double mi0B = DensityMatrixTools.MutualInformation(rho, n, new[] { 0 }, new[] { n - 1 });
            double pur = DensityMatrixTools.Purity(rho);

            // Detect crossings (both directions)
            string evt = "";
            if (prevCpsi >= 0.25 && cpsi01 < 0.25) { crossDown++; evt = $" ↓ CROSS #{crossDown + crossUp}"; }
            if (prevCpsi < 0.25 && cpsi01 >= 0.25) { crossUp++; evt = $" ↑ CROSS #{crossDown + crossUp}"; }

            // Print at crossings, at 0.5 intervals, and at start
            if (evt != "" || Math.Abs(t % 0.5) < mInterval || t < 0.2)
            {
                Console.WriteLine(FormattableString.Invariant(
                    $"{t,6:F2} | {cpsi01,8:F4} | {cpsi0B,8:F4} | {mi01,8:F6} | {mi0B,8:F6} | {pur,6:F4} | {evt}"));
            }

            prevCpsi = cpsi01;
            nextMeas += mInterval;
        }
    }

    Console.WriteLine();
    Console.WriteLine($"Total crossings: {crossDown} down + {crossUp} up = {crossDown + crossUp}");
}

// ============================================================
// SPECTRAL ANALYSIS: eigendecomposition + midpoint hypothesis test
// ============================================================
void RunSpectralAnalysis(string[] pArgs)
{
    // Usage: spectral <N> [g1,g2,...,gN]
    int n = pArgs.Length >= 2 ? int.Parse(pArgs[1]) : 3;
    double[] gammas;
    if (pArgs.Length >= 3 && pArgs[2].Contains(','))
        gammas = pArgs[2].Split(',')
            .Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray();
    else
        gammas = Enumerable.Repeat(0.05, n).ToArray();

    int d = 1 << n;
    int d2 = d * d;
    double sumGamma = gammas.Sum();

    Console.WriteLine($"=== SPECTRAL MIDPOINT ANALYSIS ===");
    Console.WriteLine($"N={n}, d={d}, d²={d2}");
    Console.WriteLine($"Gammas: [{string.Join(", ", gammas.Select(g => FormattableString.Invariant($"{g:F4}")))}]");
    Console.WriteLine($"Σγ = {sumGamma:F4} (palindromic midpoint rate)");
    Console.WriteLine();

    // Build Hamiltonian (Heisenberg chain)
    var couplings = Enumerable.Repeat(1.0, n - 1).ToArray();
    var bonds = Topology.Chain(n, couplings);
    var H = Topology.BuildHamiltonian(n, bonds);
    var prop = new LindbladPropagator(H, gammas, n);

    // Build Liouvillian matrix by probing with basis density matrices
    // Vectorization: vec(ρ)[j*d + i] = ρ[i, j] (column-major stacking)
    Console.WriteLine("Building Liouvillian matrix...");
    var L = Matrix<Complex>.Build.Dense(d2, d2);
    for (int col = 0; col < d2; col++)
    {
        var basis = Matrix<Complex>.Build.Dense(d, d);
        int bRow = col % d, bCol = col / d;
        basis[bRow, bCol] = Complex.One;
        var Lbasis = prop.EvalRHS(basis);
        for (int row = 0; row < d2; row++)
            L[row, col] = Lbasis[row % d, row / d];
    }

    // Eigendecompose
    Console.WriteLine("Eigendecomposing...");
    var evd = L.Evd();
    var eigenvalues = evd.EigenValues;
    var V = evd.EigenVectors;

    // Classify eigenvalues by decay rate bands
    double gBand = gammas[0]; // band width = one γ_base
    int nMid = 0, nSlow = 0, nFast = 0, nImmune = 0;
    Console.WriteLine($"\nEigenvalue spectrum ({d2} modes):");
    Console.WriteLine($"Band definition: SLOW (rate < {sumGamma - gBand:F3}), MID ({sumGamma - gBand:F3}-{sumGamma + gBand:F3}), FAST (> {sumGamma + gBand:F3}), IMMUNE (rate < 0.001)");

    var modeInfo = new (double rate, double freq, string band, int idx)[d2];
    for (int i = 0; i < d2; i++)
    {
        double rate = -eigenvalues[i].Real;
        double freq = eigenvalues[i].Imaginary;
        string band;
        if (rate < 0.001) { band = "IMMUNE"; nImmune++; }
        else if (Math.Abs(rate - sumGamma) < gBand) { band = "MID"; nMid++; }
        else if (rate < sumGamma) { band = "SLOW"; nSlow++; }
        else { band = "FAST"; nFast++; }
        modeInfo[i] = (rate, freq, band, i);
    }
    Console.WriteLine($"  IMMUNE: {nImmune}  SLOW: {nSlow}  MID: {nMid}  FAST: {nFast}");
    Console.WriteLine();

    // Decompose initial state |+>^N
    var psi = new Complex[d];
    double psiNorm = 1.0 / Math.Sqrt(d);
    for (int idx = 0; idx < d; idx++) psi[idx] = psiNorm;
    var rho0 = DensityMatrixTools.PureState(psi);

    // Vectorize rho0 (column-major)
    var rho0vec = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(d2);
    for (int j = 0; j < d; j++)
        for (int i = 0; i < d; i++)
            rho0vec[j * d + i] = rho0[i, j];

    // Expansion coefficients: c = V^{-1} * rho0vec
    Console.WriteLine("Computing eigenbasis decomposition of |+>^N...");
    var coeffs = V.Solve(rho0vec);

    // Initial weight distribution
    double totW0 = 0, midW0 = 0, slowW0 = 0, fastW0 = 0, immW0 = 0;
    for (int i = 0; i < d2; i++)
    {
        double w = coeffs[i].Magnitude;
        totW0 += w;
        if (modeInfo[i].band == "MID") midW0 += w;
        else if (modeInfo[i].band == "SLOW") slowW0 += w;
        else if (modeInfo[i].band == "FAST") fastW0 += w;
        else immW0 += w;
    }
    Console.WriteLine($"Initial weights: IMMUNE={immW0/totW0:P1} SLOW={slowW0/totW0:P1} MID={midW0/totW0:P1} FAST={fastW0/totW0:P1}");
    Console.WriteLine();

    // Time evolution: track spectral bands and CΨ
    Console.WriteLine("Time evolution (spectral weights at each CΨ measurement):");
    Console.WriteLine(FormattableString.Invariant(
        $"{"T",6} | {"CΨ(01)",8} | {"CΨ(0N)",8} | {"D=1-4CΨ",8} | {"IMMUNE",7} | {"SLOW",7} | {"MID",7} | {"FAST",7} | note"));

    double prevCpsi01 = 1.0;
    bool crossingFound = false;

    for (double t = 0.1; t <= 10.0; t += 0.1)
    {
        // Spectral weights at time t
        double totW = 0, midW = 0, slowW = 0, fastW = 0, immW = 0;
        for (int i = 0; i < d2; i++)
        {
            double w = coeffs[i].Magnitude * Math.Exp(-modeInfo[i].rate * t);
            totW += w;
            if (modeInfo[i].band == "MID") midW += w;
            else if (modeInfo[i].band == "SLOW") slowW += w;
            else if (modeInfo[i].band == "FAST") fastW += w;
            else immW += w;
        }

        // Reconstruct density matrix from eigenbasis for CΨ
        var rhoVec = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(d2, Complex.Zero);
        for (int i = 0; i < d2; i++)
            rhoVec += coeffs[i] * Complex.Exp(eigenvalues[i] * t) * V.Column(i);

        var rhoT = Matrix<Complex>.Build.Dense(d, d);
        for (int j = 0; j < d; j++)
            for (int i = 0; i < d; i++)
                rhoT[i, j] = rhoVec[j * d + i];

        // Enforce Hermiticity (numerical cleanup)
        for (int i = 0; i < d; i++)
        {
            rhoT[i, i] = new Complex(rhoT[i, i].Real, 0);
            for (int j = i + 1; j < d; j++)
            {
                var avg = (rhoT[i, j] + Complex.Conjugate(rhoT[j, i])) / 2.0;
                rhoT[i, j] = avg;
                rhoT[j, i] = Complex.Conjugate(avg);
            }
        }

        var rho01 = DensityMatrixTools.PartialTrace(rhoT, n, new[] { 0, 1 });
        var (cpsi01, _, _) = DensityMatrixTools.ComputeCPsi(rho01);
        var rho0N = DensityMatrixTools.PartialTrace(rhoT, n, new[] { 0, n - 1 });
        var (cpsi0N, _, _) = DensityMatrixTools.ComputeCPsi(rho0N);
        double disc = 1 - 4 * cpsi01;

        string note = "";
        if (prevCpsi01 >= 0.25 && cpsi01 < 0.25 && !crossingFound)
        {
            note = " ← CΨ(01) CROSSES 1/4";
            crossingFound = true;
        }

        // Print every 0.5 or at crossing
        if (Math.Abs(t % 0.5) < 0.05 || note != "")
        {
            Console.WriteLine(FormattableString.Invariant(
                $"{t,6:F2} | {cpsi01,8:F4} | {cpsi0N,8:F4} | {disc,8:F4} | {immW/totW,6:P0} | {slowW/totW,6:P0} | {midW/totW,6:P0} | {fastW/totW,6:P0} | {note}"));
        }

        prevCpsi01 = cpsi01;
    }
}

// ============================================================
// DISPATCH
// ============================================================
if (args.Contains("pull"))
{
    Log("The Pull Principle: Scaling, Optimization, Relay");
    Log($"Started: {DateTime.Now:yyyy-MM-dd HH:mm:ss}");
    Log();

    var sw1 = Stopwatch.StartNew();
    RunPull1();
    Log($"[Part 1 completed in {sw1.Elapsed}]");
    Log();

    var sw2 = Stopwatch.StartNew();
    RunPull2();
    Log($"[Part 2 completed in {sw2.Elapsed}]");
    Log();

    var sw3 = Stopwatch.StartNew();
    RunPull3();
    Log($"[Part 3 completed in {sw3.Elapsed}]");
    Log();
}
else
{
    RunTest0();
    RunTest1();
    RunTest2();
    RunTest3();
}

Log(new string('=', 70));
Log("COMPLETE");
Log(new string('=', 70));
Log();
Log($"Total runtime: {totalSw.Elapsed}");
Log($"Results: {outPath}");
