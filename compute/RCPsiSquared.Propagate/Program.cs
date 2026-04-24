using System.Diagnostics;
using System.Globalization;
using System.Numerics;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex;
using RCPsiSquared.Propagate;

// ============================================================
// COCKPIT MODE: scaling test of the 3-observable framework
// ============================================================
if (args.Length >= 1 && args[0] == "cockpit")
{
    try { Control.UseNativeMKL(); } catch { }
    Control.MaxDegreeOfParallelism = Environment.ProcessorCount;
    RunCockpitScaling(args);
    return;
}

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
// BRECHER MODE: variable J + variable initial state, Peak SumMI tracking
// ============================================================
if (args.Length >= 1 && args[0] == "brecher")
{
    try { Control.UseNativeMKL(); } catch { }
    Control.MaxDegreeOfParallelism = Environment.ProcessorCount;
    RunBrecherEvaluation(args);
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
// COCKPIT SCALING MODE
// ============================================================
void RunCockpitScaling(string[] cArgs)
{
    // Usage: cockpit <N> <topology> [--gamma 0.05] [--J 1.0] [--tmax 20] [--dt 0.05] [--sample 0.1]
    if (cArgs.Length < 3)
    {
        Console.Error.WriteLine("Usage: cockpit <N> <topology> [--gamma 0.05] [--J 1.0] [--tmax 20] [--dt 0.05] [--sample 0.1]");
        Environment.Exit(1);
        return;
    }

    int n = int.Parse(cArgs[1]);
    string topology = cArgs[2].ToLowerInvariant();
    if (topology != "chain" && topology != "star")
    {
        Console.Error.WriteLine($"ERROR: topology must be 'chain' or 'star', got '{topology}'");
        Environment.Exit(1);
        return;
    }

    double gamma = 0.05, J = 1.0, tMax = 20.0, dt = 0.05, sample = 0.1;
    for (int i = 3; i < cArgs.Length - 1; i++)
    {
        if (cArgs[i] == "--gamma")  gamma  = double.Parse(cArgs[i + 1], CultureInfo.InvariantCulture);
        if (cArgs[i] == "--J")      J      = double.Parse(cArgs[i + 1], CultureInfo.InvariantCulture);
        if (cArgs[i] == "--tmax")   tMax   = double.Parse(cArgs[i + 1], CultureInfo.InvariantCulture);
        if (cArgs[i] == "--dt")     dt     = double.Parse(cArgs[i + 1], CultureInfo.InvariantCulture);
        if (cArgs[i] == "--sample") sample = double.Parse(cArgs[i + 1], CultureInfo.InvariantCulture);
    }

    if (n < 5)
    {
        Console.Error.WriteLine("ERROR: V2 cockpit mode requires N >= 5");
        Environment.Exit(1);
        return;
    }

    int d = 1 << n;
    var gammas = Enumerable.Repeat(gamma, n).ToArray();

    // Build topology
    var couplings = Enumerable.Repeat(J, n - 1).ToArray();
    Bond[] bonds = topology == "chain"
        ? Topology.Chain(n, couplings)
        : Topology.Star(n, couplings);

    // V2: Bell pair at the center
    int c1 = (n - 1) / 2;
    int c2 = c1 + 1;

    // V2 pair selection with semantic labels
    (int q1, int q2, string label)[] pairs;
    if (topology == "chain")
    {
        pairs = new (int, int, string)[]
        {
            (c1, c2, $"{c1}_{c2}_center_bell"),
            (c1 - 1, c1, $"{c1 - 1}_{c1}_adjacent"),
            (0, 1, "0_1_far_edge"),
        };
    }
    else // star
    {
        pairs = new (int, int, string)[]
        {
            (c1, c2, $"{c1}_{c2}_center_bell"),
            (0, c1, $"0_{c1}_center_leaf"),
            (1, n - 1, $"1_{n - 1}_far_leaf"),
        };
    }

    // Build measure times
    int nSamples = (int)(tMax / sample) + 1;
    var measureTimes = new double[nSamples];
    for (int i = 0; i < nSamples; i++)
        measureTimes[i] = i * sample;

    // V2: Build initial state Bell+(c1,c2) ⊗ |+>^(N-2)
    Complex[] BuildInitialPsi(int nQ)
    {
        int dQ = 1 << nQ;
        int lc1 = (nQ - 1) / 2;
        int lc2 = lc1 + 1;
        int shiftC1 = nQ - 1 - lc1;
        int shiftC2 = nQ - 1 - lc2;
        var psi = new Complex[dQ];
        double norm = 1.0 / (Math.Sqrt(2) * Math.Pow(2.0, (nQ - 2) / 2.0));
        for (int idx = 0; idx < dQ; idx++)
        {
            int b1 = (idx >> shiftC1) & 1;
            int b2 = (idx >> shiftC2) & 1;
            if (b1 != b2) continue; // Bell+ requires center qubits to agree
            psi[idx] = norm;
        }
        return psi;
    }

    // CSV output setup (V2 directory)
    string csvDir = Path.GetFullPath(Path.Combine(
        Path.GetDirectoryName(AppContext.BaseDirectory)!,
        "..", "..", "..", "..", "..", "simulations", "results", "cockpit_scaling"));
    if (!Directory.Exists(csvDir))
        Directory.CreateDirectory(csvDir);
    string csvPath = Path.Combine(csvDir, $"cockpit_scaling_N{n}_{topology}.csv");

    var sw = Stopwatch.StartNew();
    Console.Error.WriteLine($"Cockpit N={n} {topology}: d={d}, pairs={pairs.Length}, samples={nSamples}, tmax={tMax}, dt={dt}");

    // Feature names for CSV header
    string[] featNames = { "phi_plus", "phi_minus", "psi_plus", "psi_minus", "purity", "svn", "concurrence", "psi_norm", "ph03" };

    using var csv = new StreamWriter(csvPath, false, System.Text.Encoding.UTF8);
    csv.WriteLine("t,pair," + string.Join(",", featNames));

    int totalSteps = nSamples;
    int heartbeatInterval = Math.Max(1, totalSteps / 10);

    if (n >= 14)
    {
        // Matrix-free path
        var mfH = new MatrixFreeHamiltonian(n, bonds);
        var psi = BuildInitialPsi(n);
        var rho0 = DensityMatrixToolsRaw.PureState(psi);
        psi = null;

        // V2 sanity check: center pair (c1,c2) must be Bell+ at t=0
        var rhoCenterT0 = DensityMatrixToolsRaw.PartialTrace(rho0, d, n, new[] { c1, c2 });
        double concT0 = DensityMatrixTools.Concurrence2Q(rhoCenterT0);
        double purT0 = DensityMatrixTools.Purity(rhoCenterT0);
        if (Math.Abs(concT0 - 1.0) > 1e-8 || Math.Abs(purT0 - 1.0) > 1e-8)
            throw new InvalidOperationException(
                $"V2 sanity fail: center pair ({c1},{c2}) at t=0 should be Bell+ " +
                $"with concurrence=1.0 purity=1.0, got conc={concT0:F6} pur={purT0:F6}");
        Console.Error.WriteLine($"  Sanity OK: center pair ({c1},{c2}) at t=0 is Bell+ (conc=1.0, pur=1.0)");

        var mfProp = new MatrixFreePropagator(mfH, gammas, n);
        int sampleIdx = 0;

        mfProp.Propagate(rho0, tMax, dt, measureTimes, (t, rho) =>
        {
            foreach (var (q1, q2, label) in pairs)
            {
                var rhoPair = DensityMatrixToolsRaw.PartialTrace(rho, d, n, new[] { q1, q2 });
                var feats = DensityMatrixTools.ExtractCockpitFeatures(rhoPair);
                csv.Write(string.Format(CultureInfo.InvariantCulture, "{0:F4},{1}",
                    t, label));
                foreach (var f in feats)
                    csv.Write(string.Format(CultureInfo.InvariantCulture, ",{0:G15}", f));
                csv.WriteLine();
            }
            csv.Flush();

            sampleIdx++;
            if (sampleIdx % heartbeatInterval == 0 || sampleIdx == totalSteps)
            {
                int pct = (int)(100.0 * sampleIdx / totalSteps);
                Console.Error.WriteLine($"Cockpit N={n} {topology}: t={t:F1} / {tMax:F1} ({pct}%) elapsed={sw.Elapsed:hh\\:mm\\:ss}");
            }
        });
    }
    else
    {
        // Dense path
        var H = Topology.BuildHamiltonian(n, bonds);
        var psi = BuildInitialPsi(n);
        var rho0 = DensityMatrixTools.PureState(psi);

        // V2 sanity check: center pair (c1,c2) must be Bell+ at t=0
        var rhoCenterT0 = DensityMatrixTools.PartialTrace(rho0, n, new[] { c1, c2 });
        double concT0 = DensityMatrixTools.Concurrence2Q(rhoCenterT0);
        double purT0 = DensityMatrixTools.Purity(rhoCenterT0);
        if (Math.Abs(concT0 - 1.0) > 1e-8 || Math.Abs(purT0 - 1.0) > 1e-8)
            throw new InvalidOperationException(
                $"V2 sanity fail: center pair ({c1},{c2}) at t=0 should be Bell+ " +
                $"with concurrence=1.0 purity=1.0, got conc={concT0:F6} pur={purT0:F6}");
        Console.Error.WriteLine($"  Sanity OK: center pair ({c1},{c2}) at t=0 is Bell+ (conc=1.0, pur=1.0)");

        var prop = new LindbladPropagator(H, gammas, n);
        int sampleIdx = 0;

        prop.Propagate(rho0, tMax, dt, measureTimes, (t, rho) =>
        {
            foreach (var (q1, q2, label) in pairs)
            {
                var rhoPair = DensityMatrixTools.PartialTrace(rho, n, new[] { q1, q2 });
                var feats = DensityMatrixTools.ExtractCockpitFeatures(rhoPair);
                csv.Write(string.Format(CultureInfo.InvariantCulture, "{0:F4},{1}",
                    t, label));
                foreach (var f in feats)
                    csv.Write(string.Format(CultureInfo.InvariantCulture, ",{0:G15}", f));
                csv.WriteLine();
            }
            csv.Flush();

            sampleIdx++;
            if (sampleIdx % heartbeatInterval == 0 || sampleIdx == totalSteps)
            {
                int pct = (int)(100.0 * sampleIdx / totalSteps);
                Console.Error.WriteLine($"Cockpit N={n} {topology}: t={t:F1} / {tMax:F1} ({pct}%) elapsed={sw.Elapsed:hh\\:mm\\:ss}");
            }
        });
    }

    sw.Stop();
    Console.WriteLine($"RESULT N={n} topology={topology} runtime={sw.Elapsed:hh\\:mm\\:ss} sample_count={nSamples} pair_count={pairs.Length} csv={csvPath}");
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

    // Heisenberg chain (J=1.0, N-1 bonds)
    var couplings = Enumerable.Repeat(1.0, n - 1).ToArray();
    var bonds = Topology.Chain(n, couplings);

    // Measure every 0.5 time units
    double mInterval = 0.5;
    int nMeas = (int)(tMax / mInterval) + 1;
    var tMeas = new double[nMeas];
    for (int i = 0; i < nMeas; i++) tMeas[i] = i * mInterval;

    double bestSumMI = -1, bestT = 0, bestPeakMI = -1;
    double cpsi01Best = 0, purBest = 0;
    double sumMI5 = 0;

    if (n >= 14)
    {
        // Matrix-free path: no Hamiltonian matrix stored.
        // Uses ~72 GB for N=15 instead of ~160 GB.
        Console.Error.WriteLine($"Matrix-free path for N={n} (d={d})");

        var mfH = new MatrixFreeHamiltonian(n, bonds);

        // Initial state: |+>^N as raw column-major Complex[]
        var psi = new Complex[d];
        double psiNorm = 1.0 / Math.Sqrt(d);
        for (int idx = 0; idx < d; idx++) psi[idx] = psiNorm;
        var rho0 = DensityMatrixToolsRaw.PureState(psi);
        psi = null; // free psi early

        var mfProp = new MatrixFreePropagator(mfH, gammas, n);

        mfProp.Propagate(rho0, tMax, dt, tMeas, (t, rho) =>
        {
            double sumMI = 0;
            for (int k = 0; k < n - 1; k++)
                sumMI += DensityMatrixToolsRaw.MutualInformation(rho, d, n,
                    new[] { k }, new[] { k + 1 });

            double mi0N = DensityMatrixToolsRaw.MutualInformation(rho, d, n,
                new[] { 0 }, new[] { n - 1 });

            if (sumMI > bestSumMI)
            {
                bestSumMI = sumMI;
                bestT = t;
                var (cpsiVal, _, _) = DensityMatrixToolsRaw.ComputeCPsi01(rho, d, n);
                cpsi01Best = cpsiVal;
                purBest = DensityMatrixToolsRaw.Purity(rho, d * d);
            }

            if (mi0N > bestPeakMI) bestPeakMI = mi0N;
            if (Math.Abs(t - 5.0) < 0.01) sumMI5 = sumMI;

            Console.Error.Write($"\r  t={t:F1} SumMI={sumMI:F4}");
        });
        Console.Error.WriteLine();
    }
    else
    {
        // Dense path: MathNet matrices, for N <= 13
        var H = Topology.BuildHamiltonian(n, bonds);

        var psi = new Complex[d];
        double psiNorm = 1.0 / Math.Sqrt(d);
        for (int idx = 0; idx < d; idx++) psi[idx] = psiNorm;
        var rho0 = DensityMatrixTools.PureState(psi);

        var prop = new LindbladPropagator(H, gammas, n);

        prop.Propagate(rho0, tMax, dt, tMeas, (t, rho) =>
        {
            double sumMI = 0;
            for (int k = 0; k < n - 1; k++)
                sumMI += DensityMatrixTools.MutualInformation(rho, n,
                    new[] { k }, new[] { k + 1 });

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
    }

    // Output: single machine-parseable line
    Console.WriteLine(FormattableString.Invariant(
        $"RESULT SumMI={bestSumMI:F6} PeakMI={bestPeakMI:F6} PeakT={bestT:F2} CPsi01={cpsi01Best:F6} Purity={purBest:F6} SumMI5={sumMI5:F6}"));
}

// ============================================================
// BRECHER EVALUATION MODE
// Variable J-profile + variable initial state, Peak SumMI over adjacent pairs.
// Designed for EQ-024 refinement Brecher-Test scaling to N >= 7.
// ============================================================
void RunBrecherEvaluation(string[] pArgs)
{
    // Usage: brecher <N> <J-profile> <initial-state> [--gamma G] [--tmax T] [--dt DT]
    //   J-profile: comma-separated N-1 doubles, e.g. "1.0,1.0,0.01,1.0,1.0"
    //   initial-state:
    //     plus                  -> |+>^N uniform superposition
    //     bits:<N-digit binary> -> Z-basis state (leftmost digit = site 0 = MSB of index)
    //     xpattern:<N chars>    -> product of |+>/|-> per char (leftmost = site 0)
    if (pArgs.Length < 4)
    {
        Console.Error.WriteLine("Usage: brecher <N> <J-profile> <initial-state> [--gamma G] [--tmax T] [--dt DT]");
        Console.Error.WriteLine("  J-profile:     comma-separated N-1 doubles");
        Console.Error.WriteLine("  initial-state: plus | bits:<N-digit bin> | xpattern:<N chars of +/->");
        Environment.Exit(1);
        return;
    }

    int n = int.Parse(pArgs[1], CultureInfo.InvariantCulture);
    int d = 1 << n;

    // Parse J-profile (N-1 values)
    var couplings = pArgs[2].Split(',')
        .Select(s => double.Parse(s, CultureInfo.InvariantCulture))
        .ToArray();

    if (couplings.Length != n - 1)
    {
        Console.Error.WriteLine($"ERROR: Expected {n - 1} J values (N-1 bonds), got {couplings.Length}");
        Environment.Exit(1);
        return;
    }

    string initialSpec = pArgs[3];

    double gammaVal = 0.05;
    double tMax = 15.0;
    double dt = 0.05;
    for (int i = 4; i < pArgs.Length - 1; i++)
    {
        if (pArgs[i] == "--gamma")
            gammaVal = double.Parse(pArgs[i + 1], CultureInfo.InvariantCulture);
        if (pArgs[i] == "--tmax")
            tMax = double.Parse(pArgs[i + 1], CultureInfo.InvariantCulture);
        if (pArgs[i] == "--dt")
            dt = double.Parse(pArgs[i + 1], CultureInfo.InvariantCulture);
    }

    var gammas = Enumerable.Repeat(gammaVal, n).ToArray();

    // RK4 stability: |lambda_max * dt| < 2.8 for RK4 stability region.
    // H-eigenvalue scale is O(max|J| * N). Conservative: dt <= 0.05 / max(1, max|J|).
    // Reduces dt automatically if user-provided dt is too large.
    double maxAbsJ = couplings.Select(j => Math.Abs(j)).DefaultIfEmpty(1.0).Max();
    double dtMax = 0.05 / Math.Max(1.0, maxAbsJ);
    if (dt > dtMax)
    {
        Console.Error.WriteLine($"Auto-dt: reducing dt from {dt} to {dtMax} for stability (max|J|={maxAbsJ:F3})");
        dt = dtMax;
    }

    Complex[] psi;
    try { psi = BuildInitialStatePsi(initialSpec, n); }
    catch (ArgumentException ex)
    {
        Console.Error.WriteLine($"ERROR: {ex.Message}");
        Environment.Exit(1);
        return;
    }

    // Measurement times: fine-grained early (catch fast transients), coarser late.
    // Python reference uses np.linspace(0.1, 15.0, 40) -> ~0.38 step; we use
    // every 0.1 up to t=2, then every 0.5 up to tMax. Guarantees we hit the
    // Peak SumMI which for SU(2)-broken initial states often sits at t ~ 0.1.
    var tMeasList = new List<double>();
    for (double t = 0.1; t < 2.0 - 1e-9; t += 0.1) tMeasList.Add(t);
    for (double t = 2.0; t < tMax + 1e-9; t += 0.5) tMeasList.Add(t);
    var tMeas = tMeasList.ToArray();

    double bestSumMI = -1, bestT = 0;
    double sumMI5 = 0;
    // End-to-end MI tracking: MI(site 0, site N-1).
    // Complements PeakSumMI (distributed adjacent-pair correlation) with a
    // direct transport measure, per Tom's MI(0, N-1) question 2026-04-24.
    double bestMI0N = -1, bestT_MI0N = 0;
    double mi0N_at_t5 = 0;
    // Mirror-pair Sum-MI (F71-symmetric multi-end transport):
    // MM = sum over k=0..(n/2-1) of MI(site k, site n-1-k).
    // At N=5 this is MI(0,4) + MI(1,3). At N=9 it is MI(0,8)+MI(1,7)+MI(2,6)+MI(3,5).
    // Captures the total F71-mirror-symmetric end-pair correlation (Tom's
    // multi-ends question 2026-04-24).
    double bestMM = -1, bestT_MM = 0;
    double mm5 = 0;
    int nMirrorPairs = n / 2; // floor(n/2): covers all non-self-mirror pairs

    var sw = Stopwatch.StartNew();

    // Brecher mode matrix-free threshold lowered to N >= 13: at N=13 dense
    // matmul costs ~8 GFLOPS^(-1) per step (d=8192, d^3 = 5e11), estimated 6h
    // per eval, versus matrix-free per-bond bit-flip application at ~15 min
    // per eval. Dense path (N <= 12) retained unchanged since faster there.
    if (n >= 13)
    {
        Console.Error.WriteLine($"Matrix-free path for N={n} (d={d})");

        var mfBonds = Topology.Chain(n, couplings);
        var mfH = new MatrixFreeHamiltonian(n, mfBonds);

        var rho0 = DensityMatrixToolsRaw.PureState(psi);
        psi = null!; // free early

        var mfProp = new MatrixFreePropagator(mfH, gammas, n);

        mfProp.Propagate(rho0, tMax, dt, tMeas, (t, rho) =>
        {
            double sumMI = 0;
            for (int k = 0; k < n - 1; k++)
                sumMI += DensityMatrixToolsRaw.MutualInformation(rho, d, n,
                    new[] { k }, new[] { k + 1 });

            double mi0N = DensityMatrixToolsRaw.MutualInformation(rho, d, n,
                new[] { 0 }, new[] { n - 1 });

            double mm = 0;
            for (int k = 0; k < nMirrorPairs; k++)
                mm += DensityMatrixToolsRaw.MutualInformation(rho, d, n,
                    new[] { k }, new[] { n - 1 - k });

            if (sumMI > bestSumMI) { bestSumMI = sumMI; bestT = t; }
            if (mi0N > bestMI0N) { bestMI0N = mi0N; bestT_MI0N = t; }
            if (mm > bestMM) { bestMM = mm; bestT_MM = t; }

            if (Math.Abs(t - 5.0) < 0.01)
            {
                sumMI5 = sumMI;
                mi0N_at_t5 = mi0N;
                mm5 = mm;
            }

            Console.Error.Write($"\r  t={t:F2} SumMI={sumMI:F4} MI0N={mi0N:F4} MM={mm:F4}");
        });
        Console.Error.WriteLine();
    }
    else
    {
        // Dense path: MathNet matrices, N <= 13
        var bonds = Topology.Chain(n, couplings);
        var H = Topology.BuildHamiltonian(n, bonds);

        var rho0 = DensityMatrixTools.PureState(psi);

        var prop = new LindbladPropagator(H, gammas, n);

        prop.Propagate(rho0, tMax, dt, tMeas, (t, rho) =>
        {
            double sumMI = 0;
            for (int k = 0; k < n - 1; k++)
                sumMI += DensityMatrixTools.MutualInformation(rho, n,
                    new[] { k }, new[] { k + 1 });

            double mi0N = DensityMatrixTools.MutualInformation(rho, n,
                new[] { 0 }, new[] { n - 1 });

            double mm = 0;
            for (int k = 0; k < nMirrorPairs; k++)
                mm += DensityMatrixTools.MutualInformation(rho, n,
                    new[] { k }, new[] { n - 1 - k });

            if (sumMI > bestSumMI) { bestSumMI = sumMI; bestT = t; }
            if (mi0N > bestMI0N) { bestMI0N = mi0N; bestT_MI0N = t; }
            if (mm > bestMM) { bestMM = mm; bestT_MM = t; }

            if (Math.Abs(t - 5.0) < 0.01)
            {
                sumMI5 = sumMI;
                mi0N_at_t5 = mi0N;
                mm5 = mm;
            }
        });
    }

    sw.Stop();

    string jStr = string.Join(",", couplings.Select(j => j.ToString("F6", CultureInfo.InvariantCulture)));
    Console.WriteLine(FormattableString.Invariant(
        $"RESULT N={n} J=[{jStr}] Initial={initialSpec} Gamma={gammaVal:F4} PeakSumMI={bestSumMI:F6} PeakT={bestT:F2} SumMI5={sumMI5:F6} PeakMI0N={bestMI0N:F6} PeakT_MI0N={bestT_MI0N:F2} MI0N5={mi0N_at_t5:F6} PeakMM={bestMM:F6} PeakT_MM={bestT_MM:F2} MM5={mm5:F6} ComputeTime={sw.Elapsed.TotalSeconds:F2}s"));
}

// Build initial state psi[d] from text spec. Throws ArgumentException on invalid input.
Complex[] BuildInitialStatePsi(string spec, int n)
{
    int d = 1 << n;

    if (spec == "plus")
    {
        // |+>^N uniform superposition
        var psi = new Complex[d];
        double norm = 1.0 / Math.Sqrt(d);
        for (int i = 0; i < d; i++) psi[i] = norm;
        return psi;
    }

    if (spec.StartsWith("bits:"))
    {
        string bits = spec.Substring(5);
        if (bits.Length != n)
            throw new ArgumentException($"bits: pattern length {bits.Length} != N {n}");
        if (!bits.All(c => c == '0' || c == '1'))
            throw new ArgumentException($"bits: pattern must contain only 0 and 1, got '{bits}'");
        // Leftmost char = site 0 = MSB of integer index
        int idx = Convert.ToInt32(bits, 2);
        var psi = new Complex[d];
        psi[idx] = 1.0;
        return psi;
    }

    if (spec.StartsWith("xpattern:"))
    {
        string pattern = spec.Substring(9);
        if (pattern.Length != n)
            throw new ArgumentException($"xpattern: pattern length {pattern.Length} != N {n}");
        // Iterative tensor product, first char = site 0 = slowest-varying = MSB of index.
        // Supported chars: '+' |+>, '-' |->, '0' |0>, '1' |1> (mixed X/Z basis products).
        double invSqrt2 = 1.0 / Math.Sqrt(2);
        Complex[] psi = new Complex[] { 1.0 };
        foreach (char c in pattern)
        {
            Complex[] ket;
            if (c == '+') ket = new Complex[] { invSqrt2, invSqrt2 };
            else if (c == '-') ket = new Complex[] { invSqrt2, -invSqrt2 };
            else if (c == '0') ket = new Complex[] { 1.0, 0.0 };
            else if (c == '1') ket = new Complex[] { 0.0, 1.0 };
            else throw new ArgumentException($"xpattern: invalid char '{c}', expected +, -, 0, or 1");

            var newPsi = new Complex[psi.Length * 2];
            for (int i = 0; i < psi.Length; i++)
            {
                newPsi[2 * i] = psi[i] * ket[0];
                newPsi[2 * i + 1] = psi[i] * ket[1];
            }
            psi = newPsi;
        }
        return psi;
    }

    if (spec.StartsWith("bonding:"))
    {
        // F67 single-excitation sinusoidal mode (uniform-J free-particle eigenstate):
        //   |psi_k> = sqrt(2/(N+1)) * sum_{j=0..N-1} sin(pi*k*(j+1)/(N+1)) |1_j>
        // where |1_j> has bit j set (site 0 = MSB of index).
        // Mirror-symmetric for any k (sin-profile under j <-> N-1-j).
        // k=1 is the longest-lived bonding mode, with decay alpha_1 = (4 gamma_0/(N+1)) sin^2(pi/(N+1)).
        string kStr = spec.Substring(8);
        int k;
        if (!int.TryParse(kStr, NumberStyles.Integer, CultureInfo.InvariantCulture, out k))
            throw new ArgumentException($"bonding: expected integer k, got '{kStr}'");
        if (k < 1 || k > n)
            throw new ArgumentException($"bonding: k={k} must be in [1, N={n}]");
        var psi = new Complex[d];
        double norm = Math.Sqrt(2.0 / (n + 1));
        for (int j = 0; j < n; j++)
        {
            double amp = norm * Math.Sin(Math.PI * k * (j + 1) / (n + 1));
            int idx = 1 << (n - 1 - j); // single-excitation basis state |1_j>
            psi[idx] = amp;
        }
        return psi;
    }

    throw new ArgumentException($"Unknown initial-state spec '{spec}'. Valid: plus | bits:<bin> | xpattern:<+/-/0/1> | bonding:<k>");
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
