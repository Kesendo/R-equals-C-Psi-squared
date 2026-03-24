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
/// Part 1: Scaling curve — MI vs N for uniform chain and hierarchical mediator.
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
/// Part 2: Pull optimization — 2:1 coupling at various levels.
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
/// Part 3: Relay protocol — time-dependent gamma.
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
