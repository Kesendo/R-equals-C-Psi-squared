using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.Foundation;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Cli.Commands;

/// <summary>The <c>mirror</c> verb: call the conductor's stand with parameters and listen to its
/// voices. Builds a <see cref="MirrorSystem"/> and exports the live F-formula readings.
///
/// <para><b>Single reading</b>: one system , its slow modes as channel-difference portfolios and
/// the F1 palindrome check.</para>
///
/// <para><b>Time</b> (<c>--evolve K</c>): unroll that spectrum over K carrier-ticks and show how much
/// of each mode is still remembered , the connection from the static spectrum to the measured decay.
/// K is the dimensionless time (t and γ cancel into K = γ₀·t); only K is readable from inside.</para>
///
/// <para><b>Memory</b> (<c>--memory</c>): the 90° rotation H ↔ M (F80). For truly H (XY, Heisenberg)
/// the mirror is perfect, M = 0; for the non-truly bond Hamiltonian (<c>--htype xybond</c>,
/// H = Σ X_lY_{l+1}) the defect carries H's spectrum rotated a quarter turn onto the imaginary memory
/// axis (Spec(M) = ±2i·Spec(H)) , the wave remembers across the turn.</para>
///
/// <para><b>Sweep</b> (<c>--sweep-J start,stop,steps</c>): move the coupling J in steps with the
/// carrier fixed, and watch the slowest (memory) mode's portfolio move , the box in motion, the
/// trajectory in classification-space. At J=0 the slowest mode is pure in the slow channel
/// (protected); as J grows it leaks into the fast neighbours. The path is the reading.</para>
///
/// Call parameters in, readings out (stdout, optional CSV via <c>--out</c>).</summary>
public static class MirrorCommand
{
    public static int Run(string[] args)
    {
        int N = GetInt(args, "--N", 0);
        if (N < 1 || N > 7) { Console.Error.WriteLine("mirror: --N required, 1..7 (spectrum via joint-popcount blocks; --memory and non-conserving H stay dense)"); return 2; }

        string htypeStr = (GetStr(args, "--htype", "XY") ?? "XY").ToLowerInvariant();
        var htype = htypeStr switch
        {
            "heisenberg" or "heis" => HamiltonianType.Heisenberg,
            _ => HamiltonianType.XY,
        };
        var topo = (GetStr(args, "--topology", "chain") ?? "chain").ToLowerInvariant() switch
        {
            "star" => TopologyKind.Star,
            "ring" => TopologyKind.Ring,
            _ => TopologyKind.Chain,
        };

        double[] gammas;
        string? gList = GetStr(args, "--gamma-list", null);
        if (gList is not null)
        {
            gammas = gList.Split(',').Select(s => double.Parse(s.Trim(), CultureInfo.InvariantCulture)).ToArray();
            if (gammas.Length != N) { Console.Error.WriteLine($"mirror: --gamma-list needs {N} values, got {gammas.Length}"); return 2; }
        }
        else
        {
            double g = GetDouble(args, "--gamma", 0.1);
            gammas = Enumerable.Repeat(g, N).ToArray();
        }
        var channels = gammas.Select((g, l) => new ChannelRate($"q{l}", g)).ToList();
        var inv = CultureInfo.InvariantCulture;
        string? outPath = GetStr(args, "--out", null);

        ComplexMatrix BuildHamiltonian(double J) => htypeStr == "xybond"
            ? NonTrulyXYBond(N, J)
            : new ChainSystem(N, J, gammas[0], htype, topo).BuildHamiltonian();
        MirrorSystem Build(double J) => new(N, BuildHamiltonian(J), channels);

        string? sweep = GetStr(args, "--sweep-J", null);
        if (sweep is not null)
            return SweepJ(sweep, N, htype, topo, gammas, channels, Build, outPath, inv);

        // ---- single reading ----
        double jSingle = GetDouble(args, "--J", 1.0);
        int top = GetInt(args, "--top", 6);
        var sys = Build(jSingle);

        // The spectrum/F1/evolve voices need the full 4^N Liouvillian eigendecomposition (a
        // sequential non-Hermitian Evd, the heavy part at large N). --memory does not: it builds M
        // by matmul and only eigendecomposes H. --no-spectrum skips the heavy reading.
        bool noSpectrum = args.Contains("--no-spectrum");

        Console.WriteLine($"# mirror: N={N} J={jSingle.ToString(inv)} H={htype} topology={topo}");
        Console.WriteLine($"# carrier gamma = [{string.Join(", ", gammas.Select(g => g.ToString("0.###", inv)))}]  sigma = {sys.TotalDephasing.ToString("0.####", inv)}");
        Console.WriteLine();
        if (!noSpectrum)
        {
            Console.WriteLine($"Spectrum (slowest {top} nonzero modes, rate + angle theta=arctan(omega/rate) (deg) + per-channel difference portfolio):");
            foreach (var grp in sys.Spectrum.Modes
                         .Where(m => m.ActualDecayRate > 1e-9)
                         .GroupBy(m => m.ActualDecayRate.ToString("0.000000", inv))
                         .OrderBy(g => g.First().ActualDecayRate)
                         .Take(top))
            {
                // A degenerate rate can hold modes at different angles (some rotate, some do
                // not); show the most-rotating one, so the slowest line matches the Rotation
                // voice's gap reading (which takes max|omega| over the modes at the gap).
                var rep = grp.OrderByDescending(m => Math.Abs(m.OscillationFrequency)).First();
                double theta = Math.Atan2(Math.Abs(rep.OscillationFrequency), rep.ActualDecayRate) * 180.0 / Math.PI;
                Console.WriteLine($"  rate {rep.ActualDecayRate.ToString("0.0000", inv),9}   angle {theta.ToString("0.0", inv),5} deg   {Portfolio(rep, inv)}");
            }
            Console.WriteLine();
            Console.WriteLine($"F1 palindrome (rate r pairs with 2*sigma - r): holds = {sys.PalindromeHolds}");
            var sample = sys.PalindromePartners.FirstOrDefault(p => p.Rate > 1e-9);
            if (sample is not null)
                Console.WriteLine($"  e.g. rate {sample.Rate.ToString("0.0000", inv)} pairs with {sample.PartnerRate.ToString("0.0000", inv)} (partner present: {sample.PartnerPresent})");

            // --evolve K: the static spectrum unrolled into time (K = carrier-ticks; t and gamma cancel into K).
            string? evStr = GetStr(args, "--evolve", null);
            if (evStr is not null && double.TryParse(evStr, NumberStyles.Float, inv, out double K))
            {
                Console.WriteLine();
                Console.WriteLine($"At K = {K.ToString(inv)} carrier-ticks  (survival = e^(-(rate/sigma)*K), sigma = {sys.TotalDephasing.ToString("0.####", inv)}):");
                Console.WriteLine("  the spectrum unrolled into time , how much of each mode is still remembered after K ticks");
                var seenK = new HashSet<string>();
                foreach (var s in sys.Evolve(K).Where(s => s.Rate > 1e-9).OrderBy(s => s.Rate))
                {
                    var key = s.Rate.ToString("0.000000", inv);
                    if (!seenK.Add(key)) continue;
                    Console.WriteLine($"  rate {s.Rate.ToString("0.0000", inv),9}   survival {s.Survival.ToString("0.0000", inv),7}   {PortfolioOf(s.Portfolio, inv)}");
                    if (seenK.Count >= top) break;
                }
            }

            // The clock, both hands on the memory mode: Takt the radial hand (decay, set by gamma),
            // Rotation the angular hand (turning omega, angle theta = arctan(omega/gap), set by J). e^(lambda t) =
            // e^(-alpha t) * e^(i omega t): the radius winds in while the angle turns round.
            var takt = sys.Takt;
            var rot = sys.Rotation;
            string tau = double.IsPositiveInfinity(takt.Tau) ? "inf" : takt.Tau.ToString("0.0000", inv);
            Console.WriteLine();
            Console.WriteLine("Clock (the two hands, both on the memory mode):");
            Console.WriteLine($"  Takt (radial)      stopped {takt.Stopped}   gap {takt.Gap.ToString("0.0000", inv)}   tau {tau}");
            Console.WriteLine($"  Rotation (angular) turning {rot.Turning}   omega {rot.Frequency.ToString("0.0000", inv)}   angle {(rot.Angle * 180.0 / Math.PI).ToString("0.0", inv)} deg  (theta = arctan(omega/gap); arctan Q at 2-level)");
        }

        // --memory: the 90 degree memory rotation (F80), H <-> M.
        if (args.Contains("--memory"))
        {
            var mr = sys.MemoryRotation;
            Console.WriteLine();
            Console.WriteLine("Memory rotation (F80, the 90 degree turn H <-> M):");
            Console.WriteLine($"  mirror-defect ||M||_F = {mr.DefectNorm.ToString("0.0000", inv)}   perfect mirror: {mr.PerfectMirror}   memory carries energy: {mr.MemoryCarriesEnergy}");
            if (mr.PerfectMirror)
                Console.WriteLine("  truly H: the palindrome holds exactly, M=0, no memory-defect (a perfect mirror)");
            else if (mr.MemoryCarriesEnergy)
            {
                Console.WriteLine("  non-truly: the defect carries H's spectrum, rotated 90 degrees onto the imaginary memory axis:");
                var seenE = new HashSet<string>();
                foreach (var p in mr.Rotation.OrderByDescending(p => Math.Abs(p.Energy)))
                {
                    var key = p.Energy.ToString("0.0000", inv);
                    if (!seenE.Add(key)) continue;
                    Console.WriteLine($"    energy {p.Energy.ToString("0.0000", inv),9}   ->   memory {p.MemoryAxisValue.ToString("0.0000", inv),9} i");
                    if (seenE.Count >= top) break;
                }
            }
        }

        if (!noSpectrum && outPath is not null)
        {
            using var w = new StreamWriter(outPath);
            w.WriteLine("rate," + string.Join(",", channels.Select(c => c.Channel)));
            foreach (var m in sys.Spectrum.Modes.OrderBy(m => m.ActualDecayRate))
                w.WriteLine(m.ActualDecayRate.ToString("0.000000", inv) + "," +
                    string.Join(",", m.Portfolio.Activity.Select(a => a.Delta.ToString("0.000000", inv))));
            Console.Error.WriteLine($"# exported {sys.Spectrum.Modes.Count} modes to {outPath}");
        }
        return 0;
    }

    private static int SweepJ(string spec, int N, HamiltonianType htype, TopologyKind topo,
        double[] gammas, IReadOnlyList<ChannelRate> channels, Func<double, MirrorSystem> build,
        string? outPath, CultureInfo inv)
    {
        var parts = spec.Split(',');
        if (parts.Length != 3) { Console.Error.WriteLine("mirror: --sweep-J needs start,stop,steps (e.g. 0,2,9)"); return 2; }
        double start = double.Parse(parts[0], NumberStyles.Float, inv);
        double stop = double.Parse(parts[1], NumberStyles.Float, inv);
        int steps = int.Parse(parts[2], inv);
        if (steps < 1) { Console.Error.WriteLine("mirror: --sweep-J steps must be >= 1"); return 2; }

        Console.WriteLine($"# mirror sweep: J in [{start.ToString(inv)}, {stop.ToString(inv)}], {steps} steps");
        Console.WriteLine($"# N={N} H={htype} topology={topo}  carrier gamma = [{string.Join(", ", gammas.Select(g => g.ToString("0.###", inv)))}]");
        Console.WriteLine($"# trajectory of the slowest (memory) mode: how its portfolio moves as the coupling grows");
        Console.WriteLine();
        Console.WriteLine($"  {"J",8}  {"rate",9}   slowest-mode portfolio");

        StreamWriter? w = outPath is not null ? new StreamWriter(outPath) : null;
        w?.WriteLine("J,rate," + string.Join(",", channels.Select(c => c.Channel)));
        try
        {
            for (int i = 0; i < steps; i++)
            {
                double J = steps == 1 ? start : start + (stop - start) * i / (steps - 1);
                var spectrum = build(J).Spectrum;
                var slow = spectrum.Modes
                    .Where(m => m.ActualDecayRate > 1e-9).OrderBy(m => m.ActualDecayRate).FirstOrDefault();
                if (slow is null) { Console.WriteLine($"  {J.ToString("0.000", inv),8}  (no nonzero mode)"); continue; }
                double res = Math.Abs(slow.ActualDecayRate - slow.Portfolio.DecayRate(spectrum.Carrier));
                Console.WriteLine($"  {J.ToString("0.000", inv),8}  {slow.ActualDecayRate.ToString("0.0000", inv),9}   {Portfolio(slow, inv)}   (law res {res.ToString("E1", inv)})");
                w?.WriteLine(J.ToString("0.000000", inv) + "," + slow.ActualDecayRate.ToString("0.000000", inv) + "," +
                    string.Join(",", slow.Portfolio.Activity.Select(a => a.Delta.ToString("0.000000", inv))));
            }
        }
        finally { w?.Dispose(); }
        if (outPath is not null) Console.Error.WriteLine($"# exported trajectory to {outPath}");
        return 0;
    }

    private static string Portfolio(CarrierMode m, CultureInfo inv) => PortfolioOf(m.Portfolio, inv);

    private static string PortfolioOf(ChannelDifferencePortfolio p, CultureInfo inv) =>
        string.Join("  ", p.Activity.Select(a => $"{a.Channel} {(100 * a.Delta).ToString("0", inv),3}%"));

    private static readonly ComplexMatrix XHalf = ComplexMatrix.Build.DenseOfArray(new Complex[,] { { 0, 1 }, { 1, 0 } });
    private static readonly ComplexMatrix YHalf = ComplexMatrix.Build.DenseOfArray(new Complex[,] { { 0, -Complex.ImaginaryOne }, { Complex.ImaginaryOne, 0 } });
    private static readonly ComplexMatrix IHalf = ComplexMatrix.Build.DenseIdentity(2);

    // Non-truly chain Pi2-odd Hamiltonian H = J * sum_l X_l Y_{l+1}: the F80 case where M != 0,
    // so --memory shows the 90 degree rotation carrying H's spectrum instead of a perfect mirror.
    private static ComplexMatrix NonTrulyXYBond(int N, double J)
    {
        var H = ComplexMatrix.Build.Dense(1 << N, 1 << N);
        for (int l = 0; l < N - 1; l++)
        {
            ComplexMatrix? term = null;
            for (int s = 0; s < N; s++)
            {
                var op = s == l ? XHalf : (s == l + 1 ? YHalf : IHalf);
                term = term is null ? op : term.KroneckerProduct(op);
            }
            H += term!.Multiply((Complex)J);
        }
        return H;
    }

    private static int GetInt(string[] a, string k, int def)
    { int i = Array.IndexOf(a, k); return (i >= 0 && i + 1 < a.Length && int.TryParse(a[i + 1], out var v)) ? v : def; }

    private static double GetDouble(string[] a, string k, double def)
    { int i = Array.IndexOf(a, k); return (i >= 0 && i + 1 < a.Length && double.TryParse(a[i + 1], NumberStyles.Float, CultureInfo.InvariantCulture, out var v)) ? v : def; }

    private static string? GetStr(string[] a, string k, string? def)
    { int i = Array.IndexOf(a, k); return (i >= 0 && i + 1 < a.Length) ? a[i + 1] : def; }
}
