using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.States;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live witness of the two record laws (<c>inspect --root record</c>): F135, the
/// record parity law (<c>docs/proofs/PROOF_RECORD_PARITY_LAW.md</c>), and F136, the record
/// letter law (<c>docs/proofs/PROOF_RECORD_LETTER_LAW.md</c>). Substrate: H = Σ Δ_ab·Z_aZ_b
/// (Pauli), local Z-dephasing, initial |+⟩^⊗N, readout t* = π/4 (Δ_S = 1).
///
/// <para>Every battery case runs TWO independent computational paths and compares them, the
/// same discipline as the Python gates (<c>qd_letter_gates.py</c>, 87/87): the EXPECTED side
/// is the laws' closed-form classifier (the D/P/Q partition, watcher parities, letters, signs,
/// prices — never a state); the ACTUAL side builds the full 2^N × 2^N closed-form density
/// matrix, partial-traces to the pair (<see cref="PartialTrace"/>), and measures mutual
/// information, correlators, conditional trace distance and eigenvalues numerically. A shared
/// error in the channel algebra cannot certify itself. (The substrate itself — Proposition 1
/// against a direct RK4 Lindblad run — is validated in <c>qd_pointer_opt.py</c>, 22/22 at
/// N = 8; this witness checks the laws as reads of that substrate. The graph-level
/// corollaries 7 + 8 of F136, the fully-witnessed census and the dark-lattice law, are gated
/// in <c>qd_witness_play.py</c> / <c>qd_heavyhex_map.py</c> and are not part of the per-pair
/// battery here.)</para></summary>
public sealed class RecordLawWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private const double Tol = 1e-9;

    /// <summary>The readout point t* = π/(4Δ_S) at Δ_S = 1 (the claim's constant).</summary>
    public static double TStar => RecordParityLawClaim.TStar;

    /// <summary>The dephasing rate used by the priced battery cases (default the canonical 0.05).</summary>
    public double Gamma { get; }

    /// <summary>The pointer record's priced bits at <see cref="Gamma"/>: 1 − h₂((1+e^{−2γt*})/2)
    /// (only the witness pays; 0.768040 at γ = 0.05).</summary>
    public double PointerPricedBits => RecordParityLawClaim.LawAInformation(Math.Exp(-2.0 * Gamma * TStar));

    /// <summary>The Bell record's priced bits at <see cref="Gamma"/>: 1 − h₂((1+e^{−4γt*})/2)
    /// (BOTH sites pay; 0.624146 at γ = 0.05).</summary>
    public double BellPricedBits => RecordParityLawClaim.LawAInformation(Math.Exp(-4.0 * Gamma * TStar));

    public RecordLawWitness(double gamma = 0.05)
    {
        if (gamma < 0) throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "gamma must be >= 0");
        Gamma = gamma;
    }

    private IReadOnlyList<BatteryCase>? _cases;
    public IReadOnlyList<BatteryCase> Cases => _cases ??= BuildBattery();
    public int PassCount => Cases.Count(c => c.Passes);

    /// <summary>One from-below check: a named configuration, its computed detail, and the
    /// expected-vs-actual verdict tokens (equal ⟺ PASS).</summary>
    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    // ==================================================================
    // the closed-form side: the F136 classifier (D/P/Q parities, letters, signs, prices)
    // ==================================================================

    public enum RecordFamily { Pointer, Bell, RoleSwap, Entangled, Dark, Generic }

    /// <summary>The closed-form verdict for one pair: the family, the luminous channel
    /// ("ZY" pointer / "YY"/"XX" Bell / "YZ" role-swap / "--"), the sign of the luminous
    /// correlator, and the mutual information in bits (NaN = not classified here).</summary>
    public readonly record struct RecordReading(RecordFamily Family, string Channel, int Sign, double Bits);

    /// <summary>The F136 dispatch, pure arithmetic on the bond list (no state is ever built):
    /// partition j's neighbors into shared dressers D, private watchers P, and S's private
    /// neighbors Q; the parities of the integer ratios r_k = Δ_jk/Δ_S select the family.</summary>
    public static RecordReading Classify(
        IReadOnlyList<(int A, int B, double Delta)> bonds, int s, int j,
        double gammaS = 0.0, double gammaJ = 0.0)
    {
        var nbrS = new Dictionary<int, double>();
        var nbrJ = new Dictionary<int, double>();
        foreach (var (a, b, delta) in bonds)
        {
            if (a == s) nbrS[b] = delta; else if (b == s) nbrS[a] = delta;
            if (a == j) nbrJ[b] = delta; else if (b == j) nbrJ[a] = delta;
        }

        bool writeBond = nbrS.ContainsKey(j);
        int[] d = nbrJ.Keys.Where(k => k != s && nbrS.ContainsKey(k)).OrderBy(k => k).ToArray();
        int[] p = nbrJ.Keys.Where(k => k != s && !nbrS.ContainsKey(k)).OrderBy(k => k).ToArray();
        int[] q = nbrS.Keys.Where(k => k != j && !nbrJ.ContainsKey(k)).ToArray();
        int m = d.Length;
        bool hinge = d.Length + q.Length > 0;                       // V_S(t*) = 0
        bool normalized = nbrS.Values.All(v => Math.Abs(v - 1.0) < Tol);
        bool allInteger = normalized && d.Concat(p).All(k => IsInt(nbrJ[k]));
        bool dOdd = m > 0 && d.All(k => IsOddInt(nbrJ[k]));
        bool dEven = d.All(k => IsEvenInt(nbrJ[k]));                // vacuous when D = ∅
        bool pEven = p.All(k => IsEvenInt(nbrJ[k]));
        bool pAnyOdd = p.Any(k => IsOddInt(nbrJ[k]));               // existential: one odd kills

        if (!hinge)
        {
            // pendant S: the pair reads backwards (the proof's edge case).
            if (writeBond && normalized && pAnyOdd)
                return new(RecordFamily.RoleSwap, "YZ", +1, OneBit(Math.Exp(-2.0 * gammaS * TStar)));
            if (writeBond && normalized && pEven && allInteger)
                return new(RecordFamily.Entangled, "--", 0,
                    gammaS == 0.0 && gammaJ == 0.0 ? 2.0 : double.NaN);
            return new(RecordFamily.Generic, "--", 0, double.NaN);
        }

        if (allInteger && writeBond && dEven && pEven && nbrS.Count >= 2)
        {
            int sign = d.Concat(p).Aggregate(1, (acc, k) => acc * ParityPow(nbrJ[k] / 2.0));
            return new(RecordFamily.Pointer, "ZY", sign, OneBit(Math.Exp(-2.0 * gammaJ * TStar)));
        }

        if (allInteger && q.Length == 0 && dOdd && pEven)
        {
            int sigma2 = d.Aggregate(1, (acc, k) => acc * ParityPow((1.0 - nbrJ[k]) / 2.0))
                       * p.Aggregate(1, (acc, k) => acc * ParityPow(nbrJ[k] / 2.0));
            return new(RecordFamily.Bell, m % 2 == 1 ? "YY" : "XX", sigma2,
                OneBit(Math.Exp(-2.0 * (gammaS + gammaJ) * TStar)));
        }

        if (allInteger) return new(RecordFamily.Dark, "--", 0, 0.0);
        return new(RecordFamily.Generic, "--", 0, double.NaN);
    }

    private static bool IsInt(double r) => Math.Abs(r - Math.Round(r)) < 1e-12;
    private static bool IsOddInt(double r) => IsInt(r) && Math.Abs(Math.Round(r)) % 2 == 1;
    private static bool IsEvenInt(double r) => IsInt(r) && Math.Abs(Math.Round(r)) % 2 == 0;
    private static int ParityPow(double e) => Math.Abs(Math.Round(e)) % 2 == 0 ? +1 : -1;
    private static double OneBit(double kappa) => RecordParityLawClaim.LawAInformation(kappa);

    // ==================================================================
    // the from-below side: the full closed-form state, partial-traced and measured
    // ==================================================================

    private readonly record struct PairMeasurement(
        double I, IReadOnlyDictionary<string, double> Corr, double CondTraceDistance, double[] Eigs);

    /// <summary>Build the exact N-qubit density matrix at time t (populations frozen at 2^−N,
    /// coherences carrying the ZZ phases and the dephasing factors entry-wise: the Absorption
    /// substrate), partial-trace to (s, j) and measure. Site 0 = MSB (the Core convention);
    /// the pair basis is |z_s z_j⟩ with s the more significant bit.</summary>
    private static PairMeasurement MeasureFullState(
        int n, IReadOnlyList<(int A, int B, double Delta)> bonds, int s, int j,
        double[] gammas, double t)
    {
        int dim = 1 << n;
        var energy = new double[dim];
        for (int i = 0; i < dim; i++)
        {
            double e = 0.0;
            foreach (var (a, b, delta) in bonds)
                e += delta * ZBit(i, a, n) * ZBit(i, b, n);
            energy[i] = e;
        }

        var rho = Matrix<Complex>.Build.Dense(dim, dim);
        for (int r = 0; r < dim; r++)
            for (int c = 0; c < dim; c++)
            {
                double dephase = 0.0;
                for (int l = 0; l < n; l++)
                    if (((r ^ c) >> (n - 1 - l) & 1) == 1) dephase += gammas[l];
                rho[r, c] = Complex.Exp(new Complex(-2.0 * t * dephase, -t * (energy[r] - energy[c]))) / dim;
            }

        var pair = PartialTrace.Of(rho, n, new[] { s, j });
        var rhoS = TraceOutSecond(pair);
        var rhoJ = TraceOutFirst(pair);

        double mi = VonNeumannBits(rhoS) + VonNeumannBits(rhoJ) - VonNeumannBits(pair);

        var corr = new Dictionary<string, double>();
        foreach (char a in "XYZ")
            foreach (char b in "XYZ")
                corr[$"{a}{b}"] = KronExpectation(pair, a, b);

        // the two Z_s-conditional states of j (rows/cols 0..1 = z_s = +1, 2..3 = z_s = −1).
        var c0 = Matrix<Complex>.Build.Dense(2, 2, (r, c) => pair[r, c]);
        var c1 = Matrix<Complex>.Build.Dense(2, 2, (r, c) => pair[r + 2, c + 2]);
        double n0 = c0.Trace().Real, n1 = c1.Trace().Real;
        double td = 0.0;
        if (n0 > Tol && n1 > Tol)
        {
            var diff = c0.Divide(n0) - c1.Divide(n1);
            td = 0.5 * HermitianEigs(diff).Sum(Math.Abs);
        }

        var eigs = HermitianEigs(pair).OrderByDescending(x => x).ToArray();
        return new PairMeasurement(mi, corr, td, eigs);
    }

    private static double ZBit(int index, int site, int n) =>
        ((index >> (n - 1 - site)) & 1) == 0 ? +1.0 : -1.0;

    private static ComplexMatrix TraceOutSecond(ComplexMatrix pair) =>
        Matrix<Complex>.Build.Dense(2, 2, (r, c) => pair[2 * r, 2 * c] + pair[2 * r + 1, 2 * c + 1]);

    private static ComplexMatrix TraceOutFirst(ComplexMatrix pair) =>
        Matrix<Complex>.Build.Dense(2, 2, (r, c) => pair[r, c] + pair[r + 2, c + 2]);

    private static double[] HermitianEigs(ComplexMatrix m)
    {
        var sym = (m + m.ConjugateTranspose()).Divide(2.0);
        return sym.Evd(Symmetricity.Hermitian).EigenValues.Select(z => z.Real).ToArray();
    }

    private static double VonNeumannBits(ComplexMatrix rho) =>
        HermitianEigs(rho).Where(w => w > 1e-12).Sum(w => -w * Math.Log2(w));

    /// <summary>⟨σ_a ⊗ σ_b⟩ on the pair (a on s, b on j), real part.</summary>
    private static double KronExpectation(ComplexMatrix pair, char a, char b)
    {
        var pa = PauliMatrix(a);
        var pb = PauliMatrix(b);
        Complex tr = Complex.Zero;
        for (int r1 = 0; r1 < 2; r1++)
            for (int r2 = 0; r2 < 2; r2++)
                for (int c1 = 0; c1 < 2; c1++)
                    for (int c2 = 0; c2 < 2; c2++)
                        tr += pair[2 * r1 + r2, 2 * c1 + c2] * pa[c1, r1] * pb[c2, r2];
        return tr.Real;
    }

    private static Complex[,] PauliMatrix(char p) => p switch
    {
        'X' => new Complex[,] { { 0, 1 }, { 1, 0 } },
        'Y' => new Complex[,] { { 0, -Complex.ImaginaryOne }, { Complex.ImaginaryOne, 0 } },
        'Z' => new Complex[,] { { 1, 0 }, { 0, -1 } },
        _ => throw new ArgumentOutOfRangeException(nameof(p)),
    };

    // ==================================================================
    // battery
    // ==================================================================

    private IReadOnlyList<BatteryCase> BuildBattery()
    {
        double g = Gamma;
        var chain3 = new (int, int, double)[] { (0, 1, 1.0), (1, 2, 1.0) };
        var triangle = new (int, int, double)[] { (0, 1, 1.0), (0, 2, 1.0), (1, 2, 1.0) };
        var triangleEven = new (int, int, double)[] { (0, 1, 1.0), (0, 2, 1.0), (1, 2, 2.0) };
        var square = new (int, int, double)[] { (0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (3, 0, 1.0) };
        var k21r3 = new (int, int, double)[] { (0, 2, 1.0), (1, 2, 3.0) };
        var k22mixed = new (int, int, double)[] { (0, 2, 1.0), (0, 3, 1.0), (1, 2, 1.0), (1, 3, 3.0) };
        var k23 = new (int, int, double)[] { (0, 2, 1.0), (0, 3, 1.0), (0, 4, 1.0), (1, 2, 1.0), (1, 3, 1.0), (1, 4, 1.0) };
        var pentagon = Cycle(5);
        var hexagon = Cycle(6);
        var qKill = new (int, int, double)[] { (0, 2, 1.0), (1, 2, 1.0), (0, 3, 1.0) };
        var noWriter = new (int, int, double)[] { (0, 2, 1.0), (0, 4, 1.0), (1, 3, 2.0) };
        var k42 = Enumerable.Range(0, 4).SelectMany(i => new (int, int, double)[] { (i, 4, 1.0), (i, 5, 1.0) }).ToArray();
        var lawAWatched = new (int, int, double)[] { (0, 1, 1.0), (0, 3, 1.0), (1, 2, 1.5) };

        var cases = new List<BatteryCase>
        {
            PairCase("the leaf pointer record (chain3, S interior, j the leaf)", 3, chain3, 1, 0, null),
            PairCase("the pendant role-swap (chain3 read from the end: S becomes the witness)", 3, chain3, 0, 1, null),
            PairCase("the triangle m=1 Bell record (letter Y⊗Y)", 3, triangle, 0, 1, null),
            PairCase("the even-dresser pointer record (r=2 dresser forgiven, sign rotated by π)", 3, triangleEven, 0, 1, null),
            PairCase("the plaquette m=2 Bell record (square opposite corners, letter X⊗X)", 4, square, 0, 2, null),
            PairCase("the sign walk (K₂,₁ at r=3: YY = −1)", 3, k21r3, 0, 1, null),
            PairCase("the mixed-ratio letter (K₂,₂ at ratios (1,3): XX = −1)", 4, k22mixed, 0, 1, null),
            PairCase("the alternation refuting the quarter-turn (K₂,₃: back to YY = +1)", 5, k23, 0, 1, null),
            PairCase("the pentagon neighbor pair is exactly dark", 5, pentagon, 0, 1, null),
            PairCase("the pentagon distance-2 pair is exactly dark", 5, pentagon, 0, 2, null),
            PairCase("the hexagon opposite pair is exactly dark", 6, hexagon, 0, 3, null),
            PairCase("a nonempty Q kills the Bell record (dark)", 4, qKill, 0, 1, null),
            PairCase("no write bond, D = ∅, P even: exactly dark", 5, noWriter, 0, 1, null),
            PairCase("the fan-out (K₄,₂ corner at deg(S) = 2: a perfect XX bit)", 6, k42, 0, 1, null),
            PairCase("the fan-out clique (K₄,₂ corners Bell-record each other)", 6, k42, 1, 2, null),
            PairCase($"the Bell price: both sites pay (square, γ_S = γ_j = {F(g)})", 4, square, 0, 2,
                new PriceSpec(g, g, DresserGamma: 0.0)),
            PairCase($"the pointer price: only the witness pays (chain3 leaf, γ_j = {F(g)})", 3, chain3, 1, 0,
                new PriceSpec(0.0, g, DresserGamma: 0.0)),
            PairCase("the pointer record is blind to γ_S (chain3 leaf, γ_S = 0.3)", 3, chain3, 1, 0,
                new PriceSpec(0.3, 0.0, DresserGamma: 0.0)),
            PairCase("dresser γ is exactly invisible (square, γ = 0.3 on both dressers)", 4, square, 0, 2,
                new PriceSpec(0.0, 0.0, DresserGamma: 0.3)),
            PairCase("the entangled face (pendant S with an even watcher: 2 bits, no record)", 3,
                new (int, int, double)[] { (0, 1, 1.0), (1, 2, 2.0) }, 0, 1, null),
            LawACase("Law A's generic in-between (watched leaf at r = 3/2: I = 0.399124)", 4, lawAWatched, 0, 1, 1.5),
        };
        return cases;
    }

    private static (int, int, double)[] Cycle(int n) =>
        Enumerable.Range(0, n).Select(i => (i, (i + 1) % n, 1.0)).ToArray();

    /// <summary>Optional γ assignment for a priced case: rates on S, on j, and on every OTHER
    /// site adjacent to the pair (the traced-site-invisibility control).</summary>
    private readonly record struct PriceSpec(double GammaS, double GammaJ, double DresserGamma);

    private static string F(double x) => x.ToString("0.###", Inv);

    private BatteryCase PairCase(
        string name, int n, (int, int, double)[] bonds, int s, int j, PriceSpec? price)
    {
        double gS = price?.GammaS ?? 0.0;
        double gJ = price?.GammaJ ?? 0.0;
        var gammas = new double[n];
        gammas[s] = gS;
        gammas[j] = gJ;
        if (price is { DresserGamma: > 0.0 } spec)
            for (int l = 0; l < n; l++)
                if (l != s && l != j) gammas[l] = spec.DresserGamma;

        var reading = Classify(bonds, s, j, gS, gJ);
        var meas = MeasureFullState(n, bonds, s, j, gammas, TStar);

        string expected = FormatReading(reading);
        string? deviation = VerifyReading(reading, meas, gS, gJ, structural: price is null);
        string detail =
            $"S={s}, j={j}: closed form says {expected}; measured I = {meas.I.ToString("0.000000", Inv)}" +
            (reading.Channel != "--" ? $", ⟨{reading.Channel}⟩ = {meas.Corr[reading.Channel].ToString("+0.000;-0.000", Inv)}" : "");
        return new BatteryCase(name, detail, expected, deviation ?? expected);
    }

    /// <summary>The F135 Law-A face: a non-integer watcher ratio's explicit in-between,
    /// checked against the measured mutual information (the classifier says Generic here;
    /// the VALUE is F135's formula).</summary>
    private BatteryCase LawACase(string name, int n, (int, int, double)[] bonds, int s, int j, double ratio)
    {
        var reading = Classify(bonds, s, j);
        double predicted = RecordParityLawClaim.LawAInformation(
            RecordParityLawClaim.RecordRadius(new[] { ratio }));
        var meas = MeasureFullState(n, bonds, s, j, new double[n], TStar);
        string expected = $"Generic I={predicted.ToString("0.000000", Inv)}";
        bool ok = reading.Family == RecordFamily.Generic && Math.Abs(meas.I - predicted) < Tol;
        string detail = $"S={s}, j={j}, watcher ratio {F(ratio)}: Law A predicts " +
                        $"I = {predicted.ToString("0.000000", Inv)}, measured {meas.I.ToString("0.000000", Inv)}";
        return new BatteryCase(name, detail, expected,
            ok ? expected : $"family={reading.Family} I={meas.I.ToString("0.000000", Inv)}");
    }

    private static string FormatReading(RecordReading r) => r.Family switch
    {
        RecordFamily.Pointer or RecordFamily.Bell or RecordFamily.RoleSwap =>
            $"{r.Family} {r.Channel} sign={(r.Sign > 0 ? "+1" : "-1")} I={r.Bits.ToString("0.000000", Inv)}",
        RecordFamily.Dark => "Dark I=0.000000",
        RecordFamily.Entangled => $"Entangled I={r.Bits.ToString("0.000000", Inv)}",
        _ => "Generic",
    };

    /// <summary>Check the measured pair against the closed-form reading; null = all checks pass.
    /// Structural (γ = 0) cases additionally pin the conditional trace distance and, for Bell,
    /// the {½, ½, 0, 0} eigenvalues; priced cases pin the κ-dressed correlator.</summary>
    private static string? VerifyReading(RecordReading r, PairMeasurement m, double gS, double gJ, bool structural)
    {
        string? Dev(string what, double got, double want) =>
            Math.Abs(got - want) < Tol ? null : $"{what}: got {got.ToString("0.000000000", Inv)}, want {want.ToString("0.000000000", Inv)}";

        switch (r.Family)
        {
            case RecordFamily.Pointer:
            {
                double kappa = Math.Exp(-2.0 * gJ * TStar);
                return Dev("I", m.I, r.Bits)
                    ?? Dev("ZY", m.Corr["ZY"], r.Sign * kappa)
                    ?? Dev("XX", m.Corr["XX"], 0.0)
                    ?? Dev("YY", m.Corr["YY"], 0.0)
                    ?? (structural ? Dev("condTD", m.CondTraceDistance, 1.0) : null);
            }
            case RecordFamily.Bell:
            {
                double kappa = Math.Exp(-2.0 * (gS + gJ) * TStar);
                string other = r.Channel == "YY" ? "XX" : "YY";
                return Dev("I", m.I, r.Bits)
                    ?? Dev(r.Channel, m.Corr[r.Channel], r.Sign * kappa)
                    ?? Dev(other, m.Corr[other], 0.0)
                    ?? Dev("ZY", m.Corr["ZY"], 0.0)
                    ?? (structural
                        ? Dev("condTD", m.CondTraceDistance, 0.0)
                          ?? Dev("eig0", m.Eigs[0], 0.5) ?? Dev("eig1", m.Eigs[1], 0.5)
                          ?? Dev("eig2", m.Eigs[2], 0.0) ?? Dev("eig3", m.Eigs[3], 0.0)
                          ?? Dev("XY", m.Corr["XY"], 0.0) ?? Dev("YX", m.Corr["YX"], 0.0)
                          ?? Dev("ZZ", m.Corr["ZZ"], 0.0)
                        : null);
            }
            case RecordFamily.RoleSwap:
            {
                double kappa = Math.Exp(-2.0 * gS * TStar);
                return Dev("I", m.I, r.Bits)
                    ?? Dev("YZ", m.Corr["YZ"], r.Sign * kappa)
                    ?? Dev("ZY", m.Corr["ZY"], 0.0);
            }
            case RecordFamily.Dark:
                return Dev("I", m.I, 0.0);
            case RecordFamily.Entangled:
                return double.IsNaN(r.Bits) ? null : Dev("I", m.I, r.Bits);
            default:
                return null;
        }
    }

    // ==================================================================
    // IInspectable
    // ==================================================================

    public string DisplayName => $"RecordLawWitness (γ={Gamma.ToString("0.###", Inv)}, t*=π/4)";

    public string Summary =>
        "the two record laws recomputed live, closed form vs full state: F135 (who records is the " +
        "parity arithmetic of the watcher ratios) and F136 (which operator is the shared-neighborhood " +
        "parity: all even + write bond → pointer, all odd → Bell with the letter alternating on the " +
        $"dresser count, mixed → dark). Battery {PassCount}/{Cases.Count} pass; every case builds the " +
        "full 2^N closed-form state and partial-traces it, never the channel formulas under test. " +
        $"Prices at γ = {Gamma.ToString("0.###", Inv)}: pointer {PointerPricedBits.ToString("0.000000", Inv)} " +
        $"(only the witness pays), Bell {BellPricedBits.ToString("0.000000", Inv)} (both pay). " +
        "Proofs: docs/proofs/PROOF_RECORD_PARITY_LAW.md + PROOF_RECORD_LETTER_LAW.md.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("the two laws",
                summary: "F135 Prop 1: every pair page is closed-form for all t (the Absorption substrate); " +
                         "Law A: watcher even → perfect, odd → blind, non-integer → 1 − h₂((1+β)/2). " +
                         "F136: private watchers must be even; the dressers' parity selects pointer/Bell; " +
                         "the letter is the dresser parity (m odd YY, m even XX); Bell pays both sites, " +
                         "pointer only the witness; a pendant S role-swaps; K_{R+1,2} fans out.");
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"),
                    provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
