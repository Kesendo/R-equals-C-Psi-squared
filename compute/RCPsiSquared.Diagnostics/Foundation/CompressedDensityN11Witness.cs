using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>An element a + b√2 + c√3 + d√6 of Q(√2, √3), with exact rational coefficients.
/// {1, √2, √3, √6} is a basis of the field over Q, so coefficient equality is equality of the
/// represented numbers, and <see cref="Sign"/> decides the sign by squaring, never by rounding.
/// Every N=11 open-chain sine mode √(2/12)·sin(kπ(z+1)/12) and every energy 4cos(kπ/12) lies in
/// this field, which is what lets <see cref="CompressedDensityN11Witness"/> compare with ==.</summary>
public readonly struct Sqrt23 : IEquatable<Sqrt23>
{
    public BigRational A { get; }
    public BigRational B { get; }
    public BigRational C { get; }
    public BigRational D { get; }

    public Sqrt23(BigRational a, BigRational b, BigRational c, BigRational d)
    {
        A = a; B = b; C = c; D = d;
    }

    public static Sqrt23 Rational(BigRational q) => new(q, BigRational.Zero, BigRational.Zero, BigRational.Zero);
    public static Sqrt23 Rational(long numerator, long denominator) =>
        Rational(new BigRational(numerator, denominator));
    public static readonly Sqrt23 Zero = Rational(0);
    public static readonly Sqrt23 One = Rational(1);
    public static readonly Sqrt23 Sqrt2 = new(0, 1, 0, 0);
    public static readonly Sqrt23 Sqrt3 = new(0, 0, 1, 0);
    public static readonly Sqrt23 Sqrt6 = new(0, 0, 0, 1);

    public bool IsZero => A.IsZero && B.IsZero && C.IsZero && D.IsZero;

    /// <summary>The exact inverse: multiply by the √2-conjugate to land in Q(√3), then by the √3-conjugate
    /// to land in Q. Throws on zero.</summary>
    public Sqrt23 Inverse()
    {
        if (IsZero) throw new DivideByZeroException("Sqrt23 zero has no inverse.");
        var conj2 = new Sqrt23(A, -B, C, -D);
        var inQ3 = this * conj2;                       // B = D = 0
        var conj3 = new Sqrt23(inQ3.A, 0, -inQ3.C, 0);
        var norm = (inQ3 * conj3).A;                    // rational, nonzero
        return (BigRational.One / norm) * (conj2 * conj3);
    }

    /// <summary>Exact rank by Gaussian elimination over the field.</summary>
    public static int Rank(Sqrt23[,] matrix)
    {
        int rows = matrix.GetLength(0), cols = matrix.GetLength(1);
        var m = (Sqrt23[,])matrix.Clone();
        int rank = 0;
        for (int col = 0; col < cols && rank < rows; col++)
        {
            int pivot = -1;
            for (int r = rank; r < rows; r++)
                if (!m[r, col].IsZero) { pivot = r; break; }
            if (pivot < 0) continue;
            for (int c = 0; c < cols; c++) (m[rank, c], m[pivot, c]) = (m[pivot, c], m[rank, c]);
            var inv = m[rank, col].Inverse();
            for (int r = 0; r < rows; r++)
            {
                if (r == rank || m[r, col].IsZero) continue;
                var factor = m[r, col] * inv;
                for (int c = col; c < cols; c++) m[r, c] = m[r, c] - factor * m[rank, c];
            }
            rank++;
        }
        return rank;
    }

    public static Sqrt23 operator +(Sqrt23 x, Sqrt23 y) => new(x.A + y.A, x.B + y.B, x.C + y.C, x.D + y.D);
    public static Sqrt23 operator -(Sqrt23 x, Sqrt23 y) => new(x.A - y.A, x.B - y.B, x.C - y.C, x.D - y.D);
    public static Sqrt23 operator -(Sqrt23 x) => new(-x.A, -x.B, -x.C, -x.D);
    public static Sqrt23 operator *(BigRational q, Sqrt23 x) => new(q * x.A, q * x.B, q * x.C, q * x.D);

    /// <summary>√2·√2 = 2, √2·√3 = √6, √2·√6 = 2√3, √3·√3 = 3, √3·√6 = 3√2, √6·√6 = 6.</summary>
    public static Sqrt23 operator *(Sqrt23 x, Sqrt23 y) => new(
        x.A * y.A + 2 * x.B * y.B + 3 * x.C * y.C + 6 * x.D * y.D,
        x.A * y.B + x.B * y.A + 3 * (x.C * y.D + x.D * y.C),
        x.A * y.C + x.C * y.A + 2 * (x.B * y.D + x.D * y.B),
        x.A * y.D + x.D * y.A + x.B * y.C + x.C * y.B);

    public static bool operator ==(Sqrt23 x, Sqrt23 y) => x.Equals(y);
    public static bool operator !=(Sqrt23 x, Sqrt23 y) => !x.Equals(y);
    public bool Equals(Sqrt23 other) => A == other.A && B == other.B && C == other.C && D == other.D;
    public override bool Equals(object? obj) => obj is Sqrt23 other && Equals(other);
    public override int GetHashCode() => HashCode.Combine(A, B, C, D);

    /// <summary>The exact sign. Write x = u + √2·v with u = A + C√3 and v = B + D√3; where u
    /// and v disagree in sign, the larger of u² and 2v² decides, and each comparison inside
    /// Q(√3) is decided the same way one level down.</summary>
    public int Sign
    {
        get
        {
            int su = SignQ3(A, C), sv = SignQ3(B, D);
            if (su == 0) return sv;
            if (sv == 0 || su == sv) return su;
            // u² − 2v² = (A² + 3C² − 2B² − 6D²) + (2AC − 4BD)√3
            int s = SignQ3(A * A + 3 * C * C - 2 * B * B - 6 * D * D, 2 * A * C - 4 * B * D);
            return su * s;
        }
    }

    private static int SignQ3(BigRational p, BigRational q)
    {
        int sp = p.Sign, sq = q.Sign;
        if (sp == 0) return sq;
        if (sq == 0 || sp == sq) return sp;
        return sp * (p * p - 3 * q * q).Sign;
    }

    /// <summary>A display value only; no decision in the witness reads it.</summary>
    public double ToDouble() =>
        Ratio(A) + Ratio(B) * Math.Sqrt(2.0) + Ratio(C) * Math.Sqrt(3.0) + Ratio(D) * Math.Sqrt(6.0);

    private static double Ratio(BigRational q) => (double)q.Numerator / (double)q.Denominator;

    public override string ToString()
    {
        var parts = new List<string>();
        if (!A.IsZero) parts.Add(A.ToString());
        if (!B.IsZero) parts.Add($"({B})√2");
        if (!C.IsZero) parts.Add($"({C})√3");
        if (!D.IsZero) parts.Add($"({D})√6");
        return parts.Count == 0 ? "0" : string.Join(" + ", parts);
    }
}

/// <summary>Live physical-cell reading of the N=11 mixed-parity obstruction to F154's
/// conditional identity and of the zero-frequency endpoint room, recomputed at inspect time in
/// exact Q(√2, √3) arithmetic (<see cref="Sqrt23"/>). Every check is an equality or an exact sign;
/// frequency membership is decided exactly as well, so no tolerance enters anywhere. The SymPy gate
/// <c>simulations/n11_compressed_density_gate.py</c> computes the same radicals independently.</summary>
public sealed class CompressedDensityN11Witness : IInspectable
{
    public const int SiteCount = 11;
    private const int ModeCount = SiteCount;
    private const int Modulus = SiteCount + 1;
    private static readonly (int Ket, int Bra)[] FrequencyDyads =
        { (1, 6), (3, 7), (5, 9), (6, 11) };
    private readonly Func<int, int, int, int> _cellIndicator;
    private readonly Lazy<Snapshot> _reading;

    public CompressedDensityN11Witness() : this(PhysicalDisagreement) { }

    /// <summary>Indicator injection permits a wrong cell rule to traverse exactly the same
    /// projection path in tests, the zero-frequency endpoint half included. Production inspection
    /// always uses the physical rule.</summary>
    public CompressedDensityN11Witness(Func<int, int, int, int> cellIndicator)
    {
        _cellIndicator = cellIndicator ?? throw new ArgumentNullException(nameof(cellIndicator));
        _reading = new(() => Reconstruct(_cellIndicator));
    }

    public sealed record Check(string Name, string Expected, string Actual, bool Passes)
    {
        public string Detail => $"{Name}: expected {Expected}; computed {Actual}";
    }

    public sealed record Snapshot(
        int FrequencyMultiplicity,
        int FrequencyMembershipMismatchCount,
        Sqrt23 LeftLocalCross,
        Sqrt23 RightLocalCross,
        Sqrt23 ContrastCross,
        Sqrt23 ContrastTraceSquare,
        Sqrt23 BalancedPhysicalCross,
        Sqrt23 BalancedConditionalPredictionCross,
        IReadOnlyList<Sqrt23> BalancedRoomCharacteristicPolynomial,
        int ZeroFrequencyMultiplicity,
        int ZeroFrequencyMembershipMismatchCount,
        Sqrt23 IdentityRayleigh,
        Sqrt23 ChiralRayleigh,
        Sqrt23 OffLocusRayleigh,
        int RoomCount,
        int RoomsReachingLowerEndpoint,
        IReadOnlyList<Check> Checks)
    {
        public bool AllPass => Checks.All(c => c.Passes);
    }

    public Snapshot Reading => _reading.Value;
    public string DisplayName => "F154 compressed density (physical N=11 live witness)";
    public string Summary
    {
        get
        {
            var r = Reading;
            bool identityFails = r.BalancedPhysicalCross != r.BalancedConditionalPredictionCross;
            return $"N=11 (1,1): C_0 cross={r.ContrastCross} ≈ {r.ContrastCross.ToDouble():G8}; " +
                   $"balanced F154 identity {(identityFails ? "fails" : "holds")}; " +
                   $"zero-frequency endpoints {r.IdentityRayleigh} and {r.ChiralRayleigh}; " +
                   $"{r.Checks.Count(c => c.Passes)}/{r.Checks.Count} exact physical-cell checks.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("model and frequency space",
                "N=11 open uniform XY, h hopping 2J, ψ_k(z)=sqrt(2/12) sin(kπ(z+1)/12). " +
                "The complete frequency Ω at E_1−E_6 has dyads (1,6),(3,7),(5,9),(6,11). " +
                $"The live one-body collision count is {Reading.FrequencyMultiplicity}, decided by exact " +
                "equality in Q(√2,√3).");
            yield return new InspectableNode("physical operator",
                "N_l|a><b| = 1[(a=l) XOR (b=l)]|a><b|; D=-2Σ_l γ_l N_l. " +
                "Every matrix element here is summed over all 121 (1,1) cells. " +
                "Rates for the identity counterexample are (2,1,1,1,1,1,1,1,1,1,0).");
            yield return new InspectableNode("zero-frequency endpoint room",
                "Every mode projector is reflection-even, so the zero-frequency room is a scalar-parity " +
                "room and F154's identity is a theorem there; its compression is −4γbar(I − G) with " +
                "F143's Gram G = (1/12)(11ᵀ + (I+R_χ)/2), R_χ the chiral mode flip k ↦ 12−k (the XY " +
                "case of PROOF_R90_FROZEN_DIVISOR Lemma 5), at every locus profile, rechecked here at the " +
                "balanced and at a signed profile. The physical-cell matrices I=Σ_k|ψ_k><ψ_k| and " +
                "Q=|ψ_1><ψ_1|−|ψ_11><ψ_11| reach its two ends, 0 and −4γbar; the lower one is F140's " +
                "frozen root at every J.");
            foreach (var check in Reading.Checks)
                yield return new InspectableNode(check.Name,
                    $"{check.Detail}; {(check.Passes ? "PASS" : "FAIL")}", provenance: NodeProvenance.Live);
            yield return new InspectableNode("scope",
                "The C_l=0 identity fails on this balanced N=11 profile. The separate " +
                "(1,1) interval theorem applies to complete frequency spaces; the off-locus " +
                "Rayleigh reading breaks its premise. The upper endpoint 0 is absent from the " +
                "nonzero-frequency room read here. Containment is a strong-coupling compression " +
                "statement; the endpoint values themselves are also exact finite-J (1,1) eigenvalues.");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    private static readonly Sqrt23[] SinTable = BuildSinTable();

    /// <summary>sin(jπ/12) for j = 0..23, exact.</summary>
    private static Sqrt23[] BuildSinTable()
    {
        var quarter = new BigRational(1, 4);
        var half = new BigRational(1, 2);
        var s = new Sqrt23[24];
        s[0] = Sqrt23.Zero;
        s[1] = quarter * (Sqrt23.Sqrt6 - Sqrt23.Sqrt2);
        s[2] = Sqrt23.Rational(1, 2);
        s[3] = half * Sqrt23.Sqrt2;
        s[4] = half * Sqrt23.Sqrt3;
        s[5] = quarter * (Sqrt23.Sqrt6 + Sqrt23.Sqrt2);
        s[6] = Sqrt23.One;
        for (int j = 7; j <= 12; j++) s[j] = s[12 - j];
        for (int j = 13; j < 24; j++) s[j] = -s[j - 12];
        return s;
    }

    private static Sqrt23 Sin(int j) => SinTable[((j % 24) + 24) % 24];

    /// <summary>E_k = 4cos(kπ/12) = 4sin((k+6)π/12).</summary>
    private static Sqrt23 Energy(int k) => new BigRational(4) * Sin(k + 6);

    private static bool Same(Sqrt23[,] x, Sqrt23[,] y)
    {
        for (int i = 0; i < x.GetLength(0); i++)
            for (int j = 0; j < x.GetLength(1); j++)
                if (x[i, j] != y[i, j]) return false;
        return true;
    }

    /// <summary>Faddeev-LeVerrier over the field: ascending coefficients of det(xI − A).</summary>
    private static Sqrt23[] CharacteristicPolynomial(Sqrt23[,] a)
    {
        int n = a.GetLength(0);
        var c = new Sqrt23[n + 1];
        c[n] = Sqrt23.One;
        var m = new Sqrt23[n, n];
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) m[i, j] = Sqrt23.Zero;
        for (int k = 1; k <= n; k++)
        {
            var next = new Sqrt23[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    var sum = i == j ? c[n - k + 1] : Sqrt23.Zero;
                    for (int l = 0; l < n; l++) sum += a[i, l] * m[l, j];
                    next[i, j] = sum;
                }
            m = next;
            var trace = Sqrt23.Zero;
            for (int i = 0; i < n; i++)
                for (int l = 0; l < n; l++) trace += a[i, l] * m[l, i];
            c[n - k] = new BigRational(-1, k) * trace;
        }
        return c;
    }

    private static Snapshot Reconstruct(Func<int, int, int, int> indicator)
    {
        var norm = new BigRational(1, 6);
        var modes = new Sqrt23[ModeCount + 1, SiteCount];
        for (int k = 1; k <= ModeCount; k++)
            for (int z = 0; z < SiteCount; z++)
                modes[k, z] = norm * (Sqrt23.Sqrt6 * Sin(k * (z + 1)));

        int dim = FrequencyDyads.Length;
        var dyads = new Sqrt23[dim, SiteCount, SiteCount];
        for (int d = 0; d < dim; d++)
        {
            var (ket, bra) = FrequencyDyads[d];
            for (int a = 0; a < SiteCount; a++)
                for (int b = 0; b < SiteCount; b++)
                    dyads[d, a, b] = modes[ket, a] * modes[bra, b];
        }

        var scannedDyads = new HashSet<(int Ket, int Bra)>();
        var scannedZeroDyads = new HashSet<(int Ket, int Bra)>();
        var frequency = Energy(1) - Energy(6);
        for (int ket = 1; ket <= ModeCount; ket++)
            for (int bra = 1; bra <= ModeCount; bra++)
            {
                var gap = Energy(ket) - Energy(bra);
                if (gap == frequency) scannedDyads.Add((ket, bra));
                if (gap.IsZero) scannedZeroDyads.Add((ket, bra));
            }
        int multiplicity = scannedDyads.Count;
        int membershipMismatches = scannedDyads.Except(FrequencyDyads).Count() +
                                   FrequencyDyads.Except(scannedDyads).Count();
        var expectedZeroDyads = Enumerable.Range(1, ModeCount).Select(k => (Ket: k, Bra: k)).ToArray();
        int zeroMembershipMismatches = scannedZeroDyads.Except(expectedZeroDyads).Count() +
                                       expectedZeroDyads.Except(scannedZeroDyads).Count();

        var local = new Sqrt23[SiteCount, dim, dim];
        for (int site = 0; site < SiteCount; site++)
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                {
                    var sum = Sqrt23.Zero;
                    for (int a = 0; a < SiteCount; a++)
                        for (int b = 0; b < SiteCount; b++)
                            if (indicator(a, b, site) != 0)
                                sum += new BigRational(indicator(a, b, site)) * (dyads[i, a, b] * dyads[j, a, b]);
                    local[site, i, j] = sum;
                }

        BigRational[] balanced = Enumerable.Repeat(BigRational.One, SiteCount).ToArray();
        balanced[0] = 2;
        balanced[SiteCount - 1] = 0;
        BigRational[] uniform = Enumerable.Repeat(BigRational.One, SiteCount).ToArray();
        BigRational[] signed = { 5, -2, 1, 1, 1, 1, 1, 1, 1, 4, -3 };
        BigRational[] offLocus = Enumerable.Repeat(BigRational.Zero, SiteCount).ToArray();
        offLocus[0] = 1;

        Sqrt23[,] Physical(BigRational[] profile)
        {
            var d = new Sqrt23[dim, dim];
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                {
                    var sum = Sqrt23.Zero;
                    for (int site = 0; site < SiteCount; site++)
                        sum += (-2 * profile[site]) * local[site, i, j];
                    d[i, j] = sum;
                }
            return d;
        }

        Sqrt23 T(BigRational[] profile, int i, int j)
        {
            var sum = Sqrt23.Zero;
            for (int site = 0; site < SiteCount; site++)
                sum += profile[site] * (dyads[i, site, site] * dyads[j, site, site]);
            return sum;
        }

        Sqrt23 DirectSize(int i, int j)
        {
            var sum = Sqrt23.Zero;
            for (int a = 0; a < SiteCount; a++)
                for (int b = 0; b < SiteCount; b++)
                    if (a != b) sum += new BigRational(2) * (dyads[i, a, b] * dyads[j, a, b]);
            return sum;
        }

        var balancedD = Physical(balanced);
        var uniformD = Physical(uniform);
        bool intervalHolds = true, uniformHolds = true;
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
            {
                var interval = (i == j ? Sqrt23.Rational(-4) : Sqrt23.Zero) + new BigRational(4) * T(balanced, i, j);
                intervalHolds &= balancedD[i, j] == interval;
                uniformHolds &= uniformD[i, j] == new BigRational(-2) * DirectSize(i, j);
            }

        var left = local[0, 2, 0];
        var right = local[SiteCount - 1, 2, 0];
        var contrast = left - right;
        int[,] crossPattern =
        {
            { 0, 1, 1, 0 },
            { 1, 0, 0, 1 },
            { 1, 0, 0, 1 },
            { 0, 1, 1, 0 },
        };
        var contrastEntry = new BigRational(-1, 72) * Sqrt23.Sqrt2;
        bool contrastMatrixHolds = true, signedAgreementHolds = true;
        var traceSquare = Sqrt23.Zero;
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
            {
                var cell = local[0, i, j] - local[SiteCount - 1, i, j];
                contrastMatrixHolds &= cell == (crossPattern[i, j] == 1 ? contrastEntry : Sqrt23.Zero);
                traceSquare += cell * cell;
                for (int site = 0; site < SiteCount; site++)
                {
                    int reflected = SiteCount - 1 - site;
                    var agreementContrast =
                        dyads[i, site, site] * dyads[j, site, site] -
                        dyads[i, reflected, reflected] * dyads[j, reflected, reflected];
                    signedAgreementHolds &= local[site, i, j] - local[reflected, i, j] ==
                                            new BigRational(-2) * agreementContrast;
                }
            }
        var balancedCross = balancedD[2, 0];
        var conditionalCross = new BigRational(-2) * DirectSize(2, 0);
        var charpoly = CharacteristicPolynomial(balancedD);
        Sqrt23[] expectedCharpoly =
        {
            Sqrt23.Rational(14392, 81), Sqrt23.Rational(15836, 81), Sqrt23.Rational(13031, 162),
            Sqrt23.Rational(44, 3), Sqrt23.One,
        };
        bool charpolyHolds = charpoly.SequenceEqual(expectedCharpoly);
        var offRayleigh = Physical(offLocus)[2, 2];
        var offExpected = new BigRational(-1, 72) * (Sqrt23.Rational(22) + new BigRational(5) * Sqrt23.Sqrt3);
        int offBelowSign = (offRayleigh + Sqrt23.Rational(4, 11)).Sign;

        // The zero-frequency room, spanned by the mode projectors P_k.
        var projector = new Sqrt23[ModeCount + 1, SiteCount, SiteCount];
        for (int k = 1; k <= ModeCount; k++)
            for (int a = 0; a < SiteCount; a++)
                for (int b = 0; b < SiteCount; b++)
                    projector[k, a, b] = modes[k, a] * modes[k, b];

        Sqrt23[,] CellRate(BigRational[] profile)
        {
            var rate = new Sqrt23[SiteCount, SiteCount];
            for (int a = 0; a < SiteCount; a++)
                for (int b = 0; b < SiteCount; b++)
                {
                    BigRational sum = 0;
                    for (int site = 0; site < SiteCount; site++)
                        sum += -2 * profile[site] * indicator(a, b, site);
                    rate[a, b] = Sqrt23.Rational(sum);
                }
            return rate;
        }

        Sqrt23[,] ZeroRoom(BigRational[] profile)
        {
            var rate = CellRate(profile);
            var d = new Sqrt23[ModeCount, ModeCount];
            for (int j = 1; j <= ModeCount; j++)
                for (int k = j; k <= ModeCount; k++)
                {
                    var sum = Sqrt23.Zero;
                    for (int a = 0; a < SiteCount; a++)
                        for (int b = 0; b < SiteCount; b++)
                            if (!rate[a, b].IsZero)
                                sum += rate[a, b] * (projector[j, a, b] * projector[k, a, b]);
                    d[j - 1, k - 1] = sum;
                    d[k - 1, j - 1] = sum;
                }
            return d;
        }

        // F143: G = (1/M)(11ᵀ + (I + R_χ)/2), R_χ the chiral involution k ↦ M − k; the zero-frequency
        // compression on the locus is −4γbar(I − G), and every profile here has γbar = 1.
        var f143Room = new Sqrt23[ModeCount, ModeCount];
        for (int j = 1; j <= ModeCount; j++)
            for (int k = 1; k <= ModeCount; k++)
            {
                var g = new BigRational(1, Modulus) *
                        (BigRational.One + new BigRational((j == k ? 1 : 0) + (j + k == Modulus ? 1 : 0), 2));
                f143Room[j - 1, k - 1] = Sqrt23.Rational(-4 * ((j == k ? BigRational.One : BigRational.Zero) - g));
            }
        var balancedZeroRoom = ZeroRoom(balanced);
        bool f143Balanced = Same(balancedZeroRoom, f143Room);
        bool f143Signed = Same(ZeroRoom(signed), f143Room);

        bool identityIsPhysicalDiagonal = true, chiralHasZeroDiagonal = true;
        var chiralNormSquared = Sqrt23.Zero;
        for (int a = 0; a < SiteCount; a++)
            for (int b = 0; b < SiteCount; b++)
            {
                var identity = Sqrt23.Zero;
                for (int k = 1; k <= ModeCount; k++) identity += projector[k, a, b];
                identityIsPhysicalDiagonal &= identity == (a == b ? Sqrt23.One : Sqrt23.Zero);
                var chiral = projector[1, a, b] - projector[ModeCount, a, b];
                if (a == b) chiralHasZeroDiagonal &= chiral.IsZero;
                chiralNormSquared += chiral * chiral;
            }

        // I and Q in the P_k coordinates; each Rayleigh value is read off the computed action,
        // and the eigenvector property is checked on every coordinate.
        var identityCoords = Enumerable.Repeat(Sqrt23.One, ModeCount).ToArray();
        var chiralCoords = new Sqrt23[ModeCount];
        for (int k = 0; k < ModeCount; k++) chiralCoords[k] = Sqrt23.Zero;
        chiralCoords[0] = Sqrt23.One;
        chiralCoords[ModeCount - 1] = -Sqrt23.One;

        (Sqrt23 Value, bool Eigen) Rayleigh(Sqrt23[] v, int pivot)
        {
            var action = new Sqrt23[ModeCount];
            for (int j = 0; j < ModeCount; j++)
            {
                var sum = Sqrt23.Zero;
                for (int k = 0; k < ModeCount; k++) sum += balancedZeroRoom[j, k] * v[k];
                action[j] = sum;
            }
            // v[pivot] is ±1, so the ratio at the pivot is the action times that sign.
            var value = v[pivot] * action[pivot];
            bool eigen = true;
            for (int j = 0; j < ModeCount; j++) eigen &= action[j] == value * v[j];
            return (value, eigen);
        }

        var (identityValue, identityEigen) = Rayleigh(identityCoords, 0);
        var (chiralValue, chiralEigen) = Rayleigh(chiralCoords, 0);

        // Every complete (1,1) frequency room at the balanced profile, exactly: the multiplicity of −4 against
        // the non-fixed orbits of the chiral transpose (a,b) ↦ (M−b, M−a), and where 0 occurs.
        var balancedRate = CellRate(balanced);
        var roomsByFrequency = new Dictionary<Sqrt23, List<(int Ket, int Bra)>>();
        for (int ket = 1; ket <= ModeCount; ket++)
            for (int bra = 1; bra <= ModeCount; bra++)
            {
                var key = Energy(ket) - Energy(bra);
                if (!roomsByFrequency.TryGetValue(key, out var list)) roomsByFrequency[key] = list = new();
                list.Add((ket, bra));
            }
        var histogram = roomsByFrequency.Values.GroupBy(r => r.Count)
            .ToDictionary(g => g.Key, g => g.Count());
        bool histogramHolds = roomsByFrequency.Count == 55 && histogram.Count == 4 &&
                              histogram.GetValueOrDefault(11) == 1 && histogram.GetValueOrDefault(4) == 6 &&
                              histogram.GetValueOrDefault(2) == 38 && histogram.GetValueOrDefault(1) == 10;
        int orbitMismatches = 0, roomsAtLower = 0, roomsAtZero = 0;
        bool zeroOnlyAtZeroFrequency = true;
        foreach (var (key, room) in roomsByFrequency)
        {
            int n = room.Count;
            var d = new Sqrt23[n, n];
            for (int i = 0; i < n; i++)
                for (int j = i; j < n; j++)
                {
                    var sum = Sqrt23.Zero;
                    for (int a = 0; a < SiteCount; a++)
                        for (int b = 0; b < SiteCount; b++)
                            if (!balancedRate[a, b].IsZero)
                                sum += balancedRate[a, b] *
                                       (modes[room[i].Ket, a] * modes[room[i].Bra, b] *
                                        (modes[room[j].Ket, a] * modes[room[j].Bra, b]));
                    d[i, j] = sum;
                    d[j, i] = sum;
                }
            var shifted = (Sqrt23[,])d.Clone();
            for (int i = 0; i < n; i++) shifted[i, i] = shifted[i, i] + Sqrt23.Rational(4);
            int lowerMultiplicity = n - Sqrt23.Rank(shifted);
            int zeroMultiplicity = n - Sqrt23.Rank(d);
            int orbits = room.Where(x => x != (Modulus - x.Bra, Modulus - x.Ket))
                .Select(x => x.Ket * 100 + x.Bra < (Modulus - x.Bra) * 100 + (Modulus - x.Ket)
                    ? (x.Ket, x.Bra) : (Modulus - x.Bra, Modulus - x.Ket))
                .Distinct().Count();
            orbitMismatches += lowerMultiplicity != orbits ? 1 : 0;
            roomsAtLower += lowerMultiplicity > 0 ? 1 : 0;
            roomsAtZero += zeroMultiplicity > 0 ? 1 : 0;
            if (zeroMultiplicity > 0 && !key.IsZero) zeroOnlyAtZeroFrequency = false;
        }

        static string Show(Sqrt23 x) => $"{x} ≈ {x.ToDouble():G10}";
        static string Verdict(bool holds) => holds ? "exact equality" : "differs";
        var checks = new List<Check>
        {
            new("complete N=11 frequency multiplicity", "4", multiplicity.ToString(), multiplicity == 4),
            new("complete N=11 frequency dyad index set", "0 mismatches", membershipMismatches.ToString(),
                membershipMismatches == 0),
            new("left physical local cross", Show(new BigRational(-1, 144) * Sqrt23.Sqrt2), Show(left),
                left == new BigRational(-1, 144) * Sqrt23.Sqrt2),
            new("right physical local cross", Show(new BigRational(1, 144) * Sqrt23.Sqrt2), Show(right),
                right == new BigRational(1, 144) * Sqrt23.Sqrt2),
            new("mirror contrast cross", Show(contrastEntry), Show(contrast), contrast == contrastEntry),
            new("complete four-dyad contrast matrix", "−√2/72 · pattern", Verdict(contrastMatrixHolds),
                contrastMatrixHolds),
            new("contrast trace square", "1/324", Show(traceSquare), traceSquare == Sqrt23.Rational(1, 324)),
            new("signed agreement contrast identity, every site", "C_l = −2P(P_ll − P_mm)P",
                Verdict(signedAgreementHolds), signedAgreementHolds),
            new("balanced physical D cross", Show(new BigRational(1, 36) * Sqrt23.Sqrt2), Show(balancedCross),
                balancedCross == new BigRational(1, 36) * Sqrt23.Sqrt2),
            new("conditional F154 prediction cross", "0", Show(conditionalCross), conditionalCross.IsZero),
            new("balanced (1,1) interval decomposition", "DΩ = −4I + 4T", Verdict(intervalHolds), intervalHolds),
            new("uniform conditional size identity", "DΩ = −2 N_XY", Verdict(uniformHolds), uniformHolds),
            new("balanced room characteristic polynomial", "(x+4)²(x² + 20x/3 + 1799/162)",
                string.Join(", ", charpoly.Select(c => c.ToString())), charpolyHolds),
            new("complete zero-frequency multiplicity", SiteCount.ToString(), scannedZeroDyads.Count.ToString(),
                scannedZeroDyads.Count == SiteCount),
            new("complete zero-frequency dyad index set", "0 mismatches", zeroMembershipMismatches.ToString(),
                zeroMembershipMismatches == 0),
            new("zero-frequency room is F143's −4γbar(I − G), balanced profile", "exact equality",
                Verdict(f143Balanced), f143Balanced),
            new("zero-frequency room is F143's −4γbar(I − G), signed profile (5,−2,1,…,1,4,−3)",
                "exact equality", Verdict(f143Signed), f143Signed),
            new("physical identity matrix", "Σ_k P_k = I on cells", Verdict(identityIsPhysicalDiagonal),
                identityIsPhysicalDiagonal),
            new("identity is a zero-frequency eigenvector", "D I ∝ I", Verdict(identityEigen), identityEigen),
            new("identity reaches the upper endpoint", "0", Show(identityValue), identityValue.IsZero),
            new("chiral projector difference has zero physical diagonal", "0", Verdict(chiralHasZeroDiagonal),
                chiralHasZeroDiagonal),
            new("chiral projector difference is nonzero", "2", Show(chiralNormSquared),
                chiralNormSquared == Sqrt23.Rational(2)),
            new("chiral difference is a zero-frequency eigenvector", "D Q ∝ Q", Verdict(chiralEigen), chiralEigen),
            new("chiral difference reaches the lower endpoint", "−4", Show(chiralValue),
                chiralValue == Sqrt23.Rational(-4)),
            new("off-locus Y Rayleigh reading", Show(offExpected), Show(offRayleigh), offRayleigh == offExpected),
            new("off-locus Rayleigh lies below −4/11", "sign −1", offBelowSign.ToString(), offBelowSign < 0),
            new("(1,1) block: 55 complete frequency rooms, dimensions {11 ×1, 4 ×6, 2 ×38, 1 ×10}",
                "55 rooms", $"{roomsByFrequency.Count} rooms; " +
                string.Join(", ", histogram.OrderByDescending(h => h.Key).Select(h => $"{h.Key} ×{h.Value}")),
                histogramHolds),
            new("every room, balanced profile: mult(−4) = non-fixed chiral-transpose orbits (exact rank)",
                "0 mismatches", orbitMismatches.ToString(), orbitMismatches == 0),
            new("rooms reaching the lower endpoint −4, balanced profile", "45 of 55", $"{roomsAtLower} of " +
                $"{roomsByFrequency.Count}", roomsAtLower == 45),
            new("0 lies in the zero-frequency room only, balanced profile", "1 room, at ω = 0",
                $"{roomsAtZero} room(s); only at ω = 0: {zeroOnlyAtZeroFrequency}",
                roomsAtZero == 1 && zeroOnlyAtZeroFrequency),
        };
        return new Snapshot(multiplicity, membershipMismatches, left, right, contrast, traceSquare,
            balancedCross, conditionalCross, charpoly, scannedZeroDyads.Count, zeroMembershipMismatches,
            identityValue, chiralValue, offRayleigh, roomsByFrequency.Count, roomsAtLower, checks.AsReadOnly());
    }

    private static int PhysicalDisagreement(int ket, int bra, int site) =>
        (ket == site) == (bra == site) ? 0 : 1;
}
