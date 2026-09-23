using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Live physical-cell reading of the N=11 mixed-parity obstruction to F154's
/// conditional identity. The exact radical equalities are owned by the companion SymPy gate;
/// this independent C# reading uses a conservative floating-point summation budget.</summary>
public sealed class CompressedDensityN11Witness : IInspectable
{
    public const int SiteCount = 11;
    private const int ModeCount = SiteCount;
    private static readonly (int Ket, int Bra)[] FrequencyDyads =
        { (1, 6), (3, 7), (5, 9), (6, 11) };
    private readonly Func<int, int, int, double> _cellIndicator;
    private readonly Lazy<Snapshot> _reading;

    public CompressedDensityN11Witness() : this(PhysicalDisagreement) { }

    /// <summary>Indicator injection permits a wrong one-sided cell rule to traverse exactly
    /// the same projection path in tests. Production inspection always uses the physical rule.</summary>
    public CompressedDensityN11Witness(Func<int, int, int, double> cellIndicator)
    {
        _cellIndicator = cellIndicator ?? throw new ArgumentNullException(nameof(cellIndicator));
        _reading = new(() => Reconstruct(_cellIndicator));
    }

    public sealed record Check(string Name, double Expected, double Actual, double ErrorBudget)
    {
        public bool Passes => Math.Abs(Actual - Expected) <= ErrorBudget;
        public string Detail => $"{Name}: expected {Expected:R}; computed {Actual:R}; " +
                                $"|difference|={Math.Abs(Actual - Expected):G5}; budget={ErrorBudget:G5}";
    }

    public sealed record Snapshot(
        int FrequencyMultiplicity,
        int FrequencyMembershipMismatchCount,
        double LeftLocalCross,
        double RightLocalCross,
        double ContrastCross,
        double ContrastMatrixResidual,
        double ContrastTraceSquare,
        double SignedAgreementContrastResidual,
        double BalancedPhysicalCross,
        double BalancedConditionalPredictionCross,
        double BalancedIntervalResidual,
        double UniformConditionalResidual,
        int ZeroFrequencyMultiplicity,
        int ZeroFrequencyMembershipMismatchCount,
        double IdentityPhysicalDiagonalResidual,
        double IdentityZeroRateResidual,
        double ChiralPhysicalDiagonalResidual,
        double ChiralNormSquared,
        double ChiralLowerEndpointResidual,
        double OffLocusRayleigh,
        double OffLocusLowerBound,
        double ErrorBudget,
        IReadOnlyList<Check> Checks);

    public Snapshot Reading => _reading.Value;
    public string DisplayName => "F154 compressed density (physical N=11 live witness)";
    public string Summary
    {
        get
        {
            var r = Reading;
            bool identityFails = Math.Abs(r.BalancedPhysicalCross - r.BalancedConditionalPredictionCross) >
                                 100 * r.ErrorBudget;
            return $"N=11 (1,1): C_0 cross={r.ContrastCross:G8}; balanced F154 identity " +
                   $"{(identityFails ? "fails" : "holds at reading precision")}; " +
                   $"interval residual={r.BalancedIntervalResidual:G5}; " +
                   $"zero-frequency endpoints residual={Math.Max(r.IdentityZeroRateResidual, r.ChiralLowerEndpointResidual):G5}; " +
                   $"{r.Checks.Count(c => c.Passes)}/{r.Checks.Count} physical-cell checks.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("model and frequency space",
                "N=11 open uniform XY, h hopping 2J, ψ_k(z)=sqrt(2/12) sin(kπ(z+1)/12). " +
                "The complete frequency Ω at E_1−E_6 has dyads (1,6),(3,7),(5,9),(6,11). " +
                $"The live one-body collision count is {Reading.FrequencyMultiplicity}.");
            yield return new InspectableNode("physical operator",
                "N_l|a><b| = 1[(a=l) XOR (b=l)]|a><b|; D=-2Σ_l γ_l N_l. " +
                "Every matrix element here is summed over all 121 (1,1) cells. " +
                "Rates for the identity counterexample are (2,1,1,1,1,1,1,1,1,1,0).");
            yield return new InspectableNode("zero-frequency endpoint preparations",
                "The physical-cell matrices I=Σ_k|ψ_k><ψ_k| and " +
                "Q=|ψ_1><ψ_1|−|ψ_11><ψ_11| are reconstructed from all 11 mode projectors. " +
                "The same cell dissipator gives D I=0 and projected D Q=−4γbar Q " +
                "on the balanced profile; these reach the two block-wide (1,1) endpoints.");
            foreach (var check in Reading.Checks)
                yield return new InspectableNode(check.Name,
                    $"{check.Detail}; {(check.Passes ? "PASS" : "FAIL")}", provenance: NodeProvenance.Live);
            yield return new InspectableNode("scope",
                "The C_l=0 identity fails on this balanced N=11 profile. The separate " +
                "(1,1) interval theorem applies to complete frequency spaces; the off-locus " +
                "Rayleigh reading breaks its premise. The zero-frequency I/Q witnesses reach both " +
                "block-wide endpoints on uniform XY; no claim puts both in every frequency room. " +
                "These are strong-coupling dissipator compressions, not finite-J Liouvillian eigenvalues. " +
                "Floating-point budgets are numerical checks, not the exact symbolic certificate.");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    private static Snapshot Reconstruct(Func<int, int, int, double> indicator)
    {
        // The numerical budget covers the direct 121-cell sums and four-dyad matrix products.
        // The exact gate separately checks the radicals with equality in SymPy.
        double budget = 8192.0 * SiteCount * (Math.BitIncrement(1.0) - 1.0);
        double[,] modes = new double[ModeCount + 1, SiteCount];
        for (int k = 1; k <= ModeCount; k++)
            for (int z = 0; z < SiteCount; z++)
                modes[k, z] = Math.Sqrt(2.0 / (SiteCount + 1)) *
                              Math.Sin(k * Math.PI * (z + 1) / (SiteCount + 1));

        double[,,] dyads = new double[FrequencyDyads.Length, SiteCount, SiteCount];
        for (int d = 0; d < FrequencyDyads.Length; d++)
        {
            var (ket, bra) = FrequencyDyads[d];
            for (int a = 0; a < SiteCount; a++)
                for (int b = 0; b < SiteCount; b++)
                    dyads[d, a, b] = modes[ket, a] * modes[bra, b];
        }

        var scannedDyads = new HashSet<(int Ket, int Bra)>();
        double frequency = Energy(1) - Energy(6);
        for (int ket = 1; ket <= ModeCount; ket++)
            for (int bra = 1; bra <= ModeCount; bra++)
                if (Math.Abs(Energy(ket) - Energy(bra) - frequency) <= budget)
                    scannedDyads.Add((ket, bra));
        int multiplicity = scannedDyads.Count;
        int membershipMismatches = scannedDyads.Except(FrequencyDyads).Count() +
                                   FrequencyDyads.Except(scannedDyads).Count();

        int dim = FrequencyDyads.Length;
        double[,,] local = new double[SiteCount, dim, dim];
        for (int site = 0; site < SiteCount; site++)
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                    for (int a = 0; a < SiteCount; a++)
                        for (int b = 0; b < SiteCount; b++)
                            local[site, i, j] += dyads[i, a, b] *
                                                 indicator(a, b, site) * dyads[j, a, b];

        double[] balanced = Enumerable.Repeat(1.0, SiteCount).ToArray();
        balanced[0] = 2.0;
        balanced[SiteCount - 1] = 0.0;
        double[] uniform = Enumerable.Repeat(1.0, SiteCount).ToArray();
        double[] offLocus = new double[SiteCount];
        offLocus[0] = 1.0;

        double Physical(double[] profile, int i, int j)
        {
            double result = 0;
            for (int site = 0; site < SiteCount; site++)
                result -= 2.0 * profile[site] * local[site, i, j];
            return result;
        }

        double T(double[] profile, int i, int j)
        {
            double result = 0;
            for (int site = 0; site < SiteCount; site++)
                result += profile[site] * dyads[i, site, site] * dyads[j, site, site];
            return result;
        }

        double directSize(int i, int j)
        {
            double result = 0;
            for (int a = 0; a < SiteCount; a++)
                for (int b = 0; b < SiteCount; b++)
                    result += (a == b ? 0.0 : 2.0) * dyads[i, a, b] * dyads[j, a, b];
            return result;
        }

        double intervalResidual = 0.0;
        double uniformResidual = 0.0;
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
            {
                intervalResidual = Math.Max(intervalResidual,
                    Math.Abs(Physical(balanced, i, j) + (i == j ? 4.0 : 0.0) -
                             4.0 * T(balanced, i, j)));
                uniformResidual = Math.Max(uniformResidual,
                    Math.Abs(Physical(uniform, i, j) + 2.0 * directSize(i, j)));
            }

        double left = local[0, 2, 0];
        double right = local[SiteCount - 1, 2, 0];
        double contrast = left - right;
        int[,] crossPattern =
        {
            { 0, 1, 1, 0 },
            { 1, 0, 0, 1 },
            { 1, 0, 0, 1 },
            { 0, 1, 1, 0 },
        };
        double contrastMatrixResidual = 0.0;
        double traceSquare = 0.0;
        double signedAgreementResidual = 0.0;
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
            {
                double cell = local[0, i, j] - local[SiteCount - 1, i, j];
                contrastMatrixResidual = Math.Max(contrastMatrixResidual,
                    Math.Abs(cell + Math.Sqrt(2.0) * crossPattern[i, j] / 72.0));
                traceSquare += cell * cell;
                for (int site = 0; site < SiteCount; site++)
                {
                    int reflected = SiteCount - 1 - site;
                    double agreementContrast =
                        dyads[i, site, site] * dyads[j, site, site] -
                        dyads[i, reflected, reflected] * dyads[j, reflected, reflected];
                    signedAgreementResidual = Math.Max(signedAgreementResidual,
                        Math.Abs(local[site, i, j] - local[reflected, i, j] +
                                 2.0 * agreementContrast));
                }
            }
        double balancedCross = Physical(balanced, 2, 0);
        double conditionalCross = -2.0 * directSize(2, 0);
        double offRayleigh = Physical(offLocus, 2, 2);
        double offLowerBound = -4.0 / SiteCount;

        var scannedZeroDyads = new HashSet<(int Ket, int Bra)>();
        for (int ket = 1; ket <= ModeCount; ket++)
            for (int bra = 1; bra <= ModeCount; bra++)
                if (Math.Abs(Energy(ket) - Energy(bra)) <= budget)
                    scannedZeroDyads.Add((ket, bra));
        var expectedZeroDyads = Enumerable.Range(1, ModeCount).Select(k => (Ket: k, Bra: k)).ToArray();
        int zeroMembershipMismatches = scannedZeroDyads.Except(expectedZeroDyads).Count() +
                                       expectedZeroDyads.Except(scannedZeroDyads).Count();

        double[,] identity = new double[SiteCount, SiteCount];
        double[,] chiral = new double[SiteCount, SiteCount];
        double[] physicalRate = new double[SiteCount * SiteCount];
        double identityDiagonalResidual = 0.0;
        double chiralDiagonalResidual = 0.0;
        double chiralNormSquared = 0.0;
        for (int a = 0; a < SiteCount; a++)
            for (int b = 0; b < SiteCount; b++)
            {
                for (int k = 1; k <= ModeCount; k++)
                    identity[a, b] += modes[k, a] * modes[k, b];
                chiral[a, b] = modes[1, a] * modes[1, b] -
                                modes[ModeCount, a] * modes[ModeCount, b];
                for (int site = 0; site < SiteCount; site++)
                    physicalRate[a * SiteCount + b] -=
                        2.0 * balanced[site] * indicator(a, b, site);
                identityDiagonalResidual = Math.Max(identityDiagonalResidual,
                    Math.Abs(identity[a, b] - (a == b ? 1.0 : 0.0)));
                if (a == b)
                    chiralDiagonalResidual = Math.Max(chiralDiagonalResidual, Math.Abs(chiral[a, a]));
                chiralNormSquared += chiral[a, b] * chiral[a, b];
            }

        double identityRateResidual = 0.0;
        double chiralRateResidual = 0.0;
        for (int k = 1; k <= ModeCount; k++)
        {
            double identityAction = 0.0;
            double chiralAction = 0.0;
            for (int a = 0; a < SiteCount; a++)
                for (int b = 0; b < SiteCount; b++)
                {
                    double projector = modes[k, a] * modes[k, b];
                    double rate = physicalRate[a * SiteCount + b];
                    identityAction += projector * rate * identity[a, b];
                    chiralAction += projector * rate * chiral[a, b];
                }
            identityRateResidual = Math.Max(identityRateResidual, Math.Abs(identityAction));
            double chiralCoefficient = k == 1 ? 1.0 : k == ModeCount ? -1.0 : 0.0;
            chiralRateResidual = Math.Max(chiralRateResidual,
                Math.Abs(chiralAction + 4.0 * chiralCoefficient));
        }

        double rootTwo = Math.Sqrt(2.0);
        var checks = new List<Check>
        {
            new("complete N=11 frequency multiplicity", 4, multiplicity, 0),
            new("complete N=11 frequency dyad index set", 0, membershipMismatches, 0),
            new("left physical local cross", -rootTwo / 144.0, left, budget),
            new("right physical local cross", rootTwo / 144.0, right, budget),
            new("mirror contrast cross", -rootTwo / 72.0, contrast, budget),
            new("complete four-dyad contrast matrix", 0.0, contrastMatrixResidual, budget),
            new("contrast trace square", 1.0 / 324.0, traceSquare, budget),
            new("signed agreement contrast identity", 0.0, signedAgreementResidual, budget),
            new("balanced physical D cross", rootTwo / 36.0, balancedCross, budget),
            new("conditional F154 prediction cross", 0.0, conditionalCross, budget),
            new("balanced (1,1) interval decomposition", 0.0, intervalResidual, budget),
            new("uniform conditional size identity", 0.0, uniformResidual, budget),
            new("complete zero-frequency multiplicity", SiteCount, scannedZeroDyads.Count, 0),
            new("complete zero-frequency dyad index set", 0, zeroMembershipMismatches, 0),
            new("physical identity matrix", 0.0, identityDiagonalResidual, budget),
            new("physical identity has zero dissipative rate", 0.0, identityRateResidual, budget),
            new("chiral projector difference has zero physical diagonal", 0.0,
                chiralDiagonalResidual, budget),
            new("chiral projector difference is nonzero", 2.0, chiralNormSquared, budget),
            new("chiral projector difference reaches lower endpoint", 0.0,
                chiralRateResidual, budget),
            new("off-locus Y Rayleigh reading", -(22.0 + 5.0 * Math.Sqrt(3.0)) / 72.0,
                offRayleigh, budget),
        };
        return new Snapshot(multiplicity, membershipMismatches, left, right, contrast,
            contrastMatrixResidual, traceSquare, signedAgreementResidual,
            balancedCross, conditionalCross,
            intervalResidual, uniformResidual,
            scannedZeroDyads.Count, zeroMembershipMismatches,
            identityDiagonalResidual, identityRateResidual, chiralDiagonalResidual,
            chiralNormSquared, chiralRateResidual,
            offRayleigh, offLowerBound, budget,
            checks.AsReadOnly());
    }

    private static double Energy(int k) => 4.0 * Math.Cos(k * Math.PI / (SiteCount + 1));

    private static double PhysicalDisagreement(int ket, int bra, int site) =>
        (ket == site) == (bra == site) ? 0.0 : 1.0;
}
