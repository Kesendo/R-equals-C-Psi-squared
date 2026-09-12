using System.Globalization;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Live finite-precision companion to the fixed-N theorem in
/// <c>docs/proofs/PROOF_MISSING_PHASE_RELAXATION_SCALE.md</c>. It independently diagonalises
/// the 49-dimensional A/(1,1) generator and its adjoint, matches the two rank-two clusters on
/// the theorem branches rooted at <c>+/-2 sqrt(2)i</c>, and reads their orthogonal-projector centre light. A separate
/// seven-dimensional Hermitian route continues the three blind Hamiltonian modes and composes
/// their centre weights into the five Kato coefficients that do not use the outer-complement
/// direction. This is a live reading of the proved local germ, not a replacement proof.</summary>
public sealed class MissingPhaseRelaxationScaleWitness : IInspectable
{
    public const int SiteCount = 7;
    public const int WatchedSeat = 3;
    public const int GeneratorDimension = SiteCount * SiteCount;
    public const int ClusterRank = 2;
    public const double DefaultEpsilon = 0.01;
    public const double DefaultGamma = 0.3;
    public const double MaxAbsEpsilon = 0.1;
    public const double MachineEpsilon = 2.2204460492503131e-16;
    public static readonly double MinimumResolvableAbsEpsilon =
        2.0 * Math.Sqrt(MachineEpsilon * SiteCount * SiteCount);

    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private readonly Lazy<Snapshot> _snapshot;

    public double Epsilon { get; }
    public double Gamma { get; }

    public MissingPhaseRelaxationScaleWitness(
        double epsilon = DefaultEpsilon,
        double gamma = DefaultGamma)
    {
        ValidateEpsilon(epsilon);
        ValidateGamma(gamma);
        ValidateResolvedSpectralScale(epsilon, gamma);
        Epsilon = epsilon;
        Gamma = gamma;
        _snapshot = new Lazy<Snapshot>(() => BuildSnapshot(epsilon, gamma));
    }

    public IReadOnlyList<ClusterReading> Clusters => _snapshot.Value.Clusters;
    public IReadOnlyList<HilbertModeReading> HilbertModes => _snapshot.Value.HilbertModes;
    public IReadOnlyList<DyadCoefficientReading> DyadCoefficients => _snapshot.Value.Dyads;
    public double GeneratorFrobeniusNorm => _snapshot.Value.GeneratorFrobeniusNorm;

    /// <summary>Roundoff scale for an orthogonal 49 by 49 projector: eps_machine times the
    /// number of scalar entries. It is a dimension law, not a fitted acceptance threshold.</summary>
    public double ProjectorErrorModel => MachineEpsilon * GeneratorDimension * GeneratorDimension;

    /// <summary>Backward-error scale for the generator and absorption reads: projector scale
    /// multiplied by max(1, ||L_A||_F). No observed residual is used to choose it.</summary>
    public double SpectralErrorModel =>
        ProjectorErrorModel * Math.Max(1.0, GeneratorFrobeniusNorm);

    public double HilbertErrorModel => MachineEpsilon * SiteCount * SiteCount;

    public string ConstructionProvenance =>
        "The 49 by 49 A/(1,1) matrix is assembled cell by cell in row-major convention from the " +
        "7 by 7 end-detuned hopping and centre-only dephasing. Evd(A) and Evd(A†) are run independently; " +
        "complete repeated-eigenvalue groups are matched to dyads of the separately continued Hilbert " +
        "blind modes, and each selected two-column invariant space is orthonormalised by its own thin QR " +
        "before its projector and light are read. The Hilbert control is a separate real-symmetric 7 by 7 EVD.";

    public string Scope =>
        "Rigid scope: fixed N = 7, the A/(1,1) generator, watched seat c = 3, " +
        "h_01 = 2(1 + epsilon), every other path hopping 2, centre-only fixed gamma > 0 (finite), and " +
        $"{MinimumResolvableAbsEpsilon.ToString("E2", Inv)} < |epsilon| <= " +
        $"{MaxAbsEpsilon.ToString("0.0", Inv)} in the fixed-gamma punctured local range. " +
        "The lower guard keeps the epsilon^2 light above the 7 by 7 roundoff model and the upper cap is " +
        "a live-witness policy guard; the combined gamma epsilon^2(1+epsilon)/2 cubic prediction must also exceed the 49D " +
        "spectral error model, and the selected decay/light must remain resolved or the live read fails closed. " +
        "This is not a punctured-neighbourhood radius proved uniformly in gamma. " +
        "This is not the short-time onset witness, not an all-N theorem, not the full 4^7 Liouvillian gap, " +
        "and not an observable lifetime.";

    public string FullSpaceOnlyFence =>
        "The 3 gamma/2 zero-cluster direction uses I_E, the four-dimensional outer-complement identity. " +
        "It is explicitly full-49D-only and is not reconstructed from the three Hilbert centre weights.";

    public sealed record SubspaceReading(
        string Source,
        int Rank,
        double ProjectorTrace,
        double SeedOverlap,
        double MeanDecay,
        double CentreLight,
        double AbsorptionResidual,
        double HermitianResidual,
        double IdempotenceResidual,
        double RitzResidual);

    public sealed record ClusterReading(
        int FrequencySign,
        SubspaceReading Right,
        SubspaceReading Adjoint,
        double IndependentProjectorDistance);

    public sealed record HilbertModeReading(
        string Label,
        double Energy,
        double CentreWeight,
        double SeedOverlap);

    public sealed record DyadCoefficientReading(
        string Label,
        double MeasuredCoefficientOverGamma,
        double ExpectedCoefficientOverGamma,
        double CubicPredictionOverGamma,
        double FourthOrderScaleOverGamma,
        double MeasuredRateCoefficient,
        double ExpectedRateCoefficient,
        bool IncludesOuterComplementIdentity);

    private sealed record Snapshot(
        IReadOnlyList<ClusterReading> Clusters,
        IReadOnlyList<HilbertModeReading> HilbertModes,
        IReadOnlyList<DyadCoefficientReading> Dyads,
        double GeneratorFrobeniusNorm);

    private sealed record SubspaceWork(SubspaceReading Reading, Matrix<Complex> Projector);

    private static Snapshot BuildSnapshot(double epsilon, double gamma)
    {
        Matrix<Complex> a = BuildAGenerator(epsilon, gamma, WatchedSeat);
        Matrix<Complex> aAdjoint = a.ConjugateTranspose();
        Matrix<Complex> deltaCentre = CentreCellProjector(WatchedSeat);
        double generatorNorm = a.FrobeniusNorm();
        double spectralErrorModel = MachineEpsilon
            * GeneratorDimension
            * GeneratorDimension
            * Math.Max(1.0, generatorNorm);

        // Two genuinely separate decompositions. In particular no A† vector or projector is
        // obtained by conjugating, inverting, or copying an Evd(A) result. The seed is instead
        // built from the three separately continued blind Hilbert modes at this epsilon; asking a
        // generic EVD for a basis of the exactly doubled 49D eigenvalue would make the match depend on
        // an arbitrary, and potentially rank-deficient, degeneracy basis.
        var evdA = a.Evd();
        var evdAdjoint = aAdjoint.Evd();

        var clusters = new List<ClusterReading>(2);
        foreach (int sign in new[] { -1, +1 })
        {
            Matrix<Complex> rightSeed = BlindDyadSeed(epsilon, sign);
            Matrix<Complex> adjointSeed = BlindDyadSeed(epsilon, sign);
            SubspaceWork right = ContinueSubspace(
                a, evdA, rightSeed, deltaCentre, gamma, "Evd(A)");
            SubspaceWork adjoint = ContinueSubspace(
                aAdjoint, evdAdjoint, adjointSeed, deltaCentre, gamma, "Evd(A†)");
            clusters.Add(new ClusterReading(
                sign,
                right.Reading,
                adjoint.Reading,
                (right.Projector - adjoint.Projector).FrobeniusNorm()));
        }

        SubspaceReading? unresolved = clusters
            .SelectMany(cluster => new[] { cluster.Right, cluster.Adjoint })
            .FirstOrDefault(reading => reading.MeanDecay <= spectralErrorModel
                                    || 2.0 * gamma * reading.CentreLight <= spectralErrorModel);
        if (unresolved is not null)
            throw new InvalidOperationException(
                $"{unresolved.Source} selected-cluster spectral rate/light is unresolved: " +
                $"decay={unresolved.MeanDecay:E3}, 2 gamma w={2.0 * gamma * unresolved.CentreLight:E3}, " +
                $"spectral error model={spectralErrorModel:E3}. The live witness fails closed in this regime.");

        IReadOnlyList<HilbertModeReading> hilbertModes = ContinueBlindHilbertModes(epsilon);
        IReadOnlyList<DyadCoefficientReading> dyads = ComposeDyadCoefficients(epsilon, gamma, hilbertModes);
        return new Snapshot(clusters, hilbertModes, dyads, generatorNorm);
    }

    private static Matrix<Complex> BlindDyadSeed(double epsilon, int frequencySign)
    {
        IReadOnlyList<MathNet.Numerics.LinearAlgebra.Vector<double>> blind =
            ContinuedBlindHilbertVectors(epsilon);
        var dMinus = blind[0];
        var dZero = blind[1];
        var dPlus = blind[2];
        var pairs = frequencySign < 0
            ? new[] { (Left: dZero, Right: dMinus), (Left: dPlus, Right: dZero) }
            : new[] { (Left: dMinus, Right: dZero), (Left: dZero, Right: dPlus) };
        Matrix<Complex> dyads = Matrix<Complex>.Build.Dense(
            GeneratorDimension,
            ClusterRank,
            (coordinate, column) =>
            {
                int row = coordinate / SiteCount;
                int col = coordinate % SiteCount;
                return pairs[column].Left[row] * pairs[column].Right[col];
            });
        return dyads.QR(QRMethod.Thin).Q;
    }

    private static SubspaceWork ContinueSubspace(
        Matrix<Complex> generator,
        MathNet.Numerics.LinearAlgebra.Factorization.Evd<Complex> evd,
        Matrix<Complex> seedBasis,
        Matrix<Complex> deltaCentre,
        double gamma,
        string source)
    {
        Matrix<Complex> seedProjector = seedBasis * seedBasis.ConjugateTranspose();
        // Rank is a property of the repeated eigenvalue, not of the two individually most
        // seed-like Schur/EVD columns. Selecting columns one by one can splice together two
        // nearby doubled branches at large gamma. Form complete numerical eigenvalue groups
        // first, retain only rank-two groups, then use the seed projector to choose one group.
        double eigenvalueGroupingScale = MachineEpsilon
            * Math.Max(1.0, generator.FrobeniusNorm())
            * GeneratorDimension
            * GeneratorDimension;
        IReadOnlyList<int[]> rankTwoGroups = CompleteEigenvalueGroups(evd.EigenValues, eigenvalueGroupingScale)
            .Where(group => group.Length == ClusterRank)
            .ToArray();
        if (rankTwoGroups.Count == 0)
            throw new InvalidOperationException(
                $"{source} exposes no complete rank-two eigenvalue cluster under the numerical error law");

        int[] selected = rankTwoGroups
            .Select(group =>
            {
                Matrix<Complex> candidateBasis = OrthonormalBasis(evd.EigenVectors, group);
                double overlap = (candidateBasis.ConjugateTranspose()
                                  * seedProjector
                                  * candidateBasis).Trace().Real;
                return (group, overlap);
            })
            .OrderByDescending(candidate => candidate.overlap)
            .First().group;

        Matrix<Complex> basis = OrthonormalBasis(evd.EigenVectors, selected);
        Matrix<Complex> projector = basis * basis.ConjugateTranspose();
        Matrix<Complex> ritz = basis.ConjugateTranspose() * generator * basis;
        Complex meanEigenvalue = ritz.Trace() / ClusterRank;
        double meanDecay = -meanEigenvalue.Real;
        double light = (projector * deltaCentre).Trace().Real / ClusterRank;
        double seedOverlap = (basis.ConjugateTranspose() * seedProjector * basis).Trace().Real;
        var reading = new SubspaceReading(
            source,
            basis.ColumnCount,
            projector.Trace().Real,
            seedOverlap,
            meanDecay,
            light,
            Math.Abs(meanDecay - 2.0 * gamma * light),
            (projector - projector.ConjugateTranspose()).FrobeniusNorm(),
            (projector * projector - projector).FrobeniusNorm(),
            (generator * basis - basis * ritz).FrobeniusNorm());
        return new SubspaceWork(reading, projector);
    }

    private static Matrix<Complex> OrthonormalBasis(Matrix<Complex> eigenvectors, int[] selected)
    {
        Matrix<Complex> candidates = Matrix<Complex>.Build.Dense(
            eigenvectors.RowCount,
            selected.Length,
            (row, column) => eigenvectors[row, selected[column]]);
        Matrix<Complex> q = candidates.QR(QRMethod.Thin).Q;
        return q.ColumnCount == selected.Length
            ? q
            : q.SubMatrix(0, q.RowCount, 0, selected.Length);
    }

    private static IReadOnlyList<int[]> CompleteEigenvalueGroups(
        MathNet.Numerics.LinearAlgebra.Vector<Complex> eigenvalues,
        double tolerance)
    {
        var remaining = new HashSet<int>(Enumerable.Range(0, eigenvalues.Count));
        var groups = new List<int[]>();
        while (remaining.Count > 0)
        {
            int seed = remaining.First();
            remaining.Remove(seed);
            var group = new List<int> { seed };
            var frontier = new Queue<int>();
            frontier.Enqueue(seed);
            while (frontier.Count > 0)
            {
                int current = frontier.Dequeue();
                int[] neighbours = remaining
                    .Where(index => Complex.Abs(eigenvalues[index] - eigenvalues[current]) <= tolerance)
                    .ToArray();
                foreach (int neighbour in neighbours)
                {
                    remaining.Remove(neighbour);
                    group.Add(neighbour);
                    frontier.Enqueue(neighbour);
                }
            }
            groups.Add(group.ToArray());
        }
        return groups;
    }

    private static IReadOnlyList<HilbertModeReading> ContinueBlindHilbertModes(double epsilon)
    {
        Matrix<double> h = BuildHopping(epsilon);
        IReadOnlyList<MathNet.Numerics.LinearAlgebra.Vector<double>> vectors =
            ContinuedBlindHilbertVectors(epsilon);
        IReadOnlyList<MathNet.Numerics.LinearAlgebra.Vector<double>> seeds = BlindHilbertSeedVectors();
        string[] labels = { "d-", "d0", "d+" };
        var readings = new List<HilbertModeReading>(3);
        for (int mode = 0; mode < vectors.Count; mode++)
        {
            MathNet.Numerics.LinearAlgebra.Vector<double> vector = vectors[mode];
            MathNet.Numerics.LinearAlgebra.Vector<double> seed = seeds[mode];
            readings.Add(new HilbertModeReading(
                labels[mode],
                vector.DotProduct(h * vector),
                vector[WatchedSeat] * vector[WatchedSeat],
                Math.Abs(seed.DotProduct(vector))));
        }
        return readings;
    }

    private static IReadOnlyList<MathNet.Numerics.LinearAlgebra.Vector<double>> BlindHilbertSeedVectors()
    {
        var seedEvd = BuildHopping(0.0).Evd(Symmetricity.Symmetric);
        return Enumerable.Range(0, SiteCount)
            .OrderBy(index => Math.Abs(seedEvd.EigenVectors[WatchedSeat, index]))
            .Take(3)
            .OrderBy(index => seedEvd.EigenValues[index].Real)
            .Select(seedEvd.EigenVectors.Column)
            .ToArray();
    }

    private static IReadOnlyList<MathNet.Numerics.LinearAlgebra.Vector<double>> ContinuedBlindHilbertVectors(
        double epsilon)
    {
        IReadOnlyList<MathNet.Numerics.LinearAlgebra.Vector<double>> seeds = BlindHilbertSeedVectors();
        var currentEvd = BuildHopping(epsilon).Evd(Symmetricity.Symmetric);
        var unused = new HashSet<int>(Enumerable.Range(0, SiteCount));
        var vectors = new List<MathNet.Numerics.LinearAlgebra.Vector<double>>(3);
        foreach (MathNet.Numerics.LinearAlgebra.Vector<double> seed in seeds)
        {
            int selected = unused
                .Select(index => (index, overlap: Math.Abs(seed.DotProduct(currentEvd.EigenVectors.Column(index)))))
                .OrderByDescending(item => item.overlap)
                .First().index;
            unused.Remove(selected);
            vectors.Add(currentEvd.EigenVectors.Column(selected));
        }
        return vectors;
    }

    private static IReadOnlyList<DyadCoefficientReading> ComposeDyadCoefficients(
        double epsilon,
        double gamma,
        IReadOnlyList<HilbertModeReading> modes)
    {
        var byLabel = modes.ToDictionary(mode => mode.Label);
        double scale = epsilon * epsilon;
        DyadCoefficientReading Read(string label, string left, string right, double expected)
        {
            double measured = 2.0 * (byLabel[left].CentreWeight + byLabel[right].CentreWeight) / scale;
            return new DyadCoefficientReading(
                label,
                measured,
                expected,
                expected * (1.0 + epsilon),
                expected * epsilon * epsilon,
                gamma * measured,
                gamma * expected,
                IncludesOuterComplementIdentity: false);
        }

        return new[]
        {
            Read("-2 sqrt(2)i: |d0><d-|", "d0", "d-", 0.5),
            Read("+2 sqrt(2)i: |d-><d0|", "d-", "d0", 0.5),
            Read("-4 sqrt(2)i: |d+><d-|", "d+", "d-", 1.0),
            Read("+4 sqrt(2)i: |d-><d+|", "d-", "d+", 1.0),
            Read("0: |d+><d+| - |d-><d-|", "d-", "d-", 1.0),
        };
    }

    internal static int RowMajorIndex(int row, int column) => row * SiteCount + column;

    internal static Matrix<double> BuildHopping(double epsilon)
    {
        var h = Matrix<double>.Build.Dense(SiteCount, SiteCount);
        for (int site = 0; site < SiteCount - 1; site++)
        {
            double hopping = site == 0 ? 2.0 * (1.0 + epsilon) : 2.0;
            h[site, site + 1] = hopping;
            h[site + 1, site] = hopping;
        }
        return h;
    }

    internal static Matrix<Complex> BuildAGenerator(double epsilon, double gamma, int watchedSeat)
    {
        if (watchedSeat < 0 || watchedSeat >= SiteCount)
            throw new ArgumentOutOfRangeException(nameof(watchedSeat));
        Matrix<double> h = BuildHopping(epsilon);
        var generator = Matrix<Complex>.Build.Dense(GeneratorDimension, GeneratorDimension);
        for (int i = 0; i < SiteCount; i++)
            for (int j = 0; j < SiteCount; j++)
            {
                int output = RowMajorIndex(i, j);
                for (int k = 0; k < SiteCount; k++)
                {
                    generator[output, RowMajorIndex(k, j)] += -Complex.ImaginaryOne * h[i, k];
                    generator[output, RowMajorIndex(i, k)] += Complex.ImaginaryOne * h[k, j];
                }
                int zi = i == watchedSeat ? -1 : 1;
                int zj = j == watchedSeat ? -1 : 1;
                generator[output, output] += gamma * (zi * zj - 1);
            }
        return generator;
    }

    private static Matrix<Complex> CentreCellProjector(int seat) =>
        Matrix<Complex>.Build.Diagonal(
            GeneratorDimension,
            GeneratorDimension,
            coordinate =>
            {
                int row = coordinate / SiteCount;
                int column = coordinate % SiteCount;
                return (row == seat) != (column == seat) ? Complex.One : Complex.Zero;
            });

    private static void ValidateEpsilon(double epsilon)
    {
        if (!double.IsFinite(epsilon)
            || epsilon == 0.0
            || 1.0 + epsilon == 1.0
            || Math.Abs(epsilon) <= MinimumResolvableAbsEpsilon
            || Math.Abs(epsilon) > MaxAbsEpsilon)
            throw new ArgumentOutOfRangeException(nameof(epsilon), epsilon,
                $"epsilon must be finite, change the represented value of 1 + epsilon, keep epsilon^2/4 " +
                $"above the Hilbert roundoff model, and satisfy {MinimumResolvableAbsEpsilon:E2} < " +
                $"|epsilon| <= {MaxAbsEpsilon} for the local N=7 witness");
    }

    private static void ValidateGamma(double gamma)
    {
        if (!double.IsFinite(gamma) || gamma <= 0.0)
            throw new ArgumentOutOfRangeException(nameof(gamma), gamma,
                "gamma must be finite and strictly positive for the fixed-gamma local theorem");
    }

    private static void ValidateResolvedSpectralScale(double epsilon, double gamma)
    {
        double generatorNorm = BuildAGenerator(epsilon, gamma, WatchedSeat).FrobeniusNorm();
        double spectralErrorModel = MachineEpsilon
            * GeneratorDimension
            * GeneratorDimension
            * Math.Max(1.0, generatorNorm);
        double predictedSmallClusterRate = gamma * epsilon * epsilon * (1.0 + epsilon) / 2.0;
        if (predictedSmallClusterRate <= spectralErrorModel)
            throw new ArgumentOutOfRangeException(nameof(gamma), gamma,
                $"gamma * epsilon^2 * (1 + epsilon) / 2 must exceed the 49D spectral error model " +
                $"({spectralErrorModel:E3}) for this live finite-precision witness");
    }

    public string DisplayName =>
        $"MissingPhaseRelaxationScaleWitness (fixed N=7, epsilon={Epsilon.ToString("G17", Inv)}, " +
        $"gamma={Gamma.ToString("G17", Inv)})";

    public string Summary
    {
        get
        {
            int spectralPasses = Clusters.SelectMany(cluster => new[] { cluster.Right, cluster.Adjoint })
                .Count(reading => reading.Rank == ClusterRank
                                  && reading.HermitianResidual <= ProjectorErrorModel
                                  && reading.IdempotenceResidual <= ProjectorErrorModel
                                  && reading.RitzResidual <= SpectralErrorModel
                                  && reading.AbsorptionResidual <= SpectralErrorModel);
            return "Fixed N=7 live spectral/geometry companion to the missing-phase relaxation-scale proof: " +
                   $"the independently built A and A† rank-two cluster reads close {spectralPasses}/4 under " +
                   "their eps_machine x norm x dimension error laws; the separate Hilbert route reads " +
                   "c0=0 and c±=epsilon^2/4+epsilon^3/4+O(epsilon^4), composing the five non-I_E Kato " +
                   "gamma-normalized rate coefficients. " + Scope;
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            foreach (var cluster in Clusters)
            {
                string frequency = cluster.FrequencySign < 0 ? "-2 sqrt(2)i" : "+2 sqrt(2)i";
                foreach (var side in new[] { cluster.Right, cluster.Adjoint })
                    yield return new InspectableNode(
                        $"{frequency}, {side.Source}, orthogonal rank-{side.Rank} projector",
                        summary: $"mean decay={side.MeanDecay.ToString("E6", Inv)}, centre light=" +
                                 $"{side.CentreLight.ToString("E6", Inv)}, |Delta-2 gamma w|=" +
                                 $"{side.AbsorptionResidual.ToString("E3", Inv)} " +
                                 $"({(side.AbsorptionResidual / SpectralErrorModel).ToString("E3", Inv)} model units); " +
                                 $"Hermitian/idempotent/Ritz residuals=" +
                                 $"{side.HermitianResidual.ToString("E2", Inv)}/" +
                                 $"{side.IdempotenceResidual.ToString("E2", Inv)}/" +
                                 $"{side.RitzResidual.ToString("E2", Inv)}; seed overlap=" +
                                 side.SeedOverlap.ToString("0.########", Inv));
            }
            yield return new InspectableNode("independent Evd(A) versus Evd(A†) construction",
                summary: ConstructionProvenance + " Projector differences are " +
                         string.Join(", ", Clusters.Select(cluster =>
                             cluster.IndependentProjectorDistance.ToString("E3", Inv))) +
                         "; a copied-left construction would make them identically zero.");
            yield return new InspectableNode("Hilbert 7x7 blind-mode centre weights",
                summary: string.Join("; ", HilbertModes.Select(mode =>
                    $"{mode.Label}: E={mode.Energy.ToString("0.######", Inv)}, " +
                    $"c={mode.CentreWeight.ToString("E7", Inv)}, overlap={mode.SeedOverlap.ToString("0.########", Inv)}")));
            yield return new InspectableNode("five non-I_E dyad coefficients",
                summary: string.Join("; ", DyadCoefficients.Select(reading =>
                    $"{reading.Label}: gamma-normalized {reading.MeasuredCoefficientOverGamma.ToString("0.########", Inv)} " +
                    $"(second order {reading.ExpectedCoefficientOverGamma.ToString("0.0", Inv)}), " +
                    $"rate coefficient {reading.MeasuredRateCoefficient.ToString("G9", Inv)} " +
                    $"at gamma={Gamma.ToString("G17", Inv)}")) +
                    ". " + FullSpaceOnlyFence);
            yield return new InspectableNode("scope and non-claims", summary: Scope);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
