using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using static RCPsiSquared.Diagnostics.Foundation.MissingPhaseReadoutAlgebra;
using Q = RCPsiSquared.Diagnostics.Foundation.ReadoutScalar;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Exact live preparation/readout companion to the N=7 local spectral theorem.
/// Reconstructs both uniform slow dyads and certifies that they exhaust the cluster (exact ranks of the
/// shifted generator and its square, left eigenoperators equal to right ones), the physical-copy partial
/// traces and decoder quadratures, the pair tiers (seven null pairs at epsilon=0, of which the four
/// even ones switch on at first order through two independent bordered solves), the three-spin reading
/// of proof section 7.1 in the full 2^7 space, and derives the leakage germ by an exact bordered
/// eigenvector-derivative solve. It consumes no saved Python result and performs no finite-epsilon
/// eigendecomposition.
/// Proof: docs/proofs/PROOF_MISSING_PHASE_SLOW_READOUT.md.
/// Root: inspect --root missingphasereadout.</summary>
public sealed class MissingPhaseSlowReadoutWitness : IInspectable
{
    public const int SiteCount = 7;
    public BigRational Gamma { get; }
    private readonly Lazy<Snapshot> _reading;
    public Snapshot Reading => _reading.Value;

    /// <summary>The pairs whose physical maps annihilate the uniform residue M_-(0).</summary>
    public static IReadOnlyList<(int A, int B)> UniformNullPairs { get; } =
        new[] { (0, 2), (0, 4), (1, 3), (1, 5), (2, 6), (3, 5), (4, 6) };

    /// <summary>The odd-site pairs, null at every epsilon of the local branch because the blind ray
    /// has no odd-site entry.</summary>
    public static IReadOnlyList<(int A, int B)> OddNullPairs { get; } = new[] { (1, 3), (1, 5), (3, 5) };

    public MissingPhaseSlowReadoutWitness() : this(new BigRational(3, 10)) { }

    public MissingPhaseSlowReadoutWitness(BigRational gamma) : this(gamma, BigRational.Zero) { }

    /// <summary>The seam a mutation needs, not a knob on the physics: a real on-site energy on the centre.
    /// It breaks C h C = -h while e+, v and e- keep their eigenvalues (none has centre weight), so the
    /// uniform dyad checks still pass, the two gates that rest on the chiral pairing (the first-order
    /// pairing and the mixed-cell antisymmetry) go red, and the first-order odd-pair null, which rests on
    /// the blind ray's even-site support alone, stays green. Every physical caller uses zero.</summary>
    internal MissingPhaseSlowReadoutWitness(BigRational gamma, BigRational centreOnsiteEnergy)
    {
        if (gamma.Sign <= 0) throw new ArgumentOutOfRangeException(nameof(gamma), "The local theorem fixes gamma > 0.");
        Gamma = gamma;
        _reading = new(() => Reconstruct(gamma, centreOnsiteEnergy));
    }

    public sealed record Check(string Name, string Expected, string Actual)
    {
        public bool Passes => Expected == Actual;
        public string Detail => $"{Name}: expected {Expected}; computed {Actual}";
    }

    public sealed record Snapshot(
        BigRational EndZzResidue, BigRational EndXxResidue, BigRational ResidueNormSquared,
        BigRational RankOneEndZzResidue, BigRational StationaryResidueNormSquared,
        BigRational SingleCopySinePairNormSquared, BigRational PairedCopySinePairNormSquared,
        BigRational AllPairSineNormSquared, BigRational DecodedHalfTraceNormSquared,
        BigRational LeakageRealCoefficient, BigRational LeakageImaginaryRootTwoCoefficient,
        int ComplexBasisPairCases, int PhysicalPairMismatchCount, int RealPartMutationMismatchCount,
        int ClusterRank, int ClusterSquareRank,
        IReadOnlyList<(int A, int B)> NullPairsAtUniformPoint,
        IReadOnlyList<(int A, int B)> NullPairsAtFirstOrder,
        IReadOnlyDictionary<(int A, int B), BigRational> FirstOrderPairNormSquared,
        BigRational ThreeSpinResidueImaginaryRootTwoCoefficient,
        BigRational ThreeSpinSineRootTwoCoefficient,
        int SeparatingThreeSitePages,
        IReadOnlyList<Check> Checks);

    public string DisplayName => "Missing-phase slow readout (exact N=7 live witness)";
    public string Summary
    {
        get
        {
            var r = Reading;
            return $"N=7 uniform complete slow residue: ZZ={r.EndZzResidue}, XX={r.EndXxResidue}, " +
                   $"norm squared={r.ResidueNormSquared}; {r.NullPairsAtUniformPoint.Count} null pairs at " +
                   $"epsilon=0, {r.NullPairsAtFirstOrder.Count} of them still null at first order; " +
                   $"Y0X1Z3 residue i sqrt(2)({r.ThreeSpinResidueImaginaryRootTwoCoefficient}); " +
                   $"leakage epsilon^2 coefficient {r.LeakageRealCoefficient} + " +
                   $"i sqrt(2)({r.LeakageImaginaryRootTwoCoefficient}); " +
                   $"{r.Checks.Count(c => c.Passes)}/{r.Checks.Count} exact checks. " +
                   "Fixed-gamma local coupling; not a d_out/d_2 lifetime.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("preparation and scope",
                "Coherent equal superposition of the end-odd one-excitation state d=(|0>-|6>)/sqrt(2) " +
                "and its global-flip copy; A_init=B_init=dd-dagger. Centre Z dephasing, XY hopping 2, " +
                $"left end hopping 2(1+epsilon), exact rational gamma={Gamma}. The witness evaluates the " +
                "uniform residue, its first epsilon-derivative and the uniform three-spin trajectory. The proof " +
                "supplies analytic continuation locally at each fixed gamma>0, not a gamma-uniform radius or " +
                "the full Liouvillian gap.");
            foreach (var check in Reading.Checks)
                yield return new InspectableNode(check.Name,
                    $"{check.Detail}; {(check.Passes ? "PASS" : "FAIL")}", provenance: NodeProvenance.Live);
            yield return new InspectableNode("nonlinear reference boundary",
                "The isolated paired decoder contribution is not a density matrix. The d_out and d_2 of the " +
                "experiment's section 7 compare each gamma>0 run with its matched gamma=0 run, a moving unitary " +
                "reference; their whole-signal lifetimes are not asserted.");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    private static Snapshot Reconstruct(BigRational gamma, BigRational centreOnsiteEnergy)
    {
        var checks = new List<Check>();
        void Equal(string name, Q actual, Q expected) => checks.Add(new(name, expected.ToString(), actual.ToString()));
        void ZeroCheck(string name, Q[,] matrix) => checks.Add(new(name, "0", NonzeroCount(matrix).ToString()));

        var b0 = OddSite(0);
        var b1 = OddSite(1);
        var b2 = OddSite(2);
        var v = Scale(Q.RootTwo / 2, Subtract(b0, b2));
        var u = Scale(Q.RootTwo / 2, Add(b0, b2));
        var ePlus = Scale(Q.One / 2, Add(Add(b0, Scale(Q.RootTwo, b1)), b2));
        var eMinus = Scale(Q.One / 2, Add(Subtract(b0, Scale(Q.RootTwo, b1)), b2));
        var dyad1 = Outer(ePlus, v);
        var dyad2 = Outer(v, eMinus);
        var initial = Lift(MissingPhaseOnsetWitness.InitialBlock(SiteCount));
        Q[,] Project(Q[,] input) => Add(Scale(Inner(dyad1, input), dyad1), Scale(Inner(dyad2, input), dyad2));
        var m = Project(initial);
        var mAdjoint = Adjoint(m);
        var cosine = Add(m, mAdjoint);
        var sine = Scale(Q.I, Subtract(mAdjoint, m));
        var h = Lift(MissingPhaseOnsetWitness.Hopping(SiteCount, BigRational.Zero));
        h[3, 3] = centreOnsiteEnergy;
        var lambda = -2 * Q.I * Q.RootTwo;

        Equal("first dyad norm squared", Inner(dyad1, dyad1), 1);
        Equal("second dyad norm squared", Inner(dyad2, dyad2), 1);
        Equal("two independent dyads are orthogonal", Inner(dyad1, dyad2), 0);
        Equal("first preparation overlap", Inner(dyad1, initial), Q.RootTwo / 4);
        Equal("second preparation overlap", Inner(dyad2, initial), Q.RootTwo / 4);
        ZeroCheck("first dyad exact A eigenoperator", Subtract(AAction(h, dyad1, gamma), Scale(lambda, dyad1)));
        ZeroCheck("second dyad exact A eigenoperator", Subtract(AAction(h, dyad2, gamma), Scale(lambda, dyad2)));
        ZeroCheck("residue exact A eigenoperator", Subtract(AAction(h, m, gamma), Scale(lambda, m)));

        // Completeness, recomputed rather than inherited: the shifted 49 x 49 generator and its square
        // both have rank 47, so the cluster is exactly the two dyads with no Jordan partner; both dyads
        // are also eigenoperators of L_A^dagger at the conjugate value, so the Riesz projector is the
        // orthogonal projection Project uses.
        var shift = Zero(49, 49);
        for (int row = 0; row < 7; row++)
            for (int col = 0; col < 7; col++)
            {
                var unit = Zero(7, 7);
                unit[row, col] = Q.One;
                var image = AAction(h, unit, gamma);
                for (int i = 0; i < 7; i++)
                    for (int j = 0; j < 7; j++) shift[7 * i + j, 7 * row + col] = image[i, j];
            }
        for (int diagonal = 0; diagonal < 49; diagonal++) shift[diagonal, diagonal] -= lambda;
        int clusterRank = Rank(shift);
        int clusterSquareRank = Rank(Multiply(shift, shift));
        Equal("slow cluster is exactly two-dimensional: rank(L_A - lambda I) = 47", clusterRank, 47);
        Equal("slow cluster has no Jordan partner: rank((L_A - lambda I)^2) = 47", clusterSquareRank, 47);
        ZeroCheck("first dyad is also a left eigenoperator (L_A^dagger at the conjugate)",
            Subtract(AAdjointAction(h, dyad1, gamma), Scale(lambda.Conjugate, dyad1)));
        ZeroCheck("second dyad is also a left eigenoperator (L_A^dagger at the conjugate)",
            Subtract(AAdjointAction(h, dyad2, gamma), Scale(lambda.Conjugate, dyad2)));

        Equal("residue trace", Trace(m), 0);
        var norm = Inner(m, m).Rational;
        Equal("complete residue norm squared", norm, new BigRational(1, 4));
        var zz = (-2 * (m[0, 0] + m[6, 6])).Rational;
        var xx = (m[0, 6] + m[6, 0]).Rational;
        Equal("physical endpoint ZZ residue", zz, new BigRational(-1, 2));
        Equal("physical endpoint XX residue", xx, new BigRational(-1, 4));
        Equal("bit-level Z0Z6 on both copies equals the matrix formula", PauliExpectation(m, null, "Z0Z6"), zz);
        Equal("bit-level X0X6 on both copies equals the matrix formula", PauliExpectation(m, null, "X0X6"), xx);
        var rankOne = Scale(Inner(dyad1, initial), dyad1);
        var rankOneZz = (-2 * (rankOne[0, 0] + rankOne[6, 6])).Rational;
        Equal("one-dyad mutation changes ZZ", rankOneZz, new BigRational(-1, 4));
        var stationary = Project(Outer(v, v));
        var stationaryNorm = Inner(stationary, stationary).Rational;
        Equal("stationary preparation has no slow residue", stationaryNorm, 0);

        // A complex basis of all traceless inputs, plus i times that basis, detects
        // premature Re(), independently of this particular physical preparation.
        int pairCases = 0, pairMismatches = 0, realPartMismatches = 0;
        for (int row = 0; row < 7; row++)
            for (int col = 0; col < 7; col++)
            {
                if (row == 6 && col == 6) continue;
                foreach (var phase in new[] { Q.One, Q.I })
                {
                    var basis = Zero(7, 7);
                    basis[row, col] = phase;
                    if (row == col) basis[6, 6] = -phase;
                    foreach (var (a, b) in Pairs())
                    {
                        pairCases++;
                        var physical = PhysicalPair(basis, a, b);
                        pairMismatches += NonzeroCount(Subtract(PairMap(basis, a, b), physical));
                        realPartMismatches += NonzeroCount(Subtract(PairMap(basis, a, b, realPartShortcut: true), physical));
                    }
                }
            }
        Equal("complex-linear map vs physical trace on full traceless basis", pairMismatches, 0);
        // Section 8's control, exercised: each of the 84 off-diagonal single-entry inputs breaks the
        // Re(M_ab) shortcut on its own pair, in both symmetric coherence cells.
        Equal("mutation: Re(M_ab) in place of (M_ab+M_ba)/2 fails the same comparison", realPartMismatches, 168);
        BigRational allPairSine = BigRational.Zero;
        foreach (var (a, b) in Pairs())
        {
            ZeroCheck($"physical pair ({a},{b}) residue map", Subtract(PairMap(m, a, b), PhysicalPair(m, a, b)));
            var image = PairMap(sine, a, b);
            allPairSine += Inner(image, image).Rational;
        }
        Equal("all physical pairs lose the uniform sine quadrature", allPairSine, 0);
        var oneCopy = PhysicalPair(sine, 0, 1, flippedCopy: false);
        var twoCopies = PhysicalPair(sine, 0, 1);
        var oneCopyNorm = Inner(oneCopy, oneCopy).Rational;
        var twoCopyNorm = Inner(twoCopies, twoCopies).Rational;
        Equal("removing the flipped copy exposes the sine quadrature", oneCopyNorm, new BigRational(1, 16));
        Equal("paired copies cancel that sine output", twoCopyNorm, 0);

        // The pair tiers at the uniform point. Only the three odd pairs are forced by the blind ray's
        // support; the four even ones are null because e+ and e- agree on even sites there.
        var uniformNull = Pairs().Where(p => PairNull(m, p.A, p.B)).ToList();
        checks.Add(new("pairs that miss the uniform residue", PairList(UniformNullPairs), PairList(uniformNull)));
        var cosineNull = Pairs().Where(p => PairNull(cosine, p.A, p.B)).ToList();
        checks.Add(new("pairs that miss the uniform cosine quadrature", PairList(UniformNullPairs), PairList(cosineNull)));
        // A mutation of the claim rather than of the object: the endpoints (0,6) put into the claimed
        // null set in place of (1,3). The set comparison above already implies the rejection.
        List<(int A, int B)> mutatedNullSet = UniformNullPairs.Select(p => p == (1, 3) ? (0, 6) : p).ToList();
        bool mutantAccepted = mutatedNullSet.All(p => PairNull(m, p.A, p.B));
        checks.Add(new("mutation of the claimed set: the endpoints (0,6) in place of (1,3) are rejected",
            "rejected", mutantAccepted ? "accepted" : "rejected"));
        // The mechanism of the odd-pair null at finite defect, computed from h_r rather than written in:
        // at r = 4/3 the kernel is one ray, it is the relaxation proof's v_r, and it has no odd-site entry,
        // so every entry of E_- X P_v + P_v X E_+^dagger on an odd pair carries a zero factor.
        checks.AddRange(BlindRayChecks(MissingPhaseOnsetWitness.Hopping(SiteCount, new BigRational(1, 3)),
            new BigRational(4, 3)));

        var dc = DecodeA(cosine);
        var ds = DecodeA(sine);
        ZeroCheck("cosine decoder vs physical controlled flips", Subtract(dc, PhysicalDecoder(cosine)));
        ZeroCheck("sine decoder vs physical controlled flips", Subtract(ds, PhysicalDecoder(sine)));
        var dc2 = Multiply(dc, dc);
        var ds2 = Multiply(ds, ds);
        ZeroCheck("decoded cosine cubic identity", Subtract(Multiply(dc2, dc), Scale(Q.One / 4, dc)));
        ZeroCheck("decoded sine cubic identity", Subtract(Multiply(ds2, ds), Scale(Q.One / 4, ds)));
        ZeroCheck("decoded c^2 s coefficient", Subtract(
            Add(Add(Multiply(dc2, ds), Multiply(Multiply(dc, ds), dc)), Multiply(ds, dc2)), Scale(Q.One / 4, ds)));
        ZeroCheck("decoded c s^2 coefficient", Subtract(
            Add(Add(Multiply(dc, ds2), Multiply(Multiply(ds, dc), ds)), Multiply(ds2, dc)), Scale(Q.One / 4, dc)));
        Equal("decoded cosine trace", Trace(dc), 0);
        Equal("decoded sine trace", Trace(ds), 0);
        Equal("decoded cosine squared trace", Trace(dc2), Q.One / 2);
        Equal("decoded sine squared trace", Trace(ds2), Q.One / 2);
        Equal("decoded mixed trace", Trace(Multiply(dc, ds)), 0);
        // Cubic identities force eigenvalues 0,+/-1/2 at c^2+s^2=1.
        // Trace-square and trace then force one of each nonzero eigenvalue.
        var decodedNormSquared = (Trace(dc2) / 2).Rational;

        var pc = Zero(7, 7);
        pc[3, 3] = 1;
        var k = Subtract(Scale(-Q.I, h), Scale(2 * (Q)gamma, pc));
        var hPrime = Subtract(Lift(MissingPhaseOnsetWitness.Hopping(7, BigRational.One)),
            Lift(MissingPhaseOnsetWitness.Hopping(7, BigRational.Zero)));
        var kPrime = Scale(-Q.I, hPrime);
        // The equivalent bordered system of proof section 6: right-hand side (lambda' I - K') e and
        // border column +e; its extra compatibility multiplier must be zero.
        Q[,] BorderedDerivative(string label, Q eigenvalue, Q[,] border, Q eigenvaluePrime)
        {
            var system = Zero(8, 8);
            var rhs = Zero(8, 1);
            var source = Subtract(Scale(eigenvaluePrime, border), Multiply(kPrime, border));
            for (int i = 0; i < 7; i++)
            {
                for (int j = 0; j < 7; j++) system[i, j] = k[i, j] - (i == j ? eigenvalue : Q.Zero);
                system[i, 7] = border[i, 0];
                system[7, i] = border[i, 0]; // transpose: K is complex symmetric
                rhs[i, 0] = source[i, 0];
            }
            var solved = Solve(system, rhs);
            ZeroCheck($"bordered eigenvector derivative residual{label}", Subtract(Multiply(system, solved), rhs));
            Equal($"bordered compatibility multiplier{label}", solved[7, 0], 0);
            var derivative = Zero(7, 1);
            for (int i = 0; i < 7; i++) derivative[i, 0] = solved[i, 0];
            return derivative;
        }
        var lambdaPrime = Inner(ePlus, Multiply(kPrime, ePlus));
        Equal("eigenvalue derivative from the left/right pairing", lambdaPrime, -Q.I * Q.RootTwo / 2);
        // e+(eps): the eigenvector of K at lambda_-(eps) continuing e+, normalized by e+^T e+(eps) = 1.
        var ePlusPrime = BorderedDerivative("", lambda, ePlus, lambdaPrime);
        // Differentiate b(r)/sqrt(b(r)^T b(r)) at r=1, directly from its entries.
        var raw = Zero(7, 1);
        var rawPrime = Zero(7, 1);
        int[] entries = { 1, 0, -1, 0, 1, 0, -1 };
        int[] derivatives = { 0, 0, -1, 0, 1, 0, -1 };
        for (int j = 0; j < 7; j++) { raw[j, 0] = entries[j]; rawPrime[j, 0] = derivatives[j]; }
        var vPrime = Subtract(Scale(Q.One / 2, rawPrime), Scale(Inner(raw, rawPrime) / 8, raw));
        var qvPrime = ReflectionEven(vPrime);
        var qePlusPrime = ReflectionEven(ePlusPrime);
        var a0 = Inner(v, b0);
        var overlap = Inner(ePlus, b0);
        var leakage = 2 * a0 * overlap * Inner(qvPrime, qePlusPrime);
        Equal("leakage epsilon^2 coefficient from derivative vectors", leakage,
            -Q.One / 16 + Q.I * Q.RootTwo * (Q)gamma / 16);

        // First order, without assuming the chiral pairing: the conjugate branch e-(eps), the eigenvector
        // of K at lambda_+(eps) continuing e-, gets its own bordered solve, and the residue derivative is
        // assembled from the projector form E_- A P_v + P_v A E_+^dagger with E_- = e+(eps) e+(eps)^T and
        // E_+ = e-(eps) e-(eps)^T over their bilinear norms (K is complex symmetric). The normalizations
        // e+^T e+(eps) = e-^T e-(eps) = 1 keep those norms stationary at first order.
        // Every derivative here is affine in gamma: gamma sits only in the centre diagonal entry of each
        // bordered matrix, whose cofactor vanishes, so both bordered determinants are gamma-free.
        var lambdaPlusPrime = Inner(eMinus, Multiply(kPrime, eMinus));
        Equal("conjugate eigenvalue derivative", lambdaPlusPrime, Q.I * Q.RootTwo / 2);
        var eMinusPrime = BorderedDerivative(" (conjugate branch)", lambda.Conjugate, eMinus, lambdaPlusPrime);
        var pairedEPlusPrime = Zero(7, 1);
        for (int i = 0; i < 7; i++) pairedEPlusPrime[i, 0] = (i % 2 == 0 ? Q.One : -Q.One) * ePlusPrime[i, 0].Conjugate;
        ZeroCheck("chiral pairing at first order: e-' = C conj(e+') from two independent solves",
            Subtract(eMinusPrime, pairedEPlusPrime));
        var eMinusProjector = OuterTranspose(ePlus, ePlus);
        var ePlusAdjoint = OuterTranspose(eMinus, eMinus);
        var blindProjector = OuterTranspose(v, v);
        ZeroCheck("projector form reproduces the uniform residue", Subtract(
            Add(Multiply(Multiply(eMinusProjector, initial), blindProjector),
                Multiply(Multiply(blindProjector, initial), ePlusAdjoint)), m));
        var eMinusPrimeConjugate = Conjugate(eMinusPrime);
        var eMinusProjectorPrime = Add(OuterTranspose(ePlusPrime, ePlus), OuterTranspose(ePlus, ePlusPrime));
        var ePlusAdjointPrime = Add(OuterTranspose(eMinusPrimeConjugate, eMinus), OuterTranspose(eMinus, eMinusPrimeConjugate));
        var blindProjectorPrime = Add(OuterTranspose(vPrime, v), OuterTranspose(v, vPrime));
        var firstOrder = Add(
            Add(Multiply(Multiply(eMinusProjectorPrime, initial), blindProjector),
                Multiply(Multiply(eMinusProjector, initial), blindProjectorPrime)),
            Add(Multiply(Multiply(blindProjectorPrime, initial), ePlusAdjoint),
                Multiply(Multiply(blindProjector, initial), ePlusAdjointPrime)));
        Equal("first-order residue is traceless", Trace(firstOrder), 0);
        Equal("first-order leakage vanishes", LeakageTrace(firstOrder), 0);
        int mixedSymmetricCells = 0;
        for (int a = 0; a < 7; a += 2)
            for (int b = 1; b < 7; b += 2)
                if (!(firstOrder[a, b] + firstOrder[b, a]).IsZero) mixedSymmetricCells++;
        Equal("first-order mixed even/odd cells are antisymmetric", mixedSymmetricCells, 0);
        var firstOrderNull = UniformNullPairs.Where(p => PairNull(firstOrder, p.A, p.B)).ToList();
        checks.Add(new("of the seven uniform null pairs, those still null at first order",
            PairList(OddNullPairs), PairList(firstOrderNull)));
        var g = (Q)gamma;
        var firstOrderNorms = new Dictionary<(int A, int B), BigRational>();
        void FirstOrderPair(int a, int b, Q diagonalSum, Q coherence)
        {
            var s = firstOrder[a, a] + firstOrder[b, b];
            var c = (firstOrder[a, b] + firstOrder[b, a]) / 2;
            Equal($"first-order pair ({a},{b}) diagonal sum", s, diagonalSum);
            Equal($"first-order pair ({a},{b}) symmetric coherence", c, coherence);
            // |s|^2 + 2|v|^2, the squared Frobenius norm of the pair map's first-order term; rational
            // on the physical chain, and left out rather than rounded if a mutation makes it irrational.
            var normSquared = s.Conjugate * s + 2 * c.Conjugate * c;
            if (normSquared.B.IsZero && normSquared.A.Im.IsZero) firstOrderNorms[(a, b)] = normSquared.A.Re;
        }
        FirstOrderPair(0, 2, Q.One / 8, -(Q)3 / 16);
        FirstOrderPair(4, 6, -Q.One / 8, Q.One / 16);
        FirstOrderPair(0, 4, -(Q.One + Q.I * Q.RootTwo * g) / 8, (Q.One - Q.I * Q.RootTwo * g) / 16);
        FirstOrderPair(2, 6, (Q.One + Q.I * Q.RootTwo * g) / 8, (Q.One + Q.I * Q.RootTwo * g) / 16);

        // Proof section 7.1 in the full 2^7 space. Bit j of a word is site j and bit 0 has Z = +1.
        var threeSpin = PauliExpectation(m, null, "Y0X1Z3");
        Equal("three-spin residue: Y0X1Z3 on both copies", threeSpin, Q.I * Q.RootTwo / 8);
        Equal("every tag Z_k, k = 2..6, reads the same residue",
            Enumerable.Range(2, 5).Count(site => !(PauliExpectation(m, null, $"Y0X1Z{site}") - threeSpin).IsZero), 0);
        Equal("control: the untagged pair Y0X1 reads zero on both copies", PauliExpectation(m, null, "Y0X1"), 0);
        Equal("control: an X tag, Y0X1X3, reads zero", PauliExpectation(m, null, "Y0X1X3"), 0);
        Equal("one copy: the untagged pair Y0X1 reads what Y0X1Z3 reads on two",
            PauliExpectation(m, null, "Y0X1", flippedCopy: false), threeSpin);
        Equal("one copy: the endpoint ZZ residue is unchanged", PauliExpectation(m, null, "Z0Z6", flippedCopy: false), zz);
        Equal("one copy: the endpoint XX residue is unchanged", PauliExpectation(m, null, "X0X6", flippedCopy: false), xx);
        var sineReading = PauliExpectation(sine, null, "Y0X1Z3");
        Equal("Y0X1Z3 reads the sine quadrature", sineReading, Q.RootTwo / 4);
        Equal("Y0X1Z3 reads nothing of the cosine quadrature", PauliExpectation(cosine, null, "Y0X1Z3"), 0);
        Equal("every one- and two-site page of the sine quadrature vanishes",
            Subsets(1).Concat(Subsets(2)).Count(sites => NonzeroCount(PhysicalPage(sine, null, sites)) > 0), 0);

        // h v = 0, h u = 2 sqrt(2) b1, h b1 = 2 sqrt(2) u: counted separately so no two failures cancel.
        // h maps the centre-free span {v, u, b1} into itself, which is why d(theta) never charges the
        // centre; that d(theta) has no centre entry is true by construction and is not re-checked.
        var motionResiduals = NonzeroCount(Multiply(h, v))
            + NonzeroCount(Subtract(Multiply(h, u), Scale(2 * Q.RootTwo, b1)))
            + NonzeroCount(Subtract(Multiply(h, b1), Scale(2 * Q.RootTwo, u)));
        Equal("d(theta) = (v + cos theta u - i sin theta b1)/sqrt(2) solves i dd/dt = h d at theta = 2 sqrt(2) t",
            motionResiduals, 0);
        var phases = new (Q C, Q S)[]
        {
            (1, 0), (0, 1), (-1, 0), (0, -1), (new BigRational(3, 5), new BigRational(4, 5)),
            (new BigRational(-5, 13), new BigRational(12, 13)), (Q.RootTwo / 2, Q.RootTwo / 2),
        };
        int mixtureMisses = 0, threeSpinMisses = 0, endpointMisses = 0, untaggedMisses = 0,
            tagMisses = 0, slopeMisses = 0;
        foreach (var (c, s) in phases)
        {
            var d = Scale(Q.RootTwo / 2, Subtract(Add(v, Scale(c, u)), Scale(Q.I * s, b1)));
            var state = Outer(d, d);
            var intercopy = Scale(new BigRational(5, 11), state); // any B weight: every reading here has at most three sites
            mixtureMisses += ConvexMixtureMismatches(d, state);
            var reading = PauliExpectation(state, intercopy, "Y0X1Z3");
            if (!(reading - s * (1 + c) * Q.RootTwo / 4).IsZero) threeSpinMisses++;
            if (!(PauliExpectation(state, intercopy, "Z0Z6") - (1 - (1 + c) * (1 + c) / 2)).IsZero) endpointMisses++;
            if (!PauliExpectation(state, intercopy, "Y0X1").IsZero) untaggedMisses++;
            for (int site = 2; site < 7; site++)
                if (!(PauliExpectation(state, intercopy, $"Y0X1Z{site}") - reading).IsZero) tagMisses++;
            var slope = 2 * PauliExpectation(state, intercopy, "Y0X1Z6") - 2 * PauliExpectation(state, intercopy, "X0Y1Z6")
                      + 2 * PauliExpectation(state, intercopy, "X5Y6Z0") - 2 * PauliExpectation(state, intercopy, "Y5X6Z0");
            // d/dt [1 - (1+c)^2/2] at theta = 2 sqrt(2) t is 2 sqrt(2) s (1+c) = 8 <Y0X1Z3>.
            if (!(slope - 8 * reading).IsZero || !(slope - 2 * Q.RootTwo * s * (1 + c)).IsZero) slopeMisses++;
        }
        // A quadratic polynomial in (cos theta, sin theta) that vanishes at five points of the unit circle
        // vanishes on all of it, so seven points certify each reading at every theta.
        Equal("the physical state is the mixture (1+eta)/2 |d_s><d_s| + (1-eta)/2 |d_a><d_a| at seven phases", mixtureMisses, 0);
        Equal("<Y0X1Z3> = sin theta (1 + cos theta)/(2 sqrt(2)) at seven phases", threeSpinMisses, 0);
        Equal("<Z0Z6> = 1 - (1 + cos theta)^2/2 at seven phases", endpointMisses, 0);
        Equal("<Y0X1> = 0 at seven phases", untaggedMisses, 0);
        Equal("<Y0X1Z_k> = <Y0X1Z3> for k = 2..6 at seven phases", tagMisses, 0);
        Equal("<i[H,Z0Z6]> = 8 <Y0X1Z3> = d<Z0Z6>/dt at seven phases", slopeMisses, 0);
        Equal("i[H,Z0Z6] = 2(Y0X1 - X0Y1)Z6 + 2(X5Y6 - Y5X6)Z0 on all 128 basis words", SlopeOperatorMismatches(), 0);

        var dPlus = Scale(Q.RootTwo / 2, Subtract(v, Scale(Q.I, b1)));  // theta = pi/2, t+ = pi/(4 sqrt 2)
        var dMinus = Scale(Q.RootTwo / 2, Add(v, Scale(Q.I, b1)));      // theta = 3 pi/2, t- = 3 pi/(4 sqrt 2)
        var aPlus = Outer(dPlus, dPlus);
        var aMinus = Outer(dMinus, dMinus);
        ZeroCheck("snapshots: A(t+) - A(t-) = 2 C_(pi/2)", Subtract(Subtract(aPlus, aMinus), Scale(2, sine)));
        var bPlus = Scale(new BigRational(2, 3), aPlus);  // unequal B weights stand in for e^(-2 gamma t+-)
        var bMinus = Scale(new BigRational(1, 7), aMinus);
        bool Separates(int[] sites) =>
            NonzeroCount(Subtract(PhysicalPage(aPlus, bPlus, sites), PhysicalPage(aMinus, bMinus, sites))) > 0;
        Equal("snapshots: all 21 pair pages agree", Subsets(2).Count(Separates), 0);
        var separating = Subsets(3).Where(Separates).Select(t => (t[0], t[1], t[2])).ToHashSet();
        Equal("snapshots: three-site pages that tell t+ from t-", separating.Count, 24);
        Equal("those pages are exactly the ones holding site 1 or 5 and an even site",
            Subsets(3).Count(t => separating.Contains((t[0], t[1], t[2])) !=
                (t.Any(site => site == 1 || site == 5) && t.Any(site => site % 2 == 0))), 0);
        Equal("snapshot t+ three-spin reading", PauliExpectation(aPlus, bPlus, "Y0X1Z3"), Q.RootTwo / 4);
        Equal("snapshot t- three-spin reading", PauliExpectation(aMinus, bMinus, "Y0X1Z3"), -Q.RootTwo / 4);
        Equal("snapshot t+ endpoint ZZ", PauliExpectation(aPlus, bPlus, "Z0Z6"), Q.One / 2);
        Equal("snapshot t- endpoint ZZ", PauliExpectation(aMinus, bMinus, "Z0Z6"), Q.One / 2);
        var noA = Zero(7, 7);
        Equal("F70: the intercopy block is invisible to every page of at most four sites",
            Enumerable.Range(1, 4).SelectMany(size => Subsets(size))
                .Count(sites => NonzeroCount(PhysicalPage(noA, aPlus, sites)) > 0), 0);
        checks.Add(new("control: the intercopy block is visible on some five-site page", "visible",
            Subsets(5).Any(sites => NonzeroCount(PhysicalPage(noA, aPlus, sites)) > 0) ? "visible" : "invisible"));

        return new(zz, xx, norm, rankOneZz, stationaryNorm, oneCopyNorm, twoCopyNorm,
            allPairSine, decodedNormSquared, leakage.A.Re, leakage.B.Im,
            pairCases, pairMismatches, realPartMismatches, clusterRank, clusterSquareRank,
            uniformNull.AsReadOnly(), firstOrderNull.AsReadOnly(), firstOrderNorms,
            threeSpin.B.Im, sineReading.B.Re, separating.Count, checks.AsReadOnly());
    }

    private static string PairList(IEnumerable<(int A, int B)> pairs) =>
        string.Join(" ", pairs.Select(p => $"({p.A},{p.B})"));

    // The pair map of section 4 vanishes exactly when its diagonal sum and symmetric coherence both do.
    private static bool PairNull(Q[,] x, int a, int b) => (x[a, a] + x[b, b]).IsZero && (x[a, b] + x[b, a]).IsZero;

    /// <summary>The three blind-ray gates, on any one-excitation hopping matrix, by exact elimination: the
    /// kernel is one ray, the ray has no odd-site entry, and it is (1,0,-r,0,r,0,-r) up to scale. The
    /// witness runs them on the physical chain at r = 4/3. The tests feed them mutated chains: a single
    /// bond between two even sites removes the kernel, and two tuned even-even bonds keep a one-dimensional
    /// kernel whose ray has odd-site entries, which the odd-site gate must reject while the dimension gate
    /// passes.</summary>
    internal static IReadOnlyList<Check> BlindRayChecks(GaussianRational[,] hopping, BigRational ratio)
    {
        var kernel = Nullspace(Lift(hopping));
        int oddSiteEntries = -1, closedMismatches = 7;
        if (kernel.Count == 1)
        {
            var ray = kernel[0];
            oddSiteEntries = Enumerable.Range(0, ray.GetLength(0)).Count(j => j % 2 == 1 && !ray[j, 0].IsZero);
            var closedRay = Zero(7, 1);
            var r = (Q)ratio;
            closedRay[0, 0] = 1; closedRay[2, 0] = -r; closedRay[4, 0] = r; closedRay[6, 0] = -r;
            closedMismatches = NonzeroCount(Subtract(ray, Scale(ray[0, 0], closedRay)));
        }
        return new Check[]
        {
            new($"finite-defect hopping at r = {ratio} has a one-dimensional kernel", "1", kernel.Count.ToString()),
            new("that blind ray has no odd-site entry", "0", oddSiteEntries.ToString()),
            new("it is the relaxation proof's (1,0,-r,0,r,0,-r) up to scale", "0", closedMismatches.ToString()),
        };
    }

    private static Q[,] OddSite(int j)
    {
        var v = Zero(7, 1);
        v[j, 0] = Q.RootTwo / 2;
        v[6 - j, 0] = -Q.RootTwo / 2;
        return v;
    }

    private static IEnumerable<(int A, int B)> Pairs()
    {
        for (int a = 0; a < 7; a++)
            for (int b = a + 1; b < 7; b++) yield return (a, b);
    }

    private static IEnumerable<int[]> Subsets(int size)
    {
        IEnumerable<int[]> Extend(int start, int[] prefix)
        {
            if (prefix.Length == size) { yield return prefix; yield break; }
            for (int site = start; site < 7; site++)
                foreach (var subset in Extend(site + 1, prefix.Append(site).ToArray())) yield return subset;
        }
        return Extend(0, Array.Empty<int>());
    }

    private static Q[,] AAction(Q[,] h, Q[,] x, BigRational gamma)
    {
        var result = Scale(-Q.I, Subtract(Multiply(h, x), Multiply(x, h)));
        for (int i = 0; i < 7; i++)
            for (int j = 0; j < 7; j++)
                if ((i == 3) != (j == 3)) result[i, j] -= 2 * (Q)gamma * x[i, j];
        return result;
    }

    // Hilbert-Schmidt adjoint: +i[h, X] for Hermitian h; the cellwise dephasing is self-adjoint.
    private static Q[,] AAdjointAction(Q[,] h, Q[,] x, BigRational gamma)
    {
        var result = Scale(Q.I, Subtract(Multiply(h, x), Multiply(x, h)));
        for (int i = 0; i < 7; i++)
            for (int j = 0; j < 7; j++)
                if ((i == 3) != (j == 3)) result[i, j] -= 2 * (Q)gamma * x[i, j];
        return result;
    }

    // With realPartShortcut the map uses Re(M_ab) instead of (M_ab+M_ba)/2: the premature shortcut
    // proof section 4 warns against, kept only as the mutation the linear-map comparison must reject.
    private static Q[,] PairMap(Q[,] m, int a, int b, bool realPartShortcut = false)
    {
        if (!Trace(m).IsZero) throw new ArgumentException("The residue map requires a traceless input.");
        var result = Zero(4, 4);
        var s = (m[a, a] + m[b, b]) / 2;
        result[0, 0] = result[3, 3] = -s;
        result[1, 1] = result[2, 2] = s;
        result[1, 2] = result[2, 1] = realPartShortcut ? RealPart(m[a, b]) : (m[a, b] + m[b, a]) / 2;
        return result;
    }

    // Independent physical partial trace of the one-particle and one-hole copies.
    private static Q[,] PhysicalPair(Q[,] m, int a, int b, bool flippedCopy = true)
    {
        var result = Zero(4, 4);
        int copies = flippedCopy ? 2 : 1;
        int rest = 127 ^ (1 << a) ^ (1 << b);
        for (int copy = 0; copy < copies; copy++)
            for (int i = 0; i < 7; i++)
                for (int j = 0; j < 7; j++)
                {
                    int row = (1 << i) ^ (copy == 1 ? 127 : 0);
                    int col = (1 << j) ^ (copy == 1 ? 127 : 0);
                    if ((row & rest) != (col & rest)) continue;
                    int r = 2 * ((row >> a) & 1) + ((row >> b) & 1);
                    int c = 2 * ((col >> a) & 1) + ((col >> b) & 1);
                    result[r, c] += m[i, j] / copies;
                }
        return result;
    }

    // The physical operator 1/2 W [[A, B], [B, A]] W^dagger as (row word, column word, value): copy one on
    // the words 1<<i, the flipped copy on 127^(1<<i), B between them. Without the flipped copy, A alone.
    private static IEnumerable<(int Row, int Col, Q Value)> PhysicalEntries(Q[,] a, Q[,]? b, bool flippedCopy)
    {
        for (int i = 0; i < 7; i++)
            for (int j = 0; j < 7; j++)
            {
                int row = 1 << i, col = 1 << j;
                if (!flippedCopy)
                {
                    if (!a[i, j].IsZero) yield return (row, col, a[i, j]);
                    continue;
                }
                var half = a[i, j] / 2;
                if (!half.IsZero)
                {
                    yield return (row, col, half);
                    yield return (row ^ 127, col ^ 127, half);
                }
                if (b is null || b[i, j].IsZero) continue;
                var cross = b[i, j] / 2;
                yield return (row, col ^ 127, cross);
                yield return (row ^ 127, col, cross);
            }
    }

    // A Pauli word written as letter-site pairs, "Y0X1Z3" for Y on site 0, X on site 1 and Z on site 3.
    private static (int Site, char Letter)[] Pauli(string word)
    {
        var ops = new (int Site, char Letter)[word.Length / 2];
        for (int t = 0; t < ops.Length; t++) ops[t] = (word[2 * t + 1] - '0', word[2 * t]);
        return ops;
    }

    // P|word> = factor |target>, with Y|0> = i|1>, Y|1> = -i|0>, Z|1> = -|1>.
    private static (int Target, Q Factor) PauliAct(int word, (int Site, char Letter)[] ops)
    {
        int target = word;
        var factor = Q.One;
        foreach (var (site, letter) in ops)
        {
            bool up = ((word >> site) & 1) == 0;
            switch (letter)
            {
                case 'X': target ^= 1 << site; break;
                case 'Y': factor *= up ? Q.I : -Q.I; target ^= 1 << site; break;
                case 'Z': if (!up) factor = -factor; break;
                default: throw new ArgumentException($"Unknown Pauli letter {letter}.");
            }
        }
        return (target, factor);
    }

    // Tr[P rho] = sum over rows of factor(row) rho[row, target(row)].
    private static Q PauliExpectation(Q[,] a, Q[,]? b, string word, bool flippedCopy = true)
    {
        var ops = Pauli(word);
        var total = Q.Zero;
        foreach (var (row, col, value) in PhysicalEntries(a, b, flippedCopy))
        {
            var (target, factor) = PauliAct(row, ops);
            if (target == col) total += factor * value;
        }
        return total;
    }

    // The reduced operator on the chosen sites, the other five or fewer traced out bit by bit.
    private static Q[,] PhysicalPage(Q[,] a, Q[,]? b, int[] sites, bool flippedCopy = true)
    {
        int mask = 0;
        foreach (var site in sites) mask |= 1 << site;
        var page = Zero(1 << sites.Length, 1 << sites.Length);
        foreach (var (row, col, value) in PhysicalEntries(a, b, flippedCopy))
        {
            if (((row ^ col) & (127 ^ mask)) != 0) continue;
            int r = 0, c = 0;
            for (int t = 0; t < sites.Length; t++)
            {
                r |= ((row >> sites[t]) & 1) << t;
                c |= ((col >> sites[t]) & 1) << t;
            }
            page[r, c] += value;
        }
        return page;
    }

    // rho(eta) = 1/2 W [[A, eta A], [eta A, A]] W^dagger against (1+eta)/2 P_s + (1-eta)/2 P_a with
    // d_s,a = (|d> +- F|d>)/sqrt(2): the eta-free part must be (P_s + P_a)/2, the eta part (P_s - P_a)/2.
    private static int ConvexMixtureMismatches(Q[,] d, Q[,] state)
    {
        var symmetric = new Dictionary<int, Q>();
        var antisymmetric = new Dictionary<int, Q>();
        for (int i = 0; i < 7; i++)
        {
            var amplitude = d[i, 0] * Q.RootTwo / 2;
            symmetric[1 << i] = amplitude;
            symmetric[127 ^ (1 << i)] = amplitude;
            antisymmetric[1 << i] = amplitude;
            antisymmetric[127 ^ (1 << i)] = -amplitude;
        }
        var still = PhysicalEntries(state, null, true).ToDictionary(e => (e.Row, e.Col), e => e.Value);
        var moving = PhysicalEntries(Zero(7, 7), state, true).ToDictionary(e => (e.Row, e.Col), e => e.Value);
        int misses = 0;
        foreach (var x in symmetric.Keys)
            foreach (var y in symmetric.Keys)
            {
                var ps = symmetric[x] * symmetric[y].Conjugate;
                var pa = antisymmetric[x] * antisymmetric[y].Conjugate;
                if (!((still.GetValueOrDefault((x, y), Q.Zero)) - (ps + pa) / 2).IsZero) misses++;
                if (!((moving.GetValueOrDefault((x, y), Q.Zero)) - (ps - pa) / 2).IsZero) misses++;
            }
        return misses;
    }

    // i[H, Z0Z6] from the hopping rule alone (an excitation crosses a bond with amplitude 2), against the
    // four Pauli terms through PauliAct: two independent routes on every basis word.
    private static int SlopeOperatorMismatches()
    {
        var terms = new (Q Coefficient, (int Site, char Letter)[] Ops)[]
        {
            (2, Pauli("Y0X1Z6")), (-2, Pauli("X0Y1Z6")), (2, Pauli("X5Y6Z0")), (-2, Pauli("Y5X6Z0")),
        };
        static int EndSign(int word) => ((word ^ (word >> 6)) & 1) == 0 ? 1 : -1;
        int mismatches = 0;
        for (int word = 0; word < 128; word++)
        {
            var commutator = new Dictionary<int, Q>();
            for (int bond = 0; bond < 6; bond++)
            {
                if (((word >> bond) & 1) == ((word >> (bond + 1)) & 1)) continue;
                int moved = word ^ (3 << bond);
                commutator[moved] = commutator.GetValueOrDefault(moved, Q.Zero)
                    + Q.I * 2 * (EndSign(word) - EndSign(moved));
            }
            var paulis = new Dictionary<int, Q>();
            foreach (var (coefficient, ops) in terms)
            {
                var (target, factor) = PauliAct(word, ops);
                paulis[target] = paulis.GetValueOrDefault(target, Q.Zero) + coefficient * factor;
            }
            foreach (var target in commutator.Keys.Union(paulis.Keys))
                if (!(commutator.GetValueOrDefault(target, Q.Zero) - paulis.GetValueOrDefault(target, Q.Zero)).IsZero)
                    mismatches++;
        }
        return mismatches;
    }

    private static Q[,] DecodeA(Q[,] m)
    {
        int[] outer = { 0, 1, 2, 4, 5, 6 };
        var result = Zero(7, 7);
        for (int i = 0; i < 6; i++)
            for (int j = 0; j < 6; j++) result[i, j] = m[outer[i], outer[j]];
        result[6, 6] = m[3, 3]; // B=0 for this isolated A contribution
        return result;
    }

    private static Q[,] PhysicalDecoder(Q[,] m)
    {
        int[] outer = { 0, 1, 2, 4, 5, 6 };
        var indices = outer.Select((site, index) => (Word: 1 << site, Index: index))
            .ToDictionary(item => item.Word, item => item.Index);
        indices[119] = 6;
        var result = Zero(7, 7);
        for (int copy = 0; copy < 2; copy++)
            for (int i = 0; i < 7; i++)
                for (int j = 0; j < 7; j++)
                {
                    int row = (1 << i) ^ (copy == 1 ? 127 : 0);
                    int col = (1 << j) ^ (copy == 1 ? 127 : 0);
                    if ((row & 8) != 0) row ^= 119; // CNOTs: centre controls every outer bit
                    if ((col & 8) != 0) col ^= 119;
                    if ((row & 8) != (col & 8)) continue; // trace the centre
                    result[indices[row & 119], indices[col & 119]] += m[i, j] / 2;
                }
        return result;
    }

    private static Q[,] ReflectionEven(Q[,] v)
    {
        var result = Zero(7, 1);
        for (int i = 0; i < 7; i++) result[i, 0] = (v[i, 0] + v[6 - i, 0]) / 2;
        return result;
    }

    // Tr(Q_R X) with Q_R = (I + R)/2 and R|j> = |6-j>.
    private static Q LeakageTrace(Q[,] x)
    {
        var total = Q.Zero;
        for (int i = 0; i < 7; i++) total += (x[i, i] + x[6 - i, i]) / 2;
        return total;
    }
}
