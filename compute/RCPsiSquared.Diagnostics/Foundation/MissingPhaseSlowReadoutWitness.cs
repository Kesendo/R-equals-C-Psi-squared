using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using static RCPsiSquared.Diagnostics.Foundation.MissingPhaseReadoutAlgebra;
using Q = RCPsiSquared.Diagnostics.Foundation.ReadoutScalar;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Exact live preparation/readout companion to the N=7 local spectral theorem.
/// Reconstructs both uniform slow dyads, physical-copy partial traces and decoder quadratures,
/// and derives the leakage germ by an exact bordered eigenvector-derivative solve.
/// It consumes no saved Python result and performs no finite-epsilon eigendecomposition.
/// Proof: docs/proofs/PROOF_MISSING_PHASE_SLOW_READOUT.md.
/// Root: inspect --root missingphasereadout.</summary>
public sealed class MissingPhaseSlowReadoutWitness : IInspectable
{
    public const int SiteCount = 7;
    public BigRational Gamma { get; }
    private readonly Lazy<Snapshot> _reading;
    public Snapshot Reading => _reading.Value;

    public MissingPhaseSlowReadoutWitness() : this(new BigRational(3, 10)) { }

    public MissingPhaseSlowReadoutWitness(BigRational gamma)
    {
        if (gamma.Sign <= 0) throw new ArgumentOutOfRangeException(nameof(gamma), "The local theorem fixes gamma > 0.");
        Gamma = gamma;
        _reading = new(() => Reconstruct(gamma));
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
        int ComplexBasisPairCases, int PhysicalPairMismatchCount, int OddPairZeros,
        IReadOnlyList<Check> Checks);

    public string DisplayName => "Missing-phase slow readout (exact N=7 live witness)";
    public string Summary
    {
        get
        {
            var r = Reading;
            return $"N=7 uniform complete slow residue: ZZ={r.EndZzResidue}, XX={r.EndXxResidue}, " +
                   $"norm squared={r.ResidueNormSquared}; leakage epsilon^2 coefficient " +
                   $"{r.LeakageRealCoefficient} + i sqrt(2)({r.LeakageImaginaryRootTwoCoefficient}); " +
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
                "uniform residue and a derivative at epsilon=0. The proof supplies analytic continuation " +
                "locally at each fixed gamma>0, not a gamma-uniform radius or the full Liouvillian gap.");
            foreach (var check in Reading.Checks)
                yield return new InspectableNode(check.Name,
                    $"{check.Detail}; {(check.Passes ? "PASS" : "FAIL")}", provenance: NodeProvenance.Live);
            yield return new InspectableNode("nonlinear reference boundary",
                "The isolated paired decoder contribution is not a density matrix. Historical d_out and " +
                "d_2 compare with a continuing unitary reference; their whole-signal lifetimes are not asserted.");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    private static Snapshot Reconstruct(BigRational gamma)
    {
        var checks = new List<Check>();
        void Equal(string name, Q actual, Q expected) => checks.Add(new(name, expected.ToString(), actual.ToString()));
        void ZeroCheck(string name, Q[,] matrix) => checks.Add(new(name, "0", NonzeroCount(matrix).ToString()));

        var b0 = OddSite(0);
        var b1 = OddSite(1);
        var b2 = OddSite(2);
        var v = Scale(Q.RootTwo / 2, Subtract(b0, b2));
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
        var lambda = -2 * Q.I * Q.RootTwo;

        Equal("first dyad norm squared", Inner(dyad1, dyad1), 1);
        Equal("second dyad norm squared", Inner(dyad2, dyad2), 1);
        Equal("two independent dyads are orthogonal", Inner(dyad1, dyad2), 0);
        Equal("first preparation overlap", Inner(dyad1, initial), Q.RootTwo / 4);
        Equal("second preparation overlap", Inner(dyad2, initial), Q.RootTwo / 4);
        ZeroCheck("first dyad exact A eigenoperator", Subtract(AAction(h, dyad1, gamma), Scale(lambda, dyad1)));
        ZeroCheck("second dyad exact A eigenoperator", Subtract(AAction(h, dyad2, gamma), Scale(lambda, dyad2)));
        ZeroCheck("complete residue exact A eigenoperator", Subtract(AAction(h, m, gamma), Scale(lambda, m)));
        Equal("residue trace", Trace(m), 0);
        var norm = Inner(m, m).Rational;
        Equal("complete residue norm squared", norm, new BigRational(1, 4));
        var zz = (-2 * (m[0, 0] + m[6, 6])).Rational;
        var xx = (m[0, 6] + m[6, 0]).Rational;
        Equal("physical endpoint ZZ residue", zz, new BigRational(-1, 2));
        Equal("physical endpoint XX residue", xx, new BigRational(-1, 4));
        var rankOne = Scale(Inner(dyad1, initial), dyad1);
        var rankOneZz = (-2 * (rankOne[0, 0] + rankOne[6, 6])).Rational;
        Equal("one-dyad mutation changes ZZ", rankOneZz, new BigRational(-1, 4));
        var stationary = Project(Outer(v, v));
        var stationaryNorm = Inner(stationary, stationary).Rational;
        Equal("stationary preparation has no slow residue", stationaryNorm, 0);

        // A complex basis of all traceless inputs, plus i times that basis, detects
        // premature Re(), independently of this particular physical preparation.
        int pairCases = 0, pairMismatches = 0;
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
                        pairMismatches += NonzeroCount(Subtract(PairMap(basis, a, b), PhysicalPair(basis, a, b)));
                    }
                }
            }
        Equal("complex-linear map vs physical trace on full traceless basis", pairMismatches, 0);
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

        int oddPairZeros = 0;
        foreach (var (a, b) in new[] { (1, 3), (1, 5), (3, 5) })
            if (NonzeroCount(PairMap(m, a, b)) == 0) oddPairZeros++;
        Equal("three odd-site pair outputs vanish", oddPairZeros, 3);
        // The same algebra beyond the uniform seed: v supported on even sites,
        // w in v^T w=0, and M=w v^T+v w^T C. This is a linear-structure control,
        // not a finite-epsilon eigenvector measurement.
        ZeroCheck("chiral pair-null structure at r=4/3", ChiralNullResidual());

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
        var hPrime = Subtract(Lift(MissingPhaseOnsetWitness.Hopping(7, BigRational.One)), h);
        var kPrime = Scale(-Q.I, hPrime);
        var lambdaPrime = Inner(ePlus, Multiply(kPrime, ePlus));
        Equal("eigenvalue derivative from the left/right pairing", lambdaPrime, -Q.I * Q.RootTwo / 2);
        var bordered = Zero(8, 8);
        var rhs = Zero(8, 1);
        var source = Subtract(Scale(lambdaPrime, ePlus), Multiply(kPrime, ePlus));
        for (int i = 0; i < 7; i++)
        {
            for (int j = 0; j < 7; j++) bordered[i, j] = k[i, j] - (i == j ? lambda : Q.Zero);
            bordered[i, 7] = ePlus[i, 0];
            bordered[7, i] = ePlus[i, 0]; // transpose: K is complex symmetric
            rhs[i, 0] = source[i, 0];
        }
        var solved = Solve(bordered, rhs);
        ZeroCheck("bordered eigenvector derivative residual", Subtract(Multiply(bordered, solved), rhs));
        Equal("bordered compatibility multiplier", solved[7, 0], 0);
        var uPrime = Zero(7, 1);
        for (int i = 0; i < 7; i++) uPrime[i, 0] = solved[i, 0];
        // Differentiate b(r)/sqrt(b(r)^T b(r)) at r=1, directly from its entries.
        var raw = Zero(7, 1);
        var rawPrime = Zero(7, 1);
        int[] entries = { 1, 0, -1, 0, 1, 0, -1 };
        int[] derivatives = { 0, 0, -1, 0, 1, 0, -1 };
        for (int j = 0; j < 7; j++) { raw[j, 0] = entries[j]; rawPrime[j, 0] = derivatives[j]; }
        var vPrime = Subtract(Scale(Q.One / 2, rawPrime), Scale(Inner(raw, rawPrime) / 8, raw));
        var qvPrime = ReflectionEven(vPrime);
        var quPrime = ReflectionEven(uPrime);
        var a0 = Inner(v, b0);
        var overlap = Inner(ePlus, b0);
        var leakage = 2 * a0 * overlap * Inner(qvPrime, quPrime);
        Equal("leakage epsilon^2 coefficient from derivative vectors", leakage,
            -Q.One / 16 + Q.I * Q.RootTwo * (Q)gamma / 16);

        return new(zz, xx, norm, rankOneZz, stationaryNorm, oneCopyNorm, twoCopyNorm,
            allPairSine, decodedNormSquared, leakage.A.Re, leakage.B.Im,
            pairCases, pairMismatches, oddPairZeros, checks.AsReadOnly());
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

    private static Q[,] AAction(Q[,] h, Q[,] x, BigRational gamma)
    {
        var result = Scale(-Q.I, Subtract(Multiply(h, x), Multiply(x, h)));
        for (int i = 0; i < 7; i++)
            for (int j = 0; j < 7; j++)
                if ((i == 3) != (j == 3)) result[i, j] -= 2 * (Q)gamma * x[i, j];
        return result;
    }

    private static Q[,] PairMap(Q[,] m, int a, int b)
    {
        if (!Trace(m).IsZero) throw new ArgumentException("The residue map requires a traceless input.");
        var result = Zero(4, 4);
        var s = (m[a, a] + m[b, b]) / 2;
        result[0, 0] = result[3, 3] = -s;
        result[1, 1] = result[2, 2] = s;
        result[1, 2] = result[2, 1] = (m[a, b] + m[b, a]) / 2;
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

    private static Q[,] ChiralNullResidual()
    {
        var r = (Q)new BigRational(4, 3);
        var v = Zero(7, 1);
        v[0, 0] = 1; v[2, 0] = -r; v[4, 0] = r; v[6, 0] = -r;
        var sum = Zero(4, 4);
        for (int j = 1; j < 7; j++)
        {
            var w = Zero(7, 1);
            w[j, 0] = Q.I;
            w[0, 0] = -Q.I * v[j, 0];
            var m = Zero(7, 7);
            for (int a = 0; a < 7; a++)
                for (int b = 0; b < 7; b++)
                    m[a, b] = w[a, 0] * v[b, 0] + v[a, 0] * w[b, 0] * (b % 2 == 0 ? 1 : -1);
            foreach (var (a, b) in new[] { (1, 3), (1, 5), (3, 5) })
            {
                var image = PairMap(m, a, b);
                // Sum entrywise squared magnitudes, so different failures cannot cancel.
                for (int row = 0; row < 4; row++)
                    for (int col = 0; col < 4; col++) sum[row, col] += image[row, col].Conjugate * image[row, col];
            }
        }
        return sum;
    }
}
