using System.Numerics;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using static RCPsiSquared.Diagnostics.Foundation.RouteBN6Modular;
using static RCPsiSquared.Diagnostics.Foundation.RouteBN6UnfoldingCertificate;

namespace RCPsiSquared.Diagnostics.Foundation;

public sealed record RouteBN6ProfileInvariants(int Alpha, int Beta, int Gamma, int Omega)
{
    public int[] ToArray() => new[] { Alpha, Beta, Gamma, Omega };
    public override string ToString() => $"({Alpha},{Beta},{Gamma},{Omega})";
}

public sealed record RouteBN6Mod101Reading(
    int Prime, int T, int A2Value, int A2Derivative, int ResidualGcdDegree,
    int LambdaCleared, int LambdaPhysical, int Nullity, int SquaredNullity,
    int ProjectorRank, int ComplementCharacteristicAtLambda,
    bool ProjectorIdempotent, bool ProjectorInEvenSector, bool ReducedResolventIdentity,
    RouteBN6ProfileInvariants OneEnd, RouteBN6ProfileInvariants EqualEnds, RouteBN6ProfileInvariants OppositeEnds,
    bool OddFirstOrderZero, bool OneEqualsEqualFirstOrder,
    int ScalarControlOmega, int CommutingControlOmega, bool EvenAdmixtureNonzero);

public sealed record RouteBN6UnfoldingReading(
    int FullDimension, int EvenDimension, int OddDimension, int ResidualDegree, int AtDegree,
    int EvenLoci, int OddLoci, bool ParityTransportExact, bool HermitianRealTExact,
    bool SourcePencilMatches, string CertificateSource, IReadOnlyList<int> CycleType, RouteBN6Mod101Reading Mod101);

/// <summary>F163's exact small-prime readings, rebuilt once per inspection instance.
/// Rebuilds both integer Core parity pencils and the full 90D XY block, joins the
/// strict 266-locus inventory, then recomputes irreducibility modulo 367 and the
/// effective end-profile matrices modulo 101. The complete bounded-CRT factorization
/// is an explicitly consumed, fingerprinted certificate input. These modular residues
/// are not real response coefficients or an experimental operating point.</summary>
public sealed class RouteBN6A2UnfoldingWitness : IInspectable
{
    private readonly Lazy<RouteBN6UnfoldingReading> _reading;
    private int _reconstructionCount;

    public RouteBN6A2UnfoldingWitness(string? certificatePath = null) =>
        _reading = new(() => Reconstruct(certificatePath), LazyThreadSafetyMode.ExecutionAndPublication);

    public RouteBN6UnfoldingReading Reading => _reading.Value;
    public int ReconstructionCount => Volatile.Read(ref _reconstructionCount);

    /// <summary>Real integer polynomial reduction uses F_p directly: 367 need not split in Z[i].
    /// A degree drop or nonsquarefree reduction is unusable, never evidence for irreducibility.</summary>
    public int[]? CycleTypeModulo367(BigInteger[] coefficients)
    {
        ArgumentNullException.ThrowIfNull(coefficients);
        if (coefficients.Length < 2 || coefficients[^1].IsZero)
            throw new ArgumentException("A nonconstant polynomial with nonzero leading coefficient is required.", nameof(coefficients));
        var reduced = coefficients.Select(c => (int)((c % 367 + 367) % 367)).ToArray();
        return reduced[^1] == 0 ? null : OcticGaloisCertificate.CycleTypeOfFpPoly(reduced, 367);
    }

    private RouteBN6UnfoldingReading Reconstruct(string? certificatePath)
    {
        Interlocked.Increment(ref _reconstructionCount);
        var even = FoldResultantCertificate.ExportRouteBN6ExactPencil(false);
        var odd = FoldResultantCertificate.ExportRouteBN6ExactPencil(true);
        var inventory = RouteBA2N6Inventory.LoadDefault();
        var certificate = Load(certificatePath, even, odd, inventory);
        bool transport = ParityTransport(even.ResidualInT, odd.ResidualInT)
            && ParityTransport(even.AtFactorInT, odd.AtFactorInT);
        Require(transport, "F163 exact t -> -t parity transport failed.");
        int[]? cycle = CycleTypeModulo367(certificate.A2);
        Require(cycle is not null && cycle.SequenceEqual(new[] { 133 }), "F163 A2 reduction is not irreducible of degree 133.");

        // All inputs are small integers (q=0 or -i, bond weights 0/1, Delta=0).
        // Complex arithmetic here is exact integer arithmetic; reject any noninteger
        // output before reducing it. No eigensolver or float rank participates.
        var d = IntegerMatrix(XxzCoherenceBlock.BuildFull(6, Complex.Zero, 0));
        var realTOne = IntegerMatrix(XxzCoherenceBlock.BuildFull(6, -Complex.ImaginaryOne, 0));
        var t = Add(realTOne, d, -1);
        int[,] Bond(int index)
        {
            var weights = new double[5];
            weights[index] = 1;
            return Add(IntegerMatrix(XxzCoherenceBlock.BuildFull(6, -Complex.ImaginaryOne, 0, weights)), d, -1);
        }
        var reflection = IntegerMatrix(XxzCoherenceBlock.ReflectionPermutation(6).ToArray());
        bool hermitian = Symmetric(d) && Symmetric(realTOne);
        Require(hermitian, "F163 real-t pencil is not symmetric.");
        var mod101 = ReadMod101(certificate.A2, even, odd, d, t, Bond(0), Bond(4), reflection);
        return new(d.GetLength(0), even.SectorDimension, odd.SectorDimension,
            even.ResidualLambdaDegree, even.AtFactorInT.Length - 1,
            inventory.Loci.Count(l => l.Parity == A2Parity.Even), inventory.Loci.Count(l => l.Parity == A2Parity.Odd),
            transport, hermitian, certificate.SourcePencilMatches, certificate.Description,
            Array.AsReadOnly((int[])cycle!.Clone()), mod101);
    }

    private static RouteBN6Mod101Reading ReadMod101(BigInteger[] a2, RouteBN6ExactPencil even, RouteBN6ExactPencil odd,
        int[,] d, int[,] hopping, int[,] left, int[,] right, int[,] reflection)
    {
        const int t = 31;
        int inv2 = (int)CrossFormCertificate.Inv(2, Prime);
        int[] a2Mod = a2.Select(c => (int)((c % Prime + Prime) % Prime)).ToArray();
        int a2Value = Evaluate(a2Mod, t), derivative = Evaluate(Derivative(a2Mod), t);
        Require(a2Value == 0 && derivative != 0, "F163 modular specialization is not a simple A2 root.");
        int[] residual = EvaluatePencil(even.ResidualInT, t);
        int[] gcd = Gcd(residual, Derivative(residual));
        Require(residual.Length == 33 && gcd.Length == 2, "F163 residual repeated root is not unique.");
        int cleared = Mod(-gcd[0]), lambda = Mod(cleared * inv2);
        int[,] l = Add(d, hopping, t), identity = Identity(d.GetLength(0));
        int[,] shifted = Add(l, identity, -lambda);
        int[] characteristic = Characteristic(l);
        // Independent full-matrix characteristic vs the Core sector factors, with
        // Lambda=2*lambda explicitly converted back to the physical eigenvalue.
        int[] SectorPolynomial(RouteBN6ExactPencil pencil)
        {
            var coefficients = MultiplyPolynomials(EvaluatePencil(pencil.ResidualInT, t), EvaluatePencil(pencil.AtFactorInT, t));
            int n = coefficients.Length - 1;
            return coefficients.Select((c, k) => Mod(c * CrossFormCertificate.PowMod(inv2, n - k, Prime))).ToArray();
        }
        Require(characteristic.SequenceEqual(MultiplyPolynomials(SectorPolynomial(even), SectorPolynomial(odd))),
            "F163 modular full block differs from the rebuilt Core sector pencils.");
        var linear = new[] { Mod(-lambda), 1 };
        var (quotient, remainder) = Divide(characteristic, MultiplyPolynomials(linear, linear));
        int complementAtLambda = Evaluate(quotient, lambda);
        Require(remainder.All(c => c == 0) && complementAtLambda != 0, "F163 full multiplicity is not exactly two.");
        var projector = Scale(MatrixPolynomial(quotient, l), (int)CrossFormCertificate.Inv(complementAtLambda, Prime));
        bool idempotent = Equal(Multiply(projector, projector), projector);
        bool evenPlane = Equal(Multiply(reflection, projector), projector);
        Require(idempotent && evenPlane && IsZero(Multiply(shifted, projector)), "F163 modular projector identity failed.");
        var complement = Add(identity, projector, -1);
        var lambdaMinusL = Add(Scale(identity, lambda), l, -1);
        var resolvent = Multiply(Multiply(complement, Inverse(Add(lambdaMinusL, projector))), complement);
        bool resolventIdentity = Equal(Multiply(Scale(shifted, -1), resolvent), complement);
        Require(resolventIdentity, "F163 reduced resolvent identity failed.");
        var a = Multiply(Multiply(projector, hopping), projector);
        int[,] First(int[,] direction) => Multiply(Multiply(projector, Scale(direction, t)), projector);
        var one = First(left);
        var equal = First(Scale(Add(left, right), inv2));
        var oddDirection = Scale(Add(left, right, -1), inv2);
        var oppositeFirst = First(oddDirection);
        var oddStep = Scale(oddDirection, t);
        var oppositeSecond = Multiply(Multiply(Multiply(Multiply(projector, oddStep), resolvent), oddStep), projector);
        var oneInvariants = Invariants(a, one);
        var equalInvariants = Invariants(a, equal);
        var oppositeInvariants = Invariants(a, oppositeSecond);
        Require(oneInvariants.Alpha != 0 && oneInvariants.Omega != 0 && equalInvariants.Omega != 0 && oppositeInvariants.Omega != 0,
            "F163 effective response fails the nonzero unfolding discriminant condition.");
        return new(Prime, t, a2Value, derivative, gcd.Length - 1, cleared, lambda,
            d.GetLength(0) - Rank(shifted), d.GetLength(0) - Rank(Multiply(shifted, shifted)), Rank(projector),
            complementAtLambda, idempotent, evenPlane, resolventIdentity,
            oneInvariants, equalInvariants, oppositeInvariants,
            IsZero(oppositeFirst), Equal(one, equal), Invariants(a, projector).Omega, Invariants(a, a).Omega,
            !IsZero(First(Add(oddDirection, Scale(Add(left, right), inv2)))));
    }

    private static int[] EvaluatePencil(BigInteger[][] rows, int t) => rows.Select(row =>
        Evaluate(row.Select(c => (int)((c % Prime + Prime) % Prime)).ToArray(), t)).ToArray();

    private static int[,] IntegerMatrix(Complex[,] source)
    {
        var result = new int[source.GetLength(0), source.GetLength(1)];
        for (int i = 0; i < result.GetLength(0); i++)
        for (int j = 0; j < result.GetLength(1); j++)
        {
            Complex z = source[i, j];
            Require(z.Imaginary == 0 && double.IsFinite(z.Real) && z.Real == Math.Truncate(z.Real),
                "F163 builder left exact integer arithmetic.");
            result[i, j] = checked((int)z.Real);
        }
        return result;
    }

    private static bool Symmetric(int[,] matrix)
    {
        for (int i = 0; i < matrix.GetLength(0); i++)
        for (int j = 0; j < i; j++)
            if (matrix[i, j] != matrix[j, i]) return false;
        return true;
    }

    public string DisplayName => "F163: N=6 Route-B A2 exact local unfolding";
    public string Summary
    {
        get
        {
            var r = Reading;
            return $"N=6 complex-q geometry: {r.FullDimension}={r.EvenDimension}+{r.OddDimension}, " +
                   $"{r.EvenLoci}+{r.OddLoci}=266 loci; local epsilon orders one/equal/opposite = 1/1/2. " +
                   "Exact pencil and modular readings live; full bounded-CRT factorization consumed. " +
                   "Not an all-N or physical measurement claim.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var r = Reading;
            var m = r.Mod101;
            yield return new InspectableNode("source and proof boundary", summary: r.CertificateSource);
            yield return new InspectableNode("exact pencil and inventory", summary:
                $"{r.FullDimension}={r.EvenDimension}+{r.OddDimension}; residual/AT degrees {r.ResidualDegree}+{r.AtDegree}; " +
                $"{r.EvenLoci}+{r.OddLoci} loci. Exact t->-t parity transport: {r.ParityTransportExact}; real-t symmetry: {r.HermitianRealTExact}.",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode("live irreducibility modulo 367", summary:
                $"Frobenius cycle [{string.Join(",", r.CycleType)}], degree-preserving real integer reduction. " +
                "Odd degree gives a real Hermitian embedding; conjugacy transports the certified semisimple double plane.",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode("live modular plane", summary:
                $"p={m.Prime}, t={m.T}, A2(t)={m.A2Value}, A2'(t)={m.A2Derivative}; Lambda={m.LambdaCleared}, lambda={m.LambdaPhysical}; " +
                $"nullities of L-lambda I and its square {m.Nullity}/{m.SquaredNullity}, projector rank {m.ProjectorRank}; " +
                $"g(lambda)={m.ComplementCharacteristicAtLambda}, reduced resolvent identity: {m.ReducedResolventIdentity}.",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode("live local response invariants", summary:
                $"(alpha,beta,gamma,Omega) mod101: one {m.OneEnd}, equal {m.EqualEnds}, opposite {m.OppositeEnds}. " +
                $"P V_- P=0: {m.OddFirstOrderZero}; one/equal first-order matrices equal: {m.OneEqualsEqualFirstOrder}. " +
                "Residues are not real response coefficients; the zero gamma residue makes no characteristic-zero zero claim.",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode("controls through the effective-matrix reading", summary:
                $"Scalar Omega={m.ScalarControlOmega}, commuting Omega={m.CommutingControlOmega}; " +
                $"even admixture reaches the plane: {m.EvenAdmixtureNonzero}. " +
                "F131 fixes parity; nonzero alpha and Omega establish the leading local orders.",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode("end-bond epsilon convention", summary:
                "L(t,epsilon)=D+t(T+epsilon V), t=i*qCSharp, Lambda=2*lambda. Interior weights stay one; " +
                "the (left,right) weights for one/equal/opposite are (1+epsilon,1), " +
                "(1+epsilon/2,1+epsilon/2), and (1+epsilon/2,1-epsilon/2). " +
                "Their two local EP2 locations move at orders epsilon/epsilon/epsilon^2.");
            yield return new InspectableNode("partner and scope", summary:
                "F89d transports the certified character to the (1,4) partner at (conj(q),-conj(lambda)-12); " +
                "this is the same 266-locus layer transported. N=6, uniform XY and Z dephasing, gamma=1, Delta=0; " +
                "complex-q analytic continuation and sufficiently small local epsilon only.");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
