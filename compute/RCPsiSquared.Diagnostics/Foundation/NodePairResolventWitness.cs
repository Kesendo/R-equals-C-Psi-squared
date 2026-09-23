using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Live exact witness for the node-pair resolvent theorem and bond criterion. It
/// recomputes a reduced resolvent by an exact bordered solve, detuned-block characteristic
/// polynomials, blind/root counts by polynomial gcd, and rational nonuniform full-chain
/// nullvectors. No floating tolerance or Python artifact is consumed.
///
/// <para>The general node-pair theorem and sufficient bond direction are not restricted to the
/// uniform centre-watched family. For any zero-free Jacobi path with fixed diagonal and
/// r not in {0,+1,-1}, the same baseline blind energy survives exactly when its original
/// eigenvector has a node at either moved-bond endpoint. The uniform-centre family additionally
/// has total-count equality; off-centre arms may acquire new blind energies.</para>
///
/// <para>Root: <c>inspect --root nodepair</c>. Proof:
/// <c>docs/proofs/PROOF_NODE_PAIR_RESOLVENT.md</c>.</para></summary>
public sealed class NodePairResolventWitness : IInspectable
{
    private readonly Lazy<NodePairResolventReading> _reading = new(Reconstruct);

    public NodePairResolventReading Reading => _reading.Value;

    public string DisplayName => "NodePairResolventWitness (exact canonical live lab)";

    public string Summary
    {
        get
        {
            var reading = Reading;
            int iffMatches = reading.CriterionCases.Count(item =>
                item.DetunedBlindCount == item.EndpointNodeCount);
            return $"exact N=7 reduced resolvent: node/node {reading.NodePairEntry}, " +
                   $"node/non-node {reading.NodeNonNodeEntry}; direct factorization " +
                   $"{(reading.FactorizationMatches ? "MATCH" : "MISMATCH")}; uniform-centre iff " +
                   $"{iffMatches}/{reading.CriterionCases.Count} canonical cases outside r=0,+1,-1; " +
                   $"off-centre N=6 new roots {reading.OffCentreBirth.BaselineBlindCount}->" +
                   $"{reading.OffCentreBirth.BornBlindCount}; nonuniform fixed-energy cases " +
                   $"{reading.NonuniformFixedEnergyCases.Count}.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var reading = Reading;
            yield return new InspectableNode(
                "node-pair reduced-resolvent zero",
                $"N=7, E=0, nodes x=1 and y=3: R[1,3]={reading.NodePairEntry}; " +
                "exact bordered solve on 0-H",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                "node/non-node nonzero control",
                $"same solve, node x=1 and non-node y=2: R[1,2]={reading.NodeNonNodeEntry} " +
                "(must be -1/8, not zero)",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                "direct determinant vs factorized polynomial",
                $"m=4, b=1, r=2: direct χ_L={reading.DirectPolynomial}; factorized=" +
                $"{reading.FactorizedPolynomial}; " +
                (reading.FactorizationMatches ? "MATCH" : "MISMATCH"),
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                "uniform-centre iff cases",
                string.Join("; ", reading.CriterionCases.Select(item =>
                    $"m={item.HalfLength} b={item.Bond} r={item.R}: " +
                    $"blind={item.DetunedBlindCount}, endpoint-nodes={item.EndpointNodeCount}")),
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                "off-centre birth at N=6",
                $"c=2, b=0, q=r^2=2: chi_L={reading.OffCentreBirth.MovedLeftPolynomial}, " +
                $"chi_R={reading.OffCentreBirth.UntouchedRightPolynomial}, " +
                $"gcd={reading.OffCentreBirth.SharedFactor}; blind baseline=" +
                $"{reading.OffCentreBirth.BaselineBlindCount}, q=2=" +
                $"{reading.OffCentreBirth.BornBlindCount}, changed opposite arm=" +
                $"{reading.OffCentreBirth.ChangedOppositeArmBlindCount}",
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                "nonuniform fixed-energy cases",
                string.Join("; ", reading.NonuniformFixedEnergyCases.Select(item =>
                    $"E={item.Energy}: baseline blind={item.BaselineBlind}, " +
                    $"endpoint node={item.HasEndpointNode}, same E blind at r=2=" +
                    item.SameEnergyBlindAfterMove)),
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                "exceptional set r=0,+1,-1",
                string.Join("; ", reading.ExceptionalCases.Select(item =>
                    $"r={item.R}: blind={item.DetunedBlindCount}, nodes={item.EndpointNodeCount}, " +
                    item.Fence)),
                provenance: NodeProvenance.Live);
            yield return new InspectableNode(
                "scope and proof",
                "PROOF_NODE_PAIR_RESOLVENT Theorem 1 is symmetry-free; Corollary C is the " +
                "same-baseline-energy iff for a zero-free Jacobi path with fixed diagonal. " +
                "Corollary D gives uniform-centre total-count equality. Arbitrary-arm new-root " +
                "classification and rate claims are not made. " +
                "Run `inspect --root nodepair`.");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    private static NodePairResolventReading Reconstruct()
    {
        var reduced = CanonicalReducedResolvent();
        var rTwo = new BigRational(2);
        var direct = DetunedPathPolynomial(4, 1, rTwo);
        var factorized = FactorizedPathPolynomial(4, 1, rTwo);

        var criterion = new[]
        {
            CriterionCase(2, 0, new BigRational(2)),
            CriterionCase(3, 0, new BigRational(2)),
            CriterionCase(4, 1, new BigRational(2, 3)),
            CriterionCase(5, 1, new BigRational(-2)),
        };
        var exceptional = new[]
        {
            ExceptionalCase(3, 0, -BigRational.One, "sign-gauge return: r^2=1"),
            ExceptionalCase(3, 0, BigRational.Zero, "zero bond: the Jacobi chain is cut"),
            ExceptionalCase(3, 0, BigRational.One, "unperturbed return: r^2=1"),
        };

        return new NodePairResolventReading(
            reduced[1, 3], reduced[1, 2], direct, factorized,
            criterion, exceptional, ReconstructOffCentreBirth(),
            ReconstructNonuniformFixedEnergyCases());
    }

    private static NodePairCriterionCase CriterionCase(int halfLength, int bond, BigRational r) =>
        new(halfLength, bond, r, EndpointNodeCount(halfLength, bond),
            DetunedBlindCount(halfLength, bond, r));

    private static NodePairExceptionalCase ExceptionalCase(
        int halfLength, int bond, BigRational r, string fence) =>
        new(halfLength, bond, r, EndpointNodeCount(halfLength, bond),
            DetunedBlindCount(halfLength, bond, r), InTheoremDomain: false, fence);

    /// <summary>Exact reduced resolvent of the N=7 uniform hopping-2 path at E=0. The upper-left
    /// block of the bordered inverse [−H v; v^T 0]^{-1} is the inverse on v^⊥.</summary>
    private static BigRational[,] CanonicalReducedResolvent()
    {
        const int n = 7;
        var bordered = ZeroMatrix(n + 1, n + 1);
        for (int site = 0; site < n - 1; site++)
        {
            bordered[site, site + 1] = -2;
            bordered[site + 1, site] = -2;
        }
        int[] kernel = { 1, 0, -1, 0, 1, 0, -1 };
        for (int site = 0; site < n; site++)
        {
            bordered[site, n] = kernel[site];
            bordered[n, site] = kernel[site];
        }

        var rhs = ZeroMatrix(n + 1, n);
        for (int site = 0; site < n; site++) rhs[site, site] = BigRational.One;
        var solved = BigRationalLinearAlgebra.Solve(bordered, rhs);
        var reduced = ZeroMatrix(n, n);
        for (int row = 0; row < n; row++)
            for (int column = 0; column < n; column++)
                reduced[row, column] = solved[row, column];
        return reduced;
    }

    private static RationalPolynomial DetunedPathPolynomial(
        int halfLength, int bond, BigRational r)
    {
        if (halfLength < 2) throw new ArgumentOutOfRangeException(nameof(halfLength));
        if (bond < 0 || bond >= halfLength - 1) throw new ArgumentOutOfRangeException(nameof(bond));
        var h = ZeroMatrix(halfLength, halfLength);
        for (int site = 0; site < halfLength - 1; site++)
        {
            BigRational hopping = site == bond ? 2 * r : new BigRational(2);
            h[site, site + 1] = hopping;
            h[site + 1, site] = hopping;
        }
        return new RationalPolynomial(BigRationalMatrixCharpoly.Characteristic(h));
    }

    private static RationalPolynomial FactorizedPathPolynomial(
        int halfLength, int bond, BigRational r) =>
        PathPolynomial(halfLength)
        - (new BigRational(4) * (r * r - BigRational.One))
          * (PathPolynomial(bond) * PathPolynomial(halfLength - bond - 2));

    private static int DetunedBlindCount(int halfLength, int bond, BigRational r) =>
        PolynomialGcd(DetunedPathPolynomial(halfLength, bond, r),
            PathPolynomial(halfLength)).Degree;

    private static int EndpointNodeCount(int halfLength, int bond) =>
        PolynomialGcd(PathPolynomial(halfLength), PathPolynomial(bond)).Degree
        + PolynomialGcd(PathPolynomial(halfLength),
            PathPolynomial(halfLength - bond - 2)).Degree;

    private static RationalPolynomial PathPolynomial(int length)
    {
        if (length < 0) throw new ArgumentOutOfRangeException(nameof(length));
        var p0 = RationalPolynomial.One;
        if (length == 0) return p0;
        var x = RationalPolynomial.Monomial(1, BigRational.One);
        var p1 = x;
        for (int size = 2; size <= length; size++)
            (p0, p1) = (p1, x * p1 - new BigRational(4) * p0);
        return p1;
    }

    private static RationalPolynomial PolynomialGcd(
        RationalPolynomial left, RationalPolynomial right)
    {
        while (!right.IsZero)
            (left, right) = (right, left.Mod(right));
        if (left.IsZero) return RationalPolynomial.Zero;
        return (BigRational.One / left[left.Degree]) * left;
    }

    private static NodePairBirthCase ReconstructOffCentreBirth()
    {
        // N=6, c=2, b=0. The moved hopping is 2*sqrt(q), so its SQUARE is rational at q=2.
        // The continuant reads q exactly without introducing an irrational matrix entry.
        BigRational z = BigRational.Zero;
        var leftBaseline = JacobiPolynomial(new[] { z, z }, new BigRational[] { 4 });
        var leftBorn = JacobiPolynomial(new[] { z, z }, new BigRational[] { 8 });
        var right = JacobiPolynomial(new[] { z, z, z }, new BigRational[] { 4, 4 });
        var changedRight = JacobiPolynomial(new[] { z, z, z }, new BigRational[] { 4, 9 });
        var shared = PolynomialGcd(leftBorn, right);
        return new NodePairBirthCase(leftBorn, right, shared,
            PolynomialGcd(leftBaseline, right).Degree, shared.Degree,
            PolynomialGcd(leftBorn, changedRight).Degree);
    }

    private static RationalPolynomial JacobiPolynomial(
        IReadOnlyList<BigRational> diagonal, IReadOnlyList<BigRational> squaredHoppings)
    {
        if (squaredHoppings.Count != Math.Max(0, diagonal.Count - 1))
            throw new ArgumentException("One squared hopping is required between adjacent sites.");
        var x = RationalPolynomial.Monomial(1, BigRational.One);
        var previous = RationalPolynomial.One;
        if (diagonal.Count == 0) return previous;
        var current = x - new RationalPolynomial(diagonal[0]);
        for (int site = 1; site < diagonal.Count; site++)
            (previous, current) = (current,
                (x - new RationalPolynomial(diagonal[site])) * current
                - squaredHoppings[site - 1] * previous);
        return current;
    }

    private static IReadOnlyList<NodePairFixedEnergyCase> ReconstructNonuniformFixedEnergyCases()
    {
        BigRational[] diagonal = { 1, 0, 1, 7, 1, 0, 1 };
        BigRational[] hoppings = { 2, 4, 5, 6, 4, 2 };
        var baseline = JacobiMatrix(diagonal, hoppings);
        hoppings[0] = 4; // r=2 on bond 0; the diagonal and all other hoppings stay fixed.
        var moved = JacobiMatrix(diagonal, hoppings);
        return new[] { new BigRational(-4), BigRational.One, new BigRational(5) }
            .Select(energy =>
            {
                var original = BigRationalLinearAlgebra.Nullspace(Shifted(baseline, energy));
                var after = BigRationalLinearAlgebra.Nullspace(Shifted(moved, energy));
                return new NodePairFixedEnergyCase(energy,
                    original.Count == 1 && original[0][3].IsZero,
                    original.Count == 1 && (original[0][0].IsZero || original[0][1].IsZero),
                    after.Count == 1 && after[0][3].IsZero);
            }).ToArray();
    }

    private static BigRational[,] JacobiMatrix(
        IReadOnlyList<BigRational> diagonal, IReadOnlyList<BigRational> hoppings)
    {
        if (hoppings.Count != diagonal.Count - 1)
            throw new ArgumentException("One hopping is required between adjacent sites.");
        var matrix = ZeroMatrix(diagonal.Count, diagonal.Count);
        for (int site = 0; site < diagonal.Count; site++)
        {
            matrix[site, site] = diagonal[site];
            if (site == hoppings.Count) continue;
            matrix[site, site + 1] = hoppings[site];
            matrix[site + 1, site] = hoppings[site];
        }
        return matrix;
    }

    private static BigRational[,] Shifted(BigRational[,] matrix, BigRational energy)
    {
        int n = matrix.GetLength(0);
        var shifted = (BigRational[,])matrix.Clone();
        for (int site = 0; site < n; site++) shifted[site, site] -= energy;
        return shifted;
    }

    private static BigRational[,] ZeroMatrix(int rows, int columns)
    {
        var matrix = new BigRational[rows, columns];
        for (int row = 0; row < rows; row++)
            for (int column = 0; column < columns; column++)
                matrix[row, column] = BigRational.Zero;
        return matrix;
    }
}

public sealed record NodePairResolventReading(
    BigRational NodePairEntry,
    BigRational NodeNonNodeEntry,
    RationalPolynomial DirectPolynomial,
    RationalPolynomial FactorizedPolynomial,
    IReadOnlyList<NodePairCriterionCase> CriterionCases,
    IReadOnlyList<NodePairExceptionalCase> ExceptionalCases,
    NodePairBirthCase OffCentreBirth,
    IReadOnlyList<NodePairFixedEnergyCase> NonuniformFixedEnergyCases)
{
    public bool FactorizationMatches => DirectPolynomial.Equals(FactorizedPolynomial);
}

public sealed record NodePairCriterionCase(
    int HalfLength,
    int Bond,
    BigRational R,
    int EndpointNodeCount,
    int DetunedBlindCount);

public sealed record NodePairExceptionalCase(
    int HalfLength,
    int Bond,
    BigRational R,
    int EndpointNodeCount,
    int DetunedBlindCount,
    bool InTheoremDomain,
    string Fence);

public sealed record NodePairBirthCase(
    RationalPolynomial MovedLeftPolynomial,
    RationalPolynomial UntouchedRightPolynomial,
    RationalPolynomial SharedFactor,
    int BaselineBlindCount,
    int BornBlindCount,
    int ChangedOppositeArmBlindCount);

public sealed record NodePairFixedEnergyCase(
    BigRational Energy,
    bool BaselineBlind,
    bool HasEndpointNode,
    bool SameEnergyBlindAfterMove);
