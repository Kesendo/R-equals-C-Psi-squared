using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The Absorption Theorem with its two mathematical objects kept distinct.
/// The Z-dephasing dissipator is diagonal on a computational-basis pair
/// <c>|A&gt;&lt;B|</c> and contributes the cell cost <c>2γ n_diff(A,B)</c>; for an
/// initially isolated cell this is its initial decay slope. A Liouvillian eigenmode is
/// generally a superposition of such cells. Its decay rate follows from the Hermitian
/// part of L and is <c>-Re λ = 2γ &lt;n_XY&gt;_v</c>, where the expectation may be
/// non-integer. Thus the coefficient 2 is exact, but the interacting eigenvalue spectrum
/// is not thereby quantized in steps of 2γ.</summary>
public sealed class AbsorptionTheoremClaim : Claim
{
    public Pi2DyadicLadderClaim Ladder { get; }

    /// <summary>The exact coefficient 2 multiplying both the basis-pair dissipator
    /// cost and the eigenmode expectation value.</summary>
    public double DissipatorCoefficient => Ladder.Term(0);

    /// <summary>Dissipator cost of one bra-ket disagreement: <c>2γ</c>. This is a
    /// basis-cell quantity, not the smallest nonzero eigenmode decay-rate step.</summary>
    public double SingleDisagreementCellCost(double gammaZero)
    {
        ValidateGamma(gammaZero);
        return DissipatorCoefficient * gammaZero;
    }

    /// <summary>Diagonal dissipator cost / isolated-cell initial decay slope of
    /// <c>|A&gt;&lt;B|</c>: <c>2γ n_diff(A,B)</c>.</summary>
    public double BasisPairDissipatorCost(int nDiff, double gammaZero)
    {
        if (nDiff < 0)
            throw new ArgumentOutOfRangeException(nameof(nDiff), nDiff, "n_diff must be >= 0.");
        ValidateGamma(gammaZero);
        return DissipatorCoefficient * gammaZero * nDiff;
    }

    /// <summary>Decay rate of a right Liouvillian eigenmode from its normalized light
    /// expectation: <c>-Re λ = 2γ &lt;n_XY&gt;_v</c>. The expectation is continuous.</summary>
    public double EigenmodeDecayRate(double averageNXy, double gammaZero)
    {
        if (!double.IsFinite(averageNXy) || averageNXy < 0)
            throw new ArgumentOutOfRangeException(nameof(averageNXy), averageNXy,
                "average n_XY must be finite and >= 0.");
        ValidateGamma(gammaZero);
        return DissipatorCoefficient * gammaZero * averageNXy;
    }

    /// <summary>Upper bound <c>2γN</c> following from <c>0 &lt;= &lt;n_XY&gt;_v &lt;= N</c>.
    /// This is a ceiling, not evidence that eigenvalues fill a quantized grid.</summary>
    public double EigenmodeDecayRateCeiling(int n, double gammaZero)
    {
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), n, "N must be >= 1.");
        ValidateGamma(gammaZero);
        return DissipatorCoefficient * gammaZero * n;
    }

    /// <summary>Recover the continuous eigenmode light expectation from its decay rate.</summary>
    public double AverageNXyFromEigenmodeDecayRate(double rate, double gammaZero)
    {
        if (!double.IsFinite(rate) || rate < 0)
            throw new ArgumentOutOfRangeException(nameof(rate), rate, "rate must be finite and >= 0.");
        if (!double.IsFinite(gammaZero) || gammaZero <= 0)
            throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "gamma must be finite and > 0.");
        return rate / (DissipatorCoefficient * gammaZero);
    }

    /// <summary>Basis-pair cell-cost sum under the bra Hamming complement:
    /// <c>2γ n_diff + 2γ(N-n_diff) = 2γN</c>.</summary>
    public double HammingComplementCellCostSum(int blockSize, double gammaZero)
    {
        if (blockSize < 1)
            throw new ArgumentOutOfRangeException(nameof(blockSize), blockSize, "block size must be >= 1.");
        ValidateGamma(gammaZero);
        return DissipatorCoefficient * gammaZero * blockSize;
    }

    public bool DissipatorCoefficientMatchesLiteral() =>
        Math.Abs(DissipatorCoefficient - 2.0) < 1e-15;

    public AbsorptionTheoremClaim(Pi2DyadicLadderClaim ladder)
        : base("Absorption Theorem: basis-pair dissipator cost 2γ n_diff; eigenmode decay -Re(λ)=2γ<n_XY>_v; coefficient 2=a_0",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_ABSORPTION_THEOREM.md + " +
               "docs/proofs/PROOF_F1_NONUNIFORM_GAMMA.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        Ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
    }

    public override string DisplayName =>
        "Absorption Theorem: basis-pair cell cost and eigenmode light expectation";

    public override string Summary =>
        $"Basis-pair dissipator diagonal/isolated initial slope is 2γ n_diff; eigenmode decay is " +
        $"-Re(λ)=2γ<n_XY>_v with generally non-integer expectation; coefficient 2=a_0; " +
        $"2γN is a ceiling, not an eigenvalue step ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("dissipator coefficient (= a_0 = 2)", DissipatorCoefficient);
            yield return new InspectableNode("basis-pair cell reading",
                summary: "D[|A><B|] contributes -2γ n_diff(A,B)|A><B|. It is the diagonal cell cost and isolated-cell initial slope, not generally an eigenvalue.");
            yield return new InspectableNode("eigenmode reading",
                summary: "For a normalized Liouvillian eigenmode v, -Re λ=2γ<n_XY>_v from the Hermitian-part Rayleigh quotient; Hamiltonian mixing permits non-integer expectations.");
            yield return new InspectableNode("Hamming-complement cell-cost sum",
                summary: $"n_diff maps to N-n_diff, so the two basis-cell costs sum to 2γN; at N=3, γ=1 the sum is {HammingComplementCellCostSum(3, 1.0):G6}.");
            yield return new InspectableNode("Pi2 anchor drift check",
                summary: $"DissipatorCoefficientMatchesLiteral={DissipatorCoefficientMatchesLiteral()} (a_0={DissipatorCoefficient}).");
        }
    }

    private static void ValidateGamma(double gammaZero)
    {
        if (!double.IsFinite(gammaZero) || gammaZero < 0)
            throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero,
                "gamma must be finite and >= 0.");
    }
}
