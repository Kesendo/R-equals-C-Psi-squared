using System.Globalization;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for the β-exotic exclusion at N = 5 (Anchor: <c>BetaExoticExcludedAtN5Claim</c>,
/// <c>inspect --root betaexotic</c>; <c>experiments/F89_SEED_EXISTENCE_REDUCTION.md</c>, section "The
/// β-exotic is excluded at N = 5 only"). It re-runs the certificate at inspect time, both R-parity
/// sectors, and reads one integer off each run: the maximum root multiplicity of the discriminant
/// D(q) = disc_Λ(F_res)(q) away from q = 0.
///
/// <para><b>Why that integer is the whole statement.</b> The order of a disc zero at a branch locus
/// reads the Puiseux exponent of the branches meeting there: 1 for a defective EP2 (exponent ½, the
/// seeds), 2 for a diabolic crossing, 2 for a cubic branch point, and 3 for the β-exotic (exponent
/// 3/2, normal form β(s) = [[0, s], [s², 0]]). Every OTHER colliding pair at the same locus contributes
/// a non-negative order, so a coincident collision can only raise the total, never mask it. Maximum
/// multiplicity 2 therefore excludes the β-exotic outright: a demand deliberately weaker than
/// squarefreeness, since the diabolic loci are genuine double roots and must survive.</para>
///
/// <para><b>What "certified" means here.</b> <see cref="FoldResultantCertificate"/> reads the layers
/// modulo one prime. The lift to ℚ(i) is one-way (reduction mod π_p ∤ lc_q(D) is a degree-preserving
/// homomorphism, so roots only merge and the mod-p maximum bounds the true one from above), and it
/// requires the prime to be good at BOTH ends of the q-axis: attaining the true deg_q D (no root
/// escaped to q = ∞) and the true q-valuation (no nonzero root collapsed onto q = 0, where the
/// q-power strip would silently discard it). Each true value is certified by its own Hadamard bound;
/// a prime attaining both is searched for. <c>DiscLayersCertified</c> is exactly that guarantee.
/// A witness that reported <c>MaxDiscMultiplicity</c> without it would be reading a diagnostic, not a
/// proof, which is what the number was before 2026-07-09.</para>
///
/// <para><b>Scope carried from the claim.</b> Per-N, not all-N: this retires N = 5. It does not
/// exclude a cubic branch point (ord disc = 2, hiding in the multiplicity-2 layer); that is ruled out
/// at a count-drop by elimination: a cubic point keeps one real branch and one conjugate pair on both
/// sides of q*, so it cannot change the real count, while an EP2 can. The all-N item (s₆ ≠ 0 at every
/// forced seed) is untouched.</para></summary>
public sealed class BetaExoticExclusionWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>The only N the certificate has been run at (n = 7 is the same call, unrun).</summary>
    public const int CertifiedN = 5;

    /// <summary>A β-exotic (Puiseux exponent 3/2) forces a disc zero of this order.</summary>
    public const int BetaExoticDiscOrder = 3;

    public int N { get; }

    private readonly FoldResultantCertificate.CompleteReport _rEven;
    private readonly FoldResultantCertificate.CompleteReport _rOdd;

    public BetaExoticExclusionWitness(int n = CertifiedN)
    {
        if (n != CertifiedN)
            throw new ArgumentOutOfRangeException(nameof(n),
                $"the β-exotic certificate has only been run at N = {CertifiedN} (n = 7 is the same call at " +
                $"residual degree 53, not yet run); got {n}");
        N = n;
        _rEven = FoldResultantCertificate.CertifyComplete(n, rOdd: false);
        _rOdd = FoldResultantCertificate.CertifyComplete(n, rOdd: true);
    }

    /// <summary>True when BOTH parity sectors are certified and carry no disc root of multiplicity ≥ 3.</summary>
    public bool BetaExoticExcluded =>
        _rEven.DiscLayersCertified && _rOdd.DiscLayersCertified &&
        _rEven.MaxDiscMultiplicity < BetaExoticDiscOrder &&
        _rOdd.MaxDiscMultiplicity < BetaExoticDiscOrder;

    public int MaxDiscMultiplicityREven => _rEven.MaxDiscMultiplicity;
    public int MaxDiscMultiplicityROdd => _rOdd.MaxDiscMultiplicity;

    public string DisplayName =>
        $"The β-exotic at N = {N}: certified disc-layer reading, both R-parities " +
        $"(max multiplicity {Math.Max(MaxDiscMultiplicityREven, MaxDiscMultiplicityROdd)} < {BetaExoticDiscOrder})";

    public string Summary =>
        BetaExoticExcluded
            ? $"EXCLUDED, exactly, over ℚ(i): disc_Λ(F_res) has no root of multiplicity ≥ {BetaExoticDiscOrder} " +
              $"off q = 0 in either sector, and a Puiseux-3/2 point would need one. Certified layers " +
              $"[{Layers(_rOdd)}] (R-odd, residual degree {_rOdd.ResidualDegree}) and [{Layers(_rEven)}] " +
              $"(R-even, degree {_rEven.ResidualDegree}); the simple roots are the √-branch (defective) loci, and " +
              "the double roots are left unidentified (an order-2 zero may be diabolic, cubic, or two coincident " +
              "defective pairs; the theorem does not need to know which). Per-N, not a law: N = 7 unrun, and the " +
              "all-N scalar s₆ ≠ 0 stays open"
            : "NOT established at inspect time: the layer prime failed certification, or a root of multiplicity " +
              $"≥ {BetaExoticDiscOrder} appeared. Read the per-parity nodes; a certificate that does not certify " +
              "is a diagnostic, and must not be reported as an exclusion";

    private static string Layers(FoldResultantCertificate.CompleteReport r) =>
        string.Join(", ", r.DiscLayerDegrees);

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return ParityNode("R-odd", _rOdd);
            yield return ParityNode("R-even", _rEven);

            yield return new InspectableNode(
                displayName: "the multiplicity reads the Puiseux exponent",
                summary: "ord of a disc zero at a branch locus: 1 for a defective EP2 (exponent ½), 2 for a " +
                         "diabolic crossing, 2 for a cubic branch point (exponent ⅓), and 3 for the β-exotic " +
                         "(exponent 3/2, β(s) = [[0,s],[s²,0]], eigenvalues ±s^{3/2}). A 2-cycle with exponent e " +
                         "gives (λ₁−λ₂)² ~ (q−q*)^{2e}. Every other colliding pair contributes non-negative " +
                         "order, so a coincident collision raises the total (3+1, 3+2) and never masks it.",
                provenance: NodeProvenance.Stored);

            yield return new InspectableNode(
                displayName: "the one-way lift (why one certified prime is a proof)",
                summary: "reduction mod a prime ideal π_p not dividing lc_q(D) is a degree-preserving ring " +
                         "homomorphism: a factor (q−α)^m survives, and distinct roots can only MERGE. So " +
                         "max-mult(D mod p) ≥ max-mult(D), and reading 2 at one certified prime proves ≤ 2 over " +
                         "ℚ(i). Certification means the prime is good at BOTH ends of the q-axis: it attains the " +
                         "true deg_q D (no root escaped to q = ∞) and the true q-valuation (no nonzero root " +
                         "collapsed onto q = 0, where the q-power strip would discard it). A Hadamard bound " +
                         "certifies each of those two true values separately; a prime attaining BOTH is then " +
                         "searched for, and the report fails closed (DiscLayersCertified = false) if none does.",
                provenance: NodeProvenance.Stored);

            yield return new InspectableNode(
                displayName: "the scope boundary (what this witness does NOT show)",
                summary: "it does not exclude a cubic branch point (3×3 Jordan): that has ord disc = 2 and hides " +
                         "in the multiplicity-2 layer. It is ruled out at a count-drop by ELIMINATION, not by " +
                         "the order of the zero: a cubic point keeps one real branch and one conjugate pair on " +
                         "both sides of q*, so it cannot change the real count, while an EP2 can. That needs " +
                         "F_res real (T commutes with the reflection R at odd N, and the AT slopes are chirally " +
                         "paired; checked at N = 5, the AT step not derived in general). The β-exclusion above " +
                         "needs only the multiplicity bound. And it does not touch the all-N item: the codim-2 " +
                         "β-exotic genericity, reduced to s₆ ≠ 0 at every forced seed, remains open. This is a " +
                         "per-N certificate. It retires N = 5.",
                provenance: NodeProvenance.Stored);
        }
    }

    private static InspectableNode ParityNode(string label, FoldResultantCertificate.CompleteReport r) =>
        new(displayName: $"{label}: max disc multiplicity {r.MaxDiscMultiplicity}, layers [{Layers(r)}], " +
                         $"certified = {r.DiscLayersCertified.ToString(Inv)}",
            summary: $"residual degree {r.ResidualDegree.ToString(Inv)} (block dim " +
                     $"{r.BlockDimension.ToString(Inv)} = AT {r.AtDegree.ToString(Inv)} ⊎ residual). Certified " +
                     $"deg_q D = {r.TrueDiscriminantDegree.ToString(Inv)} and v_q(D) = " +
                     $"{r.TrueQValuationD.ToString(Inv)} at the layer prime {r.LayerPrime.ToString(Inv)}, from " +
                     $"{r.PrimesSampled.ToString(Inv)} sampled split primes against the lc-divisor bound " +
                     $"{r.LcDivisorBoundD.ToString(Inv)} (more sampled primes than can divide lc_q(D) or the " +
                     "coefficient of q^{v_q(D)}, so both are the true values). The multiplicity-1 layer carries " +
                     "the √-branch (defective) loci, the forced seeds among them; the multiplicity-2 layer is NOT " +
                     "identified here (an order-2 zero may be a diabolic crossing, a cubic branch point, or two " +
                     "coincident defective pairs). Only the ABSENCE of a multiplicity-3 layer is load-bearing.",
            provenance: NodeProvenance.Live);

    public InspectablePayload Payload => InspectablePayload.Empty;
}
