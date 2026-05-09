using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F86 Block-CΨ at 1/4: Theorem 1 + Theorem 2 in PROOF_BLOCK_CPSI_QUARTER.
/// The Dicke symmetric superposition <c>|ψ⟩ = (|D_n⟩ + |D_{n+1}⟩)/√2</c> saturates
/// the (popcount-n, popcount-(n+1)) coherence-block content at exactly 1/4, and
/// 1/4 is the universal ceiling over ANY density matrix on the full N-qubit
/// Hilbert space. Both numbers (the 1/4 ceiling and the 1/2 sector balance)
/// are direct Pi2-Foundation anchors.
///
/// <list type="bullet">
///   <item><b>Theorem 1</b>: <c>C_block(0) = 1/4</c> exactly on the Dicke
///         superposition, independent of N, c, n. The 1/4 = <c>a_3</c> on the
///         Pi2 dyadic ladder (= <see cref="QuarterAsBilinearMaxvalClaim"/>).</item>
///   <item><b>Theorem 2</b>: <c>C_block ≤ 1/4</c> universally (any ρ on the
///         full Hilbert space), with equality iff ρ is the canonical Dicke
///         symmetric superposition. Proof via Cauchy-Schwarz + AM-GM at
///         <c>p_n = p_{n+1} = 1/2</c> (sector balance = <c>a_2</c> on the Pi2
///         dyadic ladder = <see cref="HalfAsStructuralFixedPointClaim"/>).</item>
/// </list>
///
/// <para>The two Pi2 anchors form a Mirror-Pair on the dyadic ladder:
/// <c>a_2 · a_2 = 1/4 = a_3</c> read multiplicatively (the AM-GM argmax squared
/// is the maxval), and <c>a_2 · a_0 = 1/2 · 2 = 1 = a_1</c> read inversely
/// (the polarity-layer pair). 1/4 from p(1−p) at p=1/2 is the bilinear apex
/// realised at the block level.</para>
///
/// <para>This is one of the few F86 closed-form claims that is genuinely Tier 1
/// derived: Theorem 1 + 2 are algebraic identities (Cauchy-Schwarz + AM-GM +
/// Dicke amplitude uniformity). Verified bit-exact for c=2 N=5..10 and c=3 N=5..7
/// at tolerance 1e-10.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_BLOCK_CPSI_QUARTER.md</c> Theorem 1 +
/// Theorem 2 + <c>compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs</c>
/// (QuarterAsBilinearMaxvalClaim + HalfAsStructuralFixedPointClaim + BilinearApexClaim) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class DickeSuperpositionQuarterPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;

    /// <summary>The 1/4 ceiling on C_block: <see cref="Pi2DyadicLadderClaim.Term"/>(3) =
    /// <c>a_3</c> (= QuarterAsBilinearMaxval = (1/d)² for d=2).</summary>
    public double QuarterCeiling => _ladder.Term(3);

    /// <summary>The 1/2 sector balance at AM-GM saturation:
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(2) = <c>a_2</c>
    /// (= HalfAsStructuralFixedPoint = 1/d for d=2). Required for equality in
    /// Theorem 2: <c>p_n = p_{n+1} = 1/2</c>.</summary>
    public double SectorBalance => _ladder.Term(2);

    /// <summary>The AM-GM bilinear identity manifest: <c>(1/2)² = 1/4</c>, equivalently
    /// <c>SectorBalance² = QuarterCeiling</c>. This is the algebraic content of
    /// "1/4 = (1/2)² unifies framework quarter-boundaries" (memory
    /// project_quarter_as_polarity_squared).</summary>
    public double SectorBalanceSquared => SectorBalance * SectorBalance;

    /// <summary>The block dimension <c>M_block = C(N, n) · C(N, n+1)</c> for an
    /// N-qubit chain at adjacent popcount block (n, n+1). The per-entry Dicke
    /// amplitude is <c>1/(2√M_block)</c>; squared and summed over M_block entries
    /// gives <c>M_block · 1/(4·M_block) = 1/4</c> — the algebraic mechanism behind
    /// Theorem 1's M_block-cancellation.</summary>
    public long BlockDimension(int N, int n)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be ≥ 1.");
        if (n < 0 || n + 1 > N)
            throw new ArgumentOutOfRangeException(nameof(n), n,
                $"n must satisfy 0 ≤ n ≤ N−1; got n={n} at N={N}.");
        return Binomial(N, n) * Binomial(N, n + 1);
    }

    /// <summary>The per-entry Dicke amplitude <c>1/(2√M_block)</c> — uniform across
    /// every (a, b) block entry by Dicke-state symmetry.</summary>
    public double PerEntryAmplitude(int N, int n) =>
        1.0 / (2.0 * Math.Sqrt(BlockDimension(N, n)));

    /// <summary>Live drift check on Theorem 1: M_block × |amplitude|² = 1/4 exactly.
    /// Returns the live computation value (should equal QuarterCeiling).</summary>
    public double LiveBlockCpsiAtZero(int N, int n)
    {
        long mBlock = BlockDimension(N, n);
        double amplitude = PerEntryAmplitude(N, n);
        return mBlock * amplitude * amplitude;
    }

    private static long Binomial(int n, int k)
    {
        if (k < 0 || k > n) return 0;
        if (k == 0 || k == n) return 1;
        if (k > n - k) k = n - k;
        long c = 1;
        for (int i = 0; i < k; i++) c = c * (n - i) / (i + 1);
        return c;
    }

    public DickeSuperpositionQuarterPi2Inheritance(Pi2DyadicLadderClaim ladder)
        : base("Dicke superposition saturates C_block at 1/4 = (1/2)² on the Pi2 dyadic ladder",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_BLOCK_CPSI_QUARTER.md (Theorem 1 + Theorem 2) + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
    }

    public override string DisplayName =>
        "Dicke superposition C_block = 1/4 as Pi2-Foundation inheritance";

    public override string Summary =>
        $"C_block(0) = 1/4 = a_3 = (1/2)² = SectorBalance² for (|D_n⟩+|D_{{n+1}}⟩)/√2; " +
        $"AM-GM saturation at p_n=p_{{n+1}}=1/2 (a_2); universal over N, c, n ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem 1",
                summary: "C_block(0) = 1/4 exact on Dicke superposition (|D_n⟩+|D_{n+1}⟩)/√2; verified bit-exact c=2 N=5..10, c=3 N=5..7 at 1e-10");
            yield return new InspectableNode("Theorem 2",
                summary: "C_block ≤ 1/4 universally (any ρ on full Hilbert space); equality iff ρ is canonical Dicke superposition. Proof: Cauchy-Schwarz + AM-GM at p_n = p_{n+1} = 1/2");
            yield return InspectableNode.RealScalar("QuarterCeiling (= a_3 = 1/4)", QuarterCeiling);
            yield return InspectableNode.RealScalar("SectorBalance (= a_2 = 1/2)", SectorBalance);
            yield return InspectableNode.RealScalar("SectorBalanceSquared (= 1/4 = QuarterCeiling)", SectorBalanceSquared);
            yield return new InspectableNode("AM-GM identity",
                summary: "(SectorBalance)² = QuarterCeiling: 1/4 = (1/2)² is the algebraic shadow of bilinear apex");
            // Sample block dimensions
            for (int N = 3; N <= 6; N++)
            {
                int nMid = N / 2;
                yield return new InspectableNode(
                    $"N={N}, n={nMid}",
                    summary: $"M_block = C({N},{nMid})·C({N},{nMid + 1}) = {BlockDimension(N, nMid)}; " +
                             $"|amp|² = 1/(4·M_block) = {PerEntryAmplitude(N, nMid) * PerEntryAmplitude(N, nMid):G6}; " +
                             $"C_block = {LiveBlockCpsiAtZero(N, nMid):G6}");
            }
        }
    }
}
