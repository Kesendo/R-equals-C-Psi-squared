using System;
using System.Collections.Generic;
using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>One N's reading of the bipartite gauge: how many (wKet, wBra) blocks there are, the worst residual of
/// 𝒟L(q)𝒟 = L(−q) over all of them (0.0 bit-for-bit, not machine zero), and how the blocks split into
/// gauge-commuting and gauge-anticommuting under the site reflection.</summary>
public sealed record BipartiteGaugeReading(int N, int Blocks, double WorstResidual, int Commuting, int Anticommuting);

/// <summary>THE BIPARTITE GAUGE AND ITS COMMUTATION CRITERION: why a coherence block's R-sectors are real matrix
/// families on one parity and conjugate-swapped on the other, decided by the single character (p+q)(N−1).
///
/// <para>The gauge is 𝒟 = diag((−1)^(σ(a)+σ(b))) on the (wKet, wBra) block, σ the sum of OCCUPIED SITE INDICES
/// (<see cref="WeightCoherenceBlock.BipartiteGaugeSigns"/>). Write L(q) = A + q·C with A the real Absorption-Theorem
/// diagonal and C the pure-imaginary nearest-neighbour hopping. Conjugating by 𝒟 leaves A alone and negates C,
/// because a hop moves one excitation by one site and so flips the bipartite parity on exactly one side. Hence
/// 𝒟·L(q)·𝒟 = L(−q) at every complex q, and at real q that is also conj(L), the two forms differing only by the
/// pencil reality conj L(q) = L(−q̄). Both are ±1 sign rearrangements of doubles, so this witness reports a residual
/// that is exactly 0.0 rather than a machine zero, and a nonzero one would be a fact about the construction.</para>
///
/// <para>THE CRITERION. σ(rev(c)) = w·(N−1) − σ(c), so reversing both sides multiplies every gauge sign by the SAME
/// global factor (−1)^((p+q)(N−1)). So 𝒟 commutes with the site reflection R exactly when (p+q)(N−1) is even and
/// anticommutes otherwise, with nothing in between
/// (<see cref="WeightCoherenceBlock.BipartiteGaugeCommutesWithReflection"/>). Where it commutes, 𝒟 restricts to each
/// R-parity sector, so the sector is a real matrix family and its spectrum is conjugation-closed; where it
/// anticommutes, the gauge swaps the two sectors, so what it closes is the PAIR, the spectrum of one being the conjugate of the other, and not either sector on its own. This completes the clause
/// PROOF_CODIM1_BY_ADDITIVITY §7's grading paragraph leaves open, "P, Q, 𝒟, T all commute with R up to a
/// block-global sign", which never says when that sign is −1.</para>
///
/// <para>It SUBSUMES the PARITY GATE of the odd-N effect, the sector-closure half of what the path-k story is built on. At odd N the factor N−1 is even, so
/// every block commutes and every sector is conjugation-closed regardless of p and q; the parity gate that reads
/// as "odd N" at p+q = 3 is this one character read at that weight. The criterion is SUFFICIENT and not necessary:
/// nothing here excludes some other antiunitary closing a sector on the odd parity, it is only unobserved.</para>
///
/// <para>Gauge and fold are INDEPENDENT levers, which is the reason to name the gauge as the cause. (0,3) and (2,3)
/// at N=6 are fixed by a fold leg and still anticommute; (1,3)@N=4 commutes without being fold-fixed at all. The
/// gauge makes a sector a real family; the fold pins loci to the line and doubles the ones off it.</para>
///
/// <para>SCOPE, both fences read live below: Δ = 0 and no longitudinal field. Each contributes an IMAGINARY
/// diagonal, which a diagonal of signs cannot reach, so the identity breaks by O(1) rather than degrading. Open
/// chain with nearest-neighbour hops is load-bearing too, since a next-nearest bond does not flip the bipartite
/// parity.</para>
///
/// <para>NAME COLLISION, because this repo spends the letter T three times: 𝒟 here is the bipartite gauge, called T
/// in <c>FoldResultantCertificate</c> and in the OpenArcs registry. It is NOT the fold antiunitary T = P·K, and NOT the
/// transpose leg; the bra/ket complement permutations of the cross-fold are a fourth neighbour but are called
/// P and Q, never T. Its p+q = 3 specialisation is
/// <c>BetaExoticPerNExclusionClaim</c>'s metric, which states the same commutation as "3(N+1) is even"; that has
/// the same parity as the general 3(N−1) and is not a disagreement.</para>
///
/// <para>Gate <c>BipartiteGaugeCriterionTests</c> (Category EVEN_DISC13_SCOUT). Reading
/// <c>experiments/F89_PATH_K_DIABOLIC.md</c> and the arc <c>sideways_spin_ladder</c>.</para></summary>
public sealed class BipartiteGaugeWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private static readonly int[] SweepN = { 4, 5, 6, 7 };
    private const double GenericQ = 0.7;

    /// <summary>Sweeps every (wKet, wBra) block at this N, recomputed at inspect time: the worst residual of
    /// 𝒟L(q)𝒟 = L(−q) and the commuting/anticommuting split under the reflection.</summary>
    public static BipartiteGaugeReading Read(int n, double q)
    {
        double worst = 0.0;
        int commuting = 0, anticommuting = 0, blocks = 0;
        for (int wKet = 0; wKet <= n; wKet++)
            for (int wBra = 0; wBra <= n; wBra++)
            {
                var signs = WeightCoherenceBlock.BipartiteGaugeSigns(n, wKet, wBra);
                var L = WeightCoherenceBlock.Build(n, wKet, wBra, new Complex(q, 0));
                var Lm = WeightCoherenceBlock.Build(n, wKet, wBra, new Complex(-q, 0));
                for (int i = 0; i < signs.Length; i++)
                    for (int j = 0; j < signs.Length; j++)
                        worst = Math.Max(worst, (signs[i] * signs[j] * L[i, j] - Lm[i, j]).Magnitude);

                if (WeightCoherenceBlock.BipartiteGaugeCommutesWithReflection(n, wKet, wBra)) commuting++;
                else anticommuting++;
                blocks++;
            }
        return new BipartiteGaugeReading(n, blocks, worst, commuting, anticommuting);
    }

    /// <summary>The worst |𝒟L𝒟 − conj(L)| on one block at a given Δ and optional field: 0.0 exactly inside the
    /// scope, O(1) outside it. The fence read, not a tolerance.</summary>
    public static double ScopeBreakResidual(int n, int wKet, int wBra, double delta, double[]? field = null)
    {
        var signs = WeightCoherenceBlock.BipartiteGaugeSigns(n, wKet, wBra);
        var L = WeightCoherenceBlock.Build(n, wKet, wBra, new Complex(GenericQ, 0), delta, field);
        double worst = 0.0;
        for (int i = 0; i < signs.Length; i++)
            for (int j = 0; j < signs.Length; j++)
                worst = Math.Max(worst, (signs[i] * signs[j] * L[i, j] - Complex.Conjugate(L[i, j])).Magnitude);
        return worst;
    }

    public string DisplayName =>
        "BipartiteGaugeWitness (𝒟L𝒟 = conj(L) on every block, and 𝒟 commutes with the reflection iff (p+q)(N−1) is even)";

    public string Summary
    {
        get
        {
            var r6 = Read(6, GenericQ);
            var r7 = Read(7, GenericQ);
            return "the staggered diagonal 𝒟 = diag((−1)^(σ(a)+σ(b))), σ the sum of occupied site indices, satisfies " +
                   $"𝒟·L(q)·𝒟 = L(−q) on EVERY (wKet,wBra) block and, at real q, 𝒟·L·𝒟 = conj(L): worst residual " +
                   $"{r6.WorstResidual.ToString("0.0###", Inv)} over all {r6.Blocks} blocks at N=6, exactly 0.0 rather than " +
                   "machine zero, because the whole operation is ±1 sign flips. The content is the OTHER half: since " +
                   "σ(rev(c)) = w·(N−1) − σ(c), the site reflection rescales every gauge sign by one global factor " +
                   "(−1)^((p+q)(N−1)), so 𝒟 COMMUTES with the reflection iff (p+q)(N−1) is even and anticommutes " +
                   "otherwise, never anything between. Where it commutes, 𝒟 restricts to each R-parity sector and that " +
                   "sector is a real matrix family with a conjugation-closed spectrum; where it anticommutes, the gauge " +
                   $"swaps the sectors. At odd N the factor N−1 is even, so ALL {r7.Blocks} blocks commute at N=7 and the " +
                   $"criterion is automatic there, while at N=6 the split is {r6.Commuting} commuting against " +
                   $"{r6.Anticommuting} anticommuting. That is why one character subsumes the parity-gate half of the path-k arc's odd-N effect: " +
                   "at p+q = 3 the criterion reads \"N odd\". It is SUFFICIENT and not necessary; a second closing " +
                   "antiunitary on the odd parity is unobserved, not excluded. Gauge and fold are independent: " +
                   "(0,3)@N=6 and (2,3)@N=6 are fold-fixed yet anticommute, (1,3)@N=4 commutes without being fold-fixed. " +
                   "The gauge makes the sector a real family; the fold pins loci to the line. Scope, both read live " +
                   "below: Δ = 0 and no longitudinal field, each of which adds an imaginary diagonal a sign matrix " +
                   "cannot reach. This is the clause PROOF_CODIM1 §7 leaves open when it says P, Q, 𝒟, T commute with R " +
                   "\"up to a block-global sign\" without saying when that sign is −1.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            foreach (int n in SweepN)
            {
                var r = Read(n, GenericQ);
                yield return new InspectableNode(
                    displayName: $"N={n}: {r.Blocks} blocks, {r.Commuting} commuting / {r.Anticommuting} anticommuting",
                    summary: $"worst |𝒟L(q)𝒟 − L(−q)| = {r.WorstResidual.ToString("0.0###", Inv)} over every block at " +
                             $"q={GenericQ.ToString("0.#", Inv)} (exactly 0.0, a sign rearrangement) " +
                             (n % 2 == 1
                                 ? "⟹ odd N: N−1 is even, so the character is even at every weight pair and EVERY sector is conjugation-closed."
                                 : "⟹ even N: the character alternates with p+q, so half the blocks have real-family sectors and half have swapped ones."),
                    provenance: NodeProvenance.Live);
            }

            yield return new InspectableNode(
                displayName: "the arc's four blocks, by character",
                summary: "the weight pairs the sideways-ladder arc reads, each with its criterion verdict",
                children: NamedBlocks(),
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "gauge and fold are independent levers",
                summary: FoldIndependence(),
                provenance: NodeProvenance.Live);

            double dOk = ScopeBreakResidual(6, 1, 3, 0.0);
            double dBreak = ScopeBreakResidual(6, 1, 3, 0.85);
            yield return new InspectableNode(
                displayName: "scope fence: Δ = 0",
                summary: $"(1,3)@N=6: |𝒟L𝒟 − conj(L)| = {dOk.ToString("0.0###", Inv)} at Δ=0 and " +
                         $"{dBreak.ToString("0.0###", Inv)} at Δ=0.85 ⟹ the Δ·ZZ term is an IMAGINARY diagonal, which a " +
                         "diagonal of signs cannot touch, so the identity does not degrade, it breaks. This is " +
                         "PROOF_CODIM1 §7's own Δ = 0 scope on ingredient (iv), measured rather than asserted.",
                provenance: NodeProvenance.Live);

            double fBreak = ScopeBreakResidual(5, 1, 2, 0.0, new[] { 0.3, -0.1, 0.45, 0.2, -0.35 });
            yield return new InspectableNode(
                displayName: "scope fence: no longitudinal field",
                summary: $"(1,2)@N=5 with a site-dependent Z-field: |𝒟L𝒟 − conj(L)| = {fBreak.ToString("0.0###", Inv)} " +
                         "⟹ the field is bit-flip-odd and again enters as an imaginary diagonal. Together with the Δ row " +
                         "this fixes what the gauge actually needs: a REAL diagonal, and hops that flip the bipartite parity.",
                provenance: NodeProvenance.Live);
        }
    }

    private static IEnumerable<IInspectable> NamedBlocks()
    {
        (int n, int p, int q, string note)[] rows =
        {
            (4, 1, 2, "the path-k (SE,DE) block, and the committed per-sector q-evenness read there is False"),
            (6, 1, 3, "the census block: sectors are real families, conjugation-closure residual ≤ 2e-13, q-evenness True"),
            (8, 1, 2, "the N=8 (SE,DE) block, whose no-real-defective statement is a theorem by the gcd certificate"),
            (8, 1, 4, "fold-fixed AND odd character: the case that separates the gauge from the fold"),
        };
        foreach (var (n, p, q, note) in rows)
        {
            bool commutes = WeightCoherenceBlock.BipartiteGaugeCommutesWithReflection(n, p, q);
            yield return new InspectableNode(
                displayName: $"({p},{q})@N={n}: character (p+q)(N−1) = {(p + q) * (n - 1)}",
                summary: (commutes
                             ? "EVEN ⟹ 𝒟 commutes, restricts to each R-sector, sectors are conjugation-closed real families. "
                             : "ODD ⟹ 𝒟 anticommutes, swaps the two R-sectors, and supplies nothing about either. ") + note,
                provenance: NodeProvenance.Live);
        }
    }

    private static string FoldIndependence()
    {
        static bool FoldFixed(int n, int p, int q) => n % 2 == 0 && (p == n / 2 || q == n / 2);
        bool a = WeightCoherenceBlock.BipartiteGaugeCommutesWithReflection(6, 0, 3);
        bool b = WeightCoherenceBlock.BipartiteGaugeCommutesWithReflection(6, 2, 3);
        bool c = WeightCoherenceBlock.BipartiteGaugeCommutesWithReflection(4, 1, 3);
        return $"(0,3)@N=6 fold-fixed={FoldFixed(6, 0, 3)} commutes={a}; (2,3)@N=6 fold-fixed={FoldFixed(6, 2, 3)} " +
               $"commutes={b}; (1,3)@N=4 fold-fixed={FoldFixed(4, 1, 3)} commutes={c} ⟹ self-foldedness is NEITHER " +
               "necessary NOR sufficient for conjugation-closure, so the cause of the closure is the gauge and not the " +
               "fold. The two do different work: the gauge makes a sector a real matrix family, the fold pins loci to " +
               "the line Re λ = −N and doubles the ones off it.";
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
    public NodeProvenance Provenance => NodeProvenance.Live;
}
