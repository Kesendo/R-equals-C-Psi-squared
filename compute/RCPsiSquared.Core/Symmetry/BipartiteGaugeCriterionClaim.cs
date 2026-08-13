using System;
using System.Collections.Generic;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F151 (Tier 1 derived): the BIPARTITE GAUGE and its commutation criterion, which decides whether a
/// coherence block's reflection sectors are real matrix families or conjugate-swapped, out of a single character.
///
/// <para>On the (wKet, wBra) block let 𝒟 = diag((−1)^(σ(a)+σ(b))) with σ(c) the sum of the OCCUPIED SITE INDICES
/// of a computational-basis config; this is PROOF_CODIM1_BY_ADDITIVITY §7 ingredient (iv) in block form. Writing
/// L(q) = A + q·C with A the real Absorption-Theorem diagonal and C the pure-imaginary nearest-neighbour hopping,
/// a hop moves one excitation by one site, changes σ by ±1 and so flips the gauge sign on exactly one side, while
/// 𝒟 leaves the diagonal alone. Hence 𝒟·L(q)·𝒟 = L(−q) at every complex q and, at real q, 𝒟·L·𝒟 = conj(L), the
/// second being the first composed with the pencil reality conj L(q) = L(−q̄). Both are ±1 sign rearrangements of
/// the entries, so the residual is 0.0 BIT FOR BIT and not a machine zero.</para>
///
/// <para>THE CRITERION, and this is where the content is. σ(rev(c)) = w·(N−1) − σ(c), so the site reflection
/// multiplies every gauge sign by the SAME global factor (−1)^((p+q)(N−1)). Therefore 𝒟 commutes with the
/// reflection exactly when (p+q)(N−1) is even and anticommutes otherwise, with nothing in between. Where it
/// commutes, 𝒟 restricts to each R-parity sector, that sector is a real matrix family and its spectrum is
/// conjugation-closed; where it anticommutes, the gauge swaps the two sectors, so what it closes is the PAIR, the spectrum of one being the conjugate of the other, and not either sector on its own.
/// This completes a clause §7's grading paragraph leaves open, which records that "P, Q, 𝒟, T all commute with R
/// up to a block-global sign" and never says when that sign is −1.</para>
///
/// <para>It SUBSUMES the PARITY GATE of the odd-N effect, the sector-closure half of what the path-k story is built on: at odd N the factor N−1 is even, so
/// every block commutes and every sector is conjugation-closed regardless of the weights, and the parity gate that
/// reads as "odd N" at p+q = 3 is this one character read at that weight. Note also that the BLOCK-level statement
/// is vacuous as a realness criterion, since 𝒟 makes every full block's discriminant real at any (p, q̃); all the
/// content sits at SECTOR level, in the commutation half.</para>
///
/// <para>GAUGE AND FOLD ARE INDEPENDENT LEVERS, which is why the cause of conjugation-closure is named as the gauge
/// and not the self-fold: (0,3) and (2,3) at N=6 are fixed by a fold leg and still anticommute, while (1,3)@N=4
/// commutes without being fold-fixed at all. The gauge makes a sector a real matrix family; the fold pins loci to
/// the line and doubles the ones off it.</para>
///
/// <para>SCOPE. Δ = 0 and no longitudinal field, each of which contributes an IMAGINARY diagonal that a diagonal of
/// signs cannot reach, so the identity breaks rather than degrades (measured 9.52 at Δ = 0.85 on (1,3)@N=6, 3.08
/// for a site-dependent field on (1,2)@N=5). Real J, open chain, nearest-neighbour hops, since a next-nearest bond
/// does not flip the bipartite parity. The criterion is SUFFICIENT and not necessary: on the odd parity the gauge
/// supplies nothing, which is not a proof that no other antiunitary closes those sectors, and the non-closure there
/// is measured on four probed blocks rather than derived. The descent of q-evenness from the sector charpoly to
/// F_res additionally needs the AT strands to pair under 𝒟, gated at (1,3)@N=6 and not at N=8.</para>
///
/// <para>NAME COLLISION, because the repo spends the letter T three times: 𝒟 here is the bipartite gauge, called T
/// in <c>FoldResultantCertificate</c> and in the OpenArcs registry. It is NOT the fold antiunitary T = P·K
/// (<see cref="F89BranchLocusPalindromeClaim"/>), and NOT the transpose leg; the bra/ket complement
/// permutations of <see cref="F89CrossFoldSimilarityClaim"/> are a fourth neighbour but are called P and Q,
/// never T. Its p+q = 3 specialisation is
/// <see cref="BetaExoticPerNExclusionClaim"/>'s metric, which states the same commutation as "3(N+1) is even";
/// that has the same parity as the general 3(N−1) and is not a disagreement.</para>
///
/// <para>Gate <c>BipartiteGaugeCriterionTests</c> (11 gates, Category EVEN_DISC13_SCOUT, under a second), the first
/// committed gate for the statements themselves; before it they were carried only through their consequences.
/// Live witness <c>inspect --root gauge</c> (<c>BipartiteGaugeWitness</c>). Implementation
/// <see cref="WeightCoherenceBlock.BipartiteGaugeSigns"/> and
/// <see cref="WeightCoherenceBlock.BipartiteGaugeCommutesWithReflection"/>. Anchors:
/// <c>docs/ANALYTICAL_FORMULAS.md</c> (F151) + <c>experiments/F89_PATH_K_DIABOLIC.md</c> + the arc
/// <c>sideways_spin_ladder</c>.</para></summary>
public sealed class BipartiteGaugeCriterionClaim : Claim
{
    // Parent-edge marker for Schicht-1 wiring: the sublattice gauge one level down, K·H·K = −H on the
    // single-excitation site basis. THIS is the gauge's ancestor; every committed claim that uses the
    // sublattice grading parents here (DeadSetLaw, MirrorOrderSorting, SeedExistenceCounting), and 𝒟 is
    // its Liouville-space lift, K ⊗ K in block form.
    public ChiralKClaim SublatticeGauge { get; }
    // Parent-edge marker for Schicht-1 wiring: the fold lattice, whose §7 grading clause this criterion
    // completes, and where 𝒟 already appears as ingredient (iv). The fold is where the gauge is USED, not
    // where it comes from, which is why it is the second parent and not the first.
    public F89CrossFoldSimilarityClaim FoldLattice { get; }

    /// <summary>The character (p+q)(N−1) whose parity is the whole criterion.</summary>
    public static int Character(int n, int wKet, int wBra) => (wKet + wBra) * (n - 1);

    /// <summary>Whether the gauge commutes with the site reflection on this block, so that its R-sectors are
    /// conjugation-closed real families. Delegates to the Core primitive rather than restating the rule.</summary>
    public static bool SectorsAreConjugationClosed(int n, int wKet, int wBra) =>
        WeightCoherenceBlock.BipartiteGaugeCommutesWithReflection(n, wKet, wBra);

    public BipartiteGaugeCriterionClaim(ChiralKClaim sublatticeGauge, F89CrossFoldSimilarityClaim foldLattice)
        : base("F151 bipartite gauge criterion: D = diag((-1)^(sigma(a)+sigma(b))), sigma the sum of occupied site " +
               "indices, satisfies D L(q) D = L(-q) at every complex q and D L D = conj(L) at real q, exactly (a " +
               "sign rearrangement, residual 0.0 bit-for-bit); and since sigma(rev(c)) = w(N-1) - sigma(c) the site " +
               "reflection rescales every gauge sign by ONE global factor, so D commutes with the reflection iff " +
               "(p+q)(N-1) is even and anticommutes otherwise, and where it commutes a block's R-sectors are " +
               "conjugation-closed real matrix families; this completes PROOF_CODIM1 section-7's open clause 'P, Q, D, " +
               "T all commute with R up to a block-global sign' and subsumes the parity-gate half of the path-k odd-N effect (at p+q=3 it " +
               "reads 'N odd'). Sufficient, not necessary; Delta=0, real J, open chain with NN hops",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md + " +
               "experiments/F89_PATH_K_DIABOLIC.md + " +
               "docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md")
    {
        SublatticeGauge = sublatticeGauge ?? throw new ArgumentNullException(nameof(sublatticeGauge));
        FoldLattice = foldLattice ?? throw new ArgumentNullException(nameof(foldLattice));
    }

    public override string DisplayName =>
        "F151: the bipartite gauge criterion (the gauge commutes with the reflection iff (p+q)(N-1) is even, and then the R-sectors are conjugation-closed)";

    public override string Summary =>
        "D L(q) D = L(-q) on EVERY block, exactly (0.0 bit-for-bit, a +-1 sign rearrangement), and D commutes with " +
        "the site reflection iff the character (p+q)(N-1) is even: on that parity D restricts to each R-sector, so " +
        "the sector is a real matrix family with a conjugation-closed spectrum; on the other it swaps them and " +
        "supplies nothing. Completes PROOF_CODIM1 section-7's 'up to a block-global sign' clause, and subsumes the " +
        "parity-gate half of the path-k odd-N effect, since N-1 even makes every block commute at odd N. Gauge and fold are INDEPENDENT: " +
        "(0,3)/(2,3)@N=6 are fold-fixed yet anticommute, (1,3)@N=4 commutes without being fold-fixed. Sufficient, " +
        $"not necessary; Delta=0 and no longitudinal field, both fences measured to break at O(1) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the criterion, in one line",
                summary: "sigma(rev(c)) = w(N-1) - sigma(c), so reversing both sides multiplies EVERY gauge sign by the " +
                         "same factor (-1)^((p+q)(N-1)). That is why the alternative is commute-or-anticommute with nothing " +
                         "in between, and why one character decides the whole block. Characters at the arc's blocks: " +
                         $"{Character(4, 1, 2)} at (1,2)@N=4 (odd), {Character(6, 1, 3)} at (1,3)@N=6 (even), " +
                         $"{Character(8, 1, 2)} at (1,2)@N=8 (odd), {Character(8, 1, 4)} at (1,4)@N=8 (odd).");
            yield return new InspectableNode("what it subsumes: the odd-N PARITY GATE is this character at p+q = 3",
                summary: "at odd N the factor N-1 is even, so every block commutes and every sector is conjugation-closed " +
                         "regardless of the weights (measured: 36 of 36 blocks at N=5, 64 of 64 at N=7, against 25 of 49 at " +
                         "N=6). The path-k arc's parity gate, 'real-q diabolics appear only at odd N', is this criterion read " +
                         "at p+q = 3, where (p+q)(N-1) = 3(N-1) is even exactly when N is odd.");
            yield return new InspectableNode("gauge and fold are independent levers",
                summary: "(0,3)@N=6 and (2,3)@N=6 are fixed by a fold leg and still anticommute; (1,3)@N=4 commutes without " +
                         "being fold-fixed. So self-foldedness is neither necessary nor sufficient for conjugation-closure, and " +
                         "the cause is the gauge. The two do different work: the gauge makes a sector a real matrix family, the " +
                         "fold pins loci to the line Re lambda = -N and doubles the ones off it.");
            yield return new InspectableNode("what the block-level identity does NOT give",
                summary: "D makes every FULL block's discriminant real at any (p, q-tilde), so block-level realness discriminates " +
                         "nothing at all; the content is entirely at SECTOR level and therefore entirely in the commutation half. " +
                         "The criterion is also SUFFICIENT and not necessary: on the odd parity the gauge closes only the PAIR and neither sector on its own, which " +
                         "is not a proof that no other antiunitary closes those sectors, and the non-closure there is MEASURED on " +
                         "four probed blocks.");
            yield return new InspectableNode("live lab (the witness)",
                summary: "BipartiteGaugeWitness recomputes the whole block sweep at inspect time, including both scope fences: " +
                         "inspect --root gauge. Gate BipartiteGaugeCriterionTests (Category EVEN_DISC13_SCOUT), the first " +
                         "committed gate for the statements themselves.");
        }
    }
}
