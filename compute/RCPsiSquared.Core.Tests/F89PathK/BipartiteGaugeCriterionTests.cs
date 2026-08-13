using System;
using System.Numerics;
using RCPsiSquared.Core.F89PathK;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>The BIPARTITE GAUGE and its commutation criterion, gated directly for the first time
/// (sideways_spin_ladder arc). Three documents lean on this pair of facts, and until now both were
/// carried only through their CONSEQUENCES: that the per-sector charpoly is even in q at (1,3)@N=6
/// and is not at (1,2)@N=4. The statements themselves are entry-wise sign rearrangements, so they
/// are gated here for what they are, with no eigensolver, no charpoly and no tolerance.
///
/// THE GAUGE. 𝒟 = diag((−1)^(σ(a)+σ(b))) on the (wKet, wBra) block, σ the sum of OCCUPIED SITE
/// INDICES (WeightCoherenceBlock.BipartiteGaugeSigns). Writing L(q) = A + q·C with A the real AT
/// diagonal and C the pure-imaginary nearest-neighbour hopping, a hop changes σ by ±1 on exactly one
/// side and so flips the gauge sign, giving
///
///     𝒟·L(q)·𝒟 = L(−q)   at every complex q,   and   𝒟·L(q)·𝒟 = conj(L(q))   at real q,
///
/// the second being the first composed with the pencil reality conj L(q) = L(−q̄). PROOF_CODIM1
/// §7 ingredient (iv) states the first at Δ = 0; the registry and F89_PATH_K_DIABOLIC state the
/// second. Both are ±1 sign flips on doubles, so the residual is 0.0 BIT-FOR-BIT and the assertions
/// below compare exactly. A nonzero residual here would be a fact about the construction, not a
/// rounding to be tolerated.
///
/// THE CRITERION. σ(rev(c)) = w·(n−1) − σ(c), so reversing both sides multiplies every gauge sign
/// by the SAME factor (−1)^((p+q)(n−1)). Hence 𝒟 commutes with the site reflection R iff
/// (p+q)(n−1) is even and anticommutes otherwise, with nothing in between. This completes the
/// clause PROOF_CODIM1 §7's grading paragraph leaves open ("P, Q, 𝒟, T all commute with R up to a
/// block-global sign", without saying when the sign is −1). Where it commutes, 𝒟 restricts to each
/// R-parity sector and that sector is a real matrix family; where it anticommutes, the gauge swaps
/// the sectors, so what it closes is the PAIR (one sector's spectrum is the conjugate of the other's)
/// and not either sector on its own.
///
/// NAME COLLISION, because the repo spends this letter three times: 𝒟 here is the bipartite gauge.
/// It is NOT the fold antiunitary T = P·K (F89BranchLocusPalindromeClaim) and NOT the transpose leg
/// (PROOF_CODIM1 §7 ingredient (i)); the bra/ket complement permutations are a fourth neighbour but are
/// called P and Q, never T. Its p+q = 3
/// specialisation is BetaExoticPerNExclusionClaim's metric, which states the same commutation as
/// "3(N+1) is even"; that has the same parity as the general 3(N−1) and is not a disagreement.
///
/// SCOPE, both fences gated below as negative controls: Δ = 0 and no longitudinal field. Each adds
/// an IMAGINARY diagonal, which a diagonal sign matrix cannot touch, so the identity breaks.</summary>
public sealed class BipartiteGaugeCriterionTests
{
    private readonly ITestOutputHelper _out;
    public BipartiteGaugeCriterionTests(ITestOutputHelper output) => _out = output;

    /// <summary>Every block up to N=6, both identities, exactly. The q values deliberately include a
    /// complex one, where the two identities part company: 𝒟L𝒟 = L(−q) still holds, conj(L) does not.</summary>
    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [Trait("Category", "EVEN_DISC13_SCOUT")]
    public void Gauge_ConjugatesTheBlock_ExactlyOnEveryBlock(int n)
    {
        var qs = new[] { new Complex(0.7, 0.0), new Complex(2.0, 0.0), new Complex(0.5, 1.3) };
        int blocks = 0;
        for (int wKet = 0; wKet <= n; wKet++)
            for (int wBra = 0; wBra <= n; wBra++)
            {
                var signs = WeightCoherenceBlock.BipartiteGaugeSigns(n, wKet, wBra);
                foreach (var q in qs)
                {
                    var L = WeightCoherenceBlock.Build(n, wKet, wBra, q);
                    var Lm = WeightCoherenceBlock.Build(n, wKet, wBra, -q);
                    int dim = signs.Length;
                    for (int i = 0; i < dim; i++)
                        for (int j = 0; j < dim; j++)
                        {
                            Complex gauged = signs[i] * signs[j] * L[i, j];
                            Assert.True(gauged == Lm[i, j],
                                $"𝒟L(q)𝒟 = L(−q) must be exact: N={n} ({wKet},{wBra}) q={q} at [{i},{j}], " +
                                $"{gauged} vs {Lm[i, j]}");
                            if (q.Imaginary == 0.0)
                                Assert.True(gauged == Complex.Conjugate(L[i, j]),
                                    $"𝒟L𝒟 = conj(L) must be exact at real q: N={n} ({wKet},{wBra}) at [{i},{j}]");
                        }
                }
                blocks++;
            }
        _out.WriteLine($"N={n}: {blocks} blocks, both identities exact (0.0 bit-for-bit) at 3 couplings each");
    }

    /// <summary>The criterion itself, read off the signs rather than off the formula: the reflection
    /// multiplies the gauge sign by one GLOBAL factor, and that factor is (−1)^((p+q)(n−1)).</summary>
    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [Trait("Category", "EVEN_DISC13_SCOUT")]
    public void GaugeReflectionCommutator_IsExactlyTheCharacterParity(int n)
    {
        int commuting = 0, anticommuting = 0;
        for (int wKet = 0; wKet <= n; wKet++)
            for (int wBra = 0; wBra <= n; wBra++)
            {
                var signs = WeightCoherenceBlock.BipartiteGaugeSigns(n, wKet, wBra);
                var perm = WeightCoherenceBlock.ReflectionPermutation(n, wKet, wBra);
                bool predicted = WeightCoherenceBlock.BipartiteGaugeCommutesWithReflection(n, wKet, wBra);
                int expected = ((wKet + wBra) * (n - 1)) % 2 == 0 ? 1 : -1;
                Assert.Equal(expected == 1, predicted);

                for (int t = 0; t < signs.Length; t++)
                    Assert.True(signs[perm[t]] == expected * signs[t],
                        $"the reflection must rescale the gauge by ONE global sign: N={n} " +
                        $"({wKet},{wBra}) character {(wKet + wBra) * (n - 1)} at t={t}");

                if (predicted) commuting++; else anticommuting++;
            }
        _out.WriteLine($"N={n}: {commuting} blocks commuting, {anticommuting} anticommuting, " +
                       $"every one of them global-sign exact");
    }

    /// <summary>The blocks the arc names, with their committed character values, so a future reader can
    /// find the criterion's verdict for each without recomputing: 9 at (1,2)@N=4 (odd, gauge swaps the
    /// sectors, and the committed EvenInQ read there is False), 20 at (1,3)@N=6 (even, sectors are real
    /// families, EvenInQ True), 21 at (1,2)@N=8 (odd), 35 at (1,4)@N=8 (odd, and fold-fixed, which is
    /// what makes it the case separating gauge from fold).</summary>
    [Fact]
    [Trait("Category", "EVEN_DISC13_SCOUT")]
    public void TheArcsBlocks_CarryTheRecordedCharacters()
    {
        (int n, int p, int q, int character, bool commutes)[] rows =
        {
            (4, 1, 2,  9, false),
            (6, 1, 3, 20, true),
            (8, 1, 2, 21, false),
            (8, 1, 4, 35, false),
        };
        foreach (var (n, p, q, character, commutes) in rows)
        {
            Assert.Equal(character, (p + q) * (n - 1));
            Assert.Equal(commutes, WeightCoherenceBlock.BipartiteGaugeCommutesWithReflection(n, p, q));
            _out.WriteLine($"({p},{q})@N={n}: character {character} -> " +
                           $"{(commutes ? "commutes, sectors conj-closed" : "anticommutes, sectors swapped")}");
        }
    }

    /// <summary>Self-foldedness is NEITHER necessary NOR sufficient for the criterion, which is the whole
    /// reason the cause is the gauge and not the fold. (0,3) and (2,3) at N=6 are fixed by a fold leg
    /// (bra weight 3 = N/2) and still anticommute; (1,3)@N=4 commutes without being fold-fixed at all.</summary>
    [Fact]
    [Trait("Category", "EVEN_DISC13_SCOUT")]
    public void FoldFixedness_AndTheCriterion_AreIndependent()
    {
        // Deliberately the Core predicate and not a local copy: if the repo's definition of fold-fixedness
        // moves, this control must move with it rather than stay green against a stale duplicate.
        static bool FoldFixed(int n, int p, int q) => FoldResultantCertificate.IsFoldFixedBlock(n, p, q);

        Assert.True(FoldFixed(6, 0, 3));
        Assert.False(WeightCoherenceBlock.BipartiteGaugeCommutesWithReflection(6, 0, 3));   // fold-fixed, swapped
        Assert.True(FoldFixed(6, 2, 3));
        Assert.False(WeightCoherenceBlock.BipartiteGaugeCommutesWithReflection(6, 2, 3));   // fold-fixed, swapped
        Assert.False(FoldFixed(4, 1, 3));
        Assert.True(WeightCoherenceBlock.BipartiteGaugeCommutesWithReflection(4, 1, 3));    // closed, not folded
        _out.WriteLine("(0,3)@N=6 and (2,3)@N=6 fold-fixed yet anticommuting; (1,3)@N=4 commuting yet not fold-fixed");
    }

    /// <summary>NEGATIVE CONTROL, the Δ fence. The Δ·ZZ term is a purely IMAGINARY diagonal contribution,
    /// which a diagonal of signs cannot reach, so 𝒟L𝒟 = conj(L) must FAIL at Δ ≠ 0 wherever the ZZ
    /// energies of ket and bra differ. This is PROOF_CODIM1 §7's own Δ = 0 scope on ingredient (iv),
    /// measured rather than asserted.</summary>
    [Fact]
    [Trait("Category", "EVEN_DISC13_SCOUT")]
    public void AtNonzeroDelta_TheGaugeIdentityBreaks()
    {
        const int n = 6, wKet = 1, wBra = 3;
        var q = new Complex(0.7, 0.0);
        var signs = WeightCoherenceBlock.BipartiteGaugeSigns(n, wKet, wBra);
        var L0 = WeightCoherenceBlock.Build(n, wKet, wBra, q, 0.0);
        var Ld = WeightCoherenceBlock.Build(n, wKet, wBra, q, 0.85);

        double worst0 = Worst(L0, signs), worstD = Worst(Ld, signs);
        Assert.Equal(0.0, worst0);
        // The prose says the identity BREAKS rather than degrades, so the gate is the break's exact size,
        // 2·q·Δ·(max |zz(ket) − zz(bra)|) = 2·0.7·0.85·8 = 9.52, and not a token "greater than something small".
        Assert.Equal(9.52, worstD, 10);
        _out.WriteLine($"({wKet},{wBra})@N={n}: Δ=0 residual {worst0} exactly, Δ=0.85 residual {worstD:F6}");

        static double Worst(Complex[,] L, int[] signs)
        {
            double w = 0.0;
            for (int i = 0; i < signs.Length; i++)
                for (int j = 0; j < signs.Length; j++)
                    w = Math.Max(w, (signs[i] * signs[j] * L[i, j] - Complex.Conjugate(L[i, j])).Magnitude);
            return w;
        }
    }

    /// <summary>NEGATIVE CONTROL, the field fence. A longitudinal Z-field is bit-flip-ODD and again enters
    /// as an imaginary diagonal, so it breaks the identity for the same reason Δ does. Together with the Δ
    /// row this fixes what the gauge actually needs: a REAL diagonal and hops that flip the bipartite parity.</summary>
    [Fact]
    [Trait("Category", "EVEN_DISC13_SCOUT")]
    public void UnderALongitudinalField_TheGaugeIdentityBreaks()
    {
        const int n = 5, wKet = 1, wBra = 2;
        var q = new Complex(0.7, 0.0);
        var field = new[] { 0.3, -0.1, 0.45, 0.2, -0.35 };
        var signs = WeightCoherenceBlock.BipartiteGaugeSigns(n, wKet, wBra);
        var L = WeightCoherenceBlock.Build(n, wKet, wBra, q, 0.0, field);

        double worst = 0.0;
        for (int i = 0; i < signs.Length; i++)
            for (int j = 0; j < signs.Length; j++)
                worst = Math.Max(worst, (signs[i] * signs[j] * L[i, j] - Complex.Conjugate(L[i, j])).Magnitude);

        // Again the exact break size, 2·q·(max |fe(ket) − fe(bra)|) = 1.4·2.2 = 3.08, not a token threshold.
        Assert.Equal(3.08, worst, 10);
        _out.WriteLine($"({wKet},{wBra})@N={n} with a site-dependent field: worst residual {worst:F6}");
    }
}
