using System;
using System.Linq;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The step-3 shell census engine, gated at N=5 where the full answer is known: the diamond
/// membership map (containment corollary members read σ_min ≈ 0 at their shift, everything else on the
/// probed strip excludes), the two-stage seed refinement (q* then λ_A, pushing the member-noise floor to
/// the EP pair gap), the in-window σ_min exclusion at the window-silent locus-2 (1,1) block (the datum
/// class the whole step exists for), and the Bendixson sanity invariant. The N=9 scout anchors live in
/// the SLOW_SHELLCENSUS facts; this class is the fast gate.</summary>
public class SectorShellCensusTests
{
    [Fact]
    [Trait("Category", "SHELLCENSUS")]
    public void ExpectedMembers_MatchTheContainmentCorollary()
    {
        Assert.Equal(new[]
        {
            (1, 2, "lambdaA"), (2, 3, "lambdaA"), (1, 3, "mu"), (2, 2, "mu"),
        }.ToHashSet(), SectorShellCensus.ExpectedMembers(5));

        Assert.Equal(new[]
        {
            (1, 2, "lambdaA"), (2, 3, "lambdaA"), (3, 4, "lambdaA"), (4, 5, "lambdaA"),
            (1, 7, "mu"), (2, 6, "mu"), (3, 5, "mu"), (4, 4, "mu"),
        }.ToHashSet(), SectorShellCensus.ExpectedMembers(9));
    }

    [Fact]
    [Trait("Category", "SHELLCENSUS")]
    public void N5_SeedRefinement_TightensLambdaA()
    {
        var seed = RealDefectiveSeeds.ForN(5).Single(s => Math.Abs(s.QStar - 0.620878) < 1e-6);
        var (qRefined, lambda, pairGap) = SectorShellCensus.RefineSeed(seed);
        Assert.True(Math.Abs(qRefined - seed.QStar) < 5e-4);      // stays local to the recorded locus
        Assert.True(Math.Abs(lambda.Real - (-4.618886)) < 5e-4);  // near the recorded 6-digit value
        Assert.True(pairGap < 4.8e-4);                            // at or below the recorded gap 4.771e-4
    }

    [Fact]
    [Trait("Category", "SHELLCENSUS")]
    public void N5_Locus1_MembershipMap_MatchesTheDiamond()
    {
        var seed = RealDefectiveSeeds.ForN(5).Single(s => Math.Abs(s.QStar - 0.620878) < 1e-6);
        var result = SectorShellCensus.Run(seed, new SectorShellCensus.Options());
        Assert.True(result.SeedUsable);

        var members = result.Entries.Where(e => e.Probed && e.SigmaMin < result.MemberTol)
            .Select(e => (e.P, e.W, e.Shift)).ToHashSet();
        Assert.Equal(SectorShellCensus.ExpectedMembers(5), members);

        // every probed non-member clears the exclusion threshold
        double worstNonMember = result.Entries
            .Where(e => e.Probed && !members.Contains((e.P, e.W, e.Shift)))
            .Min(e => e.SigmaMin);
        Assert.True(worstNonMember > 1e-2);

        // Bendixson SANITY invariant (near-vacuous as an accuracy check — the estimate is an upper
        // bound by construction; this catches only gross LU collapse) + convergence surfaced
        Assert.All(result.Entries.Where(e => e.Probed), e =>
        {
            Assert.True(e.SigmaMin >= e.WindowMargin - 1e-9);
            Assert.True(e.Converged);
        });

        // THE R-PARITY ALTERNATION LAW (measured by this instrument on first light, 2026-07-04, then
        // derived): R c_l R = P_Z · c_{N−1−l} (the JW string reflected from the other end), and on a
        // (p,w) block P_Z ρ P_Z = (−1)^{p+w} ρ, so R·W = (−1)^{p+w}·W·R. Every BAND step (p,p+1) has
        // p+w = 2p+1 odd ⟹ the carried value's R-parity FLIPS at every band step; the fold preserves
        // it. Member with band source index p carries the seed's parity for ODD p, the opposite for
        // EVEN p. (Diagonal climbs p+w even commute strictly — the ledger's corner-sector argument
        // is untouched.) Seed here is R-even: (1,2) even, (2,3) ODD, (1,3)=fold(1,2) even,
        // (2,2)=fold(2,3) ODD.
        foreach (var e in result.Entries.Where(e => members.Contains((e.P, e.W, e.Shift))))
        {
            // band source index: lambdaA members ARE band blocks (p,p+1); mu members are fold images
            // (p, N−p−1) with the band p preserved as the entry's P (and at odd N the transpose-
            // normalized alternative W = N−1−P has the SAME parity, so P is unambiguous)
            int pSrc = e.P;
            bool seedParityCarried = pSrc % 2 == 1;   // odd band index: seed parity; even: flipped
            bool expectEvenSector = seedParityCarried; // the seed is R-even
            if (expectEvenSector) Assert.True(e.SigmaMinEven <= e.SigmaMinOdd, $"({e.P},{e.W})x{e.Shift}: expected EVEN carrier");
            else Assert.True(e.SigmaMinOdd <= e.SigmaMinEven, $"({e.P},{e.W})x{e.Shift}: expected ODD carrier");
        }

        var summary = result.Summarize();
        Assert.Equal("PASS", summary.Verdict);
        Assert.Empty(summary.DeferredMembers);
        Assert.Empty(summary.Ambiguous);
    }

    // The in-window exclusion the whole step exists for, at the cheapest known silent case:
    // N=5 locus 2 (q*=1.077615, lambda_A=-3.7917) where the (1,1) window [-4,0] CONTAINS lambda_A
    // (the window-shell lemma is silent there — the proof's recorded negative control) yet (1,1)
    // is a certified non-member (the R4 interior-four certificates). The probe must supply the
    // exclusion the window cannot.
    [Fact]
    [Trait("Category", "SHELLCENSUS")]
    public void N5_Locus2_InWindowSilentBlock_IsExcludedByProbe()
    {
        var seed = RealDefectiveSeeds.ForN(5).Single(s => Math.Abs(s.QStar - 1.077615) < 1e-6);
        var result = SectorShellCensus.Run(seed, new SectorShellCensus.Options());
        var e11 = result.Entries.Single(e => e.P == 1 && e.W == 1 && e.Shift == "lambdaA");
        Assert.True(e11.Probed);                     // window contains Re lambda_A: no analytic exclusion
        Assert.Equal(0.0, e11.WindowMargin, 12);
        Assert.True(e11.SigmaMin > 1e-2);            // the probe supplies what the window cannot
    }
}
