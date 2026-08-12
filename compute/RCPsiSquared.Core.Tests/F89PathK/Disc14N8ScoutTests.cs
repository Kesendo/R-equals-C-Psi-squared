using RCPsiSquared.Core.F89PathK;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>The first (p, 4) block at N = 8 on this arc's exact ℤ[i] path (the builder tests
/// reach (8,3,4)/(8,4,4), ChromaticityTests; nothing on the disc/Sturm path ever has)
/// (sideways_spin_ladder arc; the gap the (1,2)@N=8 theorem's own closing paragraph names:
/// "(1,2)@N=8 absence is not N=8 absence, the (p,4) fold-fixed blocks are unexamined").
/// (1,4)@N=8 is the analogue of the census block (1,3)@N=6: the BRA fold leg f_P fixes exactly
/// the (p, N/2) blocks at even N, and at N = 6 that block carried FIFTEEN real-q defective loci
/// while its (1,2) sibling carried none. Two stages, because the sector dimension is 280 against
/// the committed 112 and no cost for it is known: the cheap pencil gate below (sub-second, it
/// pins the dimension and exercises the all-2-cycle scope and the q-linearity), then the exact
/// ℤ[i] parity scout in its own slow category.
///
/// What is DERIVED before either runs, so that a failure is a design-level finding and not a
/// code bug:
///
/// (i) The fold leg f_P: (p, q̃) ↦ (p, N−q̃) fixes this block (q̃ = 4 = N/2; PROOF_CODIM1_BY_
/// ADDITIVITY §"the fold lattice", consequence (b)). Its BARE form, F89d composed with the
/// pencil reality, is P·L(q)·P = −2N·I − L(−q), and P (the bit complement) commutes with the
/// site reflection R, so this one acts WITHIN each R-sector and gives, per sector,
/// F_res(Λ, q) = F_res(−4N − Λ, −q): a Λ-reflection about Λ₀ = −2N (Λ = 2λ, the ×2-cleared
/// sector operator, so the λ-centre −N doubles) TIED TO q → −q. That is NOT Λ-evenness at
/// fixed q, and reading it as such was this test's first draft error.
///
/// (ii) Unconditionally, A real and C entry-wise imaginary give conj-coeff(F_res)(Λ, q) =
/// F_res(Λ, −q) — the conj-antiparity the (1,2)@N=8 D-device found one level up, as Y1. Hence
/// ResidualIsReal ⟺ EvenInQ, always. Composed with (i), on a fold-fixed block all THREE reads
/// stand or fall together: EvenInShiftedLambda ⟺ EvenInQ ⟺ ResidualIsReal. Note what the gate
/// below can and cannot see: with the prediction all-False the equivalence is satisfied
/// TRIVIALLY here, as any structureless block would satisfy it; its discriminating instance is
/// the committed all-True read at (1,3)@N=6 (Disc13SturmScoutTests.N6_Disc13_SturmScout_
/// Steps0to2). This block gates the CONSISTENCY of the equivalence, not the equivalence.
///
/// (iii) The q-PRESERVING fold leg is P𝒟, and it carries the bipartite gauge 𝒟, which commutes
/// with R only up to a block-global sign — that sign is the sector swap, and it is trivial
/// exactly when the character (p+q)(N−1) is even. Here 35 is ODD, so the q-preserving leg maps
/// the R-even sector onto the R-odd one and licenses nothing within either. The fold-fixed
/// odd-character precedents agree: (0,3)/(2,3)@N=6 are conj-SWAPPED, and the fold-fixed
/// (1,2)@N=4 control — q̃ = 2 = N/2, character 9 odd — has committed EvenInShiftedLambda =
/// False (Disc13SturmScoutTests). So the PREDICTION here is all three False.
///
/// (iv) What that does and does NOT decline, because two different objects are called "real"
/// in this arc and only one of them is at stake above. COEFFICIENT-realness of F_res as a
/// bivariate over ℤ[i] needs the gauge, hence the character, and is what falls here. The
/// coefficient-realness of D = disc_Λ(F_res) is a different statement, supplied at each real q
/// by the SPECTRUM closure λ ↦ −conj(λ) − 2N that the bare fold gives within the sector
/// ([P, R] = 0), the arc's gated B3 prior — on the same AT premise as (ii), the split off the
/// AT factor having to preserve the closure, gated nowhere and empirically true at (1,2)@N=4.
/// That control is exactly this combination: Λ- and q-evenness both gated False (coefficient-
/// realness follows from q-evenness by the unconditional antiparity, logged not gated there),
/// disc real with Im = 0 at every SAMPLED prime up to the lc-divisor cap. So what declines here
/// is the G/M-descent only — no G, no composition identity — and the route that returns is the
/// ℤ one in the (1,2)@N=4 shape, NOT the ℚ(i) split retired at (1,2)@N=8 (retired on general
/// grounds, a ℤ derivative-gcd count replacing it, and CertifyDiscReImGcd is by design the
/// wrong tool on a real disc). Predicted, not measured: the mod-p discReal read at (1,4)@N=8
/// settles it, at a cost still UNPRICED and certainly not cheap — it repeats this test's exact
/// ℤ[i] bivariate build before any prime pass, so its floor is this test's 2 h 9 m — and a
/// COMPLEX answer would be a finding against B3 rather than a detail.</summary>
public class Disc14N8ScoutTests
{
    private readonly ITestOutputHelper _out;
    public Disc14N8ScoutTests(ITestOutputHelper output) => _out = output;

    /// <summary>Stage 1, cheap: the (1,4)@N=8 sector pencil exists and is well-formed. Three
    /// things are gated, all inside the construction: the all-2-cycle scope (ket weight 1 is odd
    /// at even N, so no reflection-fixed basis pair; a fixed point would THROW and would mean the
    /// weight-general path does not reach this block at all), the entry-wise R-commutation
    /// asserted per q₀ node, and the exact q-linearity of the ×2-cleared sector operator (the
    /// pencil throws unless A + 2C equals the q₀ = 2 build). The dimension is the arithmetic
    /// C(8,1)·C(8,4) = 8·70 = 560 split evenly by the 2-cycle basis, which is the claim the
    /// slow scout's atDeg + resDeg has to reproduce.</summary>
    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    [Trait("Category", "EVEN_DISC14_SCOUT")]
    public void N8_Weight14_SectorPencil_WellFormed(bool rOdd)
    {
        var (a, c) = FoldResultantCertificate.WeightSectorPencil(8, 1, 4, rOdd);
        _out.WriteLine($"N8 (1,4) rOdd={rOdd}: sector dim = {a.GetLength(0)} " +
                       $"(block C(8,1)·C(8,4) = 560, halved by the 2-cycle basis)");
        Assert.Equal(280, a.GetLength(0));
        Assert.Equal(280, a.GetLength(1));
        Assert.Equal(280, c.GetLength(0));
        Assert.Equal(280, c.GetLength(1));
        bool cNonzero = false;
        for (int i = 0; i < 280 && !cNonzero; i++)
            for (int j = 0; j < 280; j++)
                if (!c[i, j].Re.IsZero || !c[i, j].Im.IsZero) { cNonzero = true; break; }
        Assert.True(cNonzero, "the q-linear part must be nonzero: a q-free sector would carry no " +
            "coupling and the whole coalescence question would be empty at this block");
    }

    /// <summary>Stage 2: the three exact ℤ[i] structure reads on (1,4)@N=8, R-even, and the
    /// verdict is that all three are False — the (1,3)@N=6 ℤ route does NOT reach this block.
    /// Cost measured 2026-08-12 on this test's green run, 2 h 9 m quiet (against 1m7s at the
    /// (1,2)@N=8 sector dimension 112: ~116× for 2.5× the dimension, the exact ℤ[i] Berkowitz
    /// growing in dimension AND coefficient length), hence its own category; the evenness
    /// declining is what keeps even this cost finite, since the G/disc_M branch never runs.
    ///
    /// R-even only, the R-odd sector being derived by the SAME gauge-swap argument the arc used
    /// as Y2 at (1,2)@N=8 (character 21 odd), here at character 35: 𝒟 maps the R-even sector
    /// onto the R-odd one with L(q) → L(−q), so F_odd(Λ, q) = F_even(Λ, −q) and the three reads,
    /// all invariant under q → −q, are the same in both sectors. The contrast that keeps this
    /// non-vacuous: at (1,3)@N=6 the character is EVEN, there is no swap, and the two parities
    /// indeed carried DIFFERENT residual degrees, 42 and 48. The whole R-odd statement is
    /// CONDITIONAL, not just its degrees: the swap gives charpoly_odd(Λ, q) = charpoly_even(Λ,
    /// −q), and reads AND degrees descend to F_res = charpoly/AT only if 𝒟 maps AT strands to
    /// AT strands — the "AT-swap premise" this arc flagged for gating at (1,2)@N=8 and never
    /// gated. Same shape one level down: (ii)'s descent to F_res needs the AT factor to inherit
    /// the conj-antiparity, gated for the GAUGE at (1,3)@N=6 (the a+bq ↔ a−bq strand pairing)
    /// and not gated here. And it is UNRUN: falsifying it costs a second ~2 h 9 m run.
    ///
    /// One derived consequence sits inside the measured numbers already, ungated for want of an
    /// exported coefficient list: setting Λ = Λ₀ in (i)'s identity gives f(q) = f(−q), so f must
    /// be q-EVEN. deg f = 218 (even) and v_q(f) = 0 are both consistent, and an ODD deg f would
    /// have refuted (i) on the spot. Asserting every odd coefficient of f zero would turn the
    /// trivially-satisfied equivalence gate into a real one at zero compute cost.</summary>
    [Fact(DisplayName = "N=8 (1,4) scout: F_res parity reads exact over ℤ[i], R-even (the first (p,4) block read)")]
    [Trait("Category", "SLOW_EVEN_DISC14")]
    public void N8_Weight14_Scout_SturmScoutParityReads()
    {
        var r = FoldResultantCertificate.WeightSturmScout(8, 1, 4, rOdd: false, timingNodes: 0);
        _out.WriteLine($"N8 (1,4) R-even: real={r.ResidualIsReal} even={r.EvenInShiftedLambda} " +
                       $"qEven={r.EvenInQ} atDeg={r.AtDeg} resDeg={r.ResDeg} " +
                       $"mDeg={r.MDeg} fDeg={r.FDeg} vQf={r.VQf} square={r.FIsPerfectSquare} " +
                       $"uDeg={r.UDeg} uReal={r.UIsReal} compId={r.CompositionIdentityAtNodes}");
        // THE DECISION LINE, logged and deliberately NOT gated: realness decides whether the (1,3)
        // ℤ route reaches this block or whether the ℚ(i) question returns here. Both values are
        // informative and neither is a bug (see the class comment, (iii)).
        _out.WriteLine(r.ResidualIsReal
            ? "DECISION: F_res COEFFICIENT-REAL at (1,4)@N=8 => something supplies coefficient " +
              "realness at ODD character, which no committed block does; the full (1,3)@N=6 " +
              "landing/Sturm route with G and the M-descent applies to this block"
            : "DECISION: F_res COEFFICIENT-COMPLEX at (1,4)@N=8 => coefficient realness needs the " +
              "gauge (as the arc already gates for (1,3)@N=6), so the G/M-descent declines here; " +
              "the sector DISC is a different object and is predicted REAL by the B3 prior, the " +
              "(1,2)@N=4 shape, which the mod-p discReal read settles cheaply");

        // DERIVED, class comment (i)+(ii): on a fold-fixed block the three reads are EQUIVALENT.
        // Read the class comment's caveat before counting this as evidence: with all three False
        // it is satisfied trivially, and the direction that carries content is the committed
        // all-True read at (1,3)@N=6. Here it gates consistency, and a break would be real.
        Assert.True(r.EvenInShiftedLambda == r.EvenInQ && r.EvenInQ == r.ResidualIsReal,
            "the derived equivalence Λ-even ⟺ q-even ⟺ coefficient-real on a fold-fixed block " +
            "failed: that is a finding about the fold lattice or about the pencil's reality " +
            "(A real, C entry-wise imaginary), not about this code");

        // PREDICTED, class comment (iii), and deliberately labelled apart from the derivation
        // above: the character 35 is odd, so the q-preserving leg swaps the sectors and nothing
        // supplies COEFFICIENT realness, the bare fold provably being unable to do it alone (it
        // only reflects Λ composed with q → −q; the spectrum closure it DOES supply is the other
        // object, class comment (iv)). Precedent: every committed fold-fixed ODD-character block
        // behaves this way ((1,2)@N=4 control; (0,3)/(2,3)@N=6 conj-swapped).
        Assert.False(r.EvenInShiftedLambda, "PREDICTED False. A TRUE here would be the first " +
            "fold-fixed odd-character block carrying coefficient realness and would reopen the " +
            "gauge-vs-fold cause question the (1,3)@N=6 arc entry gates");

        // Sentinel reads, not gates: MDeg and CompositionIdentityAtNodes are initialised to
        // −1/false and written only under the evenness, so given the assertion above the else
        // branch is the one taken and both of its reads are definitional. Kept for shape-parity
        // with the N=4 control's
        // declining branch (Disc13SturmScoutTests) so that a design fact does not read as three
        // bugs; a future session should not count them as independent evidence.
        if (r.EvenInShiftedLambda)
        {
            Assert.Equal(r.ResDeg / 2, r.MDeg);
            Assert.True(r.CompositionIdentityAtNodes,
                "disc_Λ(F_res)(q0) = (−4)^m·f(q0)·disc_M(G)(q0)² exact at q0 = 2, 3, sign included");
        }
        else
        {
            Assert.Equal(-1, r.MDeg);
            Assert.False(r.CompositionIdentityAtNodes, "no G, no identity: the sentinel is false/n-a");
        }
        // Measured 2026-08-12, this test's green run. atDeg + resDeg = 280 is definitional
        // (resDeg is the dimension minus atDeg), so the content is the SPLIT: the AT strands
        // carry 56 of the 280, leaving a degree-224 residual — the largest residual any block
        // in this arc has produced, against 80 at (1,2)@N=8 and 42/48 at (1,3)@N=6.
        Assert.Equal(56, r.AtDeg);
        Assert.Equal(224, r.ResDeg);
        // deg f and v_q(f) ARE pinned here, deliberately and against the (1,2)@N=8 practice of
        // logging them only: there the design note's finding was that they are outputs of a
        // DECLINED M-descent and meaningless as N=8 objects. The precedent is committed and is
        // closer than that finding: the fold-fixed (1,2)@N=4 control pins deg f = 8 in exactly
        // this declining branch (Disc13SturmScoutTests). The reason, written rather than
        // done silently: f = F_res(Λ₀, q) at Λ₀ = −2N is well defined whether or not the
        // evenness holds, and on a FOLD-FIXED block Λ₀ is the true fold-line centre (at
        // (1,2)@N=8 it was not, the block not being fold-fixed, and that is what made the
        // numbers there meaningless). So these are reproducible reads of a defined object.
        Assert.Equal(218, r.FDeg);
        Assert.Equal(0, r.VQf);
        Assert.False(r.FIsPerfectSquare, "f carries no square root here; with the evenness gone " +
            "this is a read of f alone, not the (1,3) lc(f)·ũ² = f chain, which needs G");
        Assert.Equal(-1, r.UDeg);   // definitional given the line above (uDeg = isSquare ? … : −1)
    }
}
