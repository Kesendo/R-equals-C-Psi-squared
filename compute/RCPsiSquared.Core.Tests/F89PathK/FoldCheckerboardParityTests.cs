using RCPsiSquared.Core.F89PathK;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>The Λ-side SUPPORT parity of a fold-fixed block, gated exactly over ℤ[i]
/// (sideways_spin_ladder arc). Every antiparity this arc records is q-side (Y1, Y2, the pencil
/// reality); the Λ side has only ever been read as VANISHING (ExactOddTaylorAllZero, the
/// M-descent's entry condition), which dies wherever the vanishing fails. This is the weaker
/// law that survives there, and it is a determinant identity, not an antiunitary one.
///
/// THE LAW. Let the block be fixed by a fold leg: ket weight = N/2 OR bra weight = N/2 at even
/// N (PROOF_CODIM1_BY_ADDITIVITY §"the fold lattice", consequence (b) — BOTH legs, f_Q fixing
/// (N/2, q̃) exactly as f_P fixes (p, N/2)). The bare leg is the LINEAR, q-REVERSING
/// rearrangement P·S(q)·P = 2Λ₀ − S(−q) and P commutes with the site reflection, so it acts
/// inside the R-sector. Determinants give
///
///     χ_q(Λ₀ + u) = (−1)^d · χ_{−q}(Λ₀ − u),   d = deg_Λ,   Λ₀ = −2N,
///
/// at EVERY q, real or complex. Writing χ = Σ c_{j,k}·u^j·q^k that forces c_{j,k} = 0 unless
/// j + k ≡ d (mod 2): a checkerboard SUPPORT. The reality half is the separate unconditional
/// q-side antiparity, so the surviving cells are real on even j and purely imaginary on odd j.
/// ε = (−1)^d is FORCED, not measured — the residual is monic, so the cell (d, 0) settles the
/// branch before anything else is read; the gates below cross-check the sign against ResDeg
/// rather than pretending to measure it.
///
/// WHY IT IS NOT THE COMMITTED BOOLEAN EQUIVALENCE. The arc states "Λ-even ⟺ q-even ⟺
/// coefficient-real", and its transformation form ("a Λ-reflection tied to q → −q") is this
/// identity in disguise. But the equivalence is strictly WEAKER: at (1,2)@N=4 all three of its
/// booleans are False, and False ⟺ False ⟺ False implies nothing, while the checkerboard there
/// holds with its odd cells occupied. Conversely q-evenness plus Λ-evenness implies the
/// checkerboard, which is why (1,3)@N=6 is its degenerate side.
///
/// WHAT IT BUYS, and one half of it is immediate. The law says F_res(−u, −q) = ε·F_res(u, q),
/// so the Λ-roots negate under q ↦ −q while the discriminant is invariant under that negation:
/// D = disc_Λ(F_res) is EVEN in q even where the committed q-evenness read is False. Composed
/// with Y1 (conj-coeff(D)(q) = D(−q)) that forces conj-coeff(D) = D, i.e. **D is REAL on every
/// fold-fixed block** — which DERIVES the arc's B3 prior instead of measuring it, and settles
/// the disc-reality question the registry lists as the cheapest open item at (1,4)@N=8 without
/// a prime pass. Consistent with all three committed disc readings: (1,2)@N=4 fold-fixed, Im D
/// ≡ 0; (1,3)@N=6 fold-fixed, disc real; (1,2)@N=6 NOT fold-fixed, disc genuinely complex. The
/// q-evenness half is visible in the committed N=4 CLOSED FORM and was never read as a law:
/// disc(F₈) = const·q²⁴·(3q⁴+q²−1)²·P₂₀(q) with P₂₀(q) = P₁₀(q²), every factor even in q. The
/// other half is NOT a realised saving: an even-in-q D would halve the D-device's interpolation
/// nodes, but no committed device consumes a support-parity fact (they sample dense over
/// q0 = 0…dBound with no q² substitution), so that halving is designed-for, not banked.
///
/// PREMISE-FREE WHERE IT MATTERS: the identity is about det(ΛI − S(q)) and never mentions the
/// AT factor, so it is a theorem on the FULL charpoly. On F_res it additionally needs the AT
/// factor to carry the same parity, this arc's ungated AT-swap premise; both are read here, so
/// the premise is measured rather than assumed.</summary>
public class FoldCheckerboardParityTests
{
    private readonly ITestOutputHelper _out;
    public FoldCheckerboardParityTests(ITestOutputHelper output) => _out = output;

    private void Report(string tag, FoldResultantCertificate.WeightSturmScoutReport r) =>
        _out.WriteLine($"{tag}: sign={r.FoldCheckerboardSign} full={r.FullCharpolyCheckerboardSign} " +
                       $"oddCells={r.CheckerboardOddCellCount} evenCells={r.CheckerboardEvenCellCount} " +
                       $"resDeg={r.ResDeg} (reads: real={r.ResidualIsReal} qEven={r.EvenInQ} " +
                       $"ΛEven={r.EvenInShiftedLambda})");

    /// <summary>The cheap fold-fixed control at ODD character, where the odd cells must be
    /// OCCUPIED: (1,2)@N=4 has bra weight 2 = N/2, character (p+q)(N−1) = 9 odd, and committed
    /// all-False reads, so the checkerboard is the whole of its structure and cannot pass
    /// vacuously.</summary>
    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    [Trait("Category", "EVEN_DISC13_SCOUT")]
    public void N4_FoldFixedOnBra_OddCharacter_CheckerboardHoldsNonVacuously(bool rOdd)
    {
        var r = FoldResultantCertificate.WeightSturmScout(4, 1, 2, rOdd, timingNodes: 0);
        Report($"N4 (1,2) rOdd={rOdd}", r);
        Assert.Equal(1, r.FoldCheckerboardSign);
        Assert.Equal((r.ResDeg & 1) == 0 ? 1 : -1, r.FoldCheckerboardSign);   // ε = (−1)^d, forced
        Assert.Equal(1, r.FullCharpolyCheckerboardSign);   // premise-free on the full charpoly
        Assert.Equal(6, r.CheckerboardOddCellCount);
        Assert.Equal(15, r.CheckerboardEvenCellCount);
    }

    /// <summary>The control that discriminates the LAW from the way it was first mis-stated: at
    /// (2,1)@N=4 the KET weight is 2 = N/2 and the bra weight is 1, so the block is fold-fixed
    /// by the OTHER leg. A first draft of this class claimed "bra weight = N/2" and would pass
    /// its whole suite while being wrong here; §7(b) names both legs and always did.</summary>
    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    [Trait("Category", "EVEN_DISC13_SCOUT")]
    public void N4_FoldFixedOnKet_CheckerboardStillHolds(bool rOdd)
    {
        var r = FoldResultantCertificate.WeightSturmScout(4, 2, 1, rOdd, timingNodes: 0);
        Report($"N4 (2,1) rOdd={rOdd} [ket-fixed]", r);
        Assert.Equal(1, r.FoldCheckerboardSign);
        Assert.Equal(1, r.FullCharpolyCheckerboardSign);
        Assert.True(r.CheckerboardOddCellCount > 0,
            "ket-fixed at odd character must be non-degenerate, like its bra-fixed mirror");
    }

    /// <summary>The fold-fixed control at EVEN character, where the collapse is expected: at
    /// (1,3)@N=6 the gauge acts within the sector, all three committed reads are True, so every
    /// odd cell is EMPTY and the checkerboard holds degenerately. The Λ-evenness the whole
    /// (1,3) pipeline stands on IS this collapse, which is what makes the checkerboard the
    /// general case and the evenness the special one.</summary>
    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    [Trait("Category", "SLOW_EVEN_DISC13")]
    public void N6_FoldFixed_EvenCharacter_CheckerboardCollapses(bool rOdd)
    {
        var r = FoldResultantCertificate.WeightSturmScout(6, 1, 3, rOdd, timingNodes: 0);
        Report($"N6 (1,3) rOdd={rOdd}", r);
        Assert.Equal(1, r.FoldCheckerboardSign);
        Assert.Equal(1, r.FullCharpolyCheckerboardSign);
        Assert.Equal(0, r.CheckerboardOddCellCount);
        Assert.Equal(rOdd ? 325 : 247, r.CheckerboardEvenCellCount);
    }

    /// <summary>The NEGATIVE control, which keeps the law from being a property of the
    /// construction: (1,2)@N=6 is fold-fixed on NEITHER leg (1 ≠ 3, 2 ≠ 3), so the Λ side has
    /// no source and only the unconditional q-side half survives. The read must FAIL, on the
    /// full charpoly too. Its confound with the character is controlled elsewhere in this class:
    /// (1,2)@N=4 shares the odd character and passes.</summary>
    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    [Trait("Category", "EVEN_DISC13_SCOUT")]
    public void N6_NotFoldFixed_CheckerboardMustFail(bool rOdd)
    {
        var r = FoldResultantCertificate.WeightSturmScout(6, 1, 2, rOdd, timingNodes: 0);
        Report($"N6 (1,2) rOdd={rOdd} [NOT fold-fixed]", r);
        Assert.Equal(0, r.FoldCheckerboardSign);
        Assert.Equal(0, r.FullCharpolyCheckerboardSign);
    }
}
