using System;
using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class CoherenceHorizonWitnessTests
{
    [Theory]
    [InlineData(2, 1.0)]
    [InlineData(3, 1.41421)]   // √2
    [InlineData(4, 1.8785)]
    [InlineData(5, 2.37217)]
    public void Horizon_MatchesTheLiveComputedThreshold(int n, double expectedQStar)
    {
        var w = new CoherenceHorizonWitness();
        double q = w.Horizon(n);                  // computed live by bisecting Symphony.Clock.Omega
        Assert.Equal(expectedQStar, q, 2);        // the precise live Q*(N) at 2 decimals
    }

    /// <summary>The inheritance, exactly: the single-excitation block of H in the site basis IS the
    /// Hückel matrix at α=0, β=J. (J/2)(X_lX_{l+1} + Y_lY_{l+1}) = J(σ⁺_lσ⁻_{l+1} + h.c.), so restricting
    /// the number-conserving chain to the one-excitation sector leaves the tridiagonal hopping matrix and
    /// the entry-wise residual is 0.0, with no tolerance admissible. Run across couplings, not only at the
    /// unit dyadic one that cannot break it: 1/3 and 0.1 are not representable in binary, and the identity
    /// still holds bit for bit because J/2 halves exactly and the XX and YY halves resum to J.</summary>
    [Theory]
    [InlineData(2, 1.0)]
    [InlineData(3, 1.0)]
    [InlineData(4, 1.0)]
    [InlineData(5, 1.0)]
    [InlineData(6, 1.0)]
    [InlineData(3, 0.1)]
    [InlineData(4, 1.0 / 3.0)]
    [InlineData(5, 2.718281828459045)]
    [InlineData(6, 1e-7)]
    public void SingleExcitationBlock_IsExactlyTheHuckelMatrix(int n, double j)
    {
        double residual = CoherenceHorizonWitness.HuckelResidual(n, j);
        Assert.True(residual == 0.0,
            $"N={n}, J={j:R}: the single-excitation block departs from the Hückel matrix by {residual:E3}; " +
            "this residual has an exact route (both sides are built from the same J) and any non-zero " +
            "value is a finding about the construction, not a tolerance to widen");
    }

    /// <summary>The control the gate above needs, isolating the term that actually discriminates. The
    /// two HamiltonianType builders differ in their coupling convention as well as in the ZZ, so a bare
    /// XY-against-Heisenberg distance would fire on the halved hop alone and would say nothing about ZZ.
    /// The diagonal does say it: the Hückel matrix and the XY block both carry an exactly zero diagonal at
    /// every coupling, and the ZZ term puts weight there at every N.</summary>
    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void ZzPutsWeightWhereTheHuckelDiagonalIsEmpty(int n)
    {
        double xy = CoherenceHorizonWitness.XyDiagonalWeight(n);
        double zz = CoherenceHorizonWitness.ZzDiagonalWeight(n);
        Assert.True(xy == 0.0, $"N={n}: the XY single-excitation diagonal is {xy:E3}, not the exact 0.0 the identity needs");
        Assert.True(zz > 0.2, $"N={n}: the ZZ term leaves only {zz:E3} on the diagonal, so this control cannot separate it");
    }

    [Fact]
    public void N2_IsTheExceptionalPointBase()
    {
        var w = new CoherenceHorizonWitness();
        Assert.Equal(1.0, w.Horizon(2), 3);  // the EP itself: γ=J, where the ±J band mode ceases to be the gap mode
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void EpModes_CoalescerIsThe02CoherenceWithRigidityToZero(int n)
    {
        var w = new CoherenceHorizonWitness();
        var ep = w.EpModes(n);
        // the coalescing gap mode is a genuine EP: rigidity collapses as √(Q−Q*). Read at Q*(1+δ) the 2×2's own law
        // is r = √(2δ)·(1 − O(δ)) at N=2,3, and the dressed pair from N=4 sits below it (SE-block reading, gated in
        // Coalescer_RigidityObeysTheTwoByTwoLaw). On the full L the coalescer repeats N−1 times, one copy per (k,k)
        // sector with the SE block's rigidity, and the eigensolver's basis of the repeated eigenspace mixes the
        // copies; a mixture of equal-rigidity copies reads at or below that value (Cauchy-Schwarz on the dual
        // pairing), never above. So the per-vector value is basis-dependent, and the law's value plus its O(δ)
        // term bounds it from above in every basis.
        Assert.Equal(n - 1, ep.CoalescerCopies);
        double rLaw = Math.Sqrt(2.0 * CoherenceHorizonWitness.EpReadDelta) * (1.0 + 2.0 * CoherenceHorizonWitness.EpReadDelta);
        Assert.True(ep.Coalescer.Rigidity <= rLaw,
            $"N={n}: coalescer rigidity {ep.Coalescer.Rigidity:F5} exceeds the 2×2 law √(2δ)(1+2δ) = {rLaw:F5}");
        // and it is the {0,2}-coherence by SUPPORT: all of its weight sits on n_diff ∈ {0, 2}. Its WEIGHTS are
        // ½/½ (the floor Re = −2γ) at N=2,3 only; from N=4 the coalescer sits below the floor, so there is no
        // window on the mean here (the old 0.8..1.2 window, read at the handover, passed 1.07 or 1.00, whichever
        // split branch the eigensolver handed it, and called the drift noise).
        double support = ep.CoalescerHist.GetValueOrDefault(0) + ep.CoalescerHist.GetValueOrDefault(2);
        Assert.True(support > 1.0 - 1e-9, $"N={n}: coalescer weight outside n_diff ∈ {{0,2}}: {1.0 - support:E2}");
        // the Absorption Theorem on the full-L eigenvector: mean n_diff = Re λ/(−2γ), exactly (L_H anti-Hermitian,
        // D real diagonal). Read at Q = Q*(1+δ), so γ = J/Q. The residual is the 4^N-dim eigensolver's on a
        // near-defective pair (K ~ 1e3 ⟹ ~1e-11); the gate sits three decades above it.
        double gamma = CoherenceHorizonWitness.J / ep.ReadAtQ;
        double absorption = ep.CoalescerMeanNDiff - ep.Coalescer.Lambda.Real / (-2.0 * gamma);
        Assert.True(Math.Abs(absorption) < 1e-8, $"N={n}: mean n_diff − Re λ/(−2γ) = {absorption:E2}");
        // cross-dock: the full-L coalescer's n_diff = 2 weight IS the SE-block coalescer's w2 (the (1,1) block is
        // an exact sub-block, the same eigenvector read twice), to the two eigensolvers' precision.
        double w2Se = w.CoalescerAtEp(n).W2;
        Assert.True(Math.Abs(ep.CoalescerHist.GetValueOrDefault(2) - w2Se) < 1e-8,
            $"N={n}: full-L hist[2] = {ep.CoalescerHist.GetValueOrDefault(2):F10} vs SE-block w2 = {w2Se:F10}");
    }

    // ---- the excess light: the EP is on the floor Re = −2γ only where the coalescer's coherence share is ½ ----

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Coalescer_ObeysTheAbsorptionTheoremExactly(int n)
    {
        var w = new CoherenceHorizonWitness();
        var c = w.CoalescerAtEp(n);
        // Re λ = −2γ⟨n_diff⟩ for every right eigenvector (L_H anti-Hermitian, D real diagonal), so in floor
        // units Re/(−2γ) = 2·w2 exactly on the SE block (n_diff = 0 on populations, 2 on coherences). Through an
        // eigensolver the identity holds to the solver's own residual: |Re λ/(−2γ) − ⟨n_diff⟩| ≤ ‖Lv−λv‖/(2γ‖v‖)
        // (v†Lv = λ‖v‖² + v†(Lv−λv), Re v†L_Hv = 0), plus the rounding of the quotient itself, dim·2·eps in floor
        // units. That is the law gated here, not a chosen threshold; a weight outside {0,2} or a wrong index map
        // breaks it by O(1).
        double g = CoherenceHorizonWitness.J / (w.EpQ(n) * (1.0 + CoherenceHorizonWitness.EpReadDelta));
        double bound = c.EigenResidual / (2.0 * g) + n * n * 2.0 * 2.220446049250313e-16;
        Assert.True(Math.Abs(c.AbsorptionResidual) <= bound,
            $"N={n}: Re/(−2γ) − 2w2 = {c.AbsorptionResidual:E2} exceeds the eigensolver bound {bound:E2}");
        Assert.True(c.AbsIm > 1e-6, $"N={n}: the picker did not resolve the coalescer branch (|Im| = {c.AbsIm:E2})");
    }

    [Fact]
    public void TwoByTwoClosure_HoldsExactlyAtN2AndN3_AndFailsFromN4()
    {
        // The exact route to "w2 = ½ at N=2,3": P = |0⟩⟨0| − |N−1⟩⟨N−1| and C = [h,P]/J close under h,
        // [h,C] = c₂·J·P, entry-wise in integer arithmetic. Compared at 0.0: a nonzero here is a construction error.
        Assert.True(CoherenceHorizonWitness.ClosureResidual(2) == 0.0, "N=2 closure residual is not 0.0");
        Assert.True(CoherenceHorizonWitness.ClosureResidual(3) == 0.0, "N=3 closure residual is not 0.0");
        // From N=4 the leak lands on interior sites: max|[h,C] − c₂JP| = 2 exactly (two unit entries), at every N.
        for (int n = 4; n <= 12; n++)
            Assert.True(CoherenceHorizonWitness.ClosureResidual(n) == 2.0, $"N={n}: closure residual {CoherenceHorizonWitness.ClosureResidual(n)} != 2.0");
    }

    [Fact]
    public void Coalescer_IsOnTheFloorOnlyAtN2AndN3()
    {
        var w = new CoherenceHorizonWitness();
        // N=2,3: the pair is the invariant 2×2 (TwoByTwoClosure above, the exact route), which forces w2 = ½ and puts
        // the EP on the floor. Through the eigensolver: on the 2×2 the complex pair has Re λ = −2γ exactly (trace −4γ),
        // so w2 − ½ = ½·(2w2 − Re/(−2γ)) + ½·(Re/(−2γ) − 1); the first term is the identity's residual, ≤ the backward
        // error ‖Lv−λv‖/‖v‖ in floor units, the second the computed eigenvalue's error, ≤ that backward error times its
        // condition number 1/r (Wilkinson; r the phase rigidity, ~22 here), plus the quotient's rounding. A law of the
        // instrument, not a chosen number; the exact statement is the closure above.
        foreach (int n in new[] { 2, 3 })
        {
            var c = w.CoalescerAtEp(n);
            double g = CoherenceHorizonWitness.J / (w.EpQ(n) * (1.0 + CoherenceHorizonWitness.EpReadDelta));
            double bound = 0.5 * c.EigenResidual / (2.0 * g) * (1.0 + 1.0 / c.Rigidity) + n * n * 2.0 * 2.220446049250313e-16;
            Assert.True(Math.Abs(c.W2 - 0.5) <= bound,
                $"N={n}: w2 − ½ = {c.W2 - 0.5:E2} exceeds the eigensolver bound {bound:E2}");
        }
        // N≥4: strictly below the floor at Q*(1+δ), δ = 0.001, pinned from the independent Python census
        // (simulations/coherence_horizon_se_block.py L_se under numpy eig, 2026-09-05). The sixth decimal moves
        // with δ (1.4e-5 at N=4, 1.2e-4 at N=8 between δ = 0.001 and 0.002), so this pins the value AT this δ, at 5e-6.
        var census = new Dictionary<int, double> { { 4, 0.507237 }, { 5, 0.517731 }, { 6, 0.529342 }, { 8, 0.552375 } };
        foreach (var (n, w2) in census)
            Assert.True(Math.Abs(w.CoalescerAtEp(n).W2 - w2) < 5e-6, $"N={n}: w2 = {w.CoalescerAtEp(n).W2:F6}, census {w2:F6}");
        // and the strain grows with N: the coalescer's rate −4γ·w2 walks from −2γ toward the −4γ centre of the
        // resummed dispersion λ² + 8γλ + 4J²q² (PROOF_COHERENCE_HORIZON_SLOPE).
        for (int n = 3; n < CoherenceHorizonWitness.SeBlockMaxN; n++)
            Assert.True(w.CoalescerAtEp(n).W2 < w.CoalescerAtEp(n + 1).W2, $"w2 not increasing at N={n}");
    }

    [Fact]
    public void Coalescer_RigidityObeysTheTwoByTwoLaw()
    {
        var w = new CoherenceHorizonWitness();
        double sqrt2d = Math.Sqrt(2.0 * CoherenceHorizonWitness.EpReadDelta);
        // N=2,3: on the 2×2 (λ² + 4γλ + c₂J² = 0) the phase rigidity at Q*(1+δ) is √(2δ)·(1 − O(δ)); the O(δ) term
        // is the law's own next order (0.75δ at δ = 0.001), gated at 2δ.
        foreach (int n in new[] { 2, 3 })
            Assert.True(Math.Abs(w.CoalescerAtEp(n).Rigidity / sqrt2d - 1.0) <= 2.0 * CoherenceHorizonWitness.EpReadDelta,
                $"N={n}: r/√(2δ) = {w.CoalescerAtEp(n).Rigidity / sqrt2d:F5} is not 1 − O(δ)");
        // N≥4: the dressed pair sits below the 2×2 law, and the deficit grows with N (0.987, 0.968, 0.945, 0.921,
        // 0.897 at N=4..8), like the excess light: read, and pinned as monotone.
        double prev = 1.0;
        for (int n = 4; n <= CoherenceHorizonWitness.SeBlockMaxN; n++)
        {
            double ratio = w.CoalescerAtEp(n).Rigidity / sqrt2d;
            Assert.True(ratio < prev, $"N={n}: r/√(2δ) = {ratio:F5} is not below N={n - 1}'s {prev:F5}");
            prev = ratio;
        }
    }

    [Fact]
    public void EpVerdict_ReadsTheSquareRootLawAcrossTwoDecades()
    {
        // The verdict's law, on the SE block where the coalescer eigenvalue is simple: r = k·√(2δ)·(1 + aδ + O(δ²))
        // (no √δ term: the two branches are a conjugate pair and share r), so the exponent between two offsets a decade
        // apart is p = ½ + a·(δ_a − δ_b)/ln 10 + O(δ²). Read at 10δ₁, δ₁ and δ₁/10, the coarse deviation is ten times the
        // fine one to that order, and the verdict removes it by extrapolation. Gated here, at every N the SE block is
        // offered for and at two read offsets a decade apart: the verdict (the extrapolated exponent is ½ within δ₁), and
        // the law it rests on, twice, the coarse deviation over the fine one and the fine deviation at δ₁ = 10⁻³ over the
        // one at 10⁻⁴, each 10 up to the law's own next order, O(10·δ₁) relative (the data read 9.94 to 10.02 at δ₁ = 10⁻³
        // and 9.994 to 10.002 at 10⁻⁴), gated at five times that, 50·δ₁, so the gate shrinks a decade with δ₁ as the law
        // does (for the ratio across the two δ₁ the larger one sets it).
        var w = new CoherenceHorizonWitness();
        for (int n = 2; n <= CoherenceHorizonWitness.SeBlockMaxN; n++)
        {
            var at3 = w.SquareRootLaw(n, 1e-3);
            var at4 = w.SquareRootLaw(n, 1e-4);
            const double lawWidth = 50.0;
            Assert.True(at3.IsSquareRootEp, $"N={n}: p₀ − ½ = {at3.Extrapolated - 0.5:E3} at δ₁ = 1e-3 exceeds δ₁");
            Assert.True(at4.IsSquareRootEp, $"N={n}: p₀ − ½ = {at4.Extrapolated - 0.5:E3} at δ₁ = 1e-4 exceeds δ₁");
            foreach (var (label, ratio, delta1) in new[]
                     {
                         ("coarse over fine at δ₁ = 1e-3", (at3.Coarse - 0.5) / (at3.Exponent - 0.5), 1e-3),
                         ("coarse over fine at δ₁ = 1e-4", (at4.Coarse - 0.5) / (at4.Exponent - 0.5), 1e-4),
                         ("fine at 1e-3 over fine at 1e-4", (at3.Exponent - 0.5) / (at4.Exponent - 0.5), 1e-3),
                     })
                Assert.True(Math.Abs(ratio / 10.0 - 1.0) <= lawWidth * delta1,
                    $"N={n}: {label} is {ratio:F4}, not the law's 10 within 50·δ₁ = {lawWidth * delta1}");
        }
        // the static reading from a given Q* is the same computation as the instance's: identical, not close
        foreach (int n in new[] { 2, 5, 8 })
        {
            var instance = w.SquareRootLaw(n);
            var fromQ = CoherenceHorizonWitness.SquareRootLawAt(n, w.EpQ(n), CoherenceHorizonWitness.EpReadDelta);
            Assert.True(instance == fromQ, $"N={n}: SquareRootLaw {instance} and SquareRootLawAt {fromQ} differ");
        }
    }

    /// <summary>Past the SE-block fence the law's next order has grown (a ≈ −2.7 at N = 15), and the raw fine exponent
    /// misses ½ by more than δ₁ at every δ₁, so a verdict on it alone would call a genuine EP2 "no EP". The extrapolated
    /// verdict reads ½ there, which is what makes it independent of <see cref="CoherenceHorizonWitness.SeBlockMaxN"/>.</summary>
    [Theory]
    [InlineData(15)]
    [InlineData(16)]
    [InlineData(20)]
    public void EpVerdict_HoldsPastTheFence_WhereTheRawExponentMisreads(int n)
    {
        double d1 = CoherenceHorizonWitness.EpReadDelta;
        var law = CoherenceHorizonWitness.SquareRootLawAt(n, EpCharacterWitness.BisectEpQ(n), d1);
        Assert.True(Math.Abs(law.Exponent - 0.5) > d1,
            $"N={n}: the raw exponent misses ½ by only {law.Exponent - 0.5:E3}, so this N does not show why the extrapolation is needed");
        Assert.True(law.IsSquareRootEp, $"N={n}: extrapolated p₀ − ½ = {law.Extrapolated - 0.5:E3} exceeds δ₁");
        double decade = (law.Coarse - 0.5) / (law.Exponent - 0.5);
        Assert.True(Math.Abs(decade / 10.0 - 1.0) <= 50.0 * d1, $"N={n}: coarse over fine is {decade:F4}, not the law's 10 within 50·δ₁");
    }

    [Fact]
    public void EpVerdict_RejectsTheBandEdgeSurvivor()
    {
        // The control the law must reject: the band-edge survivor, read on its own (0,1) block where L = i·h − 2γ is a
        // normal matrix with simple eigenvalues, keeps r = 1 at all three offsets, so its exponent is 0 and not ½.
        var w = new CoherenceHorizonWitness();
        double d1 = CoherenceHorizonWitness.EpReadDelta, decade = CoherenceHorizonWitness.VerdictDecade;
        for (int n = 2; n <= 5; n++)
        {
            var p = CoherenceHorizonWitness.ExponentAcrossDecades(w.BandEdgeSurvivor(n, d1 * decade).Rigidity,
                w.BandEdgeSurvivor(n, d1).Rigidity, w.BandEdgeSurvivor(n, d1 / decade).Rigidity);
            Assert.False(CoherenceHorizonWitness.IsSquareRootLaw(p.Extrapolated, d1),
                $"N={n}: the survivor's extrapolated exponent {p.Extrapolated:E2} was read as the square-root law");
            Assert.True(Math.Abs(p.Extrapolated) < 0.01, $"N={n}: the survivor's exponent {p.Extrapolated:E2} is not the non-coalescing 0");
        }
    }

    /// <summary>The controls a threshold on r alone could misread, each read through the verdict's own three offsets
    /// and extrapolation. An EP3 (the companion matrix of λ³ = δ) coalesces too, with r ∝ δ^(2/3), so its exponent tends
    /// to ⅔: the verdict must reject it and read it on the EP3's side of ½. A non-normal diabolic point
    /// (S·diag(δ, −δ)·S⁻¹, eigenvectors fixed and non-orthogonal) keeps r = 0.6 at every δ; its three readings agree to
    /// the eigensolver's δ-independent error (the gap and the matrix scale together, so the eigenvectors' error is
    /// cond(S)·O(eps) at every offset), gated at 10³·eps. And the most non-normal SE mode that does not coalesce, on the
    /// coalescer's own block at N = 4 and 6, keeps a rigidity well below one that does not move.</summary>
    [Fact]
    public void EpVerdict_RejectsAnEp3ADiabolicPointAndANonCoalescingMode()
    {
        double d1 = CoherenceHorizonWitness.EpReadDelta, d0 = d1 * CoherenceHorizonWitness.VerdictDecade,
               d2 = d1 / CoherenceHorizonWitness.VerdictDecade;
        var ep3 = CoherenceHorizonWitness.ExponentAcrossDecades(CoherenceHorizonWitness.Ep3ControlRigidity(d0),
            CoherenceHorizonWitness.Ep3ControlRigidity(d1), CoherenceHorizonWitness.Ep3ControlRigidity(d2));
        Assert.False(CoherenceHorizonWitness.IsSquareRootLaw(ep3.Extrapolated, d1), $"the EP3 read as an EP2 (p₀ = {ep3.Extrapolated:F4})");
        Assert.True(ep3.Extrapolated - 0.5 > 0.1, $"the EP3 reads p₀ = {ep3.Extrapolated:F4}, not on its ⅔ side of ½");

        var rs = new[] { d0, d1, d2 }.Select(CoherenceHorizonWitness.DiabolicControlRigidity).ToArray();
        Assert.All(rs, r => Assert.True(Math.Abs(r - 0.6) <= 1e3 * 2.220446049250313e-16, $"the diabolic rigidity {r:R} is not 0.6"));
        var diabolic = CoherenceHorizonWitness.ExponentAcrossDecades(rs[0], rs[1], rs[2]);
        Assert.False(CoherenceHorizonWitness.IsSquareRootLaw(diabolic.Extrapolated, d1),
            $"the diabolic point read as an EP (p₀ = {diabolic.Extrapolated:E2})");

        foreach (int n in new[] { 4, 6 })
        {
            var mode = CoherenceHorizonWitness.NonNormalSeModeRigidities(n, EpCharacterWitness.BisectEpQ(n), d1);
            Assert.All(new[] { mode.R0, mode.R1, mode.R2 }, r => Assert.True(r > 0.2 && r < 0.99,
                $"N={n}: the followed mode's rigidity {r:F4} left the non-normal, non-coalescing band"));
            var p = CoherenceHorizonWitness.ExponentAcrossDecades(mode.R0, mode.R1, mode.R2);
            Assert.False(CoherenceHorizonWitness.IsSquareRootLaw(p.Extrapolated, d1),
                $"N={n}: the non-coalescing mode read as an EP (p₀ = {p.Extrapolated:E2})");
            Assert.True(Math.Abs(p.Extrapolated) < 0.01, $"N={n}: the non-coalescing mode's exponent {p.Extrapolated:E2} is not 0");
        }
    }

    /// <summary>The control that exercises the verdict's width: a genuine EP2 read from a Q* misplaced by ε. The offsets
    /// then read Q*(1+ε)(1+δ), so ε = 10⁻⁶ shifts the finest one (10⁻⁴) by a hundredth and moves the extrapolated
    /// exponent by about 2·10⁻³ (from below: −2.14·10⁻³ at N = 2..8, +2.16·10⁻³ for ε = −10⁻⁶), beyond δ₁, so the verdict
    /// must reject it; ε = 10⁻⁹ moves it by ~10⁻⁶ and must pass. The other controls sit ≥ 0.16 from ½ and would be rejected
    /// by any width; this one would not.</summary>
    [Theory]
    [InlineData(2)]
    [InlineData(4)]
    [InlineData(8)]
    public void EpVerdict_RejectsAnEp2ReadFromAMisplacedQStar(int n)
    {
        var w = new CoherenceHorizonWitness();
        double d1 = CoherenceHorizonWitness.EpReadDelta, eps = CoherenceHorizonWitness.MisplacedQStar;
        foreach (double sign in new[] { 1.0, -1.0 })
        {
            var shifted = CoherenceHorizonWitness.SquareRootLawAt(n, w.EpQ(n) * (1.0 + sign * eps), d1);
            Assert.False(shifted.IsSquareRootEp, $"N={n}, ε = {sign * eps:E0}: p₀ − ½ = {shifted.Extrapolated - 0.5:E3} passed the verdict");
            Assert.True(Math.Abs(shifted.Extrapolated - 0.5) < 3.0 * d1,
                $"N={n}: the misplaced read moved p₀ by {shifted.Extrapolated - 0.5:E3}, not the ~2e-3 the expansion gives");
        }
        var nudged = CoherenceHorizonWitness.SquareRootLawAt(n, w.EpQ(n) * (1.0 + 1e-9), d1);
        Assert.True(nudged.IsSquareRootEp, $"N={n}: a Q* off by 10⁻⁹ reads p₀ − ½ = {nudged.Extrapolated - 0.5:E3}");
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void FullL_CoalescerRepeatsOncePerSector(int n)
    {
        // the coalescer eigenvalue appears once in each (k,k) sector with 1 ≤ k ≤ N−1 (the free-fermion quadratic
        // operators live in all of them), N−1 copies; the count is read on a plateau: the same at 1e-11, 1e-8 and
        // 1e-5, between the copies' own spread (6.6e-14 at N = 5) and the next distinct eigenvalue (5e-2 at N = 5)
        var w = new CoherenceHorizonWitness();
        var counts = w.CoalescerCopyCounts(n, 1e-11, CoherenceHorizonWitness.CopyTolerance, 1e-5);
        Assert.All(counts, c => Assert.Equal(n - 1, c));
    }

    [Fact]
    public void Handover_EqualsTheEpExactlyWhereTheEpIsOnTheFloor()
    {
        var w = new CoherenceHorizonWitness();
        // N=2,3: w2 = ½ puts the EP on the floor, so the darker real branch reaches the floor AT the EP: Q_h = Q*.
        // Two bisections on the same 4- and 9-dim block; their agreement is bisection precision (70 and 60 steps).
        foreach (int n in new[] { 2, 3 })
            Assert.True(Math.Abs(w.EpQ(n) - w.HandoverQ(n)) < 1e-9, $"N={n}: Q* − Q_h = {w.EpQ(n) - w.HandoverQ(n):E2}");
        // N≥4: Q_h < Q* by the trace dressing ((2w2−1)/c)² to leading order; the value pinned from the
        // independent Python bisection of the same block (2026-09-05), at 1e-3 relative.
        var census = new Dictionary<int, double> { { 4, 1.969e-4 }, { 5, 1.496e-3 }, { 6, 4.972e-3 }, { 7, 1.143e-2 }, { 8, 2.132e-2 } };
        foreach (var (n, gap) in census)
        {
            double d = w.EpQ(n) - w.HandoverQ(n);
            Assert.True(Math.Abs(d - gap) / gap < 1e-3, $"N={n}: Q* − Q_h = {d:E4}, census {gap:E4}");
        }
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void Horizon_IsTheHandoverReadOnTheFullLiouvillian(int n)
    {
        var w = new CoherenceHorizonWitness();
        // Horizon(n) bisects the full L on the same criterion as the SE-block EP (the gap mode stops
        // oscillating), but the full L holds the (0,1) survivor at exactly −2γ, so between Q_h and Q* the split
        // pair's darker branch is the gap and the bisection lands on the HANDOVER. Agreement with the SE-block Q_h is
        // gated at the full-L bisection's own resolution (its γ bracket over 2^steps, mapped to Q): a law of the
        // instrument, not a chosen number.
        double res = w.HorizonResolution(n);
        Assert.True(Math.Abs(w.Horizon(n) - w.HandoverQ(n)) <= res,
            $"N={n}: Horizon {w.Horizon(n):F7} vs Q_h {w.HandoverQ(n):F7}, resolution {res:E1}");
        // from N=4 the EP sits resolvably ABOVE the handover; at N=2,3 the two coincide within the same resolution.
        double epAboveHandover = w.EpQ(n) - w.Horizon(n);
        if (n >= 4)
            Assert.True(epAboveHandover > 2.0 * res, $"N={n}: Q* − Horizon = {epAboveHandover:E2} is not above 2·res = {2.0 * res:E1}");
        else
            Assert.True(Math.Abs(epAboveHandover) <= res, $"N={n}: Q* − Horizon = {epAboveHandover:E2} exceeds res = {res:E1}");
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void EpModes_BandEdgeIsTheCoLocatedSurvivor(int n)
    {
        var w = new CoherenceHorizonWitness();
        var ep = w.EpModes(n);
        double bandIm = 2.0 * System.Math.Cos(System.Math.PI / (n + 1));
        // the band edge oscillates at 2cos(π/(N+1)) in the full L's gap band ...
        Assert.Equal(bandIm, System.Math.Abs(ep.BandEdge.Lambda.Imaginary), 2);
        // ... and does NOT coalesce. Its rigidity is read on its own (0,1) block, where L = i·h − 2γ is normal with
        // simple eigenvalues: r = 1 up to the rounding of the eigenvector norms and the inverse, a few ulps per
        // entry, gated at 64·N·eps. The full-L value (ep.BandEdgeR) reads one basis of a 2N-fold eigenvalue and is
        // not gated.
        var survivor = w.BandEdgeSurvivor(n);
        Assert.Equal(bandIm, System.Math.Abs(survivor.Lambda.Imaginary), 12);
        Assert.True(System.Math.Abs(1.0 - survivor.Rigidity) <= 64.0 * n * 2.220446049250313e-16,
            $"N={n}: the survivor's rigidity on its normal block is {survivor.Rigidity:R}, not 1");
    }

    [Fact]
    public void SqrtScalingRatio_IsConstant_AtN4_GenuineSecondOrderEp()
    {
        // Near a 2nd-order EP, Im² ∝ (Q−Q*): the ratio Im²/(Q−Q*) is the same at two offsets.
        var w = new CoherenceHorizonWitness();
        double r1 = w.SqrtScalingRatio(4, 0.03);
        double r2 = w.SqrtScalingRatio(4, 0.06);
        Assert.True(r1 > 0 && r2 > 0, $"ratios should be positive, got {r1}, {r2}");
        Assert.True(System.Math.Abs(r1 - r2) / r1 < 0.30,
            $"√-scaling ratio should be ~constant (2nd-order EP); got {r1:F3} vs {r2:F3}");
    }
}
