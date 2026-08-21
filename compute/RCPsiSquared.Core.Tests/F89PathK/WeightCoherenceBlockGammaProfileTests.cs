using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.F89PathK;
using Xunit;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>The per-site dephasing PROFILE of <see cref="WeightCoherenceBlock"/>.
///
/// <para>What these gates pin, in the order the risks actually run: the uniform route must not move at all (it
/// is load-bearing under a great many committed numbers), the profile must enter the rate and nothing else, and
/// the SITE INDEX must mean here what the rest of that file already means by it.</para>
///
/// <para>Everything here is an EXACT comparison, because an exact route exists in every case: the uniform route
/// is popcount against a sum of exact 1.0's, the profile values are dyadic, and the cross-builder gate runs at
/// H = 0 where each block IS its diagonal. Note which gate carries WHAT:
/// <see cref="TheProfileEntersTheDIAGONALAsMinusTwiceTheSumOverDisagreeingSitesOnly"/> writes the factor −2 as a
/// literal, so a wrong factor or sign WOULD fail it, but it selects the disagreeing sites with the same loop as
/// the implementation, so it could not catch a reversed site index.
/// The one that carries THAT is
/// <see cref="AgainstPerBlockLiouvillianBuilder_TheProfileMustBeREVERSED_AndAPalindromeCannotSeeIt"/>, which
/// compares against an independently written builder and has a negative control.</para></summary>
public class WeightCoherenceBlockGammaProfileTests
{
    // (n, wKet, wBra) blocks small enough to compare entry by entry, and spanning the shapes the repo builds:
    // the (0,1) vacuum edge block, the (1,2) path-k block, a same-weight block, and a cross-fold partner.
    public static IEnumerable<object[]> Blocks() => new[]
    {
        new object[] { 4, 0, 1 }, new object[] { 4, 1, 2 }, new object[] { 4, 2, 2 },
        new object[] { 5, 1, 2 }, new object[] { 5, 1, 3 }, new object[] { 6, 1, 2 },
    };

    // A non-uniform, non-palindromic, dyadic profile whose exponents lie within a few powers of two, so every
    // subset sum below is exact and the gates can use == (dyadic alone would not do it: 2^-60 beside 2.0 rounds).
    // Non-palindromic so that a site REVERSAL changes it
    // (a palindromic profile is invariant under the very error these gates exist to catch).
    private static readonly double[] DyadicTable = { 0.5, 0.25, 1.5, 0.125, 2.0, 0.75, 1.25, 0.0625, 3.0 };

    private static double[] DyadicProfile(int n) =>
        n <= DyadicTable.Length
            ? DyadicTable.Take(n).ToArray()
            : throw new ArgumentOutOfRangeException(nameof(n), $"the dyadic table holds {DyadicTable.Length} sites, not {n}");

    private static double[] Ones(int n) => Enumerable.Repeat(1.0, n).ToArray();

    [Theory]
    [MemberData(nameof(Blocks))]
    public void UniformProfileOfOnes_IsBitIdenticalToTheRouteWithNoProfile(int n, int wKet, int wBra)
    {
        // The whole committed corpus is built by the no-profile route. If the profile route at γ ≡ 1 differed
        // from it by even one ulp, every number downstream would be in question, so this is == and not ≈.
        var q = new Complex(0.7, 0.3);
        foreach (double delta in new[] { 0.0, 0.85 })
        {
            var plain = WeightCoherenceBlock.Build(n, wKet, wBra, q, delta);
            var viaProfile = WeightCoherenceBlock.Build(n, wKet, wBra, q, delta, null, Ones(n));
            int d = plain.GetLength(0);
            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                    Assert.True(plain[i, j] == viaProfile[i, j],
                        $"n={n} ({wKet},{wBra}) Δ={delta} entry [{i},{j}]: {plain[i, j]} vs {viaProfile[i, j]}");
        }
    }

    [Theory]
    [MemberData(nameof(Blocks))]
    public void TheProfileEntersTheDIAGONALAsMinusTwiceTheSumOverDisagreeingSitesOnly(int n, int wKet, int wBra)
    {
        // Isolate the dissipator: q = 0 kills the hopping and the Δ·ZZ frequency, so the block IS its
        // diagonal and the comparison needs no eigensolver and no subtraction of a Hamiltonian part.
        var gamma = DyadicProfile(n);
        var l = WeightCoherenceBlock.Build(n, wKet, wBra, Complex.Zero, 0.0, null, gamma);
        var kets = WeightCoherenceBlock.Configs(n, wKet);
        var bras = WeightCoherenceBlock.Configs(n, wBra);
        int col = 0;
        foreach (int kc in kets)
            foreach (int bc in bras)
            {
                double expected = 0.0;
                for (int s = 0; s < n; s++)
                    if ((((kc ^ bc) >> s) & 1) != 0) expected += gamma[s];       // site s is BIT s here
                Assert.True(l[col, col] == new Complex(-2.0 * expected, 0.0),
                    $"n={n} ({wKet},{wBra}) |{kc}⟩⟨{bc}|: {l[col, col]} vs {-2.0 * expected}");
                col++;
            }
    }

    [Theory]
    [MemberData(nameof(Blocks))]
    public void TheProfileNeverTouchesTheHopping(int n, int wKet, int wBra)
    {
        // Z-dephasing is diagonal in the computational basis, so a profile may move the diagonal and nothing
        // else. This is what keeps the q-linear split L(q) = A + q·C intact with the profile confined to A.
        var q = new Complex(0.7, 0.3);
        var uniform = WeightCoherenceBlock.Build(n, wKet, wBra, q, 0.85);
        var profiled = WeightCoherenceBlock.Build(n, wKet, wBra, q, 0.85, null, DyadicProfile(n));
        int d = uniform.GetLength(0);
        bool anyDiagonalMoved = false;
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                if (i != j)
                    Assert.True(uniform[i, j] == profiled[i, j], $"off-diagonal [{i},{j}] moved under the profile");
                else if (uniform[i, i] != profiled[i, i]) anyDiagonalMoved = true;
        // A checker that reports nothing has proved nothing: the negative half above is only meaningful if the
        // profile reached the block at all.
        Assert.True(anyDiagonalMoved, $"n={n} ({wKet},{wBra}): the profile changed no diagonal entry");
    }

    [Theory]
    [InlineData(5, 0, 1)]
    [InlineData(5, 1, 2)]
    public void GammaAndTheLongitudinalFieldIndexSitesTheSameWay(int n, int wKet, int wBra)
    {
        // The claim in the class doc is that gammaPerSite[s] and field[s] mean the SAME site s. Two per-site
        // knobs in one builder disagreeing about site order is silent, so it is gated rather than asserted.
        // Method: put each knob on one site alone and check they select the same basis columns.
        var selectedPerSite = new HashSet<string>();
        for (int s = 0; s < n; s++)
        {
            var oneSiteGamma = new double[n];
            oneSiteGamma[s] = 1.0;
            var oneSiteField = new double[n];
            oneSiteField[s] = 1.0;

            // The zero profile is what the GAMMA side needs: touchedByGamma is an absolute "is this entry
            // nonzero" test, so under the uniform default every disagreeing column would be nonzero for reasons
            // that have nothing to do with the one site under test. The FIELD side is a difference against a
            // reference and would survive either way; it takes the zero profile only to keep the two halves
            // symmetric.
            var zeroGamma = new double[n];
            var byGamma = WeightCoherenceBlock.Build(n, wKet, wBra, Complex.Zero, 0.0, null, oneSiteGamma);
            var byField = WeightCoherenceBlock.Build(n, wKet, wBra, Complex.One, 0.0, oneSiteField, zeroGamma);
            var reference = WeightCoherenceBlock.Build(n, wKet, wBra, Complex.One, 0.0, null, zeroGamma);

            int d = byGamma.GetLength(0);
            var touchedByGamma = new List<int>();
            var touchedByField = new List<int>();
            for (int i = 0; i < d; i++)
            {
                if (byGamma[i, i] != Complex.Zero) touchedByGamma.Add(i);
                if (byField[i, i] != reference[i, i]) touchedByField.Add(i);
            }
            Assert.True(touchedByGamma.Count > 0, $"site {s}: the one-site profile selected no coherence at all");
            Assert.Equal(touchedByField, touchedByGamma);
            selectedPerSite.Add(string.Join(",", touchedByGamma));
        }
        // Without this, an implementation that ignored s entirely and always used site 0 would satisfy every
        // assertion above, since both knobs would be wrong together and the comparison is relative.
        Assert.Equal(n, selectedPerSite.Count);
    }

    [Theory]
    [InlineData(4, 0, 1)]
    [InlineData(4, 1, 2)]
    [InlineData(5, 1, 2)]
    public void AgainstPerBlockLiouvillianBuilder_TheProfileMustBeREVERSED_AndAPalindromeCannotSeeIt(
        int n, int wKet, int wBra)
    {
        // The cross-builder convention gate, run at H = 0 so that BOTH blocks are nothing but their dephasing
        // diagonals: no Hamiltonian convention, no eigensolver, no tolerance, and the ONLY thing being compared
        // is which site each builder thinks a rate belongs to.
        //
        // WeightCoherenceBlock: site s = bit s. PerBlockLiouvillianBuilder: site l = bit N−1−l. So the same
        // physical profile is the reversed array there, and handing it over unreversed is a silent, spectrum-
        // preserving mislabel whenever H is reflection-symmetric (which the chain is). That is why it needs a
        // gate and not a comment.
        int d = 1 << n;
        ComplexMatrix zeroH = Matrix<Complex>.Build.Dense(d, d);
        var gamma = DyadicProfile(n);
        var reversed = gamma.Reverse().ToArray();

        var kets = WeightCoherenceBlock.Configs(n, wKet);
        var bras = WeightCoherenceBlock.Configs(n, wBra);
        var flat = (from kc in kets from bc in bras select kc * d + bc).ToList();

        var mine = WeightCoherenceBlock.Build(n, wKet, wBra, Complex.Zero, 0.0, null, gamma);
        var theirs = PerBlockLiouvillianBuilder.BuildBlockZ(zeroH, reversed, flat);
        for (int i = 0; i < flat.Count; i++)
            Assert.True(mine[i, i] == theirs[i, i],
                $"n={n} ({wKet},{wBra}) row {i}: reversed profile disagreed, {mine[i, i]} vs {theirs[i, i]}");

        // Negative control, so the gate above cannot pass by accident: handing the SAME array over unreversed
        // must disagree somewhere. It has to, because DyadicProfile is deliberately non-palindromic.
        var unreversed = PerBlockLiouvillianBuilder.BuildBlockZ(zeroH, gamma, flat);
        bool disagrees = Enumerable.Range(0, flat.Count).Any(i => mine[i, i] != unreversed[i, i]);
        Assert.True(disagrees, $"n={n} ({wKet},{wBra}): the unreversed profile agreed, so this gate proves nothing");

        // And the reason the mislabel can hide: at uniform γ the reversal is invisible, because a constant
        // array is its own reverse. This is not a nicety, it is why every existing consumer passed.
        var flatMine = WeightCoherenceBlock.Build(n, wKet, wBra, Complex.Zero, 0.0, null, Ones(n));
        var flatTheirs = PerBlockLiouvillianBuilder.BuildBlockZ(zeroH, Ones(n), flat);
        for (int i = 0; i < flat.Count; i++)
            Assert.True(flatMine[i, i] == flatTheirs[i, i], $"uniform γ disagreed at row {i}");
    }

    [Theory]
    [InlineData(4, 1.25, -0.25, 0.5)]
    [InlineData(4, 0.75, 1.5, 0.0)]
    [InlineData(4, 3.5, -0.125, 2.0)]
    [InlineData(5, 1.25, -0.25, 0.5)]
    [InlineData(5, 0.0625, 0.5, 0.25)]
    public void TheF89dCrossFoldSimilarityCarriesTheConstant2SumGamma_NotTwoN(
        int n, double qRe, double qIm, double delta)
    {
        // The F89d cross-fold ANTIUNITARY similarity, L(wk, n−wb)(q̄) = −P·conj(L(wk,wb)(q))·Pᵀ − 2σ·I with
        // σ = Σ_s γ_s. What this gate adds is only the CONSTANT: the identity is pinned at uniform γ by
        // WeightCoherenceBlockTests, where the constant reads 2N, and 2N turns out to be 2σ wearing γ ≡ 1. The
        // entry-level reason is one line: the bra bit-flip sends the disagreement sum S(a,b) to σ − S(a,b), so
        // the diagonal reflects about −σ. Until the profile parameter existed σ could only BE N here, so the
        // two constants could not be told apart.
        //
        // Do NOT read this as F1's palindrome. Both carry a 2σ shift and that is where the resemblance ends:
        // this map is antilinear (entry-wise conjugation and q → q̄) and pairs (wk, wb) with (wk, n−wb), while
        // F1's Π is linear; the linear form of this comparison fails by O(1), not by a rounding, at every Δ.
        // PROOF_F1_NONUNIFORM_GAMMA.md is also narrower than its title suggests: its theorem is that the
        // DISSIPATOR block M_D vanishes for any profile, not that the full M does.
        //
        // EXACT, not toleranced: with dyadic q and Δ both sides are built from exactly representable pieces and
        // the residual is 0.0 on the nose. That is a property of dyadic inputs and not of one lucky pair, which
        // is why the rows vary q and Δ rather than only n (measured 0.0 across a wider dyadic sweep). With a
        // NON-dyadic q = 1.3 − 0.2i, Δ = 0.6 it is 9e-16 instead: a fact about those decimals, not the identity.
        var gamma = DyadicProfile(n);
        double shift = 2.0 * gamma.Sum();
        var q = new Complex(qRe, qIm);
        for (int wKet = 0; wKet <= n; wKet++)
            for (int wBra = 0; wBra <= n; wBra++)
            {
                var a = WeightCoherenceBlock.Build(n, wKet, wBra, q, delta, null, gamma);
                var partner = WeightCoherenceBlock.Build(n, wKet, n - wBra, Complex.Conjugate(q), delta, null, gamma);
                var perm = WeightCoherenceBlock.BraComplementPermutation(n, wKet, wBra);
                for (int t = 0; t < perm.Length; t++)
                    for (int u = 0; u < perm.Length; u++)
                    {
                        Complex expected = -Complex.Conjugate(a[t, u]) - (t == u ? new Complex(shift, 0) : Complex.Zero);
                        Assert.True(partner[perm[t], perm[u]] == expected,
                            $"n={n} ({wKet},{wBra}) [{t},{u}]: got {partner[perm[t], perm[u]]}, want {expected}");
                    }
            }

        // Negative control, exact: swapping in the uniform constant 2N must shift every diagonal comparison by
        // exactly |2σ − 2N|. Asserted as an equality, not as "> some epsilon": the deviation has an exact route
        // (both constants are dyadic here), so a threshold would be an error code rather than a measurement.
        var a01 = WeightCoherenceBlock.Build(n, 0, 1, q, delta, null, gamma);
        var p01 = WeightCoherenceBlock.Build(n, 0, n - 1, Complex.Conjugate(q), delta, null, gamma);
        var perm01 = WeightCoherenceBlock.BraComplementPermutation(n, 0, 1);
        double gap = Math.Abs(shift - 2.0 * n);
        Assert.True(gap > 0, $"n={n}: the profile sums to N, so the two constants coincide and nothing is tested");
        for (int t = 0; t < perm01.Length; t++)
        {
            Complex expected = -Complex.Conjugate(a01[t, t]) - new Complex(2.0 * n, 0);
            Assert.True((p01[perm01[t], perm01[t]] - expected).Magnitude == gap,
                $"n={n} row {t}: the uniform shift 2N deviated by something other than |2σ − 2N| = {gap}");
        }
    }

    [Fact]
    public void TheFieldAndTheProfileCompose()
    {
        // Claimed by the field overload's doc, so gated: both knobs nonzero at once, the field entering the
        // frequency and the profile the rate, neither disturbing the other.
        const int n = 5;
        var gamma = DyadicProfile(n);
        var field = new[] { 0.5, -0.25, 0.75, 0.125, -1.0 };
        var q = new Complex(1.25, 0.0);
        var both = WeightCoherenceBlock.Build(n, 1, 2, q, 0.0, field, gamma);
        var profileOnly = WeightCoherenceBlock.Build(n, 1, 2, q, 0.0, null, gamma);
        var fieldOnly = WeightCoherenceBlock.Build(n, 1, 2, q, 0.0, field, new double[n]);
        int d = both.GetLength(0);
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
            {
                var expected = i == j ? profileOnly[i, i] + fieldOnly[i, i] : profileOnly[i, j];
                Assert.True(both[i, j] == expected, $"[{i},{j}]: {both[i, j]} vs {expected}");
            }
        // At REAL q the split is real/imaginary: the profile owns the real part, the field the imaginary.
        // (Off the real axis the field's −i·q·(fe−fe) term acquires a real part and the split closes; q is real
        // above for exactly that reason.)
        for (int i = 0; i < d; i++)
        {
            Assert.True(both[i, i].Real == profileOnly[i, i].Real, $"the field moved the rate at {i}");
            Assert.True(both[i, i].Imaginary == fieldOnly[i, i].Imaginary, $"the profile moved the frequency at {i}");
        }
    }

    [Fact]
    public void UnderAProfileTheDiagonalStopsBeingAFunctionOfNDiffAlone()
    {
        // The physics the parameter exists for, stated as a measurement rather than as prose. At uniform γ every
        // coherence with the same n_diff shares one DIAGONAL entry, which at Δ = 0 with no field is what makes
        // a block with constant n_diff normal (the Edge lemma; the Δ·ZZ frequency varies across such a block)
        // and a whole sector sit at a single Re = −2γ (the class-B family). Under
        // a profile that shared value spreads inside [−2γ_max·n_diff, −2γ_min·n_diff]. This is measured on the
        // diagonal, which is what the builder owns; the class doc points at who owns the spectral consequence.
        const int n = 5;
        var gamma = DyadicProfile(n);
        var uniformSums = new HashSet<double>();
        var profileSums = new HashSet<double>();
        var kets = WeightCoherenceBlock.Configs(n, 1);
        var bras = WeightCoherenceBlock.Configs(n, 2);
        foreach (int kc in kets)
            foreach (int bc in bras)
            {
                if (System.Numerics.BitOperations.PopCount((uint)(kc ^ bc)) != 1) continue;   // fixed n_diff
                uniformSums.Add(WeightCoherenceBlock.DisagreementSum(n, null, kc, bc));
                profileSums.Add(WeightCoherenceBlock.DisagreementSum(n, gamma, kc, bc));
            }
        Assert.Single(uniformSums);                               // one shared value: the collapse
        Assert.True(profileSums.Count > 1,
            $"the profile left {profileSums.Count} distinct value(s) at n_diff = 1; the spread is the point");
    }

    [Fact]
    public void ANullProfileIsTheUniformCaseAndAWrongLengthOneThrows()
    {
        Assert.Equal(2.0, WeightCoherenceBlock.DisagreementSum(4, null, 0b0011, 0b0000));   // popcount route
        Assert.Throws<ArgumentException>(() => WeightCoherenceBlock.DisagreementSum(4, new[] { 1.0, 1.0 }, 1, 0));
        Assert.Throws<ArgumentException>(() =>
            WeightCoherenceBlock.Build(4, 1, 2, Complex.One, 0.0, null, new[] { 1.0, 1.0 }));
    }
}
