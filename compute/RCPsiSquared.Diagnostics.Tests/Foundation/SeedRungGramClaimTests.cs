using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>From-below gate for F143 (<see cref="SeedRungGramClaim"/>, <see cref="SeedRungGramWitness"/>).
///
/// <para>The division of labour between the assertions below is the point, and it follows the repo's
/// tolerance rule. Everything that has an exact route is compared exactly, in <c>long</c> arithmetic:
/// the closed form against an independently built 2·𝟏𝟏ᵀ + I + R, the three eigenrelations against 0, the
/// completeness of the certified eigenvectors against an exact GF(p) rank. Three statements have no exact
/// route, and each is asserted on its ERROR MODEL rather than on a number: the rebuild of W from the sines
/// and the same rebuild through an eigensolver are gated on residual/(ε·M) being bounded AND not growing
/// with N, and the disorder fence is gated on SEPARATIONS (against the gap it has to beat, against the
/// uniform chain's residual) rather than on absolute sizes. A threshold that merely passes at one N would
/// catch none of it.</para></summary>
public sealed class SeedRungGramClaimTests
{
    private static IEnumerable<int> ChainLengths()
    {
        for (int n = 2; n <= SeedRungGramWitness.SweepMax; n++) yield return n;
    }

    [Fact]
    public void ScaledGram_Is_Two_Ones_Plus_Identity_Plus_The_Reflection()
    {
        foreach (int n in ChainLengths())
        {
            var g = SeedRungGramClaim.ScaledGram(n);
            int m = SeedRungGramClaim.TransformLength(n);

            for (int a = 1; a <= n; a++)
                for (int c = 1; c <= n; c++)
                {
                    // built independently of the claim: 2·𝟏𝟏ᵀ, plus I, plus the permutation matrix of R.
                    long expected = 2L
                                    + (a == c ? 1L : 0L)
                                    + (SeedRungGramClaim.ChiralPartner(n, a) == c ? 1L : 0L);
                    Assert.True(g[a - 1, c - 1] == expected,
                        $"N={n}, entry ({a},{c}): {g[a - 1, c - 1]} against {expected} (M={m})");
                }
        }
    }

    [Fact]
    public void The_Three_Eigenrelations_Are_Exactly_Zero_At_Every_N()
    {
        foreach (int n in ChainLengths())
            Assert.True(SeedRungGramWitness.ExactResidual(n) == 0L,
                $"N={n}: integer eigenrelation residual {SeedRungGramWitness.ExactResidual(n)}, expected exactly 0");
    }

    [Fact]
    public void The_Certified_Eigenvectors_Are_Independent_And_Exhaust_The_Modes()
    {
        foreach (int n in ChainLengths())
        {
            Assert.Equal(n, SeedRungGramWitness.CertifiedDimension(n));
            Assert.Equal(n, SeedRungGramWitness.CertifiedVectors(n).Count);
            Assert.True(SeedRungGramWitness.CertifiedRank(n) == n,
                $"N={n}: GF(p) rank {SeedRungGramWitness.CertifiedRank(n)} of {n}; a dependent relation " +
                "would leave a multiplicity asserted rather than shown");
        }
    }

    [Fact]
    public void The_Kernel_Is_One_Chiral_Difference_Per_Pair_And_The_Middle_Mode_Is_Not_Among_Them()
    {
        foreach (int n in ChainLengths())
        {
            var pairs = SeedRungGramWitness.ChiralPairs(n);
            Assert.Equal(SeedRungGramClaim.FrozenNullity(n), pairs.Count);

            int middle = SeedRungGramClaim.TransformLength(n) / 2;
            foreach (var (mode, partner) in pairs)
            {
                Assert.Equal(partner, SeedRungGramClaim.ChiralPartner(n, mode));
                Assert.True(mode < partner, $"N={n}: pair ({mode},{partner}) listed in the wrong order");
                if (n % 2 == 1)
                    Assert.True(mode != middle && partner != middle,
                        $"N={n}: the self-paired mode {middle} must contribute no kernel vector");
            }
        }
    }

    [Fact]
    public void The_Frozen_Count_Is_The_Site_Side_Carriers_And_Not_A_Second_Spelling()
    {
        foreach (int n in ChainLengths())
            Assert.Equal(FrozenDivisorClaim.FrozenMultiplicity(n), SeedRungGramClaim.FrozenNullity(n));
    }

    [Theory]
    [InlineData(2, 1.0)]      // the middle band is empty here, so the gap runs straight to the simple 1
    [InlineData(3, 0.25)]
    [InlineData(7, 0.125)]
    [InlineData(9, 0.1)]
    public void The_Gap_Is_One_Over_The_Transform_Length_Once_The_Middle_Band_Is_Occupied(int n, double expected)
    {
        Assert.Equal(expected, SeedRungGramClaim.Gap(n));

        // and the N = 2 case is the empty middle band, not a special-cased N: the gap runs to the simple 1
        // exactly when there is nothing between the kernel and the top.
        Assert.Equal(n == 2, SeedRungGramClaim.MiddleBandMultiplicity(n) == 0);
    }

    [Fact]
    public void The_Gap_Closes_Linearly_Which_Is_What_A_Numeric_Rank_Read_Lives_On()
    {
        // The moat is 1/M: halving it takes DOUBLING N, not a step in N. Asserted as the identity itself
        // rather than as M·gap == 1.0, which is an exact float comparison that holds only by luck of range
        // (it survives every M in 3..41 and fails first at M = 49) and would break on raising SweepMax for
        // a reason with no physics in it.
        foreach (int n in ChainLengths())
        {
            if (SeedRungGramClaim.MiddleBandMultiplicity(n) == 0) continue;
            Assert.Equal(1.0 / SeedRungGramClaim.TransformLength(n), SeedRungGramClaim.Gap(n));
        }
    }

    [Fact]
    public void The_Rebuild_From_The_Sines_Obeys_Its_Error_Model_Rather_Than_A_Threshold()
    {
        var ratios = ChainLengths().ToDictionary(n => n, SeedRungGramWitness.ResidualErrorRatio);

        // (a) bounded: the residual is the accumulation of N products of O(1/M) entries, so residual/(ε·M)
        //     is O(1). A generous constant, because the assertion that matters is (b).
        foreach (var (n, ratio) in ratios)
            Assert.True(ratio < 10.0, $"N={n}: rebuild residual is {ratio:0.##}·ε·M, outside the error model");

        // (b) NOT GROWING: if the identity were only approximately true, the ratio would climb with N.
        //     Compared across a decade of N rather than at any single one.
        double small = ratios.Where(kv => kv.Key <= 10).Max(kv => kv.Value);
        double large = ratios.Where(kv => kv.Key >= 30).Max(kv => kv.Value);
        Assert.True(large < 3.0 * small,
            $"the ratio grows with N (≤10: {small:0.##}, ≥30: {large:0.##}), which is what a residual that " +
            "is not pure rounding does");
    }

    [Fact]
    public void The_Mode_Amplitude_Matrix_Is_Doubly_Stochastic()
    {
        // Isolates the sines from the squaring: this residual would move first if the modes were wrong.
        // A row sum is N terms of size O(1/M), so the error model is ε·N and NOT a constant — the first
        // version of this assertion used a fixed multiple of ε and failed at N=16 for exactly that reason,
        // which is the rule's own point: a threshold hides the question of why there is anything to
        // tolerate.
        const double eps = 2.220446049250313e-16;
        var ratios = ChainLengths().ToDictionary(n => n, n => SeedRungGramWitness.DoublyStochasticResidual(n) / (eps * n));

        foreach (var (n, ratio) in ratios)
            Assert.True(ratio < 3.0, $"N={n}: row sums deviate by {ratio:0.##}·ε·N, outside the error model");

        double small = ratios.Where(kv => kv.Key <= 10).Max(kv => kv.Value);
        double large = ratios.Where(kv => kv.Key >= 30).Max(kv => kv.Value);
        Assert.True(large < 3.0 * small,
            $"the ratio grows with N (≤10: {small:0.##}, ≥30: {large:0.##}); the modes would not be " +
            "orthonormal if it did");
    }

    /// <summary>The fence, and what it does NOT fence off. Two earlier versions of this test were wrong in
    /// opposite directions, and both are recorded because the shape recurs. The first asserted "the
    /// involution survives" by calling <c>ModeReflectionIsAnInvolution</c>, an index map no bond profile
    /// can touch: a tautology dressed as a control. The second overcorrected and claimed survival under an
    /// ARBITRARY profile, dropping the hypothesis the parent claim carries: ±ε pairing fixes the
    /// eigenVALUES, and reading v_{M−k} ∝ K v_k off it needs a SIMPLE spectrum, which an open chain has
    /// exactly when no bond is zero (a Jacobi matrix).
    /// <see cref="A_Zero_Bond_Breaks_The_Hypothesis_And_With_It_The_Survivals"/> is the control that would
    /// have caught it.
    ///
    /// <para>The residuals are gated against ε·N/gap and not against ε·N, because an eigenvector's
    /// backward error is O(ε‖h‖/gap): at bonds [1,1,1e-10,1,1], legal under the hypothesis, the chiral
    /// residual reaches 6.7e-6, which a constant times ε·N would have called a failure.</para></summary>
    [Fact]
    public void Bond_Disorder_Splits_Only_The_Middle_Band()
    {
        const double eps = 2.220446049250313e-16;
        var moats = new Dictionary<int, (double FirstNonKernel, double UniformGap)>();

        foreach (int n in new[] { 5, 7, 9, 12 })
        {
            var reading = SeedRungGramWitness.ReadHopping(n, SeedRungGramWitness.DefaultDisorder(n));
            double model = eps * n / reading.HoppingSpectralGap;

            // (a) what survives. W R = W measured, and G annihilating every chiral difference, which reads
            //     the kernel's IDENTITY rather than the size of some low eigenvalue.
            Assert.True(reading.ChiralSymmetryResidual < 100.0 * model,
                $"N={n}: W R = W broke under disorder, {reading.ChiralSymmetryResidual:0.###e+0} against " +
                $"a model of {model:0.###e+0}");
            Assert.True(reading.ROddAnnihilationResidual < 100.0 * model,
                $"N={n}: G stopped annihilating the R-odd sector, {reading.ROddAnnihilationResidual:0.###e+0}");

            // (b) the OTHER side of the count: an annihilation test alone certifies nullity ≥ ⌊N/2⌋ and
            //     never ≤. Without this the kernel could have grown and the test would still pass. Gated
            //     on the same error model, not on a size: what has to hold is that this eigenvalue is
            //     nonzero BEYOND rounding, and nothing about how large it is.
            Assert.True(reading.FirstNonKernelEigenvalue > 100.0 * model,
                $"N={n}: the first non-kernel eigenvalue fell to {reading.FirstNonKernelEigenvalue:0.###e+0}, " +
                $"within rounding of the model {model:0.###e+0}, so the reading no longer bounds the kernel " +
                "from above");

            moats[n] = (reading.FirstNonKernelEigenvalue, SeedRungGramClaim.Gap(n));

            // (d) what breaks, measured against the scale it has to beat: the gap it is supposed to be
            //     smaller than. A spread comparable to 1/M means the middle band is no longer a band.
            Assert.True(reading.MiddleBandSpread > 0.1 * SeedRungGramClaim.Gap(n),
                $"N={n}: middle-band spread {reading.MiddleBandSpread:0.####} against gap " +
                $"{SeedRungGramClaim.Gap(n):0.####}; if this were small the fence would test nothing");

            // (e) and with it the entry-wise closed form, compared to the uniform chain rather than to a
            //     bare number: the separation is the statement, the absolute value is not.
            Assert.True(reading.ClosedFormResidual > 1e6 * SeedRungGramWitness.ClosedFormResidual(n),
                $"N={n}: disordered {reading.ClosedFormResidual:0.###e+0} against uniform " +
                $"{SeedRungGramWitness.ClosedFormResidual(n):0.###e+0}");
        }

        // (f) the MOAT is not among the survivals, and the honest form of that is the loss of the
        //     GUARANTEE rather than of the value. Under this profile the first non-kernel eigenvalue falls
        //     below the uniform chain's 1/M at some N and stays above it at others (N=5 reads 0.171 against
        //     0.167), so what disorder removes is the bound, not the margin. Asserting "< gap" at every N
        //     would have been the overclaim, and it failed here before this comment existed.
        Assert.Contains(moats, kv => kv.Value.FirstNonKernel < kv.Value.UniformGap);
        Assert.Contains(moats, kv => kv.Value.FirstNonKernel > kv.Value.UniformGap);
    }

    /// <summary>The control for the hypothesis. A zero bond is a legal double array and an illegal chain:
    /// it disconnects, the h-spectrum degenerates, and the survivals above stop holding. Without this test
    /// the scope sentence would be prose nobody could falsify.
    ///
    /// <para>What is asserted, and why not the obvious thing. The first attempt asserted that the residuals
    /// come out large at N = 6 with bonds [1,1,0,1,1]. They do, through MathNet, and that assertion was
    /// still wrong: the chain splits into two IDENTICAL trimers, so every eigenvalue of h doubles, and the
    /// equally legal basis (u^A ± u^B)/√2 of the same eigenspaces gives W R = W exactly and residuals of
    /// 0.0. Pinning either number would be pinning a LAPACK convention. So the ill-posedness is
    /// DEMONSTRATED instead: the test builds the second basis, checks it really is an orthonormal
    /// eigenbasis of the same h, and shows the two disagree. That is a statement about h, not about a
    /// solver, and it is the actual content of "the hypothesis is violated".</para></summary>
    [Fact]
    public void A_Zero_Bond_Breaks_The_Hypothesis_And_The_Reading_Stops_Being_Well_Posed()
    {
        // N = 6 into two identical trimers, and N = 4 into two dimers: both degenerate.
        foreach (var (n, bonds) in new[] { (6, new[] { 1.0, 1.0, 0.0, 1.0, 1.0 }), (4, new[] { 1.0, 0.0, 1.0 }) })
        {
            var reading = SeedRungGramWitness.ReadHopping(n, bonds);

            // (a) the mechanism, and the one part no basis convention can move: h's spectrum degenerates.
            Assert.True(reading.HoppingSpectralGap < 1e-12,
                $"N={n}: the hypothesis is the simple h-spectrum; smallest spacing read " +
                $"{reading.HoppingSpectralGap:0.###e+0}");

            // (b) and therefore W is not a function of h. Shown by scanning the rotation angle inside ONE
            //     degenerate eigenspace at a time: every choice gives an orthonormal eigenbasis of the SAME
            //     h, so if the reading were well posed they would all agree. Two earlier versions of this
            //     assertion were too narrow and both failed here. A single 45° rotation reproduces the
            //     solver's own reading at N = 4; so does rotating EVERY eigenspace by the same angle, since
            //     that keeps the chiral partners in step. Only an independent rotation breaks it.
            var (modes, spectrum) = SeedRungGramWitness.DiagonaliseChain(n, bonds);
            var readings = new List<double>();
            for (int pair = 0; pair < n / 2; pair++)
                for (int step = 0; step <= 8; step++)
                {
                    var rotated = RotateOneDegenerateSpace(n, modes, spectrum, pair, Math.PI * step / 16.0);
                    readings.Add(SeedRungGramWitness.ChiralSymmetryOfModes(
                        n, SeedRungGramWitness.SquaredAmplitudes(n, rotated)));
                }

            Assert.True(readings.Max() - readings.Min() > 0.1,
                $"N={n}: valid eigenbases of one h must disagree here, which is what makes the reading " +
                $"ill-posed; the scan spans {readings.Min():0.####} to {readings.Max():0.####}");
        }
    }

    /// <summary>Rotates by θ inside the <paramref name="target"/>-th 2-fold degenerate eigenspace and
    /// leaves the others alone. Still orthonormal, still an eigenbasis of the same h at every θ, which is
    /// exactly the point being made.</summary>
    private static double[,] RotateOneDegenerateSpace(
        int n, double[,] modes, double[] spectrum, int target, double theta)
    {
        var rotated = (double[,])modes.Clone();
        double cos = Math.Cos(theta), sin = Math.Sin(theta);
        int seen = 0;
        for (int k = 0; k + 1 < n; k++)
        {
            if (Math.Abs(spectrum[k + 1] - spectrum[k]) > 1e-10) continue;
            if (seen++ == target)
            {
                for (int l = 0; l < n; l++)
                {
                    double a = modes[l, k], b = modes[l, k + 1];
                    rotated[l, k] = cos * a + sin * b;
                    rotated[l, k + 1] = -sin * a + cos * b;
                }
            }
            k++;   // consume the pair
        }
        return rotated;
    }

    [Fact]
    public void A_Uniform_Chain_Reproduces_The_Closed_Form_Through_The_Eigensolver_Too()
    {
        // The fence above would also fire on a mistake in the eigensolver route itself, so the same route
        // is walked at zero disorder, where it must agree with the sine formula. Same object as the sine
        // rebuild, so the same law applies: a ratio to ε·M, bounded and not growing. The constant is looser
        // than the sine route's because an eigenvector carries the solver's backward error on top.
        const double eps = 2.220446049250313e-16;
        var ratios = new Dictionary<int, double>();

        foreach (int n in new[] { 5, 7, 9, 12, 20, 30, 40 })
        {
            var uniformBonds = Enumerable.Repeat(1.0, n - 1).ToArray();
            ratios[n] = SeedRungGramWitness.DisorderedResidual(n, uniformBonds) / (eps * (n + 1));
        }

        foreach (var (n, ratio) in ratios)
            Assert.True(ratio < 10.0,
                $"N={n}: the eigensolver route at uniform bonds reads {ratio:0.##}·ε·M");

        double small = ratios.Where(kv => kv.Key <= 12).Max(kv => kv.Value);
        double large = ratios.Where(kv => kv.Key >= 30).Max(kv => kv.Value);
        Assert.True(large < 3.0 * small,
            $"the ratio grows with N (≤12: {small:0.##}, ≥30: {large:0.##}), which rounding alone would not");
    }

    [Fact]
    public void The_Mode_Reflection_Is_Not_An_Involution_On_The_Heisenberg_Chain()
    {
        foreach (int n in new[] { 4, 5, 6, 7, 8, 9 })
        {
            Assert.True(SeedRungGramClaim.ModeReflectionIsAnInvolution(1, n, n + 1),
                $"N={n}: the open XY chain's modes 1..N must be closed under a ↦ (N+1) − a");
            Assert.False(SeedRungGramClaim.ModeReflectionIsAnInvolution(0, n - 1, n),
                $"N={n}: the Heisenberg chain's modes 0..N−1 must NOT be closed under k ↦ N − k, since " +
                "k = 0 is sent to N");
        }
    }

    [Fact]
    public void The_Unscaled_Entries_Are_The_Proofs_Four_Cases()
    {
        const int n = 9;                                   // odd, so the self-paired middle mode exists
        int m = SeedRungGramClaim.TransformLength(n);
        int middle = m / 2;

        Assert.Equal(1.0 / m, SeedRungGramClaim.GramEntry(n, 1, 2));                 // generic
        Assert.Equal(3.0 / (2 * m), SeedRungGramClaim.GramEntry(n, 1, 1));           // diagonal
        Assert.Equal(3.0 / (2 * m), SeedRungGramClaim.GramEntry(n, 1, m - 1));       // anti-diagonal
        Assert.Equal(2.0 / m, SeedRungGramClaim.GramEntry(n, middle, middle));       // both at once
    }

    /// <summary>The two numbers the F143 registry entry and the arc quote out of this witness. Without
    /// this test nothing in the repo holds them, and a doc number no gate touches is how the arc's own
    /// inventory found stale cells elsewhere. Bounded rather than pinned to the last digit, because the
    /// max of the ratio reproduces across summation orders and its exact digits do not.</summary>
    [Fact]
    public void The_Numbers_The_Docs_Quote_Are_The_Numbers_The_Witness_Produces()
    {
        double worstRatio = ChainLengths().Max(SeedRungGramWitness.ResidualErrorRatio);
        Assert.InRange(worstRatio, 2.0, 2.2);          // quoted as "reaching 2.09"

        int n = 9;
        double onG = SeedRungGramWitness.DisorderedResidual(n, SeedRungGramWitness.DefaultDisorder(n))
                     / (2.0 * SeedRungGramClaim.TransformLength(n));
        Assert.InRange(onG, 0.09, 0.11);               // quoted as "≈0.099 on G at N = 9"
    }

    [Fact]
    public void Guards_Reject_What_The_Statement_Does_Not_Cover()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => SeedRungGramClaim.Gap(1));
        Assert.Throws<ArgumentOutOfRangeException>(() => SeedRungGramClaim.ChiralPartner(5, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => SeedRungGramClaim.ChiralPartner(5, 6));
        Assert.Throws<ArgumentOutOfRangeException>(() => new SeedRungGramWitness(1));
        Assert.Throws<ArgumentException>(() => SeedRungGramWitness.DisorderedResidual(5, new[] { 1.0, 1.0 }));
    }

    [Fact]
    public void Claim_Is_Registered_With_Both_Typed_Parents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<SeedRungGramClaim>());

        var claim = registry.Get<SeedRungGramClaim>();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
        Assert.NotNull(claim.ChiralK);
        Assert.NotNull(claim.FrozenDivisor);

        var ancestors = registry.AncestorsOf<SeedRungGramClaim>().Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(ChiralKClaim), ancestors);
        Assert.Contains(typeof(FrozenDivisorClaim), ancestors);
    }

    [Fact]
    public void Statement_Carries_Both_Fences_And_The_Anchor_Names_The_Prior_Art()
    {
        var claim = KnowledgeRegistryFactory.BuildDefault().Get<SeedRungGramClaim>();

        Assert.Contains("UNIFORM hopping", claim.Name);
        Assert.Contains("Heisenberg", claim.Name);
        Assert.Contains("⌊N/2⌋", claim.Name);
        Assert.Contains("PROOF_R90_FROZEN_DIVISOR.md Lemma 5", claim.Anchor);
        Assert.Contains("inspect --root seedrung", claim.Anchor);
    }
}
