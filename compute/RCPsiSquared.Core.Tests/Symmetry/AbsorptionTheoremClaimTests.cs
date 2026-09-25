using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using Xunit.Abstractions;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class AbsorptionTheoremClaimTests
{
    private readonly ITestOutputHelper _out;

    public AbsorptionTheoremClaimTests(ITestOutputHelper output) => _out = output;

    private static AbsorptionTheoremClaim BuildClaim() =>
        new AbsorptionTheoremClaim(new Pi2DyadicLadderClaim());

    private static string Inv(FormattableString text) => FormattableString.Invariant(text);

    // ---- the per-site law on a real Liouvillian (Theorem 2) ----
    //
    // The object: L = −i[H, ·] + Σ_l γ_l (Z_l·Z_l − ·) from PauliDephasingDissipator.BuildZ, with a
    // random COMPLEX Hermitian H (no reality, no symmetry, no number conservation) and a γ profile
    // that is neither uniform nor palindromic, so neither the uniform formula nor a reversed site
    // order can pass by coincidence. Two routes, two kinds of gate (CLAUDE.md, "No rounding"):
    //
    //   1. Herm(L) = diag(−2·Σ_l γ_l·Δ_l) is an EXACT route. Off the diagonal every entry of
    //      (L + L†)/2 is one entry plus its conjugated mirror, and they cancel bit for bit, because
    //      (A + A†)/2 is exactly Hermitian in floating point and −i·(a+bi) is computed exactly. On
    //      the diagonal both sides add over the same sites in the SAME ascending order l = 0..N−1
    //      (the builder adds −2γ_l per disagreeing site, PerSiteEigenmodeDecayRate adds γ_l and
    //      doubles at the end), and doubling commutes with rounding, so every partial sum agrees. So the residual is
    //      compared with 0.0. The N = 3 profile (0.1, 0.2, 0.3) is chosen so that the order MATTERS:
    //      (0.1 + 0.2) + 0.3 = 0.6000000000000001 but (0.3 + 0.2) + 0.1 = 0.6, and the test shows
    //      that a descending sum on the predicted side fires the gate. At N = 2 two terms commute and
    //      no order change is visible, which the test states as well.
    //
    //   2. The per-mode law needs an eigensolver, so there is no exact route. The error model is the
    //      backward error: the law compares Re λ̂ with the Rayleigh quotient of Herm(L) on the SAME
    //      computed vector v̂, so with r = Lv̂ − λ̂v̂ the deviation is at most ‖r‖/‖v̂‖, and a
    //      backward-stable eigensolver keeps that at a modest multiple of eps·‖L‖. The eigenvalue's
    //      condition number does NOT enter, which is why a non-normal L costs nothing here. The gate
    //      is a law, not a threshold: at every coupling the ratio deviation/(eps·‖L‖_F) lies in one
    //      fixed band [1/8, 16], for both N and both sides, and within one series (one N, one side)
    //      the four readings stay within a factor 16 of each other, so the deviation follows
    //      eps·‖L‖_F to within that factor rather than drifting inside the band. The coupling runs over six decades
    //      (10⁻² to 10⁴); ‖L‖_F runs over 4.8 of them at N = 2 and 4.7 at N = 3, because at J = 10⁻² the dissipator sets
    //      it, and the test asserts that it spans more than four. The band's lower edge is what makes
    //      the deviation TRACK ‖L‖: a deviation that stayed flat while ‖L‖ grew would fall out of it.
    //      The upper edge is empirical for a dense Hessenberg-QR eigensolver at these sizes (measured
    //      0.51 to 8.1 on MKL over all sixteen readings, series spreads 2.0 to 6.0; a numpy/OpenBLAS port
    //      of the same H read 0.84 to 7.3 with one series spreading 8.7, a spread within a factor 2 of
    //      the edge, so a run off MKL sits closer to them); at J = 10⁻² the rounding of the light sums
    //      themselves, of order eps·2Σ_l γ_l, is part of it. Worst-case backward-error bounds carry a
    //      higher power of the dimension and are not what this states.
    //
    // What the second gate can and cannot catch, stated so it does not read stronger than it is:
    // given the first, it is a linear-algebra identity and cannot fail independently. What fires it
    // is a light reading that is not ⟨v|N_l|v⟩/‖v‖² (a wrong site order, the uniform formula), a
    // broken eigensolver, or an H that does not drop out. A missing normalization would NOT fire it,
    // since the eigensolver returns unit vectors; PerSiteLightProfile_IsNormalizedAndSumsToTheLightCount
    // is the test that catches that. It is the executable form of Step 3 on right AND left
    // eigenvectors, the two-sided reading of proof §2.

    private static readonly double Eps = Math.BitIncrement(1.0) - 1.0;   // 2^-52

    /// <summary>The band the eigen-route deviation stays in, in units of eps·‖L‖_F, and the largest
    /// max/min spread one series of readings may show across the couplings.</summary>
    private const double BandLower = 1.0 / 8.0, BandUpper = 16.0, SeriesSpread = 16.0;

    /// <summary>How far above the band's upper edge a wrong law must land: six orders of magnitude,
    /// the margin the proof states for the two mutations.</summary>
    private const double MutationMargin = 1e6;

    /// <summary>Six decades of coupling at fixed γ: from dissipation-dominated to H-dominated.</summary>
    private static readonly double[] Couplings = { 1e-2, 1e0, 1e2, 1e4 };

    /// <summary>Neither uniform nor palindromic (γ_0 ≠ γ_{N−1}); at N = 3 the sum of all three
    /// rates depends on the order of addition, so the exact gate can see a reordering.</summary>
    private static double[] GammaProfile(int n) => n switch
    {
        2 => new[] { 0.1, 0.2 },
        3 => new[] { 0.1, 0.2, 0.3 },
        _ => throw new ArgumentOutOfRangeException(nameof(n)),
    };

    private static double Gaussian(Random rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    /// <summary>(A + A†)/2 with complex Gaussian A: exactly Hermitian in floating point, generically
    /// complex off the diagonal, and with no conserved quantity.</summary>
    private static ComplexMatrix RandomComplexHermitian(int n, int seed)
    {
        int d = 1 << n;
        var rng = new Random(seed);
        var a = ComplexMatrix.Build.Dense(d, d);
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                a[i, j] = new Complex(Gaussian(rng), Gaussian(rng));
        return (a + a.ConjugateTranspose()) * 0.5;
    }

    private static ComplexMatrix BuildL(int n, int seed, double coupling, IReadOnlyList<double> gamma) =>
        PauliDephasingDissipator.BuildZ(RandomComplexHermitian(n, seed) * coupling, gamma);

    private static ComplexVector BasisVector(int dim, int index)
    {
        var v = ComplexVector.Build.Dense(dim);
        v[index] = Complex.One;
        return v;
    }

    private static double[] Reversed(double[] values)
    {
        var flipped = (double[])values.Clone();
        Array.Reverse(flipped);
        return flipped;
    }

    /// <summary>How the predicted diagonal is formed: through the claim (the route under test), with
    /// the site order reversed (a wrong site map), or with the right sites summed in descending order
    /// (a reordering of the same terms).</summary>
    private enum DiagonalRoute { Claim, ReversedSites, DescendingSum }

    /// <summary>max over all entries of |(L + L†)/2 − diag(predicted)|, the predicted diagonal read
    /// from each basis coherence's sharp light profile along the given route.</summary>
    private static double HermitianPartResidual(ComplexMatrix l, IReadOnlyList<double> gamma, DiagonalRoute route)
    {
        var claim = BuildClaim();
        int dim = l.RowCount;
        double worst = 0.0;
        for (int p = 0; p < dim; p++)
        {
            var light = AbsorptionTheoremClaim.PerSiteLightProfile(BasisVector(dim, p));
            double predictedDiagonal = route switch
            {
                DiagonalRoute.Claim => -claim.PerSiteEigenmodeDecayRate(gamma, light),
                DiagonalRoute.ReversedSites => -claim.PerSiteEigenmodeDecayRate(gamma, Reversed(light)),
                _ => -2.0 * DescendingSum(gamma, light),
            };
            for (int q = 0; q < dim; q++)
            {
                var herm = 0.5 * (l[p, q] + Complex.Conjugate(l[q, p]));
                if (p == q) herm -= predictedDiagonal;
                worst = Math.Max(worst, herm.Magnitude);
            }
        }
        return worst;
    }

    private static double DescendingSum(IReadOnlyList<double> gamma, double[] light)
    {
        double sum = 0.0;
        for (int l = gamma.Count - 1; l >= 0; l--) sum += gamma[l] * light[l];
        return sum;
    }

    private readonly record struct LawReading(double MaxRatio, double MaxDeviation, double Unit);

    /// <summary>Eigendecompose <paramref name="m"/> and read, for every eigenpair (μ_k, v_k), the
    /// deviation |−Re μ_k − predictedRate(⟨Δ⟩_{v_k})| in units of eps·‖m‖_F. Passing L reads the right
    /// eigenvectors; passing L† reads the left ones (L†w = μw ⟺ w†L = μ̄w†, and Re μ̄ = Re μ).</summary>
    private static LawReading ReadLaw(ComplexMatrix m, Func<double[], double> predictedRate)
    {
        var evd = m.Evd();
        double unit = Eps * m.FrobeniusNorm();
        double worst = 0.0;
        for (int k = 0; k < m.RowCount; k++)
        {
            var light = AbsorptionTheoremClaim.PerSiteLightProfile(evd.EigenVectors.Column(k));
            double deviation = Math.Abs(-evd.EigenValues[k].Real - predictedRate(light));
            worst = Math.Max(worst, deviation);
        }
        return new LawReading(worst / unit, worst, unit);
    }

    [Theory]
    [InlineData(2, 20260925)]
    [InlineData(3, 20260926)]
    public void HermitianPart_IsTheSiteWeightedDisagreementDiagonal_ExactlyForAComplexH(int n, int seed)
    {
        var gamma = GammaProfile(n);
        foreach (double coupling in Couplings)
        {
            var l = BuildL(n, seed, coupling, gamma);
            double residual = HermitianPartResidual(l, gamma, DiagonalRoute.Claim);
            Assert.True(residual == 0.0,
                Inv($"[N={n} J={coupling:G3}] Herm(L) − diag(−2·Σ_l γ_l·Δ_l) is not bit-exact 0: {residual:E3}"));

            // The mutation that must fail: the same diagonal read with the site order reversed.
            double reversed = HermitianPartResidual(l, gamma, DiagonalRoute.ReversedSites);
            Assert.True(reversed > 0.0,
                Inv($"[N={n} J={coupling:G3}] the reversed site order passed the exact gate"));

            // The gate sees a reordering of the same terms where the input allows one: at N = 3 the
            // descending sum differs in the last bit, at N = 2 two terms commute and nothing can show.
            double descending = HermitianPartResidual(l, gamma, DiagonalRoute.DescendingSum);
            if (n == 3)
                Assert.True(descending > 0.0,
                    Inv($"[N=3 J={coupling:G3}] a descending sum passed the exact gate; the profile cannot see an order change"));
            else
                Assert.True(descending == 0.0,
                    Inv($"[N=2 J={coupling:G3}] two terms should commute exactly, got {descending:E3}"));
        }
    }

    [Theory]
    [InlineData(2, 20260925)]
    [InlineData(3, 20260926)]
    public void PerSiteLaw_HoldsOnEveryRightAndLeftEigenvector_InOneBandAroundTheRoundingScale(int n, int seed)
    {
        var claim = BuildClaim();
        var gamma = GammaProfile(n);
        Func<double[], double> law = light => claim.PerSiteEigenmodeDecayRate(gamma, light);

        double normAtBottom = 0.0, normAtTop = 0.0;
        var series = new Dictionary<string, List<double>> { ["right"] = new(), ["left"] = new() };
        foreach (double coupling in Couplings)
        {
            var l = BuildL(n, seed, coupling, gamma);
            double norm = l.FrobeniusNorm();
            if (coupling == Couplings[0]) normAtBottom = norm;
            if (coupling == Couplings[^1]) normAtTop = norm;
            foreach (var (side, m) in new[] { ("right", l), ("left", l.ConjugateTranspose()) })
            {
                var reading = ReadLaw(m, law);
                series[side].Add(reading.MaxRatio);
                _out.WriteLine(Inv($"N={n} J={coupling:G3} {side,-5}: ‖L‖_F = {norm:G4}, max|dev| = {reading.MaxDeviation:E3}, ") +
                               Inv($"ratio = {reading.MaxRatio:F2} (band [{BandLower}, {BandUpper}])"));
                Assert.True(reading.MaxRatio >= BandLower && reading.MaxRatio <= BandUpper,
                    Inv($"[N={n} J={coupling:G3} {side}] per-site law deviation {reading.MaxDeviation:E3} is ") +
                    Inv($"{reading.MaxRatio:F2}·eps·‖L‖_F, outside the band [{BandLower}, {BandUpper}]"));
            }
        }

        // One series follows eps·‖L‖_F: its four readings stay within a factor SeriesSpread.
        foreach (var (side, ratios) in series)
        {
            double spread = ratios.Max() / ratios.Min();
            _out.WriteLine(Inv($"N={n} {side,-5}: max/min over the couplings = {spread:F2} (at most {SeriesSpread})"));
            Assert.True(spread <= SeriesSpread,
                Inv($"[N={n} {side}] the readings spread by {spread:F2} across the couplings, more than {SeriesSpread}"));
        }

        // The band is held over a real range of ‖L‖, not only of J.
        double decades = Math.Log10(normAtTop / normAtBottom);
        _out.WriteLine(Inv($"N={n}: ‖L‖_F spans {decades:F2} decades while J spans six"));
        Assert.True(decades > 4.0,
            Inv($"[N={n}] ‖L‖_F spans only {decades:F2} decades; the band is not exercised across scales"));
    }

    [Theory]
    [InlineData(2, 20260925)]
    [InlineData(3, 20260926)]
    public void PerSiteLaw_Gate_RejectsTheUniformFormulaAndTheReversedSiteOrder(int n, int seed)
    {
        var claim = BuildClaim();
        var gamma = GammaProfile(n);
        double gammaBar = gamma.Average();
        Func<double[], double> uniformFormula = light => claim.EigenmodeDecayRate(light.Sum(), gammaBar);
        Func<double[], double> reversedSites = light => claim.PerSiteEigenmodeDecayRate(gamma, Reversed(light));
        double floor = MutationMargin * BandUpper;

        double smallest = double.PositiveInfinity;
        foreach (double coupling in Couplings)
        {
            var l = BuildL(n, seed, coupling, gamma);
            var uniform = ReadLaw(l, uniformFormula);
            var reversed = ReadLaw(l, reversedSites);
            smallest = Math.Min(smallest, Math.Min(uniform.MaxRatio, reversed.MaxRatio));
            _out.WriteLine(Inv($"N={n} J={coupling:G3}: uniform-formula ratio {uniform.MaxRatio:E2}, ") +
                           Inv($"reversed-site ratio {reversed.MaxRatio:E2} (must exceed {floor:E1})"));
            Assert.True(uniform.MaxRatio > floor,
                Inv($"[N={n} J={coupling:G3}] the uniform formula 2γ̄·⟨n_XY⟩ lands at {uniform.MaxRatio:E2}·eps·‖L‖_F, ") +
                Inv($"not six orders above the band"));
            Assert.True(reversed.MaxRatio > floor,
                Inv($"[N={n} J={coupling:G3}] the reversed site order lands at {reversed.MaxRatio:E2}·eps·‖L‖_F, ") +
                Inv($"not six orders above the band"));
        }
        _out.WriteLine(Inv($"N={n}: smallest mutation ratio {smallest:E2} = {smallest / BandUpper:E2} × the band's upper edge"));
    }

    [Theory]
    [InlineData(HamiltonianType.XY)]
    [InlineData(HamiltonianType.Heisenberg)]
    public void PerSiteLaw_PairCoherenceAtN2_IsAnExactEigenmodeAtBothSitesPrices(HamiltonianType hType)
    {
        // Proof §2, "the recurring four": |00⟩⟨11| is an exact N = 2 eigenmode at −2(γ_0 + γ_1).
        var claim = BuildClaim();
        var gamma = GammaProfile(2);
        var h = new ChainSystem(N: 2, J: 1.0, GammaZero: 0.0, HType: hType).BuildHamiltonian();
        var l = PauliDephasingDissipator.BuildZ(h, gamma);
        const int pairIndex = 0 * 4 + 3;   // |00⟩⟨11|, row-major vec
        double rate = claim.PerSiteEigenmodeDecayRate(gamma, AbsorptionTheoremClaim.PerSiteLightProfile(BasisVector(16, pairIndex)));
        Assert.Equal(2.0 * (0.1 + 0.2), rate);

        var column = l.Column(pairIndex);
        double residual = 0.0;
        for (int p = 0; p < 16; p++)
        {
            var expected = p == pairIndex ? new Complex(-rate, 0.0) : Complex.Zero;
            residual = Math.Max(residual, (column[p] - expected).Magnitude);
        }
        Assert.True(residual == 0.0, Inv($"[{hType}] L·|00⟩⟨11| + rate·|00⟩⟨11| is not bit-exact 0: {residual:E3}"));
    }

    // ---- the per-site members, arithmetic ----

    [Fact]
    public void PerSiteLightProfile_OfABasisCoherence_IsItsDisagreementBits_SiteZeroLeftmost()
    {
        // N = 2, index x = a·4 + b for |a⟩⟨b|; site 0 is the leftmost Kronecker factor (the high bit).
        Assert.Equal(new[] { 0.0, 0.0 }, AbsorptionTheoremClaim.PerSiteLightProfile(BasisVector(16, 0 * 4 + 0)));
        Assert.Equal(new[] { 1.0, 1.0 }, AbsorptionTheoremClaim.PerSiteLightProfile(BasisVector(16, 0 * 4 + 3)));
        Assert.Equal(new[] { 0.0, 1.0 }, AbsorptionTheoremClaim.PerSiteLightProfile(BasisVector(16, 1 * 4 + 0)));
        Assert.Equal(new[] { 1.0, 0.0 }, AbsorptionTheoremClaim.PerSiteLightProfile(BasisVector(16, 2 * 4 + 0)));
        // N = 3: |000⟩⟨101| disagrees at sites 0 and 2.
        Assert.Equal(new[] { 1.0, 0.0, 1.0 }, AbsorptionTheoremClaim.PerSiteLightProfile(BasisVector(64, 0 * 8 + 5)));
    }

    [Fact]
    public void PerSiteLightProfile_IsNormalizedAndSumsToTheLightCount()
    {
        // Equal weight on |00⟩⟨00| (no light) and |00⟩⟨11| (light at both sites), unnormalized on
        // purpose (norm² = 18): this is the test a missing /‖v‖² fails, the eigen gate cannot.
        var v = ComplexVector.Build.Dense(16);
        v[0] = new Complex(0.0, 3.0);
        v[3] = new Complex(3.0, 0.0);
        var light = AbsorptionTheoremClaim.PerSiteLightProfile(v);
        Assert.Equal(new[] { 0.5, 0.5 }, light);
        Assert.Equal(1.0, light.Sum());
    }

    [Fact]
    public void PerSiteEigenmodeDecayRate_AtUniformGamma_ReducesToEigenmodeDecayRate()
    {
        // Dyadic inputs, so both routes are exact and agree bit for bit. For general inputs
        // γ·l₁ + γ·l₂ and γ·(l₁ + l₂) round differently: the reduction is an identity of exact
        // arithmetic, and this is an arithmetic check of it, not an exactness pin.
        var c = BuildClaim();
        var light = new[] { 0.5, 0.25, 1.0 };
        double perSite = c.PerSiteEigenmodeDecayRate(new[] { 0.25, 0.25, 0.25 }, light);
        Assert.Equal(0.875, perSite);
        Assert.Equal(c.EigenmodeDecayRate(light.Sum(), 0.25), perSite);
    }

    [Fact]
    public void PerSiteEigenmodeDecayRate_SharpProfile_IsTheSumOfTheDisagreeingSitesPrices()
    {
        var c = BuildClaim();
        Assert.Equal(2.0 * 0.5, c.PerSiteEigenmodeDecayRate(new[] { 0.25, 0.5, 0.125 }, new[] { 0.0, 1.0, 0.0 }));
        Assert.Equal(2.0 * (0.25 + 0.125), c.PerSiteEigenmodeDecayRate(new[] { 0.25, 0.5, 0.125 }, new[] { 1.0, 0.0, 1.0 }));
    }

    [Fact]
    public void PerSiteEigenmodeDecayRate_RejectsMalformedProfiles()
    {
        var c = BuildClaim();
        Assert.Throws<ArgumentException>(() => c.PerSiteEigenmodeDecayRate(new[] { 0.1, 0.2 }, new[] { 0.5 }));
        Assert.Throws<ArgumentException>(() => c.PerSiteEigenmodeDecayRate(Array.Empty<double>(), Array.Empty<double>()));
        Assert.Throws<ArgumentOutOfRangeException>(() => c.PerSiteEigenmodeDecayRate(new[] { -0.1 }, new[] { 0.5 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => c.PerSiteEigenmodeDecayRate(new[] { 0.1 }, new[] { 1.5 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => c.PerSiteEigenmodeDecayRate(new[] { 0.1 }, new[] { double.NaN }));
    }

    [Fact]
    public void PerSiteLightProfile_RejectsZeroAndNonPowerOfFourVectors()
    {
        Assert.Throws<ArgumentException>(() => AbsorptionTheoremClaim.PerSiteLightProfile(ComplexVector.Build.Dense(16)));
        Assert.Throws<ArgumentException>(() => AbsorptionTheoremClaim.PerSiteLightProfile(ComplexVector.Build.Dense(8, Complex.One)));
        Assert.Throws<ArgumentException>(() => AbsorptionTheoremClaim.PerSiteLightProfile(ComplexVector.Build.Dense(1, Complex.One)));
    }

    [Fact]
    public void PerSiteLawNode_IsLive_AndReadsBothSitesPricesForThePairCoherence()
    {
        var node = BuildClaim().Children.Single(ch => ch.DisplayName.StartsWith("the carrier is a vector"));
        Assert.Equal(NodeProvenance.Live, Assert.IsType<InspectableNode>(node).Provenance);
        Assert.Contains("profile (1, 1) and costs 0.3,", node.Summary);        // 2·(0.05 + 0.10)
        Assert.Contains("profile (0, 1) and cell cost 0.2,", node.Summary);    // 2·0.10, site 1 alone
    }

    [Fact]
    public void Surface_DoesNotCallTheInteractingSpectrumQuantized()
    {
        // What the guard forbids: calling the interacting spectrum quantized. 2γ₀ is a cell cost, not a
        // spacing of the interacting eigenvalues.
        var c = BuildClaim();
        foreach (var text in new[] { c.Name, c.DisplayName, c.Summary })
        {
            Assert.DoesNotContain("spectrum quantized", text, StringComparison.OrdinalIgnoreCase);
            Assert.DoesNotContain("rate-quantization", text, StringComparison.OrdinalIgnoreCase);
        }
    }

    // ---- the uniform surface ----

    [Fact]
    public void CorrelatedJumps_WithNoPauliStringJumpSet_KeepTheHermitianPartExactly()
    {
        // The condition for the integer Pauli reading is Herm(L_D) diagonal in the Pauli basis; jumps
        // that are single Pauli strings are one way to get it, not the only one. Jumps √a(Z₀ + iZ₁) and
        // √b(Z₀ − iZ₁) with a ≠ b have no Pauli-string jump set, and their dissipator differs from
        // (a+b)-rate Z-dephasing on sites 0 and 1 only by −i(a−b)(Z₀ρZ₁ − Z₁ρZ₀), which is
        // anti-Hermitian. √a = 1/2 and √b = 1/4 make every entry dyadic, so both comparisons are
        // exact (==), and the difference itself is 2|a−b| = 3/8 at its largest, not zero.
        const int d = 4;
        var z0 = ComplexMatrix.Build.DenseOfDiagonalArray(new Complex[] { 1, 1, -1, -1 });   // site 0 = high bit
        var z1 = ComplexMatrix.Build.DenseOfDiagonalArray(new Complex[] { 1, -1, 1, -1 });
        var id = ComplexMatrix.Build.DenseIdentity(d);
        var jumps = new[]
        {
            (z0 + Complex.ImaginaryOne * z1) * 0.5,
            (z0 - Complex.ImaginaryOne * z1) * 0.25,
        };
        var dissipator = ComplexMatrix.Build.Dense(d * d, d * d);
        foreach (var jump in jumps)
        {
            var jj = jump.ConjugateTranspose() * jump;
            dissipator += jump.KroneckerProduct(jump.Conjugate())
                - 0.5 * (jj.KroneckerProduct(id) + id.KroneckerProduct(jj.Transpose()));
        }
        const double rate = 0.25 + 0.0625;   // a + b = 5/16
        var zDephasing = PauliDephasingDissipator.BuildZ(ComplexMatrix.Build.Dense(d, d), new[] { rate, rate });

        var hermitianPart = 0.5 * (dissipator + dissipator.ConjugateTranspose());
        var rest = dissipator - zDephasing;
        double hermDifference = (hermitianPart - zDephasing).Enumerate().Max(c => c.Magnitude);
        double antiHermiticity = (rest + rest.ConjugateTranspose()).Enumerate().Max(c => c.Magnitude);
        double restSize = rest.Enumerate().Max(c => c.Magnitude);

        Assert.True(hermDifference == 0.0, Inv($"Herm(L_D) differs from Z-dephasing by {hermDifference:E3}"));
        Assert.True(antiHermiticity == 0.0, Inv($"the remainder is not anti-Hermitian: {antiHermiticity:E3}"));
        Assert.True(restSize == 0.375, Inv($"the remainder should be 2|a−b| = 0.375 at its largest, got {restSize:R}"));
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void DissipatorCoefficient_IsExactlyTwo()
    {
        // Exact route: Term(0) = 2^(1−0) is exactly 2.0, so compare with ==.
        Assert.True(BuildClaim().DissipatorCoefficient == 2.0);
    }

    [Fact]
    public void DissipatorCoefficientMatchesLiteral_HoldsExactly()
    {
        Assert.True(BuildClaim().DissipatorCoefficientMatchesLiteral());
    }

    // The uniform surface is arithmetic. Its inputs are dyadic, so every product and quotient below is
    // exact in floating point and the comparisons are exact.

    [Theory]
    [InlineData(0.0, 0.0)]
    [InlineData(0.0625, 0.125)]
    [InlineData(0.5, 1.0)]
    [InlineData(1.0, 2.0)]
    [InlineData(2.5, 5.0)]
    public void SingleDisagreementCellCost_IsTwoTimesGammaZero(double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().SingleDisagreementCellCost(gammaZero));
    }

    [Theory]
    [InlineData(0, 1.0, 0.0)]
    [InlineData(1, 1.0, 2.0)]
    [InlineData(3, 1.0, 6.0)]
    [InlineData(5, 1.0, 10.0)]
    [InlineData(3, 0.0625, 0.375)]
    public void EigenmodeDecayRate_IsTwoGammaTimesAverageNXY(int nXY, double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().EigenmodeDecayRate(nXY, gammaZero));
    }

    [Fact]
    public void EigenmodeDecayRate_AcceptsNonIntegerLightExpectation()
    {
        // 2·fl(4/3) = fl(8/3): doubling commutes with rounding.
        Assert.Equal(8.0 / 3.0, BuildClaim().EigenmodeDecayRate(averageNXy: 4.0 / 3.0, gammaZero: 1.0));
    }

    [Theory]
    [InlineData(0, 1.0, 0.0)]
    [InlineData(1, 1.0, 2.0)]
    [InlineData(3, 1.0, 6.0)]
    public void BasisPairDissipatorCost_EqualsEigenmodeFormulaAtSameNumericLight(int nDiff, double gammaZero, double expected)
    {
        var c = BuildClaim();
        Assert.Equal(expected, c.BasisPairDissipatorCost(nDiff, gammaZero));
        Assert.Equal(c.EigenmodeDecayRate(nDiff, gammaZero), c.BasisPairDissipatorCost(nDiff, gammaZero));
    }

    [Theory]
    [InlineData(1, 1.0, 2.0)]
    [InlineData(3, 1.0, 6.0)]
    [InlineData(5, 1.0, 10.0)]
    [InlineData(7, 0.0625, 0.875)]
    public void EigenmodeDecayRateCeiling_IsTwoGammaTimesN(int n, double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().EigenmodeDecayRateCeiling(n, gammaZero));
    }

    [Theory]
    [InlineData(0.0, 1.0, 0.0)]
    [InlineData(2.0, 1.0, 1.0)]
    [InlineData(6.0, 1.0, 3.0)]
    [InlineData(0.125, 0.0625, 1.0)]
    [InlineData(0.375, 0.0625, 3.0)]
    public void AverageNXyFromEigenmodeDecayRate_InvertsRate(double rate, double gammaZero, double expectedNXy)
    {
        Assert.Equal(expectedNXy, BuildClaim().AverageNXyFromEigenmodeDecayRate(rate, gammaZero));
    }

    [Theory]
    [InlineData(2, 1.0, 4.0)]
    [InlineData(3, 1.0, 6.0)]   // F89c path-2 anchor: pair-sum = 6γ for 3-qubit block
    [InlineData(4, 1.0, 8.0)]   // F89c path-3 prediction: pair-sum = 8γ for 4-qubit block
    [InlineData(3, 0.0625, 0.375)]
    public void HammingComplementCellCostSum_IsTwoGammaTimesBlockSize(int blockSize, double gammaZero, double expected)
    {
        Assert.Equal(expected, BuildClaim().HammingComplementCellCostSum(blockSize, gammaZero));
    }

    [Fact]
    public void HammingComplementCellCostSum_AtN3_MatchesF89cEmpiricalAnchor()
    {
        // F89c bit-exact verification: (SE,SE) + (SE,DE) eigenvalue pairs sum to 6γ at path-2.
        Assert.Equal(6.0, BuildClaim().HammingComplementCellCostSum(3, 1.0));
    }

    [Fact]
    public void EigenmodeDecayRate_NegativeAverageNXY_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().EigenmodeDecayRate(-1, 1.0));
    }

    [Fact]
    public void EigenmodeDecayRate_NegativeGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().EigenmodeDecayRate(1, -0.1));
    }

    [Fact]
    public void AverageNXyFromEigenmodeDecayRate_NegativeRate_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().AverageNXyFromEigenmodeDecayRate(-1.0, 1.0));
    }

    [Fact]
    public void AverageNXyFromEigenmodeDecayRate_ZeroGamma_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().AverageNXyFromEigenmodeDecayRate(1.0, 0.0));
    }

    [Fact]
    public void EigenmodeDecayRateCeiling_NLessThanOne_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().EigenmodeDecayRateCeiling(0, 1.0));
    }

    [Fact]
    public void HammingComplementCellCostSum_BlockSizeZero_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().HammingComplementCellCostSum(0, 1.0));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => new AbsorptionTheoremClaim(null!));
    }
}
