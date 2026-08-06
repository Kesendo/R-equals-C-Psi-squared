using System;
using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Gates for the live block-spectrum witness (<c>inspect --root blockspectrum</c>),
/// which surfaces the N=9 per-joint-popcount-sector Liouvillian spectra the SLOW_N9 test banks.
/// Gate-first: the spectrum reconstruction (Gate 2) and the Absorption floor (Gate 3) are real
/// physics checks that fire if the sector-flat-index wiring or the Hamiltonian is wrong; Gate 4
/// pins the witness's banked-headline read to the committed chain_N9.json so it cannot drift.
///
/// <para>The site-resolved half (Gates 3b..3f) covers the (0,1) block under a per-site γ PROFILE:
/// the generator identity entry-wise on four graph families (3b) and far from J=1 (3b-large-J), the
/// basis-order permutation with its mutation control (3b'), Herm(M) = −2·diag(γ) directly and
/// BIT-EXACTLY (3c, 3c-large-J), the per-mode Absorption law (3d) with the row that earns its scaled
/// threshold (3d-scale), the Bendixson and trace corollaries (3e), and the N=2 case where the floor
/// line does NOT open (3f). Two things are deliberate. The graph rows: on a uniform-J chain or ring
/// a global site-versus-bit reversal cancels, so the per-bond-J and star rows are what make the
/// ordering falsifiable. And which residual is EXACT and which is not, which is a statement about
/// the construction rather than about the identity. Herm(M) = −2·diag(γ) is bit-exact always and is
/// asserted as an exact 0.0. The per-mode residual has no exact route at all (an eigensolver on a
/// non-normal M) and is bounded by eps·‖M‖. The generator residual is the interesting one and it
/// has its own note at Gate 3b-large-J.</para></summary>
public class BlockSpectrumWitnessTests
{
    [Fact]
    public void Children_Provenance_LiveReconstructionVsBankedHeadline()
    {
        var w = new BlockSpectrumWitness(n: 5);   // small N -> live reconstruction is cheap
        var kids = ((IInspectable)w).Children.ToList();
        IInspectable Find(string namePart) => kids.Single(c => c.DisplayName.Contains(namePart));

        // Live-computed children (the badge must not contradict their "reconstructed live" prose).
        Assert.Equal(NodeProvenance.Live, Find("decomposition").Provenance);
        Assert.Equal(NodeProvenance.Live, Find("palindrome").Provenance);
        Assert.Equal(NodeProvenance.Live, Find("Absorption floor").Provenance);
        Assert.Equal(NodeProvenance.Live, Find("site-resolved band edge").Provenance);

        // Stored children: the static sector-map prose, and the banked N=9 headline (chain_N9.json).
        Assert.Equal(NodeProvenance.Stored, Find("sector map").Provenance);
        Assert.Equal(NodeProvenance.Stored, Find("banked").Provenance);
    }

    // Gate 1: the joint-popcount decomposition counts (the 100 -> 50 -> 25 story at N=9), pure
    // combinatorial. X(x)N order-2 pairing = the JSON PrimarySectorCount; the F1 Pi order-4 orbit
    // = the eig-calls the compute path actually does (Pi^2 = X(x)N).
    [Theory]
    [InlineData(8, 81, 41, 21, 4900L, 4, 4)]
    [InlineData(9, 100, 50, 25, 15876L, 4, 4)]
    public void Decomposition_CountsMatchTheKnownSectorStructure(
        int n, int sectorCount, int xnClasses, int piClasses, long maxBlock, int maxPc, int maxPr)
    {
        var d = BlockSpectrumWitness.Decomposition(n);
        Assert.Equal(sectorCount, d.SectorCount);
        Assert.Equal(xnClasses, d.XnClasses);
        Assert.Equal(piClasses, d.PiOrbitClasses);
        Assert.Equal(maxBlock, d.MaxBlock);
        Assert.Equal((maxPc, maxPr), (d.MaxPc, d.MaxPr));
    }

    // Gate 1b: the live cubic-cost speedup reproduces the banked N=9 EffectiveSpeedupOverDense
    // (645.95x) -- computed the identical (4^N)^3 / sum(block^3) way, so it must match.
    [Fact]
    public void Decomposition_CubicSpeedup_ReproducesTheBankedN9Value()
    {
        double speedup = BlockSpectrumWitness.Decomposition(9).CubicSpeedup;
        Assert.True(Math.Abs(speedup - 645.9495725858463) < 0.01,
            $"live cubic speedup {speedup} should reproduce the banked 645.95");
    }

    // Gate 2: the witness reconstructs the full spectrum sector-by-sector (every block via
    // BuildBlockZ) and it obeys F1: the multiset is symmetric about the center -sigma. A real
    // gate -- a wrong sector-flat-index extraction or a wrong H gives a spectrum that fails this.
    [Fact]
    public void ReconstructSpectrum_AtN5_IsFull_ObeysTheF1Palindrome_AboutMinusSigma()
    {
        const int n = 5;
        const double gamma = 0.5, j = 1.0;
        var (spectrum, _, skipped, full) = BlockSpectrumWitness.ReconstructSpectrum(n, gamma, j, cap: 2048);
        Assert.True(full);
        Assert.Equal(0, skipped);
        Assert.Equal(1 << (2 * n), spectrum.Length);   // 4^5 = 1024 (sector blocks sum to 4^N)

        double sigma = n * gamma;                       // 2.5
        // F1: {lambda} = {-2sigma - lambda}; the symmetric multiset has ~0 set-distance to its reflection.
        double pairing = BlockSpectrumWitness.PalindromePairingDistance(spectrum, sigma);
        Assert.True(pairing < 1e-9, $"F1 symmetry distance {pairing:E3} should be ~0");
        // the F1 floor: the maximally-dephased coherence (all N bits disagree) sits at Re = -2sigma.
        Assert.Equal(-2.0 * sigma, BlockSpectrumWitness.MinReal(spectrum), 9);
        // F4 kernel: a connected chain has kernel dimension N+1.
        Assert.Equal(n + 1, BlockSpectrumWitness.KernelDimension(spectrum));
    }

    // Gate 2b: the F1 metric must be MULTISET-aware, not set-aware. {+1 x3, -1 x1} is closed under
    // lambda -> -lambda as a SET (sigma=0) but NOT as a multiset (reflected = {-1 x3, +1 x1}). A
    // set-distance metric blindly reports ~0; the real F1 check (greedy NN with removal) must see the
    // gap -- the surplus +1 can only pair to a -1, distance 2. This is the failure mode a reconstruction
    // bug (a dropped/duplicated eigenvalue with a same-valued neighbour) would produce.
    [Fact]
    public void PalindromePairingDistance_DetectsAMultiplicityDefect()
    {
        var defective = new[]
        {
            new System.Numerics.Complex(1, 0), new System.Numerics.Complex(1, 0),
            new System.Numerics.Complex(1, 0), new System.Numerics.Complex(-1, 0),
        };
        double d = BlockSpectrumWitness.PalindromePairingDistance(defective, sigma: 0.0);
        Assert.True(d > 1.0, $"a multiplicity defect must be visible; a set-blind metric returns ~0, got {d:E3}");
    }

    // Gate 2c: the physics-review's STRONGEST counterexample, on a REAL spectrum (not a toy multiset).
    // Take a genuine F1-symmetric N=5 reconstruction, drop one eigenvalue and duplicate a different
    // (non-axis) one -- the exact defect a reconstruction-wiring bug produces. The old set/Hausdorff
    // metric gave ~3e-14 (invisible); the multiplicity-aware metric must flag it large.
    [Fact]
    public void PalindromePairingDistance_FlagsADropAndDuplicateInARealSpectrum()
    {
        const int n = 5;
        const double gamma = 0.5, j = 1.0;
        double sigma = n * gamma;
        var (clean, _, _, _) = BlockSpectrumWitness.ReconstructSpectrum(n, gamma, j, 2048);
        Assert.True(BlockSpectrumWitness.PalindromePairingDistance(clean, sigma) < 1e-9, "clean spectrum: F1 holds");

        // duplicate a non-axis eigenvalue (Re away from -sigma) over a different-valued slot:
        // removes one value, adds a copy of another -> breaks multiset F1 closure.
        int dup = System.Array.FindIndex(clean, z => System.Math.Abs(z.Real - (-sigma)) > 0.5);
        int drop = System.Array.FindIndex(clean, z => z != clean[dup]);
        Assert.True(dup >= 0 && drop >= 0 && clean[dup] != clean[drop], "test setup: distinct non-axis pair exists");
        var defective = (System.Numerics.Complex[])clean.Clone();
        defective[drop] = clean[dup];

        double d = BlockSpectrumWitness.PalindromePairingDistance(defective, sigma);
        Assert.True(d > 1e-6, $"a drop+duplicate defect must be visible; the old set-blind metric gave ~3e-14, got {d:E3}");
    }

    // Gate 3: the (0,1) band-edge sector is entirely at Re = -2gamma. On (p_c=0, p_r=1) every
    // basis coherence disagrees in exactly one bit, so L_D = -2gamma*I (scalar) and L_H restricted
    // is anti-Hermitian -> every eigenvalue has Re = -2gamma exactly (uniform gamma; Absorption).
    [Fact]
    public void BandEdgeSector_01_SitsEntirelyAtMinusTwoGamma()
    {
        const int n = 5;
        const double gamma = 0.5, j = 1.0;
        var (minRe, maxRe) = BlockSpectrumWitness.BandEdgeSectorReSpan(n, gamma, j);
        Assert.Equal(-2.0 * gamma, minRe, 9);
        Assert.Equal(-2.0 * gamma, maxRe, 9);
    }

    // The four graph families the site-resolved gates run on. Chain and ring at uniform J are
    // REVERSAL-INVARIANT, which is what lets a site<->bit mix-up hide; the per-bond J variants
    // below are not, and neither is the star.
    public static IEnumerable<object[]> BandEdgeGraphs()
    {
        foreach (int n in new[] { 3, 4, 5, 6, 7 })
        {
            yield return new object[] { n, "chain-uniform", Chain(n, _ => 1.0) };
            yield return new object[] { n, "chain-perbond", Chain(n, b => 0.4 + 0.9 * b) };   // non-palindromic
            yield return new object[] { n, "ring", Enumerable.Range(0, n).Select(i => new Bond(i, (i + 1) % n, 0.7)).ToArray() };
            yield return new object[] { n, "star", Enumerable.Range(1, n - 1).Select(i => new Bond(0, i, 1.3)).ToArray() };
        }
    }

    private static Bond[] Chain(int n, Func<int, double> j) =>
        Enumerable.Range(0, n - 1).Select(b => new Bond(b, b + 1, j(b))).ToArray();

    // The deliberately asymmetric profile: a reversed gamma cannot pass, and neither can a
    // site-ordered prediction against a block whose basis order is site n-1 down to site 0.
    private static double[] AsymmetricGamma(int n) =>
        Enumerable.Range(0, n).Select(l => 0.13 + 0.31 * l).ToArray();

    // Gate 3b: the SITE-RESOLVED (0,1) block on ANY XXX bond graph. The engine has always returned
    // it with a per-site gamma; the witness caller was uniform-only. The generator identity, in
    // this witness's convention H = sum_b (J_b/4)*(XX+YY+ZZ-I), is M = +i*(1/2)*WeightedLaplacian
    // - 2*diag(gamma), entry-wise. It fires on a wrong sign, a wrong J factor, a missing ZZ degree
    // term (which would leave the adjacency form, D10's recorded failure mode), a REVERSED gamma
    // profile (the failure PerBlockLiouvillianBuilder's own gamma-index-convention paragraph
    // records and PerBlockLiouvillianBuilderGammaOrderTests gates), and -- only because of the
    // per-bond and star rows -- a site<->bit mix-up in the Laplacian layout.
    //
    // The threshold is a SCALING LAW, not a constant: the residual is an absolute quantity built
    // from cancelling O(J*N) terms, so a fixed 1e-12 would silently become a real tolerance at
    // large J. Scale = J_max*N + 2*max(gamma). What carries that argument is Gate 3b-earns-the-scale
    // further down, which spans two decades on the NON-uniform profile; the four uniform large-J rows
    // immediately below are all exactly 0.0 and prove nothing about the scaling.
    [Theory]
    [MemberData(nameof(BandEdgeGraphs))]
    public void SiteResolvedBandEdgeBlock_IsExactlyTheLaplacianPlusDiagGammaGenerator(int n, string label, Bond[] bonds)
    {
        var gamma = AsymmetricGamma(n);
        double scale = bonds.Max(b => Math.Abs(b.Coupling)) * n + 2.0 * gamma.Max();
        double residual = BlockSpectrumWitness.GeneratorResidual(n, bonds, gamma);
        Assert.True(residual < 1e-13 * scale,
            $"[{label} N={n}] entry-wise generator residual {residual:E3} exceeds 1e-13*scale ({1e-13 * scale:E3})");
    }

    // Gate 3b-large-J: the same identity far off J=1. These four rows are all UNIFORM chains and all
    // measure exactly 0.0, so they do not exercise the scaled threshold at all -- a flat 1e-13 would
    // pass every one of them. The row that earns it is the non-uniform one below.
    //
    // THE RESIDUAL HERE IS NOT NOISE, IT IS A READING OF HOW H WAS ASSEMBLED. A site's diagonal
    // accumulates (+q, −q) for every bond it does NOT touch, and those pairs cancel around a running
    // sum. Whether they cancel bit for bit therefore depends on the ORDER the bonds arrive in, and
    // on nothing else about the physics: measured on one graph with one coupling set, changing only
    // the list order, ascending |J| leaves 1.7e-16 on the chain-end site and descending leaves 0.0.
    //
    // No ordering is exact in general: sorting the list descending by |J| makes the J~1 chain rows in
    // this file exact, but it makes the large-J one EIGHT TIMES WORSE (chain-perbond at J = pi*1e8,
    // N=5: 7.5e-09 as given, 5.96e-08 sorted descending), and on random graphs with random couplings
    // it lands nonzero about as often as ascending does. Descending is a property of particular rows,
    // not an ordering rule.
    // Note what the residual is NOT sensitive to: GeneratorResidual hands ONE list to both routes, so
    // a sort applied here moves H and the prediction together. Sorting inside the H builder alone is a
    // different measurement and gives different numbers; the two must not be quoted as one experiment.
    // So the identity is exact while the route to it is not, and
    // this gate's scaled threshold bounds the route, not the claim. Two constructions were removed
    // to get this far and both were real: the vacuum constant E0 = Σ J_b/4, which made every
    // diagonal a difference of two full sums (see HeisenbergGraph), and a prediction that halved a
    // total where the block quarters per bond (see PredictedGenerator). What is left is the order.
    [Theory]
    [InlineData(4, 1e5)]
    [InlineData(6, 1e8)]
    [InlineData(6, 1e-4)]
    [InlineData(6, 3.141592653589793e8)]
    public void SiteResolvedBandEdgeBlock_GeneratorIdentity_HoldsFarFromJEqualsOne(int n, double j)
    {
        var gamma = AsymmetricGamma(n);
        var bonds = Chain(n, _ => j);
        double scale = j * n + 2.0 * gamma.Max();
        double residual = BlockSpectrumWitness.GeneratorResidual(n, bonds, gamma);
        Assert.True(residual < 1e-13 * scale,
            $"[N={n} J={j}] residual {residual:E3} exceeds 1e-13*scale ({1e-13 * scale:E3})");
    }

    // Gate 3b-earns-the-scale: the rows a FLAT threshold would fail. The non-palindromic per-bond
    // profile at large J rounds -- at N=7, J_max ~ 1.5e9, the residual is 2.2e-08, five orders above a
    // flat 1e-13 -- while its ratio to J_max stays around 1e-18 to 2e-18. TWO decades are gated, not
    // one, because a single scale cannot tell eps*J_max*N apart from a constant or from J_max^2.
    //
    // What the upper bound is and is NOT. Over pi*10^k, k = 0..9, the ratio residual/(eps*J_max*N)
    // ranges from exactly 0.0 up to 0.14. That is an ENVELOPE, not a constant, so this is not a case-2
    // error model in the sense the rounding rule means; the residual is the file's case-3 object, a
    // reading. The bound below is therefore a ceiling on that envelope, and the second assertion is
    // the one carrying the argument: these rows must stay ABOVE a flat 1e-13, or the claim that the
    // threshold has to scale has no witness left in the file.
    [Theory]
    [InlineData(5, 3.141592653589793e7)]
    [InlineData(6, 3.141592653589793e7)]
    [InlineData(7, 3.141592653589793e7)]
    [InlineData(5, 3.141592653589793e8)]
    [InlineData(6, 3.141592653589793e8)]
    [InlineData(7, 3.141592653589793e8)]
    public void SiteResolvedBandEdgeBlock_GeneratorIdentity_ScalesWithTheCouplingNotWithAConstant(int n, double j0)
    {
        // These scales specifically. The residual is a rounding artefact, so it VANISHES at some scales
        // by luck -- pi*1e9 gives exactly 0.0 at N=5 while sitting at 2.4e-07 for N=6 and N=7. This gate
        // does not claim the residual is always nonzero; it holds two scales at which it demonstrably
        // is. A zero here means the luck changed, not that the law did.
        var gamma = AsymmetricGamma(n);
        var bonds = Chain(n, b => (0.4 + 0.9 * b) * j0);
        double jMax = bonds.Max(b => Math.Abs(b.Coupling));
        double residual = BlockSpectrumWitness.GeneratorResidual(n, bonds, gamma);

        double eps = Math.Pow(2, -52);
        Assert.True(residual < 1.0 * eps * jMax * n,
            $"[chain-perbond N={n} J_max={jMax:E3}] residual {residual:E3} exceeds eps*J_max*N " +
            $"({eps * jMax * n:E3}) -- the error model, not a tuned number, has broken");
        Assert.True(residual > 1e-13,
            $"[chain-perbond N={n}] residual {residual:E3} no longer exceeds a flat 1e-13, so this " +
            "row has stopped earning the scaled threshold and Gate 3b's scaling has no witness left");
    }

    // Gate 3b-exact: the uniform-chain generator residual is exactly 0.0, at every N and at the
    // binary-exact couplings. Three rounds of prose got the exactness map wrong in three different
    // ways while no gate held any of it -- 3b's one-sided 1e-13*scale bound sits three orders above
    // these values, so a change in the builder's accumulation order could move them to 1e-16 and
    // stale the comment again in silence. This is the falsifiable half of that map.
    //
    // The map itself, measured on the 20 Gate-3b rows: SEVENTEEN are exactly 0.0 -- every ring row,
    // every star row, chain-uniform at every N, and chain-perbond at N=3 and N=4. The three that are
    // not are chain-perbond at N=5, 6, 7 (5.6e-17, 5.6e-17, 1.7e-16). Those seventeen have an exact
    // route, so by the repository's rounding rule they are compared exactly and not tolerated: see
    // GeneratorResidual_OnTheRowsWithAnExactRoute_IsExactlyZero below.
    [Theory]
    [InlineData(3, 1.0)]
    [InlineData(4, 1.0)]
    [InlineData(5, 1.0)]
    [InlineData(6, 1.0)]
    [InlineData(7, 1.0)]
    [InlineData(6, 1e5)]
    [InlineData(6, 1e8)]
    public void GeneratorResidual_OnTheUniformChain_IsExactlyZero(int n, double j)
    {
        var gamma = AsymmetricGamma(n);
        double residual = BlockSpectrumWitness.GeneratorResidual(n, Chain(n, _ => j), gamma);
        Assert.True(residual == 0.0,
            $"[chain-uniform N={n} J={j}] expected a bit-exact 0.0, got {residual:E3} -- if this is a " +
            "deliberate builder change, the exactness map in the comments above needs remeasuring");
    }

    // The seventeen Gate-3b fixtures that currently measure exactly 0.0. This PINS those fixtures; it
    // does not claim a law, and no law is known here. Equal couplings are not sufficient: at equal J
    // the complete graph K_4 is nonzero for J = 0.4, 0.7, 1.3, -2.7 and exact for J = 1, 2, e, pi,
    // so the outcome depends on the coupling's mantissa and on the degree jointly. The mechanism is
    // that a site's H diagonal also receives (+q, -q) from bonds it does NOT touch while the
    // prediction receives only incident bonds, so exactness is whether those pairs cancel bit for bit
    // in a running sum. Chains and stars at these sizes do; K_4 and up need not.
    //
    // Why pin them anyway: Gate 3b bounds all twenty by a scaled threshold that sits three orders
    // above these values, so a builder change could move any of them off 0.0 in silence. A failure
    // here is a finding about the construction to be read, not a tolerance to widen. The three rows
    // deliberately absent (chain-perbond at N=5, 6, 7) are the ones that round today; they stay under 3b.
    public static IEnumerable<object[]> ExactRouteGraphs()
    {
        foreach (int n in new[] { 3, 4, 5, 6, 7 })
        {
            yield return new object[] { n, "chain-uniform", Chain(n, _ => 1.0) };
            yield return new object[] { n, "ring", Enumerable.Range(0, n).Select(i => new Bond(i, (i + 1) % n, 0.7)).ToArray() };
            yield return new object[] { n, "star", Enumerable.Range(1, n - 1).Select(i => new Bond(0, i, 1.3)).ToArray() };
        }
        foreach (int n in new[] { 3, 4 })
            yield return new object[] { n, "chain-perbond", Chain(n, b => 0.4 + 0.9 * b) };
    }

    [Theory]
    [MemberData(nameof(ExactRouteGraphs))]
    public void GeneratorResidual_OnTheRowsWithAnExactRoute_IsExactlyZero(int n, string label, Bond[] bonds)
    {
        var gamma = AsymmetricGamma(n);
        double residual = BlockSpectrumWitness.GeneratorResidual(n, bonds, gamma);
        Assert.True(residual == 0.0,
            $"[{label} N={n}] expected a bit-exact 0.0, got {residual:E3} -- the exact route this row " +
            "had is gone, which is a finding about the builder and not a threshold to widen");
    }

    // Gate 3b': the block's basis order is NOT site order. The sector's flat indices row*d + col
    // with row = 1<<b sort ascending, and site l sits at bit n-1-l, so the block runs site n-1 down
    // to site 0. Two halves, because the first alone is arithmetically forced by the implementation
    // and would pass even if the bit<->site semantics were inverted: (a) the permutation is
    // [n-1..0]; (b) a prediction laid out in SITE order must NOT match the block. Half (b) uses a
    // NON-PALINDROMIC per-bond J, so a global site<->bit reversal cannot cancel between the
    // Laplacian and the gamma diagonal the way it does on a uniform chain.
    [Theory]
    [InlineData(3)]
    [InlineData(5)]
    [InlineData(6)]
    public void BandEdgeSectorSiteOrder_RunsFromTheLastSiteDownToTheFirst(int n)
    {
        Assert.Equal(Enumerable.Range(0, n).Reverse().ToArray(), BlockSpectrumWitness.BandEdgeSectorSiteOrder(n));

        var gamma = AsymmetricGamma(n);
        var bonds = Chain(n, b => 0.4 + 0.9 * b);
        var m = BlockSpectrumWitness.BandEdgeSectorBlock(n, bonds, gamma);
        var lapBySite = new double[n, n];
        foreach (var b in bonds)
        {
            lapBySite[b.Site1, b.Site2] -= b.Coupling; lapBySite[b.Site2, b.Site1] -= b.Coupling;
            lapBySite[b.Site1, b.Site1] += b.Coupling; lapBySite[b.Site2, b.Site2] += b.Coupling;
        }
        double worstSiteOrdered = 0.0;
        for (int r = 0; r < n; r++)
            for (int c = 0; c < n; c++)
            {
                var predicted = System.Numerics.Complex.ImaginaryOne * 0.5 * lapBySite[r, c];
                if (r == c) predicted -= 2.0 * gamma[r];
                worstSiteOrdered = Math.Max(worstSiteOrdered, (m[r, c] - predicted).Magnitude);
            }
        Assert.True(worstSiteOrdered > 0.1,
            $"a site-ordered prediction must NOT match the block; got {worstSiteOrdered:E3} -- if this passes, " +
            "the basis-order permutation is a no-op here and Gate 3b is not testing it");
    }

    // Gate 3c: the sharp statement the profile case rests on, gated directly rather than through a
    // consequence with slack -- Herm(M) = -2*diag(gamma) on every graph. Bendixson and the trace are
    // corollaries of this; asserting only them leaves a theorem standing in for a gate.
    //
    // BIT-EXACT, and the assertion says so rather than allowing a tolerance. It is exact by
    // construction, not by luck. Setting the FORM: on this block every basis element has col = 0,
    // so an off-diagonal is exactly -i*H[row_r,row_c]. (That one entry fires is true in any
    // joint-popcount block, not just this one, so it is not what makes the result exact.) Three
    // legs carry the exactness: (a) H is real symmetric with its two mirror entries accumulated
    // from the same terms in the same order, so (m[r,c] + conj(m[c,r]))/2 cancels bit for bit;
    // (b) the bra-ket disagreement here has popcount 1, so the dephasing contribution is the SINGLE
    // term -2*gamma_site and not a sum that could round (at higher popcount the same diagonal is a
    // sum of several gammas, where a re-summation in a different order need not agree bit for bit --
    // no row here leaves the (0,1) block, so that is stated as the mechanism, not as a measurement
    // this file took); (c) the H diagonal adds exactly +0.0 to the real part.
    // (a)-(c) are why this holds on the star and per-bond rows too. Do not read that as the two
    // exactness maps being complements: the GENERATOR residual is nonzero on only three of the twenty
    // rows (chain-perbond at N=5, 6, 7), so on seventeen of them both maps say exact, for unrelated
    // reasons. This one is exact by (a)-(c) on every row, including the three. One reachable exception, pathological: a coupling large enough to overflow the H
    // diagonal sum makes the real part NaN, and NaN == 0.0 is false. Everything finite is clean. Measured 0.0 at J = 1e8 as well as at J=1,
    // so a SCALED threshold here would be a loosening with nothing behind it -- an earlier round put
    // one in, with a J=1e5 row whose stated justification ("a flat 1e-14 would fail on rounding
    // alone") was simply false. If this ever stops being exact, that is a real finding about the
    // builder, and the gate should report it rather than absorb it.
    [Theory]
    [MemberData(nameof(BandEdgeGraphs))]
    public void SiteResolvedBandEdge_HermitianPartIsExactlyMinusTwoDiagGamma(int n, string label, Bond[] bonds)
    {
        var gamma = AsymmetricGamma(n);
        double residual = BlockSpectrumWitness.HermitianPartResidual(n, bonds, gamma);
        Assert.True(residual == 0.0, $"[{label} N={n}] Herm(M) + 2*diag(gamma) is not bit-exact 0: {residual:E3}");
    }

    // Gate 3c-large-J: the same exactness at J = 1e5 and 1e8, where the entries are 13 orders above
    // the gammas they must cancel against. Not a tolerance row -- an exactness row.
    [Theory]
    [InlineData(4, 1e5)]
    [InlineData(6, 1e8)]
    [InlineData(6, 3.141592653589793e8)]   // OFF the exact grid, where the generator residual WAS 6e-8 before the gauge
    public void SiteResolvedBandEdge_HermitianPart_StaysExactAtLargeJ(int n, double j)
    {
        var gamma = AsymmetricGamma(n);
        double residual = BlockSpectrumWitness.HermitianPartResidual(n, Chain(n, _ => j), gamma);
        Assert.True(residual == 0.0, $"[N={n} J={j}] Herm(M) + 2*diag(gamma) is not bit-exact 0: {residual:E3}");
    }

    // Gate 3d: the MASTER statement, and it belongs to AbsorptionTheoremClaim, not to this witness:
    // -Re(lambda_k) = 2*sum_l gamma_l*<Delta_l>_k with <Delta_l>_k the mode's site occupancy, on the
    // RIGHT eigenvectors of a non-normal M.
    //
    // WHAT THIS CAN AND CANNOT CATCH, stated because it would otherwise read stronger than it is:
    // given Gate 3c, this is a linear-algebra identity (Re lambda is the Rayleigh quotient of
    // Herm(M) on an eigenvector), so it CANNOT fail independently -- the only things that fire it
    // are a broken eigensolver or the same gamma-index error 3c already catches. It is kept as the
    // cross-check that ties this block to its typed parent's per-channel law, and as the executable
    // form of "non-normality costs nothing here", not as independent evidence.
    //
    // This residual is inexact at every J and grows as eps*||M||. On ONE graph family so the series
    // is a series: uniform chain N=6, 1.3e-15 at J=1, 1.1e-10 at J=1e5, 3.4e-8 at J=1e8. (An earlier
    // round anchored this same series at 5.3e-15, which is the chain-perbond row -- a different
    // family, four times the J=1 baseline. Do not mix rows inside a scaling claim.) So the scaled
    // threshold belongs here, and the large-J row below is what earns it. The error model is plain
    // backward error: the residual compares an eigenvalue against the Rayleigh quotient of the SAME
    // computed eigenvector, so writing r = Mv - lambda*v it is bounded by ||r||/||v|| ~ eps*||M||.
    // The eigenvector conditioning does NOT enter -- an earlier round claimed it did and used that
    // to justify extra room; driving cond(V) to 1e8 at the N=2 defective point leaves the residual
    // at 1e-16.
    [Theory]
    [MemberData(nameof(BandEdgeGraphs))]
    public void SiteResolvedBandEdge_ObeysThePerModeAbsorptionLaw(int n, string label, Bond[] bonds)
    {
        var gamma = AsymmetricGamma(n);
        double scale = bonds.Max(b => Math.Abs(b.Coupling)) * n + 2.0 * gamma.Max();
        double residual = BlockSpectrumWitness.PerModeAbsorptionResidual(n, bonds, gamma);
        Assert.True(residual < 1e-13 * scale,
            $"[{label} N={n}] per-mode absorption residual {residual:E3} exceeds 1e-13*scale ({1e-13 * scale:E3})");
    }

    // Gate 3d-scale: the row that earns the scaling above, and the only place in this file where a
    // scaled threshold is doing real work. A flat 1e-13 fails at J=1e5 (1.1e-10) and at J=1e8
    // (3.4e-8); eps*||M|| tracks both.
    [Theory]
    [InlineData(6, 1e5)]
    [InlineData(6, 1e8)]
    public void PerModeAbsorptionResidual_GrowsLikeEpsTimesNormM(int n, double j)
    {
        var gamma = AsymmetricGamma(n);
        var bonds = Chain(n, _ => j);
        double scale = j * n + 2.0 * gamma.Max();
        double residual = BlockSpectrumWitness.PerModeAbsorptionResidual(n, bonds, gamma);
        Assert.True(residual < 1e-13 * scale,
            $"[N={n} J={j}] per-mode residual {residual:E3} exceeds 1e-13*scale ({1e-13 * scale:E3})");
        // and it really is inexact here, so the row is not silently testing an exact zero
        Assert.True(residual > 1e-12, $"[N={n} J={j}] expected a resolvable residual, got {residual:E3}");
    }

    // Gate 3e: what the profile does to the floor, as the two corollaries. Bendixson brackets every
    // Re in [-2*max gamma, -2*min gamma] and the trace pins the MEAN at -2*gamma_bar; the deep-edge
    // profile keeps sigma fixed, so the F1 center -sigma does not move either. N >= 3: at N=2 the
    // block is 2x2 and the line stays CLOSED for J >= 2*(gamma difference), so the "has opened"
    // half below is false there -- which is why the witness renders that sentence conditionally.
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void SiteResolvedBandEdge_BendixsonBracketAndTraceHold_WhileTheLineOpens(int n)
    {
        const double gamma = 0.5, j = 1.0;
        var g = BlockSpectrumWitness.DeepEdgeProfile(n, gamma);
        Assert.Equal(n * gamma, g.Sum(), 9);                       // sigma unchanged: the F1 center holds

        var (minRe, maxRe) = BlockSpectrumWitness.BandEdgeSectorReSpan(n, g, j);

        // Bendixson: the Hermitian part is exactly -2*diag(gamma), so no Re escapes its spectrum.
        Assert.True(minRe >= -2.0 * g.Max() - 1e-9, $"minRe {minRe} below the Bendixson floor {-2 * g.Max()}");
        Assert.True(maxRe <= -2.0 * g.Min() + 1e-9, $"maxRe {maxRe} above the Bendixson ceiling {-2 * g.Min()}");

        // the trace: sum of Re = -2*sigma exactly (the Absorption floor, in the mean).
        var m = BlockSpectrumWitness.BandEdgeSectorBlock(n, g, j);
        double traceRe = m.Evd().EigenValues.Sum(z => z.Real);
        Assert.Equal(-2.0 * n * gamma, traceRe, 9);

        // and the line really has opened -- otherwise the two assertions above are vacuous.
        Assert.True(maxRe - minRe > 1e-3, $"the profile must split the floor line; span width {maxRe - minRe:E3}");

        // the uniform case, same sigma, is the collapsed limit: span width 0.
        var (uMin, uMax) = BlockSpectrumWitness.BandEdgeSectorReSpan(n, gamma, j);
        Assert.Equal(-2.0 * gamma, uMin, 9);
        Assert.Equal(uMin, uMax, 9);
    }

    // Gate 3f: the N=2 exception the prose must not paper over. The 2x2 block splits by
    // 2*sqrt(d^2 - J^2/4) with d = |gamma_1 - gamma_0| the plain difference of the two rates, so once
    // |J| >= 2d the two eigenvalues have EQUAL Re and the line does not open at all. (|J|, not J:
    // the sign of the coupling never reaches Re, which is why the witness's own doc states the
    // criterion in the total coupling |c| -- BandEdgeSectorReSpan takes an H and not a J, and for a
    // multi-bond two-site graph c is the SUM of the bond couplings. This gate scans positive J.)
    // So "a profile opens the line" is
    // false in general, and the witness renders that sentence conditionally on the measured width.
    [Theory]
    [InlineData(1.0)]            // below 2d = 1.5: the line opens
    [InlineData(1.4)]
    [InlineData(1.5)]            // exactly at the coalescence
    [InlineData(2.5)]            // above: closed
    public void SiteResolvedBandEdge_AtNEqualsTwo_TheLineNeedNotOpen(double j)
    {
        const int n = 2;
        const double gamma = 0.5;
        var g = BlockSpectrumWitness.DeepEdgeProfile(n, gamma);
        Assert.Equal(n * gamma, g.Sum(), 9);
        double d = g[1] - g[0];                                    // 0.75 at gamma = 0.5
        var (minRe, maxRe) = BlockSpectrumWitness.BandEdgeSectorReSpan(n, g, j);
        Assert.Equal(-2.0 * n * gamma, minRe + maxRe, 9);          // the trace holds either way
        double expectedWidth = j < 2.0 * d ? 2.0 * Math.Sqrt(d * d - j * j / 4.0) : 0.0;
        // J = 2d is an exact coalescence (a defective 2x2), where the numerically resolved split
        // goes as sqrt(eps), not eps -- so the tolerance there is sqrt-scaled, not loosened by
        // taste. Away from it the identity is exact to 1e-9.
        double tol = Math.Abs(j - 2.0 * d) < 1e-12 ? 1e-6 : 1e-9;
        Assert.True(Math.Abs((maxRe - minRe) - expectedWidth) < tol,
            $"J={j}: width {maxRe - minRe:E6} vs predicted {expectedWidth:E6}");
    }

    // The deep-edge profile is the sibling witness's canal shape, scaled to mean gamma -- an
    // attribution, so it is gated: at gamma=1 it must reproduce SectorReductionWitness's profile
    // entry for entry, including the N=5 canal anchor [0.25, 1.5, 1.5, 1.5, 0.25].
    [Fact]
    public void DeepEdgeProfile_ReproducesTheSiblingCanalShapeAtGammaOne()
    {
        Assert.Equal(new[] { 0.25, 1.5, 1.5, 1.5, 0.25 }, BlockSpectrumWitness.DeepEdgeProfile(5, 1.0));
        for (int n = 3; n <= 8; n++)
        {
            double edge = 0.25, rest = (n - 2 * edge) / (n - 2);
            var expected = Enumerable.Range(0, n).Select(i => (i == 0 || i == n - 1) ? edge : rest).ToArray();
            Assert.Equal(expected, BlockSpectrumWitness.DeepEdgeProfile(n, 1.0));
        }
        // and the invariant the profile exists for: sum = N*gamma at every N >= 2.
        for (int n = 2; n <= 9; n++)
            Assert.Equal(n * 0.37, BlockSpectrumWitness.DeepEdgeProfile(n, 0.37).Sum(), 9);
        Assert.Throws<ArgumentOutOfRangeException>(() => BlockSpectrumWitness.DeepEdgeProfile(1, 0.5));
    }

    // Gate 4: the banked N=9 headline node reads the committed chain_N9.json and matches it. This
    // pins the witness's "stored" provenance to the artifact so it cannot silently drift.
    [Fact]
    public void ReadBankedN9_MatchesTheCommittedArtifact()
    {
        var b = BlockSpectrumWitness.ReadBankedN9();
        Assert.NotNull(b);
        Assert.Equal(262144, b!.SpectrumSize);
        Assert.Equal(-9.0, b.MinReal, 6);                 // = -2sigma, the F1 floor at N=9 (sigma=4.5)
        Assert.Equal(10, b.KernelDimension);              // = N+1, the F4 connected-chain kernel
        Assert.Equal(0.02727562511208863, b.DissipationGap, 9);
        Assert.Equal(5.736321706379341, b.MaxImag, 6);
        // The palindrome held to 3.48e-13 across the real span |MinReal| = 2*sigma = 9 (the radius
        // about the centre -sigma is 4.5), i.e. ~174 eps*2sigma. That is a BACKWARD ERROR from a
        // non-normal eigensolver: not zero, and not bit-exact. The physics behind the artifact has
        // no exact route, but THIS assertion does: the value is deserialized from a committed JSON
        // and cannot move at all unless the artifact is regenerated. So it is a PROVENANCE pin, not
        // a quality gate: it fires on any change to the artifact, better or worse, which is what is
        // wanted from a stored number. The neighbours above gate on 6 to 9 digits because their
        // values carry eigensolver noise in the last places (MinReal is -9.000000000000062); this
        // one is compared exactly because nothing recomputes it.
        Assert.True(b.MaxPairingDistance == 3.476966116100645E-13,
            $"the banked pairing distance moved: {b.MaxPairingDistance:E16}");
        Assert.Equal(0, b.OutlierPairCount);
        Assert.Equal(100, b.SectorCount);
        Assert.Equal(50, b.PrimarySectorCount);           // X(x)N order-2 classes
        Assert.Equal(15876, b.MaxBlockSize);
        Assert.True(Math.Abs(b.EffectiveSpeedup - 645.9495725858463) < 0.01);
    }

    // Render smoke: the default (N=6) witness renders six nodes without crashing, and the live
    // reconstruction node at N=6 is full (max block C(6,3)^2 = 400 < the 2048 cap).
    [Fact]
    public void Witness_RendersSixNodes_WithoutCrash()
    {
        var w = new BlockSpectrumWitness();   // default N=6, gamma=0.5, J=1
        Assert.Contains("N=6", w.DisplayName);
        Assert.False(string.IsNullOrWhiteSpace(w.Summary));

        var children = w.Children.ToList();
        Assert.Equal(6, children.Count);
        Assert.All(children, c => Assert.False(string.IsNullOrWhiteSpace(c.DisplayName)));
        Assert.All(children, c => Assert.False(string.IsNullOrWhiteSpace(c.Summary)));
    }

    // The sector map turns blockspectrum into the navigation hub: each load-bearing sector points
    // to the sector-specific witness(es) that zoom it. The five roots must all be reachable from it.
    [Fact]
    public void SectorMap_PointsToTheSectorSpecificWitnesses()
    {
        var w = new BlockSpectrumWitness();
        var map = w.Children.Single(c => c.DisplayName.Contains("sector map"));
        var text = string.Join(" | ", map.Children.Select(c => $"{c.DisplayName} {c.Summary}"));
        foreach (var root in new[] { "reduction", "ceiling", "horizon", "survivor", "secondclock" })
            Assert.Contains(root, text);
    }

    // The reverse leg of the hub: each sector-specific witness points back to blockspectrum (the
    // per-sector overview), so the navigation is bidirectional.
    [Fact]
    public void TheFiveSectorWitnesses_PointBackToBlockSpectrum()
    {
        var summaries = new[]
        {
            new SectorReductionWitness(n: 5).Summary,
            new StructuralCeilingWitness().Summary,
            new CoherenceHorizonWitness().Summary,
            new IncompletenessSurvivorWitness(6, 1.5).Summary,
            new SecondClockRegimeWitness().Summary,
        };
        Assert.All(summaries, s => Assert.Contains("blockspectrum", s));
    }
}
