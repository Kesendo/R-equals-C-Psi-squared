using System;
using RCPsiSquared.Core.F89PathK;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>Scout for the sideways_spin_ladder arc's Sturm closing route: the even-N disc layers of
/// the (1,2) residual through the D-only certificate path (<see cref="FoldResultantCertificate"/>,
/// even-N guard lifted 2026-08-10 — the D-only path never touches the corner block).
///
/// <para>N=4 splits into a CONTROL and a derived twin: the R-even layers are sourced from the
/// committed closed form, disc(F₈) = const·q²⁴·(3q⁴+q²−1)²·P₂₀(q) with P₂₀(q) = P₁₀(q²)
/// (<c>F89Path3OcticGaloisClaim</c>, which scopes itself to the R-even sector: deg 52, v_q 24,
/// factorisation (24, 8, 20) ⇒ layers [20, 4]); the R-odd row is NOT separately committed — at even
/// N the R-odd residual is the exact complex CONJUGATE of the R-even one (measured; the F89 doc's
/// sector-swap conjugacy), and conjugation fixes every layer degree, so the identical reading is
/// asserted as the derived twin, not quoted as a second source.</para>
///
/// <para>N=6 is the first data for the certificate design (residual degree 32 = the F_32 oracle's
/// symmetric-sector degree; R-odd again the conjugate). Gate: the layer aggregation identity
/// deg D − v_q = Σ (m+1)·deg(layer_m), the same consistency the N=7 disc-multiplicity test pins,
/// plus layers.Length ≤ 2 (no multiplicity-3 layer, the β-exotic point read at even N).</para>
///
/// <para>Since 2026-08-10 the class also holds the PROMOTED certificate's gates
/// (<see cref="FoldResultantCertificate.CertifyDiscReImGcd"/>): the N=6 window-free
/// no-real-nonzero-coalescence statement, both parities, and the N=4 real-disc control that must
/// fail it closed. Same day, the N=8 layer scout (Category SLOW_EVEN_DISC8, ~9 min): the dim-112
/// bivariate Berkowitz is affordable, deg_q D = 6086, v_q = 3730, layers [376, 990], so the N=8
/// certificate object is the degree-376 simple layer.</para></summary>
public class EvenNDiscLayerScoutTests
{
    private readonly ITestOutputHelper _out;
    public EvenNDiscLayerScoutTests(ITestOutputHelper o) => _out = o;

    [Fact(DisplayName = "N=4: disc layers reproduce the octic factorization [20, 4] (R-even sourced, R-odd the conjugate twin)")]
    [Trait("Category", "EVEN_DISC_SCOUT")]
    public void N4_Control_KnownOcticFactorization()
    {
        foreach (bool rOdd in new[] { false, true })
        {
            var (degD, vD, layers) = FoldResultantCertificate.DiscLayersAtNthPrime(4, rOdd, 0);
            _out.WriteLine($"N=4 {(rOdd ? "R-odd (conjugate twin)" : "R-even (control)")}: deg_q D = {degD}, v_q = {vD}, layers [{string.Join(", ", layers)}]");
            Assert.Equal(52, degD);
            Assert.Equal(24, vD);
            Assert.Equal(new[] { 20, 4 }, layers);
        }
    }

    [Fact(DisplayName = "N=4 Re/Im control: the self-fold disc is REAL, Im D ≡ 0 mod p")]
    [Trait("Category", "EVEN_DISC_SCOUT")]
    public void N4_ReImControl_DiscIsReal()
    {
        var (p, _, degRe, vRe, degIm, _, imZero, _) =
            FoldResultantCertificate.DiscReImGcdAtNthPrime(4, rOdd: false, nth: 0);
        _out.WriteLine($"N=4 R-even: p={p}, deg Re = {degRe} (v_q {vRe}), deg Im = {degIm}, ImIsZero = {imZero}");
        Assert.True(imZero, "the N=4 disc must be real (self-fold antiunitary, spec B3)");
        Assert.Equal(52, degRe);
        Assert.Equal(24, vRe);
    }

    [Fact(DisplayName = "N=6 scout: gcd(Re D, Im D) mod p = 1 — the closing certificate's live question, pinned")]
    [Trait("Category", "SLOW_EVEN_DISC")]
    public void N6_Scout_ReImGcd()
    {
        var rows = new System.Collections.Generic.List<(int DegRe, int VRe, int DegIm, int VIm, int GcdDeg)>();
        foreach (bool rOdd in new[] { false, true })
        {
            var (p, _, degRe, vRe, degIm, vIm, imZero, gcdDeg) =
                FoldResultantCertificate.DiscReImGcdAtNthPrime(6, rOdd, 0);
            _out.WriteLine($"N=6 {(rOdd ? "R-odd " : "R-even")}: p={p}, deg Re = {degRe} (v_q {vRe}), "
                           + $"deg Im = {degIm} (v_q {vIm}), gcd deg = {gcdDeg}");
            Assert.False(imZero, "the N=6 disc must be genuinely complex (disc-reality re-gate, 13-order split)");
            // PINNED, observed 2026-08-10: gcd degree 0 at the first split prime. Mod-p only — the
            // window-free certificate additionally needs the degree- and valuation-preservation
            // halves (the method's scope note); a real-q coalescence appearing at N=6 would flip
            // this to a nonzero gcd degree and fail here first.
            Assert.Equal(0, gcdDeg);
            rows.Add((degRe, vRe, degIm, vIm, gcdDeg));
        }
        // the even-N conjugation makes R-odd the exact conjugate of R-even (Re shared, Im negated):
        // one measurement, one forced twin — gate the agreement rather than presenting two sources.
        Assert.Equal(rows[0], rows[1]);
    }

    [Fact(DisplayName = "N=4 certificate control: the real self-fold disc MUST fail the complex-disc certificate closed")]
    [Trait("Category", "EVEN_DISC_SCOUT")]
    public void N4_Certificate_FailsClosedOnRealDisc()
    {
        foreach (bool rOdd in new[] { false, true })
        {
            var r = FoldResultantCertificate.CertifyDiscReImGcd(4, rOdd);
            _out.WriteLine($"N=4 {(rOdd ? "R-odd " : "R-even")}: certified = {r.Certified}, deg D = {r.TrueDiscriminantDegree}, "
                           + $"v_q = {r.TrueQValuationD}, sampled = {r.PrimesSampled} (bound {r.LcDivisorBoundD})");
            // the D-side invariants certify fine (the closed octic form), but no prime can see a
            // complex disc: Im ≡ 0, so the Re/Im gcd statement must NOT be issued: N=4 really does
            // carry four real defective loci, and a certificate that "proved" their absence would be
            // broken. The decline must come for THAT reason: the D-side device itself must succeed.
            Assert.Equal(52, r.TrueDiscriminantDegree);
            Assert.Equal(24, r.TrueQValuationD);
            Assert.True(r.DiscLayersCertified, "the D-side multi-prime device must succeed at N=4");
            Assert.False(r.Certified, "the N=4 real disc must fail the complex-disc certificate closed");
            Assert.False(r.NoRealNonzeroCoalescence);
        }
    }

    [Fact(DisplayName = "N=6 certificate: no real q ≠ 0 carries any repeated Λ-root of F_res, window-free, both parities")]
    [Trait("Category", "SLOW_EVEN_DISC")]
    public void N6_Certificate_NoRealNonzeroCoalescence()
    {
        // ~13 s per parity (26 s both), measured quiet 2026-08-10 (the D device's 499-prime
        // sampling is ~8 s; the rest is the exact bivariate builds, once in the device and once
        // per scout prime)
        foreach (bool rOdd in new[] { false, true })
        {
            var r = FoldResultantCertificate.CertifyDiscReImGcd(6, rOdd, log: _out.WriteLine);
            _out.WriteLine($"N=6 {(rOdd ? "R-odd " : "R-even")}: p={r.CertPrime}, deg D = {r.TrueDiscriminantDegree}, "
                           + $"v_q = {r.TrueQValuationD}, deg Re = {r.DegRe} (v {r.VRe}), deg Im = {r.DegIm} (v {r.VIm}), "
                           + $"gcd deg = {r.GcdDeg}, sampled = {r.PrimesSampled} (bound {r.LcDivisorBoundD})");
            Assert.True(r.Certified, "the multi-prime device must certify deg/v and find a complex-seeing prime");
            // the scout's pinned numbers, now with the degree- and valuation-preservation halves
            // closed; Re/Im pinned SIDE BY SIDE (a side swap would survive the max/min gates alone)
            Assert.Equal(926, r.TrueDiscriminantDegree);
            Assert.Equal(536, r.TrueQValuationD);
            Assert.Equal((926, 536, 925, 537), (r.DegRe, r.VRe, r.DegIm, r.VIm));
            Assert.Equal(0, r.GcdDeg);
            Assert.True(r.NoRealNonzeroCoalescence);
        }
    }

    [Fact(DisplayName = "N=8 scout: disc layers of the F_80 residual, first prime, R-even (the dim-112 cost measured)")]
    [Trait("Category", "SLOW_EVEN_DISC8")]
    public void N8_Scout_DiscLayers()
    {
        // ~9 min, measured quiet 2026-08-10 (the dim-112 bivariate Berkowitz build + one prime's
        // ~6100-node sampling); its own category so SLOW_EVEN_DISC stays a ~1 min suite. R-even
        // only for the same reason; the R-odd conjugate twin was measured identical once the same
        // day (even-N conjugation, as gated at N=6) and is not re-run here.
        var (degD, vD, layers) = FoldResultantCertificate.DiscLayersAtNthPrime(8, rOdd: false, nth: 0);
        _out.WriteLine($"N=8 R-even: deg_q D = {degD}, v_q = {vD}, layers [{string.Join(", ", layers)}]");
        // measured 2026-08-10, first split prime (mod-p reading; deg/v since certified by
        // N8_Certificate_DiscMultiplicity). The design fork this scout opened (the simple layer
        // as the only carrier of a no-defective statement, the 0.6165 diabolic assumed a real
        // root of the doubled layer) DISSOLVED 2026-08-12: the diabolic's on-axis reading was
        // a grid artifact, the full-D gcd is trivial, see N8_Certificate_NoRealNonzeroCoalescence.
        Assert.Equal(6086, degD);
        Assert.Equal(3730, vD);
        Assert.Equal(new[] { 376, 990 }, layers);
        int aggregate = 0;
        for (int m = 0; m < layers.Length; m++) aggregate += (m + 1) * layers[m];
        Assert.True(degD - vD == aggregate,
            $"layer aggregation identity broken: {degD} − {vD} != Σ (m+1)·deg = {aggregate}");
        Assert.True(layers.Length <= 2,
            $"a multiplicity-≥3 disc layer at N=8 ({layers.Length} layers): the β-exotic reading fails at even N");
    }

    [Fact(DisplayName = "N=8 scout: F_res parity reads exact over ℤ[i], R-even (realness / q-evenness / Λ-evenness)")]
    [Trait("Category", "SLOW_EVEN_DISC8")]
    public void N8_Scout_SturmScoutParityReads()
    {
        // Step 0 of the (then-open, since-closed) N=8 blocker (arc sideways_spin_ladder): the (1,3) landing runs over ℤ and
        // declines on a complex F_res, so before any ℚ(i) split design the three structure reads
        // are gated EXACTLY on the ℤ[i] bivariate (no prime, no eigensolver). Expectations are
        // DERIVED, and a failure here is a design-level finding, not a code bug: the bipartite
        // gauge acts cross-sector at (1,2)@N=8 ((p+q)(N−1) = 21 odd), no self-fold fixes (1,2)
        // (self-folded blocks are (p, N/2) = (p,4)), and the committed sibling prior is the
        // genuinely complex (1,2)@N=6 sector residual. R-even only; the R-odd sector is the
        // exact conjugate twin at even N (σ(R-even) = conj(σ(R-odd)), gated at N=6), so its
        // reads are the same by conjugation. Cost is the dim-112 exact bivariate build.
        var r = FoldResultantCertificate.WeightSturmScout(8, 1, 2, rOdd: false, timingNodes: 0);
        _out.WriteLine($"N8 (1,2) R-even: real={r.ResidualIsReal} even={r.EvenInShiftedLambda} " +
                       $"qEven={r.EvenInQ} atDeg={r.AtDeg} resDeg={r.ResDeg} " +
                       $"mDeg={r.MDeg} fDeg={r.FDeg} vQf={r.VQf} square={r.FIsPerfectSquare} " +
                       $"uDeg={r.UDeg} uReal={r.UIsReal} compId={r.CompositionIdentityAtNodes}");
        Assert.Equal(80, r.ResDeg);
        Assert.Equal(32, r.AtDeg);   // 32 + 80 = 112 = the sector dimension
        Assert.False(r.ResidualIsReal, "F_res at (1,2)@N=8 expected coefficient-complex: " +
            "coefficient realness has only the GAUGE as a source within the sector, and the character " +
            "(p+q)(N−1) = 21 is odd, so it is cross-sector; " +
            "if this fails, F_res is REAL and the whole ℚ(i) design collapses to the ℤ route");
        Assert.False(r.EvenInQ, "per-sector q-evenness needs the gauge WITHIN the sector; " +
            "at (1,2)@N=8 it swaps the sectors, so the mirrored-grid halving is not licensed");
        Assert.False(r.EvenInShiftedLambda, "no (Λ+2N)² structure here: the within-block " +
            "self-fold is the only KNOWN source of a Λ-reflection about −2N and (1,2) is not a " +
            "(p, N/2) block at N=8, so no G and no M-descent exist. Nor is the fold SUFFICIENT, " +
            "it having to be composed " +
            "with the gauge, i.e. the character (p+q)(N−1) even: measured at the fold-fixed " +
            "(1,4)@N=8, character 35 odd, where the Λ-evenness is False (Disc14N8ScoutTests)");
        Assert.Equal(-1, r.MDeg);
        Assert.Equal(-1, r.UDeg);
        Assert.False(r.CompositionIdentityAtNodes, "no G, no identity: the sentinel is false/n-a");
    }

    [Fact(DisplayName = "N=8 probe: the D-device's prime cap and degree bound, no sampling (the E0 cost probe)")]
    [Trait("Category", "SLOW_EVEN_DISC8")]
    public void N8_Probe_DiscDeviceBound()
    {
        // E0 of the N=8 design (2026-08-12): the certified D-device samples (cap + 1) primes,
        // each a full (dBound + 25)-node mod-p disc pass, so the cap decides feasibility
        // BEFORE anything runs for days. Measured on the weight path (first run of this probe,
        // ~9.5 min: the fRes build is ≤ 67 s, the q→∞ leading form + Qc SquarefreeLayers at
        // deg 80 dominate — the same one-time pre-loop step that dominates the hardwired
        // 9-min layer scout, so the per-prime marginal is NOT 8 min and is measured
        // separately): m_D = 76 on this normalization, dBound = 6168, a BOUND and not
        // attained (the measured deg_q D = 6086, deficit 82; the hardwired (1,2) path
        // normalizes differently, m_D = 117 there, and its bound 6086 is attained).
        var p8 = FoldResultantCertificate.WeightDiscDeviceBoundProbe(8, 1, 2, rOdd: false,
                                                                     timingPrimes: 2);
        _out.WriteLine($"N=8 R-even: dim={p8.BlockDim} atDeg={p8.AtDeg} resDeg={p8.ResDeg} " +
                       $"m_D={p8.MD} dBound={p8.DBound} cap={p8.LcDivisorBoundD} " +
                       $"normD bits={p8.NormDBits} rowF bits={p8.RowFBits} " +
                       $"ms/prime={p8.MsPerPrime:F0} => device sampling ~" +
                       $"{p8.MsPerPrime * (p8.LcDivisorBoundD + 1) / 60000.0:F0} min + pre-loop");
        Assert.True(p8.MsPerPrime > 0.0, "timing primes ran");
        Assert.Equal(112, p8.BlockDim);
        Assert.Equal(80, p8.ResDeg);
        Assert.Equal(76, p8.MD);
        Assert.Equal(6168, p8.DBound);
        Assert.True(p8.DBound >= 6086, "dBound must bound the measured deg_q D");
        Assert.Equal(3165, p8.LcDivisorBoundD);
        Assert.Equal(47460, p8.NormDBits);

        var p6 = FoldResultantCertificate.WeightDiscDeviceBoundProbe(6, 1, 2, rOdd: false);
        _out.WriteLine($"N=6 R-even: dim={p6.BlockDim} atDeg={p6.AtDeg} resDeg={p6.ResDeg} " +
                       $"m_D={p6.MD} dBound={p6.DBound} cap={p6.LcDivisorBoundD} " +
                       $"normD bits={p6.NormDBits}");
        Assert.True(p6.DBound >= 926, "the N=6 bound must bound the committed deg_q D = 926");
    }

    [Fact(DisplayName = "N=8 certificate: deg/v of D certified by the weight D-device, layers at a both-attaining prime, R-even")]
    [Trait("Category", "SLOW_EVEN_DISC8_DEVICE")]
    public void N8_Certificate_DiscMultiplicity()
    {
        // T1 of the N=8 design (E0-priced 2026-08-12: pre-loop ~9.5 min + 3166 primes x
        // 545 ms ~ 29 min => ~40 min; its own category so SLOW_EVEN_DISC8 stays ~20 min).
        // Upgrades the one-prime scout read to certificate grade: deg_q D and v_q certified
        // past the lc-divisor cap, layers read at a both-attaining prime, and the measured
        // MaxDiscMultiplicity bounds the true maximum from ABOVE (mod-p roots only merge),
        // so 2 here closes every mult->=3 hideout at (1,2)@N=8 outright. R-even only: one
        // parity suffices by the derived twin D_odd = conj-coeff(D_even) (gauge swap at
        // character 21 odd composed with the within-sector conj-antiparity).
        var r = FoldResultantCertificate.WeightCertifyDiscMultiplicity(8, 1, 2, rOdd: false,
                                                                       log: _out.WriteLine);
        _out.WriteLine($"N=8 R-even: deg D = {r.TrueDiscriminantDegree}, v_q = {r.TrueQValuationD}, " +
                       $"cap = {r.LcDivisorBoundD}, sampled = {r.PrimesSampled}, " +
                       $"layer prime = {r.LayerPrime}, layers [{string.Join(", ", r.DiscLayerDegrees)}], " +
                       $"max mult = {r.MaxDiscMultiplicity}");
        Assert.True(r.DiscLayersCertified, "the multi-prime device must certify deg/v past the cap");
        Assert.Equal(6086, r.TrueDiscriminantDegree);
        Assert.Equal(3730, r.TrueQValuationD);
        Assert.Equal(new[] { 376, 990 }, r.DiscLayerDegrees);
        Assert.Equal(2, r.MaxDiscMultiplicity);
    }

    [Theory(DisplayName = "N=8 certificate: no real q ≠ 0 carries any repeated Λ-root of F_res, window-free")]
    [InlineData(false)]
    [InlineData(true)]
    [Trait("Category", "SLOW_EVEN_DISC8_DEVICE")]
    public void N8_Certificate_NoRealNonzeroCoalescence(bool rOdd)
    {
        // The N=8 design's decisive experiment, first read 2026-08-12 R-even: GcdDeg = 0.
        // No real q != 0 carries ANY repeated Lambda-root of the (1,2)@N=8 RESIDUAL,
        // defective or diabolic, window-free: the N=6 theorem's exact sibling (the AT
        // factor keeps its own semisimple degeneracies, untouched by this statement).
        // The scan record's on-axis reading of the diabolic near q = 0.6165 was the
        // cell-0.002 band zoom's grid read, never an exact fact; this certificate refutes
        // it exactly, the coalescence itself lying off the axis, so the design fork's
        // simple-layer branch was never needed and the planned landing + derivative-gcd
        // count are retired unbuilt. Also pinned before the run: the Y1 conj-antiparity
        // conj-coeff(D)(q) = D(-q) forces the parity PATTERN (Re even-even, Im odd-odd);
        // the certified deg/v = 6086/3730 pin DegRe/VRe exactly, and DegIm = 6085 /
        // VIm = 3731 add the adjacent coefficients' non-vanishing (the committed N=6
        // sibling pattern: 926/536/925/537). Both parities run (the N=6 precedent; the
        // twin D_odd = conj-coeff(D_even) predicts identical numbers, Im flipping sign,
        // and the run verifies rather than relies on it).
        var r = FoldResultantCertificate.CertifyDiscReImGcd(8, rOdd, log: _out.WriteLine);
        _out.WriteLine($"N=8 {(rOdd ? "R-odd " : "R-even")}: p={r.CertPrime}, deg D = {r.TrueDiscriminantDegree}, " +
                       $"v_q = {r.TrueQValuationD}, deg Re = {r.DegRe} (v {r.VRe}), " +
                       $"deg Im = {r.DegIm} (v {r.VIm}), gcd deg = {r.GcdDeg}, " +
                       $"sampled = {r.PrimesSampled} (bound {r.LcDivisorBoundD})");
        Assert.True(r.DiscLayersCertified, "the D-side multi-prime device must certify deg/v");
        Assert.True(r.Certified, "a complex-seeing certifying prime must exist (the disc is genuinely complex)");
        Assert.Equal(6086, r.TrueDiscriminantDegree);
        Assert.Equal(3730, r.TrueQValuationD);
        Assert.Equal((6086, 3730, 6085, 3731), (r.DegRe, r.VRe, r.DegIm, r.VIm));
        Assert.Equal(0, r.GcdDeg);
        Assert.True(r.NoRealNonzeroCoalescence);
    }

    [Fact(DisplayName = "N=6 scout: disc layers of the F_32 residual, first prime, both parities (layer identity gated)")]
    [Trait("Category", "SLOW_EVEN_DISC")]
    public void N6_Scout_DiscLayers()
    {
        foreach (bool rOdd in new[] { false, true })
        {
            var (degD, vD, layers) = FoldResultantCertificate.DiscLayersAtNthPrime(6, rOdd, 0);
            _out.WriteLine($"N=6 {(rOdd ? "R-odd " : "R-even")}: deg_q D = {degD}, v_q = {vD}, layers [{string.Join(", ", layers)}]");
            int aggregate = 0;
            for (int m = 0; m < layers.Length; m++) aggregate += (m + 1) * layers[m];
            Assert.True(degD - vD == aggregate,
                $"layer aggregation identity broken: {degD} − {vD} != Σ (m+1)·deg = {aggregate}");
            Assert.True(layers.Length <= 2,
                $"a multiplicity-≥3 disc layer at N=6 ({layers.Length} layers) — the β-exotic reading fails at even N");
        }
    }
}
