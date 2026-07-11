using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the site-resolved rate cross-check (Program run mode "concentrator"), the
// continuous-Lindblad comparator for experiments/IBM_CONCENTRATOR_RELOADED.md stage 7a. It pins:
//   (1) the J convention: hopping J = 2 reproduces the design's unit-Pauli XY bond XX+YY;
//   (2) the site-resolved watching: siteGammas prices a coherence at -2 sum_l gamma_l*(bit l of i^j);
//   (3) the central prediction: the continuous slope(MP)-slope(E) matches the design's flown -0.07325/step
//       up to Trotter error, with the far-sink leakage null-consistent and coh_0(t=8) ~ 0.4554;
//   (4) that the longitudinal ZZ (Heisenberg vs XY) barely moves the payload-coherence contrast.
// All Restless additions are optional-parameter overloads; the existing constructor path is unchanged
// (guarded by the pre-existing RestlessTests, which run with siteGammas=null and zz=0).
public class ConcentratorTests
{
    static readonly World W = new();

    // ---- (1) the J convention, pinned from below ------------------------------------------------
    // The design bond is exp(-i*0.05*(XX+YY+ZZ)); its XY part is the UNIT-Pauli XX+YY. Restless builds
    // H = J*(s+_a s-_b + h.c.) = (J/2)(XX+YY), so J = 2 reproduces the design's XY hopping. Check it on the
    // cleanest witness: a single excitation on N=2, gamma=0. The single-excitation block is [[0,J],[J,0]],
    // so |01> oscillates as P(|01>) = cos^2(J t). At J=2 that is cos^2(2t) -- the design's XY frequency.
    [Fact]
    public void HoppingJ_Two_Reproduces_The_Unit_Pauli_XY_Bond()
    {
        var r = new Restless(W, 2, j: 2.0, gamma: 0.0);
        r.Seed(0b01);                                   // one excitation, |01><01|
        const double t = 0.3; const int steps = 600;
        for (int s = 0; s < steps; s++) r.Step(t / steps);
        double p01 = r[0b01, 0b01].Real;
        Assert.Equal(Math.Cos(2 * t) * Math.Cos(2 * t), p01, 6);   // cos^2(2t): J=2 = unit-Pauli XY
    }

    // ---- (2) the site-resolved watching -----------------------------------------------------------
    // siteGammas makes the watching per-site: only a coherence that DIFFERS in a watched bit is priced,
    // at -2*gamma_l per watched differing bit. A lone coherence across the sink bit decays exp(-2 gamma t);
    // a coherence across an UNWATCHED bit does not decay at all. J=0 isolates the watching from the H.
    [Fact]
    public void SiteGammas_Prices_Only_Coherences_Across_The_Watched_Bit()
    {
        var g = new double[3]; g[1] = 0.5;              // watch site 1 only
        // coherence |000><010>: differs in bit 1 (watched) -> decays exp(-2*0.5*t)
        var watched = new Restless(W, 3, j: 0.0, gamma: 0.0, siteGammas: g);
        watched.SeedCoherence(0b000, 0b010, 1.0);
        // coherence |000><001>: differs in bit 0 (UNwatched) -> no decay
        var unwatched = new Restless(W, 3, j: 0.0, gamma: 0.0, siteGammas: g);
        unwatched.SeedCoherence(0b000, 0b001, 1.0);
        const double dt = 0.05; const int steps = 10; double t = dt * steps;
        for (int s = 0; s < steps; s++) { watched.Step(dt); unwatched.Step(dt); }
        Assert.Equal(Math.Exp(-2 * 0.5 * t), watched[0b000, 0b010].Magnitude, 5);   // priced at -2 gamma
        Assert.Equal(1.0, unwatched[0b000, 0b001].Magnitude, 10);                   // far site: untouched
    }

    // ---- the arm runner (mirrors Program's "concentrator" mode) -----------------------------------
    const double JHop = 2.0, GSink = 0.5, StepT = 0.05;
    const int SubPerStep = 50;
    static readonly int[] Grid = { 1, 2, 3, 4, 6, 8 };

    static double[] Arm(int n, int payload, int? sinkSite, double zz)
    {
        var g = new double[n];
        if (sinkSite is int ss) g[ss] = GSink;
        var r = new Restless(W, n, JHop, 0.0, Topology.Chain(n), siteGammas: g, zz: zz);
        int pIdx = 1 << payload;
        r.Seed(0, 0.5); r.Seed(pIdx, 0.5); r.SeedCoherence(0, pIdx, 0.5);   // |+> on payload, |0> elsewhere
        double dt = StepT / SubPerStep;
        var coh = new List<double>();
        int gi = 0;
        for (int step = 1; step <= Grid[^1]; step++)
        {
            for (int s = 0; s < SubPerStep; s++) r.Step(dt);
            if (gi < Grid.Length && step == Grid[gi]) { coh.Add(2.0 * r[0, pIdx].Magnitude); gi++; }
        }
        return coh.ToArray();
    }

    static double SlopeLnR(double[] ca, double[] c0)
    {
        int m = Grid.Length; double sx = 0, sy = 0, sxx = 0, sxy = 0;
        for (int i = 0; i < m; i++)
        {
            double x = Grid[i], y = Math.Log(ca[i] / c0[i]);   // x = step -> slope PER STEP
            sx += x; sy += y; sxx += x * x; sxy += x * y;
        }
        return (m * sxy - sx * sy) / (m * sxx - sx * sx);
    }

    // the seed IS the design payload: coh(0) = 2|rho[0, payload]| = 1 for |+> on the payload.
    [Fact]
    public void The_Payload_Seed_Has_Unit_Initial_Coherence()
    {
        var r = new Restless(W, 5, JHop, 0.0, Topology.Chain(5), siteGammas: new double[5], zz: 1.0);
        r.Seed(0, 0.5); r.Seed(1 << 2, 0.5); r.SeedCoherence(0, 1 << 2, 0.5);
        Assert.Equal(1.0, 2.0 * r[0, 1 << 2].Magnitude, 12);
    }

    // ---- (3) the central prediction, N=5 ----------------------------------------------------------
    // The continuous Heisenberg (zz=1) slope(MP)-slope(E), independent from-below, must match the design's
    // recorded flown value -0.07325/step up to Trotter error, and reproduce coh_0(t=8) ~ 0.4554. The far
    // sink is null-consistent (|slope(E)| << the contrast), and the on-payload sink prices FASTER.
    [Fact]
    public void N5_Central_Contrast_Matches_The_Preregistered_Design()
    {
        var c0 = Arm(5, 2, null, 1.0);
        var cE = Arm(5, 2, 0, 1.0);      // far edge sink (design arm E)
        var cMP = Arm(5, 2, 2, 1.0);     // on-payload sink (design arm MP)
        double sE = SlopeLnR(cE, c0), sMP = SlopeLnR(cMP, c0), diff = sMP - sE;

        Assert.Equal(0.45543, c0[^1], 4);                       // coh_0(t=8 steps=0.4) ~ 0.4554 (design anchor)
        Assert.Equal(-0.073618, diff, 5);                       // the continuous contrast, pinned
        Assert.True(Math.Abs(diff - (-0.073249)) < 0.0025,      // within the O((J*dt)^2) Trotter bound of the flown value
            $"continuous contrast {diff:F6} too far from the flown -0.073249");
        Assert.True(Math.Abs(sE) < 1e-3, $"far-sink leakage not null-consistent: {sE:F6}");   // design -0.00029
        Assert.True(sMP < sE - 0.05, $"on-payload sink must price faster: sMP={sMP:F6} sE={sE:F6}");
    }

    // ---- (4) the ZZ (Heisenberg vs XY) contribution is tiny for this observable -------------------
    // The payload coherence is dominated by the XY hopping; the longitudinal ZZ moves the contrast by
    // only ~1e-4/step. So the XY-only Restless already gives a faithful cross-check, and Heisenberg the
    // more faithful one; either lands well inside the Trotter bound.
    [Fact]
    public void ZZ_Barely_Moves_The_Payload_Contrast()
    {
        double diffZZ = SlopeLnR(Arm(5, 2, 2, 1.0), Arm(5, 2, null, 1.0))
                      - SlopeLnR(Arm(5, 2, 0, 1.0), Arm(5, 2, null, 1.0));
        double diffXY = SlopeLnR(Arm(5, 2, 2, 0.0), Arm(5, 2, null, 0.0))
                      - SlopeLnR(Arm(5, 2, 0, 0.0), Arm(5, 2, null, 0.0));
        Assert.True(Math.Abs(diffZZ - diffXY) < 5e-4,
            $"ZZ should barely move the contrast: zz={diffZZ:F6} xy={diffXY:F6}");
    }

    // ---- the bonus: the contrast persists past the N=8 spectrum wall ------------------------------
    // At N=9 the full Liouvillian is out of reach, but a single-excitation-vs-vacuum coherence is
    // block-local (100 alive cells). The site contrast survives, slightly larger than N=5 (the far edge
    // is farther, so the leakage is even smaller), and the far-sink leakage is essentially zero.
    [Fact]
    public void N9_Site_Contrast_Persists_Past_The_Wall()
    {
        var c0 = Arm(9, 4, null, 1.0);
        var cE = Arm(9, 4, 0, 1.0);
        var cMP = Arm(9, 4, 4, 1.0);
        double sE = SlopeLnR(cE, c0), sMP = SlopeLnR(cMP, c0), diff = sMP - sE;
        Assert.True(diff < -0.07 && diff > -0.078, $"N=9 contrast off expected band: {diff:F6}");
        Assert.True(Math.Abs(sE) < 1e-4, $"N=9 far-sink leakage should be ~0: {sE:F6}");
        Assert.True(sMP < sE - 0.05, $"N=9 on-payload sink must price faster: sMP={sMP:F6} sE={sE:F6}");
    }
}
