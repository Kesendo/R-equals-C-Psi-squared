using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the pair of mirrors on the gamma axis (adopted 2026-07-21, the F134/F139
// arc's home-side move). The two involutions of the watching parameter, s: gamma -> -gamma (the
// gain turn) and s0: the anti-watch turn (agreement watched), chain through the exact generator
// identity L_anti(gamma) = L(-gamma) - 2*sigma*Id (sigma = sum_l gamma_l): the Hamiltonian leg
// never sees gamma, and per cell the turned rate is -2*sum_agree gamma_l = +2*sum_differ gamma_l
// - 2*sigma. The trajectory wears the shift as the scalar veil rho_anti(t) = e^(-2*sigma*t) *
// rho_gain(t). On rate functions s: r -> -r and s0: r -> -r - 2*sigma compose to the translation
// r -> r + 2*sigma: two mirrors make the translation by the full price, the infinite dihedral --
// F134's two-mirror shape (s: mu -> -mu, s0: mu -> 22 - mu, step 22) on the home gamma axis. A
// DIFFERENT object from the F1 fold (Pi flips the sign of L_H; this keeps H and flips only gamma).
public class GammaFoldTests
{
    static readonly World W = new();

    // the generator identity holds per cell at machine precision, for a non-uniform site profile
    // (per-site is load-bearing, uniform gamma is the special case), and the dihedral closes:
    // the turn is its own inverse, and gain-after-turn is the translation by the full price.
    [Fact]
    public void The_Mask_Identity_And_The_Dihedral_Are_Exact()
    {
        var fold = new GammaFold(W, 4);
        var ml = fold.MaskLaws();
        Assert.True(ml.WorstIdentity < 1e-12, $"L_anti(g) = L(-g) - 2*sigma*Id per cell: {ml.WorstIdentity:E1}");
        Assert.True(ml.WorstInvolution < 1e-12, $"s0 o s0 = id: {ml.WorstInvolution:E1}");
        Assert.True(ml.WorstTranslation < 1e-12, $"s o s0 = translation by 2*sigma: {ml.WorstTranslation:E1}");
        Assert.True(ml.Step > 0.5, $"the translation step must be O(1), not vacuous: {ml.Step:0.000}");
    }

    // the veil law rho_anti(t) = e^(-2*sigma*t) * rho_gain(t): twin RK4 from the same seed
    // (populations + a live coherence), worst entry residual over every probed tick. The measured
    // floor at N=3, dt=0.02, 50 ticks is ~3e-8 (the RK4 truncation of the scalar shift); pinned
    // with a 30x margin. The veil is not vacuous: anti and gain differ O(1) at the final tick,
    // and the gain world carries e^(+2*sigma*t) times the anti world's novelty (the trace is
    // blind to gamma in every world; the amplification lives in the coherences).
    [Fact]
    public void The_Trajectory_Wears_The_Price_Veil()
    {
        var fold = new GammaFold(W, 3);
        var rep = fold.Run(seed: 1, dt: 0.02, ticks: 50);
        Assert.True(rep.WorstVeil < 1e-6, $"rho_anti = e^(-2*sigma*t) * rho_gain: {rep.WorstVeil:E1}");
        Assert.True(rep.VertexSeparation > 0.1, $"anti and gain must differ O(1): {rep.VertexSeparation:0.000}");
        double expected = Math.Exp(2.0 * fold.Sigma * 1.0);
        Assert.True(Math.Abs(rep.NoveltyRatio - expected) < 1e-4,
            $"gain/anti novelty must be e^(+2*sigma*t) = {expected:0.000}: {rep.NoveltyRatio:0.000}");
    }

    // the gain flip is load-bearing: the veil against the NORMAL (+gamma) world instead of the
    // gain world misses at O(1) from below.
    [Fact]
    public void The_Gain_Flip_Is_LoadBearing()
    {
        var fold = new GammaFold(W, 3);
        var rep = fold.Run(seed: 1, dt: 0.02, ticks: 50);
        Assert.True(rep.BrokenFlip > 1e-3, $"the wrong (unflipped) world must break the veil O(1): {rep.BrokenFlip:E1}");
    }

    // the identity never touches H: the veil survives the ZZ bond unchanged (the shift is a
    // scalar, the flip acts on the dissipator alone).
    [Fact]
    public void The_Veil_Holds_With_The_ZZ_Bond()
    {
        var fold = new GammaFold(W, 3, zz: 1.0);
        var rep = fold.Run(seed: 1, dt: 0.02, ticks: 50);
        Assert.True(rep.WorstVeil < 1e-6, $"veil with ZZ: {rep.WorstVeil:E1}");
    }

    // the cross-dock with the Lattice (the scout's third find): the one-sided X^N reading of the
    // normal world equals the gain evolution of the one-sided seed in the price veil,
    // X^N rho(t) = e^(-2*sigma*t) * gain-run(X^N rho0). The committed Lattice bridge composed
    // with the veil law, checked directly.
    [Fact]
    public void The_Watched_World_Read_Through_The_Complement_Is_The_Gain_World_In_The_Veil()
    {
        var fold = new GammaFold(W, 3);
        double worst = fold.ReadThroughVeil(seed: 1, dt: 0.02, ticks: 50);
        Assert.True(worst < 1e-6, $"X^N rho(t) = veil * gain-read(t): {worst:E1}");
    }

    // ontology: the fold's own outputs.
    [Fact]
    public void Ontology_Own_Carries_The_Fold_Readings()
    {
        var fold = new GammaFold(W, 3);
        Assert.Contains("mask-identity", fold.Own);
        Assert.Contains("veil", fold.Own);
        Assert.Contains("dihedral", fold.Own);
    }
}
