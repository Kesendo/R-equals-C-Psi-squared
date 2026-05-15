using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Resonance;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F86;

/// <summary>Reconnaissance probe for direction (γ') of
/// <c>PolarityInheritanceLink.PendingDerivationNote</c>: the asymptotic structural
/// constant <c>1/g_eff_E + 1/g_eff_I ≈ 0.937</c> (witnessed but not derived).
///
/// <para>The witnessed 0.937 was computed from per-bond-CLASS g_eff means (Endpoint ≈ 1.74,
/// Interior ≈ 2.81). PolarityInheritanceLink itself flags (its N≥9 finding) that the
/// "Interior class" orbit-mixes: it is not a single F71 orbit. So 0.937 may be a real
/// sum-rule or a coarse-average artefact, and the per-class data cannot tell them apart.
/// This probe fetches the CLEAN per-F71-orbit g_eff from the <b>K_CC_pr observable</b> via
/// the existing typed <see cref="PerF71OrbitKTable"/> claim (which groups the per-bond
/// <c>C2HwhmRatio</c> witnesses by F71 orbit, c=2), with
/// <c>g_eff = 2·BareDoubledPtfXPeak / Q_peak</c> per the PolarityInheritanceLink
/// composition (<c>r_Q = BareDoubledPtfXPeak·Q_EP − 2 = 4.39382/g_eff − 2</c>, with
/// <c>Q_peak = 2 + r_Q</c>), and dumps the closure candidates so direction (γ') can be
/// decided on un-orbit-mixed data instead of guessed.</para>
///
/// <para>(An earlier draft of this probe used <c>PerF71OrbitLQPeakTable</c> by mistake —
/// the JW-track L(Q) lens, a different observable. The per-orbit Q_peak there did not
/// match the per-class K_CC_pr anchor, which caught the error. <see cref="PerF71OrbitKTable"/>
/// is the K_CC_pr-observable sibling.)</para>
///
/// <para>This is a DUMP probe, not a verification: it asserts only that the per-orbit data
/// builds. Orbits whose Q_peak escaped the default scan grid (no finite resonance) are
/// flagged and excluded from the closure sum. The (γ') verdict is read off the table.</para>
///
/// <para><b>Verdict (2026-05-14): (γ') refuted.</b> On the clean per-F71-orbit data the
/// per-orbit pairwise 1/g_eff sums drift (no constant ~0.937) and Σ(1/g_eff) scales
/// additively with orbit count, so the per-class 0.937 is an orbit-mixing artefact, not a
/// closure. Recorded in <c>PolarityInheritanceLink</c> direction (γ'); structural reason in
/// <c>docs/proofs/PROOF_F86B_OBSTRUCTION.md</c> § "The diagnosis": g_eff is not a spectral
/// primitive, so a closure among 1/g_eff values cannot hold.</para>
/// </summary>
public class F86GammaEffClosureProbeTests
{
    private readonly ITestOutputHelper _out;
    public F86GammaEffClosureProbeTests(ITestOutputHelper @out) => _out = @out;

    private static CoherenceBlock C2Block(int N) => new(N: N, n: 1, gammaZero: 0.05);

    // g_eff = 2·BareDoubledPtfXPeak / Q_peak  — PolarityInheritanceLink.ClosedFormCompositionNote:
    // r_Q = BareDoubledPtfXPeak·Q_EP − 2 = 4.39382/g_eff − 2, with Q_peak = 2 + r_Q.
    private static double GEffFromQPeak(double qPeak) =>
        2.0 * C2HwhmRatio.BareDoubledPtfXPeak / qPeak;

    [Fact]
    public void GammaEffClosure_PerF71Orbit_Recon()
    {
        var defaultGrid = ResonanceScan.DefaultQGrid();
        int[] ns = { 5, 6, 7, 8 };

        _out.WriteLine("per-F71-orbit g_eff recon — K_CC_pr observable (PerF71OrbitKTable), c=2,");
        _out.WriteLine($"default Q-grid [{defaultGrid[0]:F2}, {defaultGrid[^1]:F2}] / {defaultGrid.Length} pts.");
        _out.WriteLine($"g_eff = 2·BareDoubledPtfXPeak / Q_peak = {2.0 * C2HwhmRatio.BareDoubledPtfXPeak:F5} / Q_peak");
        _out.WriteLine("per-class anchor (PolarityInheritanceLink): g_eff_E ≈ 1.74, g_eff_I ≈ 2.81,");
        _out.WriteLine("  1/g_eff_E + 1/g_eff_I ≈ 0.937 (EmpiricalSumQPeakAsymptote = 4.12)");
        _out.WriteLine("");

        foreach (int n in ns)
        {
            var table = PerF71OrbitKTable.Build(C2Block(n));
            Assert.NotEmpty(table.OrbitWitnesses);

            _out.WriteLine($"=== N={n}  ({table.OrbitWitnesses.Count} F71 orbits) ===");
            double sumInvClean = 0.0;
            int cleanCount = 0;
            foreach (var w in table.OrbitWitnesses)
            {
                bool escaped = w.IsEscaped(defaultGrid);
                double gEff = GEffFromQPeak(w.QPeak);
                string bonds = w.Orbit.IsSelfPaired
                    ? $"{{b={w.Orbit.BondA}}}"
                    : $"{{b={w.Orbit.BondA}<->{w.Orbit.BondB}}}";
                _out.WriteLine(
                    $"  orbit {bonds,-14} Q_peak={w.QPeak:F4}  g_eff={gEff:F4}  1/g_eff={1.0 / gEff:F4}"
                    + $"  HWHM/Q*={w.HwhmLeftOverQPeak:F4}"
                    + (escaped ? "  [ESCAPED: Q_peak pinned to grid edge, no finite resonance]" : ""));
                if (!escaped) { sumInvClean += 1.0 / gEff; cleanCount++; }
            }
            _out.WriteLine($"  Sum(1/g_eff) clean ({cleanCount}/{table.OrbitWitnesses.Count})  = {sumInvClean:F4}");
            _out.WriteLine($"  mean(1/g_eff) clean             = {(cleanCount > 0 ? sumInvClean / cleanCount : 0):F4}");
            _out.WriteLine("");
        }

        _out.WriteLine("VERDICT read off the table above: is there a constant closure across N");
        _out.WriteLine("on the per-orbit data, or was 0.937 a per-class orbit-mixing artefact?");
    }
}
