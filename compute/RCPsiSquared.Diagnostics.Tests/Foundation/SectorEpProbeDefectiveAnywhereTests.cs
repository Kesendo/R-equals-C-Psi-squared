using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Guards <see cref="SectorEpProbe.ProbeDefectiveAnywhere"/>, the AT-aware "defective-anywhere" probe
/// that unblocks the Multi-Sector Monodromy Census at N ≥ 5. Where the validated
/// <see cref="SectorEpProbe.Probe"/> keys on the GLOBAL closest eigenvalue pair, at N ≥ 5 many RAW joint-popcount
/// blocks carry a PERMANENT semisimple (AT) degeneracy (a gap ~1e-15 pair at EVERY q); the global-closest probe
/// locks onto that permanent pair, reads <see cref="EpCharacter.EpKind.Diabolic"/>, and MASKS the residual
/// braid-carrying √-EP (gap ~2e-4 at a different λ). The defective-anywhere probe clusters ALL near-coalescences
/// and characterizes each, so it reaches the masked residual EP.
///
/// <para><b>Both tests probe at the sector's own precisely-scanned reference loci</b>
/// (<see cref="ReferenceDefectiveLoci.For"/>), exactly as the census does. The often-quoted "0.460212" is only the
/// 6-digit DISPLAY of the N=4 locus 0.46021217346191406; the full-precision value is what puts the residual √-EP
/// gap at ~5e-4 (below the 1e-3 cluster-link threshold), and a truncated hand-seed lands ~3e-3 off the true EP.</para>
///
/// <para><b>N=4 agreement (the strict-improvement guard).</b> At N=4 NO block has BOTH a permanent AT pair AND a
/// residual √-EP, so the defective-anywhere probe returns the SAME verdict as the global-closest
/// <see cref="SectorEpProbe.Probe"/>: the (1,2) Klein-orbit braid sectors read HasDefective (Probe reads
/// Defective), the node/core sectors do not (Probe reads Diabolic). This proves it is a strict improvement, not a
/// different instrument.</para>
///
/// <para><b>N=5 cross-fold unmask (the new capability).</b> By the F89d cross-fold
/// (<see cref="CrossFoldSimilarityWitness"/>) a (1,2) DEFECTIVE EP at q maps to a (1, N−2) = (1,3) DEFECTIVE EP at
/// q̄ with the IDENTICAL character and coalescence gap. That residual √-EP is present in the RAW (1,3) block but
/// MASKED by the permanent AT degeneracy, so <see cref="SectorEpProbe.Probe"/> under-reports it as Diabolic. The
/// defective-anywhere probe finds it, re-activating the census's N=5 cross-fold gate.</para></summary>
public class SectorEpProbeDefectiveAnywhereTests
{
    /// <summary>The (1,2) octic's real reference DEFECTIVE loci at chain length N, ascending. Precisely scanned
    /// (full double precision), so the residual √-EP gap sits below the probe's cluster-link threshold; picking
    /// the smallest real locus is deterministic. Falls back to the smallest locus overall if none is real.</summary>
    private static Complex FirstReferenceLocus(int n)
    {
        IReadOnlyList<Complex> loci = ReferenceDefectiveLoci.For(n);
        var reals = loci.Where(q => Math.Abs(q.Imaginary) < 1e-6).OrderBy(q => q.Real).ToList();
        return reals.Count > 0 ? reals[0] : loci.OrderBy(q => q.Real).ThenBy(q => q.Imaginary).First();
    }

    // N=4 must AGREE with the validated global-closest probe: braid → true, node/core → false. Probed at the
    // (1,2) octic's first real reference DEFECTIVE locus (like the census; NOT a 6-digit hand-seed).
    [Theory]
    [InlineData(1, 2, true)]    // braid sector (SE↔DE): has the residual √-EP
    [InlineData(2, 3, true)]    // braid orbit (Klein-four image of (1,2))
    [InlineData(1, 3, false)]   // node sector: only a permanent semisimple pair, no √-EP
    [InlineData(2, 2, false)]   // half-filled core: braid-free
    public void DefectiveAnywhere_N4_AgreesWithVerdict(int p, int qTilde, bool expected)
    {
        var q0 = FirstReferenceLocus(4);

        var da = SectorEpProbe.ProbeDefectiveAnywhere(4, p, qTilde, q0);
        Assert.Equal(expected, da.HasDefective);

        // The strict-improvement guard: at N=4 the AT-aware probe returns the SAME verdict as the validated
        // global-closest Probe (no N=4 block has BOTH a permanent pair AND a residual √-EP), so HasDefective ⟺
        // Probe reads Defective. At N ≥ 5 they diverge on purpose (the unmask fact below).
        var closest = SectorEpProbe.Probe(4, p, qTilde, q0);
        Assert.Equal(closest.Kind == EpCharacter.EpKind.Defective, da.HasDefective);
    }

    // The load-bearing new capability: at N=5 the (1,3) cross-fold partner of (1,2) carries a Defective EP at the
    // conjugate locus (F89d), MASKED behind the block's permanent AT degeneracy. The global-closest Probe reads
    // Diabolic there; the AT-aware probe unmasks the residual √-EP.
    [Fact(DisplayName = "N=5 (1,3) cross-fold √-EP unmasked behind the permanent AT degeneracy")]
    [Trait("Category", "SLOW_MSM")]
    public void DefectiveAnywhere_N5_Unmasks_1_3_CrossFold()
    {
        // A real (1,2) N=5 defective locus q0 (For(5) has two: q ≈ 0.6209, 1.0776). The F89d cross-fold sends a
        // (1,2) EP at q to a (1, N−2) = (1,3) EP at q̄ (= q0 for a real locus, up to ~1e-16).
        var q0 = FirstReferenceLocus(5);
        var partnerQ = Complex.Conjugate(q0);

        var da = SectorEpProbe.ProbeDefectiveAnywhere(5, 1, 3, partnerQ);
        Assert.True(da.HasDefective,
            $"cross-fold (1,3) defective EP must be found by the AT-aware probe (minDefGap={da.MinDefectiveGap:E3})");
        Assert.True(da.MinDefectiveGap < 1e-3,
            $"the unmasked cluster must be a genuine coalescence (gap={da.MinDefectiveGap:E3})");

        // The mask this probe sees through: the global-closest Probe locks onto the permanent semisimple AT pair
        // (gap ~1e-15) and reports Diabolic, NOT Defective — it MISSES the residual braid. This divergence from
        // the N=4 agreement above is exactly the new capability.
        var masked = SectorEpProbe.Probe(5, 1, 3, partnerQ);
        Assert.NotEqual(EpCharacter.EpKind.Defective, masked.Kind);
    }
}
