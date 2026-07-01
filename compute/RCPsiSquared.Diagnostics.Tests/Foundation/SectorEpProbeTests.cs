using System.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Guards <see cref="SectorEpProbe"/>: the core Multi-Sector Monodromy Census primitive
/// that, for a sector at a coupling q0, builds the RAW (p, q̃) block (Task 1's
/// <see cref="RCPsiSquared.Core.BlockSpectrum.SectorBlock"/>), finds the coalescing eigenvalue pair,
/// and classifies it defective-vs-diabolic via the trusted
/// <see cref="RCPsiSquared.Core.Numerics.EpCharacter"/>.
///
/// <para>The two q-values below are the review's decisive from-below result
/// (<c>_msm_ep_character.py</c>) on the N=4 (p=1, q̃=2) block: q≈0.658983 is the SEMISIMPLE DIABOLIC
/// (a (3q⁴+q²−1) root, SVD nullity 2, the pair merges at λ≈−4+1.318i and stays independent) and
/// q≈0.460212 is a DEFECTIVE EP (a P₁₀ real root, SVD nullity 1, a braid-carrying Jordan block).
/// Reproducing this in C# is the load-bearing validation that the Task-1 convention + EpCharacter
/// agree with the exact math.</para></summary>
public class SectorEpProbeTests
{
    [Theory]
    [InlineData(0.658983, 0.0, "Diabolic")]   // (3q^4+q^2-1) root: semisimple node, silent
    [InlineData(0.460212, 0.0, "Defective")]  // P_10 real root: braid-carrying EP
    public void Probe_ClassifiesOcticCoalescence(double qre, double qim, string kind)
    {
        var r = SectorEpProbe.Probe(N: 4, p: 1, qTilde: 2, q0: new System.Numerics.Complex(qre, qim));
        // A coalescence guard, NOT the load-bearing assertion (that is the Kind below). The threshold is
        // 1e-2, not 1e-3, because the two EP types respond DIFFERENTLY to the 6-digit rounding of these
        // reference q-values: the diabolic node splits LINEARLY (gap ≈ 4e-7 at q≈0.658983), but the
        // defective EP splits as √(q−q*), so a ~6e-7 offset of q from the exact P₁₀ root is amplified to a
        // gap ≈ 3e-3 at q≈0.460212 (the √-law itself, confirmed from-below by the census scout). Both are
        // still far below the next eigenvalue (diabolic third-eigenvalue distance ≈ 2.4e-2), so 1e-2
        // cleanly certifies "a coalescing pair, well inside the next eigenvalue" for both types.
        Assert.True(r.MinGap < 1e-2, $"expected a coalescence, gap={r.MinGap}");
        Assert.Equal(kind, r.Kind.ToString());
    }
}
