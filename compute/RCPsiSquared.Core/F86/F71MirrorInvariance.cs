using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Resonance;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 Tier-1 derived: F71 spatial-mirror invariance of per-bond Q_peak. The
/// chain-mirror symmetry F71 (R |b₀…b_{N-1}⟩ = |b_{N-1}…b₀⟩) pairs bond b with bond N−2−b.
/// Numerical observation: <c>Q_peak(b) = Q_peak(N−2−b)</c> bit-exactly for every bond pair
/// across c=2..3 N=5..7.
///
/// <para>Algebraic origin: the per-bond K-observable <c>K_b(t) = 2·Re⟨ρ(t)|S_kernel|∂ρ/∂J_b⟩</c>
/// is invariant under F71 spatial mirror because (a) the Z-dephasing dissipator is F71-symmetric,
/// (b) the Dicke probe is F71-symmetric, (c) the spatial-sum kernel S is F71-symmetric, and
/// (d) ∂L/∂J_b transforms to ∂L/∂J_{N−2−b}. Hence K_b(Q, t) = K_{N−2−b}(Q, t) as functions of
/// Q and t, and Q_peak(b) = Q_peak(N−2−b) follows.</para>
///
/// <para>Substructure refinement: even within the F71-orbit-grouped Interior bonds at
/// c=2 N=6, the central bond b=2 (=b_{N-2-2}=b_2, self-paired) gives a different Q_peak
/// (1.440) than the flanking bonds b=1, b=3 (1.648). The "Endpoint vs Interior" dichotomy
/// is a first approximation; the full structure is per-F71-orbit, with mid-chain orbits
/// showing further variation.</para>
/// </summary>
public sealed class F71MirrorInvariance : Claim
{
    public F71MirrorInvariance()
        : base("F71 spatial-mirror invariance: Q_peak(b) = Q_peak(N−2−b)",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F71.md (chain mirror) + 2026-05-03 F86 numerical verification")
    { }

    public override string DisplayName => "F71 mirror: Q_peak(b) = Q_peak(N−2−b)";

    public override string Summary =>
        $"per-bond F71 spatial-mirror pairing of Q_peak — bit-exact verified across c=2..3, N=5..7 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("algebraic origin",
                summary: "F71-symmetric L_D + S_kernel + Dicke probe → K_b(Q, t) = K_{N−2−b}(Q, t)");
            yield return new InspectableNode("substructure",
                summary: "even Interior is per-F71-orbit, not uniform — central bond differs from flanking");
            yield return new InspectableNode("numerical verification",
                summary: "5+ cases (c=2 N=5..7, c=3 N=5..6), all bond pairs match to machine precision");
        }
    }

    /// <summary>Verify F71 mirror invariance for a measured KCurve. Returns the maximum
    /// deviation |Q_peak(b) − Q_peak(N−2−b)| across all bond pairs.</summary>
    public static double MaxMirrorDeviation(KCurve curve)
    {
        int numBonds = curve.NumBonds;
        double maxDev = 0;
        for (int b = 0; b < numBonds / 2; b++)
        {
            var peakA = curve.PeakAtBond(b);
            var peakB = curve.PeakAtBond(numBonds - 1 - b);
            double dev = Math.Abs(peakA.QPeak - peakB.QPeak);
            if (dev > maxDev) maxDev = dev;
        }
        return maxDev;
    }
}
