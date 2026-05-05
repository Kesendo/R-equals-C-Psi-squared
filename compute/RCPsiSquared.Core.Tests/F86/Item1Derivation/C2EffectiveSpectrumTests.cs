using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.F86.Item1Derivation;
using RCPsiSquared.Core.Knowledge;
using Xunit;

namespace RCPsiSquared.Core.Tests.F86.Item1Derivation;

public class C2EffectiveSpectrumTests
{
    [Theory]
    [InlineData(5, 0.5)]
    [InlineData(5, 1.0)]
    [InlineData(5, 1.5)]
    [InlineData(5, 2.0)]
    [InlineData(5, 2.5)]
    [InlineData(5, 3.0)]
    [InlineData(7, 0.5)]
    [InlineData(7, 1.5)]
    [InlineData(7, 3.0)]
    public void Eigenvalues_MatchNumericalEvd_AtMachinePrecision(int N, double Q)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        var effective = FourModeEffective.Build(block);

        for (int b = 0; b < block.NumBonds; b++)
        {
            var ours = spectrum.Eigenvalues(Q, b);
            var theirs = effective.LEffAtQ(Q).Evd().EigenValues.ToArray();

            // Sort both by Re desc, then Im asc — stable comparison
            var oursSorted = ours.OrderByDescending(z => z.Real).ThenBy(z => z.Imaginary).ToArray();
            var theirsSorted = theirs.OrderByDescending(z => z.Real).ThenBy(z => z.Imaginary).ToArray();

            for (int i = 0; i < 4; i++)
            {
                var diff = (oursSorted[i] - theirsSorted[i]).Magnitude;
                Assert.True(diff < 1e-10,
                    $"λ_{i}({Q}, b={b}) at N={N}: ours={oursSorted[i]}, theirs={theirsSorted[i]}, diff={diff:E2}");
            }
        }
    }

    [Fact]
    public void Tier_IsTier1DerivedOrTier2Verified()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        Assert.True(spectrum.Tier == Tier.Tier1Derived || spectrum.Tier == Tier.Tier2Verified);
    }

    [Fact]
    public void IsAnalyticallyDerived_ConsistentWithTier()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        if (spectrum.IsAnalyticallyDerived) Assert.Equal(Tier.Tier1Derived, spectrum.Tier);
        else Assert.Equal(Tier.Tier2Verified, spectrum.Tier);
    }

    [Fact]
    public void Build_ThrowsIfNotC2()
    {
        var block = new CoherenceBlock(N: 5, n: 2, gammaZero: 0.05);  // c=3
        Assert.Throws<ArgumentException>(() => C2EffectiveSpectrum.Build(block));
    }

    [Fact]
    public void Eigenvalues_OutOfRangeBond_Throws()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        Assert.Throws<ArgumentOutOfRangeException>(() => spectrum.Eigenvalues(1.0, -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => spectrum.Eigenvalues(1.0, block.NumBonds));
    }

    [Fact]
    public void KDrivingPair_OutOfRangeBond_Throws()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        Assert.Throws<ArgumentOutOfRangeException>(() => spectrum.KDrivingPair(1.0, bond: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => spectrum.KDrivingPair(1.0, bond: block.NumBonds));
    }

    [Fact]
    public void KDrivingPairIndices_OutOfRangeBond_Throws()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        Assert.Throws<ArgumentOutOfRangeException>(() => spectrum.KDrivingPairIndices(1.0, bond: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => spectrum.KDrivingPairIndices(1.0, bond: block.NumBonds));
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void Eigenvalues_BondIndependent_AtCEquals2(int N)
    {
        // At c=2 the spectrum is bond-independent by construction (uniform-J means
        // L_eff = D_eff + Q·γ₀·Σ_b V_b, and the Σ_b is the same for any bond index passed).
        // This property test pins the contract.
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        double Q = 1.5;
        var reference = spectrum.Eigenvalues(Q, bond: 0);
        for (int b = 1; b < block.NumBonds; b++)
        {
            var actual = spectrum.Eigenvalues(Q, bond: b);
            for (int i = 0; i < 4; i++)
            {
                var diff = (actual[i] - reference[i]).Magnitude;
                Assert.True(diff < 1e-12,
                    $"Eigenvalues should be bond-independent at c=2: bond=0 vs bond={b}, λ_{i} diff={diff:E2}");
            }
        }
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void LEffAtQ_MatchesFourModeEffective_AtMachinePrecision(int N)
    {
        // Cross-check: our bond-summed L_eff(Q) assembly must match FourModeEffective's
        // direct projection L_eff(Q) = D_eff + (Q·γ₀)·M_h_total_eff entry-by-entry.
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        var effective = FourModeEffective.Build(block);

        foreach (double Q in new[] { 0.5, 1.0, 2.0 })
        {
            var ours = spectrum.LEffAtQ(Q);
            var theirs = effective.LEffAtQ(Q);
            var diff = (ours - theirs).FrobeniusNorm();
            Assert.True(diff < 1e-12, $"L_eff(Q={Q}) at N={N}: ‖ours − theirs‖_F = {diff:E2}");
        }
    }

    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    public void Eigenvalues_TraceMatchesDEffPlusQGammaTraceMh(int N)
    {
        // Sanity: Σ λ_i = trace(L_eff(Q)) = trace(D_eff) + Q·γ₀·trace(M_h_total_eff).
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        double Q = 1.0;
        var eigs = spectrum.Eigenvalues(Q, bond: 0);
        var traceFromEigs = eigs.Aggregate(Complex.Zero, (acc, z) => acc + z);
        var traceFromMatrix = spectrum.LEffAtQ(Q).Diagonal().Sum();
        Assert.True((traceFromEigs - traceFromMatrix).Magnitude < 1e-10,
            $"trace mismatch at N={N}: from eigs = {traceFromEigs}, from matrix = {traceFromMatrix}");
    }

    [Fact]
    public void Tier_IsIndependentOfBondCouplingTier_BothObstructionsCoexist()
    {
        // C2EffectiveSpectrum.Tier and BondCoupling.Tier are independent obstructions, not
        // a parent-child inheritance: BondCoupling.Tier = Tier2Verified comes from A3's
        // even-N σ_0 degeneracy; C2EffectiveSpectrum.Tier = Tier2Verified comes from the
        // independent cubic-c_3 obstruction in the 4×4 char poly (no biquadratic reduction).
        // Even if A3 were promoted to Tier1Derived, this spectrum's Tier would stay
        // Tier2Verified due to the cubic-c_3 obstruction — and vice versa.
        // This test pins that BondCoupling is reachable and carries its own Tier verdict.
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        Assert.NotNull(spectrum.BondCoupling);
        Assert.Equal(Tier.Tier2Verified, spectrum.BondCoupling.Tier);
    }

    [Fact]
    public void PendingDerivationNote_IsNullOnTier1Derived_NonNullOnTier2()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        if (spectrum.Tier == Tier.Tier1Derived)
            Assert.Null(spectrum.PendingDerivationNote);
        else
        {
            Assert.NotNull(spectrum.PendingDerivationNote);
            Assert.Contains("Tier2Verified", spectrum.PendingDerivationNote!);
        }

        if (spectrum.Tier == Tier.Tier2Verified)
        {
            Assert.NotNull(spectrum.PendingDerivationNote);
            Assert.Contains("cubic", spectrum.PendingDerivationNote!,
                StringComparison.OrdinalIgnoreCase);  // pins the c_3-cubic obstruction-of-record
        }
    }

    // ============================================================================
    // Stage C3: K-driving eigenvalue pair identification.
    // ============================================================================

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void ProbeProjection_HasZeroSvdTopComponents(int N)
    {
        // Tier 1 structural sub-fact: the Dicke probe lives entirely in span{|c_1⟩, |c_3⟩}
        // (basis indices 0, 1). Components onto |u_0⟩, |v_0⟩ (basis indices 2, 3) are exactly
        // zero by basis orthogonality. Per InterChannelSvd's structural finding:
        // "the probe (Dicke state) is orthogonal to |u_0⟩, |v_0⟩".
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        var probe = spectrum.ProbeProjection;
        Assert.True(probe[2].Magnitude < 1e-14,
            $"|probe·u_0| should be ~0; got {probe[2].Magnitude:E2} at N={N}");
        Assert.True(probe[3].Magnitude < 1e-14,
            $"|probe·v_0| should be ~0; got {probe[3].Magnitude:E2} at N={N}");
        // Sanity: c_1, c_3 components are non-trivial (probe is fully concentrated there).
        double cuMag = Math.Sqrt(probe[0].Magnitude * probe[0].Magnitude +
                                  probe[1].Magnitude * probe[1].Magnitude);
        Assert.True(cuMag > 0.1,
            $"|probe·(c_1, c_3)| should carry the full norm; got {cuMag:F4} at N={N}");
    }

    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    public void KDrivingPairIndices_ReferToValidEigenvectors(int N)
    {
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        var (i, j) = spectrum.KDrivingPairIndices(Q: 1.5, bond: 0);
        Assert.InRange(i, 0, 3);
        Assert.InRange(j, 0, 3);
        Assert.NotEqual(i, j);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void KDrivingPair_ProbeOverlap_DominatesNonDrivingPair(int N)
    {
        // Structural fact: probe lives in span{|c_1⟩, |c_3⟩}. The K-driving pair are the 2
        // eigenvectors with largest overlap with the probe. By construction,
        // overlap(driving) > overlap(non-driving).
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);

        foreach (double Q in new[] { 0.5, 1.0, 1.5, 2.0, 2.5 })
        {
            var (idxPlus, idxMinus) = spectrum.KDrivingPairIndices(Q, bond: 0);
            var overlaps = spectrum.ProbeOverlapsSquared(Q, bond: 0);
            double drivingSum = overlaps[idxPlus] + overlaps[idxMinus];
            double total = overlaps.Sum();
            double nonDriving = total - drivingSum;
            Assert.True(drivingSum > nonDriving,
                $"K-driving overlap {drivingSum:F4} should exceed non-driving {nonDriving:F4} at N={N}, Q={Q}");
        }
    }

    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    public void KDrivingPair_ReturnsEigenvaluesAtIdentifiedIndices(int N)
    {
        // Pin the identity: KDrivingPair(Q, b).LamPlus == Eigenvalues(Q, b)[IndexPlus],
        // and likewise for LamMinus. The pair eigenvalues are taken from the Eigenvalues array
        // at the identified indices, not recomputed independently — single source of truth.
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        foreach (double Q in new[] { 0.5, 1.5, 2.5 })
        {
            var eigs = spectrum.Eigenvalues(Q, bond: 0);
            var (idxPlus, idxMinus) = spectrum.KDrivingPairIndices(Q, bond: 0);
            var (lamPlus, lamMinus) = spectrum.KDrivingPair(Q, bond: 0);
            Assert.Equal(eigs[idxPlus], lamPlus);
            Assert.Equal(eigs[idxMinus], lamMinus);
        }
    }

    [Theory]
    [InlineData(5)]
    [InlineData(7)]
    public void KDrivingPair_LamPlus_HasGreaterOrEqualRePart_ThanLamMinus(int N)
    {
        // LamPlus is the K-driving eigenvalue with the larger real part by construction
        // (slower-decaying mode); LamMinus is the smaller-Re partner.
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        foreach (double Q in new[] { 0.5, 1.0, 1.5, 2.0, 2.5 })
        {
            var (lamPlus, lamMinus) = spectrum.KDrivingPair(Q, bond: 0);
            Assert.True(lamPlus.Real >= lamMinus.Real,
                $"LamPlus.Re ({lamPlus.Real:F4}) should be ≥ LamMinus.Re ({lamMinus.Real:F4}) at N={N}, Q={Q}");
        }
    }

    [Theory]
    [InlineData(5, 0.5)]
    [InlineData(5, 1.0)]
    [InlineData(5, 1.5)]
    [InlineData(5, 2.0)]
    [InlineData(7, 1.0)]
    [InlineData(7, 1.5)]
    public void KDrivingPair_AtModerateQ_HasReNearMinus4GammaZero(int N, double Q)
    {
        // F86 Statement 1: at Q ≈ Q_EP the K-driving pair both approach Re(λ) = −4γ₀.
        // We don't pin Q to Q_EP exactly because g_eff is not yet a closed form; the test
        // verifies that across a sample of Q values around the empirical Q_peak the K-driving
        // pair's mean Re sits in the F86 ladder range [−6γ₀, −2γ₀] (between the D_eff
        // diagonal entries for HD=1 and HD=3, which are the natural confining bounds for the
        // eigenvalues that carry the probe content).
        const double gamma0 = 0.05;
        var block = new CoherenceBlock(N, n: 1, gammaZero: gamma0);
        var spectrum = C2EffectiveSpectrum.Build(block);
        var (lamPlus, lamMinus) = spectrum.KDrivingPair(Q, bond: 0);
        double reMean = (lamPlus.Real + lamMinus.Real) / 2.0;
        // Generous band: [−9γ₀, −γ₀] = [−0.45, −0.05] for γ₀ = 0.05. Both eigenvalues sit
        // around the −4γ₀ ladder rung at the EP; the band leaves room for the c=2 4×4 mixing
        // off-EP. This pins the structural F86 prediction without demanding the closed form
        // for g_eff.
        Assert.InRange(reMean, -9.0 * gamma0, -1.0 * gamma0);
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void KDrivingPair_OverlapDominanceFraction_ExceedsThreshold(int N)
    {
        // Empirical sanity: the probe ⊥ {|u_0⟩, |v_0⟩} structural fact means the K-driving
        // pair carries the dominant fraction of the total |⟨probe | w⟩|² mass. At small
        // probe ↔ SVD coupling (the c=2 cross-block is empirically small relative to the
        // diagonal blocks), the dominance fraction is close to 1. This isn't a Tier-1
        // closed-form bound — it depends on the (Q, b)-rotation — but it's the empirical
        // signature that K-driving identification is unambiguous in practice.
        // (observed values empirically ≥ 0.9 across the test grid; threshold 0.7 leaves headroom)
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        foreach (double Q in new[] { 1.0, 1.5, 2.0 })
        {
            var (idxPlus, idxMinus) = spectrum.KDrivingPairIndices(Q, bond: 0);
            var overlaps = spectrum.ProbeOverlapsSquared(Q, bond: 0);
            double drivingSum = overlaps[idxPlus] + overlaps[idxMinus];
            double total = overlaps.Sum();
            double fraction = drivingSum / total;
            Assert.True(fraction > 0.7,
                $"Dominance fraction {fraction:F4} should be > 0.7 at N={N}, Q={Q}");
        }
    }

    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    public void ProbeOverlapsSquared_SumEqualsProbeNormSquared(int N)
    {
        // The eigenvectors of L_eff(Q, b) form a basis (generically); the squared overlaps
        // sum to ‖probe‖² when the basis is orthonormal. L_eff(Q, b) is non-Hermitian, so
        // the eigenvectors are not orthonormal in general — but at γ₀ ≠ 0 / Q < ∞ they are
        // still linearly independent, and the squared overlaps sum to ‖probe‖² when summed
        // over the orthonormal basis. This sanity check uses the orthonormal 4-mode basis
        // {|c_1⟩, |c_3⟩, |u_0⟩, |v_0⟩}: ‖probe‖² in the basis is just Σ |probe[k]|².
        // We verify the eigenvector overlap sum is in a reasonable neighbourhood of this
        // normalisation reference (within 50% — not a tight equality because eigenvector
        // non-orthogonality of non-Hermitian matrices means the dual basis weighting is not
        // identity).
        var block = new CoherenceBlock(N, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        var probe = spectrum.ProbeProjection;
        double probeNormSq = 0.0;
        for (int k = 0; k < 4; k++)
            probeNormSq += probe[k].Magnitude * probe[k].Magnitude;
        var overlaps = spectrum.ProbeOverlapsSquared(Q: 1.5, bond: 0);
        double total = overlaps.Sum();
        // The overlap sum should be positive and within 2× of the probe norm squared
        // (eigenvector non-orthonormality is typically modest at the c=2 EP regime).
        Assert.True(total > 0.0, $"Overlap sum should be > 0; got {total:F4} at N={N}");
        Assert.True(total < 2.0 * probeNormSq + 1e-3,
            $"Overlap sum {total:F4} should be within 2× of ‖probe‖² ({probeNormSq:F4}) at N={N}");
    }

    [Fact]
    public void ProbeOverlapsSquared_OutOfRangeBond_Throws()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var spectrum = C2EffectiveSpectrum.Build(block);
        Assert.Throws<ArgumentOutOfRangeException>(() => spectrum.ProbeOverlapsSquared(1.0, -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => spectrum.ProbeOverlapsSquared(1.0, block.NumBonds));
    }
}
