using System.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Ptf;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Tests.Ptf;

/// <summary>The s*-boundary surface as a computed object (F86 Edge 1, 2026-06-10): the slow
/// subspace's per-site light distribution, basis-free via the orthogonal projector onto the
/// range of the biorthogonal spectral projector over the slowest-rate cluster.
///
/// <para>All numbers pinned against the Python reference run (N=5 dimensionless XY chain,
/// γ-profiles summing to N, same construction as <c>simulations/light_content.py</c> /
/// <c>simulations/birth_canal_boundary_pathdependence.py</c>): uniform is sterile with the
/// distribution exactly [0.2]⁵ frozen between Q=1.5 and Q=1000 (drift ~10⁻¹⁷ in the
/// reference); flat-bulk-edge [0.25,1.5,1.5,1.5,0.25] is in the birth canal, its distribution
/// drifting [0.35034, 0.02180, 0.25572, …] → [⅓, 0, ⅓, 0, ⅓] (max per-site drift 7.76·10⁻²)
/// with the rate following 1.24829 → 4/3; the absorption cross-check rate = 2·Σ_l γ_l·light_l
/// is machine-exact even on an asymmetric degenerate cluster where the raw projector diagonal
/// (residual 7.3·10⁻²) and the PP† diagonal (3.2·10⁻²) both fail; and the projector reading
/// coincides with the per-eigenvector |v|² average wherever Degeneracy = 1.</para></summary>
public class SlowLightDistributionTests
{
    private static readonly double[] Uniform5 = { 1.0, 1.0, 1.0, 1.0, 1.0 };
    private static readonly double[] FlatBulkEdge5 = { 0.25, 1.5, 1.5, 1.5, 0.25 };
    private static readonly double[] PeakedV5 = { 0.25, 0.75, 3.0, 0.75, 0.25 };
    private static readonly double[] Asymmetric5 = { 0.3, 0.9, 1.6, 1.1, 1.1 };

    private static ComplexMatrix Liouvillian(int n, double q, double[] profile) =>
        new PostEpFlowField(n, new[] { q }, new[] { 0.0, 1.0 }, gammaProfile: profile)
            .DimensionlessLiouvillian(q);

    private static SlowLightDistribution.Reading Read(int n, double q, double[] profile) =>
        SlowLightDistribution.Compute(Liouvillian(n, q, profile), n, profile);

    // One Evd per (profile, Q), shared across the facts that read it.
    private static readonly Lazy<SlowLightDistribution.Reading> UniformLow =
        new(() => Read(5, 1.5, Uniform5));
    private static readonly Lazy<SlowLightDistribution.Reading> UniformHigh =
        new(() => Read(5, 1000.0, Uniform5));
    private static readonly Lazy<SlowLightDistribution.Reading> CanalLow =
        new(() => Read(5, 1.5, FlatBulkEdge5));
    private static readonly Lazy<SlowLightDistribution.Reading> CanalHigh =
        new(() => Read(5, 1000.0, FlatBulkEdge5));

    [Fact]
    public void Uniform_N5_Sterile_DistributionIsQFrozenAndEquipartitioned()
    {
        // Sterile zone: the slow subspace's light distribution is Q-invariant. For the uniform
        // chain it is exactly equipartitioned, [0.2]⁵, at Q=1.5 and at Q=1000 (Python reference
        // drift 5.6·10⁻¹⁷). The slowest rate is the closed form 2γ = 2 on the odd rail.
        var lo = UniformLow.Value;
        var hi = UniformHigh.Value;

        Assert.Equal(2.0, lo.Rate, 7);
        Assert.Equal(50, lo.Degeneracy);
        Assert.Equal(1.0, lo.TotalLight, 7);
        for (int l = 0; l < 5; l++)
        {
            Assert.Equal(0.2, lo.PerSiteLight[l], 9);
            Assert.True(Math.Abs(hi.PerSiteLight[l] - lo.PerSiteLight[l]) < 1e-9,
                $"sterile distribution must be Q-frozen, site {l} drifted by " +
                $"{Math.Abs(hi.PerSiteLight[l] - lo.PerSiteLight[l]):E3}");
        }
    }

    [Fact]
    public void FlatBulkEdge_N5_BirthCanal_DistributionDriftsWithQ()
    {
        // Birth canal: the Hamiltonian redistributes the slow subspace's light with Q and the
        // rate follows (Absorption). Pinned drift: center site 0.25572 → ⅓ (7.76·10⁻²), bulk
        // sites 0.02180 → 0, edges 0.35034 → ⅓; rate 1.24829 → 4/3.
        var lo = CanalLow.Value;
        var hi = CanalHigh.Value;

        Assert.Equal(10, lo.Degeneracy);
        Assert.Equal(1.2482919, lo.Rate, 6);
        Assert.Equal(0.3503416, lo.PerSiteLight[0], 6);
        Assert.Equal(0.0217967, lo.PerSiteLight[1], 6);
        Assert.Equal(0.2557233, lo.PerSiteLight[2], 6);
        Assert.Equal(0.0217967, lo.PerSiteLight[3], 6);
        Assert.Equal(0.3503416, lo.PerSiteLight[4], 6);

        Assert.Equal(4.0 / 3.0, hi.Rate, 5);
        Assert.Equal(1.0 / 3.0, hi.PerSiteLight[2], 5);

        double maxDrift = 0.0;
        for (int l = 0; l < 5; l++)
            maxDrift = Math.Max(maxDrift, Math.Abs(hi.PerSiteLight[l] - lo.PerSiteLight[l]));
        Assert.True(maxDrift > 0.05,
            $"canal distribution must drift with Q (pinned 7.76e-2), got {maxDrift:E3}");
    }

    [Fact]
    public void AbsorptionCrossCheck_MachineExact_InBothZones()
    {
        // rate = 2·Σ_l γ_l·light_l exactly (block-triangularity on the invariant subspace).
        // Python reference residuals at Q=1.5: 4.2·10⁻¹⁴ (uniform), 1.3·10⁻¹⁴ (flat-bulk-edge).
        Assert.True(UniformLow.Value.AbsorptionResidual < 1e-9,
            $"uniform absorption residual {UniformLow.Value.AbsorptionResidual:E3}");
        Assert.True(CanalLow.Value.AbsorptionResidual < 1e-9,
            $"canal absorption residual {CanalLow.Value.AbsorptionResidual:E3}");
    }

    [Fact]
    public void AbsorptionCrossCheck_MachineExact_OnAsymmetricDegenerateCluster()
    {
        // The discriminating case (degeneracy 20, no left-right symmetry): the raw biorthogonal
        // projector diagonal fails the cross-check at 7.3·10⁻² and goes negative, the PP†
        // diagonal fails at 3.2·10⁻²; the orthogonal-projector reading is exact (1.7·10⁻¹⁴ in
        // the Python reference) and stays in [0, 1] per site.
        var r = Read(5, 1.5, Asymmetric5);
        Assert.Equal(20, r.Degeneracy);
        Assert.Equal(1.6620717, r.Rate, 6);
        Assert.Equal(0.0254708, r.PerSiteLight[2], 6);
        Assert.True(r.AbsorptionResidual < 1e-9, $"absorption residual {r.AbsorptionResidual:E3}");
        foreach (double lk in r.PerSiteLight)
            Assert.InRange(lk, 0.0, 1.0);
    }

    [Fact]
    public void PeakedV_N5_LightAvoidsTheStrongCenter()
    {
        // The sterile peaked-V profile: the slow subspace puts its one light quantum on the
        // weak edges and none on the strong center, light = [¼, ¼, 0, ¼, ¼], rate 1.
        var r = Read(5, 1.5, PeakedV5);
        Assert.Equal(1.0, r.Rate, 7);
        Assert.Equal(0.25, r.PerSiteLight[0], 9);
        Assert.Equal(0.25, r.PerSiteLight[1], 9);
        Assert.Equal(0.0, r.PerSiteLight[2], 9);
        Assert.Equal(0.25, r.PerSiteLight[3], 9);
        Assert.Equal(0.25, r.PerSiteLight[4], 9);
        Assert.True(r.AbsorptionResidual < 1e-9, $"absorption residual {r.AbsorptionResidual:E3}");
    }

    [Fact]
    public void Projector_MatchesEigenvectorAverage_WhereDegeneracyIsOne()
    {
        // At Degeneracy = 1 the orthogonal projector is |m⟩⟨m|/‖m‖², so the projector light
        // must equal the per-eigenvector Rayleigh light, the regime where the old averaged
        // carrier was already exact. N=2, profile [0.7, 1.3], Q=0.2: a unique real slowest
        // mode (rate 0.16696972).
        const int n = 2;
        var profile = new[] { 0.7, 1.3 };
        var liouvillian = Liouvillian(n, 0.2, profile);
        var r = SlowLightDistribution.Compute(liouvillian, n, profile);

        Assert.Equal(1, r.Degeneracy);
        Assert.Equal(0.1669697220, r.Rate, 9);

        // The eigenvector-averaged reading, computed the pre-2026-06-10 way.
        int d = 1 << n;
        var evd = liouvillian.Evd();
        var vals = evd.EigenValues;
        double maxRe = double.NegativeInfinity;
        int slowest = -1;
        for (int k = 0; k < vals.Count; k++)
            if (vals[k].Magnitude > 1e-7 && vals[k].Real > maxRe) { maxRe = vals[k].Real; slowest = k; }
        var v = evd.EigenVectors.Column(slowest);
        var averaged = new double[n];
        double norm2 = 0.0;
        for (int x = 0; x < d * d; x++)
        {
            double w = v[x].Real * v[x].Real + v[x].Imaginary * v[x].Imaginary;
            norm2 += w;
            int diff = (x / d) ^ (x % d);
            for (int l = 0; l < n; l++)
                if (((diff >> (n - 1 - l)) & 1) != 0) averaged[l] += w;
        }
        for (int l = 0; l < n; l++)
            Assert.True(Math.Abs(r.PerSiteLight[l] - averaged[l] / norm2) < 1e-10,
                $"site {l}: projector {r.PerSiteLight[l]:E12} vs averaged {averaged[l] / norm2:E12}");
    }

    [Fact]
    public void Compute_Throws_WhenNoNonKernelMode()
    {
        // H = 0, γ = 0-free is not constructible here; instead a 1-site zero Liouvillian:
        // everything is kernel, nothing to read.
        var zero = ComplexMatrix.Build.Dense(4, 4);
        Assert.Throws<InvalidOperationException>(() =>
            SlowLightDistribution.Compute(zero, 1, new[] { 1.0 }));
    }

    [Fact]
    public void Compute_RejectsMismatchedShapes()
    {
        var liouvillian = Liouvillian(2, 1.0, new[] { 1.0, 1.0 });
        Assert.Throws<ArgumentException>(() =>
            SlowLightDistribution.Compute(liouvillian, 3, new[] { 1.0, 1.0, 1.0 }));
        Assert.Throws<ArgumentException>(() =>
            SlowLightDistribution.Compute(liouvillian, 2, new[] { 1.0 }));
    }
}
