using System.Numerics;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Knowledge;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>JDefectLightMigrationClaim (the jdefect axis's light-migration mechanism):
/// (a) the claim's own battery (identity at δJ ∈ {0, 0.1}, kernel zeros, partner pair,
/// nonzero migration); (b) an independent dense check at δJ = 0.05, a point the battery
/// does not touch, with the FULL greedy pairing and per-pair complementarity; (c) registry
/// wiring (typed parents AbsorptionTheoremClaim + F1PalindromeIdentity, Tier1Derived).
/// Python twin at N=4 and N=5: <c>simulations/jdefect_light_migration.py</c>.</summary>
public class JDefectLightMigrationClaimTests
{
    private const int N = 3;
    private const double Gamma = 0.05;
    private static readonly int Dim = 1 << N;

    private static double Light(ComplexVector v)
    {
        double num = 0.0, den = 0.0;
        for (int x = 0; x < v.Count; x++)
        {
            var c = v[x];
            double p = c.Real * c.Real + c.Imaginary * c.Imaginary;
            num += BitOperations.PopCount((uint)((x / Dim) ^ (x % Dim))) * p;
            den += p;
        }
        return num / den;
    }

    // ------------------------------------------------------------------
    // (a) the claim's battery
    // ------------------------------------------------------------------

    [Fact]
    public void Claim_BatteryAllPass()
    {
        var claim = new JDefectLightMigrationClaim();
        Assert.NotEmpty(claim.Cases);
        foreach (var c in claim.Cases)
            Assert.True(c.Passes, $"battery case '{c.Name}': expected {c.Expected}, got {c.Actual}");
    }

    // ------------------------------------------------------------------
    // (b) independent dense check at δJ = 0.05 (not a battery point)
    // ------------------------------------------------------------------

    [Fact]
    public void Identity_KernelAndComplementarity_HoldAtIntermediateDeltaJ()
    {
        var axis = DimensionAxis.JDefect(N, gamma: Gamma, defectBond: 0, deltaJMax: 0.1, points: 3);
        var L = PauliDephasingDissipator.BuildZ(axis.Hamiltonian(0.05), axis.GammaPerSite);
        var evd = L.Evd();
        var vals = evd.EigenValues.ToArray();
        var right = evd.EigenVectors;

        var light = new double[vals.Length];
        for (int k = 0; k < vals.Length; k++)
        {
            light[k] = Light(right.Column(k));
            double resid = Math.Abs(vals[k].Real + 2.0 * Gamma * light[k]);
            Assert.True(resid < 1e-9,
                $"mode {k} at δJ=0.05: |Re λ + 2γ·light| = {resid:E3} (λ = {vals[k]})");
        }

        // Kernel: exactly N+1 modes at λ ≈ 0, all dark (no migration, U(1) protection).
        int kernelCount = 0;
        for (int k = 0; k < vals.Length; k++)
            if (vals[k].Magnitude < 1e-8) // far below the smallest nonzero |λ| ≈ 0.086 at N=3
            {
                kernelCount++;
                Assert.True(light[k] < 1e-9, $"kernel mode {k} carries light {light[k]:E3}");
            }
        Assert.Equal(N + 1, kernelCount);

        // Full greedy pairing λ_s + λ_f ≈ −2Nγ with per-pair light complementarity.
        var target = new Complex(-2.0 * N * Gamma, 0.0);
        var unpaired = Enumerable.Range(0, vals.Length).ToList();
        while (unpaired.Count > 0)
        {
            int i = unpaired[0];
            unpaired.RemoveAt(0);
            int bestPos = 0;
            double best = double.MaxValue;
            for (int p = 0; p < unpaired.Count; p++)
            {
                double r = (vals[i] + vals[unpaired[p]] - target).Magnitude;
                if (r < best) { best = r; bestPos = p; }
            }
            int j = unpaired[bestPos];
            unpaired.RemoveAt(bestPos);
            Assert.True(best < 1e-8,
                $"pairing incomplete at δJ=0.05: |λ_{i} + λ_{j} + 2Nγ| = {best:E3}");
            double comp = Math.Abs(light[i] + light[j] - N);
            Assert.True(comp < 1e-8,
                $"pair ({i}, {j}) at δJ=0.05: |light_s + light_f − N| = {comp:E3}");
        }
    }

    // ------------------------------------------------------------------
    // (c) registration
    // ------------------------------------------------------------------

    [Fact]
    public void BuildDefault_ContainsJDefectLightMigrationClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<JDefectLightMigrationClaim>());
    }

    [Fact]
    public void Claim_TierIsTier1Derived()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.Equal(Tier.Tier1Derived, registry.Get<JDefectLightMigrationClaim>().Tier);
    }

    [Fact]
    public void Claim_Ancestors_ContainBothTypedParents()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var ancestors = registry.AncestorsOf<JDefectLightMigrationClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
        Assert.Contains(typeof(F1PalindromeIdentity), ancestors);
    }
}
