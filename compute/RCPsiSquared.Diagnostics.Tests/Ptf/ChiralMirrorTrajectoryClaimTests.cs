using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Propagation;
using RCPsiSquared.Core.States;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Diagnostics.Ptf;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Tests.Ptf;

/// <summary>The PTF K₁ chiral mirror (PROOF_PTF_CHIRAL_MIRROR_RATE_LAW): (a) the Step-1
/// algebra at N=4 (K₁HK₁ = −H, K₁VK₁ = −V exact, dissipator K₁-invariant, sine-mode
/// mapping); (b) the site-wise trajectory identity P_i(t; φ_k) = P_i(t; φ_{N+1−k}) through
/// the Core spectral propagator, including the U(1) minus-sign branch; (c) registry wiring
/// (typed parent ChiralKClaim) and the claim's pure-algebra battery.</summary>
public class ChiralMirrorTrajectoryClaimTests
{
    private const int N = 4;
    private const double Gamma = 0.05;
    private const double DeltaJ = 0.1;

    private static ComplexMatrix K1() =>
        PauliString.SiteOp(N, 1, PauliLetter.Z) * PauliString.SiteOp(N, 3, PauliLetter.Z);

    /// <summary>The defect V = δJ·½(X₀X₁ + Y₀Y₁) on bond (0,1).</summary>
    private static ComplexMatrix Defect() =>
        PauliHamiltonian.Bilinear(N, new[] { new Bond(0, 1, 1.0) },
            new (PauliLetter, PauliLetter, Complex)[]
            {
                (PauliLetter.X, PauliLetter.X, DeltaJ / 2.0),
                (PauliLetter.Y, PauliLetter.Y, DeltaJ / 2.0),
            }).ToMatrix();

    // ------------------------------------------------------------------
    // (a) Step-1 algebra at N = 4
    // ------------------------------------------------------------------

    [Fact]
    public void K1_AnticommutesWithChainAndDefect_AndCommutesWithDephasers()
    {
        var chain = new ChainSystem(N, J: 1.0, GammaZero: Gamma);
        var H = chain.BuildHamiltonian();
        var V = Defect();
        var k1 = K1();

        Assert.True((k1 * H * k1 + H).FrobeniusNorm() < 1e-12,
            $"K₁HK₁ ≠ −H: residual {(k1 * H * k1 + H).FrobeniusNorm():E3}");
        Assert.True((k1 * V * k1 + V).FrobeniusNorm() < 1e-12,
            $"K₁VK₁ ≠ −V: residual {(k1 * V * k1 + V).FrobeniusNorm():E3}");

        for (int l = 0; l < N; l++)
        {
            var zl = PauliString.SiteOp(N, l, PauliLetter.Z);
            Assert.True((k1 * zl - zl * k1).FrobeniusNorm() < 1e-12,
                $"[K₁, Z_{l}] ≠ 0");
        }
    }

    [Fact]
    public void K1_MapsSineModesToMirrorPartners_WithoutSign()
    {
        var k1 = K1();
        for (int k = 1; k <= N; k++)
        {
            var mapped = k1 * BondingMode.Build(N, k);
            var partner = BondingMode.Build(N, N + 1 - k);
            Assert.True((mapped - partner).L2Norm() < 1e-12,
                $"K₁ψ_{k} ≠ +ψ_{N + 1 - k}: deviation {(mapped - partner).L2Norm():E3}");
        }
    }

    // ------------------------------------------------------------------
    // (b) Trajectory identity at N = 4 (γ₀ = 0.05, δJ = 0.1 on bond (0,1))
    // ------------------------------------------------------------------

    /// <summary>Row-major vec(|ψ⟩⟨ψ|), matching PauliDephasingDissipator's convention.</summary>
    private static ComplexVector VecOuter(ComplexVector psi)
    {
        int d = psi.Count;
        var v = ComplexVector.Build.Dense(d * d);
        for (int a = 0; a < d; a++)
            for (int b = 0; b < d; b++)
                v[a * d + b] = psi[a] * Complex.Conjugate(psi[b]);
        return v;
    }

    private static double[] SitePurities(ComplexVector rhoVec)
    {
        int d = 1 << N;
        var rho = ComplexMatrix.Build.Dense(d, d, (a, b) => rhoVec[a * d + b]);
        var purities = new double[N];
        for (int i = 0; i < N; i++)
        {
            var red = PartialTrace.Of(rho, N, new[] { i });
            purities[i] = (red * red).Trace().Real;
        }
        return purities;
    }

    private static IReadOnlyList<double[]> PurityTrajectory(ComplexMatrix L, ComplexVector psi, double[] taus)
    {
        var states = SpectralPropagator.EvolveStateVectors(L, VecOuter(psi), taus);
        return states.Select(SitePurities).ToArray();
    }

    [Fact]
    public void SitePurities_Coincide_ForMirrorPairStates()
    {
        var chain = new ChainSystem(N, J: 1.0, GammaZero: Gamma);
        var L = PauliDephasingDissipator.BuildZ(chain.BuildHamiltonian() + Defect(),
            Enumerable.Repeat(Gamma, N).ToArray());
        var taus = new[] { 0.5, 2.0 };

        var pK = PurityTrajectory(L, BondingMode.PairState(N, 1), taus);
        var pPartner = PurityTrajectory(L, BondingMode.PairState(N, N), taus);

        for (int t = 0; t < taus.Length; t++)
            for (int i = 0; i < N; i++)
                Assert.True(Math.Abs(pK[t][i] - pPartner[t][i]) < 1e-9,
                    $"site {i}, t={taus[t]}: P(φ_1)={pK[t][i]:E12} vs P(φ_{N})={pPartner[t][i]:E12}");
    }

    [Fact]
    public void SitePurities_Coincide_ForUOneMinusSignBranch()
    {
        // χ = (|vac⟩ − |ψ_N⟩)/√2 is the even-sublattice image of φ_1; the U(1) phase
        // e^{iπN̂} maps it to φ_N, so its purity trajectories must match φ_N's.
        var chain = new ChainSystem(N, J: 1.0, GammaZero: Gamma);
        var L = PauliDephasingDissipator.BuildZ(chain.BuildHamiltonian() + Defect(),
            Enumerable.Repeat(Gamma, N).ToArray());
        var taus = new[] { 0.5, 2.0 };

        var chi = -BondingMode.Build(N, N);
        chi[0] += Complex.One;
        chi = chi / Math.Sqrt(chi.ConjugateDotProduct(chi).Real);

        var pChi = PurityTrajectory(L, chi, taus);
        var pPartner = PurityTrajectory(L, BondingMode.PairState(N, N), taus);

        for (int t = 0; t < taus.Length; t++)
            for (int i = 0; i < N; i++)
                Assert.True(Math.Abs(pChi[t][i] - pPartner[t][i]) < 1e-9,
                    $"site {i}, t={taus[t]}: P(χ)={pChi[t][i]:E12} vs P(φ_{N})={pPartner[t][i]:E12}");
    }

    // ------------------------------------------------------------------
    // (c) Registration + battery
    // ------------------------------------------------------------------

    [Fact]
    public void BuildDefault_ContainsChiralMirrorTrajectoryClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<ChiralMirrorTrajectoryClaim>());
    }

    [Fact]
    public void Claim_TierIsTier1Derived()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.Equal(Tier.Tier1Derived, registry.Get<ChiralMirrorTrajectoryClaim>().Tier);
    }

    [Fact]
    public void Claim_Ancestors_ContainChiralKClaim()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var ancestors = registry.AncestorsOf<ChiralMirrorTrajectoryClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(ChiralKClaim), ancestors);
    }

    [Fact]
    public void Claim_BatteryAllPass()
    {
        var claim = new ChiralMirrorTrajectoryClaim();
        Assert.NotEmpty(claim.Cases);
        foreach (var c in claim.Cases)
            Assert.True(c.Passes, $"battery case '{c.Name}': expected {c.Expected}, got {c.Actual}");
    }
}
