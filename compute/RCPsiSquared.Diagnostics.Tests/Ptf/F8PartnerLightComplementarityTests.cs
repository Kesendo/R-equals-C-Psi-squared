using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Diagnostics.Ptf;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Tests.Ptf;

/// <summary>F8 palindrome partners carry complementary light, and the standing-wave
/// coupling between them is exactly null.
///
/// <para><b>The April null (2026-04-18)</b>, banked from
/// <c>ClaudeTasks/Archiv/RESULT_PTF_STANDING_WAVE_ROLES.md</c>: for the uniform XY chain
/// under uniform Z-dephasing the Liouvillian spectrum pairs by the palindrome
/// λ_s + λ_f = −2Σγ (the F8/F1 pairing), but the hypothesized coupling of paired partners
/// through a J-defect perturbation V_L = −i[V, ·] does not exist: at N=7,
/// max |⟨W_f|V_L|M_s⟩| = 4.43·10⁻¹⁴ over all 80×80 partner pairs, the numerical floor.
/// Structural reason: V_L conserves the U(1)×U(1) joint popcount sector (a, b) of a
/// coherence |a-sector⟩⟨b-sector|, while the palindrome pairing maps (a, b) → (N−a, N−b);
/// the two never meet except at a = b = N/2.</para>
///
/// <para><b>The light-complementarity reading (2026-06-10 fresh-eyes survey)</b>: palindrome
/// partners sit on opposite banks of the PTF near/far edge by definition of being mirror
/// partners. Per right eigenmode v of the Liouvillian, the absorption identity
///
/// <code>
///   Re λ = −2γ · light(v),   light(v) = Σ_x w(x)·|v_x|² / Σ_x |v_x|²,
///   w(x) = popcount(i ⊕ j) for the row-major coherence index x = i·d + j
/// </code>
///
/// holds exactly. Derivation: λ⟨v|v⟩ = ⟨v|L|v⟩ with L = A + D, where A = −i[H, ·] is
/// anti-Hermitian in the Hilbert-Schmidt inner product for every Hermitian H (so ⟨v|A|v⟩
/// is purely imaginary) and D is diagonal real with entries −2γ·w(x); taking real parts
/// gives Re λ·⟨v|v⟩ = −2γ·Σ_x w(x)·|v_x|². Hence λ_s + λ_f = −2Nγ ⟺
/// light_s + light_f = N: the eigenvalue palindrome and the complementary light content
/// are the same statement read in two coordinates. This is the eigenmode-expectation form
/// of the Absorption Theorem (<c>docs/proofs/PROOF_ABSORPTION_THEOREM.md</c>,
/// <see cref="RCPsiSquared.Core.Symmetry.AbsorptionTheoremClaim"/>); the Hamming-complement
/// pair-sum corollary there is its sharp per-coherence special case.</para>
///
/// <para>Verified here at N=4 (256 modes, 128 pairs) with the layer's standard builders;
/// the Python convention check <c>simulations/f8_partner_light_sanity.py</c> reproduces
/// all three facts at N=3 and N=4 at the 10⁻¹⁴ floor.</para></summary>
public class F8PartnerLightComplementarityTests
{
    private const int N = 4;
    private const double Gamma = 0.05;
    private static readonly int Dim = 1 << N;

    /// <summary>One dense eigendecomposition of the N=4 uniform XY chain Liouvillian,
    /// shared by all facts. Left covectors are the rows of R⁻¹: globally biorthogonal
    /// to the right eigenvectors with ⟨W_k|M_k⟩ = 1 by construction, robust under
    /// degenerate clusters (no per-vector normalization).</summary>
    private static readonly Lazy<(Complex[] Values, ComplexMatrix Right, ComplexMatrix Left)> Decomp =
        new(() =>
        {
            var chain = new ChainSystem(N, J: 1.0, GammaZero: Gamma);
            var evd = chain.BuildLiouvillian().Evd();
            var right = evd.EigenVectors;
            return (evd.EigenValues.ToArray(), right, right.Inverse());
        });

    /// <summary>light(v) = Σ_x w(x)·|v_x|² / Σ_x |v_x|² with w(x) = popcount(i ⊕ j) at the
    /// row-major coherence index x = i·d + j (the Hamming weight of bra-ket difference),
    /// matching <see cref="PauliDephasingDissipator"/>'s vec convention.</summary>
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

    /// <summary>Greedy multiset pairing of the spectrum by λ_a + λ_b ≈ −2Nγ. Greedy
    /// (not sort-based) because sort-based pairing is unstable under degeneracy.</summary>
    private static (List<(int A, int B)> Pairs, double WorstResidual) GreedyF8Pairs(Complex[] vals)
    {
        var target = new Complex(-2.0 * N * Gamma, 0.0);
        var unpaired = Enumerable.Range(0, vals.Length).ToList();
        var pairs = new List<(int, int)>(vals.Length / 2);
        double worst = 0.0;
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
            worst = Math.Max(worst, best);
            int j = unpaired[bestPos];
            unpaired.RemoveAt(bestPos);
            pairs.Add((i, j));
        }
        return (pairs, worst);
    }

    [Fact]
    public void AbsorptionIdentity_HoldsForEveryRightEigenmode()
    {
        // Re λ = −2γ·light(v) exactly, for all 256 modes: the Hamiltonian part is
        // anti-Hermitian (purely imaginary Rayleigh quotient) and drops out of Re λ.
        var (vals, right, _) = Decomp.Value;
        for (int k = 0; k < vals.Length; k++)
        {
            double resid = Math.Abs(vals[k].Real + 2.0 * Gamma * Light(right.Column(k)));
            Assert.True(resid < 1e-9,
                $"mode {k}: |Re λ + 2γ·light| = {resid:E3} (λ = {vals[k]})");
        }
    }

    [Fact]
    public void F8Pairing_IsComplete_AndPartnersCarryComplementaryLight()
    {
        // The XY chain is bipartite/soft, so the palindrome pairing must be COMPLETE:
        // all 256 eigenvalues pair off with λ_s + λ_f = −2Nγ. By the absorption identity
        // this is the same statement as light_s + light_f = N: mirror partners sit on
        // opposite banks of the PTF near/far edge.
        var (vals, right, _) = Decomp.Value;
        var (pairs, worst) = GreedyF8Pairs(vals);

        Assert.Equal(vals.Length / 2, pairs.Count);
        Assert.True(worst < 1e-8,
            $"F8 pairing incomplete: worst |λ_s + λ_f + 2Nγ| = {worst:E3}");

        foreach (var (s, f) in pairs)
        {
            double resid = Math.Abs(Light(right.Column(s)) + Light(right.Column(f)) - N);
            Assert.True(resid < 1e-8,
                $"pair ({s}, {f}): |light_s + light_f − N| = {resid:E3} " +
                $"(λ_s = {vals[s]}, λ_f = {vals[f]})");
        }
    }

    [Fact]
    public void StandingWaveNull_CrossPartnerMatrixElementsVanish()
    {
        // The April null: the J-defect bond perturbation V_L = −i[V, ·] on bond (0, 1)
        // cannot couple F8 partners. V_L conserves the joint popcount sector (a, b) of a
        // coherence while the palindrome maps (a, b) → (N−a, N−b); the sectors are
        // disjoint except at a = b = N/2, and there the matrix element still vanishes.
        // N=7 original: max 4.43e−14 over 80×80 pairs (RESULT_PTF_STANDING_WAVE_ROLES).
        var (vals, right, left) = Decomp.Value;

        // Biorthogonality guard: rows of R⁻¹ satisfy ⟨W_k|M_l⟩ = δ_kl globally.
        Assert.True((left * right - ComplexMatrix.Build.DenseIdentity(vals.Length))
            .FrobeniusNorm() < 1e-8, "left covectors not biorthogonal to right eigenvectors");

        var vL = BondPerturbation.Build(N, siteA: 0, siteB: 1, BondPerturbation.Kind.XY);
        var elems = PerturbationMatrixElements.Compute(right, left, vL);

        var (pairs, worst) = GreedyF8Pairs(vals);
        Assert.True(worst < 1e-8, $"F8 pairing incomplete: worst residual {worst:E3}");

        double maxCross = 0.0;
        foreach (var (s, f) in pairs)
        {
            maxCross = Math.Max(maxCross, elems[f, s].Magnitude);
            maxCross = Math.Max(maxCross, elems[s, f].Magnitude);
        }
        Assert.True(maxCross < 1e-9,
            $"standing-wave coupling not null: max |⟨W_f|V_L|M_s⟩| = {maxCross:E3}");
    }
}
