using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.BlockSpectrum;

/// <summary>F71AntiPalindromicGammaSpectralInvariance (Tier 1 derived: algebraic proof
/// 2026-05-12; empirical witness 2026-05-11): for chain XY + Z-dephasing Liouvillian on
/// N qubits, the <b>F71-refined diagonal-block spectrum</b> (the multiset of eigenvalues
/// of the (F71-even, F71-odd) diagonal sub-blocks of <c>Q^T L Q</c>) depends on the
/// per-site γ-distribution only through the multiset of F71-pair-sums
/// <c>S_l := γ_l + γ_{N-1-l}</c>; cross-block entries depend on pair-differences
/// <c>D_l := γ_l − γ_{N-1-l}</c>. As a corollary, the spectrum is invariant under any
/// γ-distribution satisfying <c>γ_l + γ_{N-1-l} = 2·γ_avg</c> for all l ∈ {0..N−1},
/// where <c>γ_avg = (1/N)·Σ γ_l</c>. We call such γ-distributions <i>F71-anti-palindromic
/// around their mean</i>.
///
/// <para><b>Important caveat — this is NOT full-L spectral invariance.</b> Under
/// anti-palindromic γ, the F71-rotated Liouvillian <c>Q^T L Q</c> is generally NOT
/// block-diagonal (the (F71-even ↔ F71-odd) cross blocks carry the F71-asymmetry, see
/// <see cref="InhomogeneousGammaF71BreakingWitness"/>). Diagonalising only the diagonal
/// sub-blocks (which is what <see cref="F71MirrorBlockRefinement.ComputeSpectrumPerBlock"/>
/// does) gives the same multiset across all anti-palindromic γ-profiles, but this multiset
/// is NOT the true full-L spectrum when γ breaks F71 symmetry. The full L spectrum
/// (computed via <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock"/> or direct
/// diagonalisation of <c>L</c>) differs across anti-palindromic γ-profiles by the
/// perturbation induced by the cross blocks.</para>
///
/// <para><b>What anti-palindromy controls.</b> The diagonal blocks of <c>Q^T L Q</c> in
/// the F71-refined basis depend on γ only through the <i>per-pair sums</i>
/// <c>γ_l + γ_{N-1-l}</c>. When all pair sums equal <c>2·γ_avg</c>, the diagonal blocks
/// reduce to the same form as under uniform γ = γ_avg. The off-diagonal F71-even ↔ F71-odd
/// coupling depends on per-pair <i>differences</i> <c>γ_l − γ_{N-1-l}</c> (see
/// <see cref="InhomogeneousGammaF71BreakingWitness.F71AsymmetryNorm"/>), independent of
/// the pair-sum invariance.</para>
///
/// <para>This is a <b>strict relaxation</b> of the trivial uniform-γ case (γ_l = γ_avg)
/// and a <b>distinct condition</b> from F71 symmetry (γ_l = γ_{N-1-l}, captured by
/// <see cref="F71MirrorBlockRefinement"/>): anti-palindromic and F71-palindromic
/// γ-profiles overlap only at uniform γ. Both are <b>strictly stronger</b> than F1 alone
/// (which only requires Σγ_l invariant). For odd N, the middle site
/// <c>l = (N-1)/2</c> must equal γ_avg.</para>
///
/// <para><b>Empirical witness at N=4, 5, 6.</b> The F71-refined diagonal-sub-block
/// eigenvalue multiset is bit-identical (binned at 1e-7) across:
/// <list type="bullet">
///   <item>uniform γ = 0.45,</item>
///   <item>monotonic linear anti-palindromic γ (e.g. [0.2, 0.3, 0.4, 0.5, 0.6, 0.7] at N=6),</item>
///   <item>non-monotonic anti-palindromic γ (e.g. [0.3, 0.5, 0.4, 0.5, 0.4, 0.6] at N=6).</item>
/// </list>
/// The <i>full</i> L spectrum (joint-popcount block-eig, no F71 truncation) DIFFERS
/// across these profiles, confirming that the invariance is a property of the F71-projected
/// dynamics only. Permuted non-anti-palindromic γ
/// ([0.7, 0.2, 0.5, 0.3, 0.6, 0.4]) and concentrated γ ([0.1, 0.1, 0.1, 0.1, 0.1, 2.2])
/// produce distinct F71-refined diagonal spectra.</para>
///
/// <para>Tier outcome <see cref="Tier.Tier1Derived"/> as of 2026-05-12: the algebraic
/// proof in <c>docs/proofs/PROOF_F91_GAMMA_NINETY_DEGREES.md</c> § Algebraic proof
/// (Eqs. 1–13) shows the F71-refined diagonal-block matrix elements of L = −i[H, ·] + D
/// are <b>linear functionals of γ depending only on F71-pair-sums S_l</b>: per Step 3
/// (Eqs. 7a, 7b) <c>⟨sym|D|sym⟩ = ⟨antisym|D|antisym⟩ = −Σ_{l ∈ Δ(a, b)} S_l</c>; the
/// Hamiltonian term is γ-independent and F71-block-diagonal (Step 5, Eq. 10); the
/// pair-difference content lives entirely in the F71-cross-block entries
/// <c>⟨sym|D|antisym⟩ = −Σ_{l ∈ Δ(a, b)} D_l</c> (Step 4, Eq. 9), which do NOT enter
/// diagonal-block eigenvalues. The 90°-rotation invariance (corollary, Step 7) follows
/// because the orbit S_l = 2γ_avg ∀l is closed under R_{90}: γ_l ↦ 2γ_avg − γ_{N-1-l}.
/// Empirical witness (bit-exact across multiple anti-palindromic profiles at N=4, 5, 6)
/// independently verifies the proof.</para>
///
/// <para>Anchors: <c>compute/RCPsiSquared.Core/BlockSpectrum/JointPopcountSectors.cs</c>
/// (γ-blind parent), <c>compute/RCPsiSquared.Core/BlockSpectrum/F71MirrorBlockRefinement.cs</c>
/// (γ-symmetric parent — anti-palindromy is a distinct condition that overlaps only at uniform γ),
/// <c>compute/RCPsiSquared.Core.Tests/BlockSpectrum/F71AntiPalindromicGammaSpectralInvarianceTests.cs</c>
/// (empirical witness at N=4, 5, 6 + explicit verification that full-L spectrum DIFFERS).</para></summary>
public sealed class F71AntiPalindromicGammaSpectralInvariance : Claim
{
    private readonly JointPopcountSectors _sectors;
    private readonly F71MirrorBlockRefinement _f71;

    public F71AntiPalindromicGammaSpectralInvariance(
        JointPopcountSectors sectors,
        F71MirrorBlockRefinement f71)
        : base("F71AntiPalindromicGammaSpectralInvariance: chain XY+Z-deph F71-refined diagonal-block matrix elements depend on γ only through F71-pair-sums S_l = γ_l + γ_{N-1-l}; cross-block entries depend on pair-differences D_l = γ_l - γ_{N-1-l}; corollary: spectrum invariant under any γ-distribution satisfying γ_l + γ_{N-1-l} = 2·γ_avg (the orbit of the 90°-rotation R_{90}: γ_l ↦ 2γ_avg − γ_{N-1-l}); full L spectrum differs across the orbit (cross blocks carry the asymmetry); algebraic proof complete (PROOF_F91 § Algebraic proof, Eqs. 1–13) + bit-exact verified at N=4,5,6.",
               Tier.Tier1Derived,
               "JointPopcountSectors (γ-blind parent) + F71MirrorBlockRefinement (γ-symmetric parent; anti-palindromy distinct, both reduce to uniform-γ at intersection); algebraic proof in PROOF_F91 § Algebraic proof shows diagonal-block matrix elements are linear functionals of pair-sums S_l (Eqs. 7a, 7b, 11a, 11b), Hamiltonian γ-independent and F71-block-diagonal (Eq. 10), pair-differences confined to off-diagonal cross-block entries (Eq. 9); empirical bit-exact witness at N=4,5,6 in F71AntiPalindromicGammaSpectralInvarianceTests")
    {
        _sectors = sectors ?? throw new ArgumentNullException(nameof(sectors));
        _f71 = f71 ?? throw new ArgumentNullException(nameof(f71));
    }

    /// <summary>F71-anti-palindromic deviation norm of a per-site γ distribution:
    /// <c>sqrt(Σ_l ((γ_l + γ_{N-1-l}) − 2·γ_avg)²)</c> where
    /// <c>γ_avg = (1/N)·Σ_l γ_l</c>. Equals 0 iff every F71-pair sum equals
    /// <c>2·γ_avg</c> (γ-distribution F71-anti-palindromic around its mean),
    /// otherwise positive. The empirical scaling proxy for spectral deviation of the
    /// F71-refined diagonal sub-blocks from the uniform-γ_avg reference.</summary>
    public static double AntiPalindromicDeviation(IReadOnlyList<double> gammas)
    {
        if (gammas is null) throw new ArgumentNullException(nameof(gammas));
        int N = gammas.Count;
        if (N == 0) return 0.0;
        double avg = 0.0;
        for (int l = 0; l < N; l++) avg += gammas[l];
        avg /= N;
        double sumSq = 0.0;
        for (int l = 0; l < N; l++)
        {
            double pair = gammas[l] + gammas[N - 1 - l];
            double diff = pair - 2.0 * avg;
            sumSq += diff * diff;
        }
        return Math.Sqrt(sumSq);
    }

    /// <summary>True iff γ is F71-anti-palindromic around its mean within
    /// <paramref name="tolerance"/>. For odd N this implicitly requires
    /// <c>γ_{(N-1)/2} = γ_avg</c> (middle site self-paired).</summary>
    public static bool IsAntiPalindromic(IReadOnlyList<double> gammas, double tolerance = 1e-10)
        => AntiPalindromicDeviation(gammas) < tolerance;

    public override string DisplayName =>
        "F71AntiPalindromicGammaSpectralInvariance: F71-refined diagonal-block spectrum invariant under γ_l + γ_{N-1-l} = 2·γ_avg";

    public override string Summary =>
        $"chain XY+Z-deph F71-refined diagonal-block matrix elements depend on γ only through F71-pair-sums S_l; corollary: spectrum invariant on the anti-palindromic orbit (S_l = 2·γ_avg ∀l, closed under 90°-rotation); full L spectrum differs (cross blocks carry pair-differences D_l); algebraic proof + bit-exact N=4,5,6 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("parent-1",
                summary: "JointPopcountSectors (γ-blind: U(1)×U(1) per-side popcount conservation, holds for any γ_l)");
            yield return new InspectableNode("parent-2",
                summary: "F71MirrorBlockRefinement (γ-symmetric: exact iff γ_l = γ_{N-1-l}; distinct from anti-palindromy, overlap only at uniform γ)");
            yield return new InspectableNode("condition",
                summary: "γ_l + γ_{N-1-l} = 2·γ_avg for all l (distinct from F71 symmetry; strictly stronger than F1 alone)");
            yield return new InspectableNode("scope",
                summary: "applies to F71-refined DIAGONAL sub-block spectrum only; full L spectrum differs because (F71-even ↔ F71-odd) cross blocks are nonzero under F71-asymmetric γ");
            yield return new InspectableNode("prediction-1",
                summary: "anti-palindromic γ → F71-refined diagonal sub-block spectrum = uniform-γ_avg sub-block spectrum (bit-identical)");
            yield return new InspectableNode("prediction-2",
                summary: "non-anti-palindromic γ with same Σγ → F71-refined diagonal sub-block spectrum DIFFERS (F1-only protection insufficient)");
            yield return new InspectableNode("prediction-3",
                summary: "odd N: middle site γ_{(N-1)/2} must equal γ_avg (forced by self-pairing)");
            yield return new InspectableNode("prediction-4",
                summary: "full L spectrum DIFFERS across anti-palindromic γ-profiles (cross-block perturbation lifts the diagonal-block degeneracy)");
            yield return new InspectableNode("witness",
                summary: "bit-exact F71-refined diagonal-block spectrum equality across anti-palindromic γ-profiles at N=4,5,6 chain XY+Z-deph");
            yield return new InspectableNode("proof",
                summary: "PROOF_F91 § Algebraic proof (Eqs. 1–13, 2026-05-12): diagonal-block matrix elements ⟨sym|D|sym⟩ = ⟨antisym|D|antisym⟩ = −Σ S_l (linear in pair-sums only); cross-block ⟨sym|D|antisym⟩ = −Σ D_l (pair-differences only); H γ-independent and F71-block-diagonal; 90°-rotation invariance is corollary on the orbit S_l = 2γ_avg ∀l");
        }
    }
}
