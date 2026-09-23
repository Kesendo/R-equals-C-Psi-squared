using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>F154's conditional compression law and its distinct (1,1) interval sibling.
/// The live N=11 physical counterexample is computed by <see cref="CompressedDensityN11Witness"/>.
/// Uniform open XY has separate block-wide (1,1) endpoint witnesses at zero frequency;
/// the F154 pure-class census has its own stated scope.</summary>
public sealed class CompressedDensityLocusClaim : Claim
{
    public AbsorptionTheoremClaim Absorption { get; }
    public F71AntiPalindromicGammaSpectralInvariance Locus { get; }

    public CompressedDensityLocusClaim(
        AbsorptionTheoremClaim absorption,
        F71AntiPalindromicGammaSpectralInvariance locus)
        : base(
            "F154 compressed density on the R90 locus: for a complete ad_H eigenspace Ω, " +
            "mirror balance plus C_l = PΩ(N_l-N_mirror(l))PΩ = 0 for each site gives " +
            "DΩ = -2γbar PΩ N_XY PΩ and the corresponding size-class interval. " +
            "The independent (1,1) interval theorem holds for a Hermitian, simple, reflection-symmetric " +
            "one-excitation h and nonnegative mirror-balanced rates even when C_l is nonzero. " +
            "At N=11 XY, C_0 is nonzero and the F154 identity fails on a balanced profile; " +
            "the (1,1) interval survives. On the uniform open XY chain, zero-frequency " +
            "physical-cell witnesses I and Q reach both block-wide (1,1) endpoints. " +
            "These statements concern strong-coupling compressions, not finite-J Liouvillian eigenvalues.",
            Tier.Tier1Derived,
            "docs/ANALYTICAL_FORMULAS.md F154 + docs/proofs/PROOF_MIXED_SPACE_REFLECTION_LAW.md " +
            "+ docs/proofs/PROOF_N11_COMPRESSED_DENSITY.md + " +
            "simulations/n11_compressed_density_gate.py + " +
            "compute/RCPsiSquared.Diagnostics/Foundation/CompressedDensityN11Witness.cs " +
            "(inspect --root compresseddensity)")
    {
        Absorption = absorption ?? throw new ArgumentNullException(nameof(absorption));
        Locus = locus ?? throw new ArgumentNullException(nameof(locus));
    }

    public string ConditionalIdentity =>
        "For the complete ad_H eigenspace Ω and nonnegative mirror-balanced γ, " +
        "C_l = 0 for every l implies " +
        "DΩ = −2γbar PΩ N_XY PΩ. The physical N_l counts ket/bra disagreement at site l. " +
        "Rayleigh bounds the compressed spectrum by the supported size classes. " +
        "Pure-class vectors reach their class centres; the converse for block endpoints " +
        "has the positivity and nondegenerate-space premises stated separately.";

    public string PureClassEndpointCriterion =>
        "On the F154 rows where every relevant Ω obeys C_l=0, at γbar > 0 a multi-class " +
        "block reaches both size-class-centre endpoints iff both extreme classes have " +
        "pure-class vectors in colliding eigenspaces, provided the gated census confirms " +
        "that no nondegenerate eigenspace has a pure-class vector. At γbar = 0, physical " +
        "nonnegative mirror-balanced rates all vanish, the interval is {0}, and endpoint " +
        "attainment is automatic; the pure-vector iff does not apply.";

    public string OneExcitationInterval =>
        "For every complete (1,1) frequency space Ω of a Hermitian, simple, reflection-symmetric " +
        "one-excitation h, with γ_l ≥ 0 and γ_l+γ_mirror(l)=2γbar, " +
        "DΩ = −4γbar IΩ + 4 PΩ T PΩ, T = Σ_l γ_l |ll><ll|. " +
        "The compressed one-body rate is scalar because a shared leg in an equal-frequency " +
        "dyad pair forces the other leg to agree, and each simple reflection eigenmode has " +
        "mean rate γbar. Since T ≥ 0 and the physical D ≤ 0, " +
        "spec(DΩ) ⊂ [−4γbar, 0]. This does not require C_l = 0. " +
        "The uniform XY chain's block-wide endpoint witnesses are stated separately.";

    public string SignedAgreementContrast =>
        "In every complete (1,1) frequency space of a Hermitian, simple, reflection-symmetric " +
        "one-excitation h, the one-legged reflected density difference compresses to zero. " +
        "For m=N−1−l, C_l = −2PΩ(|ll><ll|−|mm><mm|)PΩ. This signed agreement or " +
        "double-occupancy overlap can connect opposite reflection parities; its off-diagonal " +
        "entries are coherent amplitudes, not probabilities. No F143 Gram closed form is assumed.";

    public string UniformXyEndpointWitnesses =>
        "On the uniform open XY chain at every N≥2 and every nonnegative mirror-balanced rate " +
        "profile, the zero-frequency room contains I = Σ_k P_k, the one-excitation identity, " +
        "which is physical-cell diagonal and obeys D I=0. For a distinct chiral pair, " +
        "Q = P_k−P_(N+1−k) is nonzero, has zero physical diagonal, and obeys " +
        "DΩ0 Q=−4γbar Q. Thus the union of complete (1,1) compressed spectra reaches " +
        "both block-wide endpoints [−4γbar,0], also at γbar=0 when they coincide. " +
        "This does not assert both endpoints in every frequency room.";

    public string N11Counterexample =>
        "N=11 uniform XY: Ω at frequency 4J cos(π/12) contains (1,6),(3,7),(5,9),(6,11). " +
        "For X=|ψ1><ψ6| and Y=|ψ5><ψ9|, " +
        "<Y,(N_0−N_10)X>=−sqrt(2)/72. With γ=(2,1,1,1,1,1,1,1,1,1,0), " +
        "<Y,DΩ X>=sqrt(2)/36 but the F154 conditional RHS cross entry is zero. " +
        "The profile remains mirror-balanced; it violates the C_l=0 hypothesis.";

    public string Scope =>
        "The original F154 identity is conditional, not an all-N mixed-space rule. " +
        "The new interval is only the (1,1) strong-coupling compression under a Hermitian, simple " +
        "reflection-symmetric h and nonnegative mirror-balanced rates; complete Ω is required. " +
        "Block-wide (1,1) endpoint attainment uses the uniform open XY chiral pair; " +
        "it does not extend to each N=11 frequency room or general reflection-symmetric h. " +
        "No claim about higher joint-popcount blocks, finite J, hardware, or all-N " +
        "identity restoration follows. The C# witness uses " +
        "double arithmetic with an explicit error budget; the Python gate checks the N=11 " +
        "radicals exactly.";

    public override string DisplayName => "F154: compressed density on the mirror-balanced locus";
    public override string Summary =>
        "Conditional size-class identity; N=11 physical mixed-parity obstruction; " +
        "separate (1,1) interval and uniform-XY block endpoints. Live N=11 root: compresseddensity.";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return Absorption;
            yield return Locus;
            yield return new InspectableNode("conditional F154 identity", ConditionalIdentity);
            yield return new InspectableNode("positive-rate pure-class endpoint criterion", PureClassEndpointCriterion);
            yield return new InspectableNode("independent (1,1) interval theorem", OneExcitationInterval);
            yield return new InspectableNode("signed agreement contrast", SignedAgreementContrast);
            yield return new InspectableNode("uniform XY block-wide endpoint witnesses", UniformXyEndpointWitnesses);
            yield return new InspectableNode("physical N=11 obstruction", N11Counterexample);
            yield return new InspectableNode("scope", Scope);
        }
    }
}
