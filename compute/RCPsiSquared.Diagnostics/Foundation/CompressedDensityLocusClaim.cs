using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>F154's conditional compression law and its distinct (1,1) interval sibling.
/// The live N=11 physical counterexample is computed by <see cref="CompressedDensityN11Witness"/>.
/// Endpoint attainment outside the stated pure-class census is not promoted to a theorem.</summary>
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
            "the (1,1) interval survives. These statements concern strong-coupling compressions, " +
            "not finite-J Liouvillian eigenvalues or unconditional endpoint attainment.",
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
        "For the complete ad_H eigenspace Ω and mirror-balanced γ, C_l = 0 for every l implies " +
        "DΩ = −2γbar PΩ N_XY PΩ. The physical N_l counts ket/bra disagreement at site l. " +
        "Rayleigh bounds the compressed spectrum by the supported size classes. " +
        "Endpoint attainment requires pure-class vectors in the appropriate eigenspaces; " +
        "the F154 census establishes this only on its enumerated rows.";

    public string OneExcitationInterval =>
        "For every complete (1,1) frequency space Ω of a Hermitian, simple, reflection-symmetric " +
        "one-excitation h, with γ_l ≥ 0 and γ_l+γ_mirror(l)=2γbar, " +
        "DΩ = −4γbar IΩ + 4 PΩ T PΩ, T = Σ_l γ_l |ll><ll|. " +
        "The compressed one-body rate is scalar because a shared leg in an equal-frequency " +
        "dyad pair forces the other leg to agree, and each simple reflection eigenmode has " +
        "mean rate γbar. Since T ≥ 0 and the physical D ≤ 0, " +
        "spec(DΩ) ⊂ [−4γbar, 0]. This does not require C_l = 0 or assert endpoints.";

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
        "No claim about higher joint-popcount blocks, endpoint attainment in N=11, " +
        "finite J, hardware, or all-N identity restoration follows. The C# witness uses " +
        "double arithmetic with an explicit error budget; the Python gate checks the N=11 " +
        "radicals exactly.";

    public override string DisplayName => "F154: compressed density on the mirror-balanced locus";
    public override string Summary =>
        "Conditional size-class identity; N=11 physical mixed-parity obstruction; " +
        "separate (1,1) interval theorem. Live N=11 root: compresseddensity.";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return Absorption;
            yield return Locus;
            yield return new InspectableNode("conditional F154 identity", ConditionalIdentity);
            yield return new InspectableNode("independent (1,1) interval theorem", OneExcitationInterval);
            yield return new InspectableNode("physical N=11 obstruction", N11Counterexample);
            yield return new InspectableNode("scope", Scope);
        }
    }
}
