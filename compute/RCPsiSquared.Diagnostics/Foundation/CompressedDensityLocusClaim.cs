using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>F154's conditional compression law and its independent (1,1) interval sibling.
/// The live N=11 physical counterexample and the zero-frequency endpoints are computed exactly
/// by <see cref="CompressedDensityN11Witness"/>.
///
/// <para>Five typed parents. <see cref="AbsorptionTheoremClaim"/> supplies the per-site
/// disagreement split D = −2Σ_l γ_l N_l, and <see cref="F71AntiPalindromicGammaSpectralInvariance"/>
/// the mirror-balanced locus. The (1,1) half stands on three more. <see cref="JointPopcountSectors"/>
/// makes the (1,1) block a block at all. <see cref="SeedRungGramClaim"/> (F143) is the typed home of the
/// zero-frequency room's spectrum: there the compressed dissipator is −4γ̄(I − WᵀW), and on the uniform
/// XY chain WᵀW is the Gram of PROOF_R90_FROZEN_DIVISOR Lemma 5 at M = N + 1, re-derived as F143.
/// <see cref="FrozenDivisorClaim"/> (F140) owns the lower endpoint −4γ̄ as an exact (1,1)
/// Liouvillian eigenvalue at every coupling J. The interval theorem is the part that is
/// independent of all three; the endpoints are theirs.</para></summary>
public sealed class CompressedDensityLocusClaim : Claim
{
    public AbsorptionTheoremClaim Absorption { get; }
    public F71AntiPalindromicGammaSpectralInvariance Locus { get; }

    /// <summary>Parent edge: the (1,1) block is a block of the joint-popcount grading.</summary>
    public JointPopcountSectors Sectors { get; }

    /// <summary>Parent edge: F143's Gram G (Lemma 5's on XY) gives the zero-frequency room's whole
    /// compressed spectrum on the uniform XY chain, −4γ̄(I − G), on every locus profile.</summary>
    public SeedRungGramClaim SeedRungGram { get; }

    /// <summary>Parent edge: F140's frozen root −4γ̄, the lower endpoint at every J for real h.</summary>
    public FrozenDivisorClaim FrozenDivisor { get; }

    public CompressedDensityLocusClaim(
        AbsorptionTheoremClaim absorption,
        F71AntiPalindromicGammaSpectralInvariance locus,
        JointPopcountSectors sectors,
        SeedRungGramClaim seedRungGram,
        FrozenDivisorClaim frozenDivisor)
        : base(
            "F154 compressed density on the R90 locus: for a complete ad_H eigenspace Ω, " +
            "mirror balance plus C_l = PΩ(N_l-N_mirror(l))PΩ = 0 for each site gives " +
            "DΩ = -2γbar PΩ N_XY PΩ and the corresponding size-class interval. " +
            "The independent (1,1) interval theorem holds for a Hermitian, simple, reflection-symmetric " +
            "one-excitation h and nonnegative mirror-balanced rates even when C_l is nonzero. " +
            "At N=11 XY, C_0 is nonzero and the F154 identity fails on a balanced profile; " +
            "the (1,1) interval survives. The block's two endpoints sit in the zero-frequency room, " +
            "a scalar-parity room where the F154 identity is a theorem: its compression is " +
            "-4γbar(I - WᵀW), W_lk = |ψ_k(l)|², reached by I and by ker W for every such h " +
            "(on the uniform XY chain WᵀW is the Gram of PROOF_R90_FROZEN_DIVISOR Lemma 5, " +
            "F143's G); for real h F140 makes -4γbar an exact (1,1) eigenvalue at every J.",
            Tier.Tier1Derived,
            "docs/ANALYTICAL_FORMULAS.md F154 + docs/proofs/PROOF_MIXED_SPACE_REFLECTION_LAW.md " +
            "+ docs/proofs/PROOF_N11_COMPRESSED_DENSITY.md + " +
            "simulations/n11_compressed_density_gate.py + " +
            "compute/RCPsiSquared.Diagnostics/Foundation/CompressedDensityN11Witness.cs " +
            "(inspect --root compresseddensity)")
    {
        Absorption = absorption ?? throw new ArgumentNullException(nameof(absorption));
        Locus = locus ?? throw new ArgumentNullException(nameof(locus));
        Sectors = sectors ?? throw new ArgumentNullException(nameof(sectors));
        SeedRungGram = seedRungGram ?? throw new ArgumentNullException(nameof(seedRungGram));
        FrozenDivisor = frozenDivisor ?? throw new ArgumentNullException(nameof(frozenDivisor));
    }

    public string ConditionalIdentity =>
        "For the complete ad_H eigenspace Ω and real mirror-balanced γ, " +
        "C_l = 0 for every l implies " +
        "DΩ = −2γbar PΩ N_XY PΩ. The physical N_l counts ket/bra disagreement at site l. " +
        "The identity is algebraic and needs no sign on the rates. " +
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
        "spec(DΩ) ⊂ [−4γbar, 0]. This does not require C_l = 0, and it is the part of this " +
        "claim that no parent supplies.";

    public string SignedAgreementContrast =>
        "In every complete (1,1) frequency space of a Hermitian, simple, reflection-symmetric " +
        "one-excitation h, the one-legged reflected density difference compresses to zero. " +
        "For m=N−1−l, C_l = −2PΩ(|ll><ll|−|mm><mm|)PΩ. This signed agreement or " +
        "double-occupancy overlap can connect opposite reflection parities; its off-diagonal " +
        "entries are coherent amplitudes, not probabilities. At uniform rate and zero frequency " +
        "the same agreement projector is F143's Gram.";

    public string BlockWideEndpointWitnesses =>
        "For a Hermitian, simple, reflection-symmetric one-excitation h the zero-frequency room " +
        "is Ω0 = span{P_k}, P_k = |ψ_k><ψ_k|. Every P_k is reflection-even, so Ω0 is a " +
        "scalar-parity room and C_l = 0 there by F154's parity theorem. With W_lk = |ψ_k(l)|², " +
        "whose mirror rows agree, the compression is DΩ0 = −4γbar(I − G), G = WᵀW, on every " +
        "real mirror-balanced profile (profile-blind). I = Σ_k P_k is the one-excitation " +
        "identity, physical-cell diagonal, with D I = 0: the upper endpoint. Since rank W ≤ " +
        "⌈N/2⌉, ker G has dimension at least ⌊N/2⌋ and sits at −4γbar: the lower endpoint. " +
        "On the uniform open XY chain G is Lemma 5's Gram of PROOF_R90_FROZEN_DIVISOR, F143's " +
        "(1/M)(11ᵀ + (I+R_χ)/2), M = N+1, R_χ the chiral mode flip k ↦ M−k, " +
        "so spec DΩ0 = " +
        "{−4γbar ×⌊N/2⌋, −4γbar(1−1/M) ×(⌈N/2⌉−1), 0 ×1}, and on XY a chiral difference " +
        "Q = P_k−P_(N+1−k) spans one kernel direction. Both endpoint values are also exact " +
        "(1,1) Liouvillian eigenvalues at every J: 0 through the sector projector I, which the " +
        "generator annihilates, and −4γbar with multiplicity ≥ ⌊N/2⌋ through F140 for real h. " +
        "Containment in [−4γbar, 0] is a strong-coupling compression statement; at small J the " +
        "(1,1) block can have rates below −4γbar.";

    public string PerRoomLowerEndpoint =>
        "On the uniform open XY chain the chiral transpose (a,b) ↦ (M−b, M−a) preserves a dyad's " +
        "frequency and its contact amplitudes ψ_a(z)ψ_b(z), since ψ_(M−k)(z) = (−1)^z ψ_k(z). " +
        "Each non-fixed orbit {x, τx} gives e_x − e_τx in the kernel of the room's T, so in every " +
        "room mult(−4γbar) ≥ the number of non-fixed orbits, for every balanced profile. Equality " +
        "is exact in all 55 N=11 rooms at the profile (2,1,...,1,0) (CompressedDensityN11Witness) " +
        "and measured N = 2..14 (simulations/n11_compressed_density_census.py, section A).";

    public string N11Counterexample =>
        "N=11 uniform XY: Ω at frequency 4J cos(π/12) contains (1,6),(3,7),(5,9),(6,11). " +
        "For X=|ψ1><ψ6| and Y=|ψ5><ψ9|, " +
        "<Y,(N_0−N_10)X>=−sqrt(2)/72. With γ=(2,1,1,1,1,1,1,1,1,1,0), " +
        "<Y,DΩ X>=sqrt(2)/36 but the F154 conditional RHS cross entry is zero. " +
        "The profile remains mirror-balanced; it violates the C_l=0 hypothesis. " +
        "On the uniform (1,1) block every opposite-parity pair x=(a,b), y=(c,d) in one room has " +
        "<y,C_0 x> = −4ψ_cψ_aψ_dψ_b(0) ≠ 0, so every parity-mixed (1,1) room fires; such " +
        "rooms exist, measured on N ≤ 30 (census section B), at N = 11, 14, 17, 20, 23 and 29 " +
        "and nowhere else.";

    public string Scope =>
        "The F154 identity is conditional, not an all-N mixed-space rule. " +
        "The (1,1) interval is the strong-coupling compression under a Hermitian, simple, " +
        "reflection-symmetric h and nonnegative mirror-balanced rates; complete Ω is required. " +
        "The block-wide endpoints hold for every such h; the upper endpoint is absent from the " +
        "N=11 nonzero-frequency rooms at the tested profiles. " +
        "Open: parity doors of the higher joint-popcount blocks, the per-room upper endpoint at " +
        "profiles with dark pairs, and the N=6 M0 arithmetic of the transversal certificate. " +
        "No hardware claim follows. The C# witness recomputes the N=11 readings in exact " +
        "Q(√2,√3) arithmetic, as the Python gate does in SymPy.";

    public override string DisplayName => "F154: compressed density on the mirror-balanced locus";
    public override string Summary =>
        "Conditional size-class identity; N=11 physical mixed-parity obstruction; " +
        "independent (1,1) interval; block endpoints in the scalar zero-frequency room " +
        "(Lemma 5 / F143, F140). Live N=11 root: compresseddensity.";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return Absorption;
            yield return Locus;
            yield return Sectors;
            yield return SeedRungGram;
            yield return FrozenDivisor;
            yield return new InspectableNode("conditional F154 identity", ConditionalIdentity);
            yield return new InspectableNode("positive-rate pure-class endpoint criterion", PureClassEndpointCriterion);
            yield return new InspectableNode("independent (1,1) interval theorem", OneExcitationInterval);
            yield return new InspectableNode("signed agreement contrast", SignedAgreementContrast);
            yield return new InspectableNode("block-wide endpoints in the zero-frequency room", BlockWideEndpointWitnesses);
            yield return new InspectableNode("per-room lower endpoint", PerRoomLowerEndpoint);
            yield return new InspectableNode("physical N=11 obstruction", N11Counterexample);
            yield return new InspectableNode("scope", Scope);
        }
    }
}
