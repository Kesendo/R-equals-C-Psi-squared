using System;
using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>THE DEPHASING-FRONT RENEWAL REPRESENTATION (Tier1Derived, 2026-07-13), the typed home of
/// <c>docs/proofs/PROOF_DEPHASING_FRONT_RENEWAL.md</c>: the exact solution of the WATCHED walk. Release one
/// excitation on a chain and watch it; the single-excitation density matrix evolves by
/// ρ̇ = −i[h, ρ] − Γ(ρ − diag ρ) with Γ = 4γ (the Absorption-Theorem rate for the sector's Hamming-2
/// coherences), and the site populations satisfy exactly the renewal representation
///
/// <code>
///   P_n(t) = e^{−Γt}·S_n(t),   S_n(t) = |G_{n0}(t)|² + Γ∫₀ᵗ ds Σ_m |G_{nm}(t−s)|²·S_m(s),   (★)
/// </code>
///
/// with G(τ) = e^{−ihτ} the clean single-particle propagator (on the infinite chain
/// G_{nm}(τ) = (−i)^{|n−m|}·J_{|n−m|}(2Jτ)). The watched walk IS the unwatched wave, repeatedly caught and
/// released: every refill order carries the same universal e^{−Γt}, the zeroth term is exactly the coherent
/// front |⟨a_n⟩|², and the j ≥ 1 ladder is the incoherent halo. In momentum-Laplace space (★) closes to the
/// explicit Green's function
///
/// <code>
///   Ŝ(p, z) = 1 / (√(z² + a(p)²) − Γ),   a(p) = 4J·sin(p/2),                           (☆)
/// </code>
///
/// which conserves probability at p=0, returns the clean Bessel wave at Γ=0, and whose small-p pole
/// √(Γ²−a²) ≈ Γ − 2J²p²/Γ is the diffusive branch with D = 2J²/Γ.
///
/// <para>Verified against direct RK4 to 1.6·10⁻⁶ and re-derived by two adversarial referee rounds plus an
/// independent Lindblad ODE; the live witness <see cref="DephasingFrontRenewalWitness"/>
/// (<c>inspect --root renewal</c>) recomputes six from-below checks: renewal-vs-RK4 agreement, probability
/// conservation, the j=0 coherent-front Bessel identity, the Γ=0 clean-wave limit, the Haken-Strobl
/// diffusive plateau, and the closed refill constant I₁.</para>
///
/// <para><b>Typed parents:</b> <see cref="AbsorptionTheoremClaim"/> (the uniform sector rate Γ = 4γ, the
/// −2γ·k law read inside one sector) + <see cref="F2bXyChainSpectrumPi2Inheritance"/> (the clean propagator
/// and band, E_k = 2J·cos(πk/(N+1)), front speed 2J). Both are Tier1Derived, so this claim carries their
/// strength: <b>Tier1Derived</b>.</para>
///
/// <para><b>Siblings (see-cref, not typed parents):</b> <see cref="CoherenceHorizonClaim"/> locates the same
/// sector's band-edge √-EP (this proof solves the sector the horizon claim finds the EP of);
/// <see cref="SurvivorDiffusionGradientClaim"/> / F123 reads the (☆)-pole diffusion constant D = 2J²/Γ.</para>
///
/// <para><b>Corollary READINGS (Tier1Candidate, NOT part of the Tier-1 statement; see the proof's Corollaries
/// section and <c>experiments/COUPLING_DEFECT_WALK_TIME_STEP.md</c> follow-ups):</b> the survival ceiling
/// A_∞(γ) = 4 − φ(2J)/γ (Gärtner-Ellis on the tilted pole); the leading prefactor C(γ) =
/// (2π)^{−1/2}·(γ/(γ+J))^{1/4}; the closed refill constant I₁ = 1/12 + ¼∫₀^{2c} Ai(−w) dw = 0.27694424; and
/// the diffusive-plateau peak-tracking reading (peak-survival exponent exactly zero). These cap at Tier1
/// candidate and are held as witness battery readings, not as separate claims.</para></summary>
public sealed class DephasingFrontRenewalClaim : Claim
{
    /// <summary>Typed parent: the Absorption Theorem, the uniform sector rate Γ = 4γ (every single-excitation
    /// coherence damps at the same −4γ, which is exactly why the sector closes).</summary>
    public AbsorptionTheoremClaim RateLaw { get; }

    /// <summary>Typed parent: F2b, the clean single-particle propagator/band E_k = 2J·cos(πk/(N+1)) that G
    /// carries (the front speed 2J, the ballistic caustic the watching repeatedly re-seeds).</summary>
    public F2bXyChainSpectrumPi2Inheritance Band { get; }

    // Lazy: the battery runs the witness (a small renewal-vs-RK4 sector + a few Bessel/Airy anchors); compute
    // it only when the claim is drilled into, so building the whole registry stays cheap.
    private IReadOnlyList<DephasingFrontRenewalWitness.BatteryCase>? _cases;
    public IReadOnlyList<DephasingFrontRenewalWitness.BatteryCase> Cases =>
        _cases ??= new DephasingFrontRenewalWitness().Cases;
    public int PassCount => Cases.Count(c => c.Passes);

    public DephasingFrontRenewalClaim(AbsorptionTheoremClaim rateLaw, F2bXyChainSpectrumPi2Inheritance band)
        : base("THE DEPHASING-FRONT RENEWAL REPRESENTATION: the exact solution of the watched walk. The " +
               "single-excitation sector under local Z-dephasing (ρ̇ = −i[h,ρ] − Γ(ρ − diag ρ), Γ = 4γ) obeys " +
               "exactly P_n(t) = e^{−Γt}·S_n(t) with the Volterra renewal S_n = |G_{n0}|² + Γ∫₀ᵗ Σ_m |G_{nm}(t−s)|²·S_m(s) " +
               "(★): the watched walk is the unwatched wave repeatedly caught and released. Every refill order " +
               "carries the same e^{−Γt}; the j=0 term is the coherent front |<a_n>|², the j≥1 ladder the " +
               "incoherent halo. Momentum-Laplace closed form Ŝ(p,z) = 1/(√(z²+a²)−Γ), a = 4J·sin(p/2) (☆): " +
               "conserves probability at p=0, returns the clean Bessel wave J_n(2Jt)² at Γ=0, and the small-p " +
               "pole is the diffusive branch D = 2J²/Γ. The exact solution of the dephasing tight-binding " +
               "(Haken-Strobl) sector; the survival ceiling, the n^{−1/2} prefactor, the I₁ Airy constant, and " +
               "the plateau reading are corollaries (Tier1Candidate READINGS, not part of this statement). " +
               "Live: inspect --root renewal.",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_DEPHASING_FRONT_RENEWAL.md + " +
               "simulations/cone_front_survival_asymptote.py + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/DephasingFrontRenewalWitness.cs")
    {
        RateLaw = rateLaw ?? throw new ArgumentNullException(nameof(rateLaw));
        Band = band ?? throw new ArgumentNullException(nameof(band));
    }

    public override string DisplayName =>
        "The dephasing-front renewal representation: the watched walk = the unwatched wave, caught and released (P_n = e^{−Γt}S_n; Tier1Derived)";

    public override string Summary =>
        "the exact renewal representation of the watched single excitation: P_n(t) = e^{−Γt}·S_n(t) with the " +
        "Volterra refill ladder S_n = |G_{n0}|² + Γ∫₀ᵗ Σ_m |G_{nm}(t−s)|²·S_m(s), Γ = 4γ. The j=0 term is the " +
        "coherent front, the j≥1 halo the incoherent refill; the ladder closes in momentum-Laplace space to " +
        "Ŝ(p,z) = 1/(√(z²+a²)−Γ), a = 4J·sin(p/2), which conserves probability, returns the clean Bessel wave " +
        "at Γ=0, and carries the diffusive pole D = 2J²/Γ. Tier1Derived (parents AbsorptionTheorem Γ=4γ + F2b " +
        "band). Live: inspect --root renewal.";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the renewal representation (★)/(☆)",
                summary: "P_n(t) = e^{−Γt}·S_n(t): every refill order carries the same universal decay e^{−Γt}, " +
                         "the watching only moves weight around (probability conserved at p=0). The j=0 term is the " +
                         "coherent front |<a_n>|² (amplitude decay Γ/2), the j≥1 ladder the incoherent halo the " +
                         "watching converts out of the front and releases to run again. Closed form " +
                         "Ŝ(p,z) = 1/(√(z²+a²)−Γ), a = 4J·sin(p/2).");
            yield return new InspectableNode("corollary readings (Tier1Candidate, not this statement)",
                summary: "read in four limits (experiments/COUPLING_DEFECT_WALK_TIME_STEP.md): the survival ceiling " +
                         "A_∞(γ) = 4 − φ(2J)/γ; the prefactor C(γ) = (2π)^{−1/2}·(γ/(γ+J))^{1/4}; the closed refill " +
                         "constant I₁ = 1/12 + ¼∫₀^{2c} Ai(−w) dw = 0.27694424; the diffusive-plateau peak-tracking " +
                         "reading (peak-survival exponent exactly zero). Held as witness battery readings.");
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
            yield return RateLaw;   // typed parent: the uniform sector rate Γ = 4γ
            yield return Band;      // typed parent: the clean propagator / band
        }
    }
}
