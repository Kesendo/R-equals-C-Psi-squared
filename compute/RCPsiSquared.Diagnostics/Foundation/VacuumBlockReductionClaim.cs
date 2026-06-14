using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The birth-canal boundary's slowest-rate mode is the ODD, number-changing
/// |1-excitation&gt;&lt;vacuum| (0,1) coherence (Tier1Derived; the load-bearing structure is
/// derived, the global-slowest reading is bit-exact verified). The live witness is
/// <see cref="SectorReductionWitness"/> (<c>inspect --root reduction</c>); this claim banks the
/// settled N=5 result and the analytic flat-γ blindness.
///
/// <para><b>What is DERIVED.</b> The (0,1) Liouville sector is an exact invariant sub-block of L:
/// the ket-excitation-number and bra-excitation-number (the ket#/bra# bi-grading) are conserved by
/// both H_unit (XY hopping preserves total excitation in each factor) and Z-dephasing (diagonal in
/// the coherence basis), so L is block-diagonal across the joint-popcount sectors and the (PCol=0,
/// PRow=1) block, an N-dimensional sub-block, never mixes with the rest. Its reduced generator is
/// <c>L_(1,0) = −iQ·h − 2·diag(γ)</c> with h the single-particle hopping matrix. <b>Flat-γ
/// blindness is analytic at every N</b>: at uniform γ, <c>L_(1,0) = −iQ·h − 2γ·I</c>, and −iQ·h is
/// anti-Hermitian, so Re(λ) = −2γ for every mode of the block, rate Q-invariant by uniformity
/// alone (the Absorption Theorem read inside this sector).</para>
///
/// <para><b>What is VERIFIED bit-exact.</b> That this N-dim (0,1) block carries the GLOBAL slowest
/// mode across the whole γ-profile surface at <b>N=5</b>: the reduced sub-spectrum reproduces the
/// full-4^N boundary mode of <see cref="PostEpFlowField"/> (the live <see cref="SectorReductionWitness"/>
/// pins it bit-identical at the canal profile), and the standalone scan
/// <c>simulations/birth_canal_vacuum_block_verifier.py</c> (80 surface points) finds worst gap
/// 5.9·10⁻¹² between the (0,1) sub-spectrum's slowest non-kernel rate and the full Liouvillian's.</para>
///
/// <para><b>SCOPE: N=5.</b> At N≥6 a {0,2}-coherence in a higher density block (the (2,2) sector)
/// can become the global slowest as a Q-dependent mode crossing, so the (0,1) block is no longer
/// guaranteed to carry the boundary; see the <c>birth_canal_horizon_junction</c> arc and
/// <c>simulations/birth_canal_n6_mode_crossing.py</c>. The V-Effect (Pauli-string weight w=N/2)
/// self-pair co-locates with that {0,2}-coherence, but it is a DIFFERENT decomposition
/// (RESOLVED 2026-06-14): total weight = n_diff + Z-shadow, so the dark {0,2}-coherence (peaking at w=N-1) is NOT
/// the w=N/2 self-pair. No aromaticity 4n-vs-4n+2 thesis is asserted (open; the C8 ring case breaks the
/// naive reading).</para>
///
/// <para><b>Tier1Derived honestly stated.</b> Like its sibling
/// <see cref="JDefectLightMigrationClaim"/>, this is a composition of proven structure plus a
/// bit-exact verification, scoped honestly: the invariant sub-block (from the conserved bi-grading)
/// and the flat-γ blindness (from −iQh anti-Hermitian) are DERIVED at every N; only the
/// "carries the global slowest across the surface" reading is N=5-VERIFIED (not yet generally
/// proven, because the cross-sector dominance at N≥6 is exactly the open junction).</para>
///
/// <para>Typed parent: <see cref="AbsorptionTheoremClaim"/> (the rate = −2·Σ_l γ_l·⟨Δ_l⟩ law this
/// block's −2·diag(γ) and flat-γ blindness rest on; the (0,1)-block rate is that law restricted to
/// a single conserved sector).</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_ABSORPTION_THEOREM.md</c> +
/// <c>simulations/birth_canal_vacuum_block_verifier.py</c> +
/// <c>compute/RCPsiSquared.Diagnostics/Foundation/SectorReductionWitness.cs</c>
/// (SectorReductionWitness, inspect --root reduction).</para></summary>
public sealed class VacuumBlockReductionClaim : Claim
{
    /// <summary>Typed parent: the Absorption Theorem. The (0,1) block's −2·diag(γ) decay and its
    /// flat-γ blindness are the absorption law restricted to a single conserved sector.</summary>
    public AbsorptionTheoremClaim Absorption { get; }

    /// <summary>One dense self-check tied to <see cref="SectorReductionWitness"/>.</summary>
    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    public IReadOnlyList<BatteryCase> Cases { get; }
    public int PassCount => Cases.Count(c => c.Passes);

    public VacuumBlockReductionClaim(AbsorptionTheoremClaim absorption)
        : base("Vacuum-block reduction: the birth-canal boundary's slowest mode is the odd |1-exc><vac| (0,1) coherence. " +
               "The (0,1) Liouville sector is an exact invariant sub-block (ket#/bra# bi-grading conserved by H_unit + Z-dephasing) - DERIVED; " +
               "its N-dim block L_(1,0) = -iQ·h - 2·diag(γ) carries the global slowest across the whole γ-surface at N=5 - VERIFIED bit-exact " +
               "(worst gap 5.9e-12, 80 pts; + SectorReductionWitness vs PostEpFlowField). Flat-γ blindness is analytic at every N (-iQh anti-Hermitian => Re λ = -2γ, Q-invariant). " +
               "SCOPE: N=5; at N>=6 a {0,2}-coherence (the (2,2) block) can become the global slowest (Q-dependent crossing, the birth_canal_horizon_junction arc). " +
               "The V-Effect (w=N/2) self-pair co-locates there but is a DIFFERENT decomposition (RESOLVED 2026-06-14: total weight = n_diff + Z-shadow; the dark {0,2}-coherence peaks at w=N-1, not w=N/2). No aromaticity 4n-vs-4n+2 thesis.",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_ABSORPTION_THEOREM.md + " +
               "simulations/birth_canal_vacuum_block_verifier.py + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/SectorReductionWitness.cs (SectorReductionWitness, inspect --root reduction)")
    {
        Absorption = absorption ?? throw new ArgumentNullException(nameof(absorption));
        Cases = BuildBattery();
    }

    public string InvariantSubBlock =>
        "The ket-excitation-number and bra-excitation-number are conserved by H_unit (XY hopping " +
        "preserves total excitation per factor) and by Z-dephasing (diagonal in the coherence " +
        "basis), so L is block-diagonal across the joint-popcount sectors. The (PCol=0, PRow=1) " +
        "block is an N-dim invariant sub-block: |1-exc><vac| coherences never mix with the rest, " +
        "and L_(1,0) = -iQ·h - 2·diag(γ) with h the single-particle hopping matrix.";

    public string FlatGammaBlindness =>
        "At uniform γ, L_(1,0) = -iQ·h - 2γ·I; -iQ·h is anti-Hermitian, so Re(λ) = -2γ for every " +
        "mode of the block, rate Q-invariant by uniformity alone. This is the Absorption Theorem " +
        "(-2γ per active bit) read inside the (0,1) sector; analytic at every N.";

    public string Scope =>
        "N=5. The (0,1) block IS the full boundary at N=5 (bit-exact). At N>=6 the global slowest " +
        "can cross to a {0,2}-coherence in the (2,2) density block (Q-dependent mode crossing, the " +
        "birth_canal_horizon_junction arc); the (0,1) block is then no longer guaranteed to carry " +
        "the boundary. The V-Effect (w=N/2) self-pair co-locates with that {0,2}-coherence but is a " +
        "DIFFERENT decomposition (RESOLVED 2026-06-14): n_diff = XY-weight, total Pauli weight = n_diff + " +
        "Z-shadow, so the dark {0,2}-coherence (total weight peaking at w=N-1) is NOT the w=N/2 self-pair. " +
        "No aromaticity 4n-vs-4n+2 thesis (open; C8 breaks the naive reading).";

    public override string DisplayName =>
        "Vacuum-block reduction: the N=5 birth-canal boundary = the |1-exc><vac| (0,1) block (Tier1Derived)";

    public override string Summary =>
        "the birth-canal boundary's slowest mode is the odd |1-exc><vac| (0,1) coherence; the (0,1) " +
        "sector is an exact invariant sub-block (DERIVED), carrying the global slowest at N=5 (VERIFIED " +
        "bit-exact); flat-γ blindness Re λ = -2γ analytic at every N; SCOPE N=5 (at N>=6 a {0,2}-coherence " +
        $"can win; V-Effect w=N/2 identity RESOLVED = distinct); {PassCount}/{Cases.Count} PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the invariant sub-block (DERIVED)", summary: InvariantSubBlock);
            yield return new InspectableNode("flat-γ blindness (analytic at every N)", summary: FlatGammaBlindness);
            yield return new InspectableNode("scope (N=5; the {0,2} junction; V-Effect identity RESOLVED = distinct)", summary: Scope);
            yield return new InspectableNode("verification (simulations/birth_canal_vacuum_block_verifier.py)",
                summary: "80 surface points: worst gap 5.9e-12 between the (0,1) sub-spectrum's slowest " +
                         "non-kernel rate and the full Liouvillian's; the live SectorReductionWitness pins " +
                         "the canal anchor bit-identical to PostEpFlowField at N=5.");
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
            yield return Absorption;   // typed parent edge
        }
    }

    /// <summary>Live battery through <see cref="SectorReductionWitness"/>: (a) the N=5 canal anchor
    /// reproduces the pinned boundary rates (the bit-exact reading the claim banks), (b) flat-γ
    /// blindness gives exactly 2γ at uniform γ, Q-invariant.</summary>
    private static IReadOnlyList<BatteryCase> BuildBattery()
    {
        const double tol = 1e-6;
        string Round(double v) => v.ToString("0.000000", System.Globalization.CultureInfo.InvariantCulture);

        var cases = new List<BatteryCase>();

        // (a) the N=5 canal anchor, bit-exact (the boundary the claim banks).
        var canal = new[] { 0.25, 1.5, 1.5, 1.5, 0.25 };
        double canalLo = SectorReductionWitness.VacBlockSlowest(5, 1.5, canal, TopologyKind.Chain);
        double canalHi = SectorReductionWitness.VacBlockSlowest(5, 1000.0, canal, TopologyKind.Chain);
        cases.Add(new BatteryCase(
            Name: "N=5 canal anchor (the banked boundary)",
            Detail: $"(0,1)-block slowest at canal profile, Q=1.5 -> {Round(canalLo)} and Q=1000 -> {Round(canalHi)} (matches PostEpFlowField)",
            Expected: $"{Round(1.2482918643729715)} and {Round(4.0 / 3.0)}",
            Actual: $"{Round(canalLo)} and {Round(canalHi)}"));

        // (b) flat-γ blindness: uniform γ=1 gives rate 2γ = 2, Q-invariant (analytic).
        foreach (int n in new[] { 5, 6 })
        {
            var uni = Enumerable.Repeat(1.0, n).ToArray();
            double lo = SectorReductionWitness.VacBlockSlowest(n, 1.5, uni, TopologyKind.Chain);
            double hi = SectorReductionWitness.VacBlockSlowest(n, 1000.0, uni, TopologyKind.Chain);
            bool ok = Math.Abs(lo - 2.0) < tol && Math.Abs(hi - 2.0) < tol;
            cases.Add(new BatteryCase(
                Name: $"flat-γ blindness at N={n} (Q-invariant 2γ)",
                Detail: $"uniform γ=1: rate {Round(lo)} (Q=1.5) and {Round(hi)} (Q=1000), expect 2γ=2 both",
                Expected: "2γ, Q-invariant",
                Actual: ok ? "2γ, Q-invariant" : $"{Round(lo)} / {Round(hi)}"));
        }

        return cases;
    }

    /// <summary>The shared instance with its parent built from the Pi2-Foundation root, mirroring
    /// <see cref="Core.Symmetry.ClockHandLadderClaim.Shared"/>: lets metadata be read without the
    /// full registry.</summary>
    public static VacuumBlockReductionClaim Shared { get; } =
        new VacuumBlockReductionClaim(new AbsorptionTheoremClaim(new Pi2DyadicLadderClaim()));
}
