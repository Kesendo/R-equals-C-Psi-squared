using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>F87 (Tier 1, structural): the Pauli-pair trichotomy.
///
/// <code>
///     truly  iff  ‖M‖_F < ε                                              (M ≡ 0)
///     soft   iff  every λ ∈ Spec(L) has a partner −λ − 2σ ∈ Spec(L)      (palindrome only at spectral level)
///     hard   otherwise                                                    (no partner pairing)
/// </code>
///
/// <para>The classification uses F1's palindrome residual <c>M = Π·L·Π⁻¹ + L + 2σ·I</c> as the
/// discriminator. Equivalently in Π²-class language (cf. F79, F85): a term <c>(P, Q)</c> is
/// <c>truly</c> iff #Y is even AND #Z is even, <c>pi2_odd</c> iff <c>bit_b(P) + bit_b(Q)</c>
/// is odd, <c>pi2_even_nontruly</c> iff <c>bit_b(P) + bit_b(Q)</c> is even and not truly.
/// Mixed Hamiltonians refine the 3-way trichotomy into a 4-way classification
/// (truly / pi2_odd_pure / pi2_even_nontruly / mixed).</para>
///
/// <para>F87 is the **entry point of the F-chain**: F87 → F1 (residual M as discriminator) →
/// F49/F85 (‖M‖² closed forms) → F78/F79 (M-structure) → F80 (Π²-odd Bloch sign-walk) → F81
/// (Π·M·Π⁻¹ split) → F82/F84 (T1 / thermal corrections) → F83 (anti-fraction).</para>
///
/// <para>Origin (2026-04-24 to 2026-05-03): V_EFFECT_FINE_STRUCTURE counted 14 hard / 22
/// not-hard at N=3 over the 36-element bond-pair enumeration (commit 95386cd, 2026-04-25),
/// and the 22 not-hard split into 19 soft + 3 truly under the strict operator-equation test.
/// A separate 120-element ordered enumeration extended via Π-protected-observable testing
/// (commit 96ed6da N=4, 2026-04-26; commit 6438fef N=5, 2026-04-26) gave 15 hard / 46 soft /
/// 59 truly N-stable through N=3, 4, 5. Commit 81caf67 (2026-04-27) derived the partition
/// combinatorially from Pauli-pair compatibility rules. Marrakesh hardware Δ(soft − truly)
/// = −0.722 confirmation 2026-04-26 (job d7mjnjjaq2pc73a1pk4g). Registered as F87
/// retrospectively 2026-05-03.</para>
/// </summary>
public sealed class F87TrichotomyClassification : Claim
{
    public F87TrichotomyClassification()
        : base("F87 Pauli-pair trichotomy: truly / soft / hard via F1 residual",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F87 + experiments/V_EFFECT_FINE_STRUCTURE.md + experiments/MARRAKESH_THREE_LAYERS.md")
    { }

    public override string DisplayName => "F87: Pauli-pair trichotomy (truly / soft / hard)";

    public override string Summary =>
        "every Pauli-pair Hamiltonian under single-letter dephasing classifies as truly (‖M‖=0), soft (M≠0 but spectrum pairs), or hard (no pairing); F1 residual is the discriminator";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("discriminator",
                summary: "F1 palindrome residual M = Π·L·Π⁻¹ + L + 2σ·I; ‖M‖<ε → truly, else spectrum pairing test");
            yield return new InspectableNode("Π²-class refinement",
                summary: "truly / pi2_odd_pure / pi2_even_nontruly / mixed (4-way); bit_b parity per pair");
            yield return new InspectableNode("F-chain entry point",
                summary: "F87 → F1 (M) → F49/F85 (scaling) → F78/F79 (structure) → F80 (sign-walk) → F81 (Π-split) → F82/F84 (T1/thermal) → F83 (anti-fraction)");
            yield return new InspectableNode("Klein resonance",
                summary: "F87-hardness lives in the Klein cell matching the dephase letter (Z=(0,1), X=(1,0), Y=(1,1)); SU(2)-rotation-equivalent");
            yield return new InspectableNode("origin",
                summary: "36-enum N=3 → 14/19/3 (V_EFFECT_FINE_STRUCTURE 95386cd); 120-enum N=3,4,5 → 15/46/59 N-stable (96ed6da, 6438fef); combinatorial proof 81caf67; Marrakesh Δ(soft−truly)=−0.722");
            yield return new InspectableNode("classifier",
                summary: "RCPsiSquared.Diagnostics.F87.PauliPairTrichotomy.Classify");
        }
    }
}
