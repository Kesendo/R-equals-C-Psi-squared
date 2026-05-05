using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>OOP carry of B2's per-(N, bond) cross-block numerical findings, tagged with
/// <see cref="BondClass"/>. One entry per bond; visible in the inspection tree under
/// <see cref="C2BondCoupling.CrossBlockWitnesses"/>.
///
/// <para><b>Why a witness, not a closed form:</b> the cross-block entries
/// V_b[α, j] = ⟨c_α | M_H_per_bond[b] | x_j⟩ depend on |u_0⟩, |v_0⟩, which are A3
/// Tier 2 (numerical SVD, library-tiebreaker-dependent at even N). The Frobenius norm
/// over the 2×2 cross-block is a magnitude indicator that splits Endpoint vs Interior
/// at c=2 (empirical anchor: HWHM_left/Q_peak Endpoint=0.7728 > Interior=0.7506,
/// EQ-022 (b1) 2026-05-02). Future Tier 1 promotion (per A3 PendingDerivationNote: lift
/// to projector-overlap onto the 2D top eigenspace) would replace these per-(N, b)
/// witnesses with a closed-form ratio.</para>
/// </summary>
public sealed record CrossBlockWitness(
    int Bond,
    BondClass BondClass,
    Complex EntryC1U0,    // V_b[0, 2] = ⟨c_1 | M | u_0⟩
    Complex EntryC1V0,    // V_b[0, 3] = ⟨c_1 | M | v_0⟩
    Complex EntryC3U0,    // V_b[1, 2] = ⟨c_3 | M | u_0⟩
    Complex EntryC3V0,    // V_b[1, 3] = ⟨c_3 | M | v_0⟩
    double FrobeniusNorm  // ‖cross-block 2×2‖_F
) : IInspectable
{
    public string DisplayName => $"cross-block witness b={Bond} ({BondClass})";

    public string Summary => $"‖V_b cross‖_F = {FrobeniusNorm:F4}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("bond", summary: Bond.ToString());
            yield return new InspectableNode("bond class", summary: BondClass.ToString());
            yield return InspectableNode.RealScalar("‖V_b cross‖_F", FrobeniusNorm);
            yield return new InspectableNode("V_b[0, 2] (⟨c_1 | M | u_0⟩)", summary: EntryC1U0.ToString("G6"));
            yield return new InspectableNode("V_b[0, 3] (⟨c_1 | M | v_0⟩)", summary: EntryC1V0.ToString("G6"));
            yield return new InspectableNode("V_b[1, 2] (⟨c_3 | M | u_0⟩)", summary: EntryC3U0.ToString("G6"));
            yield return new InspectableNode("V_b[1, 3] (⟨c_3 | M | v_0⟩)", summary: EntryC3V0.ToString("G6"));
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
