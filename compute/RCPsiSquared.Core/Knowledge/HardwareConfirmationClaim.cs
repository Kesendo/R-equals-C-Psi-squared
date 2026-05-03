using RCPsiSquared.Core.Confirmations;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Core.Knowledge;

/// <summary>Typed <see cref="Claim"/> wrapper around a hardware-confirmed framework
/// prediction (<see cref="Confirmation"/>). Used by F-theorem knowledge bases to surface
/// their hardware-verified predictions as <see cref="Tier.Tier2Verified"/> nodes in the
/// inspection tree, alongside the algebraic Tier-1 claims and Tier-2 empirical witnesses.
///
/// <para>The wrapper is intentionally thin: it adapts an immutable <see cref="Confirmation"/>
/// record into the <see cref="Claim"/> shape so the same hardware entry can appear under
/// multiple KBs (e.g. <c>palindrome_trichotomy</c> sits under F1 and is referenced by F87).
/// </para>
/// </summary>
public sealed class HardwareConfirmationClaim : Claim
{
    public Confirmation Confirmation { get; }

    public HardwareConfirmationClaim(Confirmation confirmation)
        : base($"hardware confirmation: {confirmation.Name}",
               Tier.Tier2Verified,
               confirmation.ExperimentDoc)
    {
        Confirmation = confirmation;
    }

    public override string DisplayName =>
        $"[HW {Confirmation.Date}] {Confirmation.Name}";

    public override string Summary =>
        $"{Confirmation.Machine} ({Confirmation.Date}); observable {Confirmation.Observable}; " +
        $"primitive {Confirmation.FrameworkPrimitive}";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("date", summary: Confirmation.Date);
            yield return new InspectableNode("machine", summary: Confirmation.Machine);
            yield return new InspectableNode("job id", summary: Confirmation.JobId);
            yield return new InspectableNode("observable", summary: Confirmation.Observable);
            yield return new InspectableNode("predicted", summary: Confirmation.PredictedValue);
            yield return new InspectableNode("measured", summary: Confirmation.MeasuredValue);
            yield return new InspectableNode("framework primitive", summary: Confirmation.FrameworkPrimitive);
            yield return new InspectableNode("hardware data", summary: Confirmation.HardwareData);
            yield return new InspectableNode("experiment doc", summary: Confirmation.ExperimentDoc);
            yield return new InspectableNode("description", summary: Confirmation.Description);
        }
    }

    /// <summary>Look up multiple confirmations by name from <see cref="ConfirmationsRegistry"/>
    /// and wrap each found entry. Names that are not in the registry are silently skipped, so
    /// the array length may be less than the input list. Use this to attach a curated set of
    /// confirmations to an F-theorem knowledge base.</summary>
    public static IReadOnlyList<HardwareConfirmationClaim> LookupAll(IEnumerable<string> names) =>
        names.Select(ConfirmationsRegistry.Lookup)
             .Where(c => c is not null)
             .Select(c => new HardwareConfirmationClaim(c!))
             .ToArray();
}
