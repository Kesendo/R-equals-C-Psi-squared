namespace RCPsiSquared.Core.Inspection;

/// <summary>Where a node's displayed value comes from, the [live]/[stored] provenance badge
/// (stranger_door's fifth door, the doormat for the outside reader).
///
/// <para><see cref="Live"/> = recomputed from state at inspect time: a witness rebuilding a
/// matrix, recounting a census, re-running an eigendecomposition. <see cref="Stored"/> = written
/// from a proof, registry, closed form, or banked artifact and surfaced as-is.</para>
///
/// <para>The rule is <b>value-origin honesty</b>: the badge marks where the NUMBER comes from,
/// per node, at the point it is computed; it is NOT the code path that surfaces it (a live
/// getter reading a banked file is Stored) and it is NOT inherited down a derived chain (a child
/// recomputed live inside a live witness is Live even though its parent's law was Stored).
/// The default direction is conservative: an unstamped frozen carrier reports Stored, which can
/// only ever under-claim liveness, never falsely claim it.</para>
/// </summary>
public enum NodeProvenance
{
    /// <summary>Recomputed from state at inspect time.</summary>
    Live,

    /// <summary>Written from a proof, registry, closed form, or banked artifact; surfaced as-is.</summary>
    Stored,
}
