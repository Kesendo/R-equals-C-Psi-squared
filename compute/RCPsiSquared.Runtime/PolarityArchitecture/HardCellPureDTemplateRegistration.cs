using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="HardCellPureDTemplate"/>
/// (F111, 8th YParity-axis Claim, Tier1Derived, promoted 2026-06-10). Standalone Claim: no ctor
/// parents in this registration extension.
///
/// <para><b>Layer-boundary note</b>: F111's empirical anchor depends on the
/// F106 enumeration tool (in RCPsiSquared.Diagnostics.Tests) and structurally
/// references F110 (<see cref="HardCellYInversionPattern"/> in Core) as the
/// parent observation; the per-pair Pure-D Template Rule sharpens F110 Aspect
/// B (228:0 Y-inversion). Aspect A of F110 in turn rests on the colouring of
/// the three non-diagonal Klein cells (PROOF_F110 §2). The Core/Runtime layer
/// registers the Claim itself; the Diagnostics-layer F106 anchor and the
/// Diagnostics-layer SLOW_F111 enumeration test are independent. The proof's
/// parents F110 and F107 are not wired as <c>b.Get&lt;...&gt;()</c> calls; the Claim
/// names them in its docstring and anchor. Separately, the F106 enumeration tool
/// lives in <c>RCPsiSquared.Diagnostics.Tests</c>, which <c>RCPsiSquared.Runtime</c>
/// does not reference (per the PolarityCubeMap architectural boundary), and the F87
/// dissipator-resonance law is a census, not a parent. Same pattern as F110's
/// <see cref="HardCellYInversionPatternRegistration"/>.</para></summary>
public static class HardCellPureDTemplateRegistration
{
    public static ClaimRegistryBuilder RegisterHardCellPureDTemplate(
        this ClaimRegistryBuilder builder) =>
        builder.Register<HardCellPureDTemplate>(b =>
            new HardCellPureDTemplate(b.Get<KleinEightCellClaim>()));
}
