using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="HardCellPureDTemplate"/>
/// (F111, 8th YParity-axis Claim, Tier1Candidate). Standalone Claim: no ctor
/// parents in this registration extension.
///
/// <para><b>Layer-boundary note</b>: F111's empirical anchor depends on the
/// F106 enumeration tool (in RCPsiSquared.Diagnostics.Tests) and structurally
/// references F110 (<see cref="HardCellYInversionPattern"/> in Core) as the
/// parent observation; the per-pair Pure-D Template Rule sharpens F110 Aspect
/// B (228:0 Y-inversion). Aspect A of F110 in turn rests on F108 Part 1+2+3
/// (Π_5bilinear palindrome family across Z, X, Y dephasing) plus the F87
/// dissipator-resonance trichotomy. The Core/Runtime layer registers the
/// Claim itself; the Diagnostics-layer F106 anchor and the Diagnostics-layer
/// SLOW_F111 enumeration test are independent. Typed parent edges to F110 /
/// F108 / F107 are documented in prose in the Claim docstring but not wired
/// as <c>b.Get&lt;...&gt;()</c> calls because <c>RCPsiSquared.Runtime</c> does
/// not reference <c>RCPsiSquared.Diagnostics</c> (per the PolarityCubeMap
/// architectural boundary). Same pattern as F110's
/// <see cref="HardCellYInversionPatternRegistration"/>.</para></summary>
public static class HardCellPureDTemplateRegistration
{
    public static ClaimRegistryBuilder RegisterHardCellPureDTemplate(
        this ClaimRegistryBuilder builder) =>
        builder.Register<HardCellPureDTemplate>(b =>
            new HardCellPureDTemplate(b.Get<KleinEightCellClaim>()));
}
