using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Tests for the Q-axis QAnchorMap (Q = J/γ₀ structural anchors). Verifies
/// the 9 canonical anchors, the J = Q·γ₀ identity, the F95 θ = arctan(Q) angle, band
/// queries, and the wave-breaking-scan subset.</summary>
public class QAnchorMapTests
{
    private readonly ITestOutputHelper _out;

    public QAnchorMapTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void CanonicalMap_HasTenAnchors()
    {
        var m = new QAnchorMap();
        Assert.Equal(10, m.Anchors.Count);
    }

    [Fact]
    public void CanonicalAnchors_IncludeAllNamedQValues_AndSqrt3()
    {
        var m = new QAnchorMap();
        var qValues = m.Anchors.Select(a => a.Q).OrderBy(q => q).ToArray();
        Assert.Equal(
            new[] { 0.2, 0.35, 1.0, 1.2, 1.5, 1.6, Math.Sqrt(3.0), 1.8, 2.0, 2.5 },
            qValues);
    }

    [Fact]
    public void JAtGamma0Point05_MatchesQTimesGamma0()
    {
        var m = new QAnchorMap();
        foreach (var a in m.Anchors)
        {
            Assert.Equal(a.Q * 0.05, a.JAtGamma0Point05, precision: 12);
        }
    }

    [Fact]
    public void JAt_SubstrateInvariant_ComputesQTimesGamma0()
    {
        var m = new QAnchorMap();
        var balance = m.AnkerAt(1.0);
        Assert.NotNull(balance);
        Assert.Equal(0.025, balance!.JAt(0.025), precision: 12);
        Assert.Equal(0.10, balance.JAt(0.10), precision: 12);
        // Substrate invariance per UniversalCarrierClaim: same Q, different γ₀, same role.
    }

    [Fact]
    public void ThetaDegrees_AtBalance_IsFortyFiveDegrees()
    {
        var m = new QAnchorMap();
        var balance = m.AnkerAt(1.0);
        Assert.NotNull(balance);
        Assert.Equal(45.0, balance!.ThetaDegrees(), precision: 10);
    }

    [Fact]
    public void ThetaDegrees_AtQ2_IsArctan2()
    {
        var m = new QAnchorMap();
        var qep = m.AnkerAt(2.0);
        Assert.NotNull(qep);
        Assert.Equal(Math.Atan(2.0) * 180.0 / Math.PI, qep!.ThetaDegrees(), precision: 10);
        // ~63.4°
    }

    [Fact]
    public void AnkerAt_FindsByQValue_WithinTolerance()
    {
        var m = new QAnchorMap();
        Assert.NotNull(m.AnkerAt(1.5));
        Assert.NotNull(m.AnkerAt(2.5));
        Assert.Null(m.AnkerAt(3.0));  // not a canonical anchor
    }

    [Fact]
    public void OnsetBand_ContainsTwoAnchors()
    {
        var m = new QAnchorMap();
        var onset = m.ByBand(QBand.Onset);
        Assert.Equal(2, onset.Count);
        Assert.Contains(onset, a => Math.Abs(a.Q - 0.2) < 1e-12);
        Assert.Contains(onset, a => Math.Abs(a.Q - 0.35) < 1e-12);
    }

    [Fact]
    public void PeakBand_ContainsFiveAnchors_IncludingSqrt3CanonicalAngleAnchor()
    {
        // 1.2 peak start + 1.5 c=2 Q_peak + 1.6 c=3 Q_peak + √3 canonical θ=60°
        // (via LindbladAbsorptionMatchAtSixtyDegreesClaim) + 1.8 c=4/5 Q_peak / peak end
        var m = new QAnchorMap();
        var peak = m.ByBand(QBand.Peak);
        Assert.Equal(5, peak.Count);
        Assert.Contains(peak, a => Math.Abs(a.Q - 1.2) < 1e-12);
        Assert.Contains(peak, a => Math.Abs(a.Q - 1.5) < 1e-12);
        Assert.Contains(peak, a => Math.Abs(a.Q - 1.6) < 1e-12);
        Assert.Contains(peak, a => Math.Abs(a.Q - Math.Sqrt(3.0)) < 1e-12);
        Assert.Contains(peak, a => Math.Abs(a.Q - 1.8) < 1e-12);
    }

    [Fact]
    public void BalanceBand_ContainsSingleAnchor_AtQOne()
    {
        var m = new QAnchorMap();
        var balance = m.ByBand(QBand.Balance);
        Assert.Single(balance);
        Assert.Equal(1.0, balance[0].Q);
    }

    [Fact]
    public void EndpointOrbitBand_ContainsTier1DerivedAtQ25()
    {
        var m = new QAnchorMap();
        var endpoint = m.ByBand(QBand.EndpointOrbit);
        Assert.Single(endpoint);
        Assert.Equal(2.5, endpoint[0].Q);
        Assert.Equal(Tier.Tier1Derived, endpoint[0].Tier);
    }

    [Fact]
    public void Tier1DerivedAnchors_AreFiveSchemaGroundedAnchors()
    {
        // The five schema-derived Q-anchors: Balance (J=γ₀), Interior pole (2−1/2),
        // Lindblad-Absorption-Match (√3 at θ=60°), Q_EP idealized (2), Endpoint pole (2+1/2).
        var m = new QAnchorMap();
        var tier1 = m.ByTier(Tier.Tier1Derived);
        Assert.Equal(5, tier1.Count);
        Assert.Contains(tier1, a => Math.Abs(a.Q - 1.0) < 1e-12);
        Assert.Contains(tier1, a => Math.Abs(a.Q - 1.5) < 1e-12);
        Assert.Contains(tier1, a => Math.Abs(a.Q - Math.Sqrt(3.0)) < 1e-12);
        Assert.Contains(tier1, a => Math.Abs(a.Q - 2.0) < 1e-12);
        Assert.Contains(tier1, a => Math.Abs(a.Q - 2.5) < 1e-12);
    }

    [Fact]
    public void WaveBreakingScanSubset_HasThreeAnchors()
    {
        // The original FractionReferenceGraph.QBasisAnkers { 1.0, 1.5, 2.0 } — the
        // anchors with documented dynamical-visibility in wave_breaking_q_anchor_scan.py.
        var m = new QAnchorMap();
        var scanSubset = m.WaveBreakingScanSubset;
        Assert.Equal(3, scanSubset.Count);
        var scanQValues = scanSubset.Select(a => a.Q).OrderBy(q => q).ToArray();
        Assert.Equal(new[] { 1.0, 1.5, 2.0 }, scanQValues);
    }

    [Fact]
    public void Render_EmitsAllAnchors()
    {
        var m = new QAnchorMap();
        var rendered = m.Render();
        Assert.Contains("Q-Anchor Map", rendered);
        Assert.Contains("Total anchors: 10", rendered);
        Assert.Contains("Wave-breaking-scan subset: 3", rendered);
        Assert.Contains("Balance", rendered);
        Assert.Contains("F86 Q_peak (c=2)", rendered);
        Assert.Contains("Endpoint orbit", rendered);
        _out.WriteLine(rendered);
    }

    [Fact]
    public void EmptyConstructor_AllowsCustomAnchors()
    {
        var custom = new[]
        {
            new QBasisAnker(0.5, 0.025, QBand.Onset, "test", Tier.OpenQuestion, "test"),
        };
        var m = new QAnchorMap(custom);
        Assert.Single(m.Anchors);
    }

    [Fact]
    public void NullAnchors_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => new QAnchorMap(null!));
    }
}
