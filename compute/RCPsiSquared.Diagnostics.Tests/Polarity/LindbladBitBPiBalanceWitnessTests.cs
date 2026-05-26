using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.Polarity;

namespace RCPsiSquared.Diagnostics.Tests.Polarity;

public class LindbladBitBPiBalanceWitnessTests
{
    private static ChainSystem MakeChain() =>
        new ChainSystem(N: 2, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.XY, Topology: TopologyKind.Chain);

    [Fact]
    public void StandardSet_HasFiveWitnesses()
    {
        var set = LindbladBitBPiBalanceWitness.StandardSet(MakeChain());
        Assert.Equal(5, set.Count);
    }

    [Fact]
    public void StandardSet_NamesAreDistinct()
    {
        var set = LindbladBitBPiBalanceWitness.StandardSet(MakeChain());
        var names = set.Select(w => w.WitnessName).ToList();
        Assert.Equal(names.Count, names.Distinct().Count());
    }

    [Fact]
    public void StandardSet_TierIsTier1Derived()
    {
        // Every standard-set witness is anchored to F112 Tier1Derived (Hermitian H +
        // bit_b-homog c). Witness 4 sits in the broader empirical envelope; witness 5
        // is the in-typed-scope-violation structural counterexample. The Claim's Tier
        // tracks the F112-anchor tier, not the per-witness scope status.
        var set = LindbladBitBPiBalanceWitness.StandardSet(MakeChain());
        foreach (var w in set)
            Assert.Equal(Tier.Tier1Derived, w.Tier);
    }

    [Fact]
    public void Witness_Heisenberg_pure_Z_balanced_Matches()
    {
        var w = LindbladBitBPiBalanceWitness.StandardSet(MakeChain())[0];
        Assert.Equal("Heisenberg_pure_Z_balanced", w.WitnessName);
        Assert.Equal("BALANCED", w.ExpectedVerdict);
        Assert.Equal("BALANCED", w.ActualVerdict);
        Assert.True(w.Matches,
            $"witness '{w.WitnessName}' expected {w.ExpectedVerdict}, got {w.ActualVerdict} " +
            $"(rel asym = {w.ActualRelativeAsymmetry:E3})");
    }

    [Fact]
    public void Witness_YZ_ZY_pi2even_balanced_Matches()
    {
        var w = LindbladBitBPiBalanceWitness.StandardSet(MakeChain())[1];
        Assert.Equal("YZ_ZY_pi2even_balanced", w.WitnessName);
        Assert.Equal("BALANCED", w.ExpectedVerdict);
        Assert.Equal("BALANCED", w.ActualVerdict);
        Assert.True(w.Matches,
            $"witness '{w.WitnessName}' expected {w.ExpectedVerdict}, got {w.ActualVerdict} " +
            $"(rel asym = {w.ActualRelativeAsymmetry:E3})");
    }

    [Fact]
    public void Witness_XY_pi2odd_balanced_Matches()
    {
        var w = LindbladBitBPiBalanceWitness.StandardSet(MakeChain())[2];
        Assert.Equal("XY_pi2odd_balanced", w.WitnessName);
        Assert.Equal("BALANCED", w.ExpectedVerdict);
        Assert.Equal("BALANCED", w.ActualVerdict);
        Assert.True(w.Matches,
            $"witness '{w.WitnessName}' expected {w.ExpectedVerdict}, got {w.ActualVerdict} " +
            $"(rel asym = {w.ActualRelativeAsymmetry:E3})");
    }

    [Fact]
    public void Witness_Heisenberg_with_T1_envelope_balanced_Matches()
    {
        var w = LindbladBitBPiBalanceWitness.StandardSet(MakeChain())[3];
        Assert.Equal("Heisenberg_with_T1_envelope_balanced", w.WitnessName);
        Assert.Equal("BALANCED", w.ExpectedVerdict);
        Assert.Equal("BALANCED", w.ActualVerdict);
        Assert.True(w.Matches,
            $"witness '{w.WitnessName}' expected {w.ExpectedVerdict}, got {w.ActualVerdict} " +
            $"(rel asym = {w.ActualRelativeAsymmetry:E3})");
    }

    [Fact]
    public void Witness_Zdrive_with_T1_envelope_BROKEN_Matches()
    {
        // The structural payload of this commit: Z-drive + σ⁻ T1 BREAKS F112 balance
        // bit-exactly at rel asym ≈ 3.85e-3 (Welle 2 hardware-fit isolation,
        // f95_angle_steering Kingston 2026-05-16). The witness is expected BROKEN.
        var w = LindbladBitBPiBalanceWitness.StandardSet(MakeChain())[4];
        Assert.Equal("Zdrive_with_T1_envelope_BROKEN", w.WitnessName);
        Assert.Equal("BROKEN", w.ExpectedVerdict);
        Assert.Equal("BROKEN", w.ActualVerdict);
        Assert.True(w.Matches,
            $"witness '{w.WitnessName}' expected {w.ExpectedVerdict}, got {w.ActualVerdict} " +
            $"(rel asym = {w.ActualRelativeAsymmetry:E3}; expected ~3.85e-3)");
        Assert.True(w.ActualRelativeAsymmetry > 1e-6,
            $"Z-drive + T1 should produce structurally non-zero rel asym; got {w.ActualRelativeAsymmetry:E3}");
    }

    [Fact]
    public void Constructor_RejectsInvalidVerdict()
    {
        var chain = MakeChain();
        var terms = new[] { new PauliPairBondTerm(PauliLetter.X, PauliLetter.X) };
        Assert.Throws<ArgumentException>(() =>
            new LindbladBitBPiBalanceWitness(
                witnessName: "test",
                chain: chain,
                bondTerms: terms,
                gammaT1: null,
                expectedVerdict: "MAYBE"));
    }

    [Fact]
    public void Constructor_RejectsEmptyWitnessName()
    {
        var chain = MakeChain();
        var terms = new[] { new PauliPairBondTerm(PauliLetter.X, PauliLetter.X) };
        Assert.Throws<ArgumentException>(() =>
            new LindbladBitBPiBalanceWitness(
                witnessName: "",
                chain: chain,
                bondTerms: terms,
                gammaT1: null,
                expectedVerdict: "BALANCED"));
    }

    [Fact]
    public void Lazy_Polarity_DoesNotComputeUntilAccessed()
    {
        // Construction is cheap; the L-build runs only on first access to .Polarity / .ActualVerdict.
        var w = LindbladBitBPiBalanceWitness.StandardSet(MakeChain())[0];
        Assert.False(w.Polarity.IsValueCreated);
        _ = w.ActualVerdict;
        Assert.True(w.Polarity.IsValueCreated);
    }
}
