using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.Polarity;

namespace RCPsiSquared.Diagnostics.Tests.Polarity;

/// <summary>F112-X witness tests: 5 BALANCED witnesses under X-dephase Π_X polarity,
/// substantive M-norm-squared checks, lazy-evaluation guarantee, and constructor
/// guards. Mirrors <see cref="LindbladBitBPiBalanceWitnessTests"/>; the F112-X
/// counterexample-BROKEN slot is intentionally absent (see witness 5 docstring in
/// <see cref="LindbladBitAPiBalanceWitness.StandardSet"/>).</summary>
public class LindbladBitAPiBalanceWitnessTests
{
    private static ChainSystem MakeChain() =>
        new ChainSystem(N: 2, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.XY, Topology: TopologyKind.Chain);

    [Fact]
    public void StandardSet_HasFiveWitnesses()
    {
        var set = LindbladBitAPiBalanceWitness.StandardSet(MakeChain());
        Assert.Equal(5, set.Count);
    }

    [Fact]
    public void StandardSet_NamesAreDistinct()
    {
        var set = LindbladBitAPiBalanceWitness.StandardSet(MakeChain());
        var names = set.Select(w => w.WitnessName).ToList();
        Assert.Equal(names.Count, names.Distinct().Count());
    }

    [Fact]
    public void StandardSet_TierIsTier1Derived()
    {
        var set = LindbladBitAPiBalanceWitness.StandardSet(MakeChain());
        foreach (var w in set)
            Assert.Equal(Tier.Tier1Derived, w.Tier);
    }

    [Fact]
    public void Witness_Heisenberg_pure_X_balanced_Matches()
    {
        var w = LindbladBitAPiBalanceWitness.StandardSet(MakeChain())[0];
        Assert.Equal("Heisenberg_pure_X_balanced", w.WitnessName);
        Assert.Equal("BALANCED", w.ExpectedVerdict);
        Assert.Equal("BALANCED", w.ActualVerdict);
        Assert.True(w.Matches,
            $"witness '{w.WitnessName}' expected {w.ExpectedVerdict}, got {w.ActualVerdict} " +
            $"(rel asym = {w.ActualRelativeAsymmetry:E3})");
    }

    [Fact]
    public void Witness_YZ_ZY_bit_a_odd_balanced_Matches()
    {
        var w = LindbladBitAPiBalanceWitness.StandardSet(MakeChain())[1];
        Assert.Equal("YZ_ZY_bit_a_odd_balanced", w.WitnessName);
        Assert.Equal("BALANCED", w.ExpectedVerdict);
        Assert.Equal("BALANCED", w.ActualVerdict);
        Assert.True(w.Matches,
            $"witness '{w.WitnessName}' expected {w.ExpectedVerdict}, got {w.ActualVerdict} " +
            $"(rel asym = {w.ActualRelativeAsymmetry:E3})");
    }

    [Fact]
    public void Witness_XY_bit_a_even_balanced_Matches()
    {
        var w = LindbladBitAPiBalanceWitness.StandardSet(MakeChain())[2];
        Assert.Equal("XY_bit_a_even_balanced", w.WitnessName);
        Assert.Equal("BALANCED", w.ExpectedVerdict);
        Assert.Equal("BALANCED", w.ActualVerdict);
        Assert.True(w.Matches,
            $"witness '{w.WitnessName}' expected {w.ExpectedVerdict}, got {w.ActualVerdict} " +
            $"(rel asym = {w.ActualRelativeAsymmetry:E3})");
    }

    [Fact]
    public void Witness_Heisenberg_with_T1_envelope_balanced_Matches()
    {
        var w = LindbladBitAPiBalanceWitness.StandardSet(MakeChain())[3];
        Assert.Equal("Heisenberg_with_T1_envelope_balanced", w.WitnessName);
        Assert.Equal("BALANCED", w.ExpectedVerdict);
        Assert.Equal("BALANCED", w.ActualVerdict);
        Assert.True(w.Matches,
            $"witness '{w.WitnessName}' expected {w.ExpectedVerdict}, got {w.ActualVerdict} " +
            $"(rel asym = {w.ActualRelativeAsymmetry:E3})");
    }

    [Fact]
    public void Witness_Xdrive_with_T1_envelope_balanced_Matches()
    {
        // X-axis analog of the BitB Zdrive+T1 BROKEN witness. Under X-deph, X-drive
        // aligns with the dephase axis (both bit_a=1) and σ⁻ is also bit_a=1, so c is
        // bit_a-homogeneous and F112-X holds. Demonstrates the F113-style break
        // mechanism is BIT_A-AXIS-FREE.
        var w = LindbladBitAPiBalanceWitness.StandardSet(MakeChain())[4];
        Assert.Equal("Xdrive_with_T1_envelope_balanced", w.WitnessName);
        Assert.Equal("BALANCED", w.ExpectedVerdict);
        Assert.Equal("BALANCED", w.ActualVerdict);
        Assert.True(w.Matches,
            $"witness '{w.WitnessName}' expected {w.ExpectedVerdict}, got {w.ActualVerdict} " +
            $"(rel asym = {w.ActualRelativeAsymmetry:E3})");
    }

    [Fact]
    public void Substantive_X_Dephase_Asymmetry_Zero_For_Bit_A_Homogeneous_Heisenberg_At_N_Equals_2()
    {
        // Hardening test: assert M is non-trivial so the bit-exact 0 assertion is substantive.
        var w = LindbladBitAPiBalanceWitness.StandardSet(MakeChain())[0];
        var pol = w.Polarity.Value;
        Assert.True(pol.MNormSquared > 1e-6,
            $"M should be non-trivial; got ‖M‖² = {pol.MNormSquared:E3}");
        Assert.True(w.ActualRelativeAsymmetry < 1e-12,
            $"F112-X in-scope rel asym should be bit-exact 0; got {w.ActualRelativeAsymmetry:E3}");
    }

    [Fact]
    public void Constructor_RejectsInvalidVerdict()
    {
        var chain = MakeChain();
        var terms = new[] { new PauliPairBondTerm(PauliLetter.X, PauliLetter.X) };
        Assert.Throws<ArgumentException>(() =>
            new LindbladBitAPiBalanceWitness(
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
            new LindbladBitAPiBalanceWitness(
                witnessName: "",
                chain: chain,
                bondTerms: terms,
                gammaT1: null,
                expectedVerdict: "BALANCED"));
    }

    [Fact]
    public void Lazy_Polarity_DoesNotComputeUntilAccessed()
    {
        var w = LindbladBitAPiBalanceWitness.StandardSet(MakeChain())[0];
        Assert.False(w.Polarity.IsValueCreated);
        _ = w.ActualVerdict;
        Assert.True(w.Polarity.IsValueCreated);
    }
}
