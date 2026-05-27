using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.Polarity;

namespace RCPsiSquared.Diagnostics.Tests.Polarity;

/// <summary>F112-Y witness tests: 5 witnesses under Y-dephase Π_Y polarity, the same
/// shape as <see cref="LindbladBitBPiBalanceWitnessTests"/> (F112-Z). Witnesses 1-4
/// expected BALANCED (in-scope and envelope), witness 5 expected BROKEN (Z-drive + T1
/// counterexample regime). Substantive M-norm-squared check + lazy-evaluation guarantee
/// + constructor guards.</summary>
public class LindbladBitBPiYBalanceWitnessTests
{
    private static ChainSystem MakeChain() =>
        new ChainSystem(N: 2, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.XY, Topology: TopologyKind.Chain);

    [Fact]
    public void StandardSet_HasFiveWitnesses()
    {
        var set = LindbladBitBPiYBalanceWitness.StandardSet(MakeChain());
        Assert.Equal(5, set.Count);
    }

    [Fact]
    public void StandardSet_NamesAreDistinct()
    {
        var set = LindbladBitBPiYBalanceWitness.StandardSet(MakeChain());
        var names = set.Select(w => w.WitnessName).ToList();
        Assert.Equal(names.Count, names.Distinct().Count());
    }

    [Fact]
    public void StandardSet_TierIsTier1Derived()
    {
        var set = LindbladBitBPiYBalanceWitness.StandardSet(MakeChain());
        foreach (var w in set)
            Assert.Equal(Tier.Tier1Derived, w.Tier);
    }

    [Fact]
    public void Witness_Heisenberg_pure_Y_balanced_Matches()
    {
        var w = LindbladBitBPiYBalanceWitness.StandardSet(MakeChain())[0];
        Assert.Equal("Heisenberg_pure_Y_balanced", w.WitnessName);
        Assert.Equal("BALANCED", w.ExpectedVerdict);
        Assert.Equal("BALANCED", w.ActualVerdict);
        Assert.True(w.Matches,
            $"witness '{w.WitnessName}' expected {w.ExpectedVerdict}, got {w.ActualVerdict} " +
            $"(rel asym = {w.ActualRelativeAsymmetry:E3})");
    }

    [Fact]
    public void Witness_YZ_ZY_pi2even_Y_balanced_Matches()
    {
        var w = LindbladBitBPiYBalanceWitness.StandardSet(MakeChain())[1];
        Assert.Equal("YZ_ZY_pi2even_Y_balanced", w.WitnessName);
        Assert.Equal("BALANCED", w.ExpectedVerdict);
        Assert.Equal("BALANCED", w.ActualVerdict);
        Assert.True(w.Matches,
            $"witness '{w.WitnessName}' expected {w.ExpectedVerdict}, got {w.ActualVerdict} " +
            $"(rel asym = {w.ActualRelativeAsymmetry:E3})");
    }

    [Fact]
    public void Witness_XY_pi2odd_Y_balanced_Matches()
    {
        var w = LindbladBitBPiYBalanceWitness.StandardSet(MakeChain())[2];
        Assert.Equal("XY_pi2odd_Y_balanced", w.WitnessName);
        Assert.Equal("BALANCED", w.ExpectedVerdict);
        Assert.Equal("BALANCED", w.ActualVerdict);
        Assert.True(w.Matches,
            $"witness '{w.WitnessName}' expected {w.ExpectedVerdict}, got {w.ActualVerdict} " +
            $"(rel asym = {w.ActualRelativeAsymmetry:E3})");
    }

    [Fact]
    public void Witness_Heisenberg_with_T1_envelope_Y_balanced_Matches()
    {
        var w = LindbladBitBPiYBalanceWitness.StandardSet(MakeChain())[3];
        Assert.Equal("Heisenberg_with_T1_envelope_Y_balanced", w.WitnessName);
        Assert.Equal("BALANCED", w.ExpectedVerdict);
        Assert.Equal("BALANCED", w.ActualVerdict);
        Assert.True(w.Matches,
            $"witness '{w.WitnessName}' expected {w.ExpectedVerdict}, got {w.ActualVerdict} " +
            $"(rel asym = {w.ActualRelativeAsymmetry:E3})");
    }

    [Fact]
    public void Witness_Zdrive_with_T1_envelope_Y_BROKEN_Matches()
    {
        // The Y-deph analog of the F112-Z f95 counterexample: Z-drive + σ⁻ T1 produces
        // non-Hermitian Π_Y-eigenspace coupling (the F113 [Z, σ⁻] = −2σ⁻ commutator
        // structure is dephase-letter independent at the H × c level).
        var w = LindbladBitBPiYBalanceWitness.StandardSet(MakeChain())[4];
        Assert.Equal("Zdrive_with_T1_envelope_Y_BROKEN", w.WitnessName);
        Assert.Equal("BROKEN", w.ExpectedVerdict);
        Assert.Equal("BROKEN", w.ActualVerdict);
        Assert.True(w.Matches,
            $"witness '{w.WitnessName}' expected {w.ExpectedVerdict}, got {w.ActualVerdict} " +
            $"(rel asym = {w.ActualRelativeAsymmetry:E3}; expected structurally > 1e-6)");
        Assert.True(w.ActualRelativeAsymmetry > 1e-6,
            $"Z-drive + T1 under Y-deph should produce structurally non-zero rel asym; " +
            $"got {w.ActualRelativeAsymmetry:E3}");
    }

    [Fact]
    public void Substantive_Y_Dephase_Asymmetry_Zero_For_Bit_B_Homogeneous_XY_At_N_Equals_2()
    {
        // Use witness 2 (XY) for the substantive check: XY is bit_b sum 0+1 = 1
        // (Π²_Y-odd / non-truly under Y-deph) so M is non-trivial, unlike Heisenberg
        // / YZ+ZY which are Π²_Y-even (truly) and have M = 0.
        var w = LindbladBitBPiYBalanceWitness.StandardSet(MakeChain())[2];
        var pol = w.Polarity.Value;
        Assert.True(pol.MNormSquared > 1e-6,
            $"M should be non-trivial; got ‖M‖² = {pol.MNormSquared:E3}");
        Assert.True(w.ActualRelativeAsymmetry < 1e-12,
            $"F112-Y in-scope rel asym should be bit-exact 0; got {w.ActualRelativeAsymmetry:E3}");
    }

    [Fact]
    public void Constructor_RejectsInvalidVerdict()
    {
        var chain = MakeChain();
        var terms = new[] { new PauliPairBondTerm(PauliLetter.X, PauliLetter.X) };
        Assert.Throws<ArgumentException>(() =>
            new LindbladBitBPiYBalanceWitness(
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
            new LindbladBitBPiYBalanceWitness(
                witnessName: "",
                chain: chain,
                bondTerms: terms,
                gammaT1: null,
                expectedVerdict: "BALANCED"));
    }

    [Fact]
    public void Lazy_Polarity_DoesNotComputeUntilAccessed()
    {
        var w = LindbladBitBPiYBalanceWitness.StandardSet(MakeChain())[0];
        Assert.False(w.Polarity.IsValueCreated);
        _ = w.ActualVerdict;
        Assert.True(w.Polarity.IsValueCreated);
    }
}
