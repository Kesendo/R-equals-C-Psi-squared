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
    public void Witness_Heisenberg_pure_X_DEGENERATE_Matches()
    {
        var w = LindbladBitAPiBalanceWitness.StandardSet(MakeChain())[0];
        Assert.Equal("Heisenberg_pure_X_DEGENERATE", w.WitnessName);
        // Structural, not a float reading: decided from the Pauli letters, so it holds
        // at every N. The float IsDegenerate only reaches exact 0.0 at N = 2.
        Assert.True(w.IsStructurallyDegenerate,
            "Heisenberg_pure_X: Pi^2-even H with homogeneous c carries no polarity content");
        Assert.Equal("DEGENERATE", w.ExpectedVerdict);
        Assert.Equal("DEGENERATE", w.ActualVerdict);
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
    public void Witness_XY_bit_a_even_DEGENERATE_Matches()
    {
        var w = LindbladBitAPiBalanceWitness.StandardSet(MakeChain())[2];
        Assert.Equal("XY_bit_a_even_DEGENERATE", w.WitnessName);
        // Structural, not a float reading: decided from the Pauli letters, so it holds
        // at every N. The float IsDegenerate only reaches exact 0.0 at N = 2.
        Assert.True(w.IsStructurallyDegenerate,
            "XY_bit_a_even: Pi^2-even H with homogeneous c carries no polarity content");
        Assert.Equal("DEGENERATE", w.ExpectedVerdict);
        Assert.Equal("DEGENERATE", w.ActualVerdict);
        Assert.True(w.Matches,
            $"witness '{w.WitnessName}' expected {w.ExpectedVerdict}, got {w.ActualVerdict} " +
            $"(rel asym = {w.ActualRelativeAsymmetry:E3})");
    }

    [Fact]
    public void Witness_Heisenberg_with_T1_envelope_DEGENERATE_Matches()
    {
        var w = LindbladBitAPiBalanceWitness.StandardSet(MakeChain())[3];
        Assert.Equal("Heisenberg_with_T1_envelope_DEGENERATE", w.WitnessName);
        // Structural, not a float reading: decided from the Pauli letters, so it holds
        // at every N. The float IsDegenerate only reaches exact 0.0 at N = 2.
        Assert.True(w.IsStructurallyDegenerate,
            "Heisenberg_with_T1_envelope: Pi^2-even H with homogeneous c carries no polarity content");
        Assert.Equal("DEGENERATE", w.ExpectedVerdict);
        Assert.Equal("DEGENERATE", w.ActualVerdict);
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
    public void Vacuous_X_Dephase_Asymmetry_Zero_For_Bit_A_Homogeneous_Heisenberg_At_N_Equals_2()
    {
        // This test used to "harden" the zero below by asserting ‖M‖² > 1e-6, on the reading
        // that a non-trivial M makes the zero substantive. Measured 2026-08-06, that reading
        // is wrong here, and the numbers are the argument: ‖M‖² = 3.200000e-1, comfortably
        // over the old bar, while ‖M_anti‖² is EXACTLY 0. The asymmetry is a difference of
        // the two halves of M_anti and can see nothing else, so ‖M‖² cannot certify that a
        // zero asymmetry means anything; here it does not. Keep the ‖M‖² line as the record
        // of that (it is a true fact about this fixture, just not a hardening), and assert
        // the honest reading beside it.
        var w = LindbladBitAPiBalanceWitness.StandardSet(MakeChain())[0];
        var pol = w.Polarity.Value;
        Assert.True(pol.MNormSquared > 1e-6,
            $"M is non-trivial on this fixture; got ‖M‖² = {pol.MNormSquared:E3}");
        Assert.True(pol.MAntiNormSquared == 0.0,
            $"Heisenberg is Π²-even ⇒ no polarity content; got ‖M_anti‖² = {pol.MAntiNormSquared:E3}");
        Assert.True(w.IsDegenerate,
            "so the zero below is 0 = 0, not a balance the tolerance earned; this witness is "
            + "NOT a substantive test of F112-X and must not be cited as one");
        Assert.True(double.IsNaN(w.ActualRelativeAsymmetry),
            "degenerate ⇒ NO ratio at all (0/0), and since 2026-08-07 the field says so "
            + $"instead of reporting a 0.0 that reads as balance; got {w.ActualRelativeAsymmetry:E3}");
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
        // Construction is cheap; the L-build runs only on first access to .Polarity /
        // .ActualVerdict. Witness [1] on purpose: it is NOT structurally degenerate, so
        // reaching its verdict genuinely needs the decomposition. Witness [0] would pass
        // this test for the wrong reason after the DEGENERATE short-circuit landed.
        var w = LindbladBitAPiBalanceWitness.StandardSet(MakeChain())[1];
        Assert.Equal("YZ_ZY_bit_a_odd_balanced", w.WitnessName);
        Assert.False(w.IsStructurallyDegenerate);
        Assert.False(w.Polarity.IsValueCreated);
        _ = w.ActualVerdict;
        Assert.True(w.Polarity.IsValueCreated);
    }

    [Fact]
    public void Degenerate_Witness_Reaches_Its_Verdict_Without_Building_L()
    {
        // The structural test is exact and cheap, so a silent configuration is known to be
        // silent before any Liouvillian is assembled. This pins that the short-circuit is
        // real and not incidental: no decomposition is forced, and the verdict is not
        // BALANCED.
        var w = LindbladBitAPiBalanceWitness.StandardSet(MakeChain())[0];
        Assert.True(w.IsStructurallyDegenerate);
        Assert.False(w.Polarity.IsValueCreated);
        Assert.Equal("DEGENERATE", w.ActualVerdict);
        Assert.False(w.Polarity.IsValueCreated);
    }
}
