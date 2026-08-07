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
    public void Witness_Heisenberg_pure_Z_DEGENERATE_Matches()
    {
        var w = LindbladBitBPiBalanceWitness.StandardSet(MakeChain())[0];
        Assert.Equal("Heisenberg_pure_Z_DEGENERATE", w.WitnessName);
        // Structural, not a float reading: decided from the Pauli letters, so it holds
        // at every N. The float IsDegenerate only reaches exact 0.0 at N = 2.
        Assert.True(w.IsStructurallyDegenerate,
            "Heisenberg_pure_Z: Pi^2-even H with homogeneous c carries no polarity content");
        Assert.Equal("DEGENERATE", w.ExpectedVerdict);
        Assert.Equal("DEGENERATE", w.ActualVerdict);
        Assert.True(w.Matches,
            $"witness '{w.WitnessName}' expected {w.ExpectedVerdict}, got {w.ActualVerdict} " +
            $"(rel asym = {w.ActualRelativeAsymmetry:E3})");

        // Heisenberg is truly, so M itself vanishes by the F83 closed form and there is no
        // polarity content: this BALANCED is 0 = 0, carrying no information about the
        // tolerance. Recorded rather than left implicit, because for two and a half months
        // the ‖M‖² denominator hid it behind a 1e-15 floor.
        Assert.True(w.IsDegenerate, "Heisenberg + pure Z ⇒ M ≡ 0 ⇒ no polarity content");
        Assert.True(w.Polarity.Value.MNormSquared == 0.0,
            $"truly ⇒ M ≡ 0; got ‖M‖² = {w.Polarity.Value.MNormSquared:E3}");
    }

    [Fact]
    public void Witness_YZ_ZY_pi2even_DEGENERATE_Matches()
    {
        var w = LindbladBitBPiBalanceWitness.StandardSet(MakeChain())[1];
        Assert.Equal("YZ_ZY_pi2even_DEGENERATE", w.WitnessName);
        // Structural, not a float reading: decided from the Pauli letters, so it holds
        // at every N. The float IsDegenerate only reaches exact 0.0 at N = 2.
        Assert.True(w.IsStructurallyDegenerate,
            "YZ_ZY_pi2even: Pi^2-even H with homogeneous c carries no polarity content");
        Assert.Equal("DEGENERATE", w.ExpectedVerdict);
        Assert.Equal("DEGENERATE", w.ActualVerdict);
        Assert.True(w.Matches,
            $"witness '{w.WitnessName}' expected {w.ExpectedVerdict}, got {w.ActualVerdict} " +
            $"(rel asym = {w.ActualRelativeAsymmetry:E3})");

        // Degenerate for a different reason than witness 1, and this is the pair that shows
        // ‖M‖² was the wrong scale: YZ+ZY is Π²-EVEN but not truly, so ‖M‖² is LARGE while
        // M_anti = L_{H_odd} is still exactly zero. The old denominator divided by that
        // large number and the vacuity was invisible.
        Assert.True(w.IsDegenerate, "Π²-even ⇒ H_odd = 0 ⇒ M_anti = 0 ⇒ no polarity content");
        Assert.True(w.Polarity.Value.MNormSquared > 1.0,
            $"non-truly Π²-even ⇒ ‖M‖² is large; got {w.Polarity.Value.MNormSquared:E3}");
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
        // The structural payload: Z-drive + σ⁻ T1 BREAKS F112 balance for these parameters
        // (ω = 0.13, γ_T1 = 0.001, N = 2). What F113 predicts is the ABSOLUTE asymmetry,
        // (4^N/2)·Σ_l ω_l·(γ_pump,l − γ_T1,l); the ratio below is that divided by the
        // polarity content. Neither is a hardware reading. The witness is expected BROKEN.
        var w = LindbladBitBPiBalanceWitness.StandardSet(MakeChain())[4];
        Assert.Equal("Zdrive_with_T1_envelope_BROKEN", w.WitnessName);
        Assert.Equal("BROKEN", w.ExpectedVerdict);
        Assert.Equal("BROKEN", w.ActualVerdict);
        Assert.True(w.Matches,
            $"witness '{w.WitnessName}' expected {w.ExpectedVerdict}, got {w.ActualVerdict} " +
            $"(rel asym = {w.ActualRelativeAsymmetry:E3})");

        // A single-site Z-drive is Π²-odd, so unlike witnesses 1 and 2 this one has real
        // polarity content and its ratio is a genuine contrast, not a 0/0.
        Assert.False(w.IsDegenerate,
            $"Z-drive is Π²-odd; ‖M_anti‖² = {w.Polarity.Value.MAntiNormSquared:E3}");
        Assert.True(w.ActualRelativeAsymmetry > 1e-6,
            $"Z-drive + T1 should produce structurally non-zero rel asym; got {w.ActualRelativeAsymmetry:E3}");

        // Value pin, so the number quoted in the witness docstring stays honest. It moved
        // when the denominator changed from max(‖M‖², 1e-15) to ‖M_anti‖²: on this fixture
        // ‖M‖² = 5.408880e-1 against ‖M_anti‖² = 2.704080e-1, so the old reading was
        // 3.845528e-3, a factor 2.000266 lower. (NOT exactly half: ‖M‖² = ‖M_zero‖² +
        // ‖M_anti‖², and the two parts merely happen to be close on this fixture.) This
        // line therefore also fails on a revert.
        Assert.Equal(7.69208011597269e-3, w.ActualRelativeAsymmetry, precision: 12);
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
        // Construction is cheap; the L-build runs only on first access to .Polarity /
        // .ActualVerdict. Witness [2] on purpose: it is NOT structurally degenerate, so
        // reaching its verdict genuinely needs the decomposition. Witness [0] would pass
        // this test for the wrong reason after the DEGENERATE short-circuit landed.
        var w = LindbladBitBPiBalanceWitness.StandardSet(MakeChain())[2];
        Assert.Equal("XY_pi2odd_balanced", w.WitnessName);
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
        var w = LindbladBitBPiBalanceWitness.StandardSet(MakeChain())[0];
        Assert.True(w.IsStructurallyDegenerate);
        Assert.False(w.Polarity.IsValueCreated);
        Assert.Equal("DEGENERATE", w.ActualVerdict);
        Assert.False(w.Polarity.IsValueCreated);
    }
}
