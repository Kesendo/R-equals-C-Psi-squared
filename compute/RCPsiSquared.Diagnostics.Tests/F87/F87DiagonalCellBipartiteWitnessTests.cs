using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F87;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>F103 §7: a diagonal-cell pair is soft iff H's hopping graph is bipartite in the
/// dephasing letter's eigenbasis (admits a chiral K with KHK=−H). These witnesses check that
/// the bipartite criterion (PredictedClass) agrees with the actual F87 trichotomy
/// (ActualClass via <see cref="PauliPairTrichotomy"/>) bit-exactly, for soft / hard-odd-cycle /
/// hard-template cases and across dephase letters.</summary>
public class F87DiagonalCellBipartiteWitnessTests
{
    private static ChainSystem MakeChainN4() => new(N: 4, J: 1.0, GammaZero: 0.05);

    private static F87DiagonalCellBipartiteWitness Named(string name) =>
        F87DiagonalCellBipartiteWitness.StandardSet(MakeChainN4()).Single(w => w.WitnessName == name);

    [Fact]
    public void StandardSet_HasFourWitnesses()
    {
        Assert.Equal(4, F87DiagonalCellBipartiteWitness.StandardSet(MakeChainN4()).Count);
    }

    [Fact]
    public void StandardSet_NamesAreDistinct()
    {
        var names = F87DiagonalCellBipartiteWitness.StandardSet(MakeChainN4())
            .Select(w => w.WitnessName).ToList();
        Assert.Equal(names.Count, names.Distinct().Count());
    }

    [Fact]
    public void StandardSet_TierIsTier1Candidate()
    {
        // bipartite ⟹ soft is derived; the converse (non-bipartite ⟹ hard) is verified, not
        // yet proved. Tier1Candidate: strong witness, missing one derivation piece.
        foreach (var w in F87DiagonalCellBipartiteWitness.StandardSet(MakeChainN4()))
            Assert.Equal(Tier.Tier1Candidate, w.Tier);
    }

    [Fact]
    public void EveryWitness_CriterionPredictsActualClass()
    {
        // The core claim: the bipartite criterion's soft/hard call equals the actual F87 verdict.
        foreach (var w in F87DiagonalCellBipartiteWitness.StandardSet(MakeChainN4()))
            Assert.True(w.CriterionAgrees,
                $"witness '{w.WitnessName}': predicted {w.PredictedClass}, actual {w.ActualClass}");
    }

    [Fact]
    public void EveryWitness_MatchesExpectedVerdict()
    {
        foreach (var w in F87DiagonalCellBipartiteWitness.StandardSet(MakeChainN4()))
            Assert.True(w.Matches,
                $"witness '{w.WitnessName}': expected {w.ExpectedClass}, actual {w.ActualClass}, bipartite={w.IsBipartite}");
    }

    [Fact]
    public void SoftWitness_IsBipartite_WithVerifiedChiralK()
    {
        var w = Named("XXZ_ZXX_Z_soft");
        Assert.Equal(PauliLetter.Z, w.DephaseLetter);
        Assert.True(w.IsBipartite);
        Assert.True(w.Criterion.Value.ChiralKVerified, "soft pair should admit a chiral K with KHK=−H");
        Assert.Equal(TrichotomyClass.Soft, w.ActualClass);
        Assert.Equal(TrichotomyClass.Soft, w.PredictedClass);
    }

    [Fact]
    public void HardOddCycleWitness_IsNotBipartite()
    {
        var w = Named("XXZ_XZX_Z_hard_oddcycle");
        Assert.False(w.IsBipartite);
        Assert.Equal(TrichotomyClass.Hard, w.ActualClass);
        Assert.Equal(TrichotomyClass.Hard, w.PredictedClass);
    }

    [Fact]
    public void HardTemplateWitness_IsNotBipartite_DiagonalLifted()
    {
        // ZZZ is a pure-D template: diagonal in the dephasing basis, so it lifts H's diagonal
        // and no chiral K can exist (rule (a)).
        var w = Named("ZZZ_XXZ_Z_hard_template");
        Assert.False(w.IsBipartite);
        Assert.Equal(TrichotomyClass.Hard, w.ActualClass);
    }

    [Fact]
    public void XDephWitness_RotationToEigenbasis_IsBipartiteAndSoft()
    {
        // ZZX+XZZ is the Z↔X mirror of the soft Z-deph witness; bipartite only in the X
        // eigenbasis, so this exercises the dephase-letter-relative rotation.
        var w = Named("ZZX_XZZ_X_soft");
        Assert.Equal(PauliLetter.X, w.DephaseLetter);
        Assert.True(w.IsBipartite);
        Assert.Equal(TrichotomyClass.Soft, w.ActualClass);
    }

    [Fact]
    public void Criterion_IsLazy_NotComputedUntilAccessed()
    {
        var w = Named("XXZ_ZXX_Z_soft");
        Assert.False(w.Criterion.IsValueCreated);
        _ = w.ActualClass;
        Assert.True(w.Criterion.IsValueCreated);
    }

    [Fact]
    public void Constructor_RejectsEmptyWitnessName()
    {
        var chain = MakeChainN4();
        var terms = new[] { new PauliTerm(new[] { PauliLetter.X, PauliLetter.X, PauliLetter.Z }, Complex.One) };
        Assert.Throws<ArgumentException>(() =>
            new F87DiagonalCellBipartiteWitness("", chain, terms, PauliLetter.Z, TrichotomyClass.Soft));
    }
}
