using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F87;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>The Pauli-basis convention guard for the DERIVED k-body hidden-Q routing soft-certifier
/// (<see cref="KBodyPalindromeRouting"/>). The per-term k-site anticommutator check is a NEW Pauli-basis
/// construction (Q_k and [T,·]_k as 4^k × 4^k matrices in the {I,X,Y,Z}^⊗k operator basis, NOT the
/// computational-basis T ⊗ I − I ⊗ Tᵀ which lives in a different basis and would silently break the check).
/// These tests pin <see cref="KBodyPalindromeRouting.Routes"/> against the spectral authority
/// <see cref="PauliPairTrichotomy.Classify(ChainSystem, IReadOnlyList{PauliTerm}, double, double, PauliLetter)"/>
/// on the 15 ground-truth cases from the verified derivation (Python, 2026-06-06):
///
/// <list type="bullet">
///   <item>the 8 discrete-routable soft sets (Routes == true AND spectral Soft);</item>
///   <item>the 4 NON-LOCAL ceiling cases (soft, NO per-site product Q, so Routes == false);</item>
///   <item>the 2 LOCAL continuous-sum cases (soft, a per-site product Q DOES exist but routes via
///     continuous-sum not per-term, so the per-term Routes still declines them, the coverage gap);</item>
///   <item>the hard XXX+XXY+YXX (no Q, spectral Hard, Routes == false).</item>
/// </list>
///
/// <para>If <see cref="KBodyPalindromeRouting.Routes"/> disagrees with the authority on ANY of these, the
/// Pauli-basis port has a bug (a convention error in Q_k, the commutator superoperator, or their shared
/// basis). The witnesses are the ground truth and are NOT to be adjusted to a buggy check.</para></summary>
public class KBodyPalindromeRoutingTests
{
    private static PauliTerm T(string label) => new(PauliLabel.Parse(label), Complex.One);
    private static List<PauliTerm> H(params string[] labels) => labels.Select(T).ToList();

    /// <summary>The spectral authority verdict at N = 4 under Z-dephasing (J = 1.0, γ₀ = 0.05).</summary>
    private static TrichotomyClass Spectral(IReadOnlyList<PauliTerm> terms) =>
        PauliPairTrichotomy.Classify(
            new ChainSystem(N: 4, J: 1.0, GammaZero: 0.05), terms, dephaseLetter: PauliLetter.Z);

    /// <summary>The 8 discrete-routable soft sets from the verified derivation: a per-site product Q (P1, P4,
    /// M2, or a period-2 mix) palindromizes every term, so Routes == true, and they are spectrally Soft.</summary>
    public static IEnumerable<object[]> RoutableSoft() => new[]
    {
        new object[] { new[] { "XIX", "XXY", "YXX" } },
        new object[] { new[] { "YIY", "XXY", "YXX" } },
        new object[] { new[] { "IYI", "XZY", "YZX" } },
        new object[] { new[] { "XIX", "IXZ", "ZXI" } },
        new object[] { new[] { "XYX", "XZY", "YZX" } },
        new object[] { new[] { "YIY", "IXZ", "ZXI" } },
        new object[] { new[] { "YYY", "XZY", "YZX" } },
        new object[] { new[] { "ZYZ", "XZY", "YZX" } },
    };

    /// <summary>The 4 NON-LOCAL ceiling cases: spectrally Soft, but NO per-site product Q palindromizes them
    /// (the mirror is genuinely non-local), so Routes == false. These stay the certifier's ceiling.</summary>
    public static IEnumerable<object[]> NonLocalCeiling() => new[]
    {
        new object[] { new[] { "XZX", "XZY", "YZX" } },
        new object[] { new[] { "YZY", "XZY", "YZX" } },
        new object[] { new[] { "IXI", "IIY", "YII" } },
        new object[] { new[] { "IYI", "IIX", "XII" } },
    };

    /// <summary>The 2 LOCAL cases the per-term router cannot reach: spectrally Soft AND a per-site product Q
    /// DOES palindromize them (a continuous-uniform M, residual ~1e-13, verified N=3,4,5,
    /// simulations/ceiling_6to4_verification.py). But that M routes via continuous-SUM cancellation, which the
    /// per-term {Q_k,[T,.]_k}=0 check does not see, AND it is an arbitrary continuous map (not a candidate), so
    /// Routes == false here too: the honest coverage gap, NOT non-locality. Not part of the 4-case ceiling.</summary>
    public static IEnumerable<object[]> LocalButNotPerTermRoutable() => new[]
    {
        new object[] { new[] { "XIX", "XIY", "YIX" } },
        new object[] { new[] { "YIY", "XIY", "YIX" } },
    };

    [Theory]
    [MemberData(nameof(RoutableSoft))]
    public void Routes_True_AndSpectralSoft_ForTheEightRoutableSets(string[] labels)
    {
        var terms = H(labels);
        Assert.True(KBodyPalindromeRouting.Routes(terms, n: 4),
            $"expected Routes == true (a per-site Q exists) for {string.Join("+", labels)}");
        Assert.Equal(TrichotomyClass.Soft, Spectral(terms));
    }

    [Theory]
    [MemberData(nameof(NonLocalCeiling))]
    public void Routes_False_ButSpectralSoft_ForTheFourNonLocalCeilingSets(string[] labels)
    {
        var terms = H(labels);
        Assert.False(KBodyPalindromeRouting.Routes(terms, n: 4),
            $"expected Routes == false (no per-site Q at all) for the non-local ceiling {string.Join("+", labels)}");
        Assert.Equal(TrichotomyClass.Soft, Spectral(terms));   // NotCertified does not imply not-soft
    }

    [Theory]
    [MemberData(nameof(LocalButNotPerTermRoutable))]
    public void Routes_False_ButLocal_ForTheTwoContinuousSumCases(string[] labels)
    {
        // These ARE local (a continuous-uniform per-site product Q exists, verified). The per-term Routes
        // declines them only because that Q routes via continuous-SUM cancellation, not per term: the gap.
        var terms = H(labels);
        Assert.False(KBodyPalindromeRouting.Routes(terms, n: 4),
            $"per-term Routes declines the continuous-sum case {string.Join("+", labels)} (the coverage gap)");
        Assert.Equal(TrichotomyClass.Soft, Spectral(terms));
    }

    [Fact]
    public void Routes_False_AndSpectralHard_ForXXXPlusXXYPlusYXX()
    {
        var terms = H("XXX", "XXY", "YXX");
        Assert.False(KBodyPalindromeRouting.Routes(terms, n: 4));
        Assert.Equal(TrichotomyClass.Hard, Spectral(terms));
    }

    /// <summary>The single load-bearing assertion: <see cref="KBodyPalindromeRouting.Routes"/> never
    /// disagrees with the spectral authority across all 15 ground-truth witnesses. A Routes == true must be
    /// Soft (constructive soundness), and the 7 Routes == false witnesses are exactly the 4 non-local ceiling
    /// cases (still Soft) plus the 2 local continuous-sum cases (still Soft) plus the hard set; the certifier
    /// is one-sided, so a Routes == false carries no claim about the spectral verdict beyond what the witness
    /// list records.</summary>
    [Fact]
    public void Routes_MatchesTheSpectralAuthority_OnAll15GroundTruthWitnesses()
    {
        var mismatches = new List<string>();

        foreach (var row in RoutableSoft())
        {
            var labels = (string[])row[0];
            var terms = H(labels);
            bool routes = KBodyPalindromeRouting.Routes(terms, n: 4);
            var spectral = Spectral(terms);
            if (!routes || spectral != TrichotomyClass.Soft)
                mismatches.Add($"{string.Join("+", labels)}: Routes={routes} spectral={spectral} (want Routes=true, Soft)");
        }

        foreach (var row in NonLocalCeiling())
        {
            var labels = (string[])row[0];
            var terms = H(labels);
            bool routes = KBodyPalindromeRouting.Routes(terms, n: 4);
            var spectral = Spectral(terms);
            if (routes || spectral != TrichotomyClass.Soft)
                mismatches.Add($"{string.Join("+", labels)}: Routes={routes} spectral={spectral} (want Routes=false, Soft)");
        }

        foreach (var row in LocalButNotPerTermRoutable())
        {
            var labels = (string[])row[0];
            var terms = H(labels);
            bool routes = KBodyPalindromeRouting.Routes(terms, n: 4);
            var spectral = Spectral(terms);
            if (routes || spectral != TrichotomyClass.Soft)
                mismatches.Add($"{string.Join("+", labels)}: Routes={routes} spectral={spectral} (want Routes=false, Soft, local)");
        }

        {
            var labels = new[] { "XXX", "XXY", "YXX" };
            var terms = H(labels);
            bool routes = KBodyPalindromeRouting.Routes(terms, n: 4);
            var spectral = Spectral(terms);
            if (routes || spectral != TrichotomyClass.Hard)
                mismatches.Add($"{string.Join("+", labels)}: Routes={routes} spectral={spectral} (want Routes=false, Hard)");
        }

        // A Routes == true must imply not-Hard (constructive soundness), the core safety property.
        Assert.True(mismatches.Count == 0,
            $"{mismatches.Count} witness(es) disagree with the spectral authority: {string.Join("; ", mismatches)}");
    }
}

/// <summary>First DIRECT coverage of the public surface of <see cref="KBodyPalindromeRouting"/>. The
/// convention-guard suite above funnels everything through <see cref="KBodyPalindromeRouting.Routes"/>;
/// these facts exercise the building blocks individually: the candidate-set enumeration
/// (<see cref="KBodyPalindromeRouting.CandidateSet"/>), the per-term anticommutator
/// (<see cref="KBodyPalindromeRouting.PerTermAnticommutes"/>, which drives BuildQk + the commutator
/// superoperator) at a verified-true point, and the reporting helper
/// (<see cref="KBodyPalindromeRouting.RoutingCandidate"/>) on both a routable set and the negative path.</summary>
public class KBodyPalindromeRoutingPublicApiTests
{
    private static PauliTerm T(string label) => new(PauliLabel.Parse(label), Complex.One);
    private static List<PauliTerm> H(params string[] labels) => labels.Select(T).ToList();

    [Fact]
    public void CandidateSet_Counts_Are_FourToThePeriod_PerPeriod()
    {
        // Period 1 only: 4^1 = 4 single-site patterns.
        Assert.Equal(4, KBodyPalindromeRouting.CandidateSet(maxPeriod: 1).Count);
        // Periods 1 and 2: 4^1 + 4^2 = 4 + 16 = 20.
        Assert.Equal(20, KBodyPalindromeRouting.CandidateSet(maxPeriod: 2).Count);
    }

    [Fact]
    public void PerTermAnticommutes_True_ForEveryTemplateAndOffset_UnderTheRoutingCandidate()
    {
        var terms = H("XIX", "XXY", "YXX");

        // The reporting helper must exhibit a certificate (non-null) for this routable set.
        string? description = KBodyPalindromeRouting.RoutingCandidate(terms, n: 4);
        Assert.NotNull(description);

        // Recover the actual candidate (Pattern + Period) the helper named, by matching its Description.
        var candidate = KBodyPalindromeRouting.CandidateSet()
            .Single(c => c.Description == description);

        // The per-term anticommutator must vanish for each template at every window-offset under that pattern.
        foreach (var term in terms)
            for (int offset = 0; offset < candidate.Period; offset++)
                Assert.True(
                    KBodyPalindromeRouting.PerTermAnticommutes(term, offset, candidate.Pattern, candidate.Period),
                    $"PerTermAnticommutes expected true for {term.Label} at offset {offset} under {description}");
    }

    [Fact]
    public void RoutingCandidate_NonNull_ForRoutableSet()
    {
        // XIX+XXY+YXX is the canonical discrete-routable soft set: a per-site Q must be exhibited.
        Assert.NotNull(KBodyPalindromeRouting.RoutingCandidate(H("XIX", "XXY", "YXX"), n: 4));
    }

    [Fact]
    public void RoutingCandidate_Null_ForNonLocalCeiling()
    {
        // XZX+XZY+YZX is spectrally soft but admits no per-site product Q (the mirror is non-local).
        Assert.Null(KBodyPalindromeRouting.RoutingCandidate(H("XZX", "XZY", "YZX"), n: 4));
    }

    [Fact]
    public void RoutingCandidate_Null_ForHardSet()
    {
        // XXX+XXY+YXX is spectrally hard: no palindromizing Q of any kind.
        Assert.Null(KBodyPalindromeRouting.RoutingCandidate(H("XXX", "XXY", "YXX"), n: 4));
    }
}
