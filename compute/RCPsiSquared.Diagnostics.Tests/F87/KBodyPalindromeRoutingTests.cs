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
///   <item>the 6 non-local ceiling cases (soft, but no per-site product Q exists, so Routes == false);</item>
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

    /// <summary>The 6 non-local ceiling cases: spectrally Soft, but NO per-site product Q palindromizes them
    /// (the mirror is non-local), so Routes == false. These stay the certifier's ceiling.</summary>
    public static IEnumerable<object[]> NonLocalCeiling() => new[]
    {
        new object[] { new[] { "XZX", "XZY", "YZX" } },
        new object[] { new[] { "XIX", "XIY", "YIX" } },
        new object[] { new[] { "YIY", "XIY", "YIX" } },
        new object[] { new[] { "YZY", "XZY", "YZX" } },
        new object[] { new[] { "IXI", "IIY", "YII" } },
        new object[] { new[] { "IYI", "IIX", "XII" } },
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
    public void Routes_False_ButSpectralSoft_ForTheSixNonLocalCeilingSets(string[] labels)
    {
        var terms = H(labels);
        Assert.False(KBodyPalindromeRouting.Routes(terms, n: 4),
            $"expected Routes == false (no per-site Q) for the non-local ceiling {string.Join("+", labels)}");
        // The spectral authority still calls these Soft: NotCertified does not imply not-soft.
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
    /// Soft (constructive soundness), and the 7 Routes == false witnesses are exactly the non-local ceiling
    /// (still Soft) plus the hard set; the certifier is one-sided, so a Routes == false carries no claim
    /// about the spectral verdict beyond what the witness list records.</summary>
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
