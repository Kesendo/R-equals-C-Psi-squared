using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F87;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>The Liouvillian-free two-term palindrome router (<see cref="TwoTermPalindromeRouting"/>)
/// reproduces both the spectral fate (truly / soft / hard) AND the hidden-symmetry Q family
/// (P1 / uniform / alternating / continuous / none) from the two bond bilinears' letters alone,
/// bit-exactly across all 36 unordered two-term combos.
///
/// <para>Ported from <c>simulations/framework/diagnostics/q_family_routing.py</c>
/// (<c>classify_two_term_palindrome</c>), whose docstring guarantees a 0-mismatch reproduction of
/// <c>fw.classify_pauli_pair</c> over the 36 combos at N = 3, 4, 5. The Python is the read-only
/// reference; these tests pin the C# port against (a) the Python's fate+family verdict table and
/// (b) the spectral authority <see cref="PauliPairTrichotomy.Classify(ChainSystem, IReadOnlyList{PauliTerm}, double, double, PauliLetter)"/>.</para></summary>
public class TwoTermPalindromeRoutingTests
{
    private static PauliTerm T(string label) => new(PauliLabel.Parse(label), Complex.One);

    /// <summary>The nine canonical two-letter bilinears (the only letter pairs with both sites lit).
    /// The 36 unordered combos are the upper triangle (a, b) with b at or after a.</summary>
    private static readonly string[] Bilinears = { "XX", "YY", "ZZ", "XY", "YX", "XZ", "ZX", "YZ", "ZY" };

    /// <summary>Expected (fate, family) verdict for each unordered pair, generated once from the Python
    /// reference <c>classify_two_term_palindrome(a, b, N=3)</c>. Keyed by "a+b" in Bilinears order.
    /// Family strings: "P1" (truly), "uniform" / "alternating" / "continuous" (soft), "None" (hard),
    /// matching the Python lowercase q_family tokens.</summary>
    private static readonly Dictionary<string, (string Fate, string Family)> Expected = new()
    {
        ["XX+XX"] = ("truly", "P1"),
        ["XX+YY"] = ("truly", "P1"),
        ["XX+ZZ"] = ("truly", "P1"),
        ["XX+XY"] = ("hard", "None"),
        ["XX+YX"] = ("hard", "None"),
        ["XX+XZ"] = ("soft", "uniform"),
        ["XX+ZX"] = ("soft", "uniform"),
        ["XX+YZ"] = ("soft", "uniform"),
        ["XX+ZY"] = ("soft", "uniform"),
        ["YY+YY"] = ("truly", "P1"),
        ["YY+ZZ"] = ("truly", "P1"),
        ["YY+XY"] = ("hard", "None"),
        ["YY+YX"] = ("hard", "None"),
        ["YY+XZ"] = ("soft", "uniform"),
        ["YY+ZX"] = ("soft", "uniform"),
        ["YY+YZ"] = ("soft", "uniform"),
        ["YY+ZY"] = ("soft", "uniform"),
        ["ZZ+ZZ"] = ("truly", "P1"),
        ["ZZ+XY"] = ("soft", "alternating"),
        ["ZZ+YX"] = ("soft", "alternating"),
        ["ZZ+XZ"] = ("soft", "uniform"),
        ["ZZ+ZX"] = ("soft", "uniform"),
        ["ZZ+YZ"] = ("soft", "uniform"),
        ["ZZ+ZY"] = ("soft", "uniform"),
        ["XY+XY"] = ("soft", "alternating"),
        ["XY+YX"] = ("soft", "alternating"),
        ["XY+XZ"] = ("hard", "None"),
        ["XY+ZX"] = ("hard", "None"),
        ["XY+YZ"] = ("hard", "None"),
        ["XY+ZY"] = ("hard", "None"),
        ["YX+YX"] = ("soft", "alternating"),
        ["YX+XZ"] = ("hard", "None"),
        ["YX+ZX"] = ("hard", "None"),
        ["YX+YZ"] = ("hard", "None"),
        ["YX+ZY"] = ("hard", "None"),
        ["XZ+XZ"] = ("soft", "uniform"),
        ["XZ+ZX"] = ("soft", "uniform"),
        ["XZ+YZ"] = ("soft", "continuous"),
        ["XZ+ZY"] = ("hard", "None"),
        ["ZX+ZX"] = ("soft", "uniform"),
        ["ZX+YZ"] = ("hard", "None"),
        ["ZX+ZY"] = ("soft", "continuous"),
        ["YZ+YZ"] = ("soft", "uniform"),
        ["YZ+ZY"] = ("soft", "uniform"),
        ["ZY+ZY"] = ("soft", "uniform"),
    };

    /// <summary>Map the C# (Fate, Family) result back to the Python lowercase tokens for comparison.</summary>
    private static string FateToken(TrichotomyClass fate) => fate switch
    {
        TrichotomyClass.Truly => "truly",
        TrichotomyClass.Soft => "soft",
        TrichotomyClass.Hard => "hard",
        _ => "?",
    };

    private static string FamilyToken(QFamily family) => family switch
    {
        QFamily.P1 => "P1",
        QFamily.Uniform => "uniform",
        QFamily.Alternating => "alternating",
        QFamily.Continuous => "continuous",
        QFamily.None_ => "None",
        _ => "?",
    };

    [Fact]
    public void Routing_ReproducesPythonFateAndFamily_OverAll36Combos()
    {
        var mismatches = new List<string>();
        for (int i = 0; i < Bilinears.Length; i++)
        {
            for (int j = i; j < Bilinears.Length; j++)
            {
                string a = Bilinears[i], b = Bilinears[j];
                var r = TwoTermPalindromeRouting.Classify(a, b);
                var (expFate, expFamily) = Expected[$"{a}+{b}"];
                string gotFate = FateToken(r.Fate);
                string gotFamily = FamilyToken(r.Family);
                if (gotFate != expFate || gotFamily != expFamily)
                    mismatches.Add($"{a}+{b}: got ({gotFate},{gotFamily}) expected ({expFate},{expFamily})");
            }
        }
        Assert.True(mismatches.Count == 0,
            $"{mismatches.Count} routing mismatch(es): {string.Join("; ", mismatches)}");
    }

    [Theory]
    [InlineData("XX", "YY", TrichotomyClass.Truly, QFamily.P1)]          // both Mother: canonical Π = P1
    [InlineData("XX", "XZ", TrichotomyClass.Soft, QFamily.Uniform)]      // shared uniform Q-family {P4}
    [InlineData("XY", "YX", TrichotomyClass.Soft, QFamily.Alternating)]  // XY/YX site-parity Q
    [InlineData("XZ", "YZ", TrichotomyClass.Soft, QFamily.Continuous)]   // same-site X&Y over dark Z
    [InlineData("XX", "XY", TrichotomyClass.Hard, QFamily.None_)]        // no Q closes it
    public void Routing_PinsTheCanonicalAnchors(string a, string b, TrichotomyClass fate, QFamily family)
    {
        var r = TwoTermPalindromeRouting.Classify(a, b);
        Assert.Equal(fate, r.Fate);
        Assert.Equal(family, r.Family);
        // Symmetry: the router is order-independent.
        var swapped = TwoTermPalindromeRouting.Classify(b, a);
        Assert.Equal(fate, swapped.Fate);
        Assert.Equal(family, swapped.Family);
    }

    [Fact]
    public void Routing_KleinIndicesMatchTheKleinConvention()
    {
        // Klein index (bit_a, bit_b) per term: XX/YY/ZZ -> (0,0); XY/YX -> (0,1); YZ/ZY -> (1,0); XZ/ZX -> (1,1).
        Assert.Equal((0, 0), TwoTermPalindromeRouting.Classify("XX", "YY").Klein1);
        Assert.Equal((0, 1), TwoTermPalindromeRouting.Classify("XY", "YX").Klein1);
        Assert.Equal((1, 0), TwoTermPalindromeRouting.Classify("YZ", "ZY").Klein1);
        Assert.Equal((1, 1), TwoTermPalindromeRouting.Classify("XZ", "ZX").Klein1);
    }

    [Fact]
    public void Routing_ParityBreak_FlagsBrokenGlobalSpinFlipSymmetry()
    {
        // XX+YY (the XY model) commutes with X^{⊗N}: no parity break.
        Assert.False(TwoTermPalindromeRouting.Classify("XX", "YY", n: 3).ParityBreak);
        // XZ carries a lone Z on one site: [Σ XZ, X^{⊗N}] != 0, parity broken.
        Assert.True(TwoTermPalindromeRouting.Classify("XX", "XZ", n: 3).ParityBreak);
    }

    [Theory]
    [InlineData("AB", "XX")]  // unknown letter 'A'/'B'
    [InlineData("X", "YY")]   // one-letter term
    [InlineData("XXX", "YY")] // three-letter term
    public void Routing_RejectsMalformedTerms(string a, string b)
    {
        Assert.ThrowsAny<System.ArgumentException>(() => TwoTermPalindromeRouting.Classify(a, b));
    }

    // -------------------------------------------------------------------------------------------
    // The substance: bit-exact agreement with the spectral authority over all 36 combos at N = 4.
    // -------------------------------------------------------------------------------------------

    [Fact]
    public void Routing_AgreesWithSpectralAuthority_OverAll36Combos_N4()
    {
        var chain = new ChainSystem(4, 1.0, 0.05);
        var mismatches = new List<string>();
        for (int i = 0; i < Bilinears.Length; i++)
        {
            for (int j = i; j < Bilinears.Length; j++)
            {
                string a = Bilinears[i], b = Bilinears[j];
                TrichotomyClass routed = TwoTermPalindromeRouting.Classify(a, b).Fate;
                TrichotomyClass authority = PauliPairTrichotomy.Classify(
                    chain, new[] { T(a), T(b) });
                if (routed != authority)
                    mismatches.Add($"{a}+{b}: routed={routed} authority={authority}");
            }
        }
        Assert.True(mismatches.Count == 0,
            $"{mismatches.Count} mismatch(es) vs spectral authority: {string.Join("; ", mismatches)}");
    }
}
