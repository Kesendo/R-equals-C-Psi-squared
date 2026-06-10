using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
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
///   <item>the 2 Z-middle golden-routed cases (soft AND per-site routable, by the period-4 golden
///     window-summed router, F116; the per-term Routes correctly returns false because the cancellation
///     is cross-template inside one window, the documented coverage gap of the per-term lens);</item>
///   <item>the 2 LOCAL continuous-sum cases (soft, a per-site product Q DOES exist but routes via
///     continuous-sum not per-term, so the per-term Routes still declines them, the coverage gap);</item>
///   <item>the 2 LOCAL single-site-field cases (the I-heavy: soft, certified by SingleSiteField, but the
///     per-term Routes declines their single-site router, so Routes == false; the 4 to 2 step);</item>
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

    /// <summary>The 2 Z-middle cases: spectrally Soft AND LOCAL after all (the period-4 golden router
    /// palindromizes them under the window-summed condition, docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md +
    /// F116; certified by the RoutingWindowSummed strategy). The old verdict "the mirror is genuinely
    /// non-local" is overturned; the per-term Routes == false stays CORRECT (each template's anticommutator
    /// alone is nonzero, the cancellation is cross-template inside one window), the documented coverage
    /// gap of the per-term lens.</summary>
    public static IEnumerable<object[]> LocalViaGoldenWindowRouting() => new[]
    {
        new object[] { new[] { "XZX", "XZY", "YZX" } },
        new object[] { new[] { "YZY", "XZY", "YZX" } },
    };

    /// <summary>The 2 I-heavy cases: spectrally Soft AND LOCAL (a site-varying per-site product Q palindromizes
    /// them, a product of single-site crossover maps, constructive, machine-precision N=4,5,6,
    /// simulations/ceiling_4to2_iheavy_local.py). PalindromeSoftCertifier.Certify certifies them via the
    /// SingleSiteField strategy. But that router is single-site (not a discrete period≤2 per-term pattern), so
    /// the per-term KBodyPalindromeRouting.Routes does not see it and returns false: the honest coverage gap of
    /// the per-term router, NOT non-locality. The 4 to 2 correction; they are no longer in the ceiling.</summary>
    public static IEnumerable<object[]> LocalViaSingleSiteFields() => new[]
    {
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
    [MemberData(nameof(LocalViaGoldenWindowRouting))]
    public void Routes_False_PerTermLens_ButSpectralSoft_ForTheTwoGoldenRoutedSets(string[] labels)
    {
        // A PER-TERM statement only: each template's anticommutator alone is nonzero, so the per-term
        // Routes correctly declines these. They ARE per-site routable (the golden window-summed router,
        // F116), certified by RoutesWindowSummed: see KBodyPalindromeRoutingWindowSummedTests below.
        var terms = H(labels);
        Assert.False(KBodyPalindromeRouting.Routes(terms, n: 4),
            $"per-term Routes declines the golden-routed Z-middle case {string.Join("+", labels)} (the per-term coverage gap)");
        Assert.Equal(TrichotomyClass.Soft, Spectral(terms));   // NotCertified-per-term does not imply not-soft
    }

    [Theory]
    [MemberData(nameof(LocalViaSingleSiteFields))]
    public void Routes_False_ButLocal_ForTheTwoSingleSiteFieldCases(string[] labels)
    {
        // These ARE local (a site-varying per-site product Q exists; PalindromeSoftCertifier certifies them via
        // SingleSiteField). The per-term Routes declines them only because the construction is single-site, not
        // a discrete period≤2 per-term pattern: the coverage gap, NOT non-locality.
        var terms = H(labels);
        Assert.False(KBodyPalindromeRouting.Routes(terms, n: 4),
            $"per-term Routes declines the single-site-field case {string.Join("+", labels)} (the coverage gap)");
        Assert.Equal(TrichotomyClass.Soft, Spectral(terms));
        Assert.Equal(PalindromeSoftCertifier.SoftStrategy.SingleSiteField,
            PalindromeSoftCertifier.Certify(terms, n: 4).Strategy);   // but Certify DOES catch them
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
    /// Soft (constructive soundness), and the 7 Routes == false witnesses are exactly the 2 Z-middle
    /// golden-routed cases (still Soft, certified by RoutingWindowSummed, F116), the 2 local
    /// single-site-field cases (still Soft, certified by SingleSiteField),
    /// the 2 local continuous-sum cases (still Soft), and the hard set; the certifier is one-sided, so a
    /// Routes == false carries no claim about the spectral verdict beyond what the witness list records.</summary>
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

        foreach (var row in LocalViaGoldenWindowRouting())
        {
            var labels = (string[])row[0];
            var terms = H(labels);
            bool routes = KBodyPalindromeRouting.Routes(terms, n: 4);
            var spectral = Spectral(terms);
            if (routes || spectral != TrichotomyClass.Soft)
                mismatches.Add($"{string.Join("+", labels)}: Routes={routes} spectral={spectral} (want Routes=false, Soft, golden-window-routed)");
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

        foreach (var row in LocalViaSingleSiteFields())
        {
            var labels = (string[])row[0];
            var terms = H(labels);
            bool routes = KBodyPalindromeRouting.Routes(terms, n: 4);
            var spectral = Spectral(terms);
            if (routes || spectral != TrichotomyClass.Soft)
                mismatches.Add($"{string.Join("+", labels)}: Routes={routes} spectral={spectral} (want Routes=false, Soft, local-single-site)");
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
    public void RoutingCandidate_Null_PerTermLens_ForTheZMiddleGoldenCase()
    {
        // XZX+XZY+YZX is spectrally soft and golden-window-routable (F116), but no PER-TERM candidate
        // routes it (the cancellation is cross-template inside one window), so the per-term reporting
        // helper correctly returns null; RoutesWindowSummed exhibits the golden certificate instead.
        Assert.Null(KBodyPalindromeRouting.RoutingCandidate(H("XZX", "XZY", "YZX"), n: 4));
    }

    [Fact]
    public void RoutingCandidate_Null_ForHardSet()
    {
        // XXX+XXY+YXX is spectrally hard: no palindromizing Q of any kind.
        Assert.Null(KBodyPalindromeRouting.RoutingCandidate(H("XXX", "XXY", "YXX"), n: 4));
    }
}

/// <summary>The WINDOW-SUMMED routing primitive (Stufe B′, F116,
/// docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md): the golden period-4 candidates
/// (<see cref="KBodyPalindromeRouting.GoldenSiteMaps"/> and the X↔Y conjugate
/// <see cref="KBodyPalindromeRouting.GoldenMirrorSiteMaps"/>), the cross-template window lemma (the
/// template-summed anticommutator vanishes at all four offsets while every per-term check fails), the
/// per-site invariants (class-swap and q² = −(2+φ)·I, the two halves of the routing derivation), and the
/// gates (hard sets and mixed-span sets are declined). The candidate values are pinned by the proof and
/// the exact-ring anchor <c>simulations/ceiling_golden_router.py</c>; these tests are the C#-port guard
/// (a sign error in any h entry flips the window lemma).</summary>
public class KBodyPalindromeRoutingWindowSummedTests
{
    private static PauliTerm T(string label) => new(PauliLabel.Parse(label), Complex.One);
    private static List<PauliTerm> H(params string[] labels) => labels.Select(T).ToList();

    [Fact]
    public void RoutesWindowSummed_ReturnsTheGoldenCandidate_ForTheXzxSet()
    {
        // The Z-middle XZX+XZY+YZX routes via the golden [a, a, b, b] pattern (a = φX + Y, b = X − φY).
        Assert.Equal("Golden[a,a,b,b] (P=4)",
            KBodyPalindromeRouting.RoutesWindowSummed(H("XZX", "XZY", "YZX"), n: 4));
    }

    [Fact]
    public void RoutesWindowSummed_ReturnsTheMirrorCandidate_ForTheYzySet()
    {
        // The X↔Y sibling YZY+XZY+YZX routes via the conjugated maps q′ = s·q·s: the X↔Y mirror is not a
        // self-equivalence of the golden router, it maps one case's routers to the other's (proof §5).
        Assert.Equal("Golden-mirror[a,a,b,b] (P=4)",
            KBodyPalindromeRouting.RoutesWindowSummed(H("YZY", "XZY", "YZX"), n: 4));
    }

    [Fact]
    public void WindowSummedAnticommutator_VanishesAtAllFourOffsets_WhileEveryPerTermCheckFails()
    {
        // The window lemma (the mechanism, proof §2): {Q_3, [XZX+XZY+YZX,·]_3} = 0 at every window offset
        // 0..3, while no single template's anticommutator vanishes (a single-template list degenerates to
        // the per-term check under the same maps). The cancellation is CROSS-TEMPLATE inside one window:
        // exactly why the per-term lens (Routes / RoutingCandidate) could never see this router.
        var terms = H("XZX", "XZY", "YZX");
        var maps = KBodyPalindromeRouting.GoldenSiteMaps;
        for (int offset = 0; offset < 4; offset++)
        {
            Assert.True(KBodyPalindromeRouting.PerWindowSummedAnticommutes(terms, offset, maps, period: 4),
                $"the template-summed anticommutator must vanish at offset {offset}");
            foreach (var term in terms)
                Assert.False(
                    KBodyPalindromeRouting.PerWindowSummedAnticommutes(new[] { term }, offset, maps, period: 4),
                    $"the per-term anticommutator for {term.Label} must NOT vanish at offset {offset}");
        }
    }

    [Fact]
    public void GoldenMaps_AreClassSwapping_WithQSquaredMinusTwoPlusPhiIdentity()
    {
        // The two per-site invariants behind the certificate, for BOTH candidate sets: class-swap (the
        // {I,Z}→{I,Z} and {X,Y}→{X,Y} blocks are zero, so the dissipator leg {W, D̂} = −2σW is automatic)
        // and q² = −(2+φ)·I (a scalar times a unitary, so W = ⊗q is invertible with condition number 1).
        double phi = (1.0 + System.Math.Sqrt(5.0)) / 2.0;
        var minusTwoPlusPhi = Matrix<Complex>.Build.DenseIdentity(4) * new Complex(-(2.0 + phi), 0.0);
        foreach (var maps in new[] { KBodyPalindromeRouting.GoldenSiteMaps, KBodyPalindromeRouting.GoldenMirrorSiteMaps })
        {
            Assert.Equal(4, maps.Count);
            foreach (var q in maps)
            {
                foreach (int r in new[] { 0, 3 })          // {I,Z} rows take nothing from {I,Z} columns
                    foreach (int c in new[] { 0, 3 })
                        Assert.Equal(Complex.Zero, q[r, c]);
                foreach (int r in new[] { 1, 2 })          // {X,Y} rows take nothing from {X,Y} columns
                    foreach (int c in new[] { 1, 2 })
                        Assert.Equal(Complex.Zero, q[r, c]);
                Assert.True((q * q - minusTwoPlusPhi).FrobeniusNorm() < 1e-12,
                    "q·q must equal −(2+φ)·Identity");
            }
        }
    }

    [Fact]
    public void RoutesWindowSummed_Null_ForAHardSet_AndForAMixedSpanSet()
    {
        // A hard set must not be window-summed-certified: a vanishing summed anticommutator plus
        // class-swap would PROVE the palindrome, so it cannot vanish on a spectrally hard set.
        Assert.Null(KBodyPalindromeRouting.RoutesWindowSummed(H("XXX", "XXY", "YXX"), n: 4));
        // Mixed spans have no shared window space: the gate declines them.
        Assert.Null(KBodyPalindromeRouting.RoutesWindowSummed(H("XX", "XZX"), n: 4));
    }
}
