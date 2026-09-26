using System;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The F89 Δ response of the (SE,DE) diabolics, certified. The (q,Δ) XxzCoherenceBlock is validated
/// against the F89 block at Δ=0; every tracker proposal is refined by the discriminant Newton on the full R-even
/// block and certified there: coincidence and seed correspondence, an isolated algebraic-2 pair, the character
/// read as geometric against algebraic multiplicity on the Riesz compression (R-3), and the order of the pair
/// discriminant's zero. At Δ=0 every sampled seed is a double zero read geo = alg = 2 (diabolic); under Δ &gt; 0 the
/// double zero splits into simple zeros, and every one the certificate claims (inside the seed's isolation disk)
/// is an EP2 read alg = 2, geo = 1. The N=5 defective control is a simple zero at Δ=0 and stays one.
///
/// <para>Reference values: the N=4 EP2 locations are bracketed to 1e-15 by a sign change of the real discriminant
/// in 40-digit arithmetic (simulations/tests/test_f89_zz_break_surface.py); the N=5..7 values come from an independent numpy port of BuildFull/BuildSym with the
/// same discriminant Newton and Riesz compression, cross-checked against the popcount-(1,2) block of the full 4^N
/// Liouvillian at N = 3, 4, 5. Location gates use the located-zero error model: an EP2 is fixed to
/// (g/s)² + u·|q| ≈ 1e-15 by its residual gap g ≈ 1e-8 and square-root slope s, a crossing to about g/s ≈ 1e-11 by
/// its linear slope; 1e-9 leaves a margin of 100 and more over both. Departure gate: the departure is computed as
/// √(‖A‖_F² − Σ|λ|²), so its rounding is k·u·‖A‖_F²/(2·dep) with k ≤ 17 (the factor measured for the same
/// cancellation in XxzCoherenceBlock.DepartureRoundingFloor), 2e-10 at the smallest sampled departure (2e-4 at
/// ‖A‖_F ≈ 6.4); the stored references add ≤ 5e-11 (ten decimals). 2e-9 sits 8 times above the sum and 7 times
/// below the closest pair of sampled departures (the two N=7 Δ = 1e-4 EP2s, 1.5e-8 apart).</para>
/// See docs/superpowers/plans/2026-06-27-f89-path4-delta-test.md.</summary>
public class XxzDeltaFlipTests
{
    private const double LocationGate = 1e-9;
    private const double DepartureGate = 2e-9;

    // (q, departure) of the reference EP2s; the Δ=0 entries are the diabolic crossings (departure 0).
    private static readonly (Complex Q, double Dep)[] N5Clean002 =
        { (new(0.638693461456209, 0.179553415463510), 0.0164610229), (new(0.642836111498104, 0.182468939434062), 0.0181502035) };
    private static readonly (Complex Q, double Dep)[] N5Clean005 =
        { (new(0.636122282920920, 0.179632001098883), 0.0402192599), (new(0.646100305976321, 0.187325655438013), 0.0450382591) };
    private static readonly (Complex Q, double Dep)[] N5NearReal002 =
        { (new(0.766885329309701, 0.026744591009047), 0.0323030377), (new(0.763972315020105, 0.020834606964548), 0.0423583888) };
    private static readonly (Complex Q, double Dep)[] N5NearReal005 =
        { (new(0.769255182403741, 0.030528760318336), 0.0829661969), (new(0.762303221483473, 0.015742622799481), 0.1051665758) };
    private static readonly (Complex Q, double Dep)[] N5Far002 =
        { (new(1.913918824946734, 1.204021447638601), 0.0541287039), (new(2.075880940635436, 1.303213838169880), 0.1514407086) };
    private static readonly (Complex Q, double Dep)[] N5Far005 =
        { (new(1.857702701098991, 1.212555379817808), 0.1166855278), (new(2.216569087055637, 1.333722902437400), 0.3839495688) };
    private static readonly (Complex Q, double Dep)[] N6RungNear002 =
        { (new(0.713172895438829, -0.218687095442581), 0.0282712386), (new(0.704443899333969, -0.211390927262037), 0.0315420326) };
    private static readonly (Complex Q, double Dep)[] N6RungNear01 =
        { (new(0.737241358192659, -0.218415884470263), 0.1570239390), (new(0.685635461519334, -0.188696187438170), 0.1597906956) };
    private static readonly (Complex Q, double Dep)[] N6RungFarA002 =
        { (new(0.759501578925560, 0.265399325546476), 0.0034587659), (new(0.766373386912108, 0.267474719592798), 0.0038068687) };
    private static readonly (Complex Q, double Dep)[] N6RungFarA01 =
        { (new(0.766649902434995, 0.283868231731283), 0.0177490800), (new(0.810051803754318, 0.300147157318646), 0.0207432524) };
    private static readonly (Complex Q, double Dep)[] N6RungFarB002 =
        { (new(1.062794886554044, 0.217082576814778), 0.0224952995), (new(1.071577824002448, 0.228859330406448), 0.0157751082) };
    private static readonly (Complex Q, double Dep)[] N6RungFarB01 =
        { (new(1.129251551545453, 0.192815022569929), 0.0642782601), (new(1.057871842141407, 0.190121782003818), 0.8455010551) };
    private static readonly (Complex Q, double Dep)[] N7Q11264At002 =
        { (new(1.127019196759119, -0.005838552257421), 0.0015491164), (new(1.124649022572506, -0.017072795701856), 0.0164956276) };
    private static readonly (Complex Q, double Dep)[] N7Q13038At002 =
        { (new(1.304563975756159, 0.011210474246574), 0.0019434409), (new(1.307107898662784, 0.024072934462780), 0.0034928013) };
    // Δ = 0.10 zero of the q = 1.3038 crossing reached from the seed; its λ lies 0.27 from the seed, beyond the
    // isolation radius 0.12, so the certificate declines it.
    private static readonly (Complex Q, double Dep)[] N7Q13038At01 =
        { (new(1.324862139742824, 0.055679503457439), 0.0136834698) };
    private static readonly (Complex Q, double Dep)[] N7Q26280At002 =
        { (new(2.305123211619291, -0.008674468298819), 0.0993297463), (new(3.063136491496817, -0.869706131304282), 0.1306314950) };
    private static readonly (Complex Q, double Dep)[] N7Q06788At00001 =
        { (new(0.678844523206904, 0.000003752503345), 0.0001976934), (new(0.678788335090143, 0.000003755920792), 0.0001977089) };
    private static readonly (Complex Q, double Dep)[] N7Q06788At002 =
        { (new(0.673252169224504, 0.000819191761567), 0.0398286617), (new(0.684501543679595, 0.000682419656153), 0.0392234259) };
    private static readonly (Complex Q, double Dep)[] N7Q06788At01 =
        { (new(0.709297010560962, 0.001995795465280), 0.1907998893), (new(0.651761333038970, 0.005396563597109), 0.2001046286) };
    // Δ = 0.10 zeros of the q = 1.1264 crossing (continuation); their λ lies 0.24 and 0.22 from the seed, beyond its
    // isolation radius 0.11, so the certificate reaches them and declines to claim them.
    private static readonly (Complex Q, double Dep)[] N7Q11264At01 =
        { (new(1.141313949402483, -0.028301984649874), 0.0091050921), (new(1.100804182618295, -0.063102057228388), 0.0772931712) };

    // A certified Jordan EP2 at one of the reference zeros: alg 2, geo 1, a simple zero of the pair discriminant.
    private static void AssertEp2AtOneOf(XxzCoherenceBlock.DeltaTrackResult r, (Complex Q, double Dep)[] refs, string label)
    {
        Console.WriteLine($"{label}: {r}");
        Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Defective, r.Verdict);
        Assert.Equal(2, r.Algebraic);
        Assert.Equal(1, r.Geometric);
        Assert.Equal(1, r.DiscriminantZeroOrder);
        Assert.False(r.Survived);
        Assert.True(r.Gap <= XxzCoherenceBlock.FullBlockCoincidenceTolerance, $"{label}: gap {r.Gap:E2}");
        var hit = refs.OrderBy(e => (e.Q - r.QCandidate).Magnitude).First();
        Assert.True((hit.Q - r.QCandidate).Magnitude <= LocationGate,
            $"{label}: q {r.QCandidate} is not a reference EP2 (nearest {hit.Q}, distance {(hit.Q - r.QCandidate).Magnitude:E2})");
        Assert.True(Math.Abs(r.Departure - hit.Dep) <= DepartureGate, $"{label}: departure {r.Departure:R} vs {hit.Dep:R}");
    }

    // A certified semisimple crossing: alg = geo = 2, a double zero of the pair discriminant.
    private static void AssertCrossingAt(XxzCoherenceBlock.DeltaTrackResult r, Complex qRef, string label)
    {
        Console.WriteLine($"{label}: {r}");
        Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Diabolic, r.Verdict);
        Assert.Equal(2, r.Algebraic);
        Assert.Equal(2, r.Geometric);
        Assert.Equal(2, r.DiscriminantZeroOrder);
        Assert.True(r.Survived);
        Assert.True(r.Gap <= XxzCoherenceBlock.FullBlockCoincidenceTolerance, $"{label}: gap {r.Gap:E2}");
        Assert.True((r.QCandidate - qRef).Magnitude <= LocationGate,
            $"{label}: q {r.QCandidate} vs the reference crossing {qRef} ({(r.QCandidate - qRef).Magnitude:E2})");
    }

    private static void AssertUncertified(XxzCoherenceBlock.DeltaTrackResult r)
    {
        Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Uncertified, r.Verdict);
        Assert.Equal(0, r.Algebraic);
        Assert.Equal(0, r.Geometric);
        Assert.Equal(0, r.DiscriminantZeroOrder);
        Assert.True(double.IsNaN(r.Departure));
        Assert.Null(r.Survived);
    }

    // For a proposal path whose starting point depends on the compressed roots: at Δ > 0 no crossing survives near
    // the seed, so the certificate is never Diabolic; what it issues is a certified EP2 (and then, if it sits at a
    // reference zero, with that zero's departure) or Uncertified when the zero it reaches has left the seed's disk.
    private static void AssertNoCrossingSurvives(XxzCoherenceBlock.DeltaTrackResult r, (Complex Q, double Dep)[] refs, string label)
    {
        Console.WriteLine($"{label}: {r}");
        Assert.NotEqual(XxzCoherenceBlock.DeltaFlipVerdict.Diabolic, r.Verdict);
        if (r.Verdict == XxzCoherenceBlock.DeltaFlipVerdict.Uncertified) { AssertUncertified(r); return; }
        Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Defective, r.Verdict);
        Assert.Equal(2, r.Algebraic);
        Assert.Equal(1, r.Geometric);
        Assert.Equal(1, r.DiscriminantZeroOrder);
        Assert.True(r.Gap <= XxzCoherenceBlock.FullBlockCoincidenceTolerance, $"{label}: gap {r.Gap:E2}");
        foreach (var e in refs.Where(e => (e.Q - r.QCandidate).Magnitude <= LocationGate))
            Assert.True(Math.Abs(r.Departure - e.Dep) <= DepartureGate, $"{label}: departure {r.Departure:R} vs {e.Dep:R}");
    }

    [Theory]
    [Trait("Category", "TASK8_DELTA_CERTIFICATION")]
    [InlineData(0.0, false, .001)]
    [InlineData(0.1, false, .001)]
    [InlineData(0.0, true, .001)]
    [InlineData(0.1, true, .001)]
    [InlineData(0.0, false, 1.0)]
    [InlineData(0.1, false, 1.0)]
    [InlineData(0.0, true, 1.0)]
    [InlineData(0.1, true, 1.0)]
    public void Path6_UntrustedGlobalOrTrackedLocator_HasNoDefinitiveCharacter(double delta, bool residualOnly, double searchTolerance)
    {
        var result = XxzCoherenceBlock.TrackDiabolicUnderDelta(7, new Complex(.6788, 0),
            new Complex(-4.557, 0), delta, residualOnly: residualOnly, coalesceTol: searchTolerance);
        Console.WriteLine($"N7 unsafe locator delta={delta:R}, residual={residualOnly}, verdict={result.Verdict}, alg={result.Algebraic}, gap={result.Gap:R}");
        AssertUncertified(result);
    }

    [Fact]
    [Trait("Category", "TASK8_DELTA_CERTIFICATION")]
    public void Path5_TrackedProposal_SearchToleranceDoesNotEnterTheCertificate()
    {
        // The legacy coalesceTol shapes nothing: the certified EP2 is the same at 1e-3 and at 1.0.
        var seed = new Complex(.7581, .260); var lam = new Complex(-5.392, 1.653);
        var tight = XxzCoherenceBlock.TrackDiabolicUnderDelta(6, seed, lam, .02, residualOnly: true, coalesceTol: .001);
        var loose = XxzCoherenceBlock.TrackDiabolicUnderDelta(6, seed, lam, .02, residualOnly: true, coalesceTol: 1.0);
        AssertEp2AtOneOf(tight, N6RungFarA002, "N6 tracked tol=1e-3");
        Assert.Equal(tight, loose);
    }

    [Fact]
    [Trait("Category", "TASK8_DELTA_CERTIFICATION")]
    public void CharacterAtDiabolicNear_RejectsUnrelatedLambdaSeed()
    {
        var result = XxzCoherenceBlock.CharacterAtDiabolicNear(4, 0,
            new Complex(GaloisMonodromyWitness.QEp, 0), new Complex(100, 0));
        AssertUncertified(result); // Unavailable character, not a borrowed unrelated pair.
    }

    [Fact]
    [Trait("Category", "TASK8_DELTA_CERTIFICATION")]
    public void DefaultTracker_UnrelatedLambdaSeed_CannotBorrowGlobalMinimumCharacter()
    {
        var result = XxzCoherenceBlock.TrackDiabolicUnderDelta(5, new Complex(.6407, .180),
            new Complex(100, 100), 0);
        AssertUncertified(result);
    }

    private static Complex[] SpectrumOf(Matrix<Complex> m)
        => m.Evd().EigenValues.ToArray();

    // multiset agreement (every a has a b within tol, and vice versa), robust to ordering, unlike sort-by-Re
    // which is unstable across the Re=-4 cluster (many eigenvalues share Re to float noise).
    private static void AssertSameSpectrum(Complex[] a, Complex[] b, double tol)
    {
        Assert.Equal(a.Length, b.Length);
        foreach (var x in a) Assert.True(b.Min(y => (x - y).Magnitude) < tol, $"{x} has no match in b");
        foreach (var y in b) Assert.True(a.Min(x => (x - y).Magnitude) < tol, $"{y} has no match in a");
    }

    // ΔTask 1: the (q,Δ) builder at Δ=0 IS the XY (SE,DE) block; at N=4 its R=+1 symmetric-sector spectrum
    // must equal F89Path3OcticBlock.BuildSeDeSymBlock(q, γ=1) to machine precision (the trusted anchor).
    [Fact]
    public void XxzBlock_Delta0_MatchesF89SymBlockSpectrum_AtQ2()
    {
        var xxz = XxzCoherenceBlock.SeDeSymSpectrum(4, new Complex(2.0, 0), 0.0);
        var f89 = SpectrumOf(F89Path3OcticBlock.BuildSeDeSymBlock(2.0, 1.0));   // (j=2, γ=1)
        AssertSameSpectrum(xxz, f89, 1e-9);
    }

    // The N=4 reading of DIABOLIC_BY_INTEGRABILITY, certified. Δ=0: the crossing at q_EP is a double zero read
    // geo = alg = 2. Δ = 0.02 and 0.10: it has split into two simple zeros, both at real q on the palindrome line
    // Re λ = −4, each read alg 2, geo 1. The references are the EP locations bracketed by a sign change of the
    // real discriminant (to 1e-15 in 40-digit arithmetic, test_f89_zz_break_surface.py) and their departures. The discriminant is geo vs alg (R-3), NOT EpCharacter.Kind, which labels
    // both Δ = 0.02 EP2s (departures 0.0176 and 0.0222 at |λ| ≈ 4.2) as "Normal".
    [Theory]
    [InlineData(0.02, 0.656935638035895896, 1.303842206475765, 0.0175505831, 0.660248940906498777, 1.30790927659822, 0.0222486857)]
    [InlineData(0.10, 0.650363027620131020, 1.243522164660881, 0.0847991865, 0.667060160649329616, 1.265942303756041, 0.1116521248)]
    public void N4_DeltaResponse_TheDiabolicSplitsIntoTwoRealAxisEp2sOnTheLine(double delta,
        double qLow, double imLow, double depLow, double qHigh, double imHigh, double depHigh)
    {
        double qEp = GaloisMonodromyWitness.QEp;                 // sqrt((-1+sqrt13)/6) ~ 0.658983
        var lamEp = new Complex(-4, 2 * qEp);
        AssertCrossingAt(XxzCoherenceBlock.CharacterAtDiabolicNear(4, 0.0, new Complex(qEp, 0), lamEp),
            new Complex(qEp, 0), "N4 Delta=0");

        // Location model at these EP2s: residual gap g ~ 1e-8, slope 0.67 or 1.47 => |dq| <~ 1e-15; the pair mean's
        // backward error is u·‖M‖ ~ 1e-15. The 1e-12 gates leave a margin of 1000.
        var refs = new[] { (Q: qLow, Im: imLow, Dep: depLow), (Q: qHigh, Im: imHigh, Dep: depHigh) };
        void AssertOnTheLine(XxzCoherenceBlock.DeltaTrackResult r, string label)
        {
            Console.WriteLine($"{label}: {r}");
            Assert.Equal(XxzCoherenceBlock.DeltaFlipVerdict.Defective, r.Verdict);
            Assert.Equal(2, r.Algebraic);
            Assert.Equal(1, r.Geometric);                        // geo<alg => DEFECTIVE (Jordan)
            Assert.Equal(1, r.DiscriminantZeroOrder);            // a simple zero: the square-root branch
            Assert.True(Math.Abs(r.QCandidate.Imaginary) <= 1e-12, $"{label}: Im q {r.QCandidate.Imaginary:E2}");
            Assert.True(Math.Abs(r.LambdaCandidate.Real + 4) <= 1e-12, $"{label}: Re lambda + 4 = {r.LambdaCandidate.Real + 4:E2}");
            var hit = refs.OrderBy(e => Math.Abs(e.Q - r.QCandidate.Real)).First();
            Assert.True(Math.Abs(hit.Q - r.QCandidate.Real) <= 1e-12, $"{label}: q {r.QCandidate.Real:R} vs {hit.Q:R}");
            Assert.True(Math.Abs(hit.Im - r.LambdaCandidate.Imaginary) <= 1e-12, $"{label}: Im lambda {r.LambdaCandidate.Imaginary:R}");
            Assert.True(Math.Abs(hit.Dep - r.Departure) <= DepartureGate, $"{label}: departure {r.Departure:R} vs {hit.Dep:R}");
        }

        AssertOnTheLine(XxzCoherenceBlock.CharacterAtDiabolicNear(4, delta, new Complex(qEp, 0), lamEp), $"N4 tracker Delta={delta}");
        var (first, second) = XxzCoherenceBlock.CertifySplitUnderDelta(4, new Complex(qEp, 0), lamEp, delta);
        AssertOnTheLine(first, $"N4 split Delta={delta} (1)");
        AssertOnTheLine(second, $"N4 split Delta={delta} (2)");
        Assert.True(Math.Abs(first.QCandidate.Real - second.QCandidate.Real) > 1e-3, "the two EP2s are distinct zeros");
    }

    // The two branch laws, read as ratios that must stay constant across decades of the offset δ from the zero:
    // at each Δ = 0.02 EP2 the gap closes as the square root (gap/√δ constant), at the Δ=0 crossing linearly
    // (gap/δ constant). Error model of a ratio: the next Puiseux/Taylor term, κ·δ relative with |κ| ≈ 150 at these
    // EP2s (1.5e-5 at δ = 1e-7), plus the eigenvalue rounding, about 2·dep·u·‖M‖/gap² relative at an EP2 (9e-8 at
    // δ = 1e-9) and u·‖M‖/gap² at the crossing. The
    // 1e-4 gate bounds both over two decades, while the wrong branch would move the ratio tenfold per two decades.
    [Fact]
    public void N4_DeltaResponse_GapLaws_SquareRootAtTheEp2s_LinearAtTheCrossing()
    {
        double qEp = GaloisMonodromyWitness.QEp;
        var lamEp = new Complex(-4, 2 * qEp);
        double[] offsets = { 1e-7, 1e-8, 1e-9 };
        double GapAt(Complex q, double delta, Complex lam)
        {
            var pair = XxzCoherenceBlock.SeDeSymSpectrum(4, q, delta).OrderBy(z => (z - lam).Magnitude).Take(2).ToArray();
            return (pair[0] - pair[1]).Magnitude;
        }
        var (first, second) = XxzCoherenceBlock.CertifySplitUnderDelta(4, new Complex(qEp, 0), lamEp, 0.02);
        foreach (var ep in new[] { first, second })
        {
            var ratios = offsets.Select(d => GapAt(ep.QCandidate + d, 0.02, ep.LambdaCandidate) / Math.Sqrt(d)).ToArray();
            Console.WriteLine($"EP2 q={ep.QCandidate.Real:R}: gap/sqrt(d) = {string.Join(", ", ratios.Select(r => r.ToString("R")))}");
            Assert.True(ratios.Max() / ratios.Min() - 1 < 1e-4, "square-root law");
            // 0.670058 at the lower EP2, 0.667211 at the upper one (independent numpy reading).
            double expected = ep.QCandidate.Real < qEp ? 0.670058 : 0.667211;
            Assert.True(Math.Abs(ratios.Last() - expected) < 1e-5, $"slope {ratios.Last():R} vs {expected}");
        }
        var linear = offsets.Select(d => GapAt(new Complex(qEp + d, 0), 0.0, lamEp) / d).ToArray();
        Console.WriteLine($"crossing: gap/d = {string.Join(", ", linear.Select(r => r.ToString("R")))}");
        Assert.True(linear.Max() / linear.Min() - 1 < 1e-4, "linear law");
        Assert.True(Math.Abs(linear.Last() - 11.5819) < 1e-3, $"crossing slope {linear.Last():R}");
    }

    // Residual-only Δ-test (the N>=6 fix): the full sym spectrum floods on AT-locked degeneracies, so the box
    // scan in TrackDiabolicUnderDelta captures AT crossings at N=6 (q_candidate jumps, false LIFTs). ResidualRootsTrackedXxz
    // uses base-label nearest-neighbour continuation from (q0=2, Δ=0), where the labels equal the F89 residual.
    // This is proposal tracking, not an invariant-set proof; full-block character remains required.
    // The base comparison pins the XxzCoherenceBlock-vs-F89 spectrum convention at path-5.
    [Fact]
    public void ResidualRootsTrackedXxz_Path5_MatchesLocatorResidual_AtBase()
    {
        var xxz = XxzCoherenceBlock.ResidualRootsTrackedXxz(5, new Complex(2, 0), 0.0)
            .OrderBy(z => z.Real).ThenBy(z => z.Imaginary).ToArray();
        var loc = PathKMonodromyScout.ResidualRootsAt(5, new Complex(2, 0))
            .OrderBy(z => z.Real).ThenBy(z => z.Imaginary).ToArray();
        Assert.Equal(32, xxz.Length);                            // F_d degree for path-5 (not 45 = full sym block)
        Assert.Equal(loc.Length, xxz.Length);
        for (int i = 0; i < xxz.Length; i++)
            Assert.True((xxz[i] - loc[i]).Magnitude < 1e-7, $"strand {i}: xxz {xxz[i]} vs locator {loc[i]}");
    }

    private static readonly (Complex Q, Complex Lam, Complex Crossing, (Complex Q, double Dep)[] At002, (Complex Q, double Dep)[] At01)[] N6Seeds =
    {
        (new(0.7090, -0.219), new(-4.151, 1.615), new(0.709029445702654, -0.218648287466866), N6RungNear002, N6RungNear01),   // rung-near
        (new(0.7581, 0.260), new(-5.392, 1.653), new(0.758135158290132, 0.260408965415563), N6RungFarA002, N6RungFarA01),   // rung-far
        (new(1.0561, 0.238), new(-5.187, 1.441), new(1.056063538421554, 0.237812167425590), N6RungFarB002, N6RungFarB01),   // rung-far
    };

    // N=6 residual-tracked proposals: each Δ=0 seed is a certified crossing at its reference q (not an AT capture,
    // which would sit elsewhere). At Δ = 0.1 the first two proposals refine to certified EP2s of their seeds. The
    // third proposal sits away from its seed's zeros; in a numpy port of this path it refines to a coalescence (pair
    // gap 5e-8) whose λ lies 0.30 from the seed, beyond the isolation radius 0.17, and is declined, so here it is held
    // only to the rule that no crossing survives.
    [Fact]
    public void Path5_ResidualProposals_CertifyTheSplitEp2s()
    {
        foreach (var (q, lam, crossing, _, _) in N6Seeds)
            AssertCrossingAt(XxzCoherenceBlock.TrackDiabolicUnderDelta(6, q, lam, 0.0, residualOnly: true), crossing, $"N6 tracked {q} Delta=0");
        AssertEp2AtOneOf(XxzCoherenceBlock.TrackDiabolicUnderDelta(6, N6Seeds[0].Q, N6Seeds[0].Lam, 0.1, residualOnly: true),
            N6Seeds[0].At01, "N6 tracked rung-near Delta=.1");
        AssertEp2AtOneOf(XxzCoherenceBlock.TrackDiabolicUnderDelta(6, N6Seeds[1].Q, N6Seeds[1].Lam, 0.1, residualOnly: true),
            N6Seeds[1].At01, "N6 tracked rung-far Delta=.1");
        AssertNoCrossingSurvives(XxzCoherenceBlock.TrackDiabolicUnderDelta(6, N6Seeds[2].Q, N6Seeds[2].Lam, 0.1, residualOnly: true),
            N6Seeds[2].At01, "N6 tracked rung-far (2) Delta=.1");
    }

    // Seed-following at N=6: each crossing splits into two EP2s at Δ = 0.02, both inside the seed's disk.
    [Fact]
    public void N6_Crossings_SplitIntoTwoEp2s_AtDelta002()
    {
        foreach (var (q, lam, _, at002, _) in N6Seeds)
        {
            var (first, second) = XxzCoherenceBlock.CertifySplitUnderDelta(6, q, lam, 0.02);
            AssertEp2AtOneOf(first, at002, $"N6 split {q} (1)");
            AssertEp2AtOneOf(second, at002, $"N6 split {q} (2)");
            Assert.True((first.QCandidate - second.QCandidate).Magnitude > 1e-3, "two distinct zeros");
        }
    }

    // N=5 default tracker: the three Δ=0 crossings are certified; at Δ = 0.05 (and the clean seed at 0.02) each
    // refines to a certified EP2 of its seed. The defective control is a simple zero at Δ=0 (alg 2, geo 1,
    // departure 0.884) and stays one at Δ = 0.02 (0.890) and 0.05 (0.897): Δ does not make a diabolic of it.
    [Fact]
    public void Path4_DefaultTracker_CertifiesTheSplitEp2s_ControlStaysDefective()
    {
        var seeds = new[]
        {
            (q: new Complex(0.6407, 0.180), lam: new Complex(-4.077, -1.115), crossing: new Complex(0.640730584145048, 0.179536271639987), at005: N5Clean005),
            (q: new Complex(0.7654, 0.024), lam: new Complex(-4.371, -2.056), crossing: new Complex(0.765442564061999, 0.024412312880001), at005: N5NearReal005),
            (q: new Complex(1.9447, 1.217), lam: new Complex(-2.455, -3.473), crossing: new Complex(1.944653817696007, 1.217453894436002), at005: N5Far005),
        };
        foreach (var (q, lam, crossing, at005) in seeds)
        {
            AssertCrossingAt(XxzCoherenceBlock.TrackDiabolicUnderDelta(5, q, lam, 0.0), crossing, $"N5 {q} Delta=0");
            AssertEp2AtOneOf(XxzCoherenceBlock.TrackDiabolicUnderDelta(5, q, lam, 0.05), at005, $"N5 {q} Delta=.05");
        }
        AssertEp2AtOneOf(XxzCoherenceBlock.TrackDiabolicUnderDelta(5, new Complex(0.6407, 0.180), new Complex(-4.077, -1.115), 0.02),
            N5Clean002, "N5 clean Delta=.02");

        var ctrlQ = new Complex(0.9938, 0.183);
        var ctrlLam = new Complex(-4.712, 0.824);
        AssertEp2AtOneOf(XxzCoherenceBlock.TrackDiabolicUnderDelta(5, ctrlQ, ctrlLam, 0),
            new[] { (new Complex(0.993804650949799, 0.182516992002082), 0.8837939194) }, "N5 control Delta=0");
        AssertEp2AtOneOf(XxzCoherenceBlock.TrackDiabolicUnderDelta(5, ctrlQ, ctrlLam, .02),
            new[] { (new Complex(1.006659475472186, 0.172059579386927), 0.8901591604) }, "N5 control Delta=.02");
        AssertEp2AtOneOf(XxzCoherenceBlock.TrackDiabolicUnderDelta(5, ctrlQ, ctrlLam, .05),
            new[] { (new Complex(1.026084821283453, 0.157520436063259), 0.8972376254) }, "N5 control Delta=.05");
    }

    // Seed-following at N=5: the clean and near-real crossings split into two certified EP2s at Δ = 0.02 and 0.05.
    // The far crossing's first zero is a certified EP2; its second is reached, but its λ has moved 0.41 (Δ = 0.02)
    // and 0.69 (Δ = 0.05) from the seed, beyond the isolation radius 0.16, so the certificate declines it.
    [Theory]
    [InlineData(0.02)]
    [InlineData(0.05)]
    public void N5_Crossings_SplitIntoEp2s_SeedFollowing(double delta)
    {
        void AssertBothEp2((XxzCoherenceBlock.DeltaTrackResult First, XxzCoherenceBlock.DeltaTrackResult Second) split,
            (Complex Q, double Dep)[] refs, string label)
        {
            AssertEp2AtOneOf(split.First, refs, $"{label} (1)");
            AssertEp2AtOneOf(split.Second, refs, $"{label} (2)");
            Assert.True((split.First.QCandidate - split.Second.QCandidate).Magnitude > 1e-3, "two distinct zeros");
        }
        AssertBothEp2(XxzCoherenceBlock.CertifySplitUnderDelta(5, new Complex(0.6407, 0.180), new Complex(-4.077, -1.115), delta),
            delta == 0.02 ? N5Clean002 : N5Clean005, $"N5 split clean Delta={delta}");
        AssertBothEp2(XxzCoherenceBlock.CertifySplitUnderDelta(5, new Complex(0.7654, 0.024), new Complex(-4.371, -2.056), delta),
            delta == 0.02 ? N5NearReal002 : N5NearReal005, $"N5 split near-real Delta={delta}");

        var farRefs = delta == 0.02 ? N5Far002 : N5Far005;
        var (farFirst, farSecond) = XxzCoherenceBlock.CertifySplitUnderDelta(5, new Complex(1.9447, 1.217), new Complex(-2.455, -3.473), delta);
        AssertEp2AtOneOf(farFirst, farRefs, $"N5 split far Delta={delta} (1)");
        AssertUncertified(farSecond);
        Assert.Contains(farRefs, e => (e.Q - farSecond.QCandidate).Magnitude <= LocationGate);
    }

    // ---- Fixed-complement compressed XXZ proposals ----
    // The tracked residual path (ResidualRootsTrackedXxz) breaks at k=6 (the F_53 AT-degeneracy flood: leaked
    // AT strands wearing residual labels, min-gap zero across the box), exactly as the XY tracker did before
    // the Delta=0 exact restriction. The XXZ port supplies compressed proposals, not residual eigenvalues:
    // M_xxz(q,Δ) = (A + qC + qΔ·G)/2 in F89's mirror basis (G the ZZ-
    // frequency generator, diagonal −2i·zzDiag), compressed onto the q-independent AT complement U_res. These
    // gates validate the matrix port; the certificate is read on the full block after the Newton refinement.

    // Foundation: the XXZ block built in F89's mirror basis (A + qC + qΔ·G) has the SAME spectrum as the
    // independently-constructed XxzCoherenceBlock.BuildSym at Δ≠0. Pins the ZZ generator (the −2i·zzDiag
    // diagonal and its ×2-cleared scaling: a reflection-invariant diagonal contributes 2·value to 2M for both
    // orbit lengths) against the trusted block, with NO compression involved.
    [Theory]
    [InlineData(5)]   // N=6
    [InlineData(6)]   // N=7
    public void AllRootsXxz_F89Basis_MatchesXxzBlockSpectrum_AtDelta(int k)
    {
        var q = new Complex(2.0, 0.3);
        const double delta = 0.1;
        var f89 = PathKMonodromyScout.AllRootsXxz(k, q, delta);
        var xxz = XxzCoherenceBlock.SeDeSymSpectrum(k + 1, q, delta);
        AssertSameSpectrum(f89, xxz, 1e-7);
    }

    // At Δ=0 the exact XXZ residual roots ARE the XY exact residual roots (G drops out): a wiring guard that
    // the port does not perturb the established Δ=0 science.
    [Theory]
    [InlineData(5)]
    [InlineData(6)]
    public void ResidualRootsCompressedXxz_Delta0_EqualsXyExactResidual(int k)
    {
        var q = new Complex(1.7, 0.0);
        var xxz = PathKMonodromyScout.ResidualRootsCompressedXxz(k, q, 0.0);
        var xy = PathKMonodromyScout.ResidualRootsExact(k, q);
        AssertSameSpectrum(xxz, xy, 1e-9);
    }

    // Compression gives 32/53 distinct proposal roots in this generic sample. It is not an invariant
    // residual subspace at Delta!=0 and does not certify a subset of the full spectrum or a perturbation bound.
    [Theory]
    [InlineData(5, 32)]
    [InlineData(6, 53)]
    public void ResidualRootsCompressedXxz_HasDistinctProposalRoots_AtDelta(int k, int fdDegree)
    {
        var q = new Complex(2.0, 0.3);
        const double delta = 0.1;
        var res = PathKMonodromyScout.ResidualRootsCompressedXxz(k, q, delta);
        Assert.Equal(fdDegree, res.Length);
        Assert.True(PathKMonodromyScout.MinGap(res) > 1e-6,
            $"the sampled compressed roots must be distinct; got min gap {PathKMonodromyScout.MinGap(res):E2}");
    }

    // N=6 compressed proposals: the same crossings at Δ=0 (where the compression is exact), and at Δ = 0.1 no
    // crossing survives; whatever the refined proposal reaches is certified or declined, never read off a gap.
    [Fact]
    public void Path5_CompressedProposals_CertifyTheCrossings_NoneSurvivesDelta()
    {
        foreach (var (q, lam, crossing, _, at01) in N6Seeds)
        {
            AssertCrossingAt(XxzCoherenceBlock.TrackDiabolicUnderDelta(6, q, lam, 0.0, exact: true), crossing, $"N6 compressed {q} Delta=0");
            AssertNoCrossingSurvives(XxzCoherenceBlock.TrackDiabolicUnderDelta(6, q, lam, 0.1, exact: true), at01, $"N6 compressed {q} Delta=.1");
        }
    }

    [Theory]
    [Trait("Category", "TASK8_DELTA_CERTIFICATION")]
    [InlineData(0.0, 0.0, true)]
    [InlineData(0.0001, 0.0, false)] // Close but split full eigenvalues.
    [InlineData(0.0, 0.01, false)] // Coincident full pair unrelated to proposed midpoint.
    public void CompressedCertification_RequiresCoincidenceAndCorrespondence(double split, double midpoint, bool certified)
    {
        var matrix = MathNet.Numerics.LinearAlgebra.Matrix<Complex>.Build.DenseDiagonal(3, 3,
            i => i == 2 ? 1 : i * split);
        var result = XxzCoherenceBlock.CertifyFullBlockProposal(matrix, Complex.Zero, new Complex(midpoint, 0));
        Assert.Equal(certified, result.IsCertifiedDiabolic);
        Assert.Equal(certified ? XxzCoherenceBlock.DeltaFlipVerdict.Diabolic : XxzCoherenceBlock.DeltaFlipVerdict.Uncertified, result.Verdict);
        Assert.True(Math.Abs(result.Gap - split) < 1e-12);
    }

    // The character is geometric against algebraic multiplicity (R-3), with the departure only guarding the
    // Jordan side against a split crossing read at the rounding level.
    [Theory]
    [Trait("Category", "TASK8_DELTA_CERTIFICATION")]
    [InlineData(2, 2, 1.5e-7, 7e-7, XxzCoherenceBlock.DeltaFlipVerdict.Diabolic)]      // departure at its rounding floor
    [InlineData(2, 1, 0.0015, 7e-7, XxzCoherenceBlock.DeltaFlipVerdict.Defective)]     // the N=7 Δ=0.02 coupling
    [InlineData(2, 1, 1.5e-7, 7e-7, XxzCoherenceBlock.DeltaFlipVerdict.Uncertified)]   // a split crossing, not a Jordan block
    [InlineData(3, 1, 0.1, 7e-7, XxzCoherenceBlock.DeltaFlipVerdict.Uncertified)]      // an EP3 is not this certificate's object
    [InlineData(1, 1, 0.0, 7e-7, XxzCoherenceBlock.DeltaFlipVerdict.Uncertified)]      // nothing coalesces
    public void CharacterVerdict_ReadsGeometricAgainstAlgebraicMultiplicity(int alg, int geo, double departure,
        double floor, XxzCoherenceBlock.DeltaFlipVerdict expected)
        => Assert.Equal(expected, XxzCoherenceBlock.CertifiedCharacterVerdict(alg, geo, departure, floor));

    // A Jordan block whose coupling is small beside |λ|: the case the relative-departure Kind misses. The
    // certificate reads it Defective with its exact coupling as departure; the same block uncoupled reads Diabolic.
    [Theory]
    [Trait("Category", "TASK8_DELTA_CERTIFICATION")]
    [InlineData(0.0015, XxzCoherenceBlock.DeltaFlipVerdict.Defective)]
    [InlineData(0.0, XxzCoherenceBlock.DeltaFlipVerdict.Diabolic)]
    public void CertifyFullBlockProposal_SmallCouplingJordanBlockAtLargeLambda(double coupling,
        XxzCoherenceBlock.DeltaFlipVerdict expected)
    {
        var lam = new Complex(-4.94, 0);
        var m = Matrix<Complex>.Build.DenseOfArray(new Complex[,]
        {
            { lam, coupling, 0 },
            { 0, lam, 0 },
            { 0, 0, lam + 1 },
        });
        var r = XxzCoherenceBlock.CertifyFullBlockProposal(m, Complex.Zero, lam);
        Console.WriteLine(r);
        Assert.Equal(expected, r.Verdict);
        Assert.Equal(2, r.Algebraic);
        Assert.Equal(expected == XxzCoherenceBlock.DeltaFlipVerdict.Defective ? 1 : 2, r.Geometric);
        // departure = the Jordan coupling exactly; its computation cancels ‖A‖² − Σ|λ|² (≈ 49), rounding ≈ u·49/c.
        if (coupling > 0) Assert.True(Math.Abs(r.Departure - coupling) < 1e-9, $"departure {r.Departure:R}");
    }

    [Fact]
    [Trait("Category", "TASK8_DELTA_CERTIFICATION")]
    public void UncertifiedSurvival_IsUnknown_NotFalse()
    {
        var unknown = new XxzCoherenceBlock.DeltaTrackResult(
            XxzCoherenceBlock.DeltaFlipVerdict.Uncertified, 0, 0, double.NaN, Complex.Zero, Complex.Zero, 0.02);
        Assert.Null(unknown.Survived);
        Assert.False(unknown.IsCertifiedDiabolic);
        Assert.True((unknown with { Verdict = XxzCoherenceBlock.DeltaFlipVerdict.Diabolic }).Survived);
        Assert.False((unknown with { Verdict = XxzCoherenceBlock.DeltaFlipVerdict.Defective }).Survived);
        Assert.False((unknown with { Verdict = XxzCoherenceBlock.DeltaFlipVerdict.Lifted }).Survived);
    }

    [Theory]
    [Trait("Category", "TASK8_DELTA_CERTIFICATION")]
    [InlineData(0.0, 0.001)]
    [InlineData(0.0, 1000.0)]
    [InlineData(0.1, 1000.0)]
    public void CompressedProposal_RejectsWrongLambdaSeed(double delta, double proposalTolerance)
    {
        var q = new Complex(0.6788, 0);
        var result = XxzCoherenceBlock.TrackDiabolicUnderDelta(7, q, new Complex(100, 0), delta,
            coalesceTol: proposalTolerance, exact: true);
        AssertUncertified(result);

        // Same public path and q with the actual seed must retain its calibrated Delta=0 pair.
        var valid = XxzCoherenceBlock.TrackDiabolicUnderDelta(7, q, new Complex(-4.557, 0), 0,
            coalesceTol: proposalTolerance, exact: true);
        AssertCrossingAt(valid, new Complex(0.678816427645286, 0), "N7 0.6788 Delta=0");
    }

    [Theory]
    [Trait("Category", "TASK8_DELTA_CERTIFICATION")]
    [InlineData(0.02)]
    [InlineData(0.0001)]
    public void CompressedProposal_SearchToleranceDoesNotEnterTheCertificate(double delta)
    {
        var q = new Complex(0.6788, 0); var lam = new Complex(-4.557, 0);
        var tight = XxzCoherenceBlock.TrackDiabolicUnderDelta(7, q, lam, delta, coalesceTol: .001, exact: true);
        var loose = XxzCoherenceBlock.TrackDiabolicUnderDelta(7, q, lam, delta, coalesceTol: 1.0, exact: true);
        AssertNoCrossingSurvives(tight, delta == 0.02 ? N7Q06788At002 : N7Q06788At00001, $"N7 0.6788 Delta={delta}");
        Assert.Equal(tight, loose);
    }

    // N=7 compressed proposals from the four real-q crossings. Δ=0: each a certified crossing (the compression is
    // exact there). Δ > 0: no crossing survives near the seed; the refined proposal is a certified EP2 or, when the
    // zero it reaches has left the seed's Δ=0 isolation disk, Uncertified. The EP2s themselves are certified
    // seed-following in N7_RealQCrossings_SplitIntoEp2s_SeedFollowing.
    [Theory]
    [Trait("Category", "TASK8_DELTA_CERTIFICATION")]
    [InlineData(1.1264, -4.942, 1.126448513252489811, 0.02)]
    [InlineData(1.3038, -5.171, 1.303811904192148, 0.02)]
    [InlineData(2.6280, -4.343, 2.628016284405854, 0.02)]
    [InlineData(0.6788, -4.557, 0.678816427645286, 0.02)]
    [InlineData(0.6788, -4.557, 0.678816427645286, 0.10)]
    [InlineData(1.1264, -4.942, 1.126448513252489811, 0.10)]
    [InlineData(1.3038, -5.171, 1.303811904192148, 0.10)]
    [InlineData(2.6280, -4.343, 2.628016284405854, 0.10)]
    public void Path6_RealQCrossings_CompressedProposals_CertifyTheCrossings_NoneSurvivesDelta(double realQ, double realLambda,
        double crossing, double delta)
    {
        var q = new Complex(realQ, 0); var lam = new Complex(realLambda, 0);
        AssertCrossingAt(XxzCoherenceBlock.TrackDiabolicUnderDelta(7, q, lam, 0.0, exact: true), new Complex(crossing, 0), $"N7 {realQ} Delta=0");
        var d = XxzCoherenceBlock.TrackDiabolicUnderDelta(7, q, lam, delta, exact: true);
        var refs = (realQ, delta) switch
        {
            (1.1264, 0.02) => N7Q11264At002,
            (1.3038, 0.02) => N7Q13038At002,
            (2.6280, 0.02) => N7Q26280At002,
            (0.6788, 0.02) => N7Q06788At002,
            (0.6788, 0.10) => N7Q06788At01,
            (1.1264, 0.10) => N7Q11264At01,
            _ => Array.Empty<(Complex Q, double Dep)>(),
        };
        AssertNoCrossingSurvives(d, refs, $"N7 compressed {realQ} Delta={delta}");
    }

    // Seed-following at N=7 for all four real-q crossings. Δ = 0.02: each splits into two certified EP2s, except that
    // the q = 2.6280 crossing's second zero (3.0631 − 0.8697i by continuation) lies too far along q to be reached
    // from the reflected start, so its second result is that EP2 or Uncertified. Δ = 0.10: the q = 0.6788 pair is
    // certified again; the q = 1.1264 zeros and the first q = 1.3038 zero are reached but lie beyond the seed's
    // isolation radius, so they are declined (Uncertified, not a death verdict). Where the locator reaches no zero
    // (the second q = 1.3038 zero, the q = 2.6280 zeros at Δ = 0.10) the test asks only that no crossing survives,
    // so a better locator that certifies them does not fail it. Δ = 1e-4: the q = 0.6788 crossing has already
    // split, 5.6e-5 apart in q, departure 1.98e-4 each.
    [Theory]
    [InlineData(1.1264, -4.942, 0.02)]
    [InlineData(1.3038, -5.171, 0.02)]
    [InlineData(2.6280, -4.343, 0.02)]
    [InlineData(0.6788, -4.557, 0.02)]
    [InlineData(0.6788, -4.557, 0.0001)]
    [InlineData(0.6788, -4.557, 0.10)]
    [InlineData(1.1264, -4.942, 0.10)]
    [InlineData(1.3038, -5.171, 0.10)]
    [InlineData(2.6280, -4.343, 0.10)]
    public void N7_RealQCrossings_SplitIntoEp2s_SeedFollowing(double realQ, double realLambda, double delta)
    {
        var (first, second) = XxzCoherenceBlock.CertifySplitUnderDelta(7, new Complex(realQ, 0), new Complex(realLambda, 0), delta);
        switch (realQ, delta)
        {
            case (1.1264, 0.10):
                AssertUncertified(first);
                AssertUncertified(second);
                foreach (var r in new[] { first, second })     // the zeros are reached: each is a reference zero
                    Assert.Contains(N7Q11264At01, e => (e.Q - r.QCandidate).Magnitude <= LocationGate);
                return;
            case (1.3038, 0.10):
                AssertUncertified(first);                          // reached, and declined by the disk rule
                Assert.Contains(N7Q13038At01, e => (e.Q - first.QCandidate).Magnitude <= LocationGate);
                AssertNoCrossingSurvives(second, N7Q13038At01, "N7 split 1.3038 Delta=0.1 (2)");
                return;
            case (2.6280, 0.10):                                   // a locator null here is not a verdict
                AssertNoCrossingSurvives(first, Array.Empty<(Complex Q, double Dep)>(), "N7 split 2.628 Delta=0.1 (1)");
                AssertNoCrossingSurvives(second, Array.Empty<(Complex Q, double Dep)>(), "N7 split 2.628 Delta=0.1 (2)");
                return;
            case (2.6280, 0.02):
                AssertEp2AtOneOf(first, N7Q26280At002, "N7 split 2.628 (1)");
                if (second.Verdict != XxzCoherenceBlock.DeltaFlipVerdict.Uncertified)
                    AssertEp2AtOneOf(second, N7Q26280At002, "N7 split 2.628 (2)");
                return;
        }
        var refs = (realQ, delta) switch
        {
            (1.1264, 0.02) => N7Q11264At002,
            (1.3038, 0.02) => N7Q13038At002,
            (0.6788, 0.02) => N7Q06788At002,
            (0.6788, 0.0001) => N7Q06788At00001,
            _ => N7Q06788At01,
        };
        AssertEp2AtOneOf(first, refs, $"N7 split {realQ} Delta={delta} (1)");
        AssertEp2AtOneOf(second, refs, $"N7 split {realQ} Delta={delta} (2)");
        Assert.True((first.QCandidate - second.QCandidate).Magnitude > 1e-5, "two distinct zeros");
    }

    // Seed-following at N=7, no box scan: the q = 1.1264 crossing is certified at its 18-digit location, and at
    // Δ = 0.02 it has split into the two EP2s 1.12701920 − 0.00583855i (departure 0.0015) and
    // 1.12464902 − 0.01707280i (0.0165), both read alg 2, geo 1. The first is the coupling EpCharacter.Kind reads
    // as Normal (relative departure 0.0002).
    [Fact]
    public void N7_RealQCrossing_SplitsIntoTwoEp2s_SeedFollowing()
    {
        var seed = new Complex(1.1264, 0); var lam = new Complex(-4.942, 0);
        var crossing = XxzCoherenceBlock.CertifyCoalescenceNear(7, 0.0, seed, lam);
        AssertCrossingAt(crossing, new Complex(1.126448513252489811, 0), "N7 1.1264 Delta=0");
        Assert.True(Math.Abs(crossing.LambdaCandidate.Real - -4.941857490410003) <= LocationGate
                    && Math.Abs(crossing.LambdaCandidate.Imaginary) <= LocationGate, $"lambda {crossing.LambdaCandidate}");
        var (first, second) = XxzCoherenceBlock.CertifySplitUnderDelta(7, seed, lam, 0.02);
        AssertEp2AtOneOf(first, N7Q11264At002, "N7 split (1)");
        AssertEp2AtOneOf(second, N7Q11264At002, "N7 split (2)");
        Assert.True((first.QCandidate - second.QCandidate).Magnitude > 1e-3, "two distinct zeros");
    }
}
