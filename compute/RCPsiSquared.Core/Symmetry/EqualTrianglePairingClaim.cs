using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Numerics;
using Bond = RCPsiSquared.Core.Symmetry.FluxAreaFifthMomentClaim.Bond;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Flux reversal on the equal triangle under cooling and heating
/// (<c>docs/proofs/PROOF_SEVENTH_TRACE_DEGREE_FOUR_TRIANGLE.md</c>, Theorem 3; <c>docs/ANALYTICAL_FORMULAS.md</c>,
/// F137, its open question on flux reversal, answered on this triangle).
///
/// <para>Three sites, every pair joined by a hopping bond with one |t| = r and one ZZ coupling Δ, cooling γ↓_l and
/// heating γ↑_l on every site, no fields, no dephasing, the loop's imaginary part Im = Im(t₀₁·t₁₂·t₂₀) not zero (a flux
/// other than 0 and π). The spectrum of L is paired, its eigenvalues in pairs λ and −c − λ, exactly when the rates have
/// one ratio γ↑/γ↓ on every damped site (cooling or heating alone included) or two sites have the same rate point.
/// Those two are F137's carriers (the one-ratio transpose and its limit at ratio 0 or ∞ through the characteristic
/// polynomial; the reflection fixing one site of the triangle). Every other rate point is ruled out by an odd centred
/// trace: the fifth where the three rate points are not on one line (<see cref="FluxAreaFifthMomentClaim"/>), the
/// seventh on a line off equal heights (<see cref="SeventhMomentDegreeFourClaim"/>, Theorem 2 with b = 0), and on equal
/// heights s_l = S, where the fifth and seventh are zero, the ninth,</para>
///
/// <para>  Tr((L + c/2)⁹) = 36864·Im·(r² + 7Δ²)·S·V(δ),  δ_l = (γ↓_l − γ↑_l)/2,</para>
///
/// <para>V(x) = (x₀ − x₁)(x₁ − x₂)(x₂ − x₀). In the integer form the witness reads, with S′ = γ↓ + γ↑ the common value
/// and a = γ↓ − γ↑: Tr((2L + c)⁹) = 1179648·Im·(r² + 7Δ²)·S′·V(a).</para>
///
/// <para>Live witness: <c>compute/RCPsiSquared.Diagnostics/Foundation/EqualTrianglePairingWitness.cs</c>,
/// <c>inspect --root equaltriangle</c>: on two equal triangles the verdict read from the rates alone against the
/// generator's odd traces up to the ninth and against its characteristic polynomial by Samuelson and Berkowitz,
/// without division and through no trace, whose coefficients of the wrong parity are all zero exactly when the spectrum
/// is paired, read on the carrier rows and on the break rows alike. Parents:
/// <see cref="FluxAreaFifthMomentClaim"/> and <see cref="SeventhMomentDegreeFourClaim"/>.</para></summary>
public sealed class EqualTrianglePairingClaim : Claim
{
    public FluxAreaFifthMomentClaim Fifth { get; }
    public SeventhMomentDegreeFourClaim Seventh { get; }

    /// <summary>What the rates decide on the equal triangle: a carrier pairs the spectrum, or the first odd centred
    /// trace that is not zero.</summary>
    public enum Verdict { PairedOneRatio, PairedSiteReflection, BrokenAtFifth, BrokenAtSeventh, BrokenAtNinth }

    public static bool IsPaired(Verdict v) => v is Verdict.PairedOneRatio or Verdict.PairedSiteReflection;

    /// <summary>Three bonds' worth of couplings on three sites: every pair hopping with one |t|² and one Δ.</summary>
    public static bool IsEqualTriangle(IReadOnlyList<Bond> bonds)
    {
        var t = FluxAreaFifthMomentClaim.HopAmplitudes(3, bonds);
        var zz = ZZ(bonds);
        static BigInteger Sq(GaussianInteger z) => z.Re * z.Re + z.Im * z.Im;
        var r01 = Sq(t[0, 1]);
        return !r01.IsZero && r01 == Sq(t[1, 2]) && r01 == Sq(t[2, 0])
               && zz[0, 1] == zz[1, 2] && zz[1, 2] == zz[0, 2];
    }

    private static BigInteger[,] ZZ(IReadOnlyList<Bond> bonds)
    {
        var zz = new BigInteger[3, 3];
        foreach (var b in bonds)
        {
            zz[b.A, b.B] += b.Delta;
            zz[b.B, b.A] += b.Delta;
        }
        return zz;
    }

    private static BigInteger LoopIm(IReadOnlyList<Bond> bonds)
    {
        var t = FluxAreaFifthMomentClaim.HopAmplitudes(3, bonds);
        return (t[0, 1] * t[1, 2] * t[2, 0]).Im;
    }

    private static void Check(IReadOnlyList<Bond> bonds, IReadOnlyList<long> cool, IReadOnlyList<long> heat)
    {
        SeventhMomentDegreeFourClaim.CheckRates(3, cool, heat);
        if (!IsEqualTriangle(bonds)) throw new ArgumentException("not a triangle with one |t| and one Δ on its three bonds", nameof(bonds));
        if (LoopIm(bonds).IsZero) throw new ArgumentException("the flux is 0 or π, where flux reversal changes nothing", nameof(bonds));
        for (int l = 0; l < 3; l++)
            if (cool[l] < 0 || heat[l] < 0) throw new ArgumentException($"a negative rate on site {l}");
    }

    /// <summary>The verdict from the rates alone.</summary>
    public static Verdict Decide(IReadOnlyList<Bond> bonds, IReadOnlyList<long> cool, IReadOnlyList<long> heat)
    {
        Check(bonds, cool, heat);
        BigInteger Omega(int a, int b) => (BigInteger)cool[a] * heat[b] - (BigInteger)cool[b] * heat[a];
        // one ratio on every damped site: rates not negative, so every pairwise cross product zero (an undamped site fits any)
        if (Omega(0, 1).IsZero && Omega(1, 2).IsZero && Omega(0, 2).IsZero) return Verdict.PairedOneRatio;
        for (int a = 0; a < 3; a++)
            for (int b = a + 1; b < 3; b++)
                if (cool[a] == cool[b] && heat[a] == heat[b]) return Verdict.PairedSiteReflection;
        if (!(Omega(0, 1) + Omega(1, 2) + Omega(2, 0)).IsZero) return Verdict.BrokenAtFifth;
        bool equalHeights = cool[0] + heat[0] == cool[1] + heat[1] && cool[1] + heat[1] == cool[2] + heat[2];
        return equalHeights ? Verdict.BrokenAtNinth : Verdict.BrokenAtSeventh;
    }

    /// <summary>Tr((2L + c)⁹) = 512·Tr((L + c/2)⁹) on equal heights: 1179648·Im·(r² + 7Δ²)·S′·V(a).</summary>
    public static BigInteger NinthOnEqualHeights(IReadOnlyList<Bond> bonds, IReadOnlyList<long> cool, IReadOnlyList<long> heat)
    {
        Check(bonds, cool, heat);
        BigInteger s = cool[0] + heat[0];
        if (cool[1] + heat[1] != s || cool[2] + heat[2] != s) throw new ArgumentException("the heights are not equal");
        var t = FluxAreaFifthMomentClaim.HopAmplitudes(3, bonds);
        var r2 = t[0, 1].Re * t[0, 1].Re + t[0, 1].Im * t[0, 1].Im;
        var delta = ZZ(bonds)[0, 1];
        BigInteger a0 = cool[0] - heat[0], a1 = cool[1] - heat[1], a2 = cool[2] - heat[2];
        var v = (a0 - a1) * (a1 - a2) * (a2 - a0);
        return 1179648 * LoopIm(bonds) * (r2 + 7 * delta * delta) * s * v;
    }

    public EqualTrianglePairingClaim(FluxAreaFifthMomentClaim fifth, SeventhMomentDegreeFourClaim seventh)
        : base("Under cooling and heating alone, on a triangle with one |t| and one Δ at a flux other than 0 and π, the spectrum is paired " +
               "exactly at one ratio γ↑/γ↓ on every damped site or two equal rate points, F137's carriers; elsewhere the fifth, the seventh or, on " +
               "equal heights, the ninth centred trace 36864·Im·(r² + 7Δ²)·S·V(δ) is not zero",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_SEVENTH_TRACE_DEGREE_FOUR_TRIANGLE.md (Theorem 3); docs/ANALYTICAL_FORMULAS.md (F137)")
    {
        Fifth = fifth ?? throw new ArgumentNullException(nameof(fifth));
        Seventh = seventh ?? throw new ArgumentNullException(nameof(seventh));
    }

    public override string DisplayName =>
        "Flux reversal on the equal triangle under cooling and heating: paired only where a carrier applies";

    public override string Summary =>
        "a triangle with one |t| and one Δ, cooling and heating per site, no fields, flux other than 0 and π: the spectrum pairs exactly at one ratio " +
        "on every damped site (cooling or heating alone included) or two equal rate points; off a line the fifth trace breaks it, on a line the " +
        "seventh off equal heights, on equal heights the ninth, 36864·Im·(r² + 7Δ²)·S·V(δ) (Tier 1 derived; live witness: inspect --root equaltriangle).";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the cases",
                summary: "the fifth trace, flux times the area of the rate triangle, leaves collinear rate points; on the triangle with b = 0 the seventh is " +
                         "−3360·e·Im·(1 + m)²·V(γ↓) on γ↑ = e + m·γ↓ and 3360·e·Im·V(γ↑) on γ↓ = e, zero at e = 0 (one ratio), at V = 0 (two equal " +
                         "rate points, the reflection) and at m = −1 (equal heights, the arrows free), where the ninth decides");
            yield return new InspectableNode("the ninth on equal heights",
                summary: "every part odd under flux reversal and an odd permutation reversing the flux, so each part is alternating in the arrows and divisible " +
                         "by V(δ); degrees zero, two, six and eight zero, degree four α·S·V(δ) with α = Im·(A·r² + B·Δ²), A = 36864 and B = 7A from the fit");
            yield return new InspectableNode("forced, not met",
                summary: "every point where the equal triangle pairs is forced by a carrier; on one triangle with equal weights and less symmetry " +
                         "(the proof's §6) the ninth leaves a surface of equal heights, and the eleventh, read in double precision at one point of it, " +
                         "breaks it there; on every other triangle F137's question stays open");
        }
    }
}
