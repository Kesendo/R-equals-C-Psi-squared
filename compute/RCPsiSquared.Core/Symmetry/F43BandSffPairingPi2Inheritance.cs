using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F43 closed form (Tier 1, proven D09): the frequency SFF of a
/// decay-rate band equals that of its palindrome-reflected band.
///
/// <code>
///   I' = { 2Nγ − d : d in I }
///   K_I(t) = K_I'(t)
/// </code>
///
/// <para>The full-spectrum palindrome sends <c>λ = −d + iω</c> to
/// <c>−λ − 2Nγ = −(2Nγ−d) − iω</c>. It therefore gives a
/// multiplicity-preserving bijection between reflected decay-rate bands and
/// complex-conjugates their frequency trace amplitudes. Their normalized and
/// unnormalized frequency SFFs agree at every real time.</para>
///
/// <para>At real Hamiltonian parameters the Absorption Theorem reads
/// <c>d = 2γ⟨n_XY⟩</c>, so the same reflection is
/// <c>⟨n_XY⟩ ↦ N−⟨n_XY⟩</c>. These are average-light bands. The Hamiltonian
/// mixes Pauli-string XY weights by ±2, so a generic Liouvillian eigenmode is
/// not assigned to an invariant fixed-integer-weight sector.</para>
///
/// <para>For a connected uniform chain with strictly positive uniform dephasing, the exact endpoint eigenspaces contain
/// the N+1 stationary modes and their N+1 palindrome partners at
/// <c>λ = −2Nγ</c>. Their frequencies are all zero. Hence their normalized
/// frequency SFF is the constant 1 for every time (unnormalized: (N+1)²), not
/// an impulse at time zero. A finite-width numerical band can contain more
/// eigenvalues and is a separate sampled object. At γ=0 the two endpoints
/// collapse into the Hamiltonian commutator kernel, whose multiplicity is not N+1.</para>
///
/// <para>Pi2-Foundation anchor: <see cref="ReflectionCoefficient"/> = 2 = a_0
/// in both the reflection span 2Nγ and the endpoint rate. The structural band
/// bijection is inherited from <see cref="F1Pi2Inheritance"/>.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F43 +
/// <c>docs/proofs/derivations/D09_SECTOR_SFF_PAIRING.md</c> +
/// <c>experiments/SPECTRAL_FORM_FACTOR.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F1Pi2Inheritance.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F43BandSffPairingPi2Inheritance : Claim, IZ2AxisClaim
{
    public Z2Axis Z2Axis => Z2Axis.BitB;

    /// <summary>No bespoke typed bit_a twin. Global Hadamard conjugation transports
    /// the spectral statement from Z-dephasing to X-dephasing.</summary>
    public Claim? BitATwin => null;

    public BitATwinClassification BitATwinStatus =>
        BitATwinClassification.CoveredByHadamardDuality;

    public Pi2DyadicLadderClaim Ladder { get; }
    public F1Pi2Inheritance F1 { get; }

    /// <summary>The coefficient 2 in the rate reflection
    /// <c>d ↦ 2Nγ−d</c>. Live from Pi2DyadicLadder a_0.</summary>
    public double ReflectionCoefficient => Ladder.Term(0);

    /// <summary>Reflect an average-light value under the palindrome:
    /// <c>⟨n_XY⟩ ↦ N−⟨n_XY⟩</c>. Fractional values are intentional because
    /// Hamiltonian-mixed eigenmodes need not have integer average light.</summary>
    public double ReflectedAverageLight(double averageLight, int N)
    {
        ValidateN(N);
        if (!double.IsFinite(averageLight) || averageLight < 0.0 || averageLight > N)
            throw new ArgumentOutOfRangeException(nameof(averageLight), averageLight,
                $"average light must be finite and in [0, {N}].");
        return N - averageLight;
    }

    /// <summary>The reflection axis in average-light coordinates.</summary>
    public double MirrorAxis(int N)
    {
        ValidateN(N);
        return N / 2.0;
    }

    /// <summary>Whether an average-light value lies on the reflection axis.</summary>
    public bool IsSelfReflected(double averageLight, int N) =>
        ReflectedAverageLight(averageLight, N) == averageLight;

    /// <summary>Reflect a decay rate in the uniform-dephasing interval
    /// <c>[0, 2Nγ]</c>.</summary>
    public double ReflectedDecayRate(double decayRate, int N, double gamma)
    {
        ValidateN(N);
        ValidateGamma(gamma);
        double span = ReflectionCoefficient * N * gamma;
        if (!double.IsFinite(decayRate) || decayRate < 0.0 || decayRate > span)
            throw new ArgumentOutOfRangeException(nameof(decayRate), decayRate,
                $"decay rate must be finite and in [0, {span}].");
        return span - decayRate;
    }

    /// <summary>The nonzero endpoint rate <c>2Nγ</c>, the palindrome partner of
    /// the stationary endpoint.</summary>
    public double EndpointPartnerRate(int N, double gamma)
    {
        ValidateN(N);
        ValidatePositiveGamma(gamma);
        return ReflectionCoefficient * N * gamma;
    }

    /// <summary>N+1 stationary modes, and therefore N+1 partner modes at
    /// <c>−2Nγ</c>, for a connected chain at strictly positive uniform dephasing.</summary>
    public int EndpointMultiplicity(int N, double gamma)
    {
        ValidateN(N);
        ValidatePositiveGamma(gamma);
        return N + 1;
    }

    /// <summary>The normalized frequency SFF of either exact zero-frequency
    /// endpoint eigenspace at strictly positive uniform dephasing. It is 1 at every finite real time.</summary>
    public double EndpointNormalizedFrequencySff(double gamma, double time)
    {
        ValidatePositiveGamma(gamma);
        if (!double.IsFinite(time))
            throw new ArgumentOutOfRangeException(nameof(time), time, "time must be finite.");
        return 1.0;
    }

    /// <summary>The unnormalized frequency SFF of either exact endpoint:
    /// <c>(N+1)²</c> at every finite real time and strictly positive uniform dephasing.</summary>
    public double EndpointUnnormalizedFrequencySff(int N, double gamma, double time)
    {
        _ = EndpointNormalizedFrequencySff(gamma, time);
        double count = EndpointMultiplicity(N, gamma);
        return count * count;
    }

    public F43BandSffPairingPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        F1Pi2Inheritance f1)
        : base("F43 palindrome-paired decay-rate-band SFF: I reflects to 2Nγ−I with multiplicity preserved; at positive uniform γ, exact endpoint multiplicity N+1 and constant normalized frequency SFF 1",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F43 + " +
               "docs/proofs/derivations/D09_SECTOR_SFF_PAIRING.md + " +
               "experiments/SPECTRAL_FORM_FACTOR.md + " +
               "compute/RCPsiSquared.Core/Symmetry/F1Pi2Inheritance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        Ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        F1 = f1 ?? throw new ArgumentNullException(nameof(f1));
    }

    public override string DisplayName =>
        "F43 palindrome-paired band SFF as Pi2-Foundation a_0 + F1 inheritance";

    public override string Summary =>
        $"decay-rate band I ↔ 2Nγ−I with equal frequency SFF; average light reflects continuously as x ↔ N−x; at positive uniform γ, connected-chain endpoint multiplicity N+1 with constant normalized K=1; coefficient 2 = a_0 (= {ReflectionCoefficient}) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F43 reflected-band identity",
                summary: "the full-spectrum palindrome gives a multiplicity-preserving bijection between reflected decay-rate bands; frequency multisets differ only by sign, so their SFFs agree for every real time");
            yield return InspectableNode.RealScalar("ReflectionCoefficient (= a_0 = 2)", ReflectionCoefficient);
            yield return new InspectableNode("average-light reading",
                summary: "under uniform dephasing, d = 2γ⟨n_XY⟩ and the reflected band has average light N−⟨n_XY⟩; fractional average light is allowed because H mixes Pauli weights by ±2");
            yield return new InspectableNode("exact endpoint bands",
                summary: $"at strictly positive uniform γ, a connected chain has N+1 stationary modes and N+1 partners at −2Nγ; all endpoint frequencies vanish, so normalized K(t) = {EndpointNormalizedFrequencySff(gamma: 1.0, time: 1.0)} constantly and unnormalized K(t) = (N+1)²; γ=0 is excluded because the endpoints collapse into the larger commutator kernel");
            yield return new InspectableNode("finite-width producer boundary",
                summary: "a numerical rate window may include modes beyond the exact endpoint eigenspace; its sampled count and SFF are not promoted to the endpoint theorem");
        }
    }

    private static void ValidateN(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F43 requires N ≥ 1.");
    }

    private static void ValidateGamma(double gamma)
    {
        if (!double.IsFinite(gamma) || gamma < 0.0)
            throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "γ must be finite and ≥ 0.");
    }

    private static void ValidatePositiveGamma(double gamma)
    {
        if (!double.IsFinite(gamma) || gamma <= 0.0)
            throw new ArgumentOutOfRangeException(nameof(gamma), gamma,
                "exact F43 endpoint multiplicity requires finite γ > 0; at γ=0 the endpoints collapse into the commutator kernel.");
    }
}
