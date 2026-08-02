using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Spectrum;

/// <summary>Closed-form dispersion of the (0,1) coherence block for an isotropic
/// Heisenberg chain (Δ = 1) of length N under uniform per-site Z-dephasing γ. The object is
/// the N-dimensional block spanned by the |0⟩⟨j| between the ferromagnet and the single
/// excitations: a single-magnon nearest-neighbour tight-binding problem with hopping 2J,
/// giving N−1 oscillation modes plus the universal palindromic decay 2γ. It sits INSIDE the
/// Pauli w=1 sector (every string with exactly one site carrying X or Y, 2N·2^(N−1)-dimensional,
/// 160 against 5 at N=5) and is far smaller than it; that sector is not L-invariant and carries
/// no spectrum of its own, so the w=1 in this type's name is historical. D10 Step 6 has the scope.
/// (Δ = 1 is required: on this block the ZZ term is exactly the −2J·deg(j) diagonal that
/// turns the adjacency matrix into the Laplacian, so at anisotropic XXZ Δ ≠ 1 the clean
/// ω_k = 4J(1 − cos(πk/N)) band does NOT survive and the formula is Heisenberg-specific.)
///
/// <list type="bullet">
///   <item><see cref="Frequencies"/>: ω_k = 4J·(1 − cos(πk/N)) for k = 1..N−1.</item>
///   <item><see cref="UniformDecayRate"/>: 2γ, k-independent, the structural constant that
///         anchors the JW track's Lorentzian-width Γ in <c>JwBondQPeakUnified</c>.</item>
///   <item><see cref="QFactors"/>: Q_k = ω_k / (2γ).</item>
///   <item><see cref="MaxQ"/>, <see cref="MinQ"/>, <see cref="MeanQ"/>, <see cref="QSpread"/>,
///         <see cref="Bandwidth"/>: closed-form algebraic consequences.</item>
/// </list>
///
/// <para>Tier 1 derived: the dispersion is proved analytically in
/// <c>docs/proofs/derivations/D10_W1_DISPERSION.md</c> via reduction to the tight-binding
/// chain; verified for N=2..6 against a full Liouvillian eigendecomposition (15/15
/// frequencies, zero error) in <c>experiments/ANALYTICAL_SPECTRUM.md</c>, and entry-wise
/// through N=10 by block closure (<c>simulations/d10_block_closure_verify.py</c>). The N≥3 ctor pre-condition rejects N≤1 (no
/// modes) and N≥0 sanity. J=0 is allowed (gives flat zero spectrum, modes exist as labels);
/// γ&gt;0 is required because Q-factors and decay rate divide by γ.</para>
///
/// <para>Anchors: <c>experiments/ANALYTICAL_SPECTRUM.md</c>,
/// <c>docs/proofs/derivations/D10_W1_DISPERSION.md</c>.</para>
/// </summary>
public sealed class W1Dispersion : Claim
{
    public const string AnchorPath =
        "experiments/ANALYTICAL_SPECTRUM.md + docs/proofs/derivations/D10_W1_DISPERSION.md";

    public int N { get; }
    public double J { get; }
    public double GammaZero { get; }

    /// <summary>The N−1 mode frequencies ω_k = 4J·(1 − cos(πk/N)), k = 1..N−1.</summary>
    public IReadOnlyList<double> Frequencies { get; }

    /// <summary>Per-mode quality factor Q_k = ω_k / (2γ). All modes of the block share the same
    /// <see cref="UniformDecayRate"/>; the Q-spread is purely from the frequency dispersion.</summary>
    public IReadOnlyList<double> QFactors { get; }

    /// <summary>Universal decay rate on the block (= 2γ) at uniform γ; k-independent. This is the anchor
    /// for the Lorentzian-width Γ used by <c>JwBondQPeakUnified</c>'s NEW-NEW regime.</summary>
    public double UniformDecayRate { get; }

    /// <summary>Q at k = N−1: 2J/γ · (1 + cos(π/N)).</summary>
    public double MaxQ { get; }

    /// <summary>Q at k = 1: 2J/γ · (1 − cos(π/N)). Approaches Jπ²/(γN²) for large N, since 1 − cos(π/N) → π²/(2N²).</summary>
    public double MinQ { get; }

    /// <summary>Mean Q across the N−1 modes: exactly 2J/γ (because Σ_{k=1}^{N−1} cos(πk/N) = 0).</summary>
    public double MeanQ { get; }

    /// <summary>Q-factor spread: cot²(π/(2N)) = MaxQ / MinQ.</summary>
    public double QSpread { get; }

    /// <summary>Frequency bandwidth ω_{N−1} − ω_1 = 8J·cos(π/N); saturates at 8J as N → ∞.</summary>
    public double Bandwidth { get; }

    public W1Dispersion(int N, double J, double gammaZero)
        : base("(0,1) coherence block dispersion", Tier.Tier1Derived, AnchorPath)
    {
        if (N < 2)
            throw new ArgumentOutOfRangeException(
                nameof(N), N, "W1Dispersion requires N ≥ 2 (the block carries no oscillating mode at N=1).");
        if (gammaZero <= 0)
            throw new ArgumentOutOfRangeException(
                nameof(gammaZero), gammaZero,
                "W1Dispersion requires gammaZero > 0 (Q-factors and decay rate divide by γ).");

        this.N = N;
        this.J = J;
        this.GammaZero = gammaZero;
        UniformDecayRate = 2.0 * gammaZero;

        var freqs = new double[N - 1];
        var qs = new double[N - 1];
        for (int k = 1; k <= N - 1; k++)
        {
            double omega = 4.0 * J * (1.0 - Math.Cos(Math.PI * k / N));
            freqs[k - 1] = omega;
            qs[k - 1] = omega / UniformDecayRate;
        }
        Frequencies = freqs;
        QFactors = qs;

        double cosPiOverN = Math.Cos(Math.PI / N);
        MaxQ = 2.0 * J / gammaZero * (1.0 + cosPiOverN);
        MinQ = 2.0 * J / gammaZero * (1.0 - cosPiOverN);
        MeanQ = 2.0 * J / gammaZero;
        double cot = Math.Cos(Math.PI / (2.0 * N)) / Math.Sin(Math.PI / (2.0 * N));
        QSpread = cot * cot;
        Bandwidth = 8.0 * J * cosPiOverN;
    }

    public override string DisplayName =>
        $"(0,1) block dispersion (N={N}, J={J}, γ={GammaZero})";

    public override string Summary =>
        $"{N - 1} modes, decay {UniformDecayRate} (= 2γ), Q ∈ [{MinQ:F4}, {MaxQ:F4}], " +
        $"spread cot²(π/(2N)) = {QSpread:F4} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("N", N);
            yield return InspectableNode.RealScalar("J", J);
            yield return InspectableNode.RealScalar("gamma_zero", GammaZero);
            yield return InspectableNode.RealScalar("UniformDecayRate (= 2γ)", UniformDecayRate);
            yield return InspectableNode.RealScalar("MaxQ", MaxQ);
            yield return InspectableNode.RealScalar("MinQ", MinQ);
            yield return InspectableNode.RealScalar("MeanQ (= 2J/γ)", MeanQ);
            yield return InspectableNode.RealScalar("QSpread (= cot²(π/(2N)))", QSpread);
            yield return InspectableNode.RealScalar("Bandwidth (= 8J·cos(π/N))", Bandwidth);
            for (int k = 0; k < Frequencies.Count; k++)
                yield return new InspectableNode(
                    $"mode k={k + 1}",
                    summary: $"ω = {Frequencies[k]:F6}, Q = {QFactors[k]:F6}");
        }
    }
}
