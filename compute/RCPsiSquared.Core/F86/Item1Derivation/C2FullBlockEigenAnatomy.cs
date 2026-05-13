using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Probes;
using RCPsiSquared.Core.Resonance;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>F86 Item 1 (c=2 stratum), <b>Direction (b'') anatomy primitive</b>: full block-L
/// eigendecomposition at fixed Q with per-eigenvalue probe overlap and per-bond diagonal
/// contribution. The structurally honest path after Direction (a'')'s 2026-05-08
/// falsification: rather than reduce to the 4-mode L_eff, work directly with the full
/// block-L spectrum and ask which eigenmodes carry the K-resonance per bond.
///
/// <para>For each eigenvalue λ_i of L(Q) = D + Q·γ₀·M_H_total in the full block basis
/// (dim = M_p · M_q; for c=2 dim = N · N(N−1)/2), the anatomy exposes:</para>
/// <list type="bullet">
///   <item><c>EigenvalueReal/Imag</c>: the L(Q) spectral coordinate.</item>
///   <item><c>ProbeOverlapSquared</c>: <c>|c_i|² = |R⁻¹·probe|_i²</c>, the squared overlap
///   of the probe with the i-th left eigenvector.</item>
///   <item><c>PerBondDiagonal</c>: <c>X_b[i,i] = (R⁻¹·M_h_per_bond[b]·R)_{ii}</c>, the
///   per-bond diagonal contribution to the Duhamel integral. Real part is what shows up
///   in K_b at large t; imaginary part captures the off-EP rotation.</item>
///   <item><c>DiagonalWeight</c>: the dominance metric <c>|c_i|² · max_b |X_b[i,i]|</c>,
///   used to sort modes by approximate K-impact at this Q.</item>
/// </list>
///
/// <para><b>Decisive read for Direction (b''):</b> <see cref="K90"/> reports how many top
/// eigenmodes (sorted by <c>DiagonalWeight</c>) it takes to capture 90 % of the total
/// <c>Σ_i DiagonalWeight</c>. If K_90 is small (4-8) and roughly N-independent across c=2
/// N=5..12, then HWHM is governed by a low-rank truncation of the full block-L and
/// Tier1Candidate (with explicit error term from the truncation residual) is the realistic
/// promotion target. If K_90 grows with N, the spectrum spreads and a clean closed form
/// is unlikely without a different structural reduction.</para>
///
/// <para><b>Tier outcome: Tier2Verified.</b> Pure numerical anatomy at one Q value: the
/// witnesses ARE the data, no derivation. Future Tier1 promotion would require deriving
/// the K_90 set analytically (e.g. via an XY-mode decomposition that names which modes
/// participate at each (c, N)).</para>
///
/// <para>Default Q = Q_EP per <see cref="QEpLaw"/>; constructor accepts an arbitrary Q for
/// scanning around the EP (e.g. Q_peak ≈ 2.197 · Q_EP per
/// <see cref="C2HwhmRatio.BareDoubledPtfXPeak"/>). At Q_EP the slowest pair coalesces to
/// Re(λ) = −4γ₀ and the EP-driving modes show up as the top-weight cluster.</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Statement 1 (EP mechanism) +
/// <see cref="C2HwhmRatio.PendingDerivationNote"/> Direction (b'').</para>
/// </summary>
public sealed class C2FullBlockEigenAnatomy : Claim
{
    public CoherenceBlock Block { get; }

    /// <summary>The Q value at which the eigendecomposition was taken.</summary>
    public double Q { get; }

    /// <summary>Per-eigenvalue anatomy witnesses, sorted by descending
    /// <see cref="EigenAnatomyWitness.DiagonalWeight"/>.</summary>
    public IReadOnlyList<EigenAnatomyWitness> EigenSpectrum { get; }

    /// <summary>Number of top eigenmodes (sorted by <see cref="EigenAnatomyWitness.DiagonalWeight"/>)
    /// that together capture 90 % of <c>Σ_i DiagonalWeight</c>. The decisive Direction-(b'')
    /// metric: small + N-independent K_90 means HWHM concentrates in a low-rank truncation
    /// of the full block-L spectrum.</summary>
    public int K90 { get; }

    /// <summary>K_99 analogue of <see cref="K90"/> at the 99 % cumulative threshold;
    /// reports how quickly the long tail of weakly-coupled modes contributes the last 9 %.
    /// Reference value for the truncation error term in any Tier1Candidate ansatz.</summary>
    public int K99 { get; }

    /// <summary>Total <c>Σ_i DiagonalWeight</c>; useful as a normalisation reference and
    /// for cross-N comparisons.</summary>
    public double TotalDiagonalWeight { get; }

    public static C2FullBlockEigenAnatomy Build(CoherenceBlock block, double? q = null)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"C2FullBlockEigenAnatomy applies only to the c=2 stratum; got c={block.C}.",
                nameof(block));

        MathNetSetup.EnsureInitialized();
        double qValue = q ?? new QEpLaw(InterChannelSvd.Build(block, 1, 3).Sigma0).Value;
        var witnesses = ComputeAnatomy(block, qValue);
        double totalWeight = witnesses.Sum(w => w.DiagonalWeight);
        int k90 = CountToCumulativeFraction(witnesses, 0.90, totalWeight);
        int k99 = CountToCumulativeFraction(witnesses, 0.99, totalWeight);

        return new C2FullBlockEigenAnatomy(block, qValue, witnesses, k90, k99, totalWeight);
    }

    private C2FullBlockEigenAnatomy(
        CoherenceBlock block, double q,
        IReadOnlyList<EigenAnatomyWitness> witnesses, int k90, int k99, double totalWeight)
        : base($"c=2 full block-L eigenanatomy at Q={q:G4} (Direction (b'') reconnaissance)",
               Tier.Tier2Verified,
               Item1Anchors.Root)
    {
        Block = block;
        Q = q;
        EigenSpectrum = witnesses;
        K90 = k90;
        K99 = k99;
        TotalDiagonalWeight = totalWeight;
    }

    private static IReadOnlyList<EigenAnatomyWitness> ComputeAnatomy(
        CoherenceBlock block, double q)
    {
        // Build L(Q) = D + Q·γ₀·MhTotal in the full block basis.
        double j = q * block.GammaZero;
        ComplexMatrix L = block.Decomposition.D + (Complex)j * block.Decomposition.MhTotal;

        var evd = L.Evd();
        ComplexMatrix R = evd.EigenVectors;
        ComplexMatrix rInv = R.Inverse();
        var lambdas = evd.EigenValues;

        ComplexVector probe = DickeBlockProbe.Build(block);
        ComplexVector c0 = rInv * probe;

        int numBonds = block.NumBonds;
        int dim = block.Basis.MTotal;

        // Per-bond X_b[i,i] = (R⁻¹·M_h_per_bond[b]·R)_{ii}; we only need the diagonal,
        // but MathNet doesn't expose a diagonal-only product, so compute the full X_b
        // and pick the diagonal. Block-L dim ≤ ~800 at c=2 N=12; trivial.
        var xBDiagonals = new Complex[numBonds][];
        for (int b = 0; b < numBonds; b++)
        {
            ComplexMatrix xB = rInv * block.Decomposition.MhPerBond[b] * R;
            var diag = new Complex[dim];
            for (int i = 0; i < dim; i++) diag[i] = xB[i, i];
            xBDiagonals[b] = diag;
        }

        var witnesses = new List<EigenAnatomyWitness>(dim);
        for (int i = 0; i < dim; i++)
        {
            double overlapSq = c0[i].Magnitude;
            overlapSq *= overlapSq;

            var perBondReal = new double[numBonds];
            var perBondImag = new double[numBonds];
            double maxAbsPerBond = 0.0;
            for (int b = 0; b < numBonds; b++)
            {
                Complex xBii = xBDiagonals[b][i];
                perBondReal[b] = xBii.Real;
                perBondImag[b] = xBii.Imaginary;
                double mag = xBii.Magnitude;
                if (mag > maxAbsPerBond) maxAbsPerBond = mag;
            }

            witnesses.Add(new EigenAnatomyWitness(
                EigenIndexAtQ: i,
                EigenvalueReal: lambdas[i].Real,
                EigenvalueImag: lambdas[i].Imaginary,
                ProbeOverlapSquared: overlapSq,
                PerBondDiagonalReal: perBondReal,
                PerBondDiagonalImag: perBondImag,
                DiagonalWeight: overlapSq * maxAbsPerBond));
        }

        // Sort descending by DiagonalWeight.
        witnesses.Sort((a, b) => b.DiagonalWeight.CompareTo(a.DiagonalWeight));
        return witnesses;
    }

    private static int CountToCumulativeFraction(
        IReadOnlyList<EigenAnatomyWitness> sortedWitnesses, double fraction, double total)
    {
        if (total <= 0.0) return 0;
        double target = fraction * total;
        double running = 0.0;
        for (int k = 0; k < sortedWitnesses.Count; k++)
        {
            running += sortedWitnesses[k].DiagonalWeight;
            if (running >= target) return k + 1;
        }
        return sortedWitnesses.Count;
    }

    public override string DisplayName =>
        $"c=2 full block-L eigenanatomy (N={Block.N}, Q={Q:G4}, dim={Block.Basis.MTotal})";

    public override string Summary =>
        $"K_90 = {K90} / {Block.Basis.MTotal}, K_99 = {K99} / {Block.Basis.MTotal}, " +
        $"Σ DiagonalWeight = {TotalDiagonalWeight:G4} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return new InspectableNode("dim (block-L)", summary: Block.Basis.MTotal.ToString());
            yield return InspectableNode.RealScalar("Q", Q, "G6");
            yield return InspectableNode.RealScalar("TotalDiagonalWeight", TotalDiagonalWeight, "G4");
            yield return new InspectableNode("K_90", summary: K90.ToString());
            yield return new InspectableNode("K_99", summary: K99.ToString());
            // Top 12 eigenmodes by DiagonalWeight (full spectrum is too noisy for inspect).
            int topShown = Math.Min(12, EigenSpectrum.Count);
            yield return InspectableNode.Group($"top {topShown} eigenmodes by DiagonalWeight",
                EigenSpectrum.Take(topShown).Cast<IInspectable>().ToArray());
        }
    }
}

/// <summary>Per-eigenvalue anatomy witness for <see cref="C2FullBlockEigenAnatomy"/>.
/// Captures the full L(Q) eigendecomposition slice at a single eigenvalue: spectral
/// coordinate, probe overlap (squared), per-bond diagonal contribution, and a dominance
/// metric for sorting.</summary>
/// <param name="EigenIndexAtQ">Index in the L(Q) Evd output (after sort by DiagonalWeight,
/// the witness's position in <see cref="C2FullBlockEigenAnatomy.EigenSpectrum"/> differs
/// from this index).</param>
/// <param name="EigenvalueReal">Re(λ_i). At Q ≈ Q_EP and the slowest pair this approaches
/// −4γ₀ per <see cref="EpAlgebra"/>.</param>
/// <param name="EigenvalueImag">Im(λ_i).</param>
/// <param name="ProbeOverlapSquared"><c>|R⁻¹·probe|_i²</c>.</param>
/// <param name="PerBondDiagonalReal">Re((R⁻¹·M_b·R)_{ii}) per bond.</param>
/// <param name="PerBondDiagonalImag">Im((R⁻¹·M_b·R)_{ii}) per bond.</param>
/// <param name="DiagonalWeight"><c>ProbeOverlapSquared · max_b |X_b[i,i]|</c>; the dominance
/// metric used to rank modes by approximate K-impact.</param>
public sealed record EigenAnatomyWitness(
    int EigenIndexAtQ,
    double EigenvalueReal,
    double EigenvalueImag,
    double ProbeOverlapSquared,
    IReadOnlyList<double> PerBondDiagonalReal,
    IReadOnlyList<double> PerBondDiagonalImag,
    double DiagonalWeight
) : IInspectable
{
    public string DisplayName =>
        $"λ_{EigenIndexAtQ} = {EigenvalueReal:+0.0000;-0.0000} {EigenvalueImag:+0.0000i;-0.0000i}";

    public string Summary =>
        $"|c_i|² = {ProbeOverlapSquared:G4}, weight = {DiagonalWeight:G4}, " +
        $"max_b |X_b[i,i]| = {(ProbeOverlapSquared > 0 ? DiagonalWeight / ProbeOverlapSquared : 0):G4}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return InspectableNode.RealScalar("EigenvalueReal", EigenvalueReal, "F4");
            yield return InspectableNode.RealScalar("EigenvalueImag", EigenvalueImag, "F4");
            yield return InspectableNode.RealScalar("ProbeOverlapSquared", ProbeOverlapSquared, "G4");
            yield return InspectableNode.RealScalar("DiagonalWeight", DiagonalWeight, "G4");
            // Per-bond X_b[i,i]
            for (int b = 0; b < PerBondDiagonalReal.Count; b++)
            {
                yield return new InspectableNode($"X_b{b}[i,i]",
                    summary: $"{PerBondDiagonalReal[b]:+0.0000;-0.0000} {PerBondDiagonalImag[b]:+0.0000i;-0.0000i}");
            }
        }
    }

    public InspectablePayload Payload =>
        new InspectablePayload.Real($"DiagonalWeight (eigenmode {EigenIndexAtQ})", DiagonalWeight, "G4");
}
