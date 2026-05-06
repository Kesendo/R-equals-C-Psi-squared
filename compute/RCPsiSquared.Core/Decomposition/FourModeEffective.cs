using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition.Views;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Probes;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Decomposition;

/// <summary>The 4-mode effective model for the F86 K_CC_pr Q-scan: all relevant operators
/// projected onto the orthonormal basis B = [|c_1⟩, |c_3⟩, |u_0⟩, |v_0⟩] from
/// <see cref="FourModeBasis"/>. Each projection is X_eff = B† · X · B.
///
/// <para>Captures the math note from the F86 structural exploration:</para>
/// <list type="bullet">
///   <item>Channel-uniform 2D subspace {|c_1⟩, |c_3⟩} = where the Dicke probe lives.</item>
///   <item>SVD-top 2D subspace {|u_0⟩, |v_0⟩} = where the EP coupling lives (under H).</item>
///   <item>The two are orthogonal in the full block; the 4-mode basis is their direct sum.</item>
///   <item>Per-bond V_b = B† · M_H_per_bond[b] · B carries the bond-position dependence
///         that splits Endpoint vs Interior shapes.</item>
/// </list>
///
/// <para><b>Structural design constraint:</b> at uniform J, L_eff(Q) = D_eff + Q·γ₀·MhTotalEff
/// is <b>bond-summed by design</b> (see <see cref="LEffAtQ"/>). Per-bond information enters
/// ONLY through the per-bond <see cref="MhPerBondEff"/>[b] in dL/dJ_b (used by K_b(Q, t)
/// Duhamel evaluation in <c>FourModeResonanceScan</c>), NOT through the L_eff spectrum or
/// eigenstates. Any direction that hopes to extract bond-class signature from the 4-mode
/// L_eff EIGENVALUES or EIGENSTATES alone (e.g. Direction (α) polarity-Bloch projection at
/// t_peak) is structurally tautological under this reduction. The bond-class signature
/// lives in the K-resonance via per-bond V_b, not in the spectrum.</para>
///
/// <para>The L_eff(Q) = D_eff + J·Σ_b V_b_eff is a 4×4 matrix; eigendecomposition is
/// instant, Duhamel integrals are tiny. The comparison <c>FourModeResonanceScan</c> vs
/// <c>ResonanceScan</c> (full block-L) reveals whether the 4-mode model captures the
/// observed universal shape — see PROOF_F86_QPEAK.md "What's missing for full Tier 1".</para>
/// </summary>
public sealed class FourModeEffective : IInspectable
{
    public CoherenceBlock Block { get; }
    public FourModeBasis Basis { get; }
    public ComplexMatrix DEff { get; }
    public IReadOnlyList<ComplexMatrix> MhPerBondEff { get; }
    public ComplexMatrix MhTotalEff { get; }
    public ComplexVector ProbeEff { get; }
    public ComplexMatrix SKernelEff { get; }

    private readonly Lazy<DiagonalRatesIn4Mode> _dEffView;
    private readonly Lazy<BlockMatrixIn4Mode> _mhTotalView;
    private readonly Lazy<IReadOnlyList<BondCouplingIn4Mode>> _mhPerBondViews;
    private readonly Lazy<ProjectedSubspaceVector> _probeView;
    private readonly Lazy<BlockMatrixIn4Mode> _sKernelView;

    /// <summary>D_eff as a typed view: per-mode rates, off-diagonal residual, trace.</summary>
    public DiagonalRatesIn4Mode DEffView => _dEffView.Value;

    /// <summary>M_h_total_eff decomposed into channel-uniform / SVD-top / cross 2×2 blocks.</summary>
    public BlockMatrixIn4Mode MhTotalView => _mhTotalView.Value;

    /// <summary>Per-bond M_h with bond-class metadata and the cross-block fingerprint that
    /// drives the F86 Endpoint vs Interior shape split.</summary>
    public IReadOnlyList<BondCouplingIn4Mode> MhPerBondViews => _mhPerBondViews.Value;

    /// <summary>Probe with channel-uniform / SVD-top fraction split (SVD-top = 0 structurally).</summary>
    public ProjectedSubspaceVector ProbeView => _probeView.Value;

    /// <summary>S_kernel decomposed in the same 4-mode block structure.</summary>
    public BlockMatrixIn4Mode SKernelView => _sKernelView.Value;

    private FourModeEffective(CoherenceBlock block, FourModeBasis basis,
        ComplexMatrix dEff, IReadOnlyList<ComplexMatrix> mhPerBondEff, ComplexMatrix mhTotalEff,
        ComplexVector probeEff, ComplexMatrix sKernelEff)
    {
        Block = block;
        Basis = basis;
        DEff = dEff;
        MhPerBondEff = mhPerBondEff;
        MhTotalEff = mhTotalEff;
        ProbeEff = probeEff;
        SKernelEff = sKernelEff;

        _dEffView = new Lazy<DiagonalRatesIn4Mode>(() => new DiagonalRatesIn4Mode(dEff));
        _mhTotalView = new Lazy<BlockMatrixIn4Mode>(() => new BlockMatrixIn4Mode("M_h_total_eff", mhTotalEff));
        _mhPerBondViews = new Lazy<IReadOnlyList<BondCouplingIn4Mode>>(() =>
        {
            int numBonds = mhPerBondEff.Count;
            var list = new BondCouplingIn4Mode[numBonds];
            for (int b = 0; b < numBonds; b++)
                list[b] = new BondCouplingIn4Mode(b, numBonds, mhPerBondEff[b]);
            return list;
        });
        _probeView = new Lazy<ProjectedSubspaceVector>(() => new ProjectedSubspaceVector("probe_eff", probeEff));
        _sKernelView = new Lazy<BlockMatrixIn4Mode>(() => new BlockMatrixIn4Mode("S_kernel_eff", sKernelEff));
    }

    public string DisplayName =>
        $"FourModeEffective (c={Block.C}, N={Block.N}, n={Block.LowerPopcount}, γ₀={Block.GammaZero:G3})";

    public string Summary =>
        $"basis ⊥-residual = {Basis.OffOrthonormalityResidual:E2}, " +
        $"‖M_h_total‖_F = {MhTotalView.Frobenius:F4}, " +
        $"probe → CU = {ProbeView.ChannelUniformFraction:P0}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("Block (CoherenceBlock)",
                summary: $"N={Block.N}, n={Block.LowerPopcount}, c={Block.C}, γ₀={Block.GammaZero:G3}");
            yield return new InspectableNode("Basis (FourModeBasis)",
                summary: $"HD={Basis.Hd1}/{Basis.Hd2}, ⊥-residual={Basis.OffOrthonormalityResidual:E2}");
            yield return DEffView;
            yield return MhTotalView;
            yield return InspectableNode.Group("M_h_per_bond_eff", MhPerBondViews, MhPerBondViews.Count);
            yield return ProbeView;
            yield return SKernelView;
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    public static FourModeEffective Build(CoherenceBlock block, int hd1 = 1, int hd2 = 3)
    {
        var basis = FourModeBasis.Build(block, hd1, hd2);
        var decomp = block.Decomposition;
        var probe = DickeBlockProbe.Build(block);
        var sKernel = SpatialSumKernel.Build(block);

        var dEff = basis.Project(decomp.D);
        var mhPerBondEff = decomp.MhPerBond.Select(m => basis.Project(m)).ToList();
        var mhTotalEff = basis.Project(decomp.MhTotal);
        var probeEff = basis.Project(probe);
        var sKernelEff = basis.Project(sKernel);

        return new FourModeEffective(block, basis, dEff, mhPerBondEff, mhTotalEff, probeEff, sKernelEff);
    }

    /// <summary>L_eff(Q) = D_eff + (Q·γ₀)·M_H_total_eff. Uniform-J across all bonds.</summary>
    public ComplexMatrix LEffAtQ(double q)
    {
        double j = q * Block.GammaZero;
        return DEff + (Complex)j * MhTotalEff;
    }
}
