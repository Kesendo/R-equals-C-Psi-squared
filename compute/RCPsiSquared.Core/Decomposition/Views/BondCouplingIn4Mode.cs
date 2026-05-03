using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Resonance;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Decomposition.Views;

/// <summary>F86 inspection wrapper for one bond's M_h_per_bond_eff[b]. Adds the per-bond
/// metadata (bond index, <see cref="BondClass"/>, total bond count) to the
/// <see cref="BlockMatrixIn4Mode"/> structural decomposition.
///
/// <para>The cross-block Frobenius is the bond-position-dependent number whose two-class
/// pattern (Endpoint vs Interior) is the F86 universal-shape fingerprint. Exposing it as
/// a typed property here means a future <c>BondCouplingDistribution</c> view can plot
/// CrossBlockFrobenius vs bond index without needing to re-derive it.</para>
/// </summary>
public sealed class BondCouplingIn4Mode : IInspectable
{
    public int BondIndex { get; }
    public int NumBonds { get; }
    public BondClass BondClass { get; }
    public BlockMatrixIn4Mode Decomposition { get; }

    public ComplexMatrix Matrix => Decomposition.Matrix;
    public double Frobenius => Decomposition.Frobenius;
    public double CrossBlockFrobenius => Decomposition.CrossBlockFrobenius;
    public double CrossFrobeniusFraction => Decomposition.CrossFrobeniusFraction;

    public BondCouplingIn4Mode(int bondIndex, int numBonds, ComplexMatrix mhPerBondEff)
    {
        if (bondIndex < 0 || bondIndex >= numBonds)
            throw new ArgumentOutOfRangeException(nameof(bondIndex), $"bondIndex {bondIndex} not in [0, {numBonds}).");
        BondIndex = bondIndex;
        NumBonds = numBonds;
        BondClass = (bondIndex == 0 || bondIndex == numBonds - 1)
            ? BondClass.Endpoint
            : BondClass.Interior;
        Decomposition = new BlockMatrixIn4Mode($"M_h_per_bond_eff[{bondIndex}]", mhPerBondEff);
    }

    public string DisplayName => $"bond {BondIndex} ({BondClass})";
    public string Summary =>
        $"‖M_b‖_F = {Frobenius:F4}, cross fraction = {CrossFrobeniusFraction:P1}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return InspectableNode.RealScalar("bond index", BondIndex);
            yield return new InspectableNode("bond class", summary: BondClass.ToString());
            foreach (var c in Decomposition.Children) yield return c;
        }
    }

    public InspectablePayload Payload => Decomposition.Payload;
}
