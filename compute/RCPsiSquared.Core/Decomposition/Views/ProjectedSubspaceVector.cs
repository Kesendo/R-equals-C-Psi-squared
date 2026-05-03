using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Decomposition.Views;

/// <summary>F86 inspection wrapper for a 4-mode-projected vector (probe, S_kernel diagonal, …).
/// Computes the support fraction in each subspace (channel-uniform vs SVD-top), the L2 norm,
/// and the maximum component magnitude.
///
/// <para>Algebraic content: the channel-uniform / SVD-top split is the structural fact used by
/// <c>FourModeEffectiveTests.Build_ProbeEff_LivesInChannelUniformPair</c>. Exposing
/// <see cref="SvdTopFraction"/> as a typed property makes the "probe is exactly perpendicular
/// to the EP partners" claim self-documenting.</para>
/// </summary>
public sealed class ProjectedSubspaceVector : IInspectable
{
    public string Label { get; }
    public ComplexVector Vector { get; }
    public Complex OnC1 => Vector[0];
    public Complex OnC3 => Vector[1];
    public Complex OnU0 => Vector[2];
    public Complex OnV0 => Vector[3];

    private readonly Lazy<double> _l2Norm;
    private readonly Lazy<double> _channelUniformNorm;
    private readonly Lazy<double> _svdTopNorm;

    public double L2Norm => _l2Norm.Value;
    public double ChannelUniformNorm => _channelUniformNorm.Value;
    public double SvdTopNorm => _svdTopNorm.Value;
    public double ChannelUniformFraction => _l2Norm.Value > 0
        ? (_channelUniformNorm.Value * _channelUniformNorm.Value) / (_l2Norm.Value * _l2Norm.Value)
        : 0;
    public double SvdTopFraction => _l2Norm.Value > 0
        ? (_svdTopNorm.Value * _svdTopNorm.Value) / (_l2Norm.Value * _l2Norm.Value)
        : 0;

    public ProjectedSubspaceVector(string label, ComplexVector vector)
    {
        if (vector.Count != 4)
            throw new ArgumentException($"ProjectedSubspaceVector expects 4 components; got {vector.Count}.");
        Label = label;
        Vector = vector;
        _l2Norm = new Lazy<double>(() => vector.L2Norm());
        _channelUniformNorm = new Lazy<double>(() => vector.SubVector(0, 2).L2Norm());
        _svdTopNorm = new Lazy<double>(() => vector.SubVector(2, 2).L2Norm());
    }

    public string DisplayName => Label;
    public string Summary =>
        $"‖v‖₂ = {L2Norm:F4}, channel-uniform = {ChannelUniformFraction:P1}, SVD-top = {SvdTopFraction:P1}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode($"on {FourModeNames.C1}",
                summary: $"|·| = {OnC1.Magnitude:F4}",
                payload: new InspectablePayload.Scalar(FourModeNames.C1, OnC1));
            yield return new InspectableNode($"on {FourModeNames.C3}",
                summary: $"|·| = {OnC3.Magnitude:F4}",
                payload: new InspectablePayload.Scalar(FourModeNames.C3, OnC3));
            yield return new InspectableNode($"on {FourModeNames.U0}",
                summary: $"|·| = {OnU0.Magnitude:E2}",
                payload: new InspectablePayload.Scalar(FourModeNames.U0, OnU0));
            yield return new InspectableNode($"on {FourModeNames.V0}",
                summary: $"|·| = {OnV0.Magnitude:E2}",
                payload: new InspectablePayload.Scalar(FourModeNames.V0, OnV0));
            yield return InspectableNode.RealScalar("‖v‖₂", L2Norm, "F6");
            yield return InspectableNode.RealScalar("channel-uniform fraction", ChannelUniformFraction, "P3");
            yield return InspectableNode.RealScalar("SVD-top fraction", SvdTopFraction, "P3");
        }
    }

    public InspectablePayload Payload =>
        new InspectablePayload.Vector(Label, Vector, FourModeNames.All);
}
