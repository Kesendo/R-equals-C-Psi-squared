using System.Collections.Generic;
using System.Numerics;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Live witness: recomputes the F89 path-3 octic EP-character at inspect time.
/// Contrast with the coherence-horizon √-EP (DEFECTIVE): this one is DIABOLIC (semisimple)
/// — eigenvalues coalesce, eigenvectors stay independent. The exact double discriminant
/// factor locates an analytic crossing after pair isolation; semisimplicity comes from the
/// twin-scalar compression and the live geometric-multiplicity/departure checks.</summary>
public sealed class F89OcticCharacterWitness : IInspectable
{
    private static double QEp => F89Path3OcticEpClaim.QEp;

    private static EpCharacter.Reading Read()
    {
        var l = F89Path3OcticBlock.BuildSeDeSymBlock(QEp, 1.0);
        return EpCharacter.Characterize(l, F89Path3OcticEpClaim.MergedEigenvalue(1.0, QEp), radius: 0.5);
    }

    public string DisplayName => "F89 path-3 octic diabolic degeneracy (EpCharacter, live)";

    public string Summary
    {
        get
        {
            var r = Read();
            return $"{r.Kind}: geo={r.Geometric}=alg={r.Algebraic}, dep={r.Departure:E2}, " +
                   $"‖P‖={r.ProjectorNorm:F2}, |cos|={r.EigenvectorMergeCos:F2} at λ_EP=−4γ+2i·q_EP";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var r = Read();
            yield return new InspectableNode("verdict",
                summary: $"{r.Kind} (geo={r.Geometric}=alg={r.Algebraic}, dep={r.Departure:E2}); " +
                         "DIABOLIC = eigenvalues coalesce, eigenvectors independent — NOT a defective EP");
            yield return InspectableNode.RealScalar("q_EP", QEp);
            yield return InspectableNode.RealScalar("projector-norm ‖P‖ (>1 ⟹ obliquely embedded / non-normal; the diabolic discriminators are geo=alg & dep≈0, not ‖P‖)", r.ProjectorNorm);
            yield return InspectableNode.RealScalar("departure-from-normality (≈0 ⟹ diabolic)", r.Departure);
            yield return new InspectableNode("grid-free anchor",
                summary: "disc(F_8) has (3q⁴+q²−1) to multiplicity 2, locating an analytic crossing after pair isolation; character comes from the twin-scalar restriction plus geo=alg and dep≈0, not from even order alone");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
