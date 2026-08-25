using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wiring of <see cref="SeatCutBlindnessClaim"/> (F157; Tier1Derived). One typed parent,
/// <see cref="F4KernelDimensionByComponentsClaim"/>, registered by <c>RegisterF1Family</c>:
/// F4's disconnected-graph extension counts the kernel with the watching on EVERY site, and its
/// open question is which seat carries the γ; this claim is the opposite corner, ONE seat in the
/// single-excitation sector, and it closes that question there. Resolution is topological, so this
/// registration may sit anywhere in the chain relative to the family that supplies the parent.</summary>
public static class SeatCutBlindnessClaimRegistration
{
    public static ClaimRegistryBuilder RegisterSeatCutBlindnessClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<SeatCutBlindnessClaim>(b =>
            new SeatCutBlindnessClaim(b.Get<F4KernelDimensionByComponentsClaim>()));
}
