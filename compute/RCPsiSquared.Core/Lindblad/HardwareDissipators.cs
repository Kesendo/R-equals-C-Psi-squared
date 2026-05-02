using System.Numerics;

namespace RCPsiSquared.Core.Lindblad;

/// <summary>Catalog of hardware-relevant dissipator classes with their Pauli decompositions
/// and per-class palindrome-residual constants (verified at N = 3, σ_offset = 0).
///
/// Each entry's contribution to ‖M‖² with rate γ_K_l per site l is
///   4^(N−1) · [c1 · Σ γ² + c2 · (Σ γ)²].
/// </summary>
public sealed record HardwareDissipator(string Name, Complex Alpha, Complex Beta, Complex Delta,
    double C1, double C2, string Description)
{
    public override string ToString() => $"{Name}: α={Alpha} β={Beta} δ={Delta}, c1={C1} c2={C2} — {Description}";
}

public static class HardwareDissipators
{
    public static readonly HardwareDissipator T1 = new(
        "T1", Alpha: 0.5, Beta: new Complex(0, -0.5), Delta: 0.0,
        C1: 3.0, C2: 4.0,
        Description: "amplitude relaxation σ⁻ = (X − iY)/2");

    public static readonly HardwareDissipator T1Pump = new(
        "T1pump", Alpha: 0.5, Beta: new Complex(0, 0.5), Delta: 0.0,
        C1: 3.0, C2: 4.0,
        Description: "amplitude pumping σ⁺ = (X + iY)/2 (thermal excitation)");

    public static readonly HardwareDissipator Tphi = new(
        "Tphi", Alpha: 0.0, Beta: 0.0, Delta: 1.0,
        C1: 0.0, C2: 16.0,
        Description: "pure dephasing σ_z (Π-respecting with σ_offset=γ; c1=c2=0 there, but c2=16 with σ_offset=0)");

    public static readonly HardwareDissipator XNoise = new(
        "Xnoise", Alpha: 1.0, Beta: 0.0, Delta: 0.0,
        C1: 16.0, C2: 16.0,
        Description: "X-axis noise / cross-talk σ_x");

    public static readonly HardwareDissipator YNoise = new(
        "Ynoise", Alpha: 0.0, Beta: 1.0, Delta: 0.0,
        C1: 0.0, C2: 16.0,
        Description: "Y-axis noise σ_y (Π-respecting with σ_offset=γ)");

    public static IReadOnlyList<HardwareDissipator> All { get; } = new[] { T1, T1Pump, Tphi, XNoise, YNoise };
}
