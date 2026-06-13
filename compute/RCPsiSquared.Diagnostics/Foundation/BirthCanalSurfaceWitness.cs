using System.Globalization;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The sterile↔birth-canal boundary as a live surface (Object-Manager root
/// <c>inspect --root surface</c>). The boundary lives in γ-profile space; this witness draws a 2D
/// slice of it (symmetric profiles by (w_edge, w_center), bulk solved by Σ=N) and reads every point
/// through six lenses plus the two-ray s* punchline. The C# home of the arc
/// <c>birth_canal_surface</c> and of <c>simulations/birth_canal_boundary_pathdependence.py</c>.
///
/// <para>Not a Claim: a live reading. The physics is the Absorption Theorem
/// (<see cref="RCPsiSquared.Core.Symmetry.AbsorptionTheoremClaim"/>): rate = 2·Σ_l γ_l·light_l, so
/// Deviation = rate(Q=1000) − rate(Q=1.5) reads the light-freeze. The reading is
/// initial-state-independent (it reads L's slow subspace, not a trajectory from ρ₀) — do not
/// "fix" it to depend on an initial state.</para>
///
/// <para>R1 (the two sterility kinds): Deviation = 0 is necessary but not sufficient for
/// light-freeze. At flat γ, rate = 2γ·(total light) is Q-invariant by uniformity alone, so the
/// probe reads sterile whether or not the light drifts. Genuine freeze (peaked-V, robust, breaks at
/// s*=0.709) vs flat-γ distribution-blindness (uniform, fragile, breaks at s*=0.105). The mechanism
/// and ray nodes tell the two kinds apart.</para>
///
/// <para>Reuses <see cref="PostEpFlowField"/> as the per-point engine, so the numbers are
/// bit-identical to the verified N=5 pins; the grid Deviation equals
/// <see cref="PostEpFlowField.BirthCanalDeviation"/> (pinned by test).</para></summary>
public sealed class BirthCanalSurfaceWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    public int N { get; }
    public int GridK { get; }

    // The N=5 admissible box: contains the verified anchors and the boundary curve (which runs
    // (w_edge, w_center) ≈ (0.25, 1.94) → (0.92, 1.05)). N=6 may need a different box (a Risk).
    private const double EdgeLo = 0.2, EdgeHi = 1.0, CenterLo = 0.5, CenterHi = 3.0;

    // M4: at the low probe Q an EP makes the eigenbasis singular and the light reading drift, which
    // shows up as an absorption-residual blow-up. Above this, mask the point as EP-suspect.
    private const double AbsorptionResidualThreshold = 1e-3;

    public BirthCanalSurfaceWitness(int n = 5, int grid = 9)
    {
        if (n < 5 || n > 6)
            throw new ArgumentOutOfRangeException(nameof(n), n,
                "the surface slice needs N in 5..6 (N=5 is the instrument; N=6 is a 4096^2 EVD per " +
                "point, coarse-grid only; N<5 has no bulk reservoir)");
        if (grid < 3) throw new ArgumentOutOfRangeException(nameof(grid), grid, "grid must be >= 3");
        N = n;
        GridK = grid;
    }

    /// <summary>One point of the surface: the coordinates, whether the bulk is positive, the
    /// Deviation = rate(Q_high) − rate(Q_low) (= <see cref="PostEpFlowField.BirthCanalDeviation"/>),
    /// the canal verdict, the EP-suspect flag (M4), the max per-site light drift between the probes
    /// (L2), and the two assembly readings (null when inadmissible).</summary>
    public sealed record SurfacePoint(
        double WEdge, double WCenter, bool Admissible,
        double Deviation, bool IsCanal, bool EpSuspect, double DriftMax,
        FlowAssemblyReading? Low, FlowAssemblyReading? High);

    /// <summary>Read one (w_edge, w_center) point through the per-point engine. Builds the profile,
    /// reads the assembly at both probe Q's via <see cref="PostEpFlowField.ReadAssembly"/>, and
    /// derives Deviation (L1), the light drift (L2), the absorption residual (L3 / M4 EP-flag).
    /// Inadmissible points return Admissible=false with null readings.</summary>
    public SurfacePoint ReadPoint(double wEdge, double wCenter)
    {
        if (!SymmetricGammaSlice.IsAdmissible(N, wEdge, wCenter))
            return new SurfacePoint(wEdge, wCenter, false, double.NaN, false, false, double.NaN, null, null);

        var profile = SymmetricGammaSlice.Profile(N, wEdge, wCenter);
        var field = new PostEpFlowField(N,
            new[] { PostEpFlowField.BirthCanalProbeQLow, PostEpFlowField.BirthCanalProbeQHigh },
            new[] { 0.0, 1.0 }, profile);
        var lo = field.ReadAssembly(PostEpFlowField.BirthCanalProbeQLow);
        var hi = field.ReadAssembly(PostEpFlowField.BirthCanalProbeQHigh);

        double dev = hi.SlowestRate - lo.SlowestRate;
        bool isCanal = Math.Abs(dev) > PostEpFlowField.BirthCanalTolerance;

        double drift = 0.0;
        for (int l = 0; l < N; l++)
            drift = Math.Max(drift, Math.Abs(hi.PerSiteLight[l] - lo.PerSiteLight[l]));

        // M4: the low-probe absorption residual is the EP detector (the high probe's residual is
        // dominated by ‖L‖ growth, not EP-conditioning, so it is not used here).
        bool epSuspect = Math.Abs(lo.SlowestRate - lo.AbsorptionRate) > AbsorptionResidualThreshold;

        return new SurfacePoint(wEdge, wCenter, true, dev, isCanal, epSuspect, drift, lo, hi);
    }

    public string DisplayName => $"BirthCanalSurfaceWitness (N={N}, chain, {GridK}x{GridK} grid)";

    public string Summary =>
        "the sterile<->birth-canal boundary as a live 2D slice of gamma-profile space " +
        "(symmetric profiles by w_edge x w_center, bulk solved by sum=N): the whole surface, not the " +
        "s*=0.709 line. Membership = PostEpFlowField.BirthCanalDeviation; mechanism = Absorption " +
        "Theorem (rate = 2*sum gamma_l*light_l). Two sterility kinds: genuine freeze (robust) vs " +
        "flat-gamma blindness (fragile).";

    public IEnumerable<IInspectable> Children => Array.Empty<IInspectable>();   // nodes added in Tasks 3-5

    public InspectablePayload Payload => InspectablePayload.Empty;
}
