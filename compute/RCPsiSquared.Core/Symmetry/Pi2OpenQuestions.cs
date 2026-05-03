using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Open theoretical items for the Π² Klein layer.</summary>
public static class Pi2OpenQuestions
{
    private const string Anchor = "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBase + Schicht-3 spectral observations 2026-05-03";

    public static IReadOnlyList<OpenQuestion> Standard { get; } = new[]
    {
        new OpenQuestion(
            "X-axis-flip mechanism for Marrakesh f83 fingerprint observables",
            "Empirical pattern (F88, locked in Pi2KleinHardwareViewTests): each H-class diagnostic " +
            "observable lives in the Klein cell that is the X-axis flip of the M-active H bilinear " +
            "cell. Why does the framework's diagnostic prescription land here algebraically? The " +
            "structural reason is open.",
            "Trace the algebraic derivation: for each H-class, compute the operator that the " +
            "Marrakesh measurement projects onto, and verify its Klein-cell membership emerges from " +
            "the F87 trichotomy + F1 commutation + dephasing-axis structure.",
            Anchor),
        new OpenQuestion(
            "2:2 truly-kernel split structural meaning",
            "Schicht 3 observation: for truly H + Z-dephasing at N=3, the 4 kernel modes split " +
            "2 in Π²_Z = +1 (Pp + Pm) and 2 in Π²_Z = −1 (Mp + Mm). Naively one might expect all " +
            "kernel mass in Pp (where the truly bilinears live). The 2:2 split reflects the kernel's " +
            "structure: span{I, Z_total, Z_total², Z_total³}, with Z_total^k of even k in Π²_Z=+1 and " +
            "odd k in Π²_Z=−1. Generalisation to k-body kernels is open.",
            "Verify the 2:2 split via direct computation of the kernel basis at N = 3, 4, 5, 6 and " +
            "match the Π²_Z parity to the power-of-Z_total structure.",
            Anchor),
        new OpenQuestion(
            "N ≥ 4 transition: slow non-kernel modes concentrate in Π²_X = −1",
            "Empirical Schicht-3 observation: at N = 2, 3, slow non-kernel modes preserve the " +
            "bilinear apex 1/2 in BOTH Π²_X axes. At N ≥ 4, the slowest 4·(N+1) non-kernel modes " +
            "concentrate in Π²_X = −1 only (Pp + Mp ≈ 0). What's the mechanism? Why does the " +
            "transition happen at N = 4? Possible connection to the 120-enum trichotomy stabilisation " +
            "at N = 4 (per project_v_effect_combinatorial: 15/46/59 N-stable from N = 4 onward).",
            "Trace the slowest-window spectral structure as N grows; identify whether this is a " +
            "Petermann-factor effect (non-normality), a Klein-cell-density effect, or a deeper " +
            "structural transition.",
            Anchor),
        new OpenQuestion(
            "k-body Klein extension (Schicht 1+2 for k ≥ 3 Pauli terms)",
            "F85 lifts F87 trichotomy to arbitrary k-body Pauli terms. The Klein decomposition " +
            "naturally extends (Π²_Z and Π²_X act diagonally on any Pauli string), but the F88 " +
            "table-of-9-bilinears would generalise to a much larger combinatorial structure at " +
            "k ≥ 3. Open: characterise the k-body Klein cells and their scaling.",
            "Enumerate k-body Pauli strings by Klein cell; check whether the bilinear apex 1/2 " +
            "still characterises slow-mode distributions or whether the structure changes.",
            "F85 + F88 + Task #53 (k-body extension)"),
        new OpenQuestion(
            "Half-integer-mirror regime and slow-mode Klein structure",
            "Tom's half-integer family w_XY = N/2 distinguishes odd N (half-integer, no modes on " +
            "mirror axis) from even N (integer, modes on axis). Schicht 3 shows the slow-mode " +
            "Klein structure changes at N = 4, but the small-N (2, 3) preservation does NOT split " +
            "cleanly along odd/even. What's the relationship — if any — between the mirror regime " +
            "and the slow-mode apex transition?",
            "Run Schicht-3 spectral scan at N = 7 (half-integer mirror) once a sparse-Krylov " +
            "method makes it feasible; compare to N = 6 (integer mirror) and N = 4 / 5 transition.",
            Anchor),
    };
}
