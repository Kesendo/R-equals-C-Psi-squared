using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Open theoretical items for the Π² Klein layer.</summary>
public static class Pi2OpenQuestions
{
    private const string Anchor = "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBase + Schicht-3 spectral observations 2026-05-03";

    public static IReadOnlyList<OpenQuestion> Standard { get; } = new[]
    {
        new OpenQuestion(
            "N ≥ 4 transition: slow non-kernel modes concentrate in Π²_X = −1",
            "Empirical Schicht-3 observation: at N = 2, 3, slow non-kernel modes preserve the " +
            "bilinear apex 1/2 in BOTH Π²_X axes. At N ≥ 4, the slowest 4·(N+1) non-kernel modes " +
            "concentrate in Π²_X = −1 only (Pp + Mp ≈ 0). What's the mechanism? Why does the " +
            "transition happen at N = 4? Possible connection to the 120-enum trichotomy stabilisation " +
            "at N = 4 (per project_v_effect_combinatorial: 15/46/59 N-stable, already at N = 3 and through N = 5).",
            "Trace the slowest-window spectral structure as N grows; identify whether this is a " +
            "Petermann-factor effect (non-normality), a Klein-cell-density effect, or a deeper " +
            "structural transition.",
            Anchor),
        new OpenQuestion(
            "k-body Klein extension (Schicht 1+2 for k ≥ 3 Pauli terms)",
            "F85 lifts F87 trichotomy to arbitrary k-body Pauli terms. The Klein decomposition " +
            "naturally extends (Π²_Z and Π²_X act diagonally on any Pauli string). The cells are " +
            "counted by the letter characters: each letter carries (bit_a, bit_b), I (0,0), X (1,0), " +
            "Y (1,1), Z (0,1), and every nontrivial character of the Klein group sums to 0 over " +
            "{I, X, Y, Z} and to −1 over {X, Y, Z}, so with identities each cell holds 4^(k−1) " +
            "strings (the all-identity one in (0,0)), and without identity letters the (0,0) cell holds (3^k + 3(−1)^k)/4 and each other cell " +
            "(3^k − (−1)^k)/4 (k = 2: 3 + 2 + 2 + 2, the F88a table of nine bilinears; " +
            "KBodyKleinCellCountTests). Open: whether the bilinear apex 1/2 still characterises the " +
            "slow-mode distributions at k ≥ 3.",
            "Read the slow modes per Klein cell at k ≥ 3 and compare with the bilinear apex 1/2.",
            "F85 + F88a + Task #53 (k-body extension)"),
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
