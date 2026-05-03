using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>The Marrakesh f83 4-class discrimination test (April 30, 2026, ibm_marrakesh
/// job d7pol1e7g7gs73cf7j90, path [4,5,6]) chose one fingerprint observable per F87 H-class.
/// Read through the framework's own 2-axis Klein lens (Π²_Z, Π²_X), each (H, observable)
/// pair sits in a specific structural relationship: the diagnostic observable lives in the
/// **X-axis-flipped Klein-cell** of the M-active bilinears of the Hamiltonian.
///
/// <para>For Mixed H = XX+XY the M-active bilinear is XY (truly XX drops by Master Lemma);
/// the rule applies to the M-active Klein-cell, not the full Hamiltonian.</para>
///
/// <para>This is an empirical structural pattern from the framework's own perspective,
/// not a derivation. The pattern is documented in ConfirmationsRegistry entry
/// "f83_pi2_class_signature_marrakesh" with measured expectation values; here we lock the
/// Klein-cell algebra alone.</para>
/// </summary>
public class Pi2KleinHardwareViewTests
{
    private static (int z, int x) KleinCell(params PauliLetter[] letters) =>
        (PiOperator.SquaredEigenvalue(letters, PauliLetter.Z),
         PiOperator.SquaredEigenvalue(letters, PauliLetter.X));

    [Fact]
    public void Marrakesh_TrulyFingerprint_LivesInXFlippedCellOfTrulyBilinears()
    {
        // Truly H: XX, YY, ZZ — bilinears in cell Pp = (+, +).
        // f83 fingerprint observable: ⟨Y₀ I Z₂⟩.
        var trulyBilinear = KleinCell(PauliLetter.X, PauliLetter.X);
        var observable = KleinCell(PauliLetter.Y, PauliLetter.I, PauliLetter.Z);

        Assert.Equal((+1, +1), trulyBilinear);    // Pp
        Assert.Equal((+1, -1), observable);       // Pm — X-axis flip of Pp
    }

    [Fact]
    public void Marrakesh_Pi2EvenNonTrulyFingerprint_LivesInXFlippedCellOfYZBilinears()
    {
        // Pi2EvenNonTruly H: YZ, ZY — bilinears in cell Pm = (+, −).
        // f83 fingerprint observable: ⟨X₀ I X₂⟩.
        var nonTrulyBilinear = KleinCell(PauliLetter.Y, PauliLetter.Z);
        var observable = KleinCell(PauliLetter.X, PauliLetter.I, PauliLetter.X);

        Assert.Equal((+1, -1), nonTrulyBilinear); // Pm
        Assert.Equal((+1, +1), observable);       // Pp — X-axis flip of Pm
    }

    [Fact]
    public void Marrakesh_Pi2OddPureFingerprint_LivesInXFlippedCellOfXYBilinears()
    {
        // Pi2OddPure (subgroup A) H: XY, YX — bilinears in cell Mp = (−, +).
        // f83 fingerprint observable: ⟨X₀ I Z₂⟩.
        var oddPureBilinear = KleinCell(PauliLetter.X, PauliLetter.Y);
        var observable = KleinCell(PauliLetter.X, PauliLetter.I, PauliLetter.Z);

        Assert.Equal((-1, +1), oddPureBilinear);  // Mp
        Assert.Equal((-1, -1), observable);       // Mm — X-axis flip of Mp
    }

    [Fact]
    public void Marrakesh_MixedFingerprint_LivesInXFlippedCellOfMActiveBilinear()
    {
        // Mixed H: XX (truly) + XY (Π²-odd). M-active bilinear is XY in Mp = (−, +)
        // (truly XX drops by Master Lemma); fingerprint observable: ⟨Z₀ I X₂⟩.
        var mActiveBilinear = KleinCell(PauliLetter.X, PauliLetter.Y);
        var observable = KleinCell(PauliLetter.Z, PauliLetter.I, PauliLetter.X);

        Assert.Equal((-1, +1), mActiveBilinear);  // Mp (M-active part of Mixed H)
        Assert.Equal((-1, -1), observable);       // Mm — X-axis flip of Mp
    }

    [Fact]
    public void Marrakesh_Soft_BreakAnchor_Observable_IsPi2Odd()
    {
        // The 2026-04-26 palindrome_trichotomy soft-break confirmation measured ⟨X₀ Z₂⟩
        // (no I in the middle — direct 2-letter spec at sites 0, 2). Klein cell:
        var softBreakObservable = KleinCell(PauliLetter.X, PauliLetter.I, PauliLetter.Z);
        // For N=3 the explicit-I form is what was measured (ConfirmationsRegistry observable
        // string "<X_0 Z_2>" implicitly identity at site 1).
        Assert.Equal((-1, -1), softBreakObservable);  // Mm — same cell as f83 mixed observable
    }

    [Fact]
    public void Klein_View_Of_All_Marrakesh_Observables_Forms_OneXFlipPattern()
    {
        // Locking the structural pattern: each (M-active-H-cell, observable-cell) pair from
        // the f83 fingerprint test is related by an X-axis flip. The pattern itself is the
        // observation; this test enforces it as a single invariant across the four classes.
        var pairs = new[]
        {
            ("truly", (+1, +1), (+1, -1)),                      // Pp → Pm
            ("pi2_even_nontruly", (+1, -1), (+1, +1)),          // Pm → Pp
            ("pi2_odd_pure", (-1, +1), (-1, -1)),               // Mp → Mm
            ("mixed (M-active part)", (-1, +1), (-1, -1)),       // Mp → Mm
        };

        foreach (var (label, hCell, obsCell) in pairs)
        {
            // X-axis flip: Π²_Z stays, Π²_X negates.
            var expected = (hCell.Item1, -hCell.Item2);
            Assert.True(expected == obsCell,
                $"{label}: expected X-flip of H-cell {hCell} = {expected}, observable cell {obsCell}");
        }
    }
}
