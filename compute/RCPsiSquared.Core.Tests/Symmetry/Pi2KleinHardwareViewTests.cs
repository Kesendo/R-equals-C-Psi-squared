using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>The Marrakesh f83 4-class discrimination test (April 30, 2026, ibm_marrakesh
/// job d7pol1e7g7gs73cf7j90, path [4,5,6]) chose one fingerprint observable per F87 H-class.
/// Read through the framework's own 2-axis Klein lens (Π²_Z, Π²_X), each (H, observable)
/// pair sits in a specific structural relationship: the diagnostic observable lives in the
/// **X-axis-flipped Klein-cell** of the M-active bilinears of the Hamiltonian.
///
/// <para>The f83 mixed H is XY + YZ, its M-active cells Mp (XY) and Pm (YZ); its fingerprint
/// Z₀IX₂ sits in Mm, the X-axis flip of Mp. Truly bilinears sit in Pp and shift no cell; the start's
/// odd-X strings supply Pm, and the y-parity turn lets the one-Y observable Y₀IZ₂ there be nonzero.</para>
///
/// <para>The reason, in docs/ANALYTICAL_FORMULAS.md's Klein view: Z dephasing keeps every string,
/// a commutator adds the term's Klein cell and turns the y-parity by y(h) + 1, and the f83 start
/// |+-+> is real with X strings in Pm (simulations/cube_old_questions_gate.py G5). The measured
/// values are in ConfirmationsRegistry entry "f83_pi2_class_signature_marrakesh"; here the
/// Klein-cell algebra alone is locked.</para>
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
        // Mixed H: XY + YZ, M-active cells Mp (XY) and Pm (YZ); fingerprint observable ⟨Z₀ I X₂⟩,
        // locked here against the XY cell.
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
        // Each (M-active-H-cell, observable-cell) pair from the f83 fingerprint test is related by an
        // X-axis flip. The reason: Z dephasing keeps strings, commutators add Klein cells and turn the
        // y-parity, and the f83 start |+-+> is real with X strings in Pm (cube_old_questions_gate.py G5).
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
