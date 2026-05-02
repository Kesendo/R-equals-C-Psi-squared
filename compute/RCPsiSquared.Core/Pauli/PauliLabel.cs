namespace RCPsiSquared.Core.Pauli;

/// <summary>String ↔ <see cref="PauliLetter"/>[] conversion. Labels read left-to-right
/// from site 0 (e.g. "IXZ" means I⊗X⊗Z, site 0 = I, site 1 = X, site 2 = Z).</summary>
public static class PauliLabel
{
    public static PauliLetter[] Parse(string label)
    {
        var letters = new PauliLetter[label.Length];
        for (int i = 0; i < label.Length; i++) letters[i] = PauliLetterExtensions.FromSymbol(label[i]);
        return letters;
    }

    public static string Format(IReadOnlyList<PauliLetter> letters)
    {
        var chars = new char[letters.Count];
        for (int i = 0; i < letters.Count; i++) chars[i] = letters[i].Symbol();
        return new string(chars);
    }
}
