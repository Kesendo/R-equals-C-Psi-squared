using System.Numerics;

namespace RCPsiSquared.Core.Pauli;

/// <summary>A single term in a Pauli-string Hamiltonian: full N-letter sequence plus a coefficient.
///
/// Structural properties (k_body, n_x, n_y, n_z, bit-parities) are derived from the letter
/// sequence and exposed as properties for diagnostic use:
/// <list type="bullet">
///   <item><see cref="Pi2Parity"/>: selects **F87** Π²-class (truly / soft / hard)</item>
///   <item><see cref="YParity"/> — independent at k≥3 (**F85** k-body generalization)</item>
///   <item><see cref="TotalBitA"/> — Z⊗N parity (commutes/anti-commutes with global Z)</item>
/// </list>
/// See docs/ANALYTICAL_FORMULAS.md (F87, F81, F85) and the bit_a/bit_b convention in
/// <see cref="PauliLetter"/>.
/// </summary>
public sealed record PauliTerm(IReadOnlyList<PauliLetter> Letters, Complex Coefficient)
{
    public int N => Letters.Count;

    public int KBody
    {
        get
        {
            int k = 0;
            foreach (var L in Letters) if (L != PauliLetter.I) k++;
            return k;
        }
    }

    public int Nx => Count(PauliLetter.X);
    public int Ny => Count(PauliLetter.Y);
    public int Nz => Count(PauliLetter.Z);

    /// <summary>Total bit_a (X+Y count). Equal to <see cref="Nx"/> + <see cref="Ny"/>.</summary>
    public int TotalBitA => PauliIndex.TotalBitA(Letters);

    /// <summary>Π²-parity (Y+Z count mod 2). 0 = Π²-even, 1 = Π²-odd.</summary>
    public int Pi2Parity => PauliIndex.TotalBitBParity(Letters);

    /// <summary>Klein-Vierergruppe Z₂×Z₂ index (bit_a parity, bit_b parity) of the term.
    /// bit_a = (#X + #Y) mod 2; bit_b = (#Y + #Z) mod 2 = <see cref="Pi2Parity"/>.
    /// 2-body bilinears: XX/YY/ZZ → (0,0) Mother (truly), XY/YX → (0,1) (Π²-odd subset),
    /// YZ/ZY → (1,0) (Π²-even non-truly), XZ/ZX → (1,1) (Π²-odd subset).</summary>
    public (int BitA, int BitB) KleinIndex => (TotalBitA & 1, Pi2Parity);

    /// <summary>Full Z₂³ structural signature (bit_a, bit_b, Y-parity). At k=2 Y-parity
    /// is determined by Klein; at k≥3 it is independent and the polarity is Z₂³ (8 sectors).</summary>
    public (int BitA, int BitB, int YParity) FullZ2Signature => (TotalBitA & 1, Pi2Parity, YParity);

    /// <summary>Y-count modulo 2. Independent from Klein at k ≥ 3-body.</summary>
    public int YParity => Ny & 1;

    public string Label => PauliLabel.Format(Letters);

    private int Count(PauliLetter letter)
    {
        int c = 0;
        foreach (var L in Letters) if (L == letter) c++;
        return c;
    }

    public static PauliTerm SingleSite(int N, int site, PauliLetter letter, Complex coefficient)
    {
        var letters = new PauliLetter[N];
        for (int i = 0; i < N; i++) letters[i] = PauliLetter.I;
        letters[site] = letter;
        return new PauliTerm(letters, coefficient);
    }

    public static PauliTerm TwoSite(int N, int site1, PauliLetter letter1, int site2, PauliLetter letter2, Complex coefficient)
    {
        var letters = new PauliLetter[N];
        for (int i = 0; i < N; i++) letters[i] = PauliLetter.I;
        letters[site1] = letter1;
        letters[site2] = letter2;
        return new PauliTerm(letters, coefficient);
    }
}
