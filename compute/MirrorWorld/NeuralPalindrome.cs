namespace MirrorWorld;

// F36/F37: the conditional matrix identity Q J Q + J + 2s I = 0, Q^2 = I.
// Adopted from docs/neural/ALGEBRAIC_PALINDROME_NEURAL.md. Q is an involutive
// permutation; Dale signs alone do not imply the identity. No eigensolver is used.
public static class NeuralPalindrome
{
    // Positive offset s; the spectral centre, if the identity holds, is -s.
    public static double Centre(double tauE, double tauI)
    {
        if (!double.IsFinite(tauE) || tauE <= 0)
            throw new ArgumentOutOfRangeException(nameof(tauE), "Time constant must be finite and positive.");
        if (!double.IsFinite(tauI) || tauI <= 0)
            throw new ArgumentOutOfRangeException(nameof(tauI), "Time constant must be finite and positive.");
        double s = 0.5 / tauE + 0.5 / tauI;
        if (!double.IsFinite(s))
            throw new OverflowException("The centre offset is not representable as a finite double.");
        return s;
    }

    // Malformed permutations throw; a valid permutation of order other than 1 or 2 returns false.
    // The empty permutation is the identity on the zero-dimensional space.
    public static bool IsInvolution(int[] permutation)
    {
        ArgumentNullException.ThrowIfNull(permutation);
        var seen = new bool[permutation.Length];
        foreach (int index in permutation)
        {
            if (index < 0 || index >= permutation.Length || seen[index])
                throw new ArgumentException("Entries must be a permutation of 0 through length minus one.", nameof(permutation));
            seen[index] = true;
        }
        for (int i = 0; i < permutation.Length; i++)
            if (permutation[permutation[i]] != i) return false;
        return true;
    }

    // Maximum absolute entry of Q J Q + J + 2s I. The empty matrix has residual zero.
    // This is a floating-point reading; zero is exact for the constructed dyadic example.
    // Diagonals are centred before summation; unrepresentable intermediates or residuals throw.
    public static double MaxResidual(double[,] j, int[] permutation, double s)
    {
        ArgumentNullException.ThrowIfNull(j);
        ArgumentNullException.ThrowIfNull(permutation);
        int n = j.GetLength(0);
        if (j.GetLength(1) != n)
            throw new ArgumentException("Matrix must be square.", nameof(j));
        if (permutation.Length != n)
            throw new ArgumentException("Permutation length must match the matrix dimension.", nameof(permutation));
        if (!IsInvolution(permutation))
            throw new ArgumentException("Permutation must be an involution.", nameof(permutation));
        if (!double.IsFinite(s))
            throw new ArgumentOutOfRangeException(nameof(s), "Centre offset must be finite.");
        foreach (double entry in j)
            if (!double.IsFinite(entry))
                throw new ArgumentException("Matrix entries must be finite.", nameof(j));

        double worst = 0;
        for (int i = 0; i < n; i++)
            for (int k = 0; k < n; k++)
            {
                double permutedEntry = j[permutation[i], permutation[k]];
                double entry = j[i, k];
                if (i == k)
                {
                    permutedEntry += s;
                    entry += s;
                    if (!double.IsFinite(permutedEntry) || !double.IsFinite(entry))
                        throw new OverflowException("A centred diagonal entry is not representable as a finite double.");
                }
                double residual = permutedEntry + entry;
                if (!double.IsFinite(residual))
                    throw new OverflowException("An entrywise residual is not representable as a finite double.");
                worst = Math.Max(worst, Math.Abs(residual));
            }
        return worst;
    }
}
