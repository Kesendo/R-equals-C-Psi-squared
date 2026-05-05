namespace RCPsiSquared.Core.Resonance;

/// <summary>Bond classification used by the F86 universal-shape statement.
///
/// <para>For a chain of N qubits with N−1 bonds (b = 0..N−2), Endpoint = {0, N−2} and
/// Interior = {1, …, N−3}. The two classes have separate universal HWHM_left/Q_peak
/// values (Interior ≈ 0.756, Endpoint ≈ 0.770) — see PROOF_F86_QPEAK Statement 2.</para>
/// </summary>
public enum BondClass
{
    Endpoint,
    Interior,
}

public static class BondClassExtensions
{
    public static IReadOnlyList<int> BondsOf(this BondClass cls, int numBonds) => cls switch
    {
        BondClass.Endpoint => new[] { 0, numBonds - 1 },
        BondClass.Interior => Enumerable.Range(1, numBonds - 2).ToArray(),
        _ => throw new ArgumentOutOfRangeException(nameof(cls)),
    };

    /// <summary>Classify a bond index. Endpoint = bonds {0, numBonds-1}; Interior otherwise.</summary>
    public static BondClass OfBond(int bond, int numBonds)
    {
        if (numBonds < 1) throw new ArgumentOutOfRangeException(nameof(numBonds));
        if (bond < 0 || bond >= numBonds) throw new ArgumentOutOfRangeException(nameof(bond));
        return (bond == 0 || bond == numBonds - 1) ? BondClass.Endpoint : BondClass.Interior;
    }
}
