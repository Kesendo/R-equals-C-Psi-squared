namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>F71-orbit-refined bond sub-classification for F86 c=2 K_b observable
/// (Tier 1, 2026-05-13). Extends the binary BondClass {Endpoint, Interior} with
/// the empirically-pinned sub-classes that correspond to distinct HWHM_left/Q_peak
/// asymptotic behaviors per the F89 AT-locked + g_eff(N, b) closed form.</summary>
public enum BondSubClass
{
    Endpoint,
    Flanking,
    Mid,
    CentralSelfPaired,
    Orbit2Escape,
    CentralEscapeOrbit3,
}

public static class BondSubClassExtensions
{
    /// <summary>Classify bond by F71 orbit + Q_peak escape threshold per
    /// PROOF_F90_F86C2_BRIDGE.md numerical verification table (N=5..8).</summary>
    public static BondSubClass Classify(int n, int bondIndex, double qPeak)
    {
        int nBonds = n - 1;
        if (bondIndex == 0 || bondIndex == nBonds - 1)
            return BondSubClass.Endpoint;
        if (qPeak > 12.0)
            return BondSubClass.CentralEscapeOrbit3;
        if (qPeak > 4.0)
            return BondSubClass.Orbit2Escape;
        if (n % 2 == 0 && bondIndex == nBonds / 2)
            return BondSubClass.CentralSelfPaired;
        if (bondIndex == 1 || bondIndex == nBonds - 2)
            return BondSubClass.Flanking;
        return BondSubClass.Mid;
    }
}
