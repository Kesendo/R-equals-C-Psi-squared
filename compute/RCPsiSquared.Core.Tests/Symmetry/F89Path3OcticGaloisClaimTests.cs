using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F89Path3OcticGaloisClaimTests
{
    private static F89Path3OcticGaloisClaim BuildClaim()
    {
        var f89 = new F89TopologyOrbitClosure(new Pi2DyadicLadderClaim());
        var atLock = new F89PathKAtLockMechanismClaim(f89);
        return new F89Path3OcticGaloisClaim(f89, atLock);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void DiscriminantPolynomialDegreeInQ_Is52()
    {
        Assert.Equal(52, F89Path3OcticGaloisClaim.DiscriminantPolynomialDegreeInQ);
    }

    [Fact]
    public void DiscriminantFactorisationStructure_IsConstantQ24SquareTimesP10()
    {
        // disc(F_8) = constant prefactor · q^24 · (3q⁴+q²-1)² · P_10(q²)
        // Tuple encodes degree contributions only (constant prefactor omitted; always present once).
        var (qPower, squareFactorDegInQ, p10DegInQ) =
            F89Path3OcticGaloisClaim.DiscriminantFactorisationStructure;
        Assert.Equal(24, qPower);
        Assert.Equal(8, squareFactorDegInQ);   // (3q⁴+q²-1)² has degree 8 in q
        Assert.Equal(20, p10DegInQ);            // P_10(q²) has degree 20 in q (= 10 in q²)
        // Total degree check: 24 (q^24) + 8 ((3q⁴+q²-1)²) + 20 (P_10) = 52
        Assert.Equal(52, qPower + squareFactorDegInQ + p10DegInQ);
    }

    [Fact]
    public void GalNotInA8_IsTier1True()
    {
        // disc non-square → Gal(F_8) ⊄ A_8 (verified at q ∈ {½, 1, 3/2, 2, 3};
        // decisively, an odd Frobenius element at q0=2).
        Assert.True(F89Path3OcticGaloisClaim.GalNotInA8);
    }

    [Fact]
    public void GalIsS8_IsTier1True()
    {
        // Gal(F_8 / Q(i)(q)) = S_8, established by specialization + Dedekind + Jordan:
        // a (5,2,1) Frobenius at 𝔭|5, q0=2, gives a 5-cycle (⇒ primitive ⇒ ⊇A_8 by
        // Jordan) and an odd permutation (⇒ ⊄A_8) ⇒ S_8. Robust to base Q(i,√5).
        Assert.True(F89Path3OcticGaloisClaim.GalIsS8);
    }

    [Fact]
    public void NonSolvableConjecture_IsNowClosed()
    {
        // Superseded: non-solvability is no longer an open conjecture — S_8 is
        // certified, so the octic is non-solvable (S_8 ⊇ A_8 simple for n ≥ 5).
        Assert.False(F89Path3OcticGaloisClaim.NonSolvableConjecture_IsOpen);
    }

    [Fact]
    public void CertifyingPrimeAtQ0Eq2_Is5()
    {
        // The split prime 𝔭|5 (Z[i]/𝔭 = F_5, i ↦ 2) factors F_8(·,2) to cycle type (5,2,1).
        Assert.Equal(5, F89Path3OcticGaloisClaim.CertifyingPrimeAtQ0Eq2);
    }
}
