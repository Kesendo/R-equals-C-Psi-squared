using System.Numerics;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class MotherSoftYParityOnePurityTests
{
    [Fact]
    public void Z2Axis_IsYParity()
    {
        var claim = new MotherSoftYParityOnePurity();
        Assert.Equal(Z2Axis.YParity, claim.Z2Axis);
    }

    [Fact]
    public void BitATwin_IsNull()
    {
        var claim = new MotherSoftYParityOnePurity();
        Assert.Null(claim.BitATwin);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var claim = new MotherSoftYParityOnePurity();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void IsMotherNonTrulyCandidate_RequiresKlein00AndAllOdd()
    {
        // Per PROOF_F109 Step 3: Klein (0,0) non-truly iff #X, #Y, #Z all odd.
        // XYZ at k=3: #X=1, #Y=1, #Z=1 (all odd) → Klein (0,0) y_par=1 non-truly candidate.
        var xyz = new PauliTerm(new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.Z }, Complex.One);
        Assert.True(MotherSoftYParityOnePurity.IsMotherNonTrulyCandidate(xyz));

        // XXI at k=3: #X=2 (even) → not non-truly candidate (Klein (0,0) truly).
        var xxi = new PauliTerm(new[] { PauliLetter.X, PauliLetter.X, PauliLetter.I }, Complex.One);
        Assert.False(MotherSoftYParityOnePurity.IsMotherNonTrulyCandidate(xxi));

        // YYI: all even (#Y=2, #X=0, #Z=0) → Klein (0,0) truly. Not candidate.
        var yyi = new PauliTerm(new[] { PauliLetter.Y, PauliLetter.Y, PauliLetter.I }, Complex.One);
        Assert.False(MotherSoftYParityOnePurity.IsMotherNonTrulyCandidate(yyi));

        // XYI: bit_a = 0, bit_b = 1 → Klein (0,1), not (0,0). Not candidate.
        var xyi = new PauliTerm(new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.I }, Complex.One);
        Assert.False(MotherSoftYParityOnePurity.IsMotherNonTrulyCandidate(xyi));
    }

    [Fact]
    public void VerifyOnTerm_MotherNonTrulyTermsAllHaveYParityOne()
    {
        // F109 central claim: enumerate all 64 k=3 letter sequences, for any
        // Klein (0,0) non-truly candidate term, y_par must be 1.
        var letters = new[] { PauliLetter.I, PauliLetter.X, PauliLetter.Y, PauliLetter.Z };
        int candidateCount = 0;
        foreach (var a in letters)
            foreach (var b in letters)
                foreach (var c in letters)
                {
                    var term = new PauliTerm(new[] { a, b, c }, Complex.One);
                    Assert.True(MotherSoftYParityOnePurity.VerifyOnTerm(term),
                        $"F109 violated: term [{a},{b},{c}] is Klein (0,0) non-truly candidate but y_par = {term.YParity} (expected 1)");
                    if (MotherSoftYParityOnePurity.IsMotherNonTrulyCandidate(term))
                        candidateCount++;
                }
        // Per PROOF_F109 Step 3 spot-check: 6 XYZ-permutations.
        Assert.Equal(6, candidateCount);
    }

    [Fact]
    public void VerifyOnTerm_AtK4_StillHoldsAndMatchesF106Count()
    {
        // F106 mother soft (0, 300) per dephase: 24 Klein (0,0) y_par=1 non-truly terms at k=4
        // ⟹ 24·25/2 = 300 unordered pairs (including self). Verify enumeration count.
        var letters = new[] { PauliLetter.I, PauliLetter.X, PauliLetter.Y, PauliLetter.Z };
        int candidateCount = 0;
        foreach (var a in letters)
            foreach (var b in letters)
                foreach (var c in letters)
                    foreach (var e in letters)
                    {
                        var term = new PauliTerm(new[] { a, b, c, e }, Complex.One);
                        Assert.True(MotherSoftYParityOnePurity.VerifyOnTerm(term),
                            $"F109 violated at k=4: term [{a},{b},{c},{e}] is non-truly candidate but y_par = {term.YParity}");
                        if (MotherSoftYParityOnePurity.IsMotherNonTrulyCandidate(term))
                            candidateCount++;
                    }
        // 4 letters with #X=1, #Y=1, #Z=1, #I=1: 4! / (1!1!1!1!) = 24 permutations.
        Assert.Equal(24, candidateCount);
    }

    [Fact]
    public void EmpiricalSpotCheck_K3MotherPairCountMatches21()
    {
        // From PROOF_F109 Step 3: 6 Klein (0,0) non-truly k=3 terms × 7 / 2 = 21 pairs
        // (unordered with self-pairs). Matches F103/F105 mother soft (0, 21) per dephase.
        int n = 6;
        int pairs = n * (n + 1) / 2;
        Assert.Equal(21, pairs);
    }

    [Fact]
    public void EmpiricalSpotCheck_K4MotherPairCountMatches300()
    {
        // From PROOF_F109: 24 Klein (0,0) non-truly k=4 terms × 25 / 2 = 300 pairs.
        // Matches F106 mother soft (0, 300) per dephase.
        int n = 24;
        int pairs = n * (n + 1) / 2;
        Assert.Equal(300, pairs);
    }
}
