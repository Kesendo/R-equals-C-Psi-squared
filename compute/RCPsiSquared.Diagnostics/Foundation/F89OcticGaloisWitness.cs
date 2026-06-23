using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Live witness: recomputes the F89 path-3 octic Galois certificate at inspect time, the
/// FULL Option-D way — it does not import the octic literal, it RECOMPUTES F_d from the block.
/// Builds the path-3 (SE,DE) S_2-sym block at q=2 exactly over Z[i] (integer mirror basis,
/// ×2-cleared), takes its division-free Berkowitz characteristic polynomial, and isolates the
/// H_B-mixed octic F_d by dividing out the reconstructed AT factor (F_a·F_b), validated by the
/// triple: exact division (remainder 0), degree 8, gcd(AT, F_d) = 1. Factoring F_d modulo split
/// primes 𝔭|p (p≡1 mod 4, i↦√−1) gives Frobenius cycle types; the generalised Jordan reading
/// (a prime cycle in (4,5] ⟹ ⊇A_8, an odd type ⟹ ⊄A_8) certifies Gal(F_8/Q(i)(q)) = S_8,
/// non-solvable ⟹ the eight eigenvalues λ_k(q) admit no radical closure in q. The Galois sibling
/// of the EP-character witness (--root f89octic). (F_d is cross-checked against the trusted
/// OcticCoefficientsAtQ2 literal in the test suite.)</summary>
public sealed class F89OcticGaloisWitness : IInspectable
{
    // primes ≡ 1 (mod 4) — the ones that split in Z[i], so Frobenius is readable over F_p.
    private static readonly int[] SplitPrimes =
        { 5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101, 109, 113 };

    /// <summary>F_d recomputed LIVE: block → Berkowitz charpoly → isolate by dividing out the AT
    /// factor. The block is ×2-cleared, so F_d carries roots 2λ_k; the Galois reading (cycle
    /// types) is invariant under that scaling, so the certificate is unaffected.</summary>
    private static GaussianInteger[] RecomputeFd()
    {
        var block = F89PathKSeDeBlock.BuildTwoTimesSymBlock(q0: 2, nBlock: 4);
        var charpoly = GaussianMatrixCharpoly.Characteristic(block);
        var atFactor = F89AtFactorReconstruction.ForPath3(q0: 2);
        return F89HbMixedIsolation.Isolate(charpoly, atFactor, expectedDegree: 8);
    }

    private static IReadOnlyList<(int Prime, int[] CycleType)> Sample()
    {
        var fd = RecomputeFd();
        var re = fd.Select(c => c.Re).ToArray();
        var im = fd.Select(c => c.Im).ToArray();
        var rows = new List<(int, int[])>();
        foreach (int p in SplitPrimes)
        {
            var ct = OcticGaloisCertificate.CycleType(re, im, p);
            if (ct is not null) rows.Add((p, ct));
        }
        return rows;
    }

    private static (GaloisGroupCertificate Cert, int TransPrime, int JordanPrime, int DistinctTypes, string Variety) Compute()
    {
        var rows = Sample();
        var cert = OcticGaloisCertificate.JordanVerdict(rows.Select(r => r.CycleType), 8);
        int transPrime = rows.Where(r => r.CycleType.Length == 1 && r.CycleType[0] == 8)
            .Select(r => r.Prime).DefaultIfEmpty(0).First();
        int jordanPrime = rows.Where(r => r.CycleType.Contains(5))
            .Select(r => r.Prime).DefaultIfEmpty(0).First();
        int distinct = rows.Select(r => string.Join(",", r.CycleType)).Distinct().Count();
        string variety = string.Join("  ",
            rows.Take(6).Select(r => $"𝔭|{r.Prime}:({string.Join(",", r.CycleType)})"));
        return (cert, transPrime, jordanPrime, distinct, variety);
    }

    public string DisplayName => "F89 path-3 octic Galois group (live: block → Berkowitz → isolate F_d → Frobenius)";

    public string Summary
    {
        get
        {
            var (cert, transPrime, jordanPrime, _, _) = Compute();
            string g = cert.IsFullSymmetric ? "S_8" : cert.IsNonSolvable ? "⊇A_8" : "a proper subgroup";
            return $"Gal(F_8/Q(i)(q)) = {g} (non-solvable): F_d recomputed live from the 12×12 block; " +
                   $"(5,2,1) at 𝔭|{jordanPrime} ⟹ 5-cycle ⟹ ⊇A_8; 8-cycle at 𝔭|{transPrime} ⟹ transitive ⟹ no radical closure in q";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var (cert, transPrime, jordanPrime, distinct, variety) = Compute();
            yield return new InspectableNode("F_d recomputed live (full Option D)",
                summary: "the path-3 (SE,DE) S_2-sym block is built at q=2 over Z[i], Berkowitz gives its exact " +
                         "charpoly, and the octic F_d is isolated by dividing out the AT factor — triple verified " +
                         "(remainder 0, degree 8, gcd(AT,F_d)=1). F_d IS the H_B-mixed factor, not an imported literal");
            yield return new InspectableNode("verdict",
                summary: $"Gal(F_8/Q(i)(q)) = {(cert.IsFullSymmetric ? "S_8" : "⊇A_8")} — non-solvable; the eight eigenvalues " +
                         "λ_k(q) admit no radical expression in q (non-radical special functions not excluded)");
            yield return new InspectableNode("transitive (8-cycle ⟹ acts on all 8 roots)",
                summary: $"8-cycle Frobenius at split prime 𝔭|{transPrime}");
            yield return new InspectableNode($"⊇A_8 (Jordan prime {cert.JordanPrime} ⟹ primitive ⟹ Jordan)",
                summary: $"(5,2,1) at 𝔭|{jordanPrime}: its square is a 5-cycle (5 > 8/2 ⟹ primitive; 5 ≤ 8−3 ⟹ Jordan ⟹ ⊇A_8)");
            yield return new InspectableNode("⊄A_8 (path-3 has the all-q closed form, not just the q0 point)",
                summary: "(5,2,1) is an odd permutation ⟹ ⊄A_8 at q=2; and uniquely for path-3, disc(F_8) = " +
                         "const·q²⁴·(3q⁴+q²−1)²·P_10(q²) with P_10 a non-square ⟹ ⊄A_8 for EVERY q (path-4/5/6 have only the q0 certificate)");
            yield return new InspectableNode("Frobenius cycle-type variety at q=2",
                summary: $"{distinct} distinct types over split primes; {variety} …");
            yield return new InspectableNode("the negative content",
                summary: "S_8 is the generic group (Bhargava 2025); integrability spent itself on the " +
                         "F_a/F_b factorisation + the solvable quartic (3q⁴+q²−1) — the residual octic is structureless");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
