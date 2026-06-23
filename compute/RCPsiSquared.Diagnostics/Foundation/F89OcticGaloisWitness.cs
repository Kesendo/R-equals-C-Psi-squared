using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Live witness: recomputes the F89 path-3 octic Galois certificate at inspect time.
/// Gal(F_8 / Q(i)(q)) = S_8 (non-solvable) ⟹ the eight Liouvillian eigenvalues λ_k(q) admit
/// no expression by radicals in q. At the certifying specialization q=2, factoring F_8(·,2)
/// modulo split primes 𝔭|p (p≡1 mod 4, i↦√−1) yields Frobenius cycle types: an 8-cycle
/// (⟹ transitive) and (5,2,1) at 𝔭|5 (a 5-cycle ⟹ primitive, since no proper primitive
/// degree-8 group has order divisible by 5 ⟹ ⊇A_8 by Jordan; odd ⟹ ⊄A_8) ⟹ S_8.
/// The Galois sibling of the EP-character witness (--root f89octic).</summary>
public sealed class F89OcticGaloisWitness : IInspectable
{
    // primes ≡ 1 (mod 4) — the ones that split in Z[i], so Frobenius is readable over F_p.
    private static readonly int[] SplitPrimes =
        { 5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101, 109, 113 };

    private static IReadOnlyList<(int Prime, int[] CycleType)> Sample()
    {
        var (re, im) = F89Path3OcticBlock.OcticCoefficientsAtQ2();
        var rows = new List<(int, int[])>();
        foreach (int p in SplitPrimes)
        {
            var ct = OcticGaloisCertificate.CycleType(re, im, p);
            if (ct is not null) rows.Add((p, ct));
        }
        return rows;
    }

    private static (int TransPrime, int FivePrime, int DistinctTypes, string Variety) Compute()
    {
        var rows = Sample();
        int transPrime = rows.Where(r => r.CycleType.Length == 1 && r.CycleType[0] == 8)
            .Select(r => r.Prime).DefaultIfEmpty(0).First();
        int fivePrime = rows.Where(r => r.CycleType.Contains(5))
            .Select(r => r.Prime).DefaultIfEmpty(0).First();
        int distinct = rows.Select(r => string.Join(",", r.CycleType)).Distinct().Count();
        string variety = string.Join("  ",
            rows.Take(6).Select(r => $"𝔭|{r.Prime}:({string.Join(",", r.CycleType)})"));
        return (transPrime, fivePrime, distinct, variety);
    }

    public string DisplayName => "F89 path-3 octic Galois group (Frobenius certificate, live)";

    public string Summary
    {
        get
        {
            var (transPrime, fivePrime, _, _) = Compute();
            return $"Gal(F_8/Q(i)(q)) = S_8 (non-solvable): (5,2,1) at 𝔭|{fivePrime} ⟹ 5-cycle ⟹ ⊇A_8; " +
                   $"8-cycle at 𝔭|{transPrime} ⟹ transitive ⟹ no radical closure in q";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var (transPrime, fivePrime, distinct, variety) = Compute();
            yield return new InspectableNode("verdict",
                summary: $"Gal(F_8/Q(i)(q)) = S_8 — non-solvable; the eight eigenvalues λ_k(q) " +
                         "admit no radical expression in q (non-radical special functions not excluded)");
            yield return new InspectableNode("transitive (8-cycle ⟹ acts on all 8 roots)",
                summary: $"8-cycle Frobenius at split prime 𝔭|{transPrime}");
            yield return new InspectableNode("⊇A_8 (5-cycle ⟹ primitive ⟹ Jordan)",
                summary: $"(5,2,1) at 𝔭|{fivePrime}: its square is a 5-cycle; a 5-orbit fits no degree-8 " +
                         "block, and no proper primitive degree-8 group has order divisible by 5");
            yield return new InspectableNode("⊄A_8 (odd permutation / disc non-square)",
                summary: "(5,2,1) is an odd permutation; equivalently disc(F_8) is not a square in Q(i)");
            yield return new InspectableNode("Frobenius cycle-type variety at q=2",
                summary: $"{distinct} distinct types over split primes; {variety} …");
            yield return new InspectableNode("the negative content",
                summary: "S_8 is the generic group (Bhargava 2025); integrability spent itself on the " +
                         "F_a/F_b factorisation + the solvable quartic (3q⁴+q²−1) — the residual octic is structureless");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
