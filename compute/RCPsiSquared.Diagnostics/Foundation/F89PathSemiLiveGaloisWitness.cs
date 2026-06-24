using System.Collections.Generic;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Semi-live witness for a path-k H_B-mixed factor (k = 4,5,6). Builds the (SE,DE) S_2-sym
/// block LIVE at q=2 over Z[i], takes its Berkowitz characteristic polynomial, and VERIFIES the
/// committed oracle F_d divides it with a coprime degree-AT complement (the triple) — proving F_d
/// IS the live block's H_B-mixed factor — then reads F_d's Frobenius cycle types: Gal(F_d) = S_d,
/// non-solvable. F_d's coefficients are imported (oracle) but live-verified against the block; the
/// path-3 octic is the fully-live sibling (F89OcticGaloisWitness). The heavy recompute lives in the
/// Children getter (so listing path-6 stays cheap until it is expanded).</summary>
public sealed class F89PathSemiLiveGaloisWitness : IInspectable
{
    private readonly int _k;

    public F89PathSemiLiveGaloisWitness(int k) => _k = k;

    private int D => F89PathKFdOracle.Degree(_k);

    public string DisplayName => $"F89 path-{_k} H_B-mixed factor F_{D} Galois group (semi-live)";

    public string Summary =>
        $"Gal(F_{D}/Q(i)(q)) = S_{D} (non-solvable) — expand to recompute: F_d verified against the live block, then Frobenius ⟹ S_{D}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            // 1) live block -> Berkowitz charpoly -> verify the oracle F_d is its factor (the triple).
            var block = F89PathKSeDeBlock.BuildTwoTimesSymBlock(q0: 2, nBlock: _k + 1);
            var charpoly = GaussianMatrixCharpoly.Characteristic(block);
            F89HbMixedIsolation.Isolate(charpoly, F89PathKFdOracle.FdScaled(_k), F89PathKFdOracle.AtDegree(_k));
            yield return new InspectableNode("F_d verified against the live block (semi-live D)",
                summary: $"the path-{_k} (SE,DE) S_2-sym block (dim {F89PathKFdOracle.SymDim(_k)}) is built at q=2 over Z[i]; " +
                         $"Berkowitz gives its charpoly; the committed F_{D} divides it exactly with a coprime degree-" +
                         $"{F89PathKFdOracle.AtDegree(_k)} AT complement (remainder 0, gcd = 1) ⟹ F_d IS the live block's H_B-mixed factor");

            // 2) Frobenius cycle types of F_d -> the generalised Jordan verdict.
            var (re, im) = F89PathKFdOracle.Fd(_k);
            var types = new List<int[]>();
            int transPrime = 0;
            GaloisGroupCertificate cert = default;
            foreach (int p in SplitPrimesUpTo(4000))
            {
                var ct = OcticGaloisCertificate.CycleType(re, im, p);
                if (ct is null) continue;
                types.Add(ct);
                if (transPrime == 0 && ct.Length == 1 && ct[0] == D) transPrime = p;
                cert = OcticGaloisCertificate.JordanVerdict(types, D);
                if (cert.IsFullSymmetric && transPrime != 0) break;
            }

            string g = cert.IsFullSymmetric ? $"S_{D}" : cert.IsNonSolvable ? $"⊇A_{D}" : "a proper subgroup";
            yield return new InspectableNode("verdict",
                summary: $"Gal(F_{D}/Q(i)(q)) = {g} — non-solvable; the {D} decay rates λ_k(q) admit no radical expression in q = J/γ");
            yield return new InspectableNode($"transitive (degree-{D} irreducible reduction)",
                summary: $"a {D}-cycle Frobenius at split prime 𝔭|{transPrime}");
            yield return new InspectableNode($"⊇A_{D} (Jordan cycle length {cert.JordanPrime})",
                summary: $"a {cert.JordanPrime}-cycle ({cert.JordanPrime} > {D}/2 ⟹ primitive; {cert.JordanPrime} ≤ {D}−3 ⟹ Jordan) ⟹ ⊇A_{D}");
            yield return new InspectableNode($"⊄A_{D} (odd Frobenius — q0 certificate only)",
                summary: $"an odd cycle type appears ⟹ ⊄A_{D}; certified at q0=2 + specialization-only-shrinks (path-k has no all-q discriminant, unlike path-3)");
            yield return new InspectableNode("the negative content",
                summary: $"S_{D} is the generic group; integrability spends itself on the AT-locked F_a/F_b factorisation, leaving the H_B-mixed F_{D} structureless");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    private static IEnumerable<int> SplitPrimesUpTo(int hi)
    {
        for (int n = 5; n < hi; n++)
            if (n % 4 == 1 && IsPrime(n)) yield return n;
    }

    private static bool IsPrime(int n)
    {
        if (n < 2) return false;
        for (int d = 2; (long)d * d <= n; d++) if (n % d == 0) return false;
        return true;
    }
}
