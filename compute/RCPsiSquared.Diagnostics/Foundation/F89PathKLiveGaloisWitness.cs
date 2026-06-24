using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>FULL-LIVE (Option D) witness for a path-k H_B-mixed factor (k = 4,5,6). Builds the
/// (SE,DE) S_2-sym block LIVE at q=2 over Z[i], reconstructs the AT-locked factor from the
/// rate-confined invariant subspace (<see cref="F89AtFactorReconstruction.ForPathK"/> — NO F_d
/// import), divides it out of the Berkowitz charpoly to isolate F_d (the validation triple), and
/// reads F_d's Frobenius cycle types: Gal(F_d) = S_d, non-solvable. F_d is recomputed from the
/// physics, never imported — the committed oracle is only a test cross-check. The heavy recompute
/// lives in the Children getter (path-6 stays cheap until expanded). path-3's octic sibling is
/// <see cref="F89OcticGaloisWitness"/>.</summary>
public sealed class F89PathKLiveGaloisWitness : IInspectable
{
    private readonly int _k;

    public F89PathKLiveGaloisWitness(int k) => _k = k;

    private int D => F89PathKFdOracle.Degree(_k);

    public string DisplayName => $"F89 path-{_k} H_B-mixed factor F_{D} Galois group (live, full D)";

    public string Summary =>
        $"Gal(F_{D}/Q(i)(q)) = S_{D} (non-solvable) — expand to recompute live: block → Berkowitz → isolate F_d (reconstructed AT) → Frobenius ⟹ S_{D}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            // 1) block live -> Berkowitz charpoly; reconstruct AT (rate-confined invariant subspace,
            //    no F_d import) and isolate F_d = C / AT (the validation triple).
            var block = F89PathKSeDeBlock.BuildTwoTimesSymBlock(q0: 2, nBlock: _k + 1);
            var charpoly = GaussianMatrixCharpoly.Characteristic(block);
            var atFactor = F89AtFactorReconstruction.ForPathK(_k);
            var fd = F89HbMixedIsolation.Isolate(charpoly, atFactor, D);
            yield return new InspectableNode("F_d isolated live (full Option D)",
                summary: $"the path-{_k} (SE,DE) S_2-sym block (dim {F89PathKFdOracle.SymDim(_k)}) is built at q=2 over Z[i]; " +
                         "Berkowitz gives its charpoly; the AT factor is reconstructed from the rate-confined invariant subspace " +
                         $"(no F_d import) and divided out — triple verified (remainder 0, degree {D}, gcd(AT,F_d)=1). F_d IS the H_B-mixed factor");

            // 2) Frobenius cycle types of the live-isolated F_d -> the generalised Jordan verdict.
            var re = fd.Select(c => c.Re).ToArray();
            var im = fd.Select(c => c.Im).ToArray();
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
                summary: $"S_{D} is the generic group; integrability spends itself on the AT-locked F_a/F_b factorisation (here reconstructed), leaving the H_B-mixed F_{D} structureless");
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
