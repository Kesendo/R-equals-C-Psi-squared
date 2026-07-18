namespace MirrorWorld;

// The witness reading (adopted 2026-07-18, F135 + F136): WHO records and WHAT it records, read off
// a pair page in closed form -- no propagator, no eigensolver, pure cosine arithmetic. The substrate
// is the Absorption Theorem's diagonal world (H = sum Delta_ab Z_a Z_b, local Z-dephasing, |+>^N):
// every pair page is a product of a write phase, the pair's own dephasing factors, and one cosine per
// third site (F135 Proposition 1). At the symmetric readout t* = pi/(4 Delta_S) the parities decide
// everything (F136, the record letter law):
//   - private watchers of the witness must be EVEN, always;
//   - the shared dressers' parity selects the family: all even + write bond -> POINTER record
//     (1 bit of Z_S in the witness's equator, the ZY channel), all odd -> BELL record
//     (1/2 Phi + 1/2 Psi, zero pointer content), mixed -> dark;
//   - the Bell LETTER alternates with the dresser count: m odd -> Y(x)Y, m even -> X(x)X,
//     sign(c1 c2) = (-1)^m, signs in closed form (the signed-coherence corollary made observable);
//   - the hinge V_S(t*) = 0 (S keeps a neighbor besides j) is global; a pure pendant S ROLE-SWAPS
//     into its neighbor's witness (1 separable YZ bit at odd watchers, 2 entangled bits at even);
//   - prices: the Bell record pays e^{-2(gamma_S+gamma_j)t*} (BOTH sites), the pointer record only
//     the witness's own gamma, and gamma on any traced site (dressers, watchers) is exactly invisible.
// At uniform coupling the statics are Hein-Eisert-Briegel 2004 graph-state marginals composed with
// the evolution's degree rotation (credited + gated in the proof); the ratio arithmetic, the signs,
// the exclusivity and the gamma dressing are the arc's own. Everything here restates
// docs/proofs/PROOF_RECORD_LETTER_LAW.md (87/87 gates) and PROOF_RECORD_PARITY_LAW.md; the J = 0
// boundary and all derivations live there, not here.
public sealed class Witness : GameObject
{
    public enum Kind { Pointer, Bell, RoleSwap, Entangled, Dark, Generic }

    static readonly double TStar = Math.PI / 4.0;   // t* at Delta_S = 1 (all couplings in Delta_S units)

    public int S { get; }
    public int J { get; }
    public Kind Family { get; }
    public char Letter { get; }                     // the luminous channel's witness letter: 'Y'/'X' (Bell), 'Y' via ZY (Pointer), 'Z' via YZ (RoleSwap), '-' else
    public int Sign { get; }                        // the luminous correlator's sign (+1/-1), 0 when dark/generic
    public double Bits { get; }                     // I(S:j) in bits, closed form (NaN = not classified here)
    public int Dressers { get; }                    // m = |D|
    public bool WriteBond { get; }

    public Witness(World world, IReadOnlyList<(int a, int b, double delta)> bonds, int s, int j,
                   double gammaS = 0.0, double gammaJ = 0.0) : base(world)
    {
        S = s; J = j;
        var nbrS = new Dictionary<int, double>();
        var nbrJ = new Dictionary<int, double>();
        foreach (var (a, b, delta) in bonds)
        {
            if (a == s) nbrS[b] = delta; else if (b == s) nbrS[a] = delta;
            if (a == j) nbrJ[b] = delta; else if (b == j) nbrJ[a] = delta;
        }
        WriteBond = nbrS.ContainsKey(j);
        var d = nbrJ.Keys.Where(k => k != s && nbrS.ContainsKey(k)).OrderBy(k => k).ToArray();   // shared dressers
        var p = nbrJ.Keys.Where(k => k != s && !nbrS.ContainsKey(k)).OrderBy(k => k).ToArray();  // private watchers of j
        var q = nbrS.Keys.Where(k => k != j && !nbrJ.ContainsKey(k)).ToArray();                  // private S-neighbors
        Dressers = d.Length;
        bool hinge = d.Length + q.Length > 0;                       // V_S(t*) = 0: S keeps a neighbor besides j

        var ratios = d.Concat(p).Select(k => nbrJ[k]).ToArray();    // r_k = Delta_jk / Delta_S (Delta_S = 1)
        bool normalized = nbrS.Values.All(v => v == 1.0);           // the class's contract: all couplings in Delta_S units
        bool allInteger = ratios.All(r => IsInt(r)) && normalized;
        bool dOdd = d.All(k => IsOddInt(nbrJ[k])) && d.Length >= 1;
        bool dEven = d.All(k => IsEvenInt(nbrJ[k]));                // vacuously true at D = empty
        bool pEven = p.All(k => IsEvenInt(nbrJ[k]));
        bool pAnyOdd = p.Any(k => IsOddInt(nbrJ[k]));               // EXISTENTIAL: one odd watcher kills beta_j AND the doubles

        if (!hinge)
        {
            // the pendant role-swap: S's only bond is the write bond; the pair reads BACKWARDS.
            // One odd-integer watcher suffices (it zeroes j's radius and both double coherences at
            // once); the other watchers may be anything. All even (or none) is the entangled face.
            if (WriteBond && normalized && pAnyOdd)
            { Family = Kind.RoleSwap; Letter = 'Z'; Sign = +1; Bits = OneBit(Math.Exp(-2.0 * gammaS * TStar)); }
            else if (WriteBond && normalized && pEven && allInteger)
            { Family = Kind.Entangled; Letter = '-'; Sign = 0;
              Bits = gammaS == 0.0 && gammaJ == 0.0 ? 2.0 : double.NaN; }   // the bare-dimer face is pinned at gamma = 0 only
            else { Family = Kind.Generic; Letter = '-'; Sign = 0; Bits = double.NaN; }
            return;
        }

        if (allInteger && WriteBond && dEven && pEven && nbrS.Count >= 2)
        {
            // POINTER record: every watcher of j even, shared and private alike; sign Prod (-1)^{r/2}.
            Family = Kind.Pointer; Letter = 'Y';
            Sign = d.Concat(p).Aggregate(1, (acc, k) => acc * ParityPow(nbrJ[k] / 2.0));
            Bits = OneBit(Math.Exp(-2.0 * gammaJ * TStar));
            return;
        }
        if (allInteger && q.Length == 0 && dOdd && pEven)
        {
            // BELL record: the letter is the dresser parity, the sign is sigma_2.
            Family = Kind.Bell; Letter = Dressers % 2 == 1 ? 'Y' : 'X';
            int sigma2 = d.Aggregate(1, (acc, k) => acc * ParityPow((1.0 - nbrJ[k]) / 2.0))
                       * p.Aggregate(1, (acc, k) => acc * ParityPow(nbrJ[k] / 2.0));
            Sign = sigma2;
            Bits = OneBit(Math.Exp(-2.0 * (gammaS + gammaJ) * TStar));
            return;
        }
        if (allInteger) { Family = Kind.Dark; Letter = '-'; Sign = 0; Bits = 0.0; return; }
        // Generic = NOT CLASSIFIED HERE, not "I < 1": the proof's fully-shared m = 1 corner (write
        // bond + one dresser + nothing else, non-integer ratio) holds a FORCED full bit with a
        // rotating channel and is deliberately left Generic by this class (proof Edge cases, gated).
        Family = Kind.Generic; Letter = '-'; Sign = 0; Bits = double.NaN;
    }

    // I = 1 - h2((1+kappa)/2): the record formula shared by both classical-quantum families.
    static double OneBit(double kappa)
    {
        double pPlus = (1.0 + kappa) / 2.0;
        return 1.0 - H2(pPlus);
    }

    static double H2(double x) =>
        x <= 0.0 || x >= 1.0 ? 0.0 : -x * Math.Log2(x) - (1.0 - x) * Math.Log2(1.0 - x);

    static bool IsInt(double r) => Math.Abs(r - Math.Round(r)) < 1e-12;
    static bool IsOddInt(double r) => IsInt(r) && ((long)Math.Round(r)) % 2 != 0;
    static bool IsEvenInt(double r) => IsInt(r) && ((long)Math.Round(r)) % 2 == 0;
    static int ParityPow(double e) => ((long)Math.Round(e)) % 2 == 0 ? +1 : -1;   // (-1)^e for integer e

    public override IReadOnlyList<string> Own => new[] { "family", "letter", "sign", "bits" };
}
