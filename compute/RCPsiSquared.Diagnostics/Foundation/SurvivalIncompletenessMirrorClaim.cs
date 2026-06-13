using System.Globalization;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The survival law and the V-Effect/incompleteness are INVERSION-MIRROR PARTNERS on the
/// Pi2 dyadic ladder, and the longest-lived dissipative mode dynamically lives in the incompleteness
/// region on dispersive matter (Tier1Candidate: the algebra is derived + live, the dynamic enactment
/// is verified at N&lt;=7, not proven). The thread (a)+(b) result of the
/// <c>survival_incompleteness_mirror</c> arc; the C# witness is
/// <see cref="IncompletenessSurvivorWitness"/> (<c>inspect --root survivor</c>).
///
/// <para><b>(a) THE ALGEBRA (derived, live).</b> On the Pi2 dyadic ladder a_n = 2^(1-n): the survival
/// law's quantum is <c>a_0 = 2</c> (<see cref="AbsorptionTheoremClaim"/>; Re(lambda) = -2*gamma*&lt;n_XY&gt;,
/// the 2 == the qubit dimension d), and the V-Effect/incompleteness baseline is <c>a_2 = 1/2</c>
/// (<see cref="HalfAsStructuralFixedPointClaim"/> Face 1, C=1/2 == 1/d). They are inversion-mirror
/// partners: <c>MirrorPartnerIndex(0)=2</c>, <c>a_0*a_2 = d*(1/d) = 1</c>, with the self-mirror a_1=1
/// sitting exactly between them, all forced by d^2-2d=0. So "the V-Effect inherits algebraically" is
/// literally typed: it shares the one ladder with the survival law, as its reciprocal mirror. Survival
/// obeys the same law model-independently (the Absorption Theorem; the XY-vs-Heisenberg reconciliation).</para>
///
/// <para><b>(b) THE DYNAMIC ENACTMENT (verified).</b> Below the handover Q the longest-lived mode lives
/// in the incompleteness region (interior filling) on DISPERSIVE EXTENDED matter - the labels are
/// physical: CHAIN = conjugated polyene chains / spin chains / the Grotthuss proton wire (dead-centre
/// half-filling survivor); RING = aromatic rings / light-harvesting macrocycles (off-centre interior
/// (2,2)/(N-2,N-2)). The HUB-LOCALIZED central-spin STAR (NV centre / quantum dot / the mediator) is the
/// COUNTEREXAMPLE: its survivor is the boundary hub coherence (1,1)/(N-1,N-1), because a hub-spoke bath
/// has no dispersion and thus no central momentum mode. Lifetime &lt;n_XY&gt; ~ c*Q^2/N^2, ring/chain -&gt; 4
/// (the cyclic-vs-open k_min^2 ratio, model-independent, the SAME 4x CHAIN_GAP reports for Heisenberg) -
/// a SEPARATE 1/N^2 magnon-admixture inheritance from the Pi2 dyadic CONSTANT ladder, the two meeting at
/// the per-mode Absorption Theorem (a_0).</para>
///
/// <para><b>STILL OPEN</b> (from <see cref="VacuumBlockReductionClaim"/>): whether the V-Effect
/// (Pauli-string weight w=N/2 self-pair) EQUALS the {0,2}-coherence (n_diff) survivor is a different
/// decomposition (weight vs bra-ket disagreement count); this claim is the constant-mirror + the
/// light-survivor, NOT that finer identity.</para>
///
/// <para>Typed parents: <see cref="AbsorptionTheoremClaim"/> (a_0, the survival law) +
/// <see cref="HalfAsStructuralFixedPointClaim"/> (a_2, the incompleteness baseline).
/// Anchors: <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> (the a_0&lt;-&gt;a_2 mirror) +
/// <c>simulations/carbon/incompleteness_survivor.py</c> + <c>simulations/carbon/_xy_vs_heisenberg_slowmode.py</c> +
/// <c>compute/RCPsiSquared.Diagnostics/Foundation/IncompletenessSurvivorWitness.cs</c>.</para></summary>
public sealed class SurvivalIncompletenessMirrorClaim : Claim
{
    /// <summary>Typed parent a_0: the survival law (Absorption Theorem), the n=0 ladder anchor.</summary>
    public AbsorptionTheoremClaim Survival { get; }

    /// <summary>Typed parent a_2: the V-Effect/incompleteness baseline (C=1/2), the n=2 ladder anchor.</summary>
    public HalfAsStructuralFixedPointClaim Incompleteness { get; }

    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    public IReadOnlyList<BatteryCase> Cases { get; }
    public int PassCount => Cases.Count(c => c.Passes);

    public SurvivalIncompletenessMirrorClaim(AbsorptionTheoremClaim survival, HalfAsStructuralFixedPointClaim incompleteness)
        : base("Survival mirrors incompleteness: the survival law (a_0=2=AbsorptionTheorem, the 2gamma quantum = qubit dim d) " +
               "and the V-Effect/incompleteness (a_2=1/2=Half, C=1/2 = 1/d) are inversion-mirror partners on the Pi2 dyadic " +
               "ladder (a_0*a_2 = d*(1/d) = 1, self-mirror a_1=1 between, forced by d^2-2d=0) - the algebra of 'der V-Effekt " +
               "vererbt sich', DERIVED+live. Dynamically VERIFIED (N<=7): below the handover Q the longest-lived mode lives in " +
               "the incompleteness region on DISPERSIVE extended matter (chain at dead-centre half-filling, ring off-centre); " +
               "the hub-localized central-spin STAR is the boundary counterexample. Lifetime <n_XY> ~ Q^2/N^2, ring/chain->4 " +
               "(model-independent), a separate 1/N^2 inheritance from the dyadic constants. OPEN: the Pauli-weight w=N/2 vs " +
               "n_diff {0,2} V-Effect identity (VacuumBlockReduction).",
               Tier.Tier1Candidate,
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "simulations/carbon/incompleteness_survivor.py + " +
               "simulations/carbon/_xy_vs_heisenberg_slowmode.py + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/IncompletenessSurvivorWitness.cs")
    {
        Survival = survival ?? throw new ArgumentNullException(nameof(survival));
        Incompleteness = incompleteness ?? throw new ArgumentNullException(nameof(incompleteness));
        Cases = BuildBattery(survival);
    }

    /// <summary>a_0 = 2 (the survival quantum, == the qubit dimension d), live from the ladder.</summary>
    public double A0 => Survival.Ladder.Term(0);

    /// <summary>a_2 = 1/2 (the incompleteness baseline, == 1/d), live from the ladder.</summary>
    public double A2 => Survival.Ladder.Term(2);

    /// <summary>The mirror partner index of a_0: 2 (so a_0's partner IS a_2).</summary>
    public int MirrorPartnerOfSurvival => Survival.Ladder.MirrorPartnerIndex(0);

    /// <summary>a_0 * a_2 = d*(1/d) = 1 (the inversion-mirror identity), live.</summary>
    public double MirrorProduct => Survival.Ladder.ProductWithMirrorPartner(0);

    public override string DisplayName =>
        "Survival mirrors incompleteness: a_0 (2gamma) and a_2 (C=1/2) are Pi2-ladder mirror partners (Tier1Candidate)";

    public override string Summary =>
        $"the survival law (a_0={A0:0.#}) and the V-Effect/incompleteness (a_2={A2:0.#}) are inversion-mirror partners " +
        $"(a_0*a_2={MirrorProduct:0.#}=d*1/d); dynamically the longest-lived mode is the interior incompleteness coherence on " +
        $"dispersive matter (chain/ring), the central-spin star the boundary counterexample; {PassCount}/{Cases.Count} PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("(a) the algebra: a_0 <-> a_2 mirror",
                summary: $"a_0={A0:0.#} (survival quantum 2gamma == qubit dim d) and a_2={A2:0.#} (incompleteness C=1/2 == 1/d) " +
                         $"are Pi2-ladder inversion-mirror partners: MirrorPartnerIndex(0)={MirrorPartnerOfSurvival}, " +
                         $"a_0*a_2={MirrorProduct:0.#}=d*(1/d), self-mirror a_1=1 between, forced by d^2-2d=0.");
            yield return new InspectableNode("(b) the dynamic enactment (the incomplete survives on dispersive matter)",
                summary: "below the handover Q the longest-lived mode is the interior incompleteness coherence on dispersive " +
                         "extended matter (chain dead-centre half-filling, ring off-centre); the hub-localized central-spin star " +
                         "is the boundary counterexample. Lifetime <n_XY> ~ Q^2/N^2, ring/chain -> 4 (model-independent). " +
                         "Live: inspect --root survivor (IncompletenessSurvivorWitness).");
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
            yield return Survival;          // typed parent a_0
            yield return Incompleteness;    // typed parent a_2
        }
    }

    private static string Classify(int n, int pc, int pr)
    {
        if (pc != pr) return "band-edge";
        if (n % 2 == 0 && pc == n / 2) return "centre";
        if (pc >= 2 && pc <= n - 2) return "off-centre interior";
        if (pc == 1 || pc == n - 1) return "boundary";
        return "?";
    }

    private static IReadOnlyList<BatteryCase> BuildBattery(AbsorptionTheoremClaim survival)
    {
        string R(double v) => v.ToString("0.######", CultureInfo.InvariantCulture);
        var L = survival.Ladder;
        var cases = new List<BatteryCase>();

        // (a) the algebra, live from the ladder.
        cases.Add(new BatteryCase("a_0 = 2 (survival quantum = qubit dim d)",
            "Pi2 ladder Term(0)", "2", R(L.Term(0))));
        cases.Add(new BatteryCase("a_2 = 1/2 (incompleteness baseline = 1/d)",
            "Pi2 ladder Term(2)", "0.5", R(L.Term(2))));
        cases.Add(new BatteryCase("a_0's mirror partner IS a_2 (index 2)",
            "MirrorPartnerIndex(0)", "2", L.MirrorPartnerIndex(0).ToString(CultureInfo.InvariantCulture)));
        cases.Add(new BatteryCase("a_0 * a_2 = d*(1/d) = 1 (the inversion mirror)",
            "ProductWithMirrorPartner(0)", "1", R(L.ProductWithMirrorPartner(0))));

        // (b) the dynamic enactment via the witness (carbon XY, N=6, Q=1.5, below the handover).
        var chain = IncompletenessSurvivorWitness.Survivor(6, 1.5, TopologyKind.Chain);
        var ring = IncompletenessSurvivorWitness.Survivor(6, 1.5, TopologyKind.Ring);
        var star = IncompletenessSurvivorWitness.Survivor(6, 1.5, TopologyKind.Star);
        cases.Add(new BatteryCase("chain survivor = interior centre (dispersive)",
            $"N=6 Q=1.5 -> sector ({chain.PCol},{chain.PRow})", "centre", Classify(6, chain.PCol, chain.PRow)));
        cases.Add(new BatteryCase("ring survivor = off-centre interior (dispersive)",
            $"N=6 Q=1.5 -> sector ({ring.PCol},{ring.PRow})", "off-centre interior", Classify(6, ring.PCol, ring.PRow)));
        cases.Add(new BatteryCase("star survivor = boundary (central-spin counterexample)",
            $"N=6 Q=1.5 -> sector ({star.PCol},{star.PRow})", "boundary", Classify(6, star.PCol, star.PRow)));
        double ratio = chain.NXy > 1e-9 ? ring.NXy / chain.NXy : 0.0;
        bool ratioOk = ratio >= 3.5 && ratio <= 4.5;
        cases.Add(new BatteryCase("ring/chain lifetime ratio ~ 4 (cyclic-vs-open k_min^2)",
            $"ring <n_XY>/chain <n_XY> = {R(ratio)}", "~4", ratioOk ? "~4" : R(ratio)));

        return cases;
    }

    /// <summary>The shared instance with its parents built from the Pi2-Foundation root, mirroring
    /// <see cref="VacuumBlockReductionClaim.Shared"/>: lets metadata be read without the full registry.</summary>
    public static SurvivalIncompletenessMirrorClaim Shared { get; } =
        new SurvivalIncompletenessMirrorClaim(
            new AbsorptionTheoremClaim(new Pi2DyadicLadderClaim()),
            new HalfAsStructuralFixedPointClaim());
}
