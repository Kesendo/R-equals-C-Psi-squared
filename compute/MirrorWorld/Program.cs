using System.Globalization;
using MirrorWorld;

CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;   // dots, not commas

// ---- run mode "grow": the diagonal protocol, step 1 (the world splits) ----
// ClaudeTasks/DIAGONAL_PROTOCOL_GAME.md, rules 1-3 running in time. No Hamiltonian yet: a field of
// possibilities under the one question, structure (the diagonal) staying while novelty (off-diagonal)
// is culled by its disagreement. The mirror is visible in the rate ladder (k decays, N-k its twin).
if (args.Length > 0 && args[0] == "grow")
{
    const int gn = 2;
    const double gg = 0.5, dt = 0.05;
    const int ticks = 20;
    var gworld = new World();
    var field = new Field(gworld, gn, gg);
    field.SeedUniform();

    Console.WriteLine("the world splits (diagonal protocol, step 1 -- rules 1-3 in time)");
    Console.WriteLine($"  source ClaudeTasks/DIAGONAL_PROTOCOL_GAME.md; N={gn}, gamma={gg}, dt={dt}");
    Console.WriteLine("  rule 2, the one question: how much do two possibilities disagree? more disagreement, faster fade.");
    Console.WriteLine("  structure = the diagonal (k=0) that stays; novelty = the off-diagonal that is culled.");
    Console.WriteLine($"  the knower's cut: {field.AliveCount} alive (upper triangle, stepped), {field.ImmortalCount} immortal (k=0, held), {field.AliveCount} mirrored (rho=rho-dagger) -- only {field.AliveCount} of {field.Dim * field.Dim} cells run.");
    Console.WriteLine($"  {"t",4} {"structure",9} {"novelty",9}   {string.Join(" ", Enumerable.Range(0, gn + 1).Select(k => $"k={k}".PadLeft(7)))}");
    for (int tick = 0; tick <= ticks; tick++)
    {
        var byK = field.WeightByDisagreement();
        Console.WriteLine($"  {field.T,4:0.00} {field.Structure,9:0.000} {field.Novelty,9:0.000}   {string.Join(" ", byK.Select(v => v.ToString("0.000").PadLeft(7)))}");
        field.Step(dt);
    }
    Console.WriteLine("  structure held flat (k=0 costs nothing); novelty fell, k=2 twice as fast as k=1 -- the mirror in the rates.");
    Console.WriteLine("  (step 2 -- run `-- live` -- adds the inner restlessness: novelty BORN from structure.)");
    return;
}

// ---- run mode "live": the diagonal protocol, step 2 (the world lives) ----
// The empty world's watching plus the inner restlessness: H (a flip-flop on the bond) makes coherence out
// of population -- novelty BORN from structure (rule 4) -- while the watching culls it. Novelty no longer
// only fades; it is born and culled, and the populations relax to a living balance.
if (args.Length > 0 && args[0] == "live")
{
    int gn = args.Length > 1 ? int.Parse(args[1]) : 2;
    string gtopo = args.Length > 2 ? args[2] : "chain";
    const double gj = 1.0, gg = 0.5, dt = 0.05;
    const int ticks = 40;
    var lworld = new World();
    var rest = new Restless(lworld, gn, gj, gg, Topology.Named(gtopo, gn));
    rest.Seed(1);                                   // |0...01>, one excitation at one end (block (1,1))

    Console.WriteLine("the world lives (diagonal protocol, step 2 -- the inner restlessness)");
    Console.WriteLine($"  source ClaudeTasks/DIAGONAL_PROTOCOL_GAME.md; N={gn}, topology={gtopo}, J={gj}, gamma={gg}, dt={dt}");
    Console.WriteLine("  rule 4, the inner motion: H (the flip-flop bond) makes coherence out of population --");
    Console.WriteLine("  novelty BORN from structure, then culled by the watching. seed = |0...01>, one excitation.");
    Console.WriteLine($"  cut (c): the loop stays in joint-popcount block (1,1) -- {rest.AliveCount} of {rest.Dim * rest.Dim} cells run, {rest.ForbiddenCount} forbidden (F63). the k=1 coherences are among them, never computed.");
    Console.WriteLine($"  {"t",4} {"structure",9} {"novelty",9}   {string.Join(" ", Enumerable.Range(0, gn + 1).Select(k => $"k={k}".PadLeft(7)))}");
    for (int tick = 0; tick <= ticks; tick++)
    {
        var byK = rest.WeightByDisagreement();
        Console.WriteLine($"  {rest.T,4:0.00} {rest.Structure,9:0.000} {rest.Novelty,9:0.000}   {string.Join(" ", byK.Select(v => v.ToString("0.000").PadLeft(7)))}");
        rest.Step(dt);
    }
    Console.WriteLine("  novelty rose from 0 (born by H from the population), churned, then damped -- the living");
    Console.WriteLine("  balance: the excitation spreads over the chain toward 1/N per site (structure persists, redistributed).");
    return;
}

// ---- run mode "concentrator": the site-resolved rate cross-check (IBM Concentrator Reloaded 7a) ----
// A standalone continuous-Lindblad comparator for experiments/IBM_CONCENTRATOR_RELOADED.md stage 7a. The
// flown design Trotterizes the isotropic Heisenberg bond exp(-i*0.05*(XX+YY+ZZ)) with a per-step Z-sink of
// off-diagonal retention e^-0.05; here we integrate the EXACT continuous Lindblad rho-dot = -i[H,rho] +
// sum_l gamma_l D_Zl[rho] to the SAME depth grid and read the payload coherence 2|<0..0|rho|0..0,1@payload>|
// (= |<X_p>+i<Y_p>|). Three arms -- 0 (no sink), E (far-edge sink), MP (on-payload sink) -- the paired-ratio
// ln-slope over the grid, and the verdict slope(MP)-slope(E) vs the pre-registered central -0.07325/step.
// Convention, read from Restless + Pair and pinned in ConcentratorTests: hopping J=2 gives the unit-Pauli XY
// bond (XX+YY), zz=1 the isotropic Heisenberg ZZ; gamma_sink=0.5 gives per-step retention e^-2*gamma*dt =
// e^-0.05 at dt=0.05 (Pair rate -2*gamma*k). Runs at N=5 (the flight) and past the N=8 wall at N=9 (bonus).
if (args.Length > 0 && args[0] == "concentrator")
{
    int cn = args.Length > 1 ? int.Parse(args[1]) : 5;
    const double jHop = 2.0, gSink = 0.5, stepT = 0.05;
    const int subPerStep = 50;                       // RK4 substeps per Trotter-step interval (continuous limit)
    int[] grid = { 1, 2, 3, 4, 6, 8 };
    int[] gridNoDeep = { 1, 2, 3, 4 };
    int payload = cn / 2;                            // center site (N=5 -> 2, N=9 -> 4)
    int edge = 0;                                    // the far-edge sink site (design arm E)
    var cworld = new World();

    double[] Arm(int? sinkSite, double zz)
    {
        var g = new double[cn];
        if (sinkSite is int ss) g[ss] = gSink;
        var r = new Restless(cworld, cn, jHop, 0.0, Topology.Chain(cn), siteGammas: g, zz: zz);
        int pIdx = 1 << payload;
        r.Seed(0, 0.5); r.Seed(pIdx, 0.5); r.SeedCoherence(0, pIdx, 0.5);   // |+> on payload, |0> elsewhere
        double dt = stepT / subPerStep;
        var coh = new List<double>();
        int gi = 0;
        for (int step = 1; step <= grid[^1]; step++)
        {
            for (int s = 0; s < subPerStep; s++) r.Step(dt);
            if (gi < grid.Length && step == grid[gi]) { coh.Add(2.0 * r[0, pIdx].Magnitude); gi++; }
        }
        return coh.ToArray();
    }

    static double SlopeLnR(double[] ca, double[] c0, int[] gr)
    {
        int m = gr.Length; double sx = 0, sy = 0, sxx = 0, sxy = 0;
        for (int i = 0; i < m; i++)
        {
            double x = gr[i], y = Math.Log(ca[i] / c0[i]);   // x = step number -> slope is PER STEP
            sx += x; sy += y; sxx += x * x; sxy += x * y;
        }
        return (m * sxy - sx * sy) / (m * sxx - sx * sx);
    }
    static double[] Pick(double[] full, int[] grid, int[] want)
    {
        var outp = new double[want.Length];
        for (int i = 0; i < want.Length; i++) outp[i] = full[Array.IndexOf(grid, want[i])];
        return outp;
    }
    static string Sg(double v) => (v < 0 ? "-" : "+") + Math.Abs(v).ToString("0.000000");   // clean forced sign

    Console.WriteLine("the site-resolved rate cross-check (IBM Concentrator Reloaded 7a, continuous-Lindblad)");
    Console.WriteLine($"  source experiments/IBM_CONCENTRATOR_RELOADED.md; N={cn}, payload |+> on site {payload}, others |0>");
    Console.WriteLine($"  H = sum_bonds (XX+YY+ZZ) [hopping J={jHop} -> unit-Pauli XY, zz=1 isotropic]; watching D_Z, gamma_sink={gSink}");
    Console.WriteLine($"  arm 0: no sink; arm E: sink on far edge site {edge}; arm MP: sink on payload site {payload}");
    Console.WriteLine($"  observable coh_a(t) = 2|rho[0, 1<<{payload}]| = |<X_{payload}>+i<Y_{payload}>|; grid (Trotter steps) = [{string.Join(",", grid)}], t = 0.05*step");
    Console.WriteLine($"  continuous limit: {subPerStep} RK4 substeps per step (dt = {stepT / subPerStep})");
    Console.WriteLine();

    var c0 = Arm(null, 1.0);
    var cE = Arm(edge, 1.0);
    var cMP = Arm(payload, 1.0);
    Console.WriteLine($"  {"step",4} {"coh_0",9} {"coh_E",9} {"coh_MP",9} {"R_E=E/0",9} {"R_MP=MP/0",10}");
    for (int i = 0; i < grid.Length; i++)
        Console.WriteLine($"  {grid[i],4} {c0[i],9:0.00000} {cE[i],9:0.00000} {cMP[i],9:0.00000} {cE[i] / c0[i],9:0.00000} {cMP[i] / c0[i],10:0.00000}");

    double sE = SlopeLnR(cE, c0, grid), sMP = SlopeLnR(cMP, c0, grid), diff = sMP - sE;
    double sE4 = SlopeLnR(Pick(cE, grid, gridNoDeep), Pick(c0, grid, gridNoDeep), gridNoDeep);
    double sMP4 = SlopeLnR(Pick(cMP, grid, gridNoDeep), Pick(c0, grid, gridNoDeep), gridNoDeep);
    Console.WriteLine();
    Console.WriteLine($"  slope(E)  = {sE,+10:0.000000}/step   slope(MP) = {sMP,+10:0.000000}/step");
    Console.WriteLine($"  slope(MP) - slope(E) = {diff,+10:0.000000}/step   (no-deep, steps 1-4: {sMP4 - sE4,+10:0.000000})");
    Console.WriteLine($"  coh_0 at step 8 (t=0.4) = {c0[^1]:0.00000}   (design clean anchor ~0.456)");

    if (cn == 5)
    {
        // isolate the ZZ contribution: rerun XY-only (zz=0) and compare the contrast
        var c0x = Arm(null, 0.0); var cEx = Arm(edge, 0.0); var cMPx = Arm(payload, 0.0);
        double diffXY = SlopeLnR(cMPx, c0x, grid) - SlopeLnR(cEx, c0x, grid);
        Console.WriteLine();
        Console.WriteLine("  cross-check vs the pre-registered design (experiments/IBM_CONCENTRATOR_RELOADED.md 7a):");
        Console.WriteLine($"    design flown Trotter (J*dt=0.05)  slope(MP)-slope(E) = -0.073249/step (recorded -0.07325)");
        Console.WriteLine($"    this continuous Heisenberg (zz=1) slope(MP)-slope(E) = {Sg(diff)}/step");
        Console.WriteLine($"    Trotter gap (continuous - flown)  = {Sg(diff - (-0.073249))}/step  (O((J*dt)^2) ~ 0.0025 bound)");
        Console.WriteLine($"    XY-only (zz=0) contrast           = {Sg(diffXY)}/step  -> ZZ shifts it by {Sg(diff - diffXY)} (tiny)");
        Console.WriteLine($"    far-sink leakage slope(E)         = {Sg(sE)}/step  (design -0.00029, null-consistency)");
    }
    else
    {
        Console.WriteLine();
        Console.WriteLine($"  N={cn} is past the N=8 spectrum wall; the site contrast is still a small block-local run (100 alive cells).");
        Console.WriteLine($"    site contrast slope(MP)-slope(E) = {Sg(diff)}/step; far-sink slope(E) = {Sg(sE)}/step (leakage).");
    }
    return;
}

// ---- run mode "mirror": the first mirror in the world of mirrors ----
// The fold-lattice legs (adopted 2026-07-03) as exact entry-wise rearrangements: no eigensolver, the
// mirror checked cell by cell. Then the play: orbits, the self-folded price, and the trajectory fold
// (the mirror shows the partner block running backward, paid at the deepest rate 2*N*gamma).
if (args.Length > 0 && args[0] == "mirror")
{
    int mn = args.Length > 1 ? int.Parse(args[1]) : 5;
    const double mj = 1.0, mg = 0.5;
    var mworld = new World();
    var mirror = new Mirror(mworld, mn, mj, mg);

    Console.WriteLine("the first mirror in the world of mirrors (the fold lattice, adopted 2026-07-03)");
    Console.WriteLine($"  source docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md section 7; N={mn}, J={mj}, gamma={mg}, price=2*N*gamma={mirror.Price}");
    Console.WriteLine("  three legs on the block lattice, each an EXACT entry-wise identity (no eigensolver):");
    Console.WriteLine("    t (swap bra/ket, keep spectrum), f_P (flip bra), f_Q (flip ket) -- the folds pay lambda -> -lambda - price.");

    if (mn > 8)
    {
        // past the wall: the full lattice is out of reach (half-filling is exponential), but the
        // memory-cut pair is not: block (1,1) is N^2 and its fold partner (1,N-1) -- one hole against
        // one excitation -- is N^2 too. The mirror walks where the spectrum died at N=8.
        Console.WriteLine();
        Console.WriteLine($"  N={mn} is past the wall (the spectrum died at N=8). The memory-cut pair stays small:");
        var (res, dim) = mirror.PastTheWallResidual();
        Console.WriteLine($"  block (1,1) and its fold partner (1,{mn - 1}) -- one hole against one excitation -- are BOTH {dim} = N^2,");
        Console.WriteLine($"  and the fold leg holds cell by cell: worst residual {res:E1}  (the mirror does not care about the wall).");
        var (wts, wnx, wnw, wworst) = mirror.PastTheWallTrajectory(dt: 0.002, ticks: 100);
        Console.WriteLine();
        Console.WriteLine($"  the trajectory fold past the wall: x forward in (1,1), w BACKWARD in (1,{mn - 1}); |w|/|x| = exp(price*t):");
        Console.WriteLine($"  {"t",6} {"|x|",12} {"|w|",14} {"|w|/|x|",14} {"exp(price*t)",14}");
        foreach (int tick in new[] { 0, 25, 50, 75, 100 })
            Console.WriteLine($"  {wts[tick],6:0.000} {wnx[tick],12:0.000000} {wnw[tick],14:0.000E+0} {wnw[tick] / wnx[tick],14:0.000E+0} {Math.Exp(mirror.Price * wts[tick]),14:0.000E+0}");
        Console.WriteLine($"  two independent runs, worst relative mismatch {wworst:E1}.");
        Console.WriteLine("  the fold partner of the memory cut is again a memory cut: the mirror maps small blocks to small");
        Console.WriteLine("  blocks, so it runs at any N the cone runs at -- the wall was a property of the spectrum, not of the mirror.");
        return;
    }

    // 1. the legs, checked cell by cell over the whole lattice
    double wT = 0, wP = 0, wQ = 0, wK = 0;
    for (int p = 0; p <= mn; p++)
        for (int q = 0; q <= mn; q++)
        {
            wT = Math.Max(wT, mirror.TransposeResidual(p, q));
            wP = Math.Max(wP, mirror.BraFoldResidual(p, q));
            wQ = Math.Max(wQ, mirror.KetFoldResidual(p, q));
            wK = Math.Max(wK, mirror.KleinResidual(p, q));
        }
    Console.WriteLine($"  all {(mn + 1) * (mn + 1)} blocks: worst residual t={wT:E1}, f_P={wP:E1}, f_Q={wQ:E1}, Klein={wK:E1}  (exact zeros: the mirror is a rearrangement, not an approximation)");

    // 2. the orbit map: the lattice folds to a fundamental domain (~1/8)
    Console.WriteLine();
    Console.WriteLine($"  the lattice folds: {(mn + 1) * (mn + 1)} blocks -> {mirror.OrbitCount()} orbits (representative per cell; '-' marks one fold, spectrum paid):");
    for (int p = 0; p <= mn; p++)
    {
        var row = new List<string>();
        for (int q = 0; q <= mn; q++)
        {
            var (rp, rq, par) = mirror.Representative(p, q);
            row.Add($"{rp}{rq}{(par == 1 ? "-" : " ")}");
        }
        Console.WriteLine($"    p={p}:  {string.Join(" ", row)}");
    }

    // 3. the self-folded price (even N): a block that is its own mirror pays out of its own trace
    if (mn % 2 == 0)
    {
        Console.WriteLine();
        Console.WriteLine($"  self-folded blocks (q=N/2): trace L = -(price/2)*dim, no computation needed:");
        for (int p = 0; p <= mn; p++)
        {
            var (tr, law) = mirror.SelfFoldedTrace(p);
            Console.WriteLine($"    (p={p}, q={mn / 2}): trace = {tr,10:0.000}  law = {law,10:0.000}  match {(Math.Abs(tr - law) < 1e-9 ? "exact" : "BROKEN")}");
        }
    }

    // 4. the trajectory fold: the mirror shows the partner running backward, at the price
    Console.WriteLine();
    var (ts, nx, nw, worst) = mirror.TrajectoryFold(1, 2, dt: 0.01, ticks: 200);
    Console.WriteLine($"  the trajectory fold on block (1,2) -> partner ({1},{mn - 2}): x runs forward under L(1,2),");
    Console.WriteLine($"  w runs the partner BACKWARD (w-dot = -L(1,{mn - 2}) w); the mirror predicts w(t) = exp(price*t) * fold(x(t)).");
    Console.WriteLine($"  {"t",5} {"|x| (decays)",13} {"|w| (grows)",13} {"|w|/|x|",12} {"exp(price*t)",13}");
    foreach (int tick in new[] { 0, 40, 80, 120, 160, 200 })
        Console.WriteLine($"  {ts[tick],5:0.0} {nx[tick],13:0.000000} {nw[tick],13:0.000} {nw[tick] / nx[tick],12:0.000} {Math.Exp(mirror.Price * ts[tick]),13:0.000}");
    Console.WriteLine($"  two independent runs, compared tick by tick: worst relative mismatch {worst:E1} (RK4 tolerance).");
    Console.WriteLine("  the mirror trades decay for growth and the exchange rate is the deepest rate in the world.");
    Console.WriteLine("  (states and their mirrors live here; the paths -- braid, monodromy -- stay in the main repo: a catalog cannot hold a way.)");
    return;
}

// ---- run mode "seed": the within-block self-dual seed (the shadow's source, as a count) ----
// Mirror gave the BETWEEN-block folds; this gives the WITHIN-block self-duality they leave untouched:
// where a state meets the mirror's null (v^T v = 0), a defective seed -- the static source the shadow
// (a large projector norm) and the i^4=1 holonomy leave behind in the main repo. Held as a COUNT, no
// eigensolver (F89's nullity surplus over GF(p)): N-1 forced seeds at odd N (the unmirrorable middle
// seat), 0 at even N (every seat mirrors, nothing forced real). Adopted 2026-07-07.
if (args.Length > 0 && args[0] == "seed")
{
    int sn = args.Length > 1 ? int.Parse(args[1]) : 9;
    const double sg = 0.5;
    var sworld = new World();
    Console.WriteLine("the within-block self-dual seed (the shadow's source, as a count -- adopted 2026-07-07)");
    Console.WriteLine("  sources experiments/F89_SEED_EXISTENCE_REDUCTION.md + docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md");
    Console.WriteLine("  Mirror holds the between-block folds; this holds the within-block self-duality they leave: where");
    Console.WriteLine("  a state meets the mirror's null (v^T v = 0), a DEFECTIVE seed. The count, ranks only, no eigensolver:");
    Console.WriteLine("  surplus = r(0+) - r(inf) = the real strands forced off the axis = the seed count (the pencil A + q*C).");
    Console.WriteLine();
    Console.WriteLine($"  {"N",3} {"dim",5} {"r(inf)",7} {"r(0+)",6} {"surplus",8} {"N-1",4} {"Z3",3} {"CJ",4} {"parity",7}   rungs (n_diff: dim, nullity)");
    for (int n = 3; n <= sn; n++)
    {
        var seed = new Seed(sworld, n, sg);
        var (rInf, r0, surplus, parts) = seed.Count();
        int dim = parts.Sum(pp => pp.Dim);
        string par = n % 2 == 1 ? "ODD" : "EVEN";
        string ps = string.Join("  ", parts.Select(pp => $"[{pp.NDiff}: {pp.Dim},{pp.Nullity}]"));
        string flag = surplus == n - 1 ? "= N-1" : surplus == 0 ? "ZERO" : "??";
        int z3 = Seed.Z3FromRInf(rInf);   // throws loudly if r(inf) ever stops being 3*Z3
        string cj = n % 2 == 1 ? (Seed.Z3ClosedFormOdd(n) == z3 ? "ok" : "??") : "-";
        Console.WriteLine($"  {n,3} {dim,5} {rInf,7} {r0,6} {surplus,8} {n - 1,4} {z3,3} {cj,4} {par,7}   {ps}  <- {flag}");
    }
    Console.WriteLine();
    Console.WriteLine("  N-1 seeds at odd N (the unmirrorable middle seat forces each half to face itself, eigenvalues");
    Console.WriteLine("  real, then collide -- Kato); 0 at even N (every seat mirrors, no real strand forced off). The between-block");
    Console.WriteLine("  fold (Mirror) holds exactly where the within-block eigenframe tears (v^T v -> 0) -- one merge, two");
    Console.WriteLine("  faces: what the way leaves behind lives here; the way itself (the shadow, the holonomy) stays outside.");
    Console.WriteLine();
    Console.WriteLine("  the resonance count, closed (adopted 2026-07-12): r(inf) = 3*Z3 at EVERY N (cyclotomic, exact);");
    Console.WriteLine("  at odd N the Conway-Jones form Z3 = (N-1)/2 + [3|n](n/3-2) + 2[15|n] (CJ column; ROT3 multiplicity");
    Console.WriteLine($"  verified, not derived). RESONANT (odd N) iff 3 | N+1 and N >= 11: next after 17 is 23, not 29. In range:");
    var res = Enumerable.Range(3, Math.Max(0, sn - 2)).Where(Seed.IsResonant);
    Console.WriteLine($"  resonant N = [{string.Join(", ", res)}]; N = 8 is the smallest EVEN N with resonances (Z3 = 2), seedless.");
    return;
}

// ---- run mode "doubleslit": the phenomenon composed from the atoms (the access layer) ----
// DoubleSlit IS Field at N=1: two places |L>,|R>, the humps = the immortal diagonal, the fringe = the
// between |L><R| (the k=1 coherence) paying -2gamma. Nothing new is computed; the atoms are assembled
// under the name so the phenomenon is recognizable where Pair + Field alone were not. The MEANING lives
// in docs/quantum/DOUBLE_SLIT_TRANSLATED.md.
if (args.Length > 0 && args[0] == "doubleslit")
{
    const double dsg = 0.05, dsdt = 0.05;   // gamma = 0.05: the canonical hardware-anchored rate
    var dsworld = new World();
    var ds = new DoubleSlit(dsworld, dsg);
    double fringe0 = ds.Fringe;
    Console.WriteLine("the double slit, composed from the atoms: Field at N=1 (Pair + Field, nothing new)");
    Console.WriteLine("  two places |L>,|R>: the humps are the immortal diagonal; the fringe is the between |L><R|,");
    Console.WriteLine($"  the k=1 coherence, rate {ds.BetweenRate:0.00} = -2*gamma -- the generator of |rho_LR(t)| = e^(-2*gamma*t).");
    Console.WriteLine("  the watching never touches the humps; it carves the between away. Meaning: DOUBLE_SLIT_TRANSLATED.md");
    Console.WriteLine();
    Console.WriteLine($"  gamma = {dsg} (canonical, hardware-anchored); the coherence 1/e time is 1/(2*gamma) = 10");
    Console.WriteLine($"  {"t",5} {"humps",8} {"fringe",8} {"fringe/f0",10} {"e^(-2gt)",9}");
    for (int s = 0; s <= 600; s++)
    {
        if (s % 120 == 0)
            Console.WriteLine($"  {ds.T,5:0.00} {ds.Humps,8:0.000} {ds.Fringe,8:0.000} {ds.Fringe / fringe0,10:0.000} {Math.Exp(ds.BetweenRate * ds.T),9:0.000}");
        ds.Watch(dsdt);
    }
    Console.WriteLine();
    Console.WriteLine("  the humps stay flat (the particle face, never watched away); the fringe fades on the exp");
    Console.WriteLine("  law (the wave face, the between paying -2gamma). The pattern is what not being watched looks like.");
    return;
}

// ---- run mode "cat": Schrodinger's cat composed as Field at N (the k=N sighting) ----
// The two definite branches (dead |0..0>, alive |1..1>) are the immortal diagonal; the coherence between
// them (k=N) pays -2*gamma*N and dies N times faster than the slit's k=1 -- the same law, the far end.
// The dynamics needs Field's 4^N, so keep N small; the N-scaling law (why a real cat never superposes) is
// pure arithmetic. Meaning in docs/quantum/SCHRODINGERS_CAT_TRANSLATED.md.
if (args.Length > 0 && args[0] == "cat")
{
    int catN = args.Length > 1 ? int.Parse(args[1]) : 4;
    const double catg = 0.05, catdt = 0.05;
    var catworld = new World();
    if (catN <= 12)                                        // the dynamics builds Field's 4^N; keep it small
    {
        var cat = new Cat(catworld, catN, catg);
        double coh0 = cat.CatCoherence;
        double tau = 1.0 / (2.0 * catg * catN);
        Console.WriteLine($"Schrodinger's cat, composed as Field at N={catN} (the k=N sighting, nothing new)");
        Console.WriteLine("  two definite branches |0..0> (dead), |1..1> (alive): the immortal diagonal, the poles that stay;");
        Console.WriteLine($"  the coherence |0..0><1..1> between them is the 'both at once', disagreement k={catN}, rate {cat.CoherenceRate:0.00}");
        Console.WriteLine($"  = -2*gamma*N. Its 1/e time is 1/(2*gamma*N) = {tau:0.0} -- it dies {catN}x faster than the slit's k=1.");
        Console.WriteLine("  Meaning: SCHRODINGERS_CAT_TRANSLATED.md");
        Console.WriteLine();
        Console.WriteLine($"  gamma = {catg}");
        Console.WriteLine($"  {"t",6} {"branches",9} {"cat-coh",8} {"coh/c0",8} {"e^(-2Ngt)",10}");
        int steps = (int)Math.Round(5 * tau / catdt);      // run to ~5 lifetimes
        int every = Math.Max(1, steps / 5);
        for (int s = 0; s <= steps; s++)
        {
            if (s % every == 0)
                Console.WriteLine($"  {cat.T,6:0.00} {cat.Branches,9:0.000} {cat.CatCoherence,8:0.000} {cat.CatCoherence / coh0,8:0.000} {Math.Exp(cat.CoherenceRate * cat.T),10:0.000}");
            cat.Watch(catdt);
        }
        Console.WriteLine();
        Console.WriteLine("  the branches stay flat (the definite poles are immortal); the 'both at once' vanishes fast.");
    }
    else
    {
        Console.WriteLine($"Schrodinger's cat at N={catN}: the dynamics builds Field's 4^N (too big past N=12), so");
        Console.WriteLine("  the run shows the N-scaling law only -- it is pure arithmetic and holds for any N (a real cat's own N).");
    }
    Console.WriteLine();
    Console.WriteLine("  the deeper reading -- why we never see a macroscopic cat both dead and alive (rate -2*N*gamma):");
    Console.WriteLine($"  {"N",4} {"rate -2Ng",10} {"1/e time",10}");
    foreach (int nn in new[] { 1, 2, 4, 8, 16, 32, 64 })
        Console.WriteLine($"  {nn,4} {-2.0 * catg * nn,10:0.00} {1.0 / (2.0 * catg * nn),10:0.0000}");
    Console.WriteLine("  the bigger the superposition, the deeper the rate and the shorter the life -- k=1 the slit,");
    Console.WriteLine("  k=N the cat, one law. A real cat is N ~ 10^23: the 'both at once' is gone before it begins.");
    return;
}

// ---- run mode "group": the mirror group and its antilinear double ----
// F118 + F119 (adopted 2026-07-04): the palindromizer factors, Pi_Z = R o D, and the two generators
// close into the dihedral D4 -- eight signed permutations of the Pauli basis, compared exactly. The
// palindrome splits along the generators (D carries the Hamiltonian sign, R carries the constant),
// the polarity cube's three axes are characters, and the antilinear triangle (theta, conj, dagger)
// docks on as the double D4 x Z2 with the transport law as its engine.
if (args.Length > 0 && args[0] == "group")
{
    int gn = args.Length > 1 ? int.Parse(args[1]) : 3;
    var gworld = new World();
    var group = new MirrorGroup(gworld, gn);
    var gammas = Enumerable.Range(0, gn).Select(l => 0.3 + 0.2 * l).ToArray();

    static string Ph(System.Numerics.Complex p)
        => p == System.Numerics.Complex.One ? "" : p == -System.Numerics.Complex.One ? "-"
         : p == System.Numerics.Complex.ImaginaryOne ? "i" : "-i";
    static string Rule(MirrorGroup.Member m)
        => string.Join("  ", "IXYZ".Select(c => { var (l, p) = m.ApplySite(c); return $"{c}->{Ph(p)}{l}"; }));

    Console.WriteLine("the mirror group and its antilinear double (F118 + F119, adopted 2026-07-04)");
    Console.WriteLine($"  sources docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md + PROOF_ANTILINEAR_TRIANGLE.md; N={gn}");
    Console.WriteLine("  the palindromizer factors: Pi_Z = R o D (D the transpose, R the reflection rho -> rho*F).");
    Console.WriteLine("  eight signed permutations of the Pauli basis -- the dihedral D4, phases in {1,-1,i,-i}:");
    var eight = new (MirrorGroup.Member M, string Form)[]
    {
        (MirrorGroup.Identity, "rho"), (MirrorGroup.PiZ, "rho^T * F"), (MirrorGroup.F, "F rho F"),
        (MirrorGroup.PiY, "F * rho^T"), (MirrorGroup.D, "rho^T"), (MirrorGroup.FD, "F rho^T F"),
        (MirrorGroup.R, "rho * F"), (MirrorGroup.FR, "F rho"),
    };
    foreach (var (m, form) in eight)
        Console.WriteLine($"    {m.Name,-4} {form,-10}  {Rule(m)}");

    Console.WriteLine($"  closure |<R,D>| = {MirrorGroup.Closure(MirrorGroup.R, MirrorGroup.D).Count}; " +
        $"Pi_Z = R o D: {MirrorGroup.Compose(MirrorGroup.R, MirrorGroup.D).Equals(MirrorGroup.PiZ)}; " +
        $"D Pi_Z D = Pi_Y (s r s = r^-1): {MirrorGroup.Compose(MirrorGroup.D, MirrorGroup.Compose(MirrorGroup.PiZ, MirrorGroup.D)).Equals(MirrorGroup.PiY)}");
    Console.WriteLine($"  all eight forms vs their operators, every string at N={gn}: worst residual {group.EightFormsWorstResidual():E1}");
    var (right, wrong) = group.PiZOnRho();
    Console.WriteLine($"  Pi_Z on a full rho: rho^T * F holds at {right:E1}; the wrong-sided F * rho^T is rejected at {wrong:0.0}");

    var (dH, dDiss, rH, rDiss) = group.PalindromeSplitResiduals(j: 1.0, delta: 0.7, gammas: gammas);
    Console.WriteLine();
    Console.WriteLine("  the palindrome splits along the generators (XXZ delta=0.7, site-dependent gamma, full basis):");
    Console.WriteLine($"    D  flips L_H (D L_H D = -L_H):              {dH:E1}      D  fixes the dissipator: {dDiss:E1}");
    Console.WriteLine($"    R  fixes L_H:                               {rH:E1}      R  reflects it, carrying the constant");
    Console.WriteLine($"       (R L_diss R = -L_diss - 2*sigma):        {rDiss:E1}   -- the same constant Mirror pays as its price.");

    Console.WriteLine();
    Console.WriteLine("  the cube of characters: bit_a = char(Ad_Z^N), bit_b = char(Ad_X^N = F), y_par = char(D);");
    int truly = PauliMode.Enumerate(gworld, gn, 0.5).Count(m =>
        MirrorGroup.D.Apply(m.Letters).Phase == System.Numerics.Complex.One
        && MirrorGroup.FD.Apply(m.Letters).Phase == System.Numerics.Complex.One);
    Console.WriteLine($"  the truly cell = the joint-fixed cell of the diagonal mirror pair (D, FD): {truly} of {1 << (2 * gn)} strings at N={gn}.");

    var triangle = new AntilinearTriangle(gworld, gn);
    Console.WriteLine();
    Console.WriteLine("  the antilinear triangle (theta transpose, conj bar, dagger adjoint; dagger = theta o conj):");
    Console.WriteLine("    graded by ell (linearity) and m (multiplicativity): theta (+,-), conj (-,+), dagger (-,-);");
    Console.WriteLine($"    the transport law mu o L_H o mu = ell*m * L_(mu H), non-Hermitian H: worst {triangle.TransportWorstResidual():E1}");
    Console.WriteLine($"    the dephasing dissipator is fixed by all three: worst {triangle.DissipatorFixedWorstResidual(gammas):E1}");
    var (herm, nonHerm) = triangle.FixedPointCollapse();
    Console.WriteLine($"    fixed-point collapse H^T = H-bar iff H = H-dagger: Hermitian {herm:E1}, non-Hermitian split {nonHerm:0.0}");
    var doubled = MirrorGroup.Closure(MirrorGroup.R, MirrorGroup.D, MirrorGroup.K);
    Console.WriteLine($"    the double <R, D, K>: order {doubled.Count}, antilinear members {doubled.Count(m => m.Antilinear)} -- D4 x Z2.");
    Console.WriteLine("  (deliberately outside, named open in F118: the letter group S3; adjoining it would assemble S3 x| D4.)");
    return;
}

// ---- run mode "klein": the parameter-side Klein V4 ----
// F91 + F92 + F93 (adopted 2026-07-04): on each parameter axis (gamma per site, J per bond, h per
// site) the palindromic mirror F71 and the anti-palindromic reshuffle R90 generate a Klein V4, and
// the proofs' sharper law is entry-wise: the F71-refined DIAGONAL blocks of L depend only on the
// pair-sums, so the whole anti-palindromic orbit shares one set of diagonal blocks, cell for cell;
// the breaking lives in the cross-blocks (the eigenvectors) only. No eigensolver anywhere.
if (args.Length > 0 && args[0] == "klein")
{
    int kn = args.Length > 1 ? int.Parse(args[1]) : 6;
    var kworld = new World();
    var klein = new ParameterKlein(kworld, kn);

    static double[] Flat(int len, double v) => Enumerable.Repeat(v, len).ToArray();
    static double[] AntiProfile(int len, double avg, double amp)      // pair-sums 2*avg by construction
    {
        var x = Flat(len, avg);
        for (int l = 0; l < len / 2; l++)
        {
            double d = amp * (l % 2 == 0 ? 1.0 : -0.6) * (l + 1);
            x[l] = avg - d;
            x[len - 1 - l] = avg + d;
        }
        return x;
    }

    Console.WriteLine("the parameter-side Klein V4 (F91 + F92 + F93, adopted 2026-07-04)");
    Console.WriteLine($"  sources docs/proofs/PROOF_F91_GAMMA_NINETY_DEGREES.md + the J/h twins; N={kn}");
    Console.WriteLine("  two commuting involutions on a parameter axis: F71 (reverse) and R90 (x -> 2 avg - reverse);");
    Console.WriteLine("  the anti-palindromic class (every pair-sum = 2 avg) is exactly the FIXED-POINT set of R90.");
    var gLadder = AntiProfile(kn, 0.45, 0.05);
    Console.WriteLine($"  witness gamma profile: [{string.Join(", ", gLadder.Select(v => v.ToString("0.000")))}]  (pair-sums all 0.900)");
    Console.WriteLine($"  R90 fixes it pointwise: {gLadder.Zip(ParameterKlein.Reshuffle(gLadder)).All(p => Math.Abs(p.First - p.Second) < 1e-12)}; the mirror does not: {gLadder.SequenceEqual(ParameterKlein.Mirror(gLadder))}");

    var jFlat = Flat(kn - 1, 1.0);
    var hZero = Flat(kn, 0.0);
    var gFlat = Flat(kn, 0.45);
    Console.WriteLine();
    Console.WriteLine("  the entry-wise law (worst refined-diagonal-block cell difference, ALL (p,q) blocks):");
    Console.WriteLine($"    F91 gamma: uniform vs anti-palindromic  {klein.DiagonalBlocksResidual(gFlat, jFlat, hZero, gLadder, jFlat, hZero):E1}   (identical, cell for cell)");
    var gBroken = (double[])gLadder.Clone();
    (gBroken[0], gBroken[1]) = (gBroken[1], gBroken[0]);               // same sum, wrong pair-sums
    Console.WriteLine($"    F91 gamma: uniform vs pair-sum-broken   {klein.DiagonalBlocksResidual(gFlat, jFlat, hZero, gBroken, jFlat, hZero):0.000}   (rejected: same total, wrong pair-sums)");
    var jAnti = AntiProfile(kn - 1, 1.0, 0.08);
    Console.WriteLine($"    F92 J:     uniform vs anti-palindromic  {klein.DiagonalBlocksResidual(gFlat, jFlat, hZero, gFlat, jAnti, hZero):E1}");
    var hAnti = AntiProfile(kn, 0.3, 0.07);
    Console.WriteLine($"    F93 h:     uniform vs anti-palindromic  {klein.DiagonalBlocksResidual(gFlat, jFlat, Flat(kn, 0.3), gFlat, jFlat, hAnti):E1}");
    Console.WriteLine();
    Console.WriteLine("  where the breaking went (the worst F71 cross-block cell):");
    Console.WriteLine($"    uniform gamma:          {klein.CrossBlockNorm(gFlat, jFlat, hZero):E1}   (F71 exact: no cross-block at all)");
    Console.WriteLine($"    anti-palindromic gamma: {klein.CrossBlockNorm(gLadder, jFlat, hZero):0.000}   (F71 broken -- but only the eigenvectors carry it)");
    Console.WriteLine("  the spectral invariance (the registry's F91-F93 statement) is the corollary: identical blocks,");
    Console.WriteLine("  identical spectra. the whole orbit shares the uniform world's decay rates; the asymmetry is a way of seeing.");
    return;
}

// ---- run mode "hardness": the hardness of the palindrome ----
// The F87 bloc (adopted 2026-07-04): the spectral trichotomy truly/soft/hard, read WITHOUT a
// spectrum -- the GF(2)[x] valuation (F115), the letter-parity purity rules (F107/F109), the cell
// rules (F110/F111), and the trace face of the all-gamma converse (F117). Two independent
// eigensolver-free certificates, shown agreeing on the K3 trio.
if (args.Length > 0 && args[0] == "hardness")
{
    var hworld = new World();
    var hardness = new Hardness(hworld);

    Console.WriteLine("the hardness of the palindrome (the F87 bloc, adopted 2026-07-04)");
    Console.WriteLine("  sources PROOF_F103_F87_Z2_CUBED_REFINEMENT.md sections 6-7 (+ 7.7 = F115), PROOF_F87_WINDOWED_");
    Console.WriteLine("  MONOMIAL_CONVERSE.md section 5 (F117), PROOF_F107/F109; the spectral classifier itself stays outside.");
    Console.WriteLine();
    Console.WriteLine("  the valuation face (F115): a diagonal-cell pair is hard iff its X/Y window masks have");
    Console.WriteLine("  DIFFERENT (1+x)-adic valuations -- the whole verdict in one subtraction. the K3 trio:");
    var trio = new (string Name, ulong Mask)[] { ("XXZ", 0b011), ("XZX", 0b101), ("ZXX", 0b110) };
    foreach (var (name, mask) in trio)
        Console.WriteLine($"    {name}: mask {Convert.ToString((long)mask, 2).PadLeft(3, '0')} -> v = {Hardness.Valuation(mask)}");
    foreach (var (i, j2) in new[] { (0, 2), (0, 1), (1, 2) })
        Console.WriteLine($"    ({trio[i].Name}, {trio[j2].Name}): {(Hardness.IsHardPair(trio[i].Mask, trio[j2].Mask) ? "HARD" : "soft")}");
    Console.WriteLine($"  hard mask-pairs (A203241) k=4,5,6: {Hardness.HardMaskPairCount(4)}, {Hardness.HardMaskPairCount(5)}, {Hardness.HardMaskPairCount(6)}; dressed: {Hardness.DressedHardCount(4)}, {Hardness.DressedHardCount(5)}, {Hardness.DressedHardCount(6)}");
    Console.WriteLine($"  obstruction ceiling min(2W-1, 2k-3): k=3,W=2 -> {Hardness.ObstructionCeiling(3, 2)} (the always-triangle)");
    Console.WriteLine();
    Console.WriteLine("  the cube face: truly forces y_par=0 under every dephase letter (F107); the mother sector's");
    Console.WriteLine("  non-truly side is all-odd, y_par=1 (F109); hard lives ONLY in the dephase letter's own Klein");
    Console.WriteLine($"  cell -- Z->{Hardness.DiagonalCell('Z')}, X->{Hardness.DiagonalCell('X')}, Y->{Hardness.DiagonalCell('Y')} -- with the Y-inversion (F110):");
    Console.WriteLine($"  k=3 split {Hardness.HardSplitK3} (N-stable, F103 = F105 bit-exact); k=N=4 fully pure 228:0 by the");
    Console.WriteLine($"  pure-D template rule (F111): 528 pairs = {Hardness.TemplateDecompositionK4} (pure-pure + pure-mixed hard, mixed-mixed soft).");
    Console.WriteLine();
    Console.WriteLine("  the trace face (F117): odd power-sums of M = A + gamma*Q -- traces, never eigenvalues.");
    double p3 = hardness.OddPowerSums(new[] { "ZII".ToCharArray() }, n: 3, j: 1.0, gamma: 0.5, upToOdd: 3)[1];
    Console.WriteLine($"    the cell-free m=3 face, H = Z_0 at N=3: p_3 = {p3:0.0} = 6*4^N*gamma = {6.0 * 64 * 0.5:0.0} (any single-site Z breaks the palindrome at every gamma)");
    var soft = hardness.OddPowerSums(new[] { "XXZ".ToCharArray(), "ZXX".ToCharArray() }, n: 4, j: 1.0, gamma: 0.5, upToOdd: 9);
    var hard = hardness.OddPowerSums(new[] { "XXZ".ToCharArray(), "XZX".ToCharArray() }, n: 4, j: 1.0, gamma: 0.5, upToOdd: 9);
    Console.WriteLine($"    soft (XXZ,ZXX) at N=4: p_1..p_9 = [{string.Join(", ", soft.Select(v => v.ToString("0.0###")))}] -- all vanish (spec symmetric)");
    Console.WriteLine($"    hard (XXZ,XZX) at N=4: p_1..p_9 = [{string.Join(", ", hard.Select(v => v.ToString("0.0###")))}]");
    Console.WriteLine($"    the deg-1 ladder is SILENT (p_7 = 0: silence is not softness); m* = 2*girth+3 = 9 fires through the");
    Console.WriteLine($"    d=3 class: p_9 = {hard[4]:0} = 2064384 * gamma^3 (the F117 CRT integer, a single positive monomial:");
    Console.WriteLine("    hard at EVERY gamma). the two certificates -- GF(2) valuation and integer traces -- agree pair for pair.");
    return;
}

// ---- run mode "router": the golden ceiling router ----
// F116 (adopted 2026-07-04): the two Z-middle ceiling cases are palindromized LOCALLY by the
// period-4 [a,a,b,b] router on the golden locus; the whole family is metallic, r(c) = the metallic
// mean. Shown here: the window lemma at all four offsets, the locus gating the frame, and the
// dense two-sided end-to-end (previously Python-only) at N=5 with site-dependent rates.
if (args.Length > 0 && args[0] == "router")
{
    var rworld = new World();
    var router = new Router(rworld);
    var phi = Formulas.F116_MetallicMean(1.0);

    Console.WriteLine("the golden ceiling router (F116, adopted 2026-07-04)");
    Console.WriteLine("  source docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md sections 1-3 + 8; the ceiling cases");
    Console.WriteLine("  XZX+XZY+YZX (and the X<->Y sibling) are soft, and LOCALLY so: W = (x)_l q_{l mod 4},");
    Console.WriteLine($"  frame [a,a,b,b] with a = phi X + Y, b = X - phi Y on the golden locus; q^2 = -(2+phi) I = {-(2 + phi):0.000} I.");
    Console.WriteLine();
    Console.WriteLine("  the window lemma ({Q_k, S} = 0 at all four offsets, the 64-dim window space):");
    Console.WriteLine($"    golden (c=1):  worst offset norm = {Router.WindowAnticommutatorNorm(1.0, phi, sibling: false):E1}");
    Console.WriteLine($"    the sibling:   worst offset norm = {Router.WindowAnticommutatorNorm(1.0, phi, sibling: true):E1}");
    Console.WriteLine("  the metallic family on the soft line t2 = t3, r(c) = (c + sqrt(c^2+4))/2 (F116_MetallicMean):");
    foreach (double c in new[] { 0.0, 1.0, 2.0, 3.0 })
    {
        double rc = Formulas.F116_MetallicMean(c);
        string name = c == 1.0 ? "golden phi" : c == 2.0 ? "silver 1+sqrt2" : c == 3.0 ? "bronze" : "the 45-degree frame";
        Console.WriteLine($"    c={c:0.0}: r = {rc:0.000000} ({name}); on the locus {Router.WindowAnticommutatorNorm(c, rc, false):E1}, off it (r+0.25) {Router.WindowAnticommutatorNorm(c, rc + 0.25, false):0.00}");
    }
    Console.WriteLine();
    Console.WriteLine("  the two-sided dense end-to-end (the proof's section 2; previously Python-only):");
    var gammas = Enumerable.Range(0, 5).Select(l => 0.3 + 0.15 * l).ToArray();
    var (chiralP, chiralQ, conjugation) = router.DenseResiduals(n: 5, c: 1.0, gammas);
    Console.WriteLine($"    N=5, site-dependent gammas: P H P^-1 = -H at {chiralP:E1}, Q H Q^-1 = -H at {chiralQ:E1},");
    Console.WriteLine($"    W L W^-1 = -L - 2 sigma on the full Pauli basis at {conjugation:E1}.");
    Console.WriteLine("  the palindrome is realized by a LOCAL period-4 product -- no non-local mirror needed: the");
    Console.WriteLine("  ceiling arc's 6 -> 4 -> 2 -> 0 closes, and the soft line's whole frame family is metallic.");
    return;
}

// ---- run mode "anti": the rules turned around ----
// The mirror's rho-level face: run the SAME living world twice, once watching disagreement (the rule,
// rate -2*gamma*k) and once watching AGREEMENT (the rule turned around, rate -2*gamma*(N-k)). The
// mirror guarantees the anti-world is the old world read through the bra complement, entry for entry:
// nothing was taken, the rule only wears its complement. The conservation law moves with it: the
// normal world keeps its trace (the diagonal is immortal), the anti-world keeps its ANTI-trace (the
// anti-diagonal, k = N, the GHZ-like maximal disagreement) while its trace dies.
if (args.Length > 0 && args[0] == "anti")
{
    int an = args.Length > 1 ? int.Parse(args[1]) : 3;
    const double aj = 1.0, ag = 0.5, adt = 0.05;
    const int aticks = 40;
    var aworld = new World();
    int s = 1, sbar = (1 << an) - 1 - s;                 // |0..01> and its complement |1..10>

    var normal = new Restless(aworld, an, aj, ag);
    normal.Seed(s, 0.5); normal.Seed(sbar, 0.5);         // structure: two populations, trace 1

    var anti = new Restless(aworld, an, aj, ag, antiWatching: true);
    anti.SeedCoherence(s, sbar, 0.5);                    // the SAME state read through X^N: |s><sbar| + twin

    Console.WriteLine("the rules turned around (the mirror's rho-level face)");
    Console.WriteLine($"  N={an}, J={aj}, gamma={ag}, dt={adt}; normal rule: rate -2*gamma*k (disagreement watched);");
    Console.WriteLine($"  turned rule: rate -2*gamma*(N-k) (AGREEMENT watched). the mirror says: anti(t) = normal(t) * X^N, exactly.");
    Console.WriteLine($"  {"t",5} {"trace(normal)",14} {"antitrace(normal)",18} {"trace(anti)",12} {"antitrace(anti)",16}");
    for (int tick = 0; tick <= aticks; tick++)
    {
        if (tick % 5 == 0)
            Console.WriteLine($"  {normal.T,5:0.00} {normal.Structure,14:0.000000} {normal.AntiStructure,18:0.000000} {anti.Structure,12:0.000000} {anti.AntiStructure,16:0.000000}");
        if (tick == aticks) break;
        normal.Step(adt); anti.Step(adt);
    }

    // the exact read-through: anti(i,j) must equal normal(i, complement j), entry for entry.
    double worstMatch = 0;
    int adim = 1 << an;
    for (int i = 0; i < adim; i++)
        for (int j = 0; j < adim; j++)
            worstMatch = Math.Max(worstMatch, (anti[i, j] - normal[i, adim - 1 - j]).Magnitude);
    Console.WriteLine($"  read-through check at t={normal.T:0.00}: worst |anti(i,j) - normal(i,~j)| = {worstMatch:E1}");

    var hn = normal.WeightByDisagreement();
    var ha = anti.WeightByDisagreement();
    Console.WriteLine($"  k-histogram normal: [{string.Join(", ", hn.Select(v => v.ToString("0.000")))}]");
    Console.WriteLine($"  k-histogram anti:   [{string.Join(", ", ha.Select(v => v.ToString("0.000")))}]  (the normal one, read backward)");
    Console.WriteLine("  the trace died in the anti-world; the anti-trace lives there instead. the law was not taken, it moved");
    Console.WriteLine("  to the anti-diagonal -- the rule turned around is the same rule wearing its complement.");
    return;
}

// ---- run mode "scale": the complexity wall and the cut ----
// Why the full-spectrum side hit the wall and a state's dynamics did not. full = 4^N (the whole Liouvillian
// an eigendecomposition tackles); single-exc = the (1,1) block a one-excitation state lives in (N^2,
// polynomial); half-fill = the largest block C(N,N/2)^2 (the spectrum's irreducible core, still exponential).
if (args.Length > 0 && args[0] == "scale")
{
    var sworld = new World();
    Console.WriteLine("the complexity wall and the cut (why the full spectrum hit it, a state's dynamics did not)");
    Console.WriteLine($"  {"N",3} {"full 4^N",13} {"single-exc",11} {"forbidden",10} {"half-fill",11} {"hf forbidden",13}");
    for (int n = 2; n <= 12; n++)
    {
        long full = 1L << (2 * n);
        long se = new Block(sworld, n, 1, 1).Size;          // N^2, the one-excitation block
        long hf = new Block(sworld, n, n / 2, n / 2).Size;  // C(N,N/2)^2, the largest block
        Console.WriteLine($"  {n,3} {full,13} {se,11} {100.0 * (full - se) / full,9:0.00}% {hf,11} {100.0 * (full - hf) / full,11:0.000}%");
    }
    Console.WriteLine("  single-exc: N^2 -- a one-excitation state runs in a polynomial corner, tractable at any N.");
    Console.WriteLine("  half-fill: ~4^N/sqrt(N) -- the full spectrum needs it, so the spectrum stays exponential.");
    Console.WriteLine("  the wall was global (all 4^N, the spectrum); a single state's dynamics is block-local.");
    return;
}

// ---- run mode "topo": play with the geometry (who can reach whom) ----
// Same handshake and seed, different bonds. All excitation-conserving, so all share the block cut (F63);
// only the dynamics -- how fast novelty is born and spreads -- differs.
if (args.Length > 0 && args[0] == "topo")
{
    int gn = args.Length > 1 ? int.Parse(args[1]) : 4;
    const double gj = 1.0, gg = 0.5, dt = 0.05;
    var tworld = new World();
    Console.WriteLine($"playing with the geometry (who can reach whom): N={gn}, seed |0...01>, J={gj}, gamma={gg}");
    Console.WriteLine("  same handshake, different bonds. all excitation-conserving -> the same block cut (F63);");
    Console.WriteLine("  only the dynamics differs. (seed site 0 is the chain's endpoint, but the star's hub.)");
    Console.WriteLine($"  {"topology",9} {"bonds",6} {"alive",6} {"forbidden",10}   novelty at t = 0.25 / 0.50 / 1.00");
    foreach (var name in new[] { "chain", "ring", "star", "complete" })
    {
        var rest = new Restless(tworld, gn, gj, gg, Topology.Named(name, gn));
        rest.Seed(1);
        var snap = new double[3];
        int[] marks = { 5, 10, 20 };
        int mi = 0;
        for (int tick = 0; tick <= 20; tick++)
        {
            if (mi < marks.Length && tick == marks[mi]) snap[mi++] = rest.Novelty;
            rest.Step(dt);
        }
        Console.WriteLine($"  {name,9} {Topology.Named(name, gn).Length,6} {rest.AliveCount,6} {rest.ForbiddenCount,10}   {snap[0],6:0.000} / {snap[1],6:0.000} / {snap[2],6:0.000}");
    }
    Console.WriteLine("  alive is identical across all four -- the cut is topology-invariant. the geometry only");
    Console.WriteLine("  steers the dynamics: more bonds from the seed (star hub, complete) -> novelty born faster.");
    return;
}

// ---- run mode "regime": how alive at each Q = J/gamma (the clock) ----
// The same world at different ratios of restlessness to watching. Low Q: the watching wins (overdamped,
// novelty barely born and quickly drained -- near death). High Q: the restlessness wins (underdamped,
// novelty churns on and on). The clock angle theta = arctan(Q) reads the regime (reused from Clock).
if (args.Length > 0 && args[0] == "regime")
{
    int gn = args.Length > 1 ? int.Parse(args[1]) : 4;
    const double gg = 0.5, dt = 0.05;
    const int ticks = 80;                            // to t = 4
    var rworld = new World();
    Console.WriteLine($"how alive at each regime Q = J/gamma (N={gn} chain, seed |0...01>, gamma={gg})");
    Console.WriteLine("  the clock: theta = arctan(Q). low Q = the watching wins (overdamped, ~death);");
    Console.WriteLine("  high Q = the restlessness wins (underdamped, novelty churns on and on).");
    Console.WriteLine($"  {"Q",5} {"theta",7} {"peak nov",9} {"t@peak",7} {"late churn (t3-4)",17}");
    foreach (double q in new[] { 0.2, 0.5, 1.0, 2.0, 5.0, 20.0 })
    {
        double j = q * gg;
        var clk = new Clock(rworld, j, gg);          // theta = arctan(Q), Q = J/gamma (the adopted clock)
        var rest = new Restless(rworld, gn, j, gg, Topology.Chain(gn));
        rest.Seed(1);
        double peak = 0, tpeak = 0, lateMin = double.MaxValue, lateMax = 0;
        for (int tick = 0; tick <= ticks; tick++)
        {
            double nov = rest.Novelty;
            if (nov > peak) { peak = nov; tpeak = rest.T; }
            if (tick >= ticks - 20) { lateMin = Math.Min(lateMin, nov); lateMax = Math.Max(lateMax, nov); }
            rest.Step(dt);
        }
        Console.WriteLine($"  {clk.Q,5:0.0} {clk.ThetaDeg,6:0.0} {peak,9:0.000} {tpeak,7:0.00} {lateMax - lateMin,17:0.000}");
    }
    Console.WriteLine("  peak rises and t@peak falls, both monotonic -- more restlessness births more novelty, sooner.");
    Console.WriteLine("  late churn is ~0 below Q=1 (overdamped: a slow bump that just decays) and switches on above it");
    Console.WriteLine("  (underdamped: novelty oscillates). It is not monotonic at high Q -- the chain is multi-mode, so");
    Console.WriteLine("  the underdamped churn beats; a single scalar cannot hold it. The clean turn is Q=1, theta=45deg.");
    return;
}

// ---- run mode "seeds": play with the seed (how many excitations) ----
// The seed picks the joint-popcount block it lives in. Empty (p=0) and full (p=N) are frozen -- no
// excitation can hop, no novelty. Half-filling is the richest: the biggest block, the most life. Block
// size and peak novelty both peak at the middle and mirror p <-> N-p (the palindrome, in the seed).
if (args.Length > 0 && args[0] == "seeds")
{
    int gn = args.Length > 1 ? int.Parse(args[1]) : 4;
    const double gj = 1.0, gg = 0.5, dt = 0.05;
    var sworld = new World();
    Console.WriteLine($"playing with the seed (how many excitations): N={gn} chain, J={gj}, gamma={gg}");
    Console.WriteLine("  the seed picks the block (p,p) it lives in. empty (p=0) and full (p=N) are frozen --");
    Console.WriteLine("  no excitation can hop, no novelty. half-filling is the richest (biggest block, most life).");
    Console.WriteLine($"  {"p (exc)",8} {"block C(N,p)^2",14} {"peak novelty",13}");
    for (int p = 0; p <= gn; p++)
    {
        int seed = (1 << p) - 1;                     // |0...0 1^p>, p excitations at one end
        var rest = new Restless(sworld, gn, gj, gg, Topology.Chain(gn));
        rest.Seed(seed);
        double peak = 0;
        for (int tick = 0; tick <= 40; tick++) { peak = Math.Max(peak, rest.Novelty); rest.Step(dt); }
        Console.WriteLine($"  {p,8} {rest.AliveCount,14} {peak,13:0.000}");
    }
    Console.WriteLine("  empty and full are frozen points (one cell, no life); the middle is where the world lives");
    Console.WriteLine("  most. block size and peak novelty both mirror p <-> N-p -- the palindrome, now in the seed.");
    return;
}

// ---- run mode "cone": the light-cone of a single excitation (the memory cut pays off) ----
// What was blocked before: large-N dynamics. The spectrum died at 4^N (N=8). Here the single excitation is
// stored on its block (N x N, not 4^N), so we watch it at an N the spectrum could never hold -- and a
// dynamical, spatial picture the static spectrum cannot show: the ballistic light-cone, decohered by the watching.
if (args.Length > 0 && args[0] == "cone")
{
    int gn = args.Length > 1 ? int.Parse(args[1]) : 60;
    const double gj = 1.0, gg = 0.05, dt = 0.1;        // canonical gamma=0.05: ballistic, then decohering
    const int rows = 28, stepsPerRow = 4;
    var cworld = new World();
    var cone = new Cone(cworld, gn, gj, gg);
    cone.Seed(gn / 2);                                  // one excitation in the middle
    Console.WriteLine($"the light-cone of a single excitation (memory cut: N={gn} sites, rho is {gn}x{gn}, not 4^{gn})");
    Console.WriteLine($"  seed in the middle, J={gj}, gamma={gg}. ballistic spreading (v~2J), the watching decoheres the front.");
    Console.WriteLine("  each row a time, each column a site, density = excitation population (per-row scaled).");
    const string ramp = " .:-=+o*#";
    for (int row = 0; row < rows; row++)
    {
        double max = 1e-12;
        for (int a = 0; a < gn; a++) max = Math.Max(max, cone.Population(a));
        var sb = new System.Text.StringBuilder();
        for (int a = 0; a < gn; a++)
        {
            int lvl = (int)Math.Round(cone.Population(a) / max * (ramp.Length - 1));
            sb.Append(ramp[Math.Clamp(lvl, 0, ramp.Length - 1)]);
        }
        Console.WriteLine($"  t={cone.T,4:0.0} |{sb}|");
        for (int s = 0; s < stepsPerRow; s++) cone.Step(dt);
    }
    Console.WriteLine($"  the cone opens at ~2J. blocked before (4^{gn} is astronomical); trivial now ({gn}x{gn} = {gn * gn} cells).");
    return;
}

// ---- run mode "spread": the ballistic-to-diffusive crossover (mining the unblocked vista) ----
// The front spread sigma(t) = sqrt(sum_a (a-a0)^2 P(a)). Ballistic (sigma ~ v t, v~2J) until the coherence
// that carries it dies (~1/4gamma), then diffusive (sigma ~ sqrt(t)). The crossover length is L ~ J/2gamma =
// Q/2. Known transport physics (ENAQT); here it is a validation of the cut, measured at N the spectrum cannot reach.
if (args.Length > 0 && args[0] == "spread")
{
    int gn = args.Length > 1 ? int.Parse(args[1]) : 80;
    const double gj = 1.0, dt = 0.05;
    int[] marks = { 20, 40, 80, 160, 320 };            // t = 1, 2, 4, 8, 16
    var sworld = new World();
    Console.WriteLine($"the front spread sigma(t): ballistic (~t) until the coherence dies (~1/4gamma), then diffusive (~sqrt t)");
    Console.WriteLine($"  N={gn} chain, seed middle, J={gj}. doubling t: sigma x2 is ballistic, sigma x1.41 is diffusive.");
    Console.WriteLine($"  {"gamma",6} {"Q",5} {"L=Q/2",6}   sigma at t = 1 / 2 / 4 / 8 / 16");
    foreach (double g in new[] { 0.0, 0.05, 0.2 })
    {
        var cone = new Cone(sworld, gn, gj, g);
        int a0 = gn / 2;
        cone.Seed(a0);
        var sig = new List<double>();
        for (int s = 0; s <= marks[^1]; s++)
        {
            if (marks.Contains(s))
            {
                double m2 = 0;
                for (int a = 0; a < gn; a++) m2 += (a - a0) * (a - a0) * cone.Population(a);
                sig.Add(Math.Sqrt(m2));
            }
            cone.Step(dt);
        }
        string q = g > 0 ? (gj / g).ToString("0.0") : "inf";
        string l = g > 0 ? (gj / (2 * g)).ToString("0.0") : "inf";
        Console.WriteLine($"  {g,6:0.00} {q,5} {l,6}   {string.Join(" / ", sig.Select(x => x.ToString("00.0")))}");
    }
    Console.WriteLine("  gamma=0 (Q=inf): pure ballistic -- sigma keeps doubling as t doubles.");
    Console.WriteLine("  gamma>0: ballistic until sigma ~ L=Q/2, then the doubling collapses to ~1.41 -- diffusive. The");
    Console.WriteLine("  cut reproduces the known crossover, at N the spectrum could never hold. The wall is broken.");
    return;
}

const double gamma = 0.5;
var world = new World();
Console.WriteLine($"empty world (Z-dephasing, no Hamiltonian), gamma={gamma}");
Console.WriteLine($"World.Own (left, the frame): {string.Join(",", world.Own)}");
Console.WriteLine();

foreach (int N in new[] { 2, 3, 4 })
{
    var modes = PauliMode.Enumerate(world, N, gamma).ToList();
    Console.WriteLine($"==== N={N}  ({modes.Count} Pauli-string modes = {1 << N}*sum C(N,k); the superposition basis) ====");

    for (int k = 0; k <= N; k++)
    {
        int count = modes.Count(m => m.K == k);
        string tag = k == 0 ? "  [populations / kernel]"
                   : k == N ? "  [drain]"
                   : k == N - k ? "  [SELF-MIRROR center]" : "";
        Console.WriteLine($"  k={k}: {count,3} modes   rate {-2.0 * gamma * k,5:0.0}   mirror k={N - k}{tag}");
    }

    Console.WriteLine(N % 2 == 0
        ? $"  fold k<->N-k about k=N/2={N / 2}: even N, the self-mirror sector EXISTS (k={N / 2})"
        : $"  fold k<->N-k about k=N/2={(N / 2.0):0.0}: odd N, axis BETWEEN rungs, NO self-mirror sector");

    int kc = (N % 2 == 0) ? N / 2 : (N - 1) / 2;   // central (even) or inner (odd)
    var cells = modes.Where(m => m.K == kc)
                     .GroupBy(m => m.Klein)
                     .OrderBy(g => g.Key.Y).ThenBy(g => g.Key.Z)
                     .Select(g => $"(nY%2={g.Key.Y},nZ%2={g.Key.Z}):{g.Count()}");
    Console.WriteLine($"  superposition / Klein quartering of sector k={kc}: {string.Join("  ", cells)}   [truly = (0,0)]");

    // Grading B (T1): the joint-popcount blocks a number-conserving Hamiltonian respects.
    Console.WriteLine($"  Grading B: {(N + 1) * (N + 1)} joint-popcount blocks (p,q), size C(N,p)*C(N,q):");
    long total = 0;
    for (int p = 0; p <= N; p++)
    {
        var sizes = new List<string>();
        for (int q = 0; q <= N; q++)
        {
            var b = new Block(world, N, p, q);
            total += b.Size;
            sizes.Add($"{b.Size,3}");
        }
        Console.WriteLine($"    p={p}: {string.Join(" ", sizes)}");
    }
    Console.WriteLine($"    sum = {total} = 4^{N}; the diagonal (p,p) carries k=0, the {N + 1} populations (the kernel)");
    Console.WriteLine();
}

// The special cases as T1 closed forms , adopted, not re-derived (all already Tier-1 in the repo).
Console.WriteLine("==== special cases (T1 closed forms, adopted) ====");
for (int n = 2; n <= 6; n++)
    Console.WriteLine($"  N={n}: V-effect 1+cos(pi/N) = {1.0 + Math.Cos(Math.PI / n):0.0000}   ceiling g2(K_N)=4/N = {4.0 / n:0.0000}");
Console.WriteLine($"  N=3: minimal EP, Q_EP = 2/g_eff (the whole rate ladder is HD=1 and HD=3)");
Console.WriteLine($"  N=4: self-fold diabolic q_EP = sqrt((-1+sqrt13)/6) = {Math.Sqrt((-1.0 + Math.Sqrt(13.0)) / 6.0):0.0000}  (real q, N=4 only); (2,2) ceiling K_4 = 2-2/sqrt3 = {2.0 - 2.0 / Math.Sqrt(3.0):0.0000}");
Console.WriteLine($"  N=5: diabolics leave the real axis into complex-conjugate pairs (only loud defective EPs remain real)");

// H on: the superposition leaves the integer grid (the on-grid d_total folds are T1, adopted).
Console.WriteLine();
Console.WriteLine("==== H on: the superposition leaves the grid (T1 d_total folds, adopted) ====");
foreach (int n in new[] { 2, 3, 4 })
{
    var bare = Redistribution.Bare(n);
    var grid = Redistribution.OnGrid(n)!;
    var diff = Enumerable.Range(0, n + 1).Select(k => grid[k] - bare[k]).ToArray();
    int off = (1 << (2 * n)) - grid.Sum();
    Console.WriteLine($"  N={n}: bare [{string.Join(",", bare)}]  ->  on-grid [{string.Join(",", grid)}]   diff [{string.Join(",", diff)}]   off-grid {off}");
}
Console.WriteLine("  edges stay N+1 (kernel + drain); even N spikes the center k=N/2; the rest bleeds off-grid (fractional <n_XY>)");

// The clock: theta = arctan(Q), Q = J/gamma (T1, F95). The two hands of one winding mode.
Console.WriteLine();
Console.WriteLine("==== the clock: theta = arctan(Q), Q = J/gamma (T1) ====");
foreach (double q in new[] { 0.0, 0.5, 1.0, 2.0, 5.0, 100.0 })
{
    var clk = new Clock(world, j: q, gamma: 1.0);   // J=q, gamma=1 => Q=q
    string mark = q == 0.0 ? "  [J=0: pure decay]" : Math.Abs(q - 1.0) < 1e-9 ? "  [theta=45deg <-> 1/4]" : q >= 100.0 ? "  [carbon, deep quantum]" : "";
    Console.WriteLine($"  Q={clk.Q,5:0.0}: theta = {clk.ThetaDeg,5:0.0} deg{mark}");
}
Console.WriteLine("  gamma=0 -> theta=90deg (pure circle, turning forever); 90deg <-> 1/2; the two hands are gamma (radial) and J (angular)");

// The survivor: the slowest non-stationary mode (T1), regime-dependent.
Console.WriteLine();
Console.WriteLine("==== the survivor: the slowest non-stationary mode (T1) ====");
foreach (int n in new[] { 2, 3, 4, 5 })
{
    var s = new Survivor(world, n);
    string lo = s.HasHalfFillingSurvivor
        ? $"low-Q: half-filling k={n / 2}, R-odd/X-odd, dark"
        : "low-Q: no (N/2,N/2) sector (odd N), no half-filling survivor";
    Console.WriteLine($"  N={n}: {lo};  hands over at Q*={s.Qstar:0.00} to the (0,1) band edge (<n_XY>=1, rate -2g)");
}

// R-parity and mod-4 (T1): the even/odd EP behaviour and the N=3-mod-4 amplification.
Console.WriteLine();
Console.WriteLine("==== R-parity & mod-4 (T1) ====");
foreach (int n in new[] { 3, 4, 5, 6 })
{
    bool epSplits = n % 2 == 1;     // odd N: EP splits by R-parity; even N: isospectral
    bool amplifies = n % 4 == 3;    // self-Pi SE mode R-odd at N=3 mod 4 -> amplifies
    Console.WriteLine($"  N={n}: {(epSplits ? "EP SPLITS by R-parity (sigma+ != sigma-)" : "EP isospectral (sigma+ = sigma-, no split)")};  " +
                      $"mod-4: {(amplifies ? "self-Pi SE mode R-odd, amplifies (~2.2x)" : "non-amplifying")}");
}
Console.WriteLine("  even N: PH partners (k,N+1-k) carry opposite R-parity -> isospectral; odd N: same R-parity + a lone R-fixed zero mode -> split");

// F-registry closed forms, adopted verbatim (J=1).
Console.WriteLine();
Console.WriteLine("==== F-registry closed forms (adopted verbatim, J=1) ====");
foreach (int n in new[] { 3, 4, 5 })
{
    Console.WriteLine($"  N={n}: F2  omega_k=4J(1-cos(pi k/N))   = [{string.Join(", ", Formulas.F2_Dispersion(n, 1.0).Select(x => x.ToString("0.00")))}]");
    Console.WriteLine($"        F2b E_k=2J cos(pi k/(N+1))       = [{string.Join(", ", Formulas.F2b_SingleExcitation(n, 1.0).Select(x => x.ToString("0.00")))}]");
}
Console.WriteLine($"  F1 residuals (N=4, gamma=0.5: sg=2, sg2=1): ||M(T1)||^2 = {Formulas.F1_T1Residual(4, 2.0, 1.0):0.0}, ||M(depol)||^2 = {Formulas.F1_DepolResidual(4, 2.0, 1.0):0.0}");
Console.WriteLine($"  F2b-corollary coherence hand omega_mem=2J cos(pi/(N+1)): N=3,4,5 = {Formulas.OmegaMem(3, 1, 0):0.000}, {Formulas.OmegaMem(4, 1, 0):0.000}, {Formulas.OmegaMem(5, 1, 0):0.000}  (sqrt2, phi, sqrt3)");
Console.WriteLine($"  coherence horizon Q*(N) exact: N=2..5 = {Formulas.Qstar(2):0.000}, {Formulas.Qstar(3):0.000}, {Formulas.Qstar(4):0.000}, {Formulas.Qstar(5):0.000}  (asymptotic 2N/pi)");
foreach (int n in new[] { 3, 4, 5 })
{
    var (mn, mx, bw) = Formulas.F3_RateBounds(n, 0.5);
    Console.WriteLine($"  N={n}: F3 min={mn:0.0}/max={mx:0.0}/bw={bw:0.0}  F4 kernel={Formulas.F4_KernelDim(n)}  F5 depol-err={Formulas.F5_DepolError(n, 0.5):0.000}  F23 XOR-frac={Formulas.F23_XorFraction(n):0.0000}  F50 w1-deg={Formulas.F50_Weight1Degeneracy(n)}");
}
Console.WriteLine($"  F33 N=3 exact rates (g=0.5) = [{string.Join(", ", Formulas.F33_N3Rates(0.5).Select(x => x.ToString("0.000")))}]  (<n_XY>=1, 4/3, 5/3)");
foreach (int n in new[] { 3, 4, 5 })
{
    var (qmax, qmin, qmean, qspread) = Formulas.F7_QSpectrum(n, 1.0, 0.5);
    Console.WriteLine($"  N={n}: F7 Q max/min/mean/spread = {qmax:0.00}/{qmin:0.00}/{qmean:0.0}/{qspread:0.00}");
}
Console.WriteLine($"  F8 2x law (N=4,g=0.5): unpaired={Formulas.F8_DecayLaw(4, 0.5).Unpaired:0.0}, paired-mean={Formulas.F8_DecayLaw(4, 0.5).PairedMean:0.0}, ratio=2");
Console.WriteLine($"  F12 single-qubit crossing t*/T2 = {Formulas.F12_CrossingFraction} (root of x^3+x=1/2)");
Console.WriteLine($"  F16 fold R=C(Psi+R)^2, boundary CPsi={Formulas.F16_FoldBoundary} (Mandelbrot u->u^2+c)");
Console.WriteLine($"  F25 CPsi Bell+ Z-deph: crossing f*={Formulas.F25_CrossingF}, K_Z={Formulas.F25_K};  F27 K_X=K_Y=ln2/8={Formulas.F27_KX:0.0000}, K_depol={Formulas.F27_KDepol}");
Console.WriteLine($"  F15 theta compass at CPsi=0.5: {Formulas.F15_ThetaDeg(0.5):0.0} deg (0 at the 1/4 crossing)");
Console.WriteLine($"  F34 qubit necessity: d^2-2d=0 -> d in {{{string.Join(", ", Formulas.F34_QubitNecessity())}}} (d=0 nothing, d=2 the qubit, the polarity root)");
Console.WriteLine($"  D1 BW=8J cos(pi/N) (N=4,J=1) = {Formulas.D1_Bandwidth(4, 1.0):0.000} -> 8J at large N");
Console.WriteLine($"  D4 crossing rhs (d-1)/2: d=2 -> {Formulas.D4_CrossingRhs(2):0.0}, d=4 -> {Formulas.D4_CrossingRhs(4):0.0}");
Console.WriteLine($"  D6 gap=2g={Formulas.D6_Gap(0.5):0.0}, mixing time (N=4) <= {Formulas.D6_MixingTime(4, 0.5):0.00};  F38 Pi^2=(-1)^(nY+nZ)=X^N (order 4)");
Console.WriteLine($"  F18 fold threshold Sg_crit/J: Bell={Formulas.F18_FoldThresholdBell}, product={Formulas.F18_FoldThresholdProduct}");
Console.WriteLine($"  F36/F37 neural (Wilson-Cowan): Q J Q + J + 2S = 0, pairing mu+mu'=-(1/tE+1/tI); C.elegans 0.013 vs 0.108 random");
Console.WriteLine($"  F61 bit_a parity Pi^2_X=Z^N=(-1)^(nX+nY)=(-1)^k (the k-parity H conserves); F63 [L,Pi^2]=0 (conserved, all N)");
Console.WriteLine($"  F39 det(Pi)=(-1)^(N 4^(N-1)): N=1->{Formulas.F39_DetPi(1)}, N>=2->{Formulas.F39_DetPi(4)};  F41 t_Pi=pi/(4J sin^2(pi/2N)) (N=4) = {Formulas.F41_PalindromicTime(4, 1.0):0.000}");
Console.WriteLine($"  F49 cross-term R(N)=sqrt((N-2)/(N 4^(N-1))): N=2,3,4 = {Formulas.F49_CrossTerm(2):0.000}, {Formulas.F49_CrossTerm(3):0.0000}, {Formulas.F49_CrossTerm(4):0.0000} (0, 1/sqrt48, 1/sqrt128);  F49b ||L_Dc||^2=g^2 4^N N (N=4,g=0.5) = {Formulas.F49b_CenteredDissipatorNormSq(4, 0.5):0.0}");
Console.WriteLine($"  F49c shadow-crossing R(N)=sqrt((N-1)/(N 4^(N-1))): N=3,4 = {Formulas.F49c_CrossTermCrossing(3):0.0000}, {Formulas.F49c_CrossTermCrossing(4):0.0000}");
Console.WriteLine($"  F55 K_death=ln10={Formulas.F55_KDeath:0.000} (99% absorption); immortal modes=N+1={Formulas.F55_ImmortalModes(4)} (N=4);  F56 cusp K(eps=1e-4,tol=1e-12) = {Formulas.F56_CriticalSlowing(1e-4, 1e-12):0.000}");
Console.WriteLine($"  F57 Bell+ dwell prefactor K_dwell/delta = {Formulas.F57_DwellPrefactorBell} (gamma-independent)");
Console.WriteLine($"  F60 GHZ_N CPsi(0)=1/(2^N-1): N=2,3,4,5 = {Formulas.F60_GhzCPsi0(2):0.000}, {Formulas.F60_GhzCPsi0(3):0.000}, {Formulas.F60_GhzCPsi0(4):0.000}, {Formulas.F60_GhzCPsi0(5):0.000} (<1/4 for N>=3)");
Console.WriteLine($"  F62 W_N CPsi(0)=2(N^2-4N+8)/(3N^3): N=2,3,4 = {Formulas.F62_WstateCPsi0(2):0.000}, {Formulas.F62_WstateCPsi0(3):0.000}, {Formulas.F62_WstateCPsi0(4):0.000};  F59 prefactor(W0=1/2,k=2,W2=0.371) = {Formulas.F59_DwellPrefactor(2, 0.5, 0.3709):0.000}");
Console.WriteLine($"  F63 [L,Pi^2]=0, 4 blocks dim 4^(N-1)={Formulas.F63_BlockDim(4)} (N=4); conserved/sector (even,odd) N=4 = {Formulas.F63_ConservedPerSector(4)}");
Console.WriteLine($"  F65 SE spectrum (4/(N+1))sin^2(k pi/(N+1)): N=3 = [{string.Join(", ", Formulas.F65_SingleExcitationRates(3).Select(x => x.ToString("0.000")))}], N=5 = [{string.Join(", ", Formulas.F65_SingleExcitationRates(5).Select(x => x.ToString("0.000")))}]");
Console.WriteLine($"  F65 Niven: rational iff N in {{0,1,2,3,5}} (N=4 golden-irrational); F66 pole multiplicity N+1 = {Formulas.F66_PoleMultiplicity(4)} (N=4)");
Console.WriteLine($"  F68 partner rate alpha_p=2g0-alpha_b (g0=0.05, alpha_b=0.0138) = {Formulas.F68_PartnerRate(0.0138, 0.05):0.0000};  F69 GHZ3+W3 optimum pair-CPsi(0) = {Formulas.F69_N3Optimum:0.000} (ratio {Formulas.F69_RatioToQuarter} to 1/4, irreducible sextic)");
Console.WriteLine($"  F70 k-local sees |dN|<=k (single-site<=1, pair<=2); F71 c1 components=floor(N/2)={Formulas.F71_C1IndependentComponents(5)} (N=5), R|psi_k>=(-1)^(k+1) [k=1->{Formulas.F71_ReflectionParity(1)}, k=2->{Formulas.F71_ReflectionParity(2)}]");
Console.WriteLine($"  F98 Dicke asymptote (N+2)/(4(N+1)): N=4,6,8 = {Formulas.F98_DickeAsymptote(4):0.000}, {Formulas.F98_DickeAsymptote(6):0.000}, {Formulas.F98_DickeAsymptote(8):0.000} (-> 1/4)");
Console.WriteLine($"  F121 qudit palindrome paired(d,N)/d^2N: d=2,N=2 = {Formulas.F121_PairedCeiling(2, 2)}/16 (full); d=3,N=2 = {Formulas.F121_PairedCeiling(3, 2)}/81 (partial, d^2-2d=0 re-seen)");
Console.WriteLine($"  F122 ceiling g2: K_5=4/N={Formulas.F122_CompleteCeiling(5):0.000}, star_6=4/(N-1)={Formulas.F122_StarCeiling(6):0.000}, K_4=2-2/sqrt3={Formulas.F122_K4Ceiling():0.000}, ring(1,1) N=5={Formulas.F122_RingCommutant(5):0.000}");
Console.WriteLine($"  F85 k-body trichotomy c(term) in {{0,1,2}}: YY->{Formulas.F85_FrobeniusFactor("YY")}, XY->{Formulas.F85_FrobeniusFactor("XY")}, YZ->{Formulas.F85_FrobeniusFactor("YZ")}, YYY->{Formulas.F85_FrobeniusFactor("YYY")} (n_YZ=3 yet c=1); Pi^2-odd count k=2,3,4 = {Formulas.F85_Pi2OddCount(2)}, {Formulas.F85_Pi2OddCount(3)}, {Formulas.F85_Pi2OddCount(4)}; ||M||^2/term = 4c ||H||^2 2^N");
Console.WriteLine($"  F97 cardioid c(phi)=z*(1-z*), |z*|=1/2: cusp c(0)={Formulas.F97_Cardioid(0).Real:0.00} (the F16 boundary), tail c(pi)={Formulas.F97_Cardioid(Math.PI).Real:0.00}, top c(pi/2)={Formulas.F97_Cardioid(Math.PI / 2).Real:0.00}+{Formulas.F97_Cardioid(Math.PI / 2).Imaginary:0.00}i");
Console.WriteLine($"  F124 band-edge invariant ||M||_F^2 + lambda_min = 2: N=4 -> {Formulas.F124_FrobeniusNormSq(4):0.0000} + {Formulas.F124_SpectralFloor(4):0.0000}; E = the F65 k=1 rung; E(N+1)^3 -> 4pi^2 = {4 * Math.PI * Math.PI:0.00}");
Console.WriteLine($"  F74 chromaticity at J=0: c(n,N)=min(n,N-1-n)+1 distinct rates, ladder 2g0*{{1,3,..,2c-1}}; N=8 row c(0..7) = {string.Join(",", Enumerable.Range(0, 8).Select(nn => Formulas.F74_Chromaticity(nn, 8)))} (pinned by direct pair enumeration in tests)");
Console.WriteLine($"  F75 mirror-pair MI = 2h(p)-h(2p) (Bell ceiling 2 bits at p=1/2); MM(0) bonding: (5,2)={Formulas.F75_MirrorPairSum(5, 2):0.000}, (7,4)={Formulas.F75_MirrorPairSum(7, 4):0.000}, (11,6)={Formulas.F75_MirrorPairSum(11, 6):0.000} (even k, noded centre, wins)");
Console.WriteLine($"  F77 MM(0) saturation: 1 + 3/(4(N+1)ln2), rescaled limit {Formulas.F77_RescaledDeviationLimit():0.0000}; closed form at N=101: {Formulas.F77_MMSaturation(101):0.00000} (exact best-k F75 sum pinned in tests)");
Console.WriteLine($"  F76 dephasing envelope lambda=e^(-4g0t): MM(t)/MM(0) at g0=0.05, t=0.1: (5,2)={Formulas.F76_Envelope(5, 2, 0.05, 0.1):0.000}, (13,4)={Formulas.F76_Envelope(13, 4, 0.05, 0.1):0.000}; the 0.93 is the g0 signature ({Formulas.F76_Envelope(5, 2, 0.025, 0.1):0.000} at g0=0.025, {Formulas.F76_Envelope(5, 2, 0.10, 0.1):0.000} at 0.10)");
Console.WriteLine($"  F95 theta-compass arctan(sqrt(c/b^2-1)): threshold b^2=1/4 at b=1/2; Lindblad face theta=arctan(Q) -- F95(g^2+J^2, g) at Q=2 = {Formulas.F95_Theta(0.25 + 1.0, 0.5) * 180 / Math.PI:0.0} deg = the clock (the compass and the clock are one)");
Console.WriteLine($"  F99 five canonical angles, alpha=sin^2/2: 0/30/45/60/90 deg -> {Formulas.F99_Alpha(0):0.000}, {Formulas.F99_Alpha(Math.PI / 6):0.000}, {Formulas.F99_Alpha(Math.PI / 4):0.000}, {Formulas.F99_Alpha(Math.PI / 3):0.000}, {Formulas.F99_Alpha(Math.PI / 2):0.000} (the Pi2 dyadic ladder); c^2(45deg) = 1+sqrt2 = {Formulas.F99_DickeWeightSq(Math.PI / 4):0.000} (the silver ratio)");
Console.WriteLine($"  F88b Pi^2-odd/memory anchors alpha: mirror(6;2,4)={Formulas.F88b_Alpha(6, 2, 4):0.0}, K-int(6;3,4)={Formulas.F88b_Alpha(6, 3, 4):0.000} (=F98: {Formulas.F98_DickeAsymptote(6):0.000}), generic(7;1,3)={Formulas.F88b_Alpha(7, 1, 3):0.0}; GHZ (HD=N) = {Formulas.F88b_Pi2OddInMemory(4, 0, 4, 4):0.0} (Pi^2-classical); Dicke alpha_total(gamma=1/2) = {Formulas.F88b_DickeAlphaTotal(0.5):0.000}");

// (the 4^N Pauli-string enumeration now lives in PauliMode.Enumerate, shared with the tests)
