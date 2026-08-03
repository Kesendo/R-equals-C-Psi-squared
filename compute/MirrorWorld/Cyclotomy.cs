namespace MirrorWorld;

// The lattice our angles are forced onto (built 2026-08-03).
//
// WHY THIS EXISTS AS AN OBJECT AND NOT AS A PARAGRAPH. The repo has a typed claim
// (NivenRationalityRootClaim, Tier1Derived) whose arithmetic is sound, and two carbon documents
// that dress the same arithmetic in physical labels that do not survive review: a "constructible
// angle" cut where Niven needs "rational multiple of pi", a completeness claim refuted by a state
// the same mechanism reaches (alpha = 1/16 at c^2 = 7 + 2*sqrt(14)), and a HOMO table computing the
// top antibonding level. Repairing the prose would trade strings for strings. So the arithmetic is
// adopted here, from the CLAIM and from the graph geometry, and never from those documents.
//
// WHAT IS OWN AND WHAT IS INHERITED, which is the whole point of putting it here:
//
//   OWN        the turn fractions themselves. A path graph with Dirichlet ends puts the single
//              excitation on theta_k = k*pi/(N+1) (Chebyshev); a ring puts it on 2*pi*k/N. That is
//              OUR choice of graph and boundary condition, and it is what hands the arithmetic
//              something to act on.
//   INHERITED  everything the arithmetic then says: which of those angles have a rational cosine,
//              and what degree the rest have. Those hold for ANY quantity of the form cos(q*pi),
//              q rational, and know nothing about dephasing, Lindblad, or carbon.
//
// So the statement worth having is not "the values are irrational" (inherited, and true of
// everything) but "we built an object whose spectrum is FORCED onto a cyclotomic lattice" (own).
//
// THE PARENT IS DELIBERATELY OPEN, and that is a finding rather than an omission. Own/Inherited are
// lists of OUTPUT NAMES. "2*cos is rational exactly at order 1, 2, 3, 4, 6" is not an output of any
// object; it is a PREDICATE OVER THE RANGE of one. The two-bucket ontology has no slot for a
// predicate, and World.Parent is null with a fixed Own, so nothing can sit below the frame either.
// Passing null here says exactly that: the arithmetic is not inherited from anything of ours,
// because nothing of ours is upstream of it. Whether the repo wants a third bucket is the question
// this object is meant to put in front of us in code rather than in prose.
//
// EXACT, NOT FLOATING POINT. Every quantity below is an integer or a set of integers: a turn
// fraction, its order, an Euler totient, a degree. No cosine is ever evaluated. That is Seed's and
// Divisor's genre and it is what lets this walk past the wall to any N.
public sealed class Cyclotomy : GameObject
{
    // Niven's theorem, in the only form this object needs: for a rational turn fraction j/n of
    // reduced order d, 2*cos(2*pi*j/n) is RATIONAL exactly when d is one of these, giving the five
    // values 2, -2, -1, 0, 1. INHERITED: this is a fact about Q and the cyclotomic fields, not
    // about us. It is stated once, here, so no document has to restate it.
    public static readonly IReadOnlyList<int> RationalOrders = new[] { 1, 2, 3, 4, 6 };

    public Cyclotomy() : base(null) { }

    public override IReadOnlyList<string> Own => new[]
    {
        "turn fractions of the path (k/(2(N+1)))",
        "turn fractions of the ring (k/N)",
    };

    public static int Gcd(int a, int b) { while (b != 0) (a, b) = (b, a % b); return a; }

    public static int Totient(int n)
    {
        int result = n;
        for (int p = 2; (long)p * p <= n; p++)
            if (n % p == 0)
            {
                while (n % p == 0) n /= p;
                result -= result / p;
            }
        if (n > 1) result -= result / n;
        return result;
    }

    // The order of a turn fraction: the smallest n' with j/n = j'/n'. Everything else keys off it.
    public static int Order(int j, int n) => n / Gcd(Math.Abs(j) % n == 0 ? n : Math.Abs(j), n);

    public static bool IsRational(int order) => RationalOrders.Contains(order);

    // [Q(2*cos(2*pi/d)) : Q]. The real subfield of the d-th cyclotomic field has degree
    // phi(d)/2 for d >= 3, and 1 for d = 1, 2. INHERITED, entirely.
    public static int Degree(int order) => order <= 2 ? 1 : Totient(order) / 2;

    // OWN: the single excitation on an open path of N sites sits at theta_k = k*pi/(N+1), which as a
    // turn fraction is k/(2(N+1)). The Dirichlet ends are ours; the Chebyshev form follows.
    public static IReadOnlyList<int> PathOrders(int n)
        => Enumerable.Range(1, n).Select(k => Order(k, 2 * (n + 1))).ToList();

    // OWN: the dissipator rates carry sin^2(theta_k), and sin^2 = (1 - cos 2*theta)/2, so their
    // arithmetic sits on the DOUBLED angle: turn fraction k/(N+1), not k/(2(N+1)). Keeping the two
    // apart is not pedantry -- the registry once justified this face on sin^2 directly and had to be
    // corrected, and the five-versus-three anchor confusion in the carbon documents is the same
    // doubling read in the other direction.
    public static IReadOnlyList<int> PathRateOrders(int n)
        => Enumerable.Range(1, n).Select(k => Order(k, n + 1)).ToList();

    // OWN: the ring closes, so the fractions are k/N.
    public static IReadOnlyList<int> RingOrders(int n)
        => Enumerable.Range(1, n).Select(k => Order(k, n)).ToList();

    public static bool AllRational(IReadOnlyList<int> orders) => orders.All(IsRational);

    public static int MaxDegree(IReadOnlyList<int> orders) => orders.Max(Degree);
}
