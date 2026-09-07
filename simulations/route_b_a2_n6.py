"""Exact N=6 exported pencils and bounded-CRT discriminant layer identities.

Works directly in t. Block reconstruction and odd-N orbit premises stay outside.
N6_TABLE is a proposal until the coefficient-difference bound proves the identity.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import multiprocessing as mp
import os
from pathlib import Path
import re
import subprocess
import threading
import time

import numpy as np
import sympy as sp

if __package__:
    from . import o2b_gcd_certificate as arithmetic
else:
    import o2b_gcd_certificate as arithmetic

Lambda, t = sp.symbols("Lambda t")
N6_TABLE = {"degD": 926, "valuation": 536, "A1": 124, "A2": 133}


@dataclass(frozen=True)
class ParityPencil:
    parity: str
    residual: sp.Poly
    at: sp.Poly


@dataclass(frozen=True)
class LayerCertificate:
    deg_d: int
    valuation: int
    a1: sp.Poly
    a2: sp.Poly
    constant: int
    proof_modulus: int
    proof_bound: int
    proof_primes: tuple[int, ...]
    discriminant_lhs: sp.Expr
    discriminant_rhs: sp.Expr


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _keys(value, expected):
    if not isinstance(value, dict) or set(value) != set(expected):
        raise ValueError(f"expected exactly the fields {sorted(expected)}")


def _parse_poly(rows, degree):
    # Each exporter row is individually trimmed; rectangular padding is invalid.
    if not isinstance(rows, list) or len(rows) != degree + 1:
        raise ValueError("wrong Lambda degree / coefficient row count")
    terms = {}
    for i, row in enumerate(rows):
        if not isinstance(row, list) or not row or len(row) > degree-i+1:
            raise ValueError("malformed trimmed coefficient row")
        for j, value in enumerate(row):
            if not isinstance(value, str) or re.fullmatch(r"0|-?[1-9][0-9]*", value) is None:
                raise ValueError("coefficients must be real decimal integer strings")
            terms[i, j] = int(value)
        if len(row) > 1 and row[-1] == "0":
            raise ValueError("noncanonical trailing zero in trimmed row")
    poly = sp.Poly.from_dict(terms, (Lambda, t), domain=sp.ZZ)
    if poly.degree(Lambda) != degree or rows[-1] != ["1"]:
        raise ValueError("polynomial must have the specified degree and be monic in Lambda")
    return poly


def load_exact_pencils(path):
    """Parse schema 1 strictly; check independently exported E/O transports."""
    data = json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=_unique_object)
    _keys(data, {"schemaVersion", "n", "lambdaConvention", "parameterConvention",
                 "coefficientOrder", "rEven", "rOdd"})
    if (type(data["schemaVersion"]) is not int or data["schemaVersion"] != 1
            or type(data["n"]) is not int or data["n"] != 6
            or data["lambdaConvention"] != "Lambda=2*lambda"
            or data["parameterConvention"] != "t=i*qCSharp"
            or data["coefficientOrder"] != "lambda-lowest-first,t-lowest-first"):
        raise ValueError("unsupported schema, N, conventions or coefficient order")
    result = {}
    for key, parity, odd in (("rEven", "E", False), ("rOdd", "O", True)):
        obj = data[key]
        _keys(obj, {"rOdd", "sectorDimension", "residualInT", "atFactorInT"})
        if obj["rOdd"] is not odd or type(obj["sectorDimension"]) is not int or obj["sectorDimension"] != 45:
            raise ValueError("parity or sector dimension mismatch")
        result[parity] = ParityPencil(parity, _parse_poly(obj["residualInT"],32),
                                    _parse_poly(obj["atFactorInT"],13))
    for field in ("residual", "at"):
        e, o = getattr(result["E"],field), getattr(result["O"],field)
        expected = sp.Poly.from_dict({(i,j): (-1)**j*c for (i,j),c in e.terms()},
                                     (Lambda,t), domain=sp.ZZ)
        if o != expected:
            raise ValueError(f"E/O transport fails for {field}")
    return result


def _matrix(f):
    if not isinstance(f, sp.Poly) or f.gens != (Lambda,t) or f.domain != sp.ZZ:
        raise ValueError("expected Poly(Lambda,t) over ZZ")
    m, d = f.degree(Lambda), max(0,f.degree(t))
    if m < 2 or any(i == m and (j != 0 or c != 1) for (i,j),c in f.terms()):
        raise ValueError("expected a monic polynomial of Lambda degree at least two")
    rows = [[0]*(d+1) for _ in range(m+1)]
    for (i,j),c in f.terms():
        rows[i][j] = int(c)
    return rows


def discriminant_degree_bound(f):
    """Sylvester support bound independent of proposed discriminant degree.

    Each F row shifted by Lambda^s has total degree <= D+s, and each F'
    row <= D-1+s. A determinant term uses all columns 0..2m-2 once.
    Subtracting their Lambda degrees gives (2m-1)D-m^2. The ordinary
    row bound (2m-1)*deg_t(F) is also valid; use the smaller one.
    """
    _matrix(f)
    m = f.degree(Lambda)
    d = max(i+j for (i,j),c in f.terms())
    return min((2*m-1)*max(0,f.degree(t)), (2*m-1)*d-m*m)


def _det_bound(f):
    rows = _matrix(f)
    m = len(rows)-1
    norms = [sum(abs(c) for c in row) for row in rows]
    # Polynomial l1 is submultiplicative: det <= permanent <= row-sum product.
    return sum(norms)**(m-1)*sum(i*v for i,v in enumerate(norms))**m


def coefficient_difference_bound(f, constant, a1, a2):
    """Bound EVERY coefficient of disc(F)-C*t^v*A1*A2^2."""
    return (_det_bound(f)+abs(int(constant))*sum(abs(int(c)) for c in a1.all_coeffs())
            *sum(abs(int(c)) for c in a2.all_coeffs())**2)


def _disc_stream(f, primes, workers):
    m = f.degree(Lambda)
    sign = (-1)**(m*(m-1)//2)
    nodes = list(range(discriminant_degree_bound(f)+1))
    stream = arithmetic.disc_stream(_matrix(f), nodes, primes, workers)
    try:
        for p, resultant in stream:
            yield p, (resultant*sign)%p
    finally:
        stream.close()


def _canonical_layer(poly):
    if (not isinstance(poly,sp.Poly) or poly.gens != (t,) or poly.domain != sp.ZZ
            or poly.is_zero or poly.LC() <= 0 or poly.content() != 1 or poly.eval(0) == 0):
        raise ValueError("layers must be primitive, positive-leading ZZ[t], nonzero at zero")


def certify_layer_identity(f, valuation, a1, a2, constant=None, workers=1):
    """Public fail-closed judge: finite-prime agreement is never a certificate.

    Missing C is proposed by centered CRT of leading ratios past 2*B_D.
    Final equality with modulus > 2*(B_D+B_R) proves the integer identity.
    Every proof prime preserves degree/valuation. Squarefree coprime reduction
    with leading coefficients intact proves shape over Q by Gauss's lemma.
    """
    _matrix(f)
    if type(valuation) is not int or valuation < 0 or type(workers) is not int or workers < 1:
        raise ValueError("invalid valuation or worker count")
    _canonical_layer(a1)
    _canonical_layer(a2)
    if constant is not None and (type(constant) is not int or constant == 0):
        raise ValueError("C must be a nonzero integer")
    degree = valuation+a1.degree()+2*a2.degree()
    if degree > discriminant_degree_bound(f):
        raise ValueError("proposed identity exceeds discriminant support")
    lc = int(a1.LC()*a2.LC()**2)
    bound_d = _det_bound(f)
    pool, pending = [], []
    modulus, rc, mc = 1, 0, 1
    discriminant_residues = [0]*(degree+1)
    shape = False
    started = time.perf_counter()
    stream = _disc_stream(f, arithmetic.prime_stream(), workers)
    try:
        for p, disc in stream:
            if p in pool:
                raise ValueError("duplicate proof prime")
            if len(disc)-1 != degree or arithmetic.valuation_at_zero(disc) != valuation or lc%p == 0:
                raise ValueError(f"proof prime {p} does not preserve degree/valuation")
            pool.append(p)
            pending.append((p,disc))
            discriminant_residues = [arithmetic.crt_pair(old,modulus,int(new),p)
                                     for old,new in zip(discriminant_residues,disc)]
            modulus *= p
            if constant is None:
                ratio = int(disc[0])*pow(lc%p,-1,p)%p
                rc = arithmetic.crt_pair(rc,mc,ratio,p)
                mc *= p
                if mc <= 2*bound_d:
                    if len(pool)%32 == 0:
                        print(f"proof C: {len(pool)} primes, {modulus.bit_length()} bits, {time.perf_counter()-started:.1f}s",flush=True)
                    continue
                constant = rc if rc <= mc//2 else rc-mc
                if constant == 0:
                    raise ValueError("zero reconstructed constant")
            for q, dp in pending:
                x = arithmetic.sympy_poly_desc_modp(a1,q)
                y = arithmetic.sympy_poly_desc_modp(a2,q)
                rhs = arithmetic.polymul(arithmetic.polymul(y,y,q),x,q)
                rhs = np.concatenate([(rhs*(constant%q))%q,np.zeros(valuation,dtype=np.int64)])
                if not np.array_equal(arithmetic.polytrim(rhs),dp):
                    raise ValueError(f"integer layer identity fails modulo {q}")
                if not shape:
                    shape = (len(arithmetic.polygcd(x,arithmetic.polyderiv(x,q),q)) == 1
                             and len(arithmetic.polygcd(y,arithmetic.polyderiv(y,q),q)) == 1
                             and len(arithmetic.polygcd(x,y,q)) == 1)
            pending.clear()
            bound = coefficient_difference_bound(f,constant,a1,a2)
            if modulus > 2*bound:
                if not shape:
                    raise ValueError("no squarefree coprime layer shape in the completed proof pool")
                # These coefficients come from the determinant residues, never RHS.
                # modulus > 2*bound >= 2*B_D proves the centered lift over ZZ.
                lhs = sp.Poly.from_list([c if c <= modulus//2 else c-modulus
                                         for c in discriminant_residues],t,domain=sp.ZZ)
                rhs = constant*sp.Poly(t**valuation,t)*a1*a2**2
                if lhs != rhs:
                    raise ValueError("independently reconstructed discriminant differs from RHS")
                return LayerCertificate(degree,valuation,a1,a2,constant,modulus,bound,
                                        tuple(pool),lhs.as_expr(),rhs.as_expr())
    finally:
        stream.close()
    raise ValueError("prime stream exhausted before the integer proof")


def _lift_layers(f, workers):
    table = {"v":N6_TABLE["valuation"], "A1":N6_TABLE["A1"], "A2":N6_TABLE["A2"]}
    residues, modulus, candidate, stable = None,1,None,0
    started = time.perf_counter()
    stream = _disc_stream(f,arithmetic.prime_stream(),workers)
    try:
        for used,(p,disc) in enumerate(stream,1):
            if len(disc)-1 != N6_TABLE["degD"]:
                raise ValueError("modular discriminant degree differs from N6 proposal")
            try:
                first,second,_ = arithmetic.layers_from_disc(disc,p,table)
            except AssertionError as exc:
                raise ValueError(f"N6 layer proposal fails at prime {p}") from exc
            incoming = [list(map(int,first)),list(map(int,second))]
            if residues is None:
                residues = incoming
            else:
                residues = [[arithmetic.crt_pair(a,modulus,b,p) for a,b in zip(old,new)]
                            for old,new in zip(residues,incoming)]
            modulus *= p
            if used%8:
                continue
            print(f"lift: {used} primes, {modulus.bit_length()} bits, {time.perf_counter()-started:.1f}s",flush=True)
            proposed = []
            for layer in residues:
                rationals = [arithmetic.rational_reconstruct(c,modulus) for c in layer]
                if any(c is None for c in rationals):
                    break
                poly = sp.Poly.from_list([sp.Rational(a,b) for a,b in rationals],t,domain=sp.QQ)
                poly = poly.clear_denoms(convert=True)[1].primitive()[1]
                proposed.append(poly if poly.LC()>0 else -poly)
            if len(proposed) != 2:
                candidate,stable = None,0
            elif proposed == candidate:
                stable += 1
                if stable == 2:
                    return tuple(proposed)
            else:
                candidate,stable = proposed,0
    finally:
        stream.close()


def prove_layer_identity(parity_pencil, workers=1):
    """Lift E's proposed layers, then prove the full integer identity."""
    if parity_pencil.parity != "E":
        raise ValueError("prove E; O is transported only by the strict loader")
    a1,a2 = _lift_layers(parity_pencil.residual,workers)
    if (a1.degree(),a2.degree()) != (124,133):
        raise ValueError("lifted layer degree mismatch")
    return certify_layer_identity(parity_pencil.residual,536,a1,a2,workers=workers)


class _ProcessTreeSampler:
    """Sample simultaneous process-tree RSS and OS-reported per-process peaks.

    psutil is a benchmark-only dependency. Polling can miss short-lived children;
    the JSON states the sampling interval instead of claiming an exact tree peak.
    """
    def __init__(self, pid):
        import psutil
        self.psutil, self.pid = psutil, pid
        self.peak_tree = 0
        self.peaks = {}
        self.stop_event = threading.Event()
        self.thread = threading.Thread(target=self._run, daemon=True)

    def _sample(self):
        try:
            root = self.psutil.Process(self.pid)
            processes = [root]+root.children(recursive=True)
        except self.psutil.Error:
            return
        current = 0
        for process in processes:
            try:
                info = process.memory_info()
                current += info.rss
                key = f"{process.pid}:{process.create_time()}"
                peak = getattr(info,"peak_wset",info.rss)
                self.peaks[key] = max(self.peaks.get(key,0),peak)
            except self.psutil.Error:
                continue
        self.peak_tree = max(self.peak_tree,current)

    def _run(self):
        while not self.stop_event.is_set():
            self._sample()
            self.stop_event.wait(0.02)

    def __enter__(self):
        self._sample()
        self.thread.start()
        return self

    def __exit__(self, *args):
        self._sample()
        self.stop_event.set()
        self.thread.join()

    def report(self):
        return {"peak_tree_rss_bytes":self.peak_tree,
                "observed_process_peak_wset_bytes":self.peaks,
                "memory_sampling_seconds":0.02}


def recommend_workers(samples, available_memory):
    """Choose measured projected speed among runs using <= half available RAM."""
    eligible = [s for s in samples if s["peak_tree_rss_bytes"] <= available_memory//2]
    if not eligible:
        raise ValueError("no measured worker configuration fits the memory budget")
    return min(eligible,key=lambda s:(s.get("projected_proof_seconds",s["seconds_per_prime"]),
                                      s["workers"]))["workers"]


def _benchmark_probe(kind, poly, queue):
    try:
        if kind == "s1":
            sequence = sp.subresultants(poly.as_expr(),poly.diff(Lambda).as_expr(),Lambda)
            result = {"linear_entries":sum(sp.degree(s,Lambda)==1 for s in sequence)}
        elif kind == "real":
            intervals = poly.intervals(eps=sp.Rational(1,10**6))
            if not intervals:
                raise ValueError("N6 A2 candidate has no real root to probe")
            interval,multiplicity = intervals[0]
            refined = poly.refine_root(*interval,eps=sp.Rational(1,10**20))
            result = {"isolation_interval":str(interval),"refined_interval":str(refined),
                      "multiplicity":int(multiplicity)}
        else:
            # Floating roots propose boxes only; exact argument-principle counts
            # decide whether a box isolates a root of this integer polynomial.
            coefficients = np.array([float(c/poly.LC()) for c in poly.all_coeffs()])
            guesses = sorted((z for z in np.roots(coefficients) if z.imag > 0.001),key=abs)
            result = None
            for z in guesses[:8]:
                re = sp.Rational(str(round(float(z.real),10)))
                im = sp.Rational(str(round(float(z.imag),10)))
                radius = sp.Rational(1,10**4)
                lower,upper = re-radius+sp.I*(im-radius),re+radius+sp.I*(im+radius)
                if poly.count_roots(lower,upper) != 1:
                    continue
                tight = sp.Rational(1,10**8)
                lo,hi = re-tight+sp.I*(im-tight),re+tight+sp.I*(im+tight)
                if poly.count_roots(lo,hi) == 1:
                    result = {"isolation_box":[str(lower),str(upper)],
                              "refined_box":[str(lo),str(hi)],"exact_root_count":1}
                    break
            if result is None:
                raise ValueError("numeric proposals did not yield an exact isolated refined root")
        queue.put({"status":"complete",**result})
    except Exception as exc:
        queue.put({"status":"error","error":str(exc)})


def _timed_probe(kind, poly, budget=30):
    queue = mp.Queue()
    process = mp.Process(target=_benchmark_probe,args=(kind,poly,queue))
    started = time.perf_counter()
    process.start()
    with _ProcessTreeSampler(process.pid) as memory:
        process.join(budget)
        if process.is_alive():
            process.terminate()
            process.join()
            report = {"status":"timeout","budget_seconds":budget}
        else:
            try:
                report = queue.get(timeout=2)
            except Exception as exc:
                report = {"status":"error","error":str(exc),"exitcode":process.exitcode}
    queue.close()
    return {**report,"seconds":time.perf_counter()-started,**memory.report()}


def benchmark(path):
    """Measure actual export, N6 arithmetic, child memory and unresolved S1 work.

    The fresh A2 lift is only a timing candidate. This benchmark never promotes
    it to an exact N6 layer; prove_layer_identity is the separate certificate.
    """
    import psutil
    root = Path(__file__).resolve().parents[1]
    output = root/"simulations/results/route_b_a2_n6_benchmark.txt"
    fresh = root/"simulations/results/route_b_a2_n6_benchmark_residual.tmp.json"
    report = {"schemaVersion":1}
    command = ["dotnet","run","--project","compute/RCPsiSquared.Cli","-c","Release",
               "--","route-b-n6-residual","--out",str(fresh)]
    started = time.perf_counter()
    exporter = subprocess.Popen(command,cwd=root,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    with _ProcessTreeSampler(exporter.pid) as memory:
        stdout,stderr = exporter.communicate()
    report["export"] = {"seconds":time.perf_counter()-started,"returncode":exporter.returncode,
                        "stdout":stdout,"stderr":stderr,**memory.report()}
    if exporter.returncode:
        output.write_text(json.dumps(report,indent=2),encoding="utf-8")
        raise RuntimeError(f"fresh C# export failed: {stderr}")
    print(f"benchmark export: {report['export']['seconds']:.3f}s",flush=True)
    started = time.perf_counter()
    pencils = load_exact_pencils(fresh)
    if pencils != load_exact_pencils(path):
        raise ValueError("fresh export differs from the supplied exact pencils")
    e = pencils["E"]
    report["load_seconds"] = time.perf_counter()-started
    started = time.perf_counter()
    with _ProcessTreeSampler(os.getpid()) as memory:
        stream = _disc_stream(e.residual,arithmetic.prime_stream(),1)
        p,disc = next(stream)
        stream.close()
        layers = arithmetic.layers_from_disc(disc,p,{"v":536,"A1":124,"A2":133})
    report["one_prime"] = {"seconds":time.perf_counter()-started,"prime":p,
                           "degree":len(disc)-1,"valuation":layers[2],**memory.report()}
    samples = []
    for workers in (1,2,4,8):
        if workers > (psutil.cpu_count(logical=False) or os.cpu_count() or 1):
            continue
        started = time.perf_counter()
        with _ProcessTreeSampler(os.getpid()) as memory:
            stream = _disc_stream(e.residual,arithmetic.prime_stream(),workers)
            next(stream)
            startup = time.perf_counter()-started
            steady = time.perf_counter()
            for _ in range(32):
                next(stream)
            per_prime = (time.perf_counter()-steady)/32
            stream.close()
        # The row-sum determinant bound estimates the actual proof pool size.
        count = ((_det_bound(e.residual)*2).bit_length()+24)//25
        samples.append({"workers":workers,"startup_seconds":startup,
                        "seconds_per_prime":per_prime,"sample_primes":32,
                        "projected_proof_seconds":startup+count*per_prime,**memory.report()})
        print(f"benchmark workers={workers}: steady {per_prime:.3f}s/prime",flush=True)
    available = psutil.virtual_memory().available
    workers = recommend_workers(samples,available)
    report.update(worker_samples=samples,recommended_workers=workers,available_memory_bytes=available)
    report["s1"] = _timed_probe("s1",e.residual)
    print(f"benchmark exact S1: {report['s1']['status']}",flush=True)
    unresolved = report["s1"]["status"] != "complete"
    report["task4_checkpoint"] = {
        "status":"REQUIRES_REVIEWED_ALTERNATIVE" if unresolved else "REQUIRES_S1_VALIDATION",
        "required_route":"CRT_WITH_RIGOROUS_COEFFICIENT_BOUND" if unresolved else "EXACT_S1_VALIDATION",
        "s1_proved":False,
        "reason":"A timing probe is not the Task4 uniqueness/nonvanishing certificate."}
    started = time.perf_counter()
    with _ProcessTreeSampler(os.getpid()) as memory:
        _,a2 = _lift_layers(e.residual,workers)
    report["roots"] = {"source":"fresh N6 E A2 CRT candidate; timing only",
                       "candidate_degree":a2.degree(),"lift_seconds":time.perf_counter()-started,
                       "lift_memory":memory.report(),"probes":{}}
    for kind in ("real","nonreal"):
        report["roots"]["probes"][kind] = _timed_probe(kind,a2)
        print(f"benchmark N6 A2 {kind}: {report['roots']['probes'][kind]['status']}",flush=True)
    output.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2),flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--residual-json",type=Path,required=True)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--benchmark",action="store_true")
    action.add_argument("--prove-layers",action="store_true")
    parser.add_argument("--workers",type=int,default=1)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    if args.benchmark:
        benchmark(args.residual_json)
        return
    started = time.perf_counter()
    cert = prove_layer_identity(load_exact_pencils(args.residual_json)["E"],args.workers)
    print(f"EXACT: degD={cert.deg_d} valuation={cert.valuation} A1={cert.a1.degree()} A2={cert.a2.degree()}")
    print(f"C={cert.constant}")
    print(f"proof_primes={len(cert.proof_primes)} proof_modulus={cert.proof_modulus} proof_bound={cert.proof_bound}")
    print(f"A1={cert.a1.as_expr()}\nA2={cert.a2.as_expr()}")
    print(f"runtime_seconds={time.perf_counter()-started:.3f}")


if __name__ == "__main__":
    main()
