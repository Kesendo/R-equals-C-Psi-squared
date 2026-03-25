#!/usr/bin/env python3
"""
IBM Quantum Historical Calibration Analysis — Corrected Physics
================================================================
KEY INSIGHT: For single transmon qubits, purity is NOT monotonically decreasing.
T1 relaxation drives state toward |0> (pure), competing with T2 dephasing.
Purity passes through a MINIMUM then rises back to 1.

Parameter r = T2/(2*T1) determines minimum:
  r < r* ~ 0.213: Purity crosses C*Psi = 1/4
  r > r* ~ 0.213: NEVER reaches 1/4

IBM transmons: r ~ 0.4-0.7 -> NEVER cross in purity sense.

Usage:
    python ibm_history_analysis.py --mode synthetic
    python ibm_history_analysis.py --mode full --token YOUR_TOKEN
"""
import os, csv, json
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from datetime import datetime, timedelta, timezone
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from pathlib import Path
from dotenv import load_dotenv

# Load .env from script directory
load_dotenv(Path(__file__).parent / ".env")

IBM_TOKEN = os.environ.get("IBM_QUANTUM_TOKEN", "YOUR_TOKEN_HERE")
BACKEND_NAME = "ibm_torino"
OUTPUT_DIR = Path(__file__).parent / "results"
OUTPUT_DIR.mkdir(exist_ok=True)
HISTORY_DAYS = 180
SAMPLE_EVERY_N_DAYS = 1
R_STAR = 0.212755
LN4 = np.log(4)

def purity_single_qubit(t, T1, T2):
    if t <= 0: return 1.0
    e1 = np.exp(-t/T1) if T1 > 0 else 0.0
    e2 = np.exp(-t/T2) if T2 > 0 else 0.0
    rho00 = 1.0 - 0.5*e1
    rho11 = 0.5*e1
    return rho00**2 + rho11**2 + 0.5*e2**2

def cpsi_from_purity(P, d=2):
    return (P - 1.0/d) / (1.0 - 1.0/d)

def cpsi(t, T1, T2):
    return cpsi_from_purity(purity_single_qubit(t, T1, T2))

def find_purity_minimum(T1, T2):
    if T1 <= 0 or T2 <= 0: return None, None, None
    result = minimize_scalar(lambda t: purity_single_qubit(t, T1, T2),
                             bounds=(0.01, max(5*T2, 5*T1)), method='bounded')
    t_min, P_min = result.x, result.fun
    return t_min, P_min, cpsi_from_purity(P_min)

def find_quarter_crossing(T1, T2):
    t_min, P_min, c_min = find_purity_minimum(T1, T2)
    if c_min is None or c_min >= 0.25: return None, None
    try:
        t_star = brentq(lambda t: cpsi(t, T1, T2) - 0.25, 0.01, t_min, xtol=0.01)
        return t_star, t_star/T2
    except: return None, None

def coherence_crossing_time(T2): return T2 * LN4
def r_parameter(T1, T2): return T2/(2.0*T1) if T1 > 0 else None

def n_qubit_ghz_min_cpsi(T1, T2, N):
    T2_eff = T2 / N
    _, _, c_min = find_purity_minimum(T1, T2_eff)
    return c_min

def compute_qubit_record(date_str, qubit_id, T1_us, T2_us, freq_ghz=None):
    r = r_parameter(T1_us, T2_us)
    t_min, P_min, cpsi_min = find_purity_minimum(T1_us, T2_us)
    t_star, t_star_ratio = find_quarter_crossing(T1_us, T2_us)
    t_coh = coherence_crossing_time(T2_us)
    dist = (cpsi_min - 0.25) if cpsi_min is not None else None
    crosses = cpsi_min < 0.25 if cpsi_min is not None else False
    return {
        'date': date_str, 'qubit': qubit_id,
        'T1_us': round(T1_us,3), 'T2_us': round(T2_us,3),
        'frequency_GHz': round(freq_ghz,6) if freq_ghz else None,
        'r_param': round(r,6) if r else None,
        't_min_us': round(t_min,3) if t_min else None,
        't_min_over_T2': round(t_min/T2_us,6) if t_min else None,
        'purity_min': round(P_min,6) if P_min else None,
        'cpsi_min': round(cpsi_min,6) if cpsi_min else None,
        'distance_from_quarter': round(dist,6) if dist is not None else None,
        'crosses_quarter': crosses,
        't_star_us': round(t_star,3) if t_star else None,
        't_star_over_T2': round(t_star_ratio,6) if t_star_ratio else None,
        't_coherence_crossing_us': round(t_coh,3),
        't_coh_over_T2': round(LN4,6),
    }

def collect_historical_data(token, days=180):
    from qiskit_ibm_runtime import QiskitRuntimeService
    print(f"Connecting to {BACKEND_NAME}...")
    service = QiskitRuntimeService(channel="ibm_quantum_platform", token=token)
    backend = service.backend(BACKEND_NAME)
    num_qubits = backend.num_qubits
    print(f"Backend: {BACKEND_NAME}, {num_qubits} qubits")
    end_date = datetime.now(timezone.utc)
    start_date = end_date - timedelta(days=days)
    print(f"History: {start_date.date()} to {end_date.date()}")
    all_records, errors = [], 0
    current_date = start_date
    day_count = 0
    while current_date <= end_date:
        day_count += 1
        date_str = current_date.strftime("%Y-%m-%d")
        try:
            target_hist = backend.target_history(datetime=current_date)
            for q in range(num_qubits):
                try:
                    qp = target_hist.qubit_properties[q]
                    if qp is None: continue
                    t1, t2 = qp.t1, qp.t2
                    if t1 is None or t2 is None or t1 <= 0 or t2 <= 0: continue
                    freq = qp.frequency
                    rec = compute_qubit_record(date_str, q, t1*1e6, t2*1e6,
                                               freq/1e9 if freq else None)
                    all_records.append(rec)
                except: continue
            if day_count % 10 == 0:
                print(f"  {date_str}: {len(all_records)} records")
        except Exception as e:
            errors += 1
            if errors <= 5: print(f"  {date_str} failed: {e}")
        current_date += timedelta(days=SAMPLE_EVERY_N_DAYS)
    print(f"Done: {len(all_records)} records, {errors} errors")
    return all_records

def generate_synthetic_data(n_qubits=133, n_days=180):
    print(f"Synthetic: {n_qubits} qubits x {n_days} days...")
    np.random.seed(42)
    records = []
    start = datetime(2025, 8, 1)
    base_T1 = np.random.lognormal(np.log(170), 0.3, n_qubits)
    base_T2 = np.random.lognormal(np.log(200), 0.35, n_qubits)
    base_T2 = np.minimum(base_T2, 1.95*base_T1)
    base_freq = np.random.uniform(4.5, 5.2, n_qubits)
    for day in range(n_days):
        ds = (start + timedelta(days=day)).strftime("%Y-%m-%d")
        for q in range(n_qubits):
            T1 = base_T1[q] * np.random.lognormal(0, 0.08)
            T2 = base_T2[q] * np.random.lognormal(0, 0.10)
            T2 = min(T2, 1.98*T1)
            if np.random.random() < 0.02: T2 *= np.random.uniform(0.15, 0.5)
            T1, T2 = max(T1, 5.0), max(T2, 2.0)
            records.append(compute_qubit_record(ds, q, T1, T2, base_freq[q]))
        if (day+1) % 30 == 0: print(f"  Day {day+1}/{n_days}")
    print(f"Generated {len(records)} records")
    return records

def save_csv(records, fn="ibm_torino_history.csv"):
    fp = OUTPUT_DIR / fn
    if not records: return fp
    with open(fp, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=records[0].keys())
        w.writeheader(); w.writerows(records)
    print(f"Saved {len(records)} -> {fp}")
    return fp

def load_csv(fn="ibm_torino_history.csv"):
    fp = OUTPUT_DIR / fn
    records = []
    with open(fp) as f:
        for row in csv.DictReader(f):
            for k in ['T1_us','T2_us','frequency_GHz','r_param','t_min_us',
                       't_min_over_T2','purity_min','cpsi_min','distance_from_quarter',
                       't_star_us','t_star_over_T2','t_coherence_crossing_us','t_coh_over_T2']:
                row[k] = float(row[k]) if row.get(k) and row[k] not in ('None','') else None
            row['qubit'] = int(row['qubit'])
            row['crosses_quarter'] = row['crosses_quarter'] == 'True'
            records.append(row)
    print(f"Loaded {len(records)} from {fp}")
    return records

def analyze(records):
    valid = [r for r in records if r['cpsi_min'] is not None]
    r_vals = np.array([r['r_param'] for r in valid])
    cpsi_mins = np.array([r['cpsi_min'] for r in valid])
    dists = np.array([r['distance_from_quarter'] for r in valid])
    T1s = np.array([r['T1_us'] for r in valid])
    T2s = np.array([r['T2_us'] for r in valid])
    crosses = np.array([r['crosses_quarter'] for r in valid])
    t_min_ratios = np.array([r['t_min_over_T2'] for r in valid if r['t_min_over_T2']])
    n_cross = crosses.sum()

    print("\n" + "="*65)
    print("  IBM TORINO — ¼-GRENZ-ANALYSE (R = CΨ²)")
    print("="*65)
    print(f"  Datenpunkte:  {len(valid):,}")
    print(f"  Qubits:       {len(set(r['qubit'] for r in valid))}")
    print(f"  Zeitraum:     {valid[0]['date']} — {valid[-1]['date']}")
    print(f"\n  r = T₂/(2T₁): mean={r_vals.mean():.4f}, median={np.median(r_vals):.4f}")
    print(f"  r range: [{r_vals.min():.4f}, {r_vals.max():.4f}]")
    print(f"  r* (kritisch) = {R_STAR:.4f}")
    print(f"  Anteil r < r*: {(r_vals<R_STAR).sum()}/{len(valid)} ({100*(r_vals<R_STAR).mean():.1f}%)")
    print(f"\n  C·Ψ_min: mean={cpsi_mins.mean():.4f}, median={np.median(cpsi_mins):.4f}")
    print(f"  Abstand ¼: mean={dists.mean():.4f}")
    print(f"  Kreuzt ¼: {n_cross:,}/{len(valid):,} ({100*n_cross/len(valid):.2f}%)")
    if n_cross > 0:
        cq = set(r['qubit'] for r in valid if r['crosses_quarter'])
        cr = [r['t_star_over_T2'] for r in valid if r['crosses_quarter'] and r['t_star_over_T2']]
        print(f"  Crossing t*/T₂: {np.mean(cr):.4f} ± {np.std(cr):.4f}")
        print(f"  Qubits: {sorted(cq)[:15]}{'...' if len(cq)>15 else ''}")
    print(f"\n  Multi-Qubit GHZ (T1={np.median(T1s):.0f}, T2={np.median(T2s):.0f} μs):")
    for N in [1,2,3,4,5,8,10,20]:
        cm = n_qubit_ghz_min_cpsi(np.median(T1s), np.median(T2s), N)
        x = "✓" if cm is not None and cm < 0.25 else "✗"
        cv = f"{cm:.4f}" if cm is not None else "n/a"
        print(f"    N={N:>2}: C·Ψ_min={cv} {x}")
    return {'records':valid,'r_vals':r_vals,'cpsi_mins':cpsi_mins,'dists':dists,
            'T1s':T1s,'T2s':T2s,'crosses':crosses,'t_min_ratios':t_min_ratios}

def gen_theory():
    rs = np.linspace(0.005, 0.99, 300)
    cms = []
    for r in rs:
        T2=100; T1=T2/(2*r)
        _,_,cm = find_purity_minimum(T1, T2)
        cms.append(cm)
    return rs, np.array(cms)

def plot_all(data):
    rv, cm, ds = data['r_vals'], data['cpsi_mins'], data['dists']
    T1s, T2s, crosses = data['T1s'], data['T2s'], data['crosses']
    records = data['records']
    rs_th, cm_th = gen_theory()

    fig, axes = plt.subplots(2, 2, figsize=(16,14))
    fig.suptitle('IBM Torino — Wie nah kommt jeder Qubit an C·Ψ = ¼?', fontsize=15, fontweight='bold')

    ax = axes[0][0]
    ax.scatter(rv, cm, s=2, alpha=0.15, c='steelblue', rasterized=True)
    ax.plot(rs_th, cm_th, 'r-', lw=3, label='Theorie', zorder=10)
    ax.axhline(0.25, color='darkred', lw=2, ls='--', alpha=0.7, label='C·Ψ = ¼')
    ax.axvline(R_STAR, color='orange', lw=2, ls=':', alpha=0.7, label=f'r* = {R_STAR:.3f}')
    ax.fill_between([0,R_STAR], 0, 0.25, alpha=0.08, color='green')
    ax.set_xlabel('r = T₂/(2T₁)'); ax.set_ylabel('C·Ψ_min')
    ax.set_title('Minimale Quantizität vs Dekohärenzcharakter')
    ax.legend(fontsize=8); ax.set_xlim(0,1); ax.set_ylim(0,1); ax.grid(True, alpha=0.3)

    ax = axes[0][1]
    ax.hist(ds, bins=80, color='coral', alpha=0.7, edgecolor='black', lw=0.3, density=True)
    ax.axvline(0, color='darkred', lw=2.5)
    nb = (ds<0).sum(); na = (ds>=0).sum()
    ax.text(0.95, 0.95, f'Über ¼: {na:,}\nUnter ¼: {nb:,}', transform=ax.transAxes,
            fontsize=10, va='top', ha='right', bbox=dict(boxstyle='round',facecolor='lightyellow',alpha=0.8))
    ax.set_xlabel('C·Ψ_min − 0.25'); ax.set_ylabel('Dichte')
    ax.set_title('Abstand zur ¼-Grenze'); ax.grid(True, alpha=0.3)

    ax = axes[1][0]
    ax.hist(rv, bins=80, color='teal', alpha=0.7, edgecolor='black', lw=0.3, density=True)
    ax.axvline(R_STAR, color='orange', lw=2.5, label=f'r*={R_STAR:.3f}')
    ax.axvline(rv.mean(), color='red', lw=1.5, ls='--', label=f'Mean={rv.mean():.3f}')
    pct = 100*(rv<R_STAR).mean()
    ax.text(0.95, 0.95, f'{pct:.1f}% unter r*', transform=ax.transAxes, fontsize=10,
            va='top', ha='right', bbox=dict(boxstyle='round',facecolor='lightyellow',alpha=0.8))
    ax.set_xlabel('r = T₂/(2T₁)'); ax.set_ylabel('Dichte')
    ax.set_title('Verteilung r'); ax.legend(fontsize=10); ax.grid(True, alpha=0.3)

    ax = axes[1][1]
    sc = ax.scatter(T1s, T2s, c=ds, s=2, alpha=0.2, cmap='RdYlGn_r', vmin=-0.1, vmax=0.6, rasterized=True)
    t1l = np.linspace(0, T1s.max(), 100)
    ax.plot(t1l, 2*t1l, 'k--', lw=1.5, alpha=0.3, label='T₂=2T₁')
    ax.plot(t1l, 2*R_STAR*t1l, 'orange', lw=2, alpha=0.7, label=f'r=r*')
    ax.set_xlabel('T₁ (μs)'); ax.set_ylabel('T₂ (μs)')
    ax.set_title('T₁-T₂ Landschaft'); ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
    fig.colorbar(sc, ax=ax, shrink=0.8).set_label('C·Ψ_min − 0.25')
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR/'quarter_boundary_analysis.png', dpi=150, bbox_inches='tight')
    plt.close(fig); print(f"  -> quarter_boundary_analysis.png")

    # Time series
    fig2, axes2 = plt.subplots(3, 1, figsize=(16,14))
    fig2.suptitle('IBM Torino — Zeitreihen', fontsize=14, fontweight='bold')
    dd = {}
    for r in records:
        d=r['date']
        if d not in dd: dd[d]={'c':[],'r':[],'d':[]}
        dd[d]['c'].append(r['cpsi_min']); dd[d]['r'].append(r['r_param']); dd[d]['d'].append(r['distance_from_quarter'])
    sd = sorted(dd.keys()); xd = [datetime.strptime(d,'%Y-%m-%d') for d in sd]

    ax=axes2[0]
    ax.fill_between(xd, [np.percentile(dd[d]['c'],5) for d in sd], [np.percentile(dd[d]['c'],95) for d in sd], alpha=0.15, color='steelblue')
    ax.plot(xd, [np.median(dd[d]['c']) for d in sd], 'b-', lw=1.5, label='Median')
    ax.plot(xd, [np.min(dd[d]['c']) for d in sd], 'g-', lw=1, alpha=0.6, label='Min')
    ax.axhline(0.25, color='darkred', lw=2, ls='--', label='¼')
    ax.set_ylabel('C·Ψ_min'); ax.legend(fontsize=8, ncol=4); ax.grid(True, alpha=0.3)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

    ax=axes2[1]
    ax.fill_between(xd, [np.percentile(dd[d]['r'],5) for d in sd], [np.percentile(dd[d]['r'],95) for d in sd], alpha=0.15, color='teal')
    ax.plot(xd, [np.median(dd[d]['r']) for d in sd], '-', color='teal', lw=1.5)
    ax.plot(xd, [np.min(dd[d]['r']) for d in sd], 'g-', lw=1, alpha=0.6)
    ax.axhline(R_STAR, color='orange', lw=2, ls='--', label=f'r*={R_STAR:.3f}')
    ax.set_ylabel('r'); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')

    ax=axes2[2]
    cpd = [sum(1 for r in records if r['date']==d and r['crosses_quarter']) for d in sd]
    ax.bar(xd, cpd, color='coral', alpha=0.7)
    ax.set_ylabel('Crossings/Tag'); ax.grid(True, alpha=0.3, axis='y')
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    fig2.tight_layout()
    fig2.savefig(OUTPUT_DIR/'time_series.png', dpi=150, bbox_inches='tight')
    plt.close(fig2); print(f"  -> time_series.png")

    # Multi-qubit
    fig3, axes3 = plt.subplots(1, 2, figsize=(16,7))
    fig3.suptitle('Wann kreuzt GHZ die ¼-Grenze?', fontsize=14, fontweight='bold')
    ax=axes3[0]; T2r=np.median(T2s); Nr=np.arange(1,21)
    for rv2,c,l in [(0.05,'green','r=0.05'),(R_STAR,'orange',f'r*={R_STAR:.3f}'),
                     (0.4,'steelblue','r=0.4'),(np.median(rv),'red',f'IBM median r={np.median(rv):.3f}')]:
        T1v=T2r/(2*rv2)
        ax.plot(Nr, [n_qubit_ghz_min_cpsi(T1v,T2r,int(n)) for n in Nr], 'o-', color=c, lw=2, ms=5, label=l)
    ax.axhline(0.25, color='darkred', lw=2, ls='--'); ax.set_xlabel('N (GHZ)'); ax.set_ylabel('C·Ψ_min')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3); ax.set_ylim(0,1)

    ax=axes3[1]
    cN=[]
    for r2 in rv:
        T1v=T2r/(2*r2) if r2>0 else 1e10
        found=False
        for N in range(1,30):
            cm2=n_qubit_ghz_min_cpsi(T1v,T2r,N)
            if cm2 is not None and cm2<0.25: cN.append(N); found=True; break
        if not found: cN.append(30)
    cN=np.array(cN)
    ax.hist(cN, bins=range(1,32), color='purple', alpha=0.7, edgecolor='black', lw=0.5, density=True)
    ax.axvline(np.median(cN), color='red', lw=2, ls='--', label=f'Median N={np.median(cN):.0f}')
    ax.set_xlabel('Kritisches N'); ax.set_ylabel('Dichte'); ax.legend(fontsize=10); ax.grid(True, alpha=0.3)
    fig3.tight_layout()
    fig3.savefig(OUTPUT_DIR/'multiqubit_prediction.png', dpi=150, bbox_inches='tight')
    plt.close(fig3); print(f"  -> multiqubit_prediction.png")
    print(f"\nAlle Plots in: {OUTPUT_DIR}/")

def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--mode', choices=['collect','analyze','synthetic','full'], default='synthetic')
    p.add_argument('--days', type=int, default=HISTORY_DAYS)
    p.add_argument('--token', type=str, default=None)
    a = p.parse_args()
    tok = a.token or IBM_TOKEN
    if a.mode in ('collect','full'):
        if tok == "YOUR_TOKEN_HERE":
            print("No token! Running synthetic...\n")
            recs = generate_synthetic_data()
            save_csv(recs, "ibm_torino_synthetic.csv")
        else:
            recs = collect_historical_data(tok, a.days)
            save_csv(recs)
        d = analyze(recs)
        if d: plot_all(d)
    elif a.mode == 'analyze':
        recs = load_csv()
        d = analyze(recs)
        if d: plot_all(d)
    elif a.mode == 'synthetic':
        recs = generate_synthetic_data()
        save_csv(recs, "ibm_torino_synthetic.csv")
        d = analyze(recs)
        if d: plot_all(d)
    print("\nFertig!")

if __name__ == '__main__':
    main()
