"""Render the validated branch and its two cross-parity coincidences."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]


def main():
    first=json.loads((ROOT/'simulations/results/route_b_n4_range_even.json').read_text())['paths']['0_-1']['rows']
    second=json.loads((ROOT/'simulations/results/route_b_n4_fold_path.json').read_text())
    from route_b_n4_bridge_certificate import audit_connection
    from flint import ctx
    ctx.prec=384
    bridge=json.loads((ROOT/'simulations/results/route_b_n4_bridge_certificate.json').read_text())
    assert bridge['success'] and audit_connection(bridge['tiles'])
    assert not audit_connection(bridge['tiles'],wrong_germ=True)
    from route_b_n4_incoming_certificate import audit,close_crossings
    incoming=json.loads((ROOT/'simulations/results/route_b_n4_incoming_certificate.json').read_text())
    assert incoming['success'] and audit(incoming['tiles'],incoming['seed_anchor_dyadic'],incoming['fold_anchor_dyadic'])
    close_crossings(incoming)
    crossings=json.loads((ROOT/'simulations/results/route_b_n4_incoming_crossing.json').read_text())['roots']
    assert all(r['jordan_block_sizes']==[2,1] for r in crossings)
    finite=second['rows'];tail=second['compactification']['rows']
    q0=np.sqrt((np.sqrt(13)-1)/6);ec=np.sqrt(2)-2;ei=np.sqrt(3)-2
    fig,(ax,bx)=plt.subplots(1,2,figsize=(11.5,4.8),gridspec_kw={'width_ratios':[1.3,1]})
    ax.plot([1/q0]+[1/float(r['q']) for r in first]+[.5],[0]+[float(r['epsilon']) for r in first]+[ec],color='#167a9b',label='Seed to turn (sampled curve)')
    ax.plot([.5]+[1/float(r['q']) for r in finite]+[float(r['z']) for r in tail],[ec]+[float(r['epsilon']) for r in finite]+[float(r['epsilon']) for r in tail],color='#b75b18',label='Second arm (sampled curve)')
    ax.axvspan(0,1/q0,color='#287544',alpha=.10,label='Validated parameter range')
    ax.scatter([1/float(r['center'][1]) for r in crossings],[float(r['center'][2]) for r in crossings],
               c='#943946',marker='D',s=38,zorder=6,label='Cross-parity coincidence: J2 + J1')
    ax.scatter([1/q0,.5,0],[0,ec,ei],c=['#167a9b','#423c83','#287544'],marker='o',zorder=5)
    ax.annotate('Diabolic seed\nq = 0.658983',(1/q0,0),xytext=(-80,-40),textcoords='offset points')
    ax.annotate('Exact EP2 turn\nq = 2',(.5,ec),xytext=(15,17),textcoords='offset points')
    ax.annotate('Exact limit\nepsilon = sqrt(3) - 2',(0,ei),xytext=(13,3),textcoords='offset points')
    ax.set_xlabel('Inverse coupling z = 1/q (infinity at the left)');ax.set_ylabel('Equal-end perturbation epsilon')
    ax.set_title('One path through the turn');ax.legend(loc='upper left',fontsize=8)
    qq=[float(r['q']) for r in finite]+[float(r['q']) for r in tail]
    vv=[float(r['omega'])-float(r['q']) for r in finite]+[float(r['v']) for r in tail]
    bx.plot(qq,vv,color='#b75b18',label='Tracked frequency offset')
    bx.axhline(-np.sqrt(5)/16,color='#287544',ls='--',label='Selected limit: -sqrt(5)/16')
    bx.axhline(-np.sqrt(2)/4,color='gray',ls=':',label='Other negative germ: -sqrt(2)/4')
    bx.axvspan(2,max(qq),color='#287544',alpha=.10,label='Validated arm: q >= 2')
    bx.set_xscale('log');bx.set_xlabel('Coupling q');bx.set_ylabel('Frequency offset omega - q')
    bx.set_title('Identifying the germ at infinity');bx.legend(fontsize=8,loc='lower right')
    for a in (ax,bx):a.grid(alpha=.2)
    fig.suptitle('N=4: a connected reflection-even EP branch from the diabolic seed to infinity',fontsize=13)
    fig.text(.5,.025,f'Validated by {incoming["tile_count"]} incoming + {bridge["tile_count"]} outgoing tiles. Lines: numerical samples.\n'
             'Full-block character: J2 + J1 at the diamonds, EP2 elsewhere on the open branch; the seed is diabolic.',ha='center',fontsize=8.5)
    fig.tight_layout(rect=(0,.085,1,.94))
    for ext in ['png','svg']:fig.savefig(ROOT/f'visualizations/route_b_n4_second_arm.{ext}',dpi=170)


if __name__=='__main__':main()
