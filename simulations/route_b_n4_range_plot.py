"""Render tracked N4 real EP paths and independently certified boundaries."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from flint import arb

ROOT=Path(__file__).resolve().parents[1]
load=lambda name:json.loads((ROOT/f'simulations/results/route_b_n4_range_{name}.json').read_text())


def main():
    fig,axes=plt.subplots(1,3,figsize=(13,4.5))
    colors=['#167a9b','#b75b18']
    one=load('one');even=load('even');odd=load('precise')
    for track in one['tracks']:
        rows=track['points'];axes[0].plot([float(r['epsilon']) for r in rows],[np.sqrt(float(r['y'])) for r in rows],color=colors[track['branch']])
    for track in even['paths'].values():
        rows=track['rows'];axes[1].plot([float(r['epsilon']) for r in rows],[float(r['q']) for r in rows],color=colors[track['branch']])
    for track in odd.values():
        rows=track['rows']
        for sign in (-1,1):axes[2].plot([sign*float(r['epsilon']) for r in rows],[float(r['q']) for r in rows],color=colors[track['branch']])
    for row in load('certificate')['roots']:
        ax=axes[['one','even','odd'].index(row['profile'])]
        e=float(arb(row['epsilon']));q=float(arb(row['q']))
        for sign in ((-1,1) if row['profile']=='odd' else (1,)):
            ax.scatter(sign*e,q,marker='D' if row['kind']=='triple' else 's',s=45,color='#aa2233' if row['kind']=='triple' else '#423c83',zorder=5)
    for ax,title in zip(axes,['One end','Equal ends','Opposite ends']):
        ax.set_title(title);ax.set_xscale('symlog',linthresh=.02)
        ax.set_xlabel('End perturbation epsilon (symlog)');ax.set_ylabel('Real coupling q')
        ax.axhline(np.sqrt((np.sqrt(13)-1)/6),color='gray',lw=.6,ls=':')
        ax.axvline(0,color='gray',lw=.6);ax.grid(alpha=.2)
    axes[0].set_xticks([-.01,0,.01,.1,1],labels=['-0.01','0','0.01','0.1','1'])
    fig.suptitle('N=4 self-mirrored EP paths: two different boundaries',fontsize=14)
    legend=[Line2D([0],[0],color=c,label=f'Initial coefficient {s}') for c,s in zip(colors,['negative','positive'])]
    legend += [Line2D([0],[0],color=c,marker=m,ls='',label=l) for c,m,l in [('#aa2233','D','Certified EP3'),('#423c83','s','Certified EP2 turn')]]
    fig.legend(handles=legend,loc='lower center',ncol=4,frameon=False)
    fig.text(.5,.10,'Lines are numerical continuation; boundary certificates do not certify the whole connecting interval.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.14,1,.92))
    for ext in ['png','svg']:fig.savefig(ROOT/f'visualizations/route_b_n4_range.{ext}',dpi=170)


if __name__=='__main__':main()
