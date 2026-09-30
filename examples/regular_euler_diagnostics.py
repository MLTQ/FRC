"""Matched physical comparison of Gauss-r and regular-coordinate fluid states."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
from .mixed_euler_study import initialize as old_initialize
from .regular_euler_study import initialize as new_initialize
from .mixed_euler_sampling import sample as old_sample,grid,differences
from .regular_euler_sampling import sample as new_sample
from .cylindrical_entropy_reference import physical
from .radial_work_records import sha


def collect(root,regular):
    initializer,sampler=(new_initialize,new_sample) if regular else (old_initialize,old_sample)
    records=[json.loads((root/f'fluid_{i:02d}.json').read_text()) for i in range(14)];states={};profiles=[];representations=[]
    for i,record in enumerate(records):
        ref,d,_=initializer(record['spec']);a=np.load(root/record['arrays_file']);r,z,volume=grid(ref.radius,2*ref.half_length,32,256,24)
        states[i]=sampler(d,a['state'],r,len(z));f=physical(states[i],r[:,None],ref)
        if i in (4,5,6,7):
            profiles.append(dict(case=i,elements=d.elements,maximum_velocity_over_sound=float(np.max(abs(f['velocity']))/ref.sound),minimum_density_over_background=float(np.min(f['n'])/ref.number_density)))
        if record['spec']['manufactured']:
            exact=ref.state(r[:,None],z[None,:],record['duration'])
            target=(ref.regular_state(d.x,d.z,record['duration']) if regular else ref.state(d.r,d.z,record['duration']))
            if regular:target[0,0,:,1:3]=0.
            represented=sampler(d,target,r,len(z))
            representations.append(dict(case=i,exact_state_interpolation=differences(represented,exact,r,volume,ref,ref.radius/8),
                evolved_difference_from_exact_nodal_interpolant=differences(states[i],represented,r,volume,ref,ref.radius/8)))
    comparisons=[]
    for a,b in ((4,5),(5,6),(6,7),(5,10),(10,11),(1,8),(8,9)):
        comparisons.append(dict(cases=[a,b],errors=differences(states[a],states[b],r,volume,ref,ref.radius/8)))
    return dict(records=records,comparisons=comparisons,profiles=profiles,representation_split=representations)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('old',type=Path);parser.add_argument('new',type=Path);parser.add_argument('output',type=Path);args=parser.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    old=collect(args.old,False);new=collect(args.new,True)
    summary=dict(old_source=json.loads((args.old/'provenance.json').read_text())['source_sha256'],new_source=json.loads((args.new/'provenance.json').read_text())['source_sha256'],old=old,new=new)
    (out/'comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold'})
    fig,axes=plt.subplots(2,3,figsize=(14,8.4),layout='constrained');x=np.array([2,4,8,16]);colors=['#c3563c','#157f8e']
    for data,label,color in ((old,'Previous Gauss-r representation',colors[0]),(new,'Regular-coordinate formulation',colors[1])):
        axes[0,0].loglog(x,[abs(r['nodal_continuous_energy_difference']) for r in data['records'][4:8]],'-o',label=label,color=color)
        axes[0,1].loglog(x,[r['maximum_velocity_over_sound'] for r in data['profiles']],'-o',label=label,color=color)
        axes[0,2].loglog(x,[r['manufactured_errors']['n']['relative_L2'] for r in data['records'][:4]],'-o',label=label,color=color)
        axes[1,0].loglog(x,[r['manufactured_errors']['velocity']['axis_relative_L2'] for r in data['records'][:4]],'-o',label=label,color=color)
        axes[1,1].loglog(x[:-1],[r['errors']['velocity']['relative_L2'] for r in data['comparisons'][:3]],'-o',label=label,color=color)
        chosen=[data['records'][i] for i in (5,10,11)];axes[1,2].loglog([1,.5,.25],[abs(r['physical_energy_work_remainder']) for r in chosen],'-o',label=label,color=color)
    axes[0,0].set(title='Actual continuous-energy discrepancy',ylabel='Absolute grid / continuous energy difference',xlabel='Radial elements')
    axes[0,1].set(title='Between-point velocity spikes removed',ylabel='Maximum physical component speed / sound speed',xlabel='Radial elements')
    axes[0,2].set(title='Exact mixed manufactured solution',ylabel='Continuous density relative L2 error',xlabel='Radial elements (axial points = 8 × elements)')
    axes[1,0].set(title='Fixed central region r < R/8',ylabel='RMS velocity error / sound speed',xlabel='Radial elements')
    axes[1,1].set(title='Unforced mixed-flow mesh convergence',ylabel='RMS velocity difference / sound speed',xlabel='Coarse elements (compared with twice as many)')
    axes[1,2].set(title='Unforced time-step control at 4 elements',ylabel='Absolute relative grid-energy drift',xlabel='Time step / each formulation’s base step')
    for ax in axes.flat:ax.grid(True,which='major',alpha=.2)
    for ax in [*axes[0],axes[1,0]]:ax.set_xticks(x,[str(v) for v in x]);ax.xaxis.set_minor_locator(NullLocator())
    axes[1,1].set_xticks(x[:-1],[str(v) for v in x[:-1]]);axes[1,1].xaxis.set_minor_locator(NullLocator())
    axes[1,2].set_xticks([.25,.5,1],['0.25','0.5','1']);axes[1,2].xaxis.set_minor_locator(NullLocator())
    axes[0,0].legend(fontsize=8,frameon=False)
    fig.suptitle('Mixed cylindrical flow: resolving the physical axis representation',fontsize=18,fontweight='bold')
    fig.supxlabel('Same physical initial fields and duration • no numerical atmosphere or energy repair • magnetic coupling remains open',fontsize=10)
    fig.savefig(out/'report.png',dpi=180);fig.savefig(out/'report.pdf');plt.close(fig)
    (out/'report_source.json').write_text(json.dumps(dict(script_sha256=sha(__file__),old_audit_sha256=sha(args.old/'audit.json'),new_audit_sha256=sha(args.new/'audit.json')),indent=2)+'\n')
    print('Saved matched physical comparison and scientific figure')


if __name__=='__main__':main()
