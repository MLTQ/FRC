"""Independent species-rate, physical conservation and convergence report."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .two_entropy_reference import physical,parent_rate
from .two_entropy_study import initialize
from .two_entropy_records import sha


def species_rates(u,rate,ref):
    f=physical(u,ref);gamma=np.array([ref.gamma_e,ref.gamma_i]);nt=rate[:,0]
    st=(rate[:,4:]-f['entropy']*nt[:,None])/f['n'][:,None]
    pressure=f['pressure']*(st+gamma*nt[:,None]/f['n'][:,None])
    temperature=f['temperature']*(st+(gamma-1)*nt[:,None]/f['n'][:,None])
    return pressure,temperature


def main():
    root=Path(__file__).resolve().parents[1];directory=root/'outputs/two_entropy_verification';output=root/'outputs/two_entropy_diagnostics';output.mkdir(exist_ok=True)
    rates=json.loads((directory/'rates.json').read_text());trajectories=json.loads((directory/'trajectories.json').read_text());diagnostics=[]
    for record in rates:
        path=directory/record['arrays_file'];assert sha(path)==record['arrays_sha256'];ref,solver,_=initialize(record['spec'])
        with np.load(path) as a:
            actual=species_rates(a['state'],a['rate'],ref);exact=species_rates(a['state'],parent_rate(ref,solver.z,True)[0],ref)
            errors={f'{key}_{species}':float(np.linalg.norm(value[:,j]-target[:,j])/np.linalg.norm(target[:,j]))
                    for key,value,target in zip(('pressure','temperature'),actual,exact) for j,species in enumerate(('electron','ion'))}
        diagnostics.append(dict(spec=record['spec'],species_rate_errors=errors))
    runs=[]
    for record in trajectories:
        path=directory/record['arrays_file'];assert sha(path)==record['arrays_sha256'];ref,solver,_=initialize(record['spec'])
        with np.load(path) as a:
            history=a['history'];scale=solver.dx*np.sum(abs(a['initial']),axis=0)
            drift=np.max(abs(history[:,1:7]-history[0,1:7]),axis=0)/np.maximum(scale,ref.number_density*solver.length*1e-15)
            f=physical(a['final'],ref);v=ref.fields(np.zeros_like(solver.z),solver.z,record['end_time'])
            species_errors={key:float(np.linalg.norm(f[field][:,j]-v[label])/np.linalg.norm(v[label]))
                for field,base,labels in [('pressure','p',('pe','pi')),('temperature','T',('Te','Ti'))]
                for j,(label,key) in enumerate(zip(labels,(base+'_electron',base+'_ion')))}
            entropy_error=[float(np.max(abs(f['entropy'][:,j]-v[key]))) for j,key in enumerate(('electron','ion'))]
        runs.append(dict(file=path.name,spec=record['spec'],maximum_normalized_state_drift=drift.tolist(),
                         nodal_species_errors=species_errors,maximum_specific_entropy_errors=entropy_error))
    (output/'measurements.json').write_text(json.dumps(dict(rates=diagnostics,trajectories=runs),indent=2)+'\n')
    plt.rcParams.update({'font.size':11,'axes.titlesize':13,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12,9));fig.subplots_adjust(left=.10,right=.97,top=.82,bottom=.10,hspace=.43,wspace=.32)
    fig.suptitle('Joint transport removes the spatial energy defect',x=.10,y=.955,ha='left',fontsize=20,weight='bold')
    fig.text(.10,.894,'Periodic axial verification • physical two-species entropy • no energy reset',fontsize=12,color='#44546a')
    for order,color in [(2,'#9297a5'),(4,'#9674aa'),(6,'#126d9c')]:
        selected=[r for r in rates if r['spec']['order']==order]
        axes[0,0].loglog([r['spec']['points'] for r in selected],[r['energy_rate_relative'] for r in selected],'o-',label=f'Order {order}',color=color)
    axes[0,0].set(title='Original Euler physical-energy rate',ylabel='Relative L2 error');axes[0,0].legend()
    for key,label,color in [('pressure_electron','Electron pressure','#126d9c'),('pressure_ion','Ion pressure','#be8530')]:
        selected=[r for r in diagnostics if r['spec']['order']==6]
        axes[0,1].loglog([r['spec']['points'] for r in selected],[r['species_rate_errors'][key] for r in selected],'o-',label=label,color=color)
    axes[0,1].set(title='Both species retain high-order accuracy',ylabel='Relative L2 rate error');axes[0,1].legend()
    for kind,label,color in [('contact','Full contact transit','#126d9c'),('acoustic','Nonlinear acoustic wave','#2b8d74')]:
        selected=[r for r in trajectories if r['spec']['kind']==kind and 'step_factor' not in r['spec'] and 'direction' not in r['spec']]
        axes[1,0].loglog([r['spec']['points'] for r in selected],[r['norms'][-1]['n'] for r in selected],'o-',label=label,color=color)
    axes[1,0].set(title='Independent continuous density errors',ylabel='Relative L2 error');axes[1,0].legend()
    selected=[trajectories[i] for i in (1,8,9)];factor=[1,2,4]
    axes[1,1].loglog(factor,[abs(r['nodal_energy_drift']) for r in selected],'o-',label='Nodal physical energy',color='#126d9c')
    axes[1,1].loglog(factor,[abs(r['continuous_energy_drift']) for r in selected],'x--',label='Independent continuous energy',color='#be8530')
    axes[1,1].set(title='Fixed-grid contact: halve the time step',ylabel='Absolute relative energy change',xlabel='Time-step refinement factor');axes[1,1].legend(fontsize=9)
    for i,a in enumerate(axes.flat):
        a.grid(True,alpha=.2);a.xaxis.set_minor_locator(matplotlib.ticker.NullLocator());a.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        a.set_xticks([16,32,64,128] if i<3 else factor)
        if i<3:a.set_xlabel('Periodic axial points')
    fig.text(.10,.025,'Nodal spatial conservation is exact to roundoff; RK4 time error remains. Cylindrical, magnetic/current and physical-edge coupling are still required.',fontsize=9,color='#44546a')
    fig.savefig(output/'report.png',dpi=190);plt.close(fig)
    source=Path(__file__);(output/'report_source.py').write_bytes(source.read_bytes());(output/'report_source.md').write_bytes(source.with_suffix('.md').read_bytes())


if __name__=='__main__':main()
