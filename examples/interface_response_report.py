"""Measured local conservation progress and the failed response-refinement gate."""
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
from .implicit_plasma_records import provenance,freeze,sha


def run(destination):
    destination=Path(destination)
    if destination.exists():raise ValueError('Refusing overwrite')
    inputs={}
    def read(name):
        path=Path('outputs')/name/'summary.json';inputs[str(path)]=sha(path)
        data=json.loads(path.read_text())
        if not data['terminal']:raise ValueError('Complete measured input required')
        return data
    budgets=read('matched_interface_budgets')['cases'];fine=read('matched_interface_axial256')['cases']
    audits=[read('matched_mode_guarded_'+s+'_audit') for s in ('ne4_z32_t30','ne4_z32_t60','ne4_z32_t120','ne4_z64_t60')]
    history=read('matched_mode_response_history')
    comparisons=[read('matched_mode_time_comparison'),read('matched_mode_time120_comparison')]
    energies=read('matched_mode_signal_energy')['cases']+read('matched_mode_time120_signal_energy')['cases']+read('matched_mode_medium_signal_energy')['cases']
    continuum=read('matched_mode_extracted_continuum');probe=read('matched_mode_extracted_response')
    read('matched_mode_extracted_fluid_refinement');read('matched_mode_phase_audit')
    identity=provenance(['examples.interface_response_report']);destination.mkdir(parents=True);freeze(destination,identity)
    plt.rcParams.update({'font.size':9,'axes.titlesize':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(3,2,figsize=(14,12),layout='constrained')
    colors=['#386D99','#659B47','#8F74AD','#B74938'];names=['density','radial_velocity','radial_B','azimuthal_B','axial_B','current','electric','electron_pressure','ion_pressure']
    labels=['n','vᵣ','Bᵣ','Bθ','Bz','J','E','pₑ','pᵢ'];x=np.arange(len(names))
    ax=axes[0,0];selected=[budgets[0],budgets[1],budgets[3],fine[1]]
    for key,label in [('energy','Total energy'),('electron_entropy','Electron entropy'),('current','Current continuity'),('particles','Particles')]:
        ax.loglog([c['spec']['field_nz'] for c in selected],[c['quadrature'][1]['channels'][key]['cell_balance_relative_L2'] for c in selected],'-o',label=label,ms=4)
    ax.set(title='A  Local physical balances improve under refinement',xlabel='Axial cells / fluid points (radial mesh also refined)',ylabel='Relative cell-balance L2 defect')
    ax.set_xticks([32,64,128,256],labels=['32','64','128','256']);ax.xaxis.set_minor_locator(NullLocator());ax.legend(fontsize=8);ax.grid(alpha=.2,which='both')
    ax=axes[0,1]
    for audit,label,color in zip(audits,['Coarse: 30 steps','Coarse: 60 steps','Coarse: 120 steps','Larger grid: 60 steps'],colors):
        e=audit['quadrature'][1]['physical']['final']['errors']
        ax.semilogy(x,[e[n]['relative_L2'] for n in names],'-o',label=label,color=color,ms=4)
    ax.set(title='B  Full 100 μs response: current accuracy regresses',xticks=x,xticklabels=labels,ylabel='Relative physical-field L2 error')
    ax.legend(fontsize=8);ax.grid(alpha=.2,which='both')
    ax=axes[1,0]
    for index,label,color in [(1,'Coarse: 60 steps',colors[1]),(2,'Coarse: 120 steps',colors[2]),(3,'Larger grid: 60 steps',colors[3])]:
        case=history['cases'][index];path=Path('outputs/matched_mode_response_history')/case['arrays_file']
        if sha(path)!=case['arrays_sha256']:raise ValueError('Original history arrays required')
        inputs[str(path)]=sha(path);raw=np.load(path)
        ax.semilogy(raw['time']*1e6,raw['q16_diagnostics'][:,11,1],label=label,color=color)
    ax.set(title='C  Unwanted current grows on the larger grid',xlabel='Time (μs)',ylabel='Orthogonal current / initial reference norm')
    ax.legend(fontsize=8);ax.grid(alpha=.2,which='both')
    ax=axes[1,1];case=continuum['cases'][-1];power=np.array(case['radial_mode_power_fractions'])
    for k,label in [(0,'Density'),(8,'Axial B')]:
        ax.plot(np.arange(len(power)),power[:,k]*case['captured_direction_fractions'][k],'-o',label=label,ms=3)
    ax.axvline(16,color='#777777',ls=':',label='16 radial field cells')
    ax.set(title='D  Extracted response lies near the radial grid limit',xlabel='Independent Bessel radial-mode index',ylabel='Fraction of original component squared norm',xlim=(0,32))
    ax.legend(fontsize=8);ax.grid(alpha=.2)
    ax=axes[2,0]
    for j,case in enumerate(energies):
        drift=case['quadrature'][1]['signed_energy_drift_over_initial_quadratic']
        ax.semilogy([j-.09,j+.09],[abs(drift['minus']),abs(drift['plus'])],'o',color=colors[3] if j==4 else colors[0])
    ax.set(title='E  Energy drift measured against the small signal',xticks=np.arange(5),xticklabels=['30 steps','60 steps','Double ε','120 steps','Larger grid'],ylabel='|Actual total-energy drift| / initial quadratic scale')
    ax.grid(alpha=.2,which='both')
    ax=axes[2,1]
    for case,label,color in zip(comparisons,['60 versus 30 steps','120 versus 60 steps'],colors[:2]):
        errors=case['quadrature'][1]['differences']
        ax.semilogy(x,[errors[n]['relative_L2'] for n in names],'-o',label=label,color=color,ms=4)
    ax.set(title='F  Time controls do not yet show fourth-order convergence',xticks=x,xticklabels=labels,ylabel='Difference / earlier-resolution response norm')
    ax.legend(fontsize=8);ax.grid(alpha=.2,which='both')
    fig.suptitle('FRC sandbox verification: local-balance progress, open spatial/time accuracy gate\nFilled periodic test cylinder; physical masses; no formation or laboratory stability claim',fontsize=14)
    fig.savefig(destination/'verification.png',dpi=170);fig.savefig(destination/'verification.pdf');plt.close(fig)
    result=dict(terminal=True,source_sha256=identity['source_sha256'],input_sha256=inputs,
        quadrature=20,history_quadrature=16,local_budget_cases=8,new_complete_signed_trajectories=4,
        medium_current_relative_error=audits[-1]['quadrature'][1]['physical']['final']['errors']['current']['relative_L2'],
        extracted_current_growth=probe['quadrature'][1]['channels']['current']['inferred_growth'],
        strong_ampere_relative=continuum['cases'][-1]['strong_ampere_relative'],
        limitations=['Full growing-response mesh and time convergence not established','Quadratic energy is only an initial normalization, not a separate invariant',
            'Projected continuum comparison retains uncaptured radial norm','No heated moving vacuum, transport, ends/coils or formation adoption'],
        files={name:sha(destination/name) for name in ('verification.png','verification.pdf')})
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
