"""Scientific diagnostic figure from preserved baseline and radial-work evidence."""
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator


def main():
    root=Path('outputs');dest=root/'log_fluid_diagnostics';dest.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.18,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,3,figsize=(17,10));blue='#176b99';orange='#ce6622';green='#16876b'
    for folder,name,color in [('log_fluid_pilot','Original pressure stencil',orange),('radial_work_time_controls','Radial work pair, half step',blue)]:
        index=0 if folder=='log_fluid_pilot' else 14
        r=json.loads((root/folder/f'fluid_{index:02d}.json').read_text());a=np.load(root/folder/r['arrays_file']);h=a['history']
        axes[0,0].semilogy(h[:,0]/r['end_time'],np.maximum(abs(h[:,11]/h[0,11]-1),1e-16),label=name,color=color)
    axes[0,0].set(xlabel='Fraction of full contact transit',ylabel='Peak speed departure / initial speed',title='Long contact: failure and a completed control')
    axes[0,0].legend(fontsize=8);axes[0,0].annotate('Original run rejects at 0.24 transit',(.24,10),(.30,2),arrowprops={'arrowstyle':'->'},fontsize=8)
    a=np.load(root/'log_fluid_radial_diagnosis/radial_05.npz')
    r=json.loads((root/'log_fluid_radial_diagnosis/records.json').read_text())['records'][5];scale=.035/r['sound']
    for key,label,color,marker in [('eigenvalues','Original pressure',orange,'x'),('adjoint_eigenvalues','Energy-adjoint linear control',blue,'.')]:
        v=a[key]*scale;axes[0,1].scatter(v.real,v.imag,s=17,label=label,color=color,marker=marker)
    axes[0,1].axvline(0,color='gray',lw=.7);axes[0,1].set(xlabel='Growth rate × R/c',ylabel='Frequency × R/c',title='Radial acoustic spectrum · nr32')
    axes[0,1].legend(fontsize=8)
    r=json.loads((root/'radial_work_time_controls/fluid_14.json').read_text());a=np.load(root/'radial_work_time_controls/fluid_14.npz');h=a['history'];e0=h[0,5]
    axes[0,2].plot(h[:,0]/r['end_time'],(h[:,5]-e0)/e0*1e9,label='Physical energy change',color=blue)
    axes[0,2].plot(h[:,0]/r['end_time'],h[:,12]/e0*1e9,'--',label='Integrated spatial defect',color=orange)
    axes[0,2].plot(h[:,0]/r['end_time'],(h[:,4]-h[0,4])/e0*1e9,label='Conservative energy account',color=green)
    axes[0,2].set(xlabel='Fraction of full contact transit',ylabel='Relative energy change × 10⁹',title='Energy accounting remains explicit · nr8')
    axes[0,2].legend(fontsize=8)
    rates=json.loads((root/'radial_work_rates/rates.json').read_text());audit=json.loads((root/'radial_work_rates/rate_audit.json').read_text())['records']
    n=np.array([8,16,32,64]);whole=[r['norms']['velocity']['relative_L2'] for r in rates[4:8]]
    axes[1,0].loglog(n,whole,'o-',label='Whole-volume acceleration',color=blue)
    for key,color in [('radial',orange),('axial',green)]:
        axes[1,0].loglog(n,[r['velocity_components'][key]['axis_relative_L2'] for r in audit[4:8]],'s--',label='First-annulus '+key,color=color)
    axes[1,0].set(xlabel='Radial cells',ylabel='Relative physical L₂ error',title='Remaining failure: local radial force accuracy')
    axes[1,0].set_xticks(n,labels=n);axes[1,0].xaxis.set_minor_locator(NullLocator());axes[1,0].legend(fontsize=8)
    trans=json.loads((root/'log_fluid_translation_probe/records.json').read_text())['records']
    for key,color in [('electron',blue),('ion',green)]:
        axes[1,1].loglog(n,[r['species'][key]['relative_chain_rule_defect'] for r in trans[:4]],'o-',label=key.capitalize(),color=color)
    axes[1,1].loglog(n,6e-7*(8/n)**6,':',label='Sixth-order guide',color='gray')
    axes[1,1].set(xlabel='Radial cells',ylabel='Thermal defect / uncancelled work',title='Uniform translation: force-only obstruction')
    axes[1,1].set_xticks(n,labels=n);axes[1,1].xaxis.set_minor_locator(NullLocator());axes[1,1].legend(fontsize=8)
    from examples.log_fluid_reference import LogFluidReference
    from examples.log_fluid_norms import FluidSampler
    from plasma_bench.mhd_mesh import MHDMesh
    from plasma_bench.log_fluid import LogFluidState
    from examples.radial_work_study import STATE_KEYS
    r=json.loads((root/'radial_work_pilot/fluid_03.json').read_text());a=np.load(root/'radial_work_pilot/fluid_03.npz')
    ref=LogFluidReference(**r['reference']);grid=FluidSampler(MHDMesh(ref.radius,ref.half_length,8,16),5,14)
    state=LogFluidState(*(a[k+'_final'] for k in STATE_KEYS));f=grid.physical(state,ref);exact=ref.fields(grid.r,grid.z,r['end_time']);row=4*14
    axes[1,2].plot(grid.z[0]*1000,(exact['n'][row]/ref.number_density-1)*100,color='black',label='Independent characteristic solution')
    axes[1,2].plot(grid.z[0]*1000,(f['n'][row]/ref.number_density-1)*100,'--',color=blue,label='Self-consistent nonlinear wave')
    axes[1,2].set(xlabel='Axial position (mm)',ylabel='Density perturbation (%)',title='Unforced finite-amplitude wave · nr8')
    axes[1,2].legend(fontsize=8)
    fig.suptitle('Moving-fluid diagnostics: progress, preserved failures and the next accuracy barrier',fontsize=17,y=.98)
    fig.text(.5,.018,'Bounded axisymmetric hydrodynamic verification · no field/current evolution, physical edge or formation claim · experimental radial-work variant',ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.045,1,.95));fig.savefig(dest/'report.png',dpi=160);plt.close(fig)


if __name__=='__main__':main()
