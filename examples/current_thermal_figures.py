"""Scientific current/thermal convergence and retained-limitation figures."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.constants import electron_mass as ME
from .current_thermal_norms import physical_energy
from .toroidal_thermal_symbol import spectrum
from .electron_reference import ElectronReference


def figures(directory,records,archives,meshes,comparisons,symbol):
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.23,'figure.dpi':150})
    fig,axes=plt.subplots(2,3,figsize=(14,8.5),constrained_layout=True)
    ax=axes[0,0];sizes=(8,16,32,64,128)
    for key,label in (('local_pressure_defect_rms_pa','volume RMS'),('local_pressure_defect_max_pa','maximum'),('axis_pressure_defect_max_pa','axis maximum')):
        ax.loglog(sizes,[records[f'original_{n}']['metrics'][key] for n in sizes],'o-',label=label)
    ax.set(title='Original fixture: local energy mismatch',xlabel='radial cells',ylabel='ion-pressure equivalent [Pa]');ax.legend()
    ax=axes[0,1]
    for prefix,label in (('wave_h1','physical mass, n=10¹⁸'),('wave_target','physical mass, n=10²⁰'),('wave_strong','strong wave, n=10²⁰')):
        ax.loglog((8,16,32),[records[f'{prefix}_{n}']['metrics']['magnetic_relative_L2'] for n in (8,16,32)],'o-',label=label)
    ax.set(title='Exact nonlinear wave: physical B error',xlabel='radial cells',ylabel='L2 / magnetic perturbation');ax.legend()
    ax=axes[0,2]
    for key,label in (('current_relative_L2','current volume L2'),('axis_current_relative_L2','current first annulus')):
        ax.loglog((8,16,32),[records[f'wave_h1_{n}']['metrics'][key] for n in (8,16,32)],'o-',label=label)
    ax.set(title='Actual current error is independently checked',xlabel='radial cells',ylabel='relative L2 error');ax.legend()
    ax=axes[1,0]
    for prefix,label in (('wave_h1','base'),('wave_inertial_strong','strong inertia fixture')):
        for key,style in (('radial','o-'),('axial','s--')):
            ax.loglog((8,16,32),[records[f'{prefix}_{n}']['metrics'][key+'_power_relative_L2'] for n in (8,16,32)],style,label=label+', '+key)
    ax.set(title='Full continuum face energy flux',xlabel='radial cells',ylabel='relative face-power L2 error');ax.legend(fontsize=8)
    ax=axes[1,1]
    for family in ('original','wave'):
        values=comparisons[family+'_time']['errors'];steps=sorted(map(int,values))
        ax.loglog(steps,[values[str(n)]['magnetic']['rms_difference'] for n in steps],'o-',label=family)
    ax.set(title='Fixed mesh, separate time refinement',xlabel='steps to common final time',ylabel='physical B difference from 320 steps [T]');ax.legend()
    ax=axes[1,2]
    for name,label in (('original_128','physical mass, 80 steps'),('original_fine_half_time','physical mass, 160 steps'),('massless_original_128','massless joint, 80 steps')):
        record=records[name];h=record['history_columns']
        ax.semilogy(np.asarray(h['time_s'])[1:]*1e6,h['local_pressure_defect_max_pa'][1:],label=label)
        if not record['completed']:ax.plot(h['time_s'][-1]*1e6,h['local_pressure_defect_max_pa'][-1],'rx',ms=9)
    ax.set(title='Maximum local defect at nr=128',xlabel='accepted time [µs]',ylabel='pressure-equivalent mismatch [Pa]');ax.legend(fontsize=8)
    fig.suptitle('Constrained fixed-density current + two temperatures: measured accuracy, not formation',fontsize=14)
    one=directory/'current_thermal_convergence.png';fig.savefig(one);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(14,4.8),constrained_layout=True)
    ax=axes[0]
    for name,label in (('original_128','physical mass, 80 steps'),('original_fine_half_time','physical mass, 160 steps'),('massless_original_128','massless joint, 80 steps')):
        rec=records[name];h=rec['history_columns'];ax.plot(np.asarray(h['time_s'])*1e6,h['min_inferred_ion_pressure_pa'],label=label)
        if not rec['completed']:
            f=rec['failure']
            if 'diagnostic' in f:ax.plot(f['trial_time_s']*1e6,f['diagnostic']['min_inferred_ion_pressure_pa'],'rx',ms=9)
    ax.axhline(0,color='black',lw=.8);ax.set(title='Independent conservative ion partition',xlabel='time [µs]',ylabel='minimum inferred ion pressure [Pa]');ax.legend(fontsize=8)
    ax=axes[1];name='original_128';a=archives[name];mesh=meshes[name];cfg=records[name]['configuration'];mu=ME*cfg['ion_mass']/(ME+cfg['ion_mass'])
    energies=physical_energy(mesh,a['final_psi'],a['final_g'],a['final_current'],a['number_density'],a['final_electron_entropy'],a['final_ion_entropy'],reduced_mass=mu,order=8)
    defect=(cfg['gamma_i']-1)*(a['final_total_energy']-sum(energies.values()));limit=abs(defect).max()
    plot=ax.pcolormesh(mesh.zf*1e3,mesh.rf*1e3,defect,cmap='RdBu_r',vmin=-limit,vmax=limit,shading='flat',rasterized=True)
    ax.grid(False);ax.set(title='Physical mass, nr=128, t=1 µs',xlabel='z [mm]',ylabel='r [mm]');fig.colorbar(plot,ax=ax,label='local energy mismatch [Pa equivalent]')
    ax=axes[2];point=symbol['moving_maximum'];angle=np.linspace(0,np.pi,361)
    values=spectrum(ElectronReference(),point['r_m'],point['z_m'],angle)
    for label,value in (('pressure only',values['pressure_growth_per_k']),('joint, fixed n',values['joint_growth_per_k']),('moving massless',np.max(abs(values['moving_eigenvalues'].imag),axis=-1))):
        ax.plot(angle*180/np.pi,value,label=label)
    ax.set(title='Massless toroidal frozen-symbol warning',xlabel='wavevector angle from r [degrees]',ylabel='growth rate / |k| [m/s]');ax.legend(fontsize=8)
    fig.suptitle('Local energy and continuum-model limits remain explicit',fontsize=14)
    two=directory/'current_thermal_limits.png';fig.savefig(two);plt.close(fig)
    return [one,two]
