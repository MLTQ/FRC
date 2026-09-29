"""Scientific figures exposing coupled pressure errors and the slow-mode cause."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def plot(report,directory):
    rates={row['configuration']['name']:row for row in report['rate_cases']}
    runs={row['configuration']['name']:row for row in report['evolution_cases']}
    fig,axes=plt.subplots(2,2,figsize=(11.5,8.4),layout='constrained')
    ax=axes[0,0]
    for prefix,label,nrs in (('moving_','Manufactured, n = 10¹⁸',(8,16,32,64,128)),
                              ('target_','Manufactured, n = 10²⁰',(8,16,32,64)),
                              ('wave_0_','Slow wave',(8,16,32,64)),('wave_2_','Fast wave',(8,16,32))):
        ax.loglog(nrs,[rates[prefix+str(n)]['metrics']['ion_pressure']['relative_L2'] for n in nrs],'o-',label=label)
    ax.axhline(1,color='.5',ls=':',lw=1);ax.set_title('Instantaneous ion-pressure rate error')
    ax.set_ylabel('Error / physical pressure-rate signal');ax.legend(fontsize=8)
    ax=axes[0,1];nrs=(8,16,32)
    for prefix,color,label in (('moving_','C0','n = 10¹⁸'),('target_','C1','n = 10²⁰')):
        for metric,style,detail in (('relative_L2','o-','full pressure'),('error_in_change_relative_L2','s--','physical change')):
            ax.loglog(nrs,[runs[prefix+str(n)]['metrics']['ion_pressure'][metric] for n in nrs],style,color=color,label=label+', '+detail)
    ax.axhline(1,color='.5',ls=':',lw=1);ax.set_title('Manufactured evolution: two normalizations')
    ax.set_ylabel('Ion-pressure error / stated signal');ax.legend(fontsize=8)
    ax=axes[1,0]
    for key,label in (('magnetic','Magnetic field'),('electron_pressure','Electron pressure'),('ion_pressure','Ion pressure')):
        ax.loglog(nrs,[runs['wave_'+str(n)]['metrics'][key]['relative_L2'] for n in nrs],'o-',label=label)
    ax.axhline(1,color='.5',ls=':',lw=1);ax.set_title('Unforced fast wave, one quarter period')
    ax.set_ylabel('Error / small oscillation amplitude');ax.legend(fontsize=8)
    ax=axes[1,1];diagnosis=[row for row in report['flux_diagnosis']['cases'] if row['number_density_m3']==1e18]
    for key,label,style in (('central','Current central-flux error','o-'),('diffusive','Diffusive contribution','s-'),
                             ('cubic_linear','Candidate on linear momentum','^-'),('cubic_finite_amplitude','Candidate on finite-amplitude state','x--')):
        ax.loglog([row['nr'] for row in diagnosis],[row['mass_rate_relative_L2'][key] for row in diagnosis],style,label=label)
    ax.set_title('Slow-mode continuity: reconstruction dominates')
    ax.set_ylabel('Rate error or contribution / physical rate');ax.legend(fontsize=8)
    for ax in axes.flat:
        ax.set_xlabel('Radial cells (axial cells = 2 × radial cells)');ax.grid(True,which='both',alpha=.22)
        ax.set_xticks([8,16,32,64,128],labels=['8','16','32','64','128'])
    fig.suptitle('Finite-mass moving reference: conservation passes, pressure accuracy remains limited',fontsize=12)
    path=directory/'accuracy.png';fig.savefig(path,dpi=170);plt.close(fig)
    return [path]
