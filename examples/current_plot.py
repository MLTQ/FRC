"""Scientific plots of current-component accuracy and conditional physical scales."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from .inertial_screening_reference import InertialScreening


def accuracy(directory,rows,arrays):
    fig,axes=plt.subplots(2,3,figsize=(16,9),layout='constrained')
    keys=('field_relative_volume_L2','current_relative_volume_L2','local_heat_relative_volume_L2')
    for de in (0.,.05,.2,1.):
        selected=[r for r in rows if r['name'] in [f'decay_d{de:g}_{n}' for n in (8,16,32,64)]]
        for ax,key in zip(axes[0],keys):
            ax.loglog([r['configuration']['radial_cells'] for r in selected],[r['final_errors'][key] for r in selected],'o-',label=f'de/R={de:g}')
    for de in (.05,.2,1.,4.):
        selected=[r for r in rows if r['name'] in [f'screen_d{de:g}_{n}' for n in (8,16,32,64)]]
        axes[1,0].loglog([r['configuration']['radial_cells'] for r in selected],[r['final_errors']['current_relative_volume_L2'] for r in selected],'o-',label=f'de/R={de:g}')
    de=.05;ref=InertialScreening(electron_length=de);radius=np.linspace(0,1,500)
    axes[1,1].plot(radius,ref.fields(radius,np.zeros_like(radius))['Jtheta'],'k--',label='continuum')
    for n in (8,16,32,64):
        name=f'screen_d0.05_{n}'
        if name not in arrays:continue
        a=arrays[name];ids=np.r_[0,1+np.arange(n)*(2*n+1)+n]
        r=a['final_coordinates'][ids,0]
        axes[1,1].plot(r,r*a['final_current_over_radius'][ids],'.-',label=f'nr={n}')
    for n in (8,16,32):
        values=[]
        for de in (.05,.2,1.,4.):
            a=arrays[f'screen_d{de:g}_{n}'];values.append(2*a['final_material_a'][0]/1e-8)
        axes[1,2].semilogx([.05,.2,1.,4.],values,'o-',label=f'nr={n}')
    lengths=np.logspace(np.log10(.05),np.log10(4.),200)
    axes[1,2].semilogx(lengths,[float(InertialScreening(electron_length=d).fields(0.,0.)['Bz']/1e-8) for d in lengths],'k--',label='continuum')
    titles=['Decay: magnetic field','Decay: current','Decay: local collisional heat','Collisionless screening: current error','Strong screening: actual equatorial current','Field penetration as carriers decrease']
    for ax,title in zip(axes.ravel(),titles):ax.set_title(title);ax.grid(alpha=.25,which='both');ax.legend(fontsize=8)
    for ax in (*axes[0],axes[1,0]):
        ax.set(xlabel='Radial cells',ylabel='Relative volume L2 error');ax.set_xticks([8,16,32,64],labels=[8,16,32,64]);ax.xaxis.set_minor_formatter(NullFormatter())
    axes[1,1].set(xlabel='Radius [m]',ylabel='Jθ [A/m²]',xlim=(.7,1.005))
    axes[1,2].set(xlabel='Electron inertial length / sphere radius',ylabel='Central Bz / applied bias')
    fig.suptitle('Fixed-material electron-current response: true magnetic vacuum, explicit inertia\nComponent verification; no moving heated-edge or formation claim',fontsize=15)
    fig.savefig(Path(directory)/'current_accuracy.png',dpi=140);plt.close(fig)


def regimes(directory,transport):
    fig,axes=plt.subplots(1,3,figsize=(16,5),layout='constrained')
    rows=[r for r in transport['cases'] if r['name'].startswith('Ar_') and r['electron_temperature_ev']==5. and r['structure_length_m']==.035]
    density=[r['electron_density_m3'] for r in rows]
    for key,label in [('hall_length_ratio','di / R'),('electron_inertia_over_length','de / R'),('electron_mean_free_path_over_length','electron mean free path / R')]:
        axes[0].loglog(density,[r[key] for r in rows],'o-',label=label)
    edge=transport['edge']['samples'];x=[r['distance_to_edge_m'] for r in edge]
    for key,label in [('hall_length_ratio','di / Lne'),('electron_inertia_over_length','de / Lne'),('debye_over_length','Debye / Lne')]:
        axes[1].loglog(x,[r[key] for r in edge],label=label)
    for key,label in [('drift_over_electron_thermal_speed','drift / electron thermal speed'),('drift_over_ion_thermal_speed','drift / ion thermal speed')]:
        axes[2].loglog(x,[r[key] for r in edge],label=label)
    axes[0].set(title='Conditional Ar+ bulk scales',xlabel='Electron density [m⁻³]',ylabel='Scale ratio')
    axes[1].set(title='Prescribed density edge: local scales',xlabel='Distance to edge [m]',ylabel='Scale / density-gradient length')
    axes[2].set(title='Fixed current scale toward the edge',xlabel='Distance to edge [m]',ylabel='Required carrier drift ratio')
    for ax in axes:
        ax.axhline(.1,color='gray',linestyle=':',label='0.1 review threshold');ax.grid(alpha=.25,which='both');ax.legend(fontsize=8)
    axes[1].invert_xaxis();axes[2].invert_xaxis()
    fig.suptitle('Explicit diagnostic assumptions: R=3.5 cm, Te=5 eV, Ti=1 eV, B=5 mT\nDensity edge and current are prescribed diagnostics, not a solved plasma boundary',fontsize=14)
    fig.savefig(Path(directory)/'transport_scales.png',dpi=140);plt.close(fig)
