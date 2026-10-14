"""Show saved free-boundary trajectories, convergence and physical rotation."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.constants import electron_mass as ME


def load(path):
    path=Path(path);data=json.loads((path/'summary.json').read_text())
    if not data['terminal'] or not data['success']:raise ValueError('Complete successful trajectories required')
    return path,data


def run(destination):
    destination=Path(destination);destination.mkdir(exist_ok=False,parents=True)
    sphere=[load(f'outputs/material_vacuum_sphere_r{n}_t800'+('_v2' if n==4 else '')) for n in (4,8,16)]
    rotating=[load(f'outputs/material_vacuum_rotating_r{n}_t800') for n in (8,16)]
    fig,axes=plt.subplots(2,2,figsize=(11.8,9),layout='constrained')
    for ax,(path,d),title in zip(axes[0],(sphere[-1],rotating[-1]),('Spherical expansion','Rotating elongated expansion')):
        initial=np.asarray(d['initial']['boundary_coordinates']);last=d['records'][-1];actual=np.asarray(last['boundary_coordinates'])
        angle=np.linspace(0,2*np.pi,300);a=last['radial_axis_reference_m'];c=last['axial_axis_reference_m']
        for points,color,label in ((initial,'#87929a','Initial free surface'),(actual,'#007f86','Evolved free surface')):
            transverse=np.r_[points[:,0],-points[::-1,0]];z=np.r_[points[:,1],points[::-1,1]]
            ax.plot(transverse*100,z*100,color=color,label=label,lw=2)
        ax.plot(100*a*np.sin(angle),100*c*np.cos(angle),'--',color='#cb6b28',label='Independent affine solution',lw=1.1)
        ax.set(xlabel='Transverse position (cm)',ylabel='Axial position (cm)',title=title,aspect='equal');ax.legend(fontsize=8)
        ax.text(.02,.02,'Exact empty exterior',transform=ax.transAxes,color='#555',fontsize=9)
    colors=['#9ba6ac','#607d8b','#007f86']
    for (path,d),color in zip(sphere,colors):
        rows=d['records'];t=np.array([v['time'] for v in rows])*1e6
        error=np.array([v['maximum_boundary_position_error_over_initial_radius'] for v in rows])*d['initial_radial_axis_m']*1e6
        axes[1,0].semilogy(t,error,color=color,label=f"{d['radial_cells']}×{d['angular_cells']} material cells")
    axes[1,0].set(xlabel='Time (µs)',ylabel='Maximum surface-position error (µm)',title='Observed spatial convergence');axes[1,0].legend(fontsize=8)
    for (path,d),color in zip(rotating,colors[1:]):
        initial=np.load(path/'initial.npz');particles=initial['particles'];ids=initial['ids'];q,_=np.polynomial.legendre.leggauss(d['quadrature'])
        xi,tau=np.meshgrid((q+1)/2,(q+1)/2,indexing='ij');xi=xi.ravel();tau=tau.ravel()
        shape=np.stack(((1-xi)*(1-tau),xi*(1-tau),xi*tau,(1-xi)*tau),axis=-1)
        values=[];times=[];exact=[]
        for row in [d['initial'],*d['records']]:
            state=initial if row['time']==0 else np.load(path/row['arrays_file'])
            radius=np.einsum('qn,cn->cq',shape,state['coordinates'][ids,0]);mass=d['ion_mass_kg']+ME
            inertia=float(np.sum(mass*particles*radius**2));values.append(float(state['angular_momentum'].sum()/inertia)/d['rotation_s'])
            times.append(row['time']*1e6);exact.append((d['initial_radial_axis_m']/row['radial_axis_reference_m'])**2)
        axes[1,1].plot(times,values,color=color,label=f"{d['radial_cells']}×{d['angular_cells']} cells")
    axes[1,1].plot(times,exact,'--',color='#cb6b28',label='Independent angular-momentum law')
    axes[1,1].set(xlabel='Time (µs)',ylabel='Inertia-weighted angular speed / initial',title='Rotation slows during expansion');axes[1,1].legend(fontsize=8)
    for ax in axes.ravel():ax.grid(alpha=.18)
    fig.suptitle('Moving two-species fluid • physical-vacuum expansion',fontsize=15)
    fig.supxlabel('Zero-current Euler limit; density and pressure vanish at the surface. Magnetic coupling and edge transport remain unfinished.',fontsize=10)
    fig.savefig(destination/'material_vacuum_expansion.png',dpi=180);fig.savefig(destination/'material_vacuum_expansion.pdf');plt.close(fig)
    (destination/'inputs.json').write_text(json.dumps([dict(path=str(p),source_sha256=d['source_sha256']) for p,d in [*sphere,*rotating]],indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('destination');run(p.parse_args().destination)
