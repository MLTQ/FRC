"""Actual material-fluid wall shock versus independent Rankine-Hugoniot jumps."""
import json
from pathlib import Path
from dataclasses import replace
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import atomic_mass, elementary_charge as E
from scipy.sparse.linalg import splu
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.material_cylinder import MaterialCylinder
from plasma_bench.material_two_species import MaterialTwoSpecies
from plasma_bench.material_current_space import MaterialCurrentSpace
from plasma_bench.material_dual_entropy import MaterialDualEntropy
from plasma_bench.material_bulk_viscosity import bulk_viscosity
from .implicit_plasma_records import provenance,freeze


def run(destination, cells=(50,100,200), time_refinement=1):
    out=Path(destination)
    if out.exists(): raise ValueError('Fresh shock output directory required')
    out.mkdir();identity=provenance(['examples.material_wall_shock']);freeze(out,identity)
    gamma=5/3;incoming=1e4;length=1.;duration=.2*length/incoming
    mass=39.948*atomic_mass;density_number=1e18;results=[]
    fig,axes=plt.subplots(1,3,figsize=(14,4),constrained_layout=True)
    if time_refinement < 1 or min(cells) < 4:
        raise ValueError('Positive time refinement and at least four axial cells required')
    for nz in cells:
        mesh=MHDMesh(.05,.5*length,4,nz);domain=MaterialCylinder(mesh,quadrature=4)
        shape=(4*nz,16);n=np.full(shape,density_number)
        # p0/(rho0*u0²)=0.01, split equally between the two initial species.
        from scipy.constants import electron_mass
        rho0=(mass+electron_mass)*density_number;p0=.01*rho0*incoming**2
        temperature=np.full(shape+(2,),p0/(2*E*density_number))
        velocity=np.zeros_like(domain.fluid_reference);velocity[:,1]=incoming
        right=domain.ids[:,-1];left=domain.ids[:,0];velocity[right,1]=0.
        fluid=MaterialTwoSpecies(domain,domain.fluid_reference,n,temperature,mass,velocity=velocity)
        species=MaterialDualEntropy(fluid,MaterialCurrentSpace(mesh,degree=1,quadrature=4))
        state=fluid.initial;matrix=fluid.gram(fluid.mass*fluid.particles)
        fixed=np.r_[left,right];free=np.setdiff1d(np.arange(fluid.nodes),fixed);factor=splu(matrix[free][:,free])
        def derivative(s):
            fields=fluid.fields(s);rhs=fluid.derivative(s);visc=bulk_viscosity(fluid,s,fields)
            force=rhs.momentum+visc['force'];acc=np.zeros(fluid.nodes);acc[free]=factor.solve(force[free,1])
            axial=matrix@acc;reaction=axial-force[:,1]
            heat=np.zeros(fields['temperature'].shape);heat[...,1]=visc['heating']
            entropy=species.heating_rate(dict(fields=fields,volume=fields['volume'].ravel()),heat.reshape(-1,2))
            rhs=replace(rhs,coordinates=np.column_stack((np.zeros(fluid.nodes),fields['velocity'][:,1])),
                momentum=np.column_stack((np.zeros(fluid.nodes),axial)),entropy=entropy)
            return rhs,float(reaction@fields['velocity'][:,1]),visc['power']
        initial_energy=fluid.fields(state)['energy'];steps=8*nz*time_refinement;dt=duration/steps;work=0.;heat=0.
        for step in range(steps):
            k1,w1,h1=derivative(state);k2,w2,h2=derivative(state.added(k1,.5*dt))
            k3,w3,h3=derivative(state.added(k2,.5*dt));k4,w4,h4=derivative(state.added(k3,dt))
            state=state.added(k1,dt/6).added(k2,dt/3).added(k3,dt/3).added(k4,dt/6)
            work+=dt*(w1+2*w2+2*w3+w4)/6;heat+=dt*(h1+2*h2+2*h3+h4)/6
        fields=fluid.fields(state);z=fields['points'][...,1];volume=fields['volume']
        # RH with stationary downstream, incoming u0, shock speed -s.
        a=.5*(gamma-1)*incoming;cs2=gamma*p0/rho0
        speed=.5*(-(incoming-a)+np.sqrt((incoming+a)**2+4*cs2))
        compression=(incoming+speed)/speed;pressure2=p0+rho0*incoming*(incoming+speed)
        front=.5*length-speed*duration;behind=z>=front
        exact=(np.where(behind,compression,1.),np.where(behind,0.,1.),np.where(behind,pressure2/p0,1.))
        actual=(fields['n']/density_number,fields['sampled_velocity'][...,1]/incoming,np.sum(fields['pressure'],axis=-1)/p0)
        errors=[float(np.sum(volume*abs(x-y))/np.sum(volume*abs(y))) for x,y in zip(actual,exact)]
        result=dict(cells=nz,steps=steps,duration_s=duration,dt_s=dt,shock_position_m=front,
            exact_compression=compression,exact_pressure_ratio=pressure2/p0,
            relative_volume_l1_errors=dict(zip(('density','velocity','pressure'),errors)),
            relative_energy_work_residual=(fields['energy']-initial_energy-work)/initial_energy,
            numerical_viscous_heat_j=heat,boundary_work_j=work)
        results.append(result);print(result,flush=True)
        np.savez_compressed(out/f'wall_shock_{nz}.npz',**{k:getattr(state,k) for k in state.__dataclass_fields__},particles=fluid.particles,
            points=fields['points'],density=fields['n'],pressure=fields['pressure'],velocity=fields['sampled_velocity'])
        for ax,value in zip(axes,actual):
            zz=np.sum(volume*z,axis=1)/volume.sum(axis=1);vv=np.sum(volume*value,axis=1)/volume.sum(axis=1)
            ax.plot(zz[:nz],vv[:nz],label=f'{nz} axial cells')
    for ax,values,title in zip(axes,((1,compression),(1,0),(1,pressure2/p0)),('Density / upstream','Velocity / incoming','Pressure / upstream')):
        ax.plot([-.3,front,front,.5],[values[0],values[0],values[1],values[1]],'k--',label='Exact shock')
        ax.set(xlim=(.32,.5),xlabel='Axial position [m]',title=title);ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Moving-material wall shock · numerical bulk viscosity · independent jump conditions')
    fig.savefig(out/'wall_shock.png',dpi=160);fig.savefig(out/'wall_shock.pdf');plt.close(fig)
    (out/'summary.json').write_text(json.dumps(dict(source_sha256=identity['source_sha256'],results=results,
       scope='One-dimensional rigid-side wall-shock fixture; no magnetic field/current or plasma-vacuum formation claim.'),indent=2)+'\n')


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('destination')
    p.add_argument('--cells',type=int,nargs='+',default=[50,100,200]);p.add_argument('--time-refinement',type=int,default=1)
    a=p.parse_args();run(a.destination,a.cells,a.time_refinement)
