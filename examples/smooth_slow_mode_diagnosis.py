"""Independent linear Ohm cancellation controls for the unresolved slow mode."""
from dataclasses import asdict
import numpy as np
from scipy.constants import elementary_charge as E,atomic_mass
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.current_thermal import CurrentThermal
from plasma_bench.kernel_induction import KernelInduction
from plasma_bench.primal_velocity import regular_velocity
from plasma_bench.dual_kinematics import velocity
from .finite_mass_wave_reference import FiniteMassWave
from .finite_mass_manufactured import cell_average
from .current_thermal_norms import vectors


def measure():
    cases=[];arrays={}
    ref=FiniteMassWave(.035,.24,1e18,39.948*atomic_mass,.005,6e18*E,gamma_i=1.4,branch=0,amplitude=.01)
    for nr in (8,16,32,64):
        mesh=MHDMesh(ref.radius,ref.length/2,nr,2*nr)
        model=CurrentThermal(mesh,np.full((nr+1,2*nr),ref.number_density),ref.ion_mass,gamma_i=ref.gamma_i)
        geo=model.geometry
        def linear_momenta(r,z):
            v=ref.fields(r,z)['velocity'];rho=np.full_like(r,ref.density)
            return np.stack((rho,rho*v[...,0],rho*v[...,2],np.zeros_like(r),rho*r*v[...,1]),axis=-1)
        matter=cell_average(mesh,linear_momenta)
        # Perturbation-only magnetic input excludes stationary-bias current noise
        # in this explicit linear diagnostic; full production inputs stay intact.
        psi=ref.fields(*mesh.corners)['psi'];psi[:,-1]=psi[:,0]
        g=cell_average(mesh,lambda r,z:ref.fields(r,z)['magnetic'][...,1]/r)
        packed=geo.pack(psi,g);j=model.current(packed)
        jnumerical=geo.sample_electric(j);vnumerical=velocity(geo,regular_velocity(model.mapping,matter))[0]
        exact=ref.fields(geo.r,geo.z);bias=np.array([0.,0.,ref.bias])
        fields={};summary={}
        for label,v,jj in (('both_numerical',vnumerical,jnumerical),('exact_velocity',exact['velocity'],jnumerical),
                           ('exact_current',vnumerical,exact['current']),('both_exact',exact['velocity'],exact['current'])):
            # Uniform-n pressure gradient has identically zero curl at linear
            # order and is omitted only in this explicit magnetic-rate control.
            electric=np.cross(-v+model.delta*jj/(ref.number_density*E),bias)
            forcing=model.mass@geo.project_electric(electric)
            solved=KernelInduction(model,tolerance=1e-10,max_iterations=2000).solve(forcing)
            fields[label]=solved.magnetic
            sampled=geo.magnetic_samples(*geo.unpack(solved.magnetic));target=ref.derivatives(geo.r,geo.z)[2]['magnetic'].real
            weight=geo.weight[...,None]
            summary[label]=float(np.sqrt(np.sum(weight*(sampled-target)**2)/np.sum(weight*target**2)))
        norm=lambda a:float(np.sqrt(np.sum(geo.weight[...,None]*a*a)))
        ideal=np.cross(-exact['velocity'],bias);hall=np.cross(model.delta*exact['current']/(ref.number_density*E),bias)
        cases.append(dict(nr=nr,reference=asdict(ref),magnetic_rate_relative_L2=summary,
            ideal_to_net_electric_norm=norm(ideal)/norm(ideal+hall),hall_to_net_electric_norm=norm(hall)/norm(ideal+hall)))
        for label,value in fields.items():arrays[f'nr{nr}_{label}']=value
        print('linear slow controls',nr,summary,flush=True)
    return dict(cases=cases,scope='Linear constant-density Ohm cancellation diagnosis. Exact substitutions are oracle controls, not solver fixes or new physical initial conditions.'),arrays
