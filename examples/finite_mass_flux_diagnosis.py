"""Separate reconstruction and dissipation errors in the slow-mode mass flux."""
from dataclasses import replace
import numpy as np
from scipy.constants import atomic_mass, elementary_charge as E
from plasma_bench.mhd_mesh import MHDMesh
from plasma_bench.finite_mass_state import evaluate
from plasma_bench.finite_mass_flux import scalar_faces, physical, rusanov
from plasma_bench.toroidal_reconstruction import ToroidalReconstruction
from plasma_bench.toroidal_state import radius_squared
from plasma_bench.annular_polynomial import normal_momentum_flux
from .finite_mass_wave_reference import FiniteMassWave
from .finite_mass_wave_case import initialize
from .finite_mass_manufactured import cell_average
from .moving_species_reference import MovingSpeciesReference, parent_sources


def measure():
    records=[];arrays={}
    for density in (1e18,1e20):
        ref=FiniteMassWave(.035,.24,density,39.948*atomic_mass,.005,6*density*E,
                          gamma_i=1.4,amplitude=1.,branch=0)
        for nr in (8,16,32,64,128):
            mesh=MHDMesh(ref.radius,ref.length/2,nr,2*nr)
            state,par=initialize(mesh,ref);ev=evaluate(mesh,state,par)
            faces=ToroidalReconstruction(mesh,ev.primitive,order=par.spatial_order).faces()
            fractions=scalar_faces(mesh,(par.gamma_e-1)*ev.components['electron']/ev.primitive[...,3],par.spatial_order,upper=1.)
            kinetic=scalar_faces(mesh,ev.components['current'],par.spatial_order)
            fields=mesh.face_fields(state.psi);central=[];total=[]
            for axis,radius in ((0,mesh.rf[:,None]),(1,np.sqrt(radius_squared(mesh)))):
                a,b=(faces[2*axis+i].copy() for i in (0,1))
                a[...,4+axis]=b[...,4+axis]=fields[axis]
                ul,fl,al=physical(a,fractions[2*axis],kinetic[2*axis],radius,axis,par)
                ur,fr,ar=physical(b,fractions[2*axis+1],kinetic[2*axis+1],radius,axis,par)
                area=2*np.pi*(mesh.rf[:,None]*mesh.dz if axis==0 else mesh.radial_measure[:,None])
                central.append(.5*(fl[...,0]+fr[...,0])*area)
                total.append(rusanov(a,b,*fractions[2*axis:2*axis+2],*kinetic[2*axis:2*axis+2],fields[axis],radius,axis,par)[...,0]*area)
            divergence=lambda flux:-(np.diff(flux[0],axis=0)+np.diff(flux[1],axis=1))/mesh.volumes
            exact=cell_average(mesh,lambda r,z:ref.derivatives(r,z)[2]['density'].real)
            def linear(r,z):
                v=ref.fields(r,z)['velocity'];zero=np.zeros_like(r)
                return np.stack((zero,ref.density*v[...,0],ref.density*v[...,2],zero,ref.density*r*v[...,1]),axis=-1)
            linear_state=cell_average(mesh,linear)
            values=dict(total=divergence(total),central=divergence(central),
                diffusive=divergence([a-b for a,b in zip(total,central)]),
                cubic_finite_amplitude=divergence(normal_momentum_flux(mesh,state.conserved)),
                cubic_linear=divergence(normal_momentum_flux(mesh,linear_state)))
            norm=lambda value:float(np.sqrt(np.sum(mesh.volumes*value*value)))
            errors={key:norm(value-(0. if key=='diffusive' else exact))/norm(exact) for key,value in values.items()}
            record=dict(number_density_m3=density,nr=nr,omega_s=ref.omega,mass_rate_relative_L2=errors)
            records.append(record);prefix=f'n{int(np.log10(density))}_{nr}_'
            arrays.update({prefix+key:value for key,value in dict(exact=exact,**values).items()})
            print('flux diagnosis',density,nr,errors,flush=True)
    # The source's regular-axis limit is independent of the production solver.
    controls=[];base=MovingSpeciesReference()
    for density in (1e18,1e20):
        ref=replace(base,density=density);z=np.linspace(-ref.length/2,ref.length/2,129,endpoint=False)
        for phase in (.23,.48,.73):
            heats=[parent_sources(ref,np.full_like(z,fraction*ref.radius),z,phase/ref.frequency)[1]['e']
                   for fraction in (1e-6,1e-8,1e-10,1e-12)]
            scale=max(float(abs(heats[-1]).max()),1e-300)
            differences=[float(abs(value-heats[-1]).max()/scale) for value in heats[:-1]]
            assert differences[-1]<1e-12
            controls.append(dict(number_density_m3=density,phase=phase,
                radius_fractions=[1e-6,1e-8,1e-10],reference_radius_fraction=1e-12,
                relative_maximum_heat_difference=differences))
    return dict(cases=records,axis_heat_limit_controls=controls,
        scope='Diagnostic continuity reconstruction only; cubic flux is not adopted by the coupled solver. Diffusive entry is its magnitude, other entries are errors relative to the exact linear mass rate.'),arrays
