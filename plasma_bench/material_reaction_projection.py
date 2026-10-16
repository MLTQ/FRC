"""Conservative reaction moments to variable-mass material support blocks."""
import numpy as np
from scipy.constants import elementary_charge as EV
from scipy.linalg import solve_banded
from .species_material import MaterialSpecies,SpeciesMaterialState
from .material_species_set import MaterialSpeciesSet,MaterialSetState
from .material_subcells import MaterialSubcells
from .material_reaction_partition import select_partition,field_projection_error


def project_reactions(overlap,state,*,cfl=.15,maximum_plasma_phase=.04,numerical_viscosity=(.5,2.),viscosity_limiter='minmod',coalesce_ulps=0,velocity_cells=None,minimum_density_width=0.):
    """Recover material velocity and thermal energy from full evolved moments.

    Masses and charges are cell-averaged on a shared interval partition. Nodal
    momentum is an L2 projection; thermal is actual total energy minus resolved
    kinetic energy, with no clipping or redistribution if inadmissible.
    """
    source=overlap.source;network=source.network;values=state.moments;density=overlap.density
    numbers=values[...,overlap.number]*density;mass=source.mass
    # Keep the source's fixed coordinates through the material return mapping.
    momentum=values[...,source.momentum]*source.momentum_scale*density
    energy=values[...,source.energy]*EV*density
    amounts=np.sum(state.weights[...,None]*numbers,axis=1)
    original_amounts=amounts.copy()
    total_energy=np.sum(state.weights[...,None]*energy,axis=1)
    # Never remove an endpoint of any occupied species support, even when that
    # support is extremely thin. Only redundant internal partition cuts qualify.
    present=np.any(values[...,overlap.number]>0,axis=1) if velocity_cells is not None else amounts>0
    kept,partition_record=select_partition(state.edges,present,roundoff_ulps=coalesce_ulps,minimum_width=minimum_density_width)
    edges=state.edges[kept];groups=np.searchsorted(kept,np.arange(len(state.edges)-1),side='right')-1
    count=len(edges)-1
    right=(state.points-edges[groups,None])/np.diff(edges)[groups,None]
    weighted=state.weights[...,None,None]*momentum
    load_left=np.zeros((count,len(mass),3));load_right=np.zeros_like(load_left)
    np.add.at(load_left,groups,np.sum(weighted*(1-right[...,None,None]),axis=1))
    np.add.at(load_right,groups,np.sum(weighted*right[...,None,None],axis=1))
    merged_n=np.zeros((count,len(mass)));merged_e=np.zeros_like(merged_n)
    np.add.at(merged_n,groups,amounts);np.add.at(merged_e,groups,total_energy)
    amounts,total_energy=merged_n,merged_e
    models=[];states=[];ids=[];records=[]
    for species in range(len(mass)):
        if velocity_cells is not None:
            yscale=float(values[...,overlap.number[species]].max())
            if yscale==0:continue
            population_scale=density*yscale
            nn=np.zeros(count);ll=np.zeros((count,3));lr=np.zeros_like(ll);ee=np.zeros(count)
            np.add.at(nn,groups,np.sum(state.weights*(values[...,overlap.number[species]]/yscale),axis=1))
            pm=(values[...,source.momentum[species]]/yscale)*source.momentum_scale
            em=(values[...,source.energy[species]]/yscale)*EV
            np.add.at(ll,groups,np.sum(state.weights[...,None]*pm*(1-right[...,None]),axis=1))
            np.add.at(lr,groups,np.sum(state.weights[...,None]*pm*right[...,None],axis=1))
            np.add.at(ee,groups,np.sum(state.weights*em,axis=1))
            occupied=nn>0
            expected=np.zeros(count,bool);np.logical_or.at(expected,groups,present[:,species])
            if np.any(expected&~occupied):raise ValueError('Positive normalized population is unrepresentable; no discarded species')
        else:occupied=amounts[:,species]>0
        changes=np.flatnonzero(np.diff(np.r_[False,occupied,False]));pairs=changes.reshape(-1,2)
        for first,last in pairs:
            faces=edges[first:last+1];width=np.diff(faces);n=amounts[first:last,species]
            if velocity_cells is not None:
                model=MaterialSubcells(nn[first:last],mass[species],overlap.charge[species],(state.edges[0],state.edges[-1]),faces,velocity_cells,number_scale=population_scale,
                    left_electric=state.left_electric,cfl=cfl,maximum_plasma_phase=maximum_plasma_phase,
                    numerical_viscosity=numerical_viscosity,viscosity_limiter=viscosity_limiter)
                load=np.zeros((len(faces),3));load[:-1]+=ll[first:last];load[1:]+=lr[first:last]
                load=model.restrict_load(load);faces=np.asarray(faces,np.longdouble)
            else:
                model=MaterialSpecies(n[None,:],[mass[species]],[overlap.charge[species]],(state.edges[0],state.edges[-1]),
                    left_electric=state.left_electric,cfl=cfl,maximum_plasma_phase=maximum_plasma_phase,
                    numerical_viscosity=numerical_viscosity,viscosity_limiter=viscosity_limiter)
                load=np.zeros((len(faces),3));load[:-1]+=load_left[first:last,species];load[1:]+=load_right[first:last,species]
            velocity=solve_banded((1,1),model.matrix[0],load,check_finite=False)
            vl,vr=model.fine_velocity(velocity[None,...]) if velocity_cells is not None else (velocity[:-1],velocity[1:])
            scale=model.number_scale if velocity_cells is not None else 1.
            relative_n=model.relative_numbers if velocity_cells is not None else n
            kinetic=mass[species]*relative_n/6*np.sum(vl*vl+vl*vr+vr*vr,axis=1)
            thermal=(ee[first:last] if velocity_cells is not None else total_energy[first:last,species])-kinetic
            result=SpeciesMaterialState(faces[None,:].copy(),velocity[None,...],thermal[None,:]);model.validate(result)
            input_p=(np.sum(ll[first:last]+lr[first:last],axis=0)*scale if velocity_cells is not None else np.sum(load_left[first:last,species]+load_right[first:last,species],axis=0))
            output_p=np.sum(mass[species]*relative_n[:,None]*(vl+vr)/2,axis=0)*scale
            records.append(dict(species=int(species),cells=len(n),support_m=[float(faces[0]),float(faces[-1])],
                velocity_cells=getattr(model,'velocity_cells',len(n)),
                minimum_thermal_j_m2=float(thermal.min()*scale),thermal_coordinate_scale=scale,momentum_difference_kg_per_m_s=(output_p-input_p).tolist(),
                minimum_cell_width_m=float(width.min())))
            models.append(model);states.append(result);ids.append(species)
    if not models:raise ValueError('No occupied material support; retain exact vacuum without a matter solver')
    model=MaterialSpeciesSet(models,ids,reference_velocity=overlap.frame);result=MaterialSetState(tuple(states));model.validate(result)
    reservoirs=dict(photon_energy_j_m2=float(EV*density*np.sum(state.weights*values[...,network.photons])),
                    input_heat_j_m2=float(EV*density*np.sum(state.weights*values[...,network.heating])),
                    binding_energy_j_m2=float(EV*np.sum(amounts[:,:network.count]*network.energy_ev)))
    field_record=field_projection_error(state.edges,original_amounts,overlap.charge,state.left_electric,model,result)
    return model,result,dict(blocks=records,**partition_record,**field_record,**reservoirs)
