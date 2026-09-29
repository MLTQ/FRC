"""Positive magnetic-space elimination of the fixed-n toroidal rate system."""
import numpy as np
from scipy.sparse.linalg import LinearOperator, cg
from plasma_bench.hall_matrices import transport_matrix


def discrete(model,psi,g,ue,ui):
    geo=model.geometry;mesh=model.mesh;nodes=(mesh.nr+1)*mesh.nz
    if np.any(psi):raise ValueError('This auxiliary rate solve requires a purely toroidal initial magnetic field')
    packed=geo.pack(psi,g);j=model.current(packed)
    star=model.star_field(packed,j);correction=star-geo.magnetic_samples(psi,g)
    pressure=model.pressure(ue,ui)[0]
    sources=(model.delta*transport_matrix(geo,star,model.density)@j+model.mass@pressure,
             model.delta*transport_matrix(geo,correction,model.density)@j)
    curl=geo.faraday[nodes:];root=np.sqrt(model.hessian.diagonal()[nodes:])
    def mass_inverse(v):
        return np.concatenate([solve(a) for solve,a in zip(geo.solves,geo.split_electric(v))])
    def operator(y):
        current=mass_inverse(curl.T@(root*y))
        return y+root*(curl@mass_inverse(model.inertia@current))
    system=LinearOperator((len(root),len(root)),matvec=operator,dtype=float)
    values=[];residuals=[]
    for source in sources:
        # Check that the omitted poloidal magnetic-rate subspace has no forcing.
        electric_source=mass_inverse(source)
        if np.any((geo.faraday@electric_source)[:nodes]):
            raise ValueError('Non-toroidal forcing is outside the auxiliary Schur solve')
        rhs=root*(curl@electric_source)
        solution,info=cg(system,rhs,rtol=2e-13,atol=0.,maxiter=1000)
        if info!=0:raise RuntimeError(f'Positive magnetic-rate solve did not converge: {info}')
        rate=solution/root
        current_rate=-mass_inverse(curl.T@(root*solution))
        electric=mass_inverse(source+model.inertia@current_rate)
        faraday=curl@electric
        # Reconstruct the eliminated variables and check the original equations.
        ohm=model.mass@electric-model.inertia@current_rate-source
        ampere=model.mass@current_rate+curl.T@(root*solution)
        norms=(np.linalg.norm(operator(solution)-rhs)/max(np.linalg.norm(rhs),np.finfo(float).tiny),
               np.linalg.norm(root*(rate-faraday))/max(np.linalg.norm(root*rate),np.finfo(float).tiny),
               np.linalg.norm(ohm)/max(np.linalg.norm(source)+np.linalg.norm(model.inertia@current_rate),np.finfo(float).tiny),
               np.linalg.norm(ampere)/max(np.linalg.norm(curl.T@(root*solution)),np.finfo(float).tiny))
        if max(norms)>1e-10:raise RuntimeError(f'Reconstructed mixed-equation residual failed: {norms}')
        values.append(rate.reshape(mesh.nr,mesh.nz));residuals.append(float(max(norms)))
    return values,residuals
