"""Independent exponential of the semidiscrete helicity equations."""
import numpy as np
from scipy.linalg import expm
from scipy.constants import mu_0 as MU0


def exact_spatial_evolution(solver, initial, time):
    """Separate temporal error from the FE spatial error; no midpoint update used."""
    p = solver.parameters
    final = np.zeros_like(initial)
    for radial, eigenvalue in enumerate(solver.radial_eigenvalues):
        for axial, k in enumerate(solver.kz):
            total2 = eigenvalue+k*k
            inertia = 1+p.current_length_squared*total2
            for h in (-1, 1):
                hall = h*p.hall_coefficient*p.bias*k*np.sqrt(total2)/MU0
                matrix = np.array([[(-p.resistivity/MU0*total2-1j*hall)/inertia, 1j*k*p.alfven/inertia],
                                   [1j*k*p.alfven, 0.]])
                pair = np.array([initial[radial, axial, 0]+h*initial[radial, axial, 1],
                                 initial[radial, axial, 2]+h*initial[radial, axial, 3]])/2
                magnetic, velocity = expm(time*matrix)@pair
                final[radial, axial] += np.array([magnetic, h*magnetic, velocity, h*velocity])
    return final
