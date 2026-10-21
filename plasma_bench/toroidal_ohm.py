"""Mimetic toroidal currents, resistive E_r/E_z and a skew-adjoint Hall field for the 2.5D conservative path."""
from dataclasses import dataclass
from functools import lru_cache
import numpy as np
from scipy.constants import mu_0 as MU0, e as QE


@dataclass(frozen=True)
class ToroidalGeometry:
    """Dual volumes [m^3] that pair each discrete current with its electric field.

    cell (nr, nz); radial faces (nr+1, nz) at r = rf; axial faces (nr, nz+1) at z = zf; nodes (nr+1, nz+1).
    The axis face holds the disk r < r_1/sqrt(2) (pi <r^2>_0 dz), so J_z(0) = 2 g_0/mu0; wall, inner-wall
    and end faces hold the half cell inside the gas. Node volumes 2 pi r A_n match the solver's
    poloidal Ohmic quadrature; axis nodes (no E_theta) have zero volume.
    """
    r2: np.ndarray
    cell: np.ndarray
    radial: np.ndarray
    axial: np.ndarray
    node: np.ndarray


@lru_cache(maxsize=8)
def geometry(mesh):
    dr, dz = mesh.dr, mesh.dz
    rf, rc = mesh.rf, mesh.r
    r2 = .5*(rf[:-1]**2+rf[1:]**2)[:, None]
    radial = np.broadcast_to((2*np.pi*rf*dr*dz)[:, None], (mesh.nr+1, mesh.nz)).copy()
    radial[0] = np.pi*rf[0]*dr*dz if mesh.annular else np.pi*r2[0, 0]*dz
    radial[-1] *= .5
    axial = np.broadcast_to((2*np.pi*rc*dr*dz)[:, None], (mesh.nr, mesh.nz+1)).copy()
    axial[:, [0, -1]] *= .5
    area = np.full((mesh.nr+1, mesh.nz+1), dr*dz)
    area[[0, -1]] *= .5
    area[:, [0, -1]] *= .5
    node = 2*np.pi*rf[:, None]*area
    for array in (r2, radial, axial, node):
        array.setflags(write=False)
    return ToroidalGeometry(r2, np.asarray(mesh.volumes), radial, axial, node)


def currents(mesh, g, boundary=0.):
    """a = <r^2> g (cell r B_theta) and the currents conjugate to the B_theta/r update.

    With the cell energy <r^2> g^2/(2 mu0) and Faraday's law on each r-z cell, summation by parts gives
    dW_theta/dt = -sum(E_z J_z V_radial) - sum(E_r J_r V_axial) exactly, for
        J_z = 2 pi dz (a_i - a_{i-1}) / (mu0 V_radial),   J_r = -2 pi dr (a_j - a_{j-1}) / (mu0 V_axial),
    with a = 0 beyond the axis (r B_theta = 0) and a = `boundary` beyond the walls and ends: the r B_theta
    of the exact vacuum outside an insulating vessel, mu0 I_link/(2 pi) for a current I_link threading it
    (zero without one; only an annular chamber can be threaded)."""
    geo = geometry(mesh)
    a = geo.r2*g
    radial = np.zeros((mesh.nr+2, mesh.nz))
    radial[1:-1] = a
    axial = np.zeros((mesh.nr, mesh.nz+2))
    axial[:, 1:-1] = a
    if boundary:
        radial[-1] = boundary
        if mesh.annular:
            radial[0] = boundary
        axial[:, [0, -1]] = boundary
    jz = 2*np.pi*mesh.dz*np.diff(radial, axis=0)/(MU0*geo.radial)
    jr = -2*np.pi*mesh.dr*np.diff(axial, axis=1)/(MU0*geo.axial)
    return a, jz, jr


def face_toroidal_field(mesh, a, boundary=0.):
    """B_theta on radial and axial faces for the Poynting flux: (a_{k-1}+a_k)/2 divided by rf or r_c.

    Zero on the axis; `boundary`/r (the vacuum r B_theta, normally zero) on every gas boundary. With
    S_r = -E_z B_theta/mu0 and S_z = E_r B_theta/mu0, each cell's gas then receives exactly half the J.E of
    each interior face it touches and all of it on its boundary faces, so a current-free region exchanges
    no energy with the gas; the boundary faces carry the vacuum's Poynting flux."""
    radial = np.zeros((mesh.nr+1, mesh.nz))
    radial[1:-1] = .5*(a[:-1]+a[1:])/mesh.rf[1:-1, None]
    axial = np.zeros((mesh.nr, mesh.nz+1))
    axial[:, 1:-1] = .5*(a[:, :-1]+a[:, 1:])/mesh.r[:, None]
    if boundary:
        radial[-1] = boundary/mesh.rf[-1]
        if mesh.annular:
            radial[0] = boundary/mesh.rf[0]
        axial[:, [0, -1]] = boundary/mesh.r[:, None]
    return radial, axial


def node_currents(jz, jr):
    """|J_z|, |J_r| averaged to nodes (for the Chodura drift only), shape (nr+1, nz+1)."""
    nr1, nz = jz.shape
    z = np.zeros((nr1, nz+1))
    z[:, :-1] += .5*jz; z[:, 1:] += .5*jz
    z[:, [0, -1]] *= 2
    r = np.zeros((nr1, nz+1))
    r[:-1] += .5*jr; r[1:] += .5*jr
    r[[0, -1]] *= 2
    return z, r


def face_average(cell, axis):
    """Cell values to the faces normal to `axis` (0: radial faces, 1: axial faces); boundaries one-sided."""
    pad = np.pad(cell, ((1, 1), (0, 0)) if axis == 0 else ((0, 0), (1, 1)), mode="edge")
    return .5*(pad[:-1]+pad[1:]) if axis == 0 else .5*(pad[:, :-1]+pad[:, 1:])


def ohmic(mesh, eta_nodes, jz, jr):
    """Resistive E_z, E_r with node eta averaged to the faces, and the exact paired heat rate [W]."""
    geo = geometry(mesh)
    ez = .5*(eta_nodes[:, :-1]+eta_nodes[:, 1:])*jz
    er = .5*(eta_nodes[:-1]+eta_nodes[1:])*jr
    return ez, er, float(np.sum(ez*jz*geo.radial)+np.sum(er*jr*geo.axial))


def face_heat_to_cells(mesh, radial_power, axial_power):
    """Per-cell heating density [W/m^3] from face J.E [W]: half to each neighbour, all of a boundary face
    to its single cell."""
    geo = geometry(mesh)
    heat = np.zeros((mesh.nr, mesh.nz))
    share = radial_power.copy()
    share[1:-1] *= .5
    heat += share[:-1]+share[1:]
    share = axial_power.copy()
    share[:, 1:-1] *= .5
    heat += share[:, :-1]+share[:, 1:]
    return heat/geo.cell


def hall_fields(mesh, g, jz, jr, jtheta, br_faces, bz_faces, coefficient):
    """Hall electric field E_H = c J x B, c = w(n)/(n e) per cell, built skew-adjoint.

    E_theta (nodes) = c (J_z B_r - J_r B_z); E_z (radial faces) = c (J_r B_theta - J_theta B_r);
    E_r (axial faces) = c (J_theta B_z - J_z B_theta). Each product is formed where its normal field lives
    (B_r on radial faces, B_z on axial faces, B_theta in cells) and moved to the other location with
    the adjoint of the averaging that brings the current there. Then
        sum(J_theta E_theta V_node) + sum(J_z E_z V_radial) + sum(J_r E_r V_axial) = 0
    to roundoff for any state: the discrete Hall field does no work.
    Returns (E_theta nodes (nr+1, nz+1), E_z (nr+1, nz), E_r (nr, nz+1)).
    """
    geo = geometry(mesh)
    c_radial, c_axial = face_average(coefficient, 0), face_average(coefficient, 1)
    # (1) J_theta J_z B_r on radial faces
    x = geo.radial*c_radial*br_faces*jz
    # (2) J_theta J_r B_z on axial faces
    y = geo.axial*c_axial*bz_faces*jr
    etheta = np.zeros((mesh.nr+1, mesh.nz+1))
    etheta[:, :-1] += .5*x; etheta[:, 1:] += .5*x
    etheta[:-1] -= .5*y; etheta[1:] -= .5*y
    etheta = np.divide(etheta, geo.node, out=np.zeros_like(etheta), where=geo.node > 0)
    jtheta_radial = .5*(jtheta[:, :-1]+jtheta[:, 1:])
    jtheta_axial = .5*(jtheta[:-1]+jtheta[1:])
    ez = -c_radial*br_faces*jtheta_radial
    er = c_axial*bz_faces*jtheta_axial
    # (3) J_z J_r B_theta in cells
    z = geo.cell*coefficient*np.sqrt(geo.r2)*g
    jr_cell, jz_cell = .5*(jr[:, :-1]+jr[:, 1:]), .5*(jz[:-1]+jz[1:])
    zr, zz = .5*z*jr_cell, .5*z*jz_cell
    ez[:-1] += zr/geo.radial[:-1]; ez[1:] += zr/geo.radial[1:]
    er[:, :-1] -= zz/geo.axial[:, :-1]; er[:, 1:] -= zz/geo.axial[:, 1:]
    return etheta, ez, er


def hall_coefficient(rho, ion_mass, cutoff_density=None):
    """c = w/(n e) [m^3/C] with n = rho/ion_mass and the DECLARED cutoff w = 1/(1+(rho_c/rho)^4).

    The cutoff removes the Hall term from near-vacuum cells (e.g. the 1%-of-fill gap floor), whose
    whistler step dx^2/d_i would otherwise set the time step. Its maximum effect, c <= 0.57 ion_mass/
    (rho_c e) at rho = 1.32 rho_c, bounds the whistler rate. Returns (c, w)."""
    weight = np.ones_like(rho) if cutoff_density is None else 1/(1+(cutoff_density/rho)**4)
    return weight*ion_mass/(rho*QE), weight


def work(mesh, jtheta, etheta, jz, ez, jr, er):
    """sum J.E over the paired dual volumes [W]."""
    geo = geometry(mesh)
    return float(np.sum(jtheta*etheta*geo.node)+np.sum(jz*ez*geo.radial)+np.sum(jr*er*geo.axial))
