"""Independent analytic fixtures for the developmental MHD kernel."""
import numpy as np
from scipy.special import j0, j1, jn_zeros
from scipy.optimize import brentq
from scipy.constants import mu_0 as MU0
from .mhd_state import from_primitive


def radial_average(mesh, function):
    """Gauss quadrature of a smooth function over each radial annulus."""
    x, weights = np.polynomial.legendre.leggauss(8)
    r = mesh.r[:, None]+mesh.dr*x[None, :]/2
    return np.sum(function(r)*r*weights[None, :], axis=1)*mesh.dr/(2*mesh.radial_measure)


def radial_fast_mode(mesh, *, epsilon=1e-5, gamma=5/3, time=0.0):
    """Linear cylindrical fast mode at rho=p=1, vA=1; exact to O(epsilon²)."""
    k = jn_zeros(1, 1)[0]/mesh.radius
    speed = np.sqrt(gamma+1)
    omega = k*speed
    r, z = mesh.corners
    b0 = np.sqrt(MU0)
    psi = .5*b0*r*r+epsilon*b0*r*j1(k*r)/k*np.cos(omega*time)
    psi[-1] = .5*b0*mesh.radius**2
    mode = radial_average(mesh, lambda rr: j0(k*rr))[:, None]
    velocity = epsilon*speed*radial_average(mesh, lambda rr: j1(k*rr))[:, None]*np.sin(omega*time)
    state = from_primitive(mesh, psi, 1+epsilon*mode*np.cos(omega*time),
                           1+gamma*epsilon*mode*np.cos(omega*time), vr=velocity, gamma=gamma)
    return state, omega


def resistive_bessel_mode(mesh, *, amplitude=1e-5, gamma=5/3):
    k = jn_zeros(1, 1)[0]/mesh.radius
    r, z = mesh.corners
    psi = amplitude*r*j1(k*r)/k
    psi[[0, -1]] = 0
    return from_primitive(mesh, psi, 1.0, .1, gamma=gamma), k


def advected_flux(mesh, *, amplitude=1e-7, velocity=1.0, time=0.0):
    r, z = mesh.corners
    psi = amplitude*r*r*(1-r*r/mesh.radius**2)**3*(1+.5*np.cos(np.pi*(z-velocity*time)/mesh.half_length))
    return from_primitive(mesh, psi, 1.0, 1.0, vz=velocity)


def sod_initial(mesh):
    r, z = mesh.centers
    return from_primitive(mesh, np.zeros((mesh.nr+1, mesh.nz+1)),
                           np.where(z < 0, 1, .125), np.where(z < 0, 1, .1), gamma=1.4)


def radial_pressure_balance(mesh, *, balanced=True):
    """Straight axial field with known J×B and optional exact pressure balance."""
    r, _ = mesh.corners
    b0, alpha = np.sqrt(MU0), .5/mesh.radius**2
    psi = b0*(r*r/2-alpha*r**4/4)
    _, bz = mesh.cell_fields(psi)
    pressure = 2-bz*bz/(2*MU0) if balanced else 2.0
    return from_primitive(mesh, psi, 1.0, pressure), b0, alpha


def sod_exact(position, time, gamma=1.4):
    """Exact self-similar Sod solution before periodic copies can interfere."""
    if time <= 0:
        raise ValueError("Exact Sod solution requires positive time")
    rho_l, p_l, rho_r, p_r = 1.0, 1.0, .125, .1
    a_l, a_r = np.sqrt(gamma*p_l/rho_l), np.sqrt(gamma*p_r/rho_r)
    def pressure_function(p, rho, pk):
        a = np.sqrt(gamma*pk/rho)
        if p > pk:
            return (p-pk)*np.sqrt(2/((gamma+1)*rho)/(p+(gamma-1)/(gamma+1)*pk))
        return 2*a/(gamma-1)*((p/pk)**((gamma-1)/(2*gamma))-1)
    p_star = brentq(lambda p: pressure_function(p, rho_l, p_l)+pressure_function(p, rho_r, p_r), p_r, p_l, xtol=1e-14)
    v_star = .5*(pressure_function(p_star, rho_r, p_r)-pressure_function(p_star, rho_l, p_l))
    a_star = a_l*(p_star/p_l)**((gamma-1)/(2*gamma))
    shock = a_r*np.sqrt((gamma+1)/(2*gamma)*p_star/p_r+(gamma-1)/(2*gamma))
    xi = np.asarray(position)/time
    rho, v, p = np.full_like(xi, rho_r), np.zeros_like(xi), np.full_like(xi, p_r)
    left = xi < -a_l
    fan = (xi >= -a_l) & (xi < v_star-a_star)
    left_star = (xi >= v_star-a_star) & (xi < v_star)
    right_star = (xi >= v_star) & (xi < shock)
    rho[left], p[left] = rho_l, p_l
    sound = 2/(gamma+1)*(a_l-(gamma-1)*xi[fan]/2)
    rho[fan] = rho_l*(sound/a_l)**(2/(gamma-1))
    v[fan] = 2/(gamma+1)*(a_l+xi[fan])
    p[fan] = p_l*(sound/a_l)**(2*gamma/(gamma-1))
    rho[left_star] = rho_l*(p_star/p_l)**(1/gamma)
    p[left_star], v[left_star] = p_star, v_star
    ratio = (gamma-1)/(gamma+1)
    rho[right_star] = rho_r*(p_star/p_r+ratio)/(ratio*p_star/p_r+1)
    p[right_star], v[right_star] = p_star, v_star
    return rho, v, p
