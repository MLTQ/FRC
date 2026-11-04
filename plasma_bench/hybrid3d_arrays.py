"""Float64 loop fusion for the existing electron-energy difference operator."""
import numpy as np
from numba import njit, prange, literal_unroll


@njit(cache=True, parallel=True)
def compression(fraction, velocity, d, out):
    nx, ny, nz = fraction.shape
    for index in prange(nx):
        i = np.int64(index)
        ip, im = (i+1 if i+1 < nx else 0), (i-1 if i > 0 else nx-1)
        for j in range(ny):
            jp, jm = (j+1 if j+1 < ny else 0), (j-1 if j > 0 else ny-1)
            for k in range(nz):
                kp, km = (k+1 if k+1 < nz else 0), (k-1 if k > 0 else nz-1)
                value = 0.
                value += (fraction[ip,j,k]*velocity[0,ip,j,k]
                          -fraction[im,j,k]*velocity[0,im,j,k])/(2*d[0])
                value += (fraction[i,jp,k]*velocity[1,i,jp,k]
                          -fraction[i,jm,k]*velocity[1,i,jm,k])/(2*d[1])
                value += (fraction[i,j,kp]*velocity[2,i,j,kp]
                          -fraction[i,j,km]*velocity[2,i,j,km])/(2*d[2])
                out[i,j,k] = value


@njit(cache=True, parallel=True)
def transport(u, velocity, links, plus_wall, minus_wall, d, fields, boundary, outgoing):
    nx, ny, nz = u.shape
    for index in prange(nx):
        i = np.int64(index)
        ip, im = (i+1 if i+1 < nx else 0), (i-1 if i > 0 else nx-1)
        for j in range(ny):
            jp, jm = (j+1 if j+1 < ny else 0), (j-1 if j > 0 else ny-1)
            for k in range(nz):
                kp, km = (k+1 if k+1 < nz else 0), (k-1 if k > 0 else nz-1)
                loss, wall = 0., 0.
                for a in literal_unroll((0, 1, 2)):
                    p = (ip,j,k) if a == 0 else ((i,jp,k) if a == 1 else (i,j,kp))
                    m = (im,j,k) if a == 0 else ((i,jm,k) if a == 1 else (i,j,km))
                    v = velocity[a,i,j,k]
                    speed = .5*(v+velocity[a,p[0],p[1],p[2]])
                    previous = .5*(velocity[a,m[0],m[1],m[2]]+v)
                    upstream = u[i,j,k] if speed >= 0. else u[p[0],p[1],p[2]]
                    fields[a,i,j,k] = links[a,i,j,k]*speed*upstream
                    loss += (links[a,i,j,k]*max(speed,0.)
                             +links[a,m[0],m[1],m[2]]*max(-previous,0.))/d[a]
                    wall += plus_wall[a,i,j,k]*max(v,0.)+minus_wall[a,i,j,k]*max(-v,0.)
                outgoing[i,j,k], boundary[i,j,k] = loss, wall


@njit(cache=True, parallel=True)
def hall_power(u, n, current, d, gamma, charge, floor, out):
    nx, ny, nz = n.shape
    for index in prange(nx):
        i = np.int64(index)
        ip, im = (i+1 if i+1 < nx else 0), (i-1 if i > 0 else nx-1)
        for j in range(ny):
            jp, jm = (j+1 if j+1 < ny else 0), (j-1 if j > 0 else ny-1)
            for k in range(nz):
                kp, km = (k+1 if k+1 < nz else 0), (k-1 if k > 0 else nz-1)
                denominator = charge*max(n[i,j,k],floor)
                dx = ((gamma-1)*u[ip,j,k]-(gamma-1)*u[im,j,k])/(2*d[0])
                dy = ((gamma-1)*u[i,jp,k]-(gamma-1)*u[i,jm,k])/(2*d[1])
                dz = ((gamma-1)*u[i,j,kp]-(gamma-1)*u[i,j,km])/(2*d[2])
                value = 0.
                value += current[0,i,j,k]*(-dx/denominator)
                value += current[1,i,j,k]*(-dy/denominator)
                value += current[2,i,j,k]*(-dz/denominator)
                out[i,j,k] = value


@njit(cache=True, parallel=True)
def eta_to_edges(eta, mask_x, mask_y, mask_z, out):
    nx, ny, nz = eta.shape
    for index in prange(nx):
        i = np.int64(index)
        ip = i+1 if i+1 < nx else 0
        for j in range(ny):
            jp = j+1 if j+1 < ny else 0
            for k in range(nz):
                kp = k+1 if k+1 < nz else 0
                out[0,i,j,k] = .5*(eta[i,j,k]+eta[ip,j,k])*mask_x[i,j,k]
                out[1,i,j,k] = .5*(eta[i,j,k]+eta[i,jp,k])*mask_y[i,j,k]
                out[2,i,j,k] = .5*(eta[i,j,k]+eta[i,j,kp])*mask_z[i,j,k]


@njit(cache=True, parallel=True)
def joule_heat(currents, eta, out):
    nx, ny, nz = out.shape
    for index in prange(nx):
        i = np.int64(index)
        im = i-1 if i > 0 else nx-1
        for j in range(ny):
            jm = j-1 if j > 0 else ny-1
            for k in range(nz):
                km = k-1 if k > 0 else nz-1
                value = 0.
                for a in literal_unroll((0, 1, 2)):
                    m = (im,j,k) if a == 0 else ((i,jm,k) if a == 1 else (i,j,km))
                    here = eta[a][i,j,k]*currents[a][i,j,k]**2
                    previous = eta[a][m[0],m[1],m[2]]*currents[a][m[0],m[1],m[2]]**2
                    value += .5*(here+previous)
                out[i,j,k] = value


@njit(cache=True, parallel=True)
def energy_rhs(u, fraction, fields, boundary, compression_value, heat, h, gamma, d,
               full, work, outflow):
    nx, ny, nz = u.shape
    for index in prange(nx):
        i = np.int64(index)
        im = i-1 if i > 0 else nx-1
        for j in range(ny):
            jm = j-1 if j > 0 else ny-1
            for k in range(nz):
                km = k-1 if k > 0 else nz-1
                divergence = 0.
                divergence += (fields[0,i,j,k]-fields[0,im,j,k])/d[0]
                divergence += (fields[1,i,j,k]-fields[1,i,jm,k])/d[1]
                divergence += (fields[2,i,j,k]-fields[2,i,j,km])/d[2]
                escaped = gamma*u[i,j,k]*boundary[i,j,k]
                power = -((gamma-1)*u[i,j,k])*compression_value[i,j,k]
                outflow[i,j,k], work[i,j,k] = escaped, power
                full[i,j,k] = u[i,j,k]*fraction[i,j,k]+h*(
                    -divergence-escaped+power+heat[i,j,k])


class EnergyArrays:
    """Per-closure scratch space; results are borrowed until the next call."""
    def __init__(self, electrons):
        shape = electrons.u.shape
        self.d = electrons.d
        self.links = np.asarray(electrons.link_area)
        self.plus_wall = np.asarray(electrons.boundary_plus)
        self.minus_wall = np.asarray(electrons.boundary_minus)
        self.compression_value = np.empty(shape)
        self.fields = np.empty((3,)+shape)
        self.boundary, self.outgoing = np.empty(shape), np.empty(shape)
        self.full, self.work, self.outflow = (np.empty(shape) for _ in range(3))
        self.edge_resistivity = np.empty((3,)+shape)

    def pressure_power(self, u, n, current, gamma, charge, floor):
        # RK stages retain all three work arrays until their weighted sum.
        out = np.empty_like(n)
        hall_power(u, n, current, self.d, gamma, charge, floor, out)
        return out

    def edge_eta(self, eta, masks):
        eta_to_edges(eta, *masks, self.edge_resistivity)
        return self.edge_resistivity

    def heating(self, currents, eta):
        # As with pressure work, earlier stage values must remain owned.
        out = np.empty_like(currents[0])
        joule_heat(tuple(currents), tuple(eta), out)
        return out

    def divergence_velocity(self, fraction, velocity):
        compression(fraction, velocity, self.d, self.compression_value)
        return self.compression_value

    def fluxes(self, u, velocity):
        transport(u, velocity, self.links, self.plus_wall, self.minus_wall, self.d,
                  self.fields, self.boundary, self.outgoing)
        return self.fields, self.boundary, self.outgoing

    def rhs(self, u, fraction, fields, boundary, compression_value, heat, h, gamma):
        energy_rhs(u, fraction, np.asarray(fields), boundary, compression_value, heat, h, gamma,
                   self.d, self.full, self.work, self.outflow)
        return self.full, self.work, self.outflow
