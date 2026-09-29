"""Independent exact common-velocity contact and nonlinear Euler simple wave."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import atomic_mass, elementary_charge as E, mu_0 as MU0


@dataclass(frozen=True)
class WarmReference:
    kind: str = 'contact'
    length: float = .24
    number_density: float = 1e20
    ion_mass: float = 39.948*atomic_mass
    electron_ev: float = 5.
    ion_ev: float = 1.
    gamma: float = 5/3
    bias: float = .01
    density_amplitude: float = .25
    temperature_amplitude: float = .12
    velocity: float = 2400.
    wave_mach: float = .05

    @property
    def rho0(self):
        return self.ion_mass*self.number_density

    @property
    def pe0(self):
        return self.number_density*E*self.electron_ev

    @property
    def p0(self):
        return self.number_density*E*(self.electron_ev+self.ion_ev)

    @property
    def sound_speed(self):
        return np.sqrt(self.gamma*self.p0/self.rho0)

    @property
    def period(self):
        return self.length/(abs(self.velocity) if self.kind == 'contact' else self.sound_speed)

    def fields(self, z, time):
        z = np.asarray(z)
        k = 2*np.pi/self.length
        if self.kind == 'contact':
            phase = k*(z-self.velocity*time)
            rho = self.rho0*(1+self.density_amplitude*np.cos(phase))
            pe = self.pe0*(1+self.temperature_amplitude*np.sin(phase+.23))
            p = np.full_like(rho, self.p0)
            v = np.full_like(rho, self.velocity)
        elif self.kind == 'simple_wave':
            c0, amplitude = self.sound_speed, self.wave_mach*self.sound_speed
            a = .5*(self.gamma+1)
            if abs(time)*abs(a*amplitude*k) >= 1:
                raise ValueError('Exact simple wave is single valued only before characteristic crossing')
            foot = z-c0*time
            for _ in range(20):
                step = (foot+c0*time+a*amplitude*time*np.sin(k*foot)-z)/(1+a*amplitude*time*k*np.cos(k*foot))
                foot = foot-step
                if np.max(abs(step)) < 2e-16*self.length:
                    break
            if np.max(abs(foot+c0*time+a*amplitude*time*np.sin(k*foot)-z)) > 1e-13*self.length:
                raise RuntimeError('Characteristic inversion did not converge')
            v = amplitude*np.sin(k*foot)
            c = c0+.5*(self.gamma-1)*v
            rho = self.rho0*(c/c0)**(2/(self.gamma-1))
            p = self.p0*(rho/self.rho0)**self.gamma
            pe = self.pe0*(rho/self.rho0)**self.gamma
        else:
            raise ValueError('Reference kind must be contact or simple_wave')
        energy = p/(self.gamma-1)+.5*rho*v*v+self.bias**2/(2*MU0)
        return dict(density=rho, velocity=v, pressure=p, electron_pressure=pe,
                    ion_pressure=p-pe, momentum=rho*v, energy=energy)

    def cell_means(self, mesh, time=0., order=12):
        nodes, weights = np.polynomial.legendre.leggauss(order)
        fields = self.fields(mesh.z[:, None]+.5*mesh.dz*nodes, time)
        return {name: np.broadcast_to(.5*np.sum(value*weights, axis=-1)[None], (mesh.nr, mesh.nz)).copy()
                for name, value in fields.items()}
