"""Conservative axial transport of the full reacting species moment state."""
import numpy as np
from .species_euler_flux import hllc


class AtomicTransport:
    """Periodic axial finite volumes, fixed Nref, and explicit local source."""
    def __init__(self, source, cells, length, reference_density, *, cfl=.3):
        if (not isinstance(cells,int) or cells<4 or not np.isfinite([length,reference_density,cfl]).all()
                or length<=0 or reference_density<=0 or not 0<cfl<=.4):
            raise ValueError('At least four cells, positive physical scales and 0<CFL<=0.4 required')
        self.source=source;self.cells=cells;self.length=float(length);self.density=float(reference_density)
        self.dx=length/cells;self.x=(np.arange(cells)+.5)*self.dx;self.cfl=cfl;self.gamma=5/3
        self.number=np.r_[np.arange(source.network.count),source.network.electrons]
        self.mass_ratio=source.mass/source.neutral_mass
        self.frame=np.asarray(getattr(source,'reference_velocity',np.zeros(3)))/source.speed_scale

    def primitive(self,state):
        if state.shape!=(self.cells,self.source.size):raise ValueError('Complete spatial species state required')
        fields=[self.source.fields(row) for row in state]
        return np.stack([np.column_stack((f['numbers'],(self.gamma-1)*f['thermal_ev'],
                          f['velocity']/self.source.speed_scale)) for f in fields])

    def timestep(self,state):
        p=self.primitive(state);number=p[...,0]
        sound=np.sqrt(self.gamma*np.divide(p[...,1],number,out=np.zeros_like(number),where=number>0)/self.mass_ratio)
        speed=self.source.speed_scale*np.max(np.where(number>0,abs(p[...,4]+self.frame[2])+sound,0.))
        return self.cfl*self.dx/speed if speed else np.inf

    def derivative(self,state):
        p=self.primitive(state)
        backward=p-np.roll(p,1,axis=0);forward=np.roll(p,-1,axis=0)-p
        slope=np.where(backward*forward>0,np.sign(backward)*np.minimum(abs(backward),abs(forward)),0.)
        left=p+.5*slope;right=np.roll(p-.5*slope,-1,axis=0)
        flux=hllc(left,right,self.mass_ratio,self.frame[2],gamma=self.gamma)*self.source.speed_scale
        rate=-(flux-np.roll(flux,1,axis=0))/self.dx
        result=np.zeros_like(state);result[:,self.number]=rate[...,0]
        result[:,self.source.momentum]=rate[...,1:4];result[:,self.source.energy]=rate[...,4]
        return result

    def transport_step(self,state,dt):
        if not np.isfinite(dt) or dt<=0:raise ValueError('Positive finite transport step required')
        first=state+dt*self.derivative(state)
        result=.5*state+.5*(first+dt*self.derivative(first))
        self.primitive(result)
        return result

    def step(self,state,time,dt,*,reactions=True,heat=None):
        if not reactions:return self.transport_step(state,dt)
        def source_step(values,start,duration):
            initial=np.zeros(self.cells) if heat is None else np.broadcast_to(heat(self.x,start),(self.cells,))
            midpoint=np.zeros(self.cells) if heat is None else np.broadcast_to(heat(self.x,start+duration/2),(self.cells,))
            return np.asarray([self.source.step(row,duration,self.density,initial_heat=initial[k],midpoint_heat=midpoint[k])
                               for k,row in enumerate(values)])
        first=source_step(state,time,dt/2)
        transported=self.transport_step(first,dt)
        return source_step(transported,time+dt/2,dt/2)

    def totals(self,state):
        source=self.source;n=source.network
        # Integrals per unit area in reference-number/eV/momentum units.
        sums=self.dx*state.sum(axis=0);number=sums[self.number]
        return dict(nuclei=float(number[:-1].sum()),charge=float(n.charge@number[:-1]-number[-1]),
            mass_over_neutral=float(source.mass@number/source.neutral_mass),
            momentum_scaled=sums[source.momentum].sum(axis=0),
            energy_ev=float(sums[source.energy].sum()+n.energy_ev@sums[:n.count]+sums[n.photons]-sums[n.heating]))
