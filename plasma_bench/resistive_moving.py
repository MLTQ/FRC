"""Physical resistivity and Ohmic heat composed with actual material dynamics."""
from copy import copy
import numpy as np
from .moving_plasma import MovingPlasma
from .material_resistivity import MaterialResistivity
from .resistive_substeps import integrate_diffusion


class ResistiveMovingPlasma(MovingPlasma):
    """Scalar-eta Strang reference; still no physical transport closure or remapping."""
    def __init__(self,domain,initial_a,*,resistivity,diffusion_tolerance=1e-4,**kwargs):
        if not np.isscalar(resistivity) or not np.isfinite(resistivity) or resistivity < 0:
            raise ValueError("Finite nonnegative scalar physical resistivity required")
        self.resistivity=float(resistivity)
        if not np.isfinite(diffusion_tolerance) or not 0 < diffusion_tolerance < 1:
            raise ValueError("Diffusion accuracy tolerance must be finite and in (0,1)")
        self.diffusion_tolerance=float(diffusion_tolerance)
        self.diffusion_substeps=self.diffusion_rejected_trials=0
        self.maximum_diffusion_error=0.
        self.cell_ohmic_heat=np.zeros(domain.fluid_cell_count)
        self.maximum_resistive_balance_residual=0.
        self.maximum_resistive_linear_residual=0.
        super().__init__(domain,initial_a,**kwargs)

    def _resist(self,dt):
        operator=MaterialResistivity(self.domain,self.state.coordinates,self.resistivity)
        a=self.flux.coefficients(self.state.coordinates);bias=float(self.bias.value(self.time))
        result=integrate_diffusion(operator,a,dt,bias_start=bias,tolerance=self.diffusion_tolerance)
        labels=self.flux.labels_for_coefficients(self.state.coordinates,result.material_a)
        volume=self.domain.volumes(operator.elements)
        invariant=self.pressure_invariant+(self.gamma-1)*result.cell_heat*volume**(self.gamma-1)
        if not np.all(np.isfinite(invariant)) or not np.all(np.isfinite(invariant/volume**self.gamma)):
            raise ValueError("Nonfinite thermal proposal; accepted state retained")
        self.flux.labels=labels;self.pressure_invariant=invariant
        self.cell_ohmic_heat+=result.cell_heat
        self.maximum_resistive_balance_residual=max(self.maximum_resistive_balance_residual,result.maximum_balance_residual)
        self.maximum_resistive_linear_residual=max(self.maximum_resistive_linear_residual,result.maximum_linear_residual)
        self.diffusion_substeps+=result.accepted_substeps
        self.diffusion_rejected_trials+=result.rejected_trials
        self.maximum_diffusion_error=max(self.maximum_diffusion_error,result.maximum_accepted_error)

    def step(self,dt):
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError("Finite positive material step required")
        if self.resistivity == 0:
            return super().step(dt)
        draft=copy(self)
        draft.state=self.state.copy();draft.flux=copy(self.flux)
        draft.flux.labels=self.flux.labels.copy()
        draft.pressure_invariant=self.pressure_invariant.copy()
        draft.cell_ohmic_heat=self.cell_ohmic_heat.copy()
        # In extended time, resistivity holds t/bias fixed; ideal evolution advances t.
        draft._resist(.5*dt)
        MovingPlasma.step(draft,dt)
        draft._resist(.5*dt)
        for name in ("state","flux","pressure_invariant","time","steps","source_work",
                     "cell_ohmic_heat","maximum_resistive_balance_residual","maximum_resistive_linear_residual",
                     "diffusion_substeps","diffusion_rejected_trials","maximum_diffusion_error"):
            setattr(self,name,getattr(draft,name))

    def diagnostics(self):
        result=super().diagnostics()
        result.update(physical_resistivity_ohm_m=self.resistivity,
                      ohmic_heat_j=float(self.cell_ohmic_heat.sum()),
                      minimum_cell_ohmic_heat_j=float(self.cell_ohmic_heat.min()),
                      maximum_resistive_substep_balance_residual_j=self.maximum_resistive_balance_residual,
                      maximum_resistive_linear_relative_residual=self.maximum_resistive_linear_residual,
                      diffusion_tolerance=self.diffusion_tolerance,
                      diffusion_accepted_substeps=self.diffusion_substeps,
                      diffusion_rejected_trials=self.diffusion_rejected_trials,
                      maximum_accepted_diffusion_error=self.maximum_diffusion_error)
        return result
