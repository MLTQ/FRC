"""Causal positive-cage pulse control from observed ion-cloud motion."""
import numpy as np


class CageDrive:
    def __init__(self, config):
        self.config = config
        self.mode = config.get('drive_mode', 'legacy')
        if self.mode not in ('legacy', 'dc', 'fixed', 'hops'):
            raise ValueError('Unknown cage drive')
        if not (0 < config.get('release_fraction', .8) < config.get('return_fraction', 1.)):
            raise ValueError('Feedback release radius must be below return radius')
        for key in ('response_ns', 'sensor_ns', 'min_dwell_ns', 'warmup_ns'):
            if config.get(key, 0) < 0:
                raise ValueError('Control times must be nonnegative')
        if config['bias_v'] < 0 or config.get('pulse_period_ns', 1000) <= 0:
            raise ValueError('Positive voltage and positive fixed period required')
        if not 0 < config.get('pulse_duty', .5) < 1:
            raise ValueError('Pulse duty must be between zero and one')
        self.value = 0.
        self.gate = True
        self.reference = None
        self.filtered = None
        self.last_switch = 0.
        self.last_on = None
        self.last_period = None
        self.pulses = 1
        self.switches = []
        self.last_observation = None

    def legacy(self, t):
        c = self.config
        ramp = 1. if c['rise_ns'] == 0 else 1-np.exp(-t/(c['rise_ns']*1e-9))
        return c['bias_v']*ramp*(1+c['modulation']*np.sin(2*np.pi*c['frequency_mhz']*1e6*t))

    def initial(self):
        self.value = self.legacy(0.) if self.mode == 'legacy' else 0.
        return self.value

    def advance(self, t, dt, radius, radial_speed, ion_count):
        """Observe the old state, choose a target, then advance voltage to t+dt."""
        if self.mode == 'legacy':
            self.value = self.legacy(t+dt)
            return self.value
        c = self.config
        if radius is not None:
            if self.filtered is None:
                self.filtered = radius
            tau = c['sensor_ns']*1e-9
            blend = 1. if tau == 0 else -np.expm1(-dt/tau)
            self.filtered += blend*(radius-self.filtered)
        warmup = c['warmup_ns']*1e-9
        desired = True
        if t >= warmup:
            if self.reference is None:
                self.reference = self.filtered
            if self.mode == 'fixed':
                phase = ((t-warmup)/(c['pulse_period_ns']*1e-9)) % 1.
                desired = phase < c['pulse_duty']
            elif self.mode == 'hops':
                desired = self.gate
                valid = ion_count >= c['min_control_ions'] and self.reference is not None
                ready = t-self.last_switch >= c['min_dwell_ns']*1e-9
                if valid and ready:
                    if self.gate and self.filtered <= c['release_fraction']*self.reference and radial_speed < 0:
                        desired = False
                    elif not self.gate and self.filtered >= c['return_fraction']*self.reference and radial_speed > 0:
                        desired = True
        if desired != self.gate:
            self.gate = desired
            self.last_switch = t
            if desired:
                self.pulses += 1
                if self.last_on is not None:
                    self.last_period = t-self.last_on
                self.last_on = t
            self.switches.append(dict(time_us=t*1e6, on=desired, radius_m=radius,
                                      filtered_radius_m=self.filtered, radial_speed_m_s=radial_speed))
        target = c['bias_v'] if self.gate else 0.
        tau = c['response_ns']*1e-9
        self.value += (target-self.value)*(1. if tau == 0 else -np.expm1(-dt/tau))
        return self.value

    def diagnostics(self):
        return dict(mode=self.mode, gate_on=self.gate, voltage_v=float(self.value),
                    filtered_radius_m=self.filtered, reference_radius_m=self.reference,
                    completed_pulse_intervals=max(0, self.pulses-2),
                    measured_period_us=None if self.last_period is None else self.last_period*1e6,
                    measured_frequency_mhz=None if self.last_period is None else 1e-6/self.last_period,
                    switches=len(self.switches))
