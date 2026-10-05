"""Limited-memory inverse proposals from accepted full-residual secants."""
import numpy as np


class GaussSecant:
    def __init__(self, inverse, memory=16):
        if type(memory) is not int or memory < 1:
            raise ValueError('Positive secant memory required')
        self.inverse = inverse; self.memory = memory; self.steps = []; self.images = []

    def append(self, step, residual_change):
        image = self.inverse(residual_change)
        length = np.linalg.norm(image)
        if length == 0 or not np.isfinite(length): return
        self.steps.append(step.ravel()/length); self.images.append(image.ravel()/length)
        self.steps = self.steps[-self.memory:]; self.images = self.images[-self.memory:]

    def correction(self, error):
        value = self.inverse(error)
        if not self.images: return -value, 0
        images = np.column_stack(self.images); steps = np.column_stack(self.steps)
        coefficients, _, rank, _ = np.linalg.lstsq(images, value.ravel(), rcond=1e-10)
        return (-value.ravel()-(steps-images)@coefficients).reshape(error.shape), int(rank)
