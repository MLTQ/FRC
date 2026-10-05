"""Original centered Newton directions with reusable full-Jacobian products."""
import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres


def direction(residual, change, error, approximate_inverse, memory):
    shape = change.shape; size = change.size
    increment = np.cbrt(np.finfo(float).eps)*(1+np.max(abs(change)))
    history = []

    def product(value):
        value = value.reshape(shape); norm = np.max(abs(value))
        if norm == 0: return np.zeros(size)
        h = increment/norm
        plus = residual(change+h*value)[0]; minus = residual(change-h*value)[0]
        image = (plus-minus)/(2*h)
        memory.append(value, image)
        return image.ravel()

    jacobian = LinearOperator((size, size), matvec=product, dtype=float)
    inverse = LinearOperator((size, size), dtype=float,
        matvec=lambda value: approximate_inverse(value.reshape(shape)).ravel())
    correction, info = gmres(jacobian, -error.ravel(), M=inverse,
        rtol=1e-4, atol=0., restart=30, maxiter=5,
        callback=lambda norm: history.append(float(norm)), callback_type='pr_norm')
    if info != 0 or not np.isfinite(correction).all():
        raise RuntimeError(f'Newton correction did not converge (GMRES info={info}); no input modified')
    return correction.reshape(shape), dict(krylov_iterations=len(history), linear_residual_history=history)
