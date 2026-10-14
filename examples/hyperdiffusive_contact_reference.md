# Independent continuum effect of numerical damping

Reuse the original Fourier/Bessel contact operator and its frozen initial mode.
For its analytic curl matrix C, the added magnetic block before scaling is
-eta6/mu0 (I+C mu/(e^2 n) C/mu0)^-1 C^6. This follows directly from
E6=eta6 curl^4 J, J=curl B/mu0 and Faraday, including variable-density current
inertia. No production mesh, mass, curl, thermal projection or eigenmode fit
enters this calculation. Q6 is quadratic, so gives no linear entropy source.

Evaluate the matrix exponential at the declared duration, starting from the
unchanged collisionless continuum mode. Compare modified and original fields
component by component; the shared radial Bessel norm cancels in each ratio.
Also expose the original exponential versus exp(lambda t) consistency error.
This measures resolved-mode distortion in the continuum equations, not actual
discrete stability, nonlinear accuracy or convergence. Save both operators,
initial/final vectors and existing reference arrays. No recursive replay or
new verification framework is introduced.
