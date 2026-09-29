# Axisymmetric linear Hall cylinder

`HallCylinderMesh` defines a uniform radial piecewise-linear mesh and periodic
odd axial Fourier grid. `LinearHallCylinder` assembles weighted radial mass and
vector-potential stiffness, numerically diagonalizes them, and delegates all
coupled evolution to `LinearHallModes`. No analytic Bessel eigenvalue is used.
All radial degrees of freedom and all resolved axial frequencies are retained.
The input spectrum must be resolved: the sampler cannot diagnose aliasing.

`from_nodal` takes real axis-regular a=psi/r, Btheta, velocity stream u and
vtheta, with exact zero traces at both radial endpoints. Fields are interpolated,
not silently projected onto a known reference eigenmode. `sample` reconstructs
Br/Btheta/Bz, vr/vtheta/vz or weak Jr/Jtheta/Jz in SI, including the regular
axis limit. Sampling chooses the outer element at an interior radial node;
first-derivative fields are discontinuous there. `quadrature_points` returns
physical 2 pi r dr dz weights for independent field/current/heat checks.
`angular_momentum` integrates rho r vtheta; it is conserved in this reference.
`step` and `energies` check real-field Fourier conjugacy and leave caller states
unchanged. Returned heat is total positive collisional dissipation, with no
thermal feedback. State coefficients are not existing compressible MHD states.

Boundary traces produce a restricted zero-flux perturbation sector with zero
energy flux. They are neither actual free vacuum nor a claimed perfect conductor
with eta>0. See LINEAR_HALL_NUMERICS for the precise conditions. This linear
incompressible component cannot be installed as a nonlinear formation model.
