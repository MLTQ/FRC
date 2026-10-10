# Original term capture without changing a solve

Loads an accepted signed full-RHS archive, verifies its raw hash and bitwise
reconstructed initial fixture, then evaluates the original bulk, Lorentz,
current transport, species pressure/entropy and wall reaction terms at those
exact states. Their original sum must reproduce the saved complete fluid rate
bitwise. The full original Ohm equation is checked against saved accepted Jt/E
under the unchanged tolerance. No new induction solve or fitted state is used.

Saves all term coefficients: fluid contributions, weak electric forces and
mass-inverted electric coefficients, including physical inertia from saved Jt.
Full and reduced pressure plus its exact potential are distinct. Retains all
coefficients of the archived states and rates for physical reconstruction.
The wall reaction is separately visible; no individual force is silently clipped.

Outputs include source/input hashes, full-equation checks and immutable raw
arrays. Fresh replay is required. A decomposition is diagnostic and does not
modify production operators, accepted rates or any evolution source. Potential
curl cancellation and reconstruction roundoff must be reported explicitly.
