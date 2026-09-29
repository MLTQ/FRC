# Independent common-coordinate sampler checks

Two distinct radial and angular grids carry a known affine deformation, affine
velocity, linear scalar magnetic coefficient and linear current-over-radius.
At common physical reference points, compare positions, velocity, analytic B
and current against closed formulas. All reference vertices and the axis also
reconstruct exactly without an inverse center Jacobian. Outside-polygon points
raise instead of silently extrapolating. Constant P0 heat checks units and cell
lookup. These tests check the comparison machinery, not the MHD solution.
The affine cylindrical-volume ratio is independently .9²*1.1 and is also checked.
