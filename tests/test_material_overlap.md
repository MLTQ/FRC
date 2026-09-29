# Independent overlap integration checks

Checks a triangle/rectangle clipped area and first radial moment against exact
rational integrals. Integrates a discontinuous indicator whose jump lies on an
unrelated mesh edge, against its exact cylindrical-volume integral. Two distinct
graded/radial/angular material grids then cover the entire smaller polygon;
their common integral is checked against independent conical-frustum volume.
These tests expose cell-edge sampling errors that uniform Gauss nodes can miss.
