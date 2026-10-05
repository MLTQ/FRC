# Independent exact-face reconstruction checks

Checks current, divergence and magnetic traces against smooth analytic
polynomials on both sides of radial/periodic axial faces. The histopolant
current is supplied as exact cell moments. A degree-six poloidal flux has
known distinct one-sided degree-five derivative errors: subtract the derivative
of the exact nodal root product from the true derivative. This verifies the
chosen side without copying the sampler's construction. A two-element density
jump must retain both exact values, with no displacement or averaging.
