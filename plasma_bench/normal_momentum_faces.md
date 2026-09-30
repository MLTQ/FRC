# Centered normal momentum faces

Shared normal faces need a stencil centered on the FACE. Reusing an arbitrary-
point cell reconstruction that averages two cell-containing stencils biases the
normal trace and can create an exponentially growing periodic advection mode.
This module uses one centered degree-3/5 true-moment fit at each normal face.

Periodic axial faces are addressed by integer index, avoiding floating-point
floor decisions at exact cell boundaries. The resulting constant-translation
divergence is skew-symmetric in the uniform finite-volume metric. Degree 5 gives
the classical sixth-order centered derivative, with purely imaginary Fourier
frequencies. RK4 still has its explicit imaginary-axis stability restriction;
this is not a nonlinear positivity or shock-stability theorem.

Radial faces use r^power P(x), x=r²/2, and annular moments, clipping stencils to
the physical interior only. The caller enforces its impermeable axis/wall flux.
This module changes normal momentum traces, not the mean-preserving density
polynomials or their continuous normal-flux primitive. No diffusion is added.
