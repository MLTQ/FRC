# Original nonsingular Euler reference

`fields_x` independently states the same manufactured fields using x=r²/2,
u=vr/r and Omega=vtheta/r, without division by zero or square-root derivatives
at the axis. Complex-step derivatives enter the original material equations:
div(v)=2u+2x u_x+vz_z; Dt u=−u²+Omega²−p_x/(Mn);
Dt Omega=−2u Omega; Dt vz=−p_z/(Mn); Dt sa=0. These supply `parent_rate`.

`regular_source` is the exact time derivative minus that continuous parent
rate; it is verification forcing only. `initial` zeroes the explicitly inactive
axis radial/angular slots, which are not independent physical DOFs.
`energy` is an independent original EOS/kinetic formula supporting complex
differentiation; inactive slots contribute exactly zero because x=0.
