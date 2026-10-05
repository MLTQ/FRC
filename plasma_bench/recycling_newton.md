# Reusable full Gauss Jacobian products

Uses exactly the centered physical-residual perturbation size and GMRES
settings of the frozen original Newton fallback. Each already computed Jv
is also passed to a limited-memory inverse proposal as (v,Jv). No extra RHS,
modified Jacobian, altered physical equation or relaxed tolerance is introduced.
Unlike tiny accepted-iterate differences, these products use the original
cuberoot(machine epsilon) finite-difference scale. Their truncation error
remains controlled by acceptance against the actual nonlinear residual.

The caller owns local proposal memory and must commit it only after successful
physical stage/endpoint acceptance. GMRES uses the original fixed approximate
inverse; changing memory during matvec does not alter its operator. This module
only proposes directions and preserves the original failure contract.
