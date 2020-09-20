The subdirectories of this directory contain plots from PyClaw simulations.
The directory name encodes the parameters used.  The first part is:

- Claw1: Classic Clawpack with no second-order corrections
- Claw2: Classic Clawpack with second-order corrections (Lax-Wendroff-LeVeque)
- Sharpclaw: Sharpclaw with WENO5 and SSP(10,4)

The second part is the name of the Riemann solver used.
The third part is the resolution.
Finally, the values of any other parameters that were changed
are given.  "epsilon" refers to the parameter in Kemm's HLLEMC Riemann solver.
