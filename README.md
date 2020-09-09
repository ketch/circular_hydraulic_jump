Code and manuscript for numerical solution of the circular hydraulic jump problem.

## Riemann solvers:

- `rpn2_sw_hllemcc.f90`: Kemm's HLLEMCC solver, meant to avoid carbuncles while
  being less dissipative than HLLE.  It doesn't seem to work at all if the carbuncle
  fix is off (it should be a Roe solver in that case).  This requires more testing.
- `rpn2swq-hllemcc.f90`: Same as above, but this version works on quadrilateral
  (mapped) grids.  Not tested yet.
