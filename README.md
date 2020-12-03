Code and manuscript for numerical solution of the circular hydraulic jump problem.

## Riemann solvers:

- `rpn2_sw_hllemcc.f90`: Kemm's HLLEMCC solver, meant to avoid carbuncles while
  being less dissipative than HLLE.  It doesn't seem to work at all if the carbuncle
  fix is off (it should be a Roe solver in that case).  This requires more testing.
- `rpn2_swq_hllemcc.f90`: Same as above, but this version works on quadrilateral
  (mapped) grids.
- `rpn2_shallow_es.f90`: New entropy-dissipative solver.  Works on mapped grids.

### Transverse solvers:

- `rpt2_dummy.f90`: Copy of generic dummy solver from Clawpack, only to be used with dimensional splitting.
- `rpt2_shallow_roe_with_efix.f90`: Standard Clawpack transverse solver for Cartesian grids.
- `rpt2_shallow_roe_mapped`: Updated version of Clawpack 4 transverse solver for mapped grids.  We should add this to clawpack/riemann after it is fully tested.

### Old ones that we're not actually using:
- `rpn2_shallow_hllc.f90`
- `rpn2swq-hllemcc.f90`
- `rpn2swq_hllemcc.f90`
- `rpn2swq_hllemccroef.f90`
- `rpt2_swq.f90`
- `rpt2_shallow_roe_with_efix_annulus.f90`



## To-do list

- [X] Test mapped grid with circular boundaries
- [ ] Test mapped grid with perturbed nodes
- Run high-resolution simulations with combinations of the following options
    - [X] Riemann solvers (Roe, HLLE, HLLEMCC, HLLC, GeoClaw, others?)
    - [ ] Jump forcing conditions (bathymetry, fixed outflow, friction)
    - [X] Grids (Cartesian, circle-ish, perturbed)
    - [ ] 2nd-order LW / High-order WENO
