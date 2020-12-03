***************************
***** TO RUN THE CODE *****
***************************
* Set use_petsc=True
* Change the parameters from main (not on the definition of setup)
* Select test_case to be 1, ..., 5. In the paper we use 1 and 5.
* Set rand_inflow (which is a percentage) to introduce a perturbation. In the paper we use either rand_inflow=0.0 or 1.0
* Choose riemann_solver to be either 'es' or 'hllemcc'
* If riemann_solver='es' then
   * To use the blended solver set: use_dmin_blended=1.0 and set_Ri=None
   * To use Roe's solver set: use_dmin_blended=0.0 and set_Ri=0.0
   * To use Rusanov's solver set: use_dmin_blended=0.0 and set_Ri=1.0

*******************************
***** TO PLOT THE RESULTS *****
*******************************
* Run plot.py
* Set use_petsc=True
* Set plot_schlieren to either True or False
* If desired set xlim and ylim
