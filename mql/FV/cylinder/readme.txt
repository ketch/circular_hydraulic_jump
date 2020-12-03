***************************
***** TO RUN THE CODE *****
***************************
* Set use_petsc to either True or False
* Change the parameters from main (not on the definition of setup)
* Choose riemann_solver to be either 'es' or 'hllemcc'
* If riemann_solver='es' then
   * To use the blended solver set: use_dmin_blended=1.0 and set_Ri=None
   * To use Roe's solver set: use_dmin_blended=0.0 and set_Ri=0.0
   * To use Rusanov's solver set: use_dmin_blended=0.0 and set_Ri=1.0

*******************************
***** TO PLOT THE RESULTS *****
*******************************
* Run plot.py
* Set use_petsc and plot_schlieren to either True or False
