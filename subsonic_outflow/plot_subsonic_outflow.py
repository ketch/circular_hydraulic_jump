import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
from matplotlib import rc
import numpy as np
import os
from scipy.integrate import ode
from scipy import integrate
from scipy.interpolate import interp1d
from subsonic_outflow import *

def plot_q(frame,
           file_format='petsc',
           file_prefix='claw',
           path='./_output/'):
    import sys
    sys.path.append('.')

    sol=Solution(frame,file_format=file_format,read_aux=False,path=path,file_prefix=file_prefix)
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers
    mx=len(x); my=len(y)

    h=sol.state.q[0,:,:]
    hu=sol.state.q[1,:,:]

    if frame < 10:
        str_frame = "000"+str(frame)
    elif frame < 100:
        str_frame = "00"+str(frame)
    elif frame < 1000:
        str_frame = "0"+str(frame)
    else:
        str_frame = str(frame)
    #
    
    # Plot numerical solution
    pl.plot(x,h[:,my/2],'-k',lw=2)
    pl.plot(r0+dx/2.0,h[0,my/2],'ko',lw=10)

    # plot exact solution in the inner ghost cells
    pl.plot(r_backward,hh_backward,'--r')
    pl.plot(r0-dx/2.0,hm1,'go',lw=10)
    pl.plot(r0-3*dx/2.0,hm2,'go',lw=10)

    # plot exact solution in the domain and outer ghost cells
    pl.plot(r_forward,hh_forward,'--b')

    #
    pl.suptitle("t= "+str(sol.state.t),fontsize=20)        
    pl.savefig('./_plots/h_'+str_frame+'_slices.png')
    pl.close()
    
    # compute L1 error
    return dx*np.sum(np.abs(hInterp_forward(x)-h[:,my/2]))
#

if __name__== "__main__":
    use_petsc = False

    if use_petsc:
        from clawpack.petclaw.solution import Solution
        file_format = 'petsc'
        file_prefix = 'claw'
    else:
        from clawpack.pyclaw.solution import Solution
        file_format = 'ascii'
        file_prefix = 'fort'
    #
    
    if not os.path.exists('./_plots'): os.mkdir('_plots')
    from_frame = 0
    to_frame   = 100
    frames=xrange(from_frame,to_frame+1)
    Ncells_interp = 10000
    
    # ********************************************************* #
    # ***** COMPUTE EXACT SOLUTION (INCLUDING BOUNDARIES) ***** #
    # ********************************************************* #
    # read some data from the solution (at any frame)
    sol=Solution(0,file_format=file_format,read_aux=False,path='./_output/',file_prefix=file_prefix)
    x=sol.state.grid.x.centers
    dx=x[1]-x[0]
    h0 =  sol.state.problem_data['h0']
    r0 =  sol.state.problem_data['r0']
    u0 =  sol.state.problem_data['u0']
    g =  sol.state.problem_data['grav']
    rOutflow =  sol.state.problem_data['rOutflow']
    beta = r0*h0*u0
    
    # inner boundary
    dh=ode(steady_rhs_backward,
           steady_rhs_backward_jacobian).set_integrator('vode',
                                                  method='bdf',
                                                  nsteps=1E5)
    dh.set_initial_value(h0,r0).set_f_params(beta).set_jac_params(beta)
    rEnd=r0-3*dx/2.0 # for 2 ghost cells
    dr=1E-6*dx
    r_backward=np.linspace(r0-dr,rEnd,Ncells_interp)
    hh_backward = r_backward*0
    for i in range(len(r_backward)):
        hh_backward[i] = dh.integrate(r_backward[i])[0]
    #    
    hInterp_backward = interp1d(r_backward,hh_backward)
    hm1 = hInterp_backward(r0-dx/2.0)
    hm2 = hInterp_backward(r0-3*dx/2.0)
    
    # domain and outer boundary
    r_forward = np.linspace(r0,rOutflow + 2*dx,Ncells_interp)
    hh_forward = integrate.odeint(steady_rhs,h0,r_forward,args=(beta,g))[:,0]
    hInterp_forward = interp1d(r_forward,hh_forward[:])
    hp1 = hInterp_forward(rOutflow+dx/2.0)
    hp2 = hInterp_forward(rOutflow+3*dx/2.0)
    # ***** END OF COMPUTATION OF EXACT SOLUTION ***** #
    
    ##########################################
    print('**********************')
    print('**********************')
    print('Plotting solution ...')
    for i in frames:
        L1_error = plot_q(frame=i,file_format=file_format,file_prefix=file_prefix)
        print ('frame '+str(i)+' plotted'+'... L1 error: '+str(L1_error))
    #
#


