from clawpack.petclaw.solution import Solution
#from petclaw.io.petsc import read_petsc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
from matplotlib import rc
import numpy as np
import os
from ent_residual import *
from scipy.integrate import ode
from scipy import integrate
from scipy.interpolate import interp1d

def plot_q(frame,
           with_jump=True,
           file_prefix='claw',
           path='./_output/'):
    import sys
    sys.path.append('.')

    sol=Solution(frame,file_format='petsc',read_aux=False,path=path,file_prefix=file_prefix)
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers
    yy,xx = np.meshgrid(y,x)
    mx=len(x); my=len(y)

    h=sol.state.q[0,:,:]
    hu=sol.state.q[1,:,:]
    q1=sol.state.q[1,:,:]
    Ri = get_ent_residual(sol.state)

    if frame < 10:
        str_frame = "000"+str(frame)
    elif frame < 100:
        str_frame = "00"+str(frame)
    elif frame < 1000:
        str_frame = "0"+str(frame)
    else:
        str_frame = str(frame)
    #

    #fig, axs = pl.subplots(2,1)
    #axs[0].plot(x,h[:,my/2],'--k',lw=2)
    #axs[1].plot(x,Ri[:,my/2],'-k',lw=2)
    
    ###################33
    # ***** COMPUTE EXACT SOLUTION TO USE AT BOUNDARIES ***** #    
    dx=x[1]-x[0]
    h0 =  sol.state.problem_data['h0']
    r0 =  sol.state.problem_data['r0']
    h0 =  sol.state.problem_data['h0']
    u0 =  sol.state.problem_data['u0']
    g =  sol.state.problem_data['grav']
    rOutflow =  sol.state.problem_data['rOutflow']
    beta = r0*h0*u0

    #import pdb; pdb.set_trace()
    #sol.state.get_qbc_from_q(2,sol.state.qbc)
    pl.plot(x,h[:,my/2],'-k',lw=2)
    #pl.plot(r0+dx/2.0,h[0,my/2],'ko',lw=10)
    
    # inner boundary
    dh=ode(steady_rhs_v2,
           steady_rhs_v2_jacobian).set_integrator('vode',
                                                  method='bdf',
                                                  nsteps=1E5)
    dh.set_initial_value(h0,r0).set_f_params(beta).set_jac_params(beta)
    rEnd=r0-3*dx/2.0 # for 2 ghost cells
    dr=1E-3*dx
    r_backward=np.linspace(r0-dr,rEnd,1000)
    hh_backward = r_backward*0
    for i in range(len(r_backward)):
        hh_backward[i] = dh.integrate(r_backward[i])[0]
    #    
    hInterp = interp1d(r_backward,hh_backward)
    hm1 = hInterp(r0-dx/2.0)
    hm2 = hInterp(r0-3*dx/2.0)

    #pl.plot(r_backward,hh_backward,'--r')
    #
    # outer boundary
    r_forward = np.linspace(r0,1.1*rOutflow,1000)
    hh_forward = integrate.odeint(steady_rhs,h0,r_forward,args=(beta,g))[:,0]
    hInterp = interp1d(r_forward,hh_forward[:])
    hp1 = hInterp(rOutflow+dx/2.0)
    hp2 = hInterp(rOutflow+3*dx/2.0)

    #pl.plot(r_forward,hh_forward,'--b')
    
    #import pdb; pdb.set_trace()
    #pl.plot(r0+dx/2.0,h0,'go',lw=10)
    #pl.plot(r0-dx/2.0,hm1,'go',lw=10)
    #pl.plot(r0-3*dx/2.0,hm2,'go',lw=10)
    # END OF COMPUTATION OF EXACT SOLUTION AT BOUNDARIES #    
    #######################
    
    pl.suptitle("t= "+str(sol.state.t),fontsize=20)
    pl.savefig('./_plots_chj_1D/h_'+str_frame+'_slices.png')
    pl.close()
#

def plot_frame(frame,
               file_prefix='claw',
               path='./_output/'):
    import sys
    sys.path.append('.')

    sol=Solution(frame,file_format='petsc',read_aux=False,path=path,file_prefix=file_prefix)
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers
    yy,xx = np.meshgrid(y,x)
    mx=len(x); my=len(y)

    h=sol.state.q[0,:,:]

    if frame < 10:
        str_frame = "000"+str(frame)
    elif frame < 100:
        str_frame = "00"+str(frame)
    elif frame < 1000:
        str_frame = "0"+str(frame)
    else:
        str_frame = str(frame)
    #
    pl.plot(x,h[:,my/2],lw=2)
#

def steady_rhs_v2(r,h,beta):
    g=1.0
    a=g/beta**2 * r**3
    b=r
    return h/(a * h**3 - b)
#

def steady_rhs_v2_jacobian(r,h,beta):
    g=1.0
    a=g/beta**2 * r**3
    b=r
    return (-2*a*h**3-b)/(a*h**3-b)**2
#


def steady_rhs(h,r,beta,g=1.):
    return h/(g/beta**2 * r**3 * h**3 - r)
#

def get_L1_error(frame,
                 color='-r',
                 plot_exact=True,
                 with_jump=True,
                 file_prefix='claw',
                 path='./_output/'):
    sol=Solution(frame,file_format='petsc',read_aux=False,path=path,file_prefix=file_prefix)
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers
    mx=len(x); my=len(y)

    h=sol.state.q[0,:,:]
    hu=sol.state.q[1,:,:]
    
    from scipy import integrate
    r0 = sol.state.problem_data['r0']
    h0 = sol.state.problem_data['h0']
    u0 = sol.state.problem_data['u0']
    beta = r0*h0*u0
    g = sol.state.problem_data['grav']

    rJump = 0.3
    rOutflow = 1.0
    import hydraulic_jump_1D
    rStab,hStab,uStab,hL = initial_and_boundary_data(r_jump = rJump,
                                                     r_inner = r0,
                                                     r_outer = rOutflow,
                                                     num_cells = 10001,
                                                     g=g,
                                                     h_in=h0,
                                                     u_in=u0)
    h_exact = np.zeros_like(h)
    if with_jump:
        from scipy.interpolate import interp1d
        hInterp = interp1d(rStab,hStab)
        rcInterp = np.minimum(np.maximum(x,r0),rOutflow)
        #hh  = (x<r0)*h0 + (x>=r0)*hInterp(rcInterp)
        hh  = (x<r0)*h0 + (x>=r0)*hInterp(x)
        for i in range(my):
            h_exact[:,i] = hh
        #
    else:
        x_exact = np.linspace(0.1,1.0,10000)
        hh = integrate.odeint(steady_rhs,h0,x_exact,args=(beta,g))
        from scipy.interpolate import interp1d        
        hInterp = interp1d(x_exact,hh[:,0])
        for i in range(my):
            h_exact[:,i] = hInterp(x)
    #

    pl.xticks(fontsize=20)
    pl.yticks(fontsize=20)
    pl.plot(x,h[:,my/2],color,lw=3)
    if plot_exact:
        pl.plot(x,h_exact[:,my/2],'--k',lw=2)    
    #
    pl.savefig('num_vs_exact.png')

    pl.clf()
    Fro = hu/(h*np.sqrt(g*h));    
    pl.plot(x,Fro[:,my/2],'-b',lw=3)
    pl.savefig('Froude_number.png')
        
    dx=x[1]-x[0]; dy=y[1]-y[0]
    L1_error = dx*dy*np.sum(np.abs(h - h_exact))
    print("**********... L1-error: ", L1_error)
#
   
def plot_Ri(frame,
            file_prefix='claw',
            path='./_output/'):
    sol=Solution(frame,file_format='petsc',read_aux=False,path=path,file_prefix=file_prefix)
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers
    mx=len(x); my=len(y)
    Ri = get_ent_residual(sol.state)

    pl.xticks(fontsize=20)
    pl.yticks(fontsize=20)
    pl.plot(x,Ri[:,my/2],'-k',lw=3)
    

#
   
if __name__== "__main__":
    if not os.path.exists('./_plots_chj_1D'): os.mkdir('_plots_chj_1D')
    from_frame = 0
    to_frame   = 100
    frames=range(from_frame,to_frame+1)

    with_jump = False
    time = 100 if with_jump else 50

    time=50
    get_L1_error(time,with_jump=with_jump,plot_exact=True)
    
    if False:
        pl.figure(figsize=(8,4))
        get_L1_error(5,with_jump=with_jump,path='./_output_roe_N400',plot_exact=False,color='-c')
        get_L1_error(5,with_jump=with_jump,path='./_output_llf_N400',plot_exact=False,color='-r')
        get_L1_error(5,with_jump=with_jump,path='./_output_blended_N400',plot_exact=False,color='--b')
        pl.title('t=0.5',fontsize=20)
        pl.yticks([0.0,0.1,0.2,0.3])
        pl.ylim([0.0,0.3])
        pl.savefig('./_plots_chj_1D/outflow_1D_t0p5.png')
        
        pl.figure(figsize=(8,4))
        get_L1_error(13,with_jump=with_jump,path='./_output_roe_N400',plot_exact=False,color='-c')
        get_L1_error(13,with_jump=with_jump,path='./_output_llf_N400',plot_exact=False,color='-r')
        get_L1_error(13,with_jump=with_jump,path='./_output_blended_N400',plot_exact=False,color='--b')
        pl.title('t=1.3',fontsize=20)
        pl.yticks([0.0,0.1,0.2,0.3])
        pl.ylim([0.0,0.3])
        pl.savefig('./_plots_chj_1D/outflow_1D_t1p3.png')

        pl.figure(figsize=(8,4))
        get_L1_error(50,with_jump=with_jump,path='./_output_roe_N400',plot_exact=False,color='-c')
        get_L1_error(50,with_jump=with_jump,path='./_output_llf_N400',plot_exact=False,color='-r')
        get_L1_error(50,with_jump=with_jump,path='./_output_blended_N400',plot_exact=False,color='--b')
        pl.title('t=5',fontsize=20)
        pl.yticks([0.0,0.1,0.2,0.3])
        pl.ylim([0.0,0.3])
        pl.savefig('./_plots_chj_1D/outflow_1D_t10p0.png')

        # plot entropy residual
        pl.figure(figsize=(8,4))
        plot_Ri(5,path='./_output_blended_N400')
        pl.yticks([0.0,0.2,0.4,0.6,0.8,1.0])
        pl.savefig('./_plots_chj_1D/outflow_1D_Ri_t0p5.png')
        
        pl.figure(figsize=(8,4))    
        plot_Ri(13,path='./_output_blended_N400')
        pl.yticks([0.0,0.2,0.4,0.6,0.8,1.0])
        pl.savefig('./_plots_chj_1D/outflow_1D_Ri_t1p3.png')
        
        pl.figure(figsize=(8,4))
        plot_Ri(50,path='./_output_blended_N400')
        pl.yticks([0.0,0.2,0.4,0.6,0.8,1.0])
        pl.savefig('./_plots_chj_1D/outflow_1D_Ri_t10p0.png')
    #

    if False:
        pl.figure(figsize=(8,4))
        get_L1_error(5,with_jump=with_jump,path='./_output_roe_N400',plot_exact=False,color='-c')
        get_L1_error(5,with_jump=with_jump,path='./_output_llf_N400',plot_exact=False,color='-r')
        get_L1_error(5,with_jump=with_jump,path='./_output_blended_N400',plot_exact=False,color='-b')
        pl.title('t=0.5',fontsize=20)
        pl.yticks([0.1,0.2,0.3,0.4])
        pl.ylim([0.05,0.4])
        pl.savefig('./_plots_chj_1D/chj_1D_t0p5.png')
        
        pl.figure(figsize=(8,4))
        get_L1_error(13,with_jump=with_jump,path='./_output_roe_N400',plot_exact=False,color='-c')
        get_L1_error(13,with_jump=with_jump,path='./_output_llf_N400',plot_exact=False,color='-r')
        get_L1_error(13,with_jump=with_jump,path='./_output_blended_N400',plot_exact=False,color='-b')
        pl.title('t=1.3',fontsize=20)
        pl.yticks([0.1,0.2,0.3,0.4])
        pl.ylim([0.05,0.4])
        pl.savefig('./_plots_chj_1D/chj_1D_t1p3.png')

        pl.figure(figsize=(8,4))
        get_L1_error(100,with_jump=with_jump,path='./_output_roe_N400',plot_exact=False,color='-c')
        get_L1_error(100,with_jump=with_jump,path='./_output_llf_N400',plot_exact=False,color='-r')
        get_L1_error(100,with_jump=with_jump,path='./_output_blended_N400',plot_exact=False,color='-b')
        pl.title('t=10',fontsize=20)
        pl.yticks([0.1,0.2,0.3,0.4])
        pl.ylim([0.05,0.4])
        pl.savefig('./_plots_chj_1D/chj_1D_t10p0.png')

        # plot entropy residual
        pl.figure(figsize=(8,4))
        plot_Ri(5,path='./_output_blended_N400')
        pl.yticks([0.0,0.1,0.2,0.3,0.4])
        pl.savefig('./_plots_chj_1D/chj_1D_Ri_t0p5.png')
        
        pl.figure(figsize=(8,4))    
        plot_Ri(13,path='./_output_blended_N400')
        pl.yticks([0.0,0.1,0.2,0.3,0.4])
        pl.savefig('./_plots_chj_1D/chj_1D_Ri_t1p3.png')
        
        pl.figure(figsize=(8,4))
        plot_Ri(100,path='./_output_blended_N400')
        pl.yticks([0.0,0.1,0.2,0.3,0.4])
        pl.savefig('./_plots_chj_1D/chj_1D_Ri_t10p0.png')
        #

    pl.clf()
    if True:
        ##########################################
        if not os.path.exists('./_plots'): os.mkdir('./_plots')
        print('**********************')
        print('**********************')
        print('Plotting solution ...')
        for i in frames:
            plot_q(frame=i,
                   with_jump=with_jump)
            print ('frame '+str(i)+' plotted')
        #
    #


