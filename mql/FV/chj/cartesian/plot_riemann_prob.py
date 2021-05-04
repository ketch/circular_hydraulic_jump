from clawpack.petclaw.solution import Solution
#from petclaw.io.petsc import read_petsc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
from matplotlib import rc
#rc('text', usetex=True)
import numpy as np
import os
from ent_residual import * 

bed_dry=False
hl = 0.005
hr = 0.001
g  = 1.0
x0 = 5.0
from scipy.optimize import fsolve

def fun(x):
    return -8*g*hr*x**2*(np.sqrt(g*hl)-x)**2+(x**2-g*hr)**2*(x**2+g*hr)
#
def exact_solution(x,t):
    if bed_dry:
        xA=x0-t*np.sqrt(g*hl)
        xB=x0+2*t*np.sqrt(g*hl)
        h = (x<=xA) * hl + (x>xA)*(x<xB) * 4.0/(9*g)*(np.sqrt(g*hl)-(x-x0)/(2*t))**2
        u = (x>xA)*(x<xB) * 2.0/3*((x-x0)/t+np.sqrt(g*hl))
    #
    else:
        cg=np.sqrt(g*hl)
        cd=np.sqrt(g*hr)    
        cm = fsolve(fun,0.5*(cd+cg))[0]
        hm = cm**2/g

        xA=x0-t*np.sqrt(g*hl)
        xB=x0+t*(2*np.sqrt(g*hl)-3*cm)
        xC=x0+t*(2*cm**2*(np.sqrt(g*hl)-cm))/(cm**2-g*hr)

        h = (x<=xA) * hl + (x>xA)*(x<xB) * 4.0/(9*g)*(np.sqrt(g*hl)-(x-x0)/(2*t))**2 + (x>=xB)*(x<=xC) * hm + (x>xC)*hr
        u = (x>xA)*(x<xB) * 2.0/3*((x-x0)/t+np.sqrt(g*hl)) + (x>=xB)*(x<=xC) * 2 * (np.sqrt(g*hl)-cm)
    #
    return h,h*u
#

def plot_q(frame,
           file_prefix='claw',
           path='./_output/'):
    import sys
    sys.path.append('.')

    sol=Solution(frame,file_format='petsc',read_aux=False,path=path,file_prefix=file_prefix)
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers
    mx=len(x); my=len(y)

    h=sol.state.q[0,:,:]
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

    # compute error
    h_exact, hu_exact = exact_solution(x,sol.state.t)

    if True:
        fig, axs = pl.subplots(2,1)
        axs[0].plot(x,h[:,my/2],'-r',lw=2)
        axs[0].plot(x,h_exact,'--k',lw=2)
        axs[0].set_ylim([0,0.0055])    
        axs[1].plot(x,Ri[:,my/2],'-r',lw=1)    
        pl.suptitle("t= "+str(sol.state.t),fontsize=20)
    else:
        pl.figure(figsize=(10,4))        
        pl.plot(x,h[:,my/2],'-r',lw=3)
        pl.plot(x,h_exact,'--k',lw=2)
        pl.xticks(fontsize=20)
        pl.yticks(fontsize=20)
        pl.ylim([-0.0005,0.0055])
    #
    pl.savefig('./_plots_riemann_prob/h_'+str_frame+'_slices.png')
    pl.close()

def get_L1_error(frame,
                 color='-r',
                 file_prefix='claw',
                 path='./_output/'):
    sol=Solution(frame,file_format='petsc',read_aux=False,path=path,file_prefix=file_prefix)
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers
    mx=len(x); my=len(y)

    h=sol.state.q[0,:,:]

    # compute error
    h_exact, hu_exact = exact_solution(x,sol.state.t)
    dx = x[1]-x[0]
    error = dx*np.sum(np.abs(h_exact - h[:,my/2]))
    print "********** L1 error: ", error
    
    pl.plot(x,h[:,my/2],color,lw=3)
    pl.plot(x,h_exact,'--k',lw=2)
    pl.xticks(fontsize=20)
    pl.yticks(fontsize=20)
    pl.ylim([-0.0005,0.0055])

#

if __name__== "__main__":
    if not os.path.exists('./_plots_riemann_prob'): os.mkdir('_plots_riemann_prob')
    from_frame = 0
    to_frame   = 10
    frames=xrange(from_frame,to_frame+1)

    get_L1_error(10)
    #pl.savefig('./_plots_riemann_prob/plot.png')
    
    if False:
        pl.figure(figsize=(10,4))        
        get_L1_error(10,path='./_output_roe_N400/',color='-c')
        get_L1_error(10,path='./_output_roe_N800/',color='-r')
        get_L1_error(10,path='./_output_roe_N1600/',color='-b')
        pl.savefig('./_plots_riemann_prob/rp_dry_roe.png')
        pl.close()

        pl.figure(figsize=(10,4))        
        get_L1_error(10,path='./_output_llf_N400/',color='-c')
        get_L1_error(10,path='./_output_llf_N800/',color='-r')
        get_L1_error(10,path='./_output_llf_N1600/',color='-b')
        pl.savefig('./_plots_riemann_prob/rp_dry_llf.png')
        pl.close()

        pl.figure(figsize=(10,4))        
        get_L1_error(10,path='./_output_blended_N400/',color='-c')
        get_L1_error(10,path='./_output_blended_N800/',color='-r')
        get_L1_error(10,path='./_output_blended_N1600/',color='-b')
        pl.savefig('./_plots_riemann_prob/rp_dry_blended.png')
        pl.close()

        pl.figure(figsize=(10,4))        
        get_L1_error(10,path='./_output_roe_with_lmin_N400/',color='-c')
        get_L1_error(10,path='./_output_roe_with_lmin_N800/',color='-r')
        get_L1_error(10,path='./_output_roe_with_lmin_N1600/',color='-b')
        pl.savefig('./_plots_riemann_prob/rp_dry_roe_with_lmin.png')
        pl.close()
    #

    if False:
        pl.figure(figsize=(8,4))
        get_L1_error(10,path='./_output_roe_N400/',color='-c')
        get_L1_error(10,path='./_output_roe_N800/',color='-r')
        get_L1_error(10,path='./_output_roe_N1600/',color='-b')
        pl.xlim([3,7])
        pl.savefig('./_plots_riemann_prob/rp_wet_roe.png')
        pl.close()

        pl.figure(figsize=(8,4))        
        get_L1_error(10,path='./_output_llf_N400/',color='-c')
        get_L1_error(10,path='./_output_llf_N800/',color='-r')
        get_L1_error(10,path='./_output_llf_N1600/',color='-b')
        pl.xlim([3,7])
        pl.savefig('./_plots_riemann_prob/rp_wet_llf.png')
        pl.close()

        pl.figure(figsize=(8,4))        
        get_L1_error(10,path='./_output_blended_N400/',color='-c')
        get_L1_error(10,path='./_output_blended_N800/',color='-r')
        get_L1_error(10,path='./_output_blended_N1600/',color='-b')
        pl.xlim([3,7])
        pl.savefig('./_plots_riemann_prob/rp_wet_blended.png')
        pl.close()
    #
    
    if True:
        ##########################################
        if not os.path.exists('./_plots'): os.mkdir('./_plots')
        print('**********************')
        print('**********************')
        print('Plotting solution ...')
        for i in frames:
            plot_q(frame=i)
            print ('frame '+str(i)+' plotted')
        #
    #
