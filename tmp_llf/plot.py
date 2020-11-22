from clawpack.pyclaw.solution import Solution
#from petclaw.io.petsc import read_petsc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
from matplotlib import rc
#rc('text', usetex=True)
import numpy as np
import os

# *********************************************************** #
# ********** FUNCTIONS TO COMPUTE ENTROPY RESIDUAL ********** #
# *********************************************************** #
tol=1E-7
def reconstruct_solution(h):
    x_hL = np.zeros_like(h)
    x_hR = np.zeros_like(h)
    x_hL[2:-2,:] = (-3.*h[:-4,:] + 27*h[1:-3,:] + 47*h[2:-2,:] - 13*h[3:-1,:] + 2*h[4:,:])/60.
    x_hR[2:-2,:] = (2.*h[:-4,:] - 13*h[1:-3,:] + 47*h[2:-2,:] + 27*h[3:-1,:] - 3*h[4:,:])/60.

    y_hL = np.zeros_like(h)
    y_hR = np.zeros_like(h)
    y_hL[:,2:-2] = (-3.*h[:,:-4] + 27*h[:,1:-3] + 47*h[:,2:-2] - 13*h[:,3:-1] + 2*h[:,4:])/60.
    y_hR[:,2:-2] = (2.*h[:,:-4] - 13*h[:,1:-3] + 47*h[:,2:-2] + 27*h[:,3:-1] - 3*h[:,4:])/60.

    return x_hL, x_hR, y_hL, y_hR
#
# ENTROPY AND ITS DERIVATIVES
def get_entropy(g,h,hu,hv):
    return 0.5*g*h**2 + 0.5*(hu**2+hv**2)/(h+tol)
def dEntropy_dh(g,h,hu,hv):
    return g*h - 0.5*(hu**2+hv**2)/(h**2+tol)
def dEntropy_dhu(g,h,hu,hv):
    return hu/(h+tol)
def dEntropy_dhv(g,h,hu,hv):
    return hv/(h+tol)
# FLUXES OF THE SHALLOW WATER EQUATIONS
def xFluxTerm1(g,h,hu,hv):
    return hu
def xFluxTerm2(g,h,hu,hv):
    return hu**2/(h+tol) + 0.5*g*h**2
def xFluxTerm3(g,h,hu,hv):
    return hu*hv/(h+tol)
def yFluxTerm1(g,h,hu,hv):
    return hv
def yFluxTerm2(g,h,hu,hv):
    return hu*hv/(h+tol)
def yFluxTerm3(g,h,hu,hv):
    return hv**2/(h+tol) + 0.5*g*h**2
# DERIVATIVES OF ENTROPY FLUX 
def xEntFlux(g,h,hu,hv):
    return (0.5*(hu**2+hv**2)/(h+tol) + g*h**2)*hu/(h+tol) 
def yEntFlux(g,h,hu,hv):
    return (0.5*(hu**2+hv**2)/(h+tol) + g*h**2)*hv/(h+tol)
# SUM INTEGRALS OVER PATCH
def sumIntOverPatch(I):
    I[1:-1,1:-1] = (I[1:-1,1:-1] +
                    I[:-2,1:-1] + I[2:,1:-1] +
                    I[1:-1,:-2] + I[1:-1,2:] +
                    I[2:,2:] + I[:-2,2:] + I[2:,:-2] +I[:-2,:-2])
    
# *********************************************************** #

def plot_q(frame,
           bathymetry=False,
           path='./_output/',
           plot_pcolor=True,
           plot_ent_residual=False,
           hLim=None,
           momMagLim=None,
           plot_slices=True):
    import sys
    sys.path.append('.')

    sol=Solution(frame,read_aux=False,path=path)
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers
    mx=len(x); my=len(y)

    h=sol.state.q[0,:,:]
    if bathymetry:
        b=sol.aux[0,:,:]        
        eta=h+b
    else:
        eta=h
    #
    q1=sol.state.q[1,:,:]
    q2=sol.state.q[2,:,:]
    qMag = np.sqrt(q1**2 + q2**2)

    # old solution 
    if plot_ent_residual:
        old_sol = Solution(frame,read_aux=False,path=path,file_prefix='old_soln')
        old_h = old_sol.state.q[0,:,:]
        if bathymetry:
            old_eta = old_h + b
        else:
            old_eta = old_h
        #
        old_q1 = old_sol.state.q[1,:,:]
        old_q2 = old_sol.state.q[2,:,:]
        old_qMag = np.sqrt(old_q1**2 + old_q2**2)
    #    
    yy,xx = np.meshgrid(y,x)

    if frame < 10:
        str_frame = "000"+str(frame)
    elif frame < 100:
        str_frame = "00"+str(frame)
    elif frame < 1000:
        str_frame = "0"+str(frame)
    else:
        str_frame = str(frame)

    if plot_pcolor:
        fig = pl.figure(figsize=(24,8))
        # ***** Water depth ***** #
        pl.subplot(121)
        pl.pcolormesh(xx,yy,eta)
        cb = pl.colorbar();
        if hLim is not None:
            pl.clim(hLim[0],hLim[1])
        #
        # ***** Magnitude of momentum ***** #
        pl.subplot(122)
        pl.pcolormesh(xx,yy,qMag)
        cb = pl.colorbar();
        if momMagLim is not None:
            pl.clim(momMagLim[0],momMagLim[1])
        #

        pl.suptitle("t= "+str(sol.state.t),fontsize=20)
        pl.savefig('./_plots/pcolor/'+str_frame+'.png')
        pl.close()
    if plot_ent_residual:
        g=sol.state.problem_data['grav']
        dx = sol.state.delta[0]
        dy = sol.state.delta[1]
        
        # reconstruct solution at the interfaces
        x_hL, x_hR, y_hL, y_hR = reconstruct_solution(h)
        x_huL, x_huR, y_huL, y_huR = reconstruct_solution(q1)
        x_hvL, x_hvR, y_hvL, y_hvR = reconstruct_solution(q2)

        int_div_ent_flux = ( (xEntFlux(g,x_hR,x_huR,x_hvR) - xEntFlux(g,x_hL,x_huL,x_hvL))*dy +
                             (yEntFlux(g,y_hR,y_huR,y_hvR) - yEntFlux(g,y_hL,y_huL,y_hvL))*dx )
        sumIntOverPatch(int_div_ent_flux)
        # derivative of the entropy 
        eta_prime_term1 = dEntropy_dh(g,h,q1,q2)
        eta_prime_term2 = dEntropy_dhu(g,h,q1,q2)
        eta_prime_term3 = dEntropy_dhv(g,h,q1,q2)
        
        # integral of the divergence of the flux 
        int_div_flux_term1 = ( (xFluxTerm1(g,x_hR,x_huR,x_hvR) - xFluxTerm1(g,x_hL,x_huL,x_hvL))*dy +
                               (yFluxTerm1(g,y_hR,y_huR,y_hvR) - yFluxTerm1(g,y_hL,y_huL,y_hvL))*dx )
        int_div_flux_term2 = ( (xFluxTerm2(g,x_hR,x_huR,x_hvR) - xFluxTerm2(g,x_hL,x_huL,x_hvL))*dy +
                               (yFluxTerm2(g,y_hR,y_huR,y_hvR) - yFluxTerm2(g,y_hL,y_huL,y_hvL))*dx )
        int_div_flux_term3 = ( (xFluxTerm3(g,x_hR,x_huR,x_hvR) - xFluxTerm3(g,x_hL,x_huL,x_hvL))*dy +
                               (yFluxTerm3(g,y_hR,y_huR,y_hvR) - yFluxTerm3(g,y_hL,y_huL,y_hvL))*dx )
        sumIntOverPatch(int_div_flux_term1)
        sumIntOverPatch(int_div_flux_term2)
        sumIntOverPatch(int_div_flux_term3)
        # get dot product of the derivative of entropy and the int of the div of the flux
        eta_prime_times_int_div_flux = (eta_prime_term1 * int_div_flux_term1 +
                                        eta_prime_term2 * int_div_flux_term2 +
                                        eta_prime_term3 * int_div_flux_term3)
        
        # get l2 norm of the derivative of the entropy
        #abs_eta_prime = np.sqrt(eta_prime_term1**2 + eta_prime_term2**2 + eta_prime_term3**2) 
        
        # get the l2 norm of the integral of the divergence of the flux
        #abs_int_div_flux = np.sqrt(int_div_flux_term1**2 + int_div_flux_term2**2 + int_div_flux_term3**2)
        
        # entropy residual 
        Ri = ( np.abs(int_div_ent_flux - eta_prime_times_int_div_flux) /
               #(np.abs(int_div_ent_flux) + abs_eta_prime * abs_int_div_flux + tol) )
               (np.abs(int_div_ent_flux) + np.abs(eta_prime_times_int_div_flux) + tol) )

        # entropy residual based on time derivative 
        dt = sol.t-old_sol.t
        entropy_tnm1 = get_entropy(g,old_h,old_q1,old_q2)
        entropy_tn = get_entropy(g,h,q1,q2)
        int_dt_ent = (entropy_tn - entropy_tnm1)/dt*dx*dy
        sumIntOverPatch(int_dt_ent)
        ent_residual = int_dt_ent + int_div_ent_flux
        Ri_v2 = (np.abs(ent_residual)/
                 (np.abs(int_dt_ent) + np.abs(int_div_ent_flux)+tol))
        #Ri_v2 = ( ent_residual*(ent_residual>0.) / 
        #          (np.abs(dt_entropy) + np.abs(int_div_ent_flux)+tol))        
        # plot
        fig = pl.figure(figsize=(24,8))
        pl.subplot(121)
        # entropy dissipation
        pl.pcolormesh(xx,yy, Ri)
        #pl.pcolormesh(xx,yy,neg_ent_residual)
        #pl.xlim([0,0.3])
        #pl.ylim([0,0.3])
        cb = pl.colorbar();
        # entropy production
        pl.subplot(122)
        pl.pcolormesh(xx,yy, Ri_v2)
        #pl.xlim([0,0.3])
        #pl.ylim([0,0.3])
        cb = pl.colorbar();        

        pl.suptitle("t= "+str(sol.state.t),fontsize=20)        
        pl.savefig('./_plots/pcolor/ent_residual_'+str_frame+'.png')
        pl.close()
    if plot_slices:
        pl.figure(figsize=(8,3))
        pl.plot(x,eta[:,my//2],'-r',lw=1)
        #pl.plot(x,qMag[:,my/2],'-r',lw=1)
        pl.title("t= "+str(sol.state.t),fontsize=20)
        pl.xlabel('x',fontsize=20)
        pl.ylabel('Surface',fontsize=20)
        pl.xticks(size=20); pl.yticks(size=20)
        pl.savefig('./_plots/slices/qMag_'+str_frame+'_slices.png')
        pl.close()

if __name__== "__main__":
    if not os.path.exists('./_plots'): os.mkdir('_plots')
    if not os.path.exists('./_plots/pcolor'): os.mkdir('_plots/pcolor')
    if not os.path.exists('./_plots/pcolor/momMag'): os.mkdir('_plots/pcolor/momMag')
    if not os.path.exists('./_plots/pcolor/eta'): os.mkdir('_plots/pcolor/eta')
    if not os.path.exists('./_plots/slices'): os.mkdir('_plots/slices')
    from_frame = 0
    to_frame   = 1000
    frames=range(from_frame,to_frame+1)

    ##########################################
    # Get min and max values for color plots #
    sol=Solution(0,read_aux=False,path='./_output/')
    hMin = sol.state.q[0,:,:].min(); hMax = sol.state.q[0,:,:].max()
    print (hMin, hMax)
    q1=sol.state.q[1,:,:]; q2=sol.state.q[2,:,:]
    momMag = np.sqrt(q1**2 + q2**2)
    momMagMin = momMag.min(); momMagMax = momMag.max()
    ###########################################
    
    if not os.path.exists('./_plots'): os.mkdir('./_plots')
    print('**********************')
    print('**********************')
    print('Plotting solution ...')
    for i in frames:
        plot_q(frame=i,
               plot_pcolor=True,
               plot_ent_residual=False,
               #hLim=[hMin,hMax],
               #momMagLim=[momMagMin,momMagMax],
               plot_slices=True)
        print ('frame '+str(i)+' plotted')

