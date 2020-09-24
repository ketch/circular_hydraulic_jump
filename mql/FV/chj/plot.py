from clawpack.petclaw.solution import Solution
#from petclaw.io.petsc import read_petsc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
from matplotlib import rc
#rc('text', usetex=True)
import numpy as np
import os

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
    
def get_entropy(g,h,hu,hv):
    return 0.5*g*h**2 + 0.5*(hu**2+hv**2)/(h+1E-7)
#

def xFlux(g,h,hu,hv):
    tol=1E-7
    return (0.5*(hu**2+hv**2)/(h+tol) + g*h**2)*hu/(h+tol) 
#
def yFlux(g,h,hu,hv):
    tol=1E-7
    return (0.5*(hu**2+hv**2)/(h+tol) + g*h**2)*hv/(h+tol) 
#

def plot_q(frame,
           bathymetry=False,
           file_prefix='claw',
           path='./_output/',
           plot_pcolor=True,
           plot_ent_residual=False,
           hLim=None,
           momMagLim=None,
           plot_slices=True):
    import sys
    sys.path.append('.')

    sol=Solution(frame,file_format='petsc',read_aux=False,path=path,file_prefix=file_prefix)
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
        old_sol = Solution(frame,file_format='petsc',read_aux=False,path=path,file_prefix='old_soln')
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
        fig = pl.figure(figsize=(12,4))
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
        entropy_tnm1 = get_entropy(g,old_h,old_q1,old_q2)
        entropy_tn = get_entropy(g,h,q1,q2)

        dt = sol.t-old_sol.t
        dx = sol.state.delta[0]
        dy = sol.state.delta[1]
        tol = 1E-7
        
        # reconstruct solution at the interfaces
        x_hL, x_hR, y_hL, y_hR = reconstruct_solution(h)
        x_huL, x_huR, y_huL, y_huR = reconstruct_solution(q1)
        x_hvL, x_hvR, y_hvL, y_hvR = reconstruct_solution(q2)

        int_div_ent_flux = (xFlux(g,x_hR,x_huR,x_hvR) - xFlux(g,x_hL,x_huL,x_hvL))*dy + (yFlux(g,y_hR,y_huR,y_hvR) - yFlux(g,y_hL,y_huL,y_hvL))*dx
        
        ent_residual = (entropy_tn - entropy_tnm1)/dt*dx*dy + int_div_ent_flux
        xc, yc = sol.state.p_centers
        rc = np.sqrt(xc**2 + yc**2)
        in_solid = np.where((rc <= 0.1))
        ent_residual[in_solid[0],in_solid[1]] = 0
        
        sv=100
        #ent_residual = sv*(ent_residual >= sv) - sv*(ent_residual <= -sv) + ent_residual*(ent_residual < sv)*(ent_residual>-sv)

        pos_ent_residual = ent_residual * (ent_residual > 0.0)
        neg_ent_residual = ent_residual * (ent_residual <= 0.0)
        
        # plot
        fig = pl.figure(figsize=(12,4))
        pl.subplot(121)
        # entropy dissipation
        pl.pcolormesh(xx,yy,neg_ent_residual)
        cb = pl.colorbar();
        # entropy production
        pl.subplot(122)
        pl.pcolormesh(xx,yy,pos_ent_residual)
        cb = pl.colorbar();        

        pl.suptitle("t= "+str(sol.state.t),fontsize=20)        
        pl.savefig('./_plots/pcolor/ent_residual_'+str_frame+'.png')
        pl.close()
    if plot_slices:
        pl.figure(figsize=(8,3))
        pl.plot(x,eta[:,my/2],'-r',lw=1)
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
    frames=xrange(from_frame,to_frame+1)

    ##########################################
    # Get min and max values for color plots #
    sol=Solution(0,file_format='petsc',read_aux=False,path='./_output/',file_prefix='claw')
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
               plot_ent_residual=True,
               #hLim=[hMin,hMax],
               #momMagLim=[momMagMin,momMagMax],
               plot_slices=False)
        print ('frame '+str(i)+' plotted')

