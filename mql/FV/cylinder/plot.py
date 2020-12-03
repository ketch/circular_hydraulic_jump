from clawpack.petclaw.solution import Solution as petclaw_Solution
from clawpack.pyclaw.solution import Solution as pyclaw_Solution
#from petclaw.io.petsc import read_petsc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
from matplotlib import rc
#rc('text', usetex=True)
import numpy as np
import os

def plot_q(frame,
           plot_schlieren=True,
           use_petsc=True,
           bathymetry=False,
           file_prefix='claw',
           path='./_output/',
           plot_ent_residual=False,
           plot_pcolor=True,
           momMagLim=None,
           plot_slices=True):
    import sys
    sys.path.append('.')

    if use_petsc:
        sol=petclaw_Solution(frame,file_format='petsc',read_aux=False,path=path,file_prefix=file_prefix)
    else:
        file_prefix='fort'
        sol=pyclaw_Solution(frame,read_aux=False,path=path,file_prefix=file_prefix)
    #    
    x=sol.state.grid.x.centers; y=sol.state.grid.y.centers
    mx=len(x); my=len(y)

    h=sol.state.q[0,:,:]
    q1=sol.state.q[1,:,:]
    q2=sol.state.q[2,:,:]
    qMag = np.sqrt(q1**2 + q2**2)

    yy,xx = np.meshgrid(y,x)

    if frame < 10:
        str_frame = "000"+str(frame)
    elif frame < 100:
        str_frame = "00"+str(frame)
    elif frame < 1000:
        str_frame = "0"+str(frame)
    else:
        str_frame = str(frame)
    #

    grad_h = np.gradient(h)
    abs_grad_h = np.sqrt(grad_h[0]**2 + grad_h[1]**2)
    
    if plot_pcolor:
        fig = pl.figure(figsize=(4,8))
        if plot_schlieren:
            import clawpack.visclaw.colormaps as cm
            scm = cm.schlieren_colormap()
            pl.pcolormesh(xx,yy,abs_grad_h,cmap=scm)                
        #
        else:
            pl.pcolormesh(xx,yy,qMag)
            if momMagLim is not None:
                pl.clim(momMagLim[0],momMagLim[1])
            #
        #
        pl.axis('equal')
        pl.xticks(fontsize=20)
        pl.yticks(fontsize=20)
        pl.tight_layout()
        pl.xlim([0,40])
        pl.ylim([0,100])
        pl.title("t= "+str(sol.state.t),fontsize=20)        
        pl.savefig('./_plots/'+str_frame+'.png')
        pl.close()
    if plot_ent_residual:
        raise NotImplementedError
    if plot_slices:
        raise NotImplementedError
    #
#

if __name__== "__main__":
    if not os.path.exists('./_plots'): os.mkdir('_plots')
    from_frame = 0
    to_frame   = 1000
    frames=xrange(from_frame,to_frame+1)

    use_petsc = True
    plot_schlieren = True

    ##########################################
    # Get min and max values for color plots #
    if use_petsc:
        sol=petclaw_Solution(0,file_format='petsc',read_aux=False,path='./_output/',file_prefix='claw')
    else:
        sol=pyclaw_Solution(0,read_aux=False,path='./_output/',file_prefix='fort')
    #
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
               plot_schlieren=plot_schlieren,
               use_petsc=use_petsc,
               plot_pcolor=True,
               plot_ent_residual=False,
               #momMagLim=[momMagMin,momMagMax],
               plot_slices=False)
        print ('frame '+str(i)+' plotted')

