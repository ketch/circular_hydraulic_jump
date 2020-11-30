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
           use_petsc=True,
           bathymetry=False,
           file_prefix='claw',
           path='./_output/',
           plot_ent_residual=False,
           plot_pcolor=True,
           hLim=None,
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
        # ***** Water depth ***** #
        #pl.subplot(121)
        #pl.pcolormesh(xx,yy,h)
        #pl.axis('equal')
        #cb = pl.colorbar();
        #if hLim is not None:
        #    pl.clim(hLim[0],hLim[1])
        #
        # ***** Magnitude of momentum ***** #
        #pl.subplot(122)
        import clawpack.visclaw.colormaps as cm
        scm = cm.schlieren_colormap()
        pl.pcolormesh(xx,yy,abs_grad_h,cmap=scm)
        #pl.pcolormesh(xx,yy,qMag)
        pl.axis('equal')
        #cb = pl.colorbar();
        #pl.yticks([])
        #if momMagLim is not None:
        #    pl.clim(momMagLim[0],momMagLim[1])
        #
        pl.xticks(fontsize=20)
        pl.yticks(fontsize=20)
        pl.tight_layout()
        pl.xlim([0,40])
        pl.ylim([0,100])
        pl.title("t= "+str(sol.state.t),fontsize=20)

        
        pl.savefig('./_plots/pcolor/'+str_frame+'.png')
        pl.close()
    if plot_ent_residual:
        raise NotImplementedError
    if plot_slices:
        pl.figure(figsize=(8,3))
        pl.plot(x,h[:,my/2],'-r',lw=1)
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

    use_petsc = False
    
    ##########################################
    # Get min and max values for color plots #
    #sol=Solution(0,file_format='petsc',read_aux=False,path='./_output/',file_prefix='claw')
    #hMin = sol.state.q[0,:,:].min(); hMax = sol.state.q[0,:,:].max()
    #print (hMin, hMax)
    #q1=sol.state.q[1,:,:]; q2=sol.state.q[2,:,:]
    #momMag = np.sqrt(q1**2 + q2**2)
    #momMagMin = momMag.min(); momMagMax = momMag.max()
    ###########################################
    
    if not os.path.exists('./_plots'): os.mkdir('./_plots')
    print('**********************')
    print('**********************')
    print('Plotting solution ...')
    for i in frames:
        plot_q(frame=i,
               use_petsc=use_petsc,
               plot_pcolor=True,
               plot_ent_residual=False,
               #hLim=[hMin,hMax],
               #momMagLim=[momMagMin,momMagMax],
               plot_slices=False)
        print ('frame '+str(i)+' plotted')

