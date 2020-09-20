from clawpack.petclaw.solution import Solution
#from petclaw.io.petsc import read_petsc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
from matplotlib import rc
#rc('text', usetex=True)
import numpy as np
import os

def plot_q(frame,
           bathymetry=False,
           file_prefix='claw',
           path='./_output/',
           plot_pcolor=True,
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
               #hLim=[hMin,hMax],
               #momMagLim=[momMagMin,momMagMax],
               plot_slices=False)
        print ('frame '+str(i)+' plotted')

