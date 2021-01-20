from clawpack.petclaw.solution import Solution
#from petclaw.io.petsc import read_petsc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
from matplotlib import rc
#rc('text', usetex=True)
import numpy as np
import os

tol=1E-7

def plot_q(frame,
           plot_jump=False,
           plot_boundaries=False,
           plot_schlieren=True,
           plot_mesh=False,
           clim=None,
           xlim=None,
           ylim=None,
           file_prefix='claw',
           path='./_output/',
           plot_pcolor=True,
           plot_ent_residual=False,
           plot_slices=True):
    import sys
    sys.path.append('.')

    import mapc2p
    sol=Solution(frame,file_format='petsc',read_aux=False,path=path,file_prefix=file_prefix)
    
    xc_centers, yc_centers = sol.state.grid.c_centers    
    xc = np.zeros((xc_centers.shape[0],xc_centers.shape[1]+1))
    yc = np.zeros((yc_centers.shape[0],yc_centers.shape[1]+1))
    xc[:,:-1] = xc_centers
    xc[:,-1]  = xc_centers[:,0]
    yc[:,:-1] = yc_centers
    yc[:,-1]  = yc_centers[:,0]    
    
    xp,yp = mapc2p.mapc2p(xc,yc)
    
    # plot mesh
    if plot_mesh:
        pl.figure(figsize=(20,20))
        xc_edges, yc_edges = sol.state.grid.c_nodes
        xp_edges, yp_edges = mapc2p.mapc2p(xc_edges,yc_edges)
        pl.plot(xp_edges,yp_edges,'-k',lw=1.25)
        pl.plot(xp_edges.transpose(),yp_edges.transpose(),'-k',lw=1.25)
        pl.axis('equal')
        pl.xticks(fontsize=20)
        pl.yticks(fontsize=20)
        pl.xlim([-1,1])
        pl.ylim([-1,1])
        pl.savefig('mesh.png')
    #
    
    h=sol.state.q[0,:,:]
    q1=sol.state.q[1,:,:]
    q2=sol.state.q[2,:,:]
    qMag = np.sqrt(q1**2 + q2**2)

    grad_h = np.gradient(h)
    abs_grad_h = np.sqrt(grad_h[0]**2 + grad_h[1]**2)
    
    if frame < 10:
        str_frame = "000"+str(frame)
    elif frame < 100:
        str_frame = "00"+str(frame)
    elif frame < 1000:
        str_frame = "0"+str(frame)
    else:
        str_frame = str(frame)
    #

    print ('min and max of h: ', h.min(), h.max())
    # get the location of the jump
    mx = h.shape[0]
    my = h.shape[1]
    radius = np.zeros(my)
    rp = np.sqrt(xp**2 + yp**2)
    for j in range(my):
        for i in range(mx):
            if h[i,j] > 1.1*h[0,j]:
                radius[j] = rp[i,j]
                break
            #
        #
    #
    
    if plot_pcolor:
        fig = pl.figure(figsize=(8,8))
        
        if plot_schlieren:
            import clawpack.visclaw.colormaps as cm
            scm = cm.schlieren_colormap()
            pl.pcolormesh(xp,yp,abs_grad_h,cmap=scm)
        else:
            #pl.pcolormesh(xp,yp,qMag,cmap='Blues')
            #pl.pcolormesh(xp,yp,qMag)
            pl.pcolormesh(xp,yp,h)
            if clim is not None:
                pl.clim([clim[0],clim[1]])
            #
            pl.colorbar()
            if plot_jump:
                angles = np.linspace(0,2*np.pi,len(radius))
                xr = radius*np.cos(angles)
                yr = radius*np.sin(angles)
                pl.plot(xr,yr,'-k',lw=1)
            #
        #

        # plot boundaries
        if plot_boundaries:
            xo=np.linspace(-1,1,1000)
            yo=np.sqrt(1.0-xo**2)
            pl.plot(xo,yo,'-k',lw=1)
            pl.plot(xo,-yo,'-k',lw=1)
            
            xi=np.linspace(-0.1,0.1,1000)
            yi=np.sqrt(0.1**2-xi**2)
            pl.plot(xi,yi,'-k',lw=1)
            pl.plot(xi,-yi,'-k',lw=1)
        #
        pl.axis('equal')
        pl.xticks([])
        pl.yticks([])
        pl.axis('off')

        if xlim is not None:
            pl.xlim(xlim)
        #
        if ylim is not None:
            pl.ylim(ylim)
        #
        #pl.suptitle("t= "+str(sol.state.t),fontsize=20)
        pl.savefig('./_plots/'+str_frame+'.png')
        pl.close()
    #
    if plot_slices:
        raise NotImplementedError
    if plot_ent_residual:
        raise NotImplementedError
    #
#

if __name__== "__main__":
    if not os.path.exists('./_plots'): os.mkdir('_plots')
    from_frame = 0
    to_frame   = 1000
    frames=xrange(from_frame,to_frame+1)

    plot_schlieren = False
    xlim=None #[0,0.4]
    ylim=None #[0,0.4]
    
    if not os.path.exists('./_plots'): os.mkdir('./_plots')
    print('**********************')
    print('**********************')
    print('Plotting solution ...')
    if True:
        sol=Solution(0,file_format='petsc',read_aux=False,path='./_output/',file_prefix='claw')
        q1=sol.state.q[1,:,:]
        q2=sol.state.q[2,:,:]
        qMag = np.sqrt(q1**2 + q2**2)
        clim=[qMag.min(),qMag.max()]
        #for i in [50]:
        for i in frames:
            plot_q(frame=i,
                   #clim=clim,
                   xlim=xlim,
                   ylim=ylim,
                   plot_schlieren=plot_schlieren,
                   plot_pcolor=True,
                   plot_mesh=False,
                   plot_ent_residual=False,
                   plot_slices=False)
            print ('frame '+str(i)+' plotted')
        #
    #


