"""
Code to simulate 2D circular hydraulic jump in PyClaw.

This version uses a mapped (annular) grid.

In order to get a jump, you must either introduce a bathymetric barrier (circular ring),
include bottom friction, or impose a subcritical flow at the outer boundary.
Bathymetry is not enabled for the mapped grid solver at present.
"""

from matplotlib import animation
import matplotlib.pyplot as plt
from clawpack import pyclaw
import numpy as np

# If plotting in a Jupyter notebook, do this:
#from clawpack.visclaw.JSAnimation import IPython_display

# Mapped grid parameters
r_lower = 0.2
r_upper = 1.0
theta_lower = 0.0
theta_upper = 2*np.pi

def compute_geometric_values(grid, aux):
    """
    Determine the normal vector and the length ratio for each edge.
    Store them in aux.
    Inputs:
        - xp, yp: physical coordinates of nodes (corners)
    """
    xcorners, ycorners = grid.p_edges
    v = np.empty(2)
    w = np.empty(2)
    xpcorn = np.empty(5)
    ypcorn = np.empty(5)

    dx, dy = grid.delta

    for i in range(xcorners.shape[0]):
        for j in range(xcorners.shape[1]):

            xpcorn[0] = xcorners[i,j]
            xpcorn[1] = xcorners[i+1,j]
            xpcorn[2] = xcorners[i+1,j+1]
            xpcorn[3] = xcorners[i,j+1]
            xpcorn[4] = xcorners[i,j]

            ypcorn[0] = ycorners[i,j]
            ypcorn[1] = ycorners[i+1,j]
            ypcorn[2] = ycorners[i+1,j+1]
            ypcorn[3] = ycorners[i,j+1]
            ypcorn[4] = ycorners[i,j]

            area = 0.5*sum((ypcorn[:-1]+ypcorn[1:])*(xpcorn[1:]-xpcorn[:-1]))

            # Compute bottom/left edge vectors
            v[0] = xcorners[i+1,j] - xcorners[i,j]
            v[1] = ycorners[i+1,j] - ycorners[i,j]
            w[0] = xcorners[i,j+1] - xcorners[i,j]
            w[1] = ycorners[i,j+1] - ycorners[i,j]

            xlength = np.linalg.norm(v)
            ylength = np.linalg.norm(y)

            # Normal to x edge
            aux[0,i,j] =  v[1]/xlength
            aux[1,i,j] = -v[2]/xlength

            aux[2,i,j] = xlength/dx # x edge length ratio

            # Normal to y edge
            aux[3,i,j] = -w[1]/ylength
            aux[4,i,j] =  w[2]/ylength

            aux[5,i,j] = ylength/dy # y edge length ratio

            aux[6,i,j] = area

def mapc2p_annulus(xc, yc):
    """
    Specifies the mapping to curvilinear coordinates.

    Inputs: c_centers = Computational cell centers
                 [array ([Xc1, Xc2, ...]), array([Yc1, Yc2, ...])]

    Output: p_centers = Physical cell centers
                 [array ([Xp1, Xp2, ...]), array([Yp1, Yp2, ...])]
    """  
    p_centers = []

    # Polar coordinates (first coordinate = radius,  second coordinate = theta)
    p_centers.append(xc[:]*np.cos(yc[:]))
    p_centers.append(xc[:]*np.sin(yc[:]))
    
    return p_centers

def plot_surface(claw, make_anim=True, save_plots=False, frames=101, val='surface',
                 vmin=0., vmax=0.5, clim=None, bathymetry=False, plotdir='./_plots'):
    """
    Plot results of 2D shallow water simulation as a pcolor plot and animate
    (intended for Jupyter notebook).

    If variable bathymetry is used, set bathymetry=True.
    """
    print('plotting')
    fig = plt.figure(figsize=[10,10])
    ax1 = fig.add_subplot(111)
    try:
        # Use solution in memory
        frame = claw.frames[0]
    except:
        # Read solution from disk
        frame = pyclaw.Solution(0,path=claw,read_aux=bathymetry)
    if bathymetry:
        b = frame.aux[0,:,:]
    else:
        b = 0.
    h = frame.q[0,:,:]
    if val == 'surface':
        qq = np.maximum(b,h+b)
    elif val=='u':
        hu = frame.q[1,:,:]
        qq = hu #/h
        

    x, y = frame.state.grid.p_centers    

    im = ax1.imshow(qq.T, #cmap='Blues',
                    extent=[x.min(), x.max(), y.min(), y.max()],
                    vmin=vmin, vmax=vmax,
                    interpolation='nearest', origin='lower')
                    
    plt.colorbar(im)
    def fplot(frame_number):
        try:
            frame = claw.frames[frame_number]
        except:
            frame = pyclaw.Solution(frame_number,path=claw,read_aux=bathymetry)
        h = frame.q[0,:,:]
        if val == 'surface':
            qq = np.maximum(b,h+b)
        elif val=='u':
            hu = frame.q[1,:,:]
            qq = hu #/h
        im.set_data(qq.T)
        if clim:
            im.set_clim(*clim)
        if save_plots:
            fname = plotdir+'/frame'+str(frame_number).zfill(4)+'.eps'
            fig.savefig(fname)   
        return im,

    if make_anim:
        return animation.FuncAnimation(fig, fplot, frames=frames, interval=40, repeat=True)
    elif save_plots:
        import os
        if not os.path.exists(plotdir):
            os.makedirs(plotdir)
        for i in range(frames):
            fplot(i)


def jet(state, dim, _, qbc, __, num_ghost):
    "Jet inflow BC at inner boundary."
    rc, thetac = state.grid.p_centers_with_ghost(num_ghost)
    r0 = r_lower
    h0 =  state.problem_data['h0']
    u0 =  state.problem_data['u0']
    
    if dim.name == 'r':
        qbc[0,:num_ghost,:] = h0[:num_ghost,:]
        qbc[1,:num_ghost,:] = h0*u0
        qbc[2,:num_ghost,:] = 0.

def subsonic_boundary_upper(state, dim, _, qbc, __, num_ghost):
    "Subsonic outflow BC at fixed Froude number."
    xc, yc = state.grid.p_centers_with_ghost(num_ghost)
    rc = np.sqrt(xc**2 + yc**2)
    r0 = r_lower
    h0 =  state.problem_data['h0']
    u0 =  state.problem_data['u0']
    beta = r0*h0*u0
    F = state.problem_data['F_bdy'] # Froude number at boundary
    g = state.problem_data['grav']
    
    h = (beta/(rc*F*np.sqrt(g)))**(2./3)
    unorm = beta / (rc*h)
    
    if dim.name == 'r':
        qbc[0,-num_ghost:,:] = h[-num_ghost:,:]
        qbc[1,-num_ghost:,:] = h[-num_ghost:,:]*unorm[-num_ghost:,:]
        qbc[2,-num_ghost:,:] = 0.
    else:
        raise Exception(dim)
        
def step_friction(solver, state, dt):
    "Friction source term:  -cf u / h.  This version is for Classic."
    cf = state.problem_data['cf']
    q = state.q
    h = q[0,:,:]
    u = q[1,:,:]/h
    v = q[2,:,:]/h

    q[1,:,:] = q[1,:,:] - dt*cf*u/h
    q[2,:,:] = q[2,:,:] - dt*cf*v/h
        
def setup(h0=0.5, u0=0.75, h_inf=0.15, g=1., num_cells_r=100,
          num_cells_theta=100, tfinal=1,
          solver_type='classic', num_output_times=10,
          boundary='subcritical', outdir='./_output', friction=False,
          friction_coeff=0.01, F_bdy=0.1, use_petsc=False, 
          kalpha=1./3, kbeta=1.3, kepsilon=1.e-3):
    
    from clawpack import riemann
    if use_petsc:
        from clawpack import petclaw as pyclaw
    else:
        from clawpack import pyclaw
    import shallow_quad_hllemcc_2D

    if solver_type == 'classic':
        solver = pyclaw.ClawSolver2D(riemann.shallow_quad_hllemcc_2D)
        solver.cfl_max     = 0.9
        solver.cfl_desired = 0.8
        solver.transverse_waves = 2
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver2D(riemann.shallow_quad_hllemcc_2D)

    solver.num_eqn = 3
    solver.num_waves = 3
    solver.fwave = False

    # Periodic BCs in theta direction
    solver.bc_lower[1] = pyclaw.BC.periodic
    solver.bc_upper[1] = pyclaw.BC.periodic
    solver.aux_bc_lower[1] = pyclaw.BC.periodic
    solver.aux_bc_upper[1] = pyclaw.BC.periodic

    # Jet at inner boundary
    solver.bc_lower[0] = pyclaw.BC.custom
    solver.user_bc_lower = jet

    if boundary == 'outflow':
        if friction == False:
            raise Exception('Either friction or subcritical outflow is required to produce a jump.')
        solver.bc_upper[0] = pyclaw.BC.extrap
    elif boundary == 'subcritical':  # subcritical boundary condition
        solver.bc_upper[0] = pyclaw.BC.custom
        solver.user_bc_upper = subsonic_boundary_upper
    elif boundary == 'wall':
        solver.bc_upper[0] = pyclaw.BC.wall

    if friction:
        solver.step_source = step_friction
        solver.source_split = 1

    r = pyclaw.Dimension(r_lower,r_upper,num_cells_r,name='r')
    y = pyclaw.Dimension(ylower,yupper,num_cells_theta,name='theta')
    domain = pyclaw.Domain([r,theta])

    state = pyclaw.State(domain,3)
    
    rc, thetac = state.p_centers
    
    state.problem_data['r0'] = r_lower
    state.problem_data['h0'] = h0
    state.problem_data['u0'] = u0
    state.problem_data['grav'] = g   # Gravitational force
    state.problem_data['kalpha'] = kalpha   # Kemm's alpha
    state.problem_data['kbeta'] = kbeta   # Kemm's beta
    state.problem_data['kepsilon'] = kepsilon   # Kemm's epsilon
    state.problem_data['F_bdy'] = F_bdy
    state.problem_data['cf'] = friction_coeff

    state.q[0,:,:] = h0
    state.q[1,:,:] = 0.
    state.q[2,:,:] = 0.

    #===========================================================================
    # Set up controller and controller parameters
    #===========================================================================
    claw = pyclaw.Controller()
    claw.tfinal = tfinal
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.num_output_times = num_output_times
    claw.outdir = outdir
    if num_cells < 400:
        claw.keep_copy = True
        #claw.output_format = None
    else:
        claw.keep_copy = False

    return claw

if __name__ == "__main__":
    from clawpack.pyclaw.util import run_app_from_main
    claw = run_app_from_main(setup)
