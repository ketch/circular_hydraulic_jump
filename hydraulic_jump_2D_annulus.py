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
from mapc2p import mapc2p, mapc2p_annulus

# If plotting in a Jupyter notebook, do this:
#from clawpack.visclaw.JSAnimation import IPython_display

# Mapped grid parameters
r_lower = 0.1
r_upper = 1.0
theta_lower = 0.0
theta_upper = 2*np.pi

def compute_geometry(grid):
    r"""Computes
        a_x
        a_y
        length_ratio_left
        b_x
        b_y
        length_ratio_bottom
        cell_area
    """

    dx, dy = grid.delta
    area_min = 1.e6
    area_max = 0.0

    x_corners, y_corners = grid.p_nodes

    lower_left_y, lower_left_x = y_corners[:-1,:-1], x_corners[:-1,:-1]
    upper_left_y, upper_left_x = y_corners[:-1,1: ], x_corners[:-1,1: ]
    lower_right_y, lower_right_x = y_corners[1:,:-1], x_corners[1:,:-1]
    upper_right_y, upper_right_x = y_corners[1:,1: ], x_corners[1:,1: ]

    a_x =   upper_left_y - lower_left_y  #upper left and lower left
    a_y = -(upper_left_x - lower_left_x)
    anorm = np.sqrt(a_x**2 + a_y**2)
    a_x, a_y = a_x/anorm, a_y/anorm
    length_ratio_left = anorm/dy

    b_x = -(lower_right_y - lower_left_y)  #lower right and lower left
    b_y =   lower_right_x - lower_left_x
    bnorm = np.sqrt(b_x**2 + b_y**2)
    b_x, b_y = b_x/bnorm, b_y/bnorm
    length_ratio_bottom = bnorm/dx

    area = 0*grid.c_centers[0]
    area += 0.5 * (lower_left_y+upper_left_y)*(upper_left_x-lower_left_x)
    area += 0.5 * (upper_left_y+upper_right_y)*(upper_right_x-upper_left_x)
    area += 0.5 * (upper_right_y+lower_right_y)*(lower_right_x-upper_right_x)
    area += 0.5 * (lower_right_y+lower_left_y)*(lower_left_x-lower_right_x)
    area = area/(dx*dy)
    area_min = min(area_min, np.min(area))
    area_max = max(area_max, np.max(area))

    return a_x, a_y, length_ratio_left, b_x, b_y, length_ratio_bottom, area


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
    h0 =  state.problem_data['h0']
    u0 =  state.problem_data['u0']

    xc, yc = state.grid.p_centers_with_ghost(num_ghost)
    rc = np.sqrt(xc**2 + yc**2)
    
    if dim.name == 'r':
        qbc[0,:num_ghost,:] = h0
        qbc[1,:num_ghost,:] = h0*u0*xc[:num_ghost,:]/r_lower
        qbc[2,:num_ghost,:] = h0*u0*yc[:num_ghost,:]/r_lower

def subsonic_boundary_upper(state, dim, _, qbc, __, num_ghost):
    "Subsonic outflow BC at fixed Froude number."
    xc, yc = state.grid.p_centers_with_ghost(num_ghost)
    rc = np.sqrt(xc**2 + yc**2)
    #rc, thetac = state.grid.c_centers_with_ghost(num_ghost)
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
        qbc[1,-num_ghost:,:] = h[-num_ghost:,:]*unorm[-num_ghost:,:]*xc[-num_ghost:,:]/(rc[-num_ghost:,:]+1.e-7)
        qbc[2,-num_ghost:,:] = h[-num_ghost:,:]*unorm[-num_ghost:,:]*yc[-num_ghost:,:]/(rc[-num_ghost:,:]+1.e-7)
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
          num_cells_theta=100, tfinal=10,
          solver_type='classic', num_output_times=100,
          boundary='subcritical', outdir='./_output', friction=False,
          friction_coeff=0.01, F_bdy=0.1, use_petsc=True, 
          kalpha=1./3, kbeta=1./3, kepsilon=1.e-3):
    
    import shallow_quad_hllemcc_2D
    if use_petsc:
        from clawpack import petclaw as pyclaw
    else:
        from clawpack import pyclaw
    import shallow_quad_hllemcc_2D

    if solver_type == 'classic':
        solver = pyclaw.ClawSolver2D(shallow_quad_hllemcc_2D)
        solver.cfl_max     = 0.5
        solver.cfl_desired = 0.45
        solver.num_eqn = 3
        solver.num_waves = 3
        solver.fwave = False
        solver.dimensional_split=False
        solver.transverse_waves = 1 #2
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver2D(shallow_quad_hllemcc_2D)

    solver.num_eqn = 3
    solver.num_waves = 3
    solver.fwave = False

    # Periodic BCs in theta direction
    solver.bc_lower[1] = pyclaw.BC.periodic
    solver.bc_upper[1] = pyclaw.BC.periodic
    solver.aux_bc_lower[1] = pyclaw.BC.periodic
    solver.aux_bc_upper[1] = pyclaw.BC.periodic

    solver.aux_bc_lower[0] = pyclaw.BC.extrap
    solver.aux_bc_upper[0] = pyclaw.BC.extrap

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
    theta = pyclaw.Dimension(theta_lower,theta_upper,num_cells_theta,name='theta')
    domain = pyclaw.Domain([r,theta])
    domain.grid.mapc2p = mapc2p_annulus

    state = pyclaw.State(domain,3,7)
    
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

    a_x, a_y, length_left, b_x, b_y, length_bottom, area = compute_geometry(state.grid)

    state.aux[0,:,:] = a_x
    state.aux[1,:,:] = a_y
    state.aux[2,:,:] = length_left
    state.aux[3,:,:] = b_x
    state.aux[4,:,:] = b_y
    state.aux[5,:,:] = length_bottom
    state.aux[6,:,:] = area
    state.index_capa = 6 # aux[6,:,:] holds the capacity function

    
    state.q[0,:,:] = 0.15 + 0.1*np.random.rand(*state.q[0,:,:].shape)
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
    if num_cells_r < 400:
        claw.keep_copy = True
        #claw.output_format = None
    else:
        claw.keep_copy = False
    claw.setplot = setplot

    return claw

def setplot(plotdata):
    """ 
    Plot solution using VisClaw.
    """
    from mapc2p import mapc2p
    import numpy as np
    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.mapc2p = mapc2p
    
    # Figure for contour plot
    plotfigure = plotdata.new_plotfigure(name='contour', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'q[0]'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = 0
    plotitem.contour_levels = np.linspace(0., 1., 10)
    plotitem.contour_colors = 'k'
    plotitem.patchedges_show = 1
    plotitem.MappedGrid = True

    # Figure for pcolor plot
    plotfigure = plotdata.new_plotfigure(name='q[0]', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'q[0]'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colormaps.red_yellow_blue
    plotitem.pcolor_cmin = 0.
    plotitem.pcolor_cmax = 1.
    plotitem.add_colorbar = True
    plotitem.MappedGrid = True

    # Figure for pcolor y-momentum plot
    plotfigure = plotdata.new_plotfigure(name='rho*v', figno=2)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'q[2]'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 2
    plotitem.pcolor_cmap = colormaps.red_yellow_blue
    plotitem.pcolor_cmin = -1.
    plotitem.pcolor_cmax = 1.
    plotitem.add_colorbar = True
    plotitem.MappedGrid = True


    return plotdata

if __name__ == "__main__":
    from clawpack.pyclaw.util import run_app_from_main
    claw = run_app_from_main(setup,setplot)
