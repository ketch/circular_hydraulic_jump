"""
Code to simulate 2D circular hydraulic jump in PyClaw.

Currently uses a Cartesian grid.

In order to get a jump, you must either introduce a bathymetric barrier (circular ring),
include bottom friction, or impose a subcritical flow at the outer boundary.
To work with bathymetry, you have to use a Riemann solver that handles variable
bathymetry; for now this means the GeoClaw solver.

We could also get a jump by introducing bottom friction; we haven't tried this yet.

Options:
    - Riemann solver: GeoClaw, Roe, HLLE, HLLEMCC
"""

from matplotlib import animation
import matplotlib.pyplot as plt
from clawpack import pyclaw
import numpy as np

# If plotting in a Jupyter notebook, do this:
#from clawpack.visclaw.JSAnimation import IPython_display

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

def set_jet_values(_, state):
    """
    Set solution value inside the centeral jet region.
    This should really be a boundary condition, but it's
    simpler to impose this way.
    """
    xc, yc = state.p_centers
    rc = np.sqrt(xc**2 + yc**2)
    r0 = state.problem_data['r0']
    h0 =  state.problem_data['h0']
    u0 =  state.problem_data['u0']

    in_jet = np.where(rc<r0)
    
    state.q[0,in_jet[0],in_jet[1]] = h0
    state.q[1,in_jet[0],in_jet[1]] = h0*u0*xc[in_jet]/(rc[in_jet]+1.e-7)
    state.q[2,in_jet[0],in_jet[1]] = h0*u0*yc[in_jet]/(rc[in_jet]+1.e-7)

def subsonic_boundary_lower(state, dim, _, qbc, __, num_ghost):
    "Subsonic outflow BC at fixed Froude number.  Handles bottom/left boundaries."
    xc, yc = state.grid.p_centers_with_ghost(num_ghost)
    rc = np.sqrt(xc**2 + yc**2)
    r0 = state.problem_data['r0']
    h0 =  state.problem_data['h0']
    u0 =  state.problem_data['u0']
    beta = r0*h0*u0
    F = state.problem_data['F_bdy'] # Froude number at boundary
    g = state.problem_data['grav']
    
    h = (beta/(rc*F*np.sqrt(g)))**(2./3)
    unorm = beta / (rc*h)
    
    if dim.name == 'x':
        qbc[0,:num_ghost,:] = h[:num_ghost,:]
        qbc[1,:num_ghost,:] = h[:num_ghost,:]*unorm[:num_ghost,:]*xc[:num_ghost,:]/rc[:num_ghost,:]
        qbc[2,:num_ghost,:] = h[:num_ghost,:]*unorm[:num_ghost,:]*yc[:num_ghost,:]/rc[:num_ghost,:]       
    elif dim.name == 'y':
        qbc[0,:,:num_ghost] = h[:,:num_ghost]
        qbc[1,:,:num_ghost] = h[:,:num_ghost]*unorm[:,:num_ghost]*xc[:,:num_ghost]/rc[:,:num_ghost]
        qbc[2,:,:num_ghost] = h[:,:num_ghost]*unorm[:,:num_ghost]*yc[:,:num_ghost]/rc[:,:num_ghost]

def subsonic_boundary_upper(state, dim, _, qbc, __, num_ghost):
    "Subsonic outflow at BC fixed Froude number.  Handles top/right boundaries."
    xc, yc = state.grid.p_centers_with_ghost(num_ghost)
    rc = np.sqrt(xc**2 + yc**2)
    r0 = state.problem_data['r0']
    h0 =  state.problem_data['h0']
    u0 =  state.problem_data['u0']
    beta = r0*h0*u0
    F = state.problem_data['F_bdy'] # Froude number at boundary
    g = state.problem_data['grav']
    
    h = (beta/(rc*F*np.sqrt(g)))**(2./3)
    unorm = beta / (rc*h)
    
    if dim.name == 'x':
        qbc[0,-num_ghost:,:] = h[-num_ghost:,:]
        qbc[1,-num_ghost:,:] = h[-num_ghost:,:]*unorm[-num_ghost:,:]*xc[-num_ghost:,:]/rc[-num_ghost:,:]
        qbc[2,-num_ghost:,:] = h[-num_ghost:,:]*unorm[-num_ghost:,:]*yc[-num_ghost:,:]/rc[-num_ghost:,:]     
    elif dim.name == 'y':
        qbc[0,:,-num_ghost:] = h[:,-num_ghost:]
        qbc[1,:,-num_ghost:] = h[:,-num_ghost:]*unorm[:,-num_ghost:]*xc[:,-num_ghost:]/rc[:,-num_ghost:]
        qbc[2,:,-num_ghost:] = h[:,-num_ghost:]*unorm[:,-num_ghost:]*yc[:,-num_ghost:]/rc[:,-num_ghost:]
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
        
def setup(h0=0.5, u0=0.75, r0=0.1, h_inf=0.15, g=1., num_cells=100, tfinal=1,
          solver_type='classic', num_output_times=10, riemann_solver='hlle',
          boundary='subcritical', outdir='./_output', friction=False,
          friction_coeff=0.01, F_bdy=0.1, use_petsc=False, 
          kalpha=1./3, kbeta=1.3, kepsilon=1.e-3):
    
    from clawpack import riemann
    if use_petsc:
        from clawpack import petclaw as pyclaw
    else:
        from clawpack import pyclaw
    #import shallow_hllc_2D
    import shallow_hllemcc_2D
    import shallow_hllemccRoEF_2D

    riemann_solver = riemann_solver.lower()
    print(riemann_solver)

    if solver_type == 'classic':
        if riemann_solver == 'hlle':
            solver = pyclaw.ClawSolver2D(riemann.shallow_hlle_2D)
            solver.fwave = False
        elif riemann_solver == 'roe':
            solver = pyclaw.ClawSolver2D(riemann.shallow_roe_with_efix_2D)
            solver.fwave = False
        elif riemann_solver == 'hllc':
            solver = pyclaw.ClawSolver2D(shallow_hllc_2D)
            solver.num_eqn = 3
            solver.num_waves = 3
            solver.fwave = False
        elif riemann_solver == 'hllemcc':
            solver = pyclaw.ClawSolver2D(shallow_hllemcc_2D)
            solver.num_eqn = 3
            solver.num_waves = 3
            solver.fwave = False
        elif riemann_solver == 'hllemccroef':
            solver = pyclaw.ClawSolver2D(shallow_hllemccRoEF_2D)
            solver.num_eqn = 3
            solver.num_waves = 3
            solver.fwave = False
        elif riemann_solver == 'hllem':
            # Not yet implemented
            solver = pyclaw.ClawSolver2D(riemann.shallow_hllem_2D)
            solver.fwave = False
        elif riemann_solver == 'geoclaw':
            solver = pyclaw.ClawSolver2D(riemann.sw_aug_2D)
            solver.fwave = True
        else:
            raise Exception('Unrecognized Riemann solver') 
        #solver.dimensional_split=True
        #solver.limiters = pyclaw.limiters.tvd.minmod
        solver.cfl_max     = 0.9
        solver.cfl_desired = 0.8
        solver.transverse_waves = 2
    elif solver_type == 'sharpclaw':
        if riemann_solver == 'hlle':
            solver = pyclaw.SharpClawSolver2D(riemann.shallow_hlle_2D)
        elif riemann_solver == 'hllem':
            # not yet implemented
            solver = pyclaw.SharpClawSolver2D(riemann.shallow_hllem_2D)
        elif riemann_solver == 'geoclaw':
            solver = pyclaw.SharpClawSolver2D(riemann.sw_aug_2D)
        else:
            raise Exception('Riemann solver must be hlle or geoclaw')

    bathymetry = False

    if boundary == 'outflow':
        # Use bathymetry to create hydrualic jump
        if friction == False:
            assert riemann_solver == 'geoclaw'
            bathymetry = True  # Otherwise, no jump

        solver.bc_lower[0] = pyclaw.BC.extrap
        solver.bc_upper[0] = pyclaw.BC.extrap
        solver.bc_lower[1] = pyclaw.BC.extrap
        solver.bc_upper[1] = pyclaw.BC.extrap
    elif boundary == 'subcritical':  # subcritical boundary condition
        solver.bc_lower[0] = pyclaw.BC.custom
        solver.bc_upper[0] = pyclaw.BC.custom
        solver.bc_lower[1] = pyclaw.BC.custom
        solver.bc_upper[1] = pyclaw.BC.custom
        solver.user_bc_lower = subsonic_boundary_lower
        solver.user_bc_upper = subsonic_boundary_upper
    elif boundary == 'wall':
        solver.bc_lower[0] = pyclaw.BC.wall
        solver.bc_upper[0] = pyclaw.BC.wall
        solver.bc_lower[1] = pyclaw.BC.wall
        solver.bc_upper[1] = pyclaw.BC.wall

    solver.before_step = set_jet_values

    if friction:
        solver.step_source = step_friction
        solver.source_split = 1

    xlower = -1;  xupper =  1.
    ylower = -1;  yupper =  1.

    mx = num_cells
    my = num_cells

    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    y = pyclaw.Dimension(ylower,yupper,my,name='y')
    domain = pyclaw.Domain([x,y])

    if bathymetry:
        state = pyclaw.State(domain,3,1)
        solver.aux_bc_lower[0] = pyclaw.BC.extrap
        solver.aux_bc_upper[0] = pyclaw.BC.extrap
        solver.aux_bc_lower[1] = pyclaw.BC.extrap
        solver.aux_bc_upper[1] = pyclaw.BC.extrap
    else:
        state = pyclaw.State(domain,3)
    
    xc, yc = state.p_centers
    
    rc = np.sqrt(xc**2 + yc**2)
    if bathymetry:
        state.aux[0,:,:] = 0.1*np.exp(-300*(rc-0.8)**2) - 0.2*np.exp(-1000*(rc-0.75)**2)

    state.problem_data['r0'] = r0
    state.problem_data['h0'] = h0
    state.problem_data['u0'] = u0
    state.problem_data['grav'] = g   # Gravitational force
    state.problem_data['kalpha'] = kalpha   # Kemm's alpha
    state.problem_data['kbeta'] = kbeta   # Kemm's beta
    state.problem_data['kepsilon'] = kepsilon   # Kemm's epsilon
    state.problem_data['F_bdy'] = F_bdy
    state.problem_data['cf'] = friction_coeff

    state.q[0,:,:] = (rc<r0)*h0 + (rc>=r0)*0.15
    state.q[1,:,:] = (rc<r0)*h0*u0*xc/(rc+1.e-7)
    state.q[2,:,:] = (rc<r0)*h0*u0*yc/(rc+1.e-7)

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
    #claw = setup(h0=0.5, u0=0.75,riemann_solver='hllemccroef',num_output_times=100,tfinal=20,num_cells=100,friction=True)
    #claw = setup(h0=0.5, u0=0.75,riemann_solver='hlle',num_output_times=100,tfinal=20,num_cells=100,friction=True)    
    #claw.run()
    
    # If plotting in a Jupyter notebook
    #anim=plot_surface(claw)
    #from IPython.display import HTML
    #HTML(anim.to_html5_video())
    
    # If plotting in the terminal    
    #from clawpack.pyclaw import plot
    #plot.interactive_plot()

    # To write plots to disk:
    #import sys
    #print(sys.argv)
    #plotdir = sys.argv[1]
    #plot_surface(claw, make_anim=False, save_plots=True, frames=len(claw.frames), plotdir=plotdir)

#-----------
