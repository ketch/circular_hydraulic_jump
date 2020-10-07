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

tol = 1E-7

def get_FrOutflow(hInflow,rInflow,velInflow,rOutflow,hOutflow,grav):
    beta = hInflow * rInflow * velInflow
    velOutflow = beta/(rOutflow*hOutflow)
    return velOutflow/np.sqrt(grav*hOutflow)
#
    
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

def set_jet_values(self, state, old_state=None):
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

    # in jet 
    in_jet = np.where(rc<r0)    
    state.q[0,in_jet[0],in_jet[1]] = h0
    state.q[1,in_jet[0],in_jet[1]] = h0*u0*xc[in_jet]/(rc[in_jet]+1.e-7)
    state.q[2,in_jet[0],in_jet[1]] = h0*u0*yc[in_jet]/(rc[in_jet]+1.e-7)

    #out of disk
    beta = r0*h0*u0
    g = state.problem_data['grav']
    Fr = state.problem_data['FrOutflow']
    rOutflow = state.problem_data['rOutflow']
    hOutflow = (beta/(rOutflow*Fr*np.sqrt(g)))**(2./3)
    velOutflow = beta / (rOutflow*hOutflow)
    out_disk = np.where(rc>=rOutflow)
    state.q[0,out_disk[0],out_disk[1]] = hOutflow
    state.q[1,out_disk[0],out_disk[1]] = hOutflow*velOutflow*xc[out_disk]/(rc[out_disk]+1.e-7)
    state.q[2,out_disk[0],out_disk[1]] = hOutflow*velOutflow*yc[out_disk]/(rc[out_disk]+1.e-7)
#

def subsonic_boundary_lower(state, dim, _, qbc, __, num_ghost):
    "Subsonic outflow BC at fixed Froude number.  Handles bottom/left boundaries."
    xc, yc = state.grid.p_centers_with_ghost(num_ghost)
    rc = np.sqrt(xc**2 + yc**2)
    r0 = state.problem_data['r0']
    h0 =  state.problem_data['h0']
    u0 =  state.problem_data['u0']
    beta = r0*h0*u0
    Fr = state.problem_data['FrOutflow'] # Froude number at boundary
    g = state.problem_data['grav']

    h = (beta/(rc*Fr*np.sqrt(g)))**(2./3)
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
    Fr = state.problem_data['FrOutflow'] # Froude number at boundary
    g = state.problem_data['grav']
    
    h = (beta/(rc*Fr*np.sqrt(g)))**(2./3)
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
#

def steady_rhs(h,r,beta,g=1.):
    return h/(g/beta**2 * r**3 * h**3 - r)
#
    
def initial_and_boundary_data(r_jump = 1.,
                              r_inner = 0.1,
                              r_outer = 4.0,
                              num_cells = 501,
                              g=1.,
                              h_in=0.5,
                              u_in=0.75):
    from scipy import integrate
    rc = np.linspace(r_inner,r_outer,num_cells)
    i_jump = np.argmin(np.abs(rc-r_jump))
    
    # Left boundary
    h_inner = h_in
    u_inner = u_in
    beta_inner = r_inner*h_inner*u_inner
    
    h = 0*rc
    u = 0*rc
    
    rvals = rc[:i_jump+1]
    beta = rvals[0]*h_inner*u_inner
    hh = integrate.odeint(steady_rhs,h_inner,rvals,args=(beta,g))
    hh = np.squeeze(hh)
    uu = beta/(hh*rvals)
    h[:i_jump+1] = hh[:]
    u[:i_jump+1] = uu[:]
    
    # Jump in h
    # Left side of jump
    h_m = h[i_jump]; u_m = u[i_jump]
    aleph = (-3*h_m+np.sqrt(h_m**2+8*h_m*u_m**2/g))/2.
    # Right side of jump
    h_p = h_m + aleph; u_p = h_m*u_m/h_p
    h[i_jump+1] = h_p; u[i_jump+1] = u_p
    
    # Outer part of solution
    beta_outer = rc[i_jump+1]*h[i_jump+1]*u[i_jump+1]
    rvals = rc[i_jump+1:]
    #print (beta_inner,beta_outer)
    hh = integrate.odeint(steady_rhs,h_p,rvals,args=(beta_outer,g))
    hh = np.squeeze(hh)
    uu = beta_outer/(rvals*hh)
    h[i_jump+1:] = hh[:]
    u[i_jump+1:] = uu[:]

    return rc, h, u, h_m
#

def setup(
        # about initial condition
        initialConditionType=0,
        rJump=0.25,
        # inflow parameters
        h0=0.5,
        u0=0.75,
        r0=0.1,
        # outflow parameters
        rOutflow=1.0,
        FrOutflow=0.1,
        # physical parameters
        g=1.,
        tfinal=1,
        num_output_times=10,
        friction=False,
        friction_coeff=0.01,
        boundary='subcritical',
        # numerical parameters
        num_cells=100, 
        solver_type='classic',
        riemann_solver='hlle',
        outdir='./_output',
        use_petsc=True,
        kalpha=1./3, kbeta=1.3, kepsilon=1.e-3):
    
    from clawpack import riemann
    if use_petsc:
        from clawpack import petclaw as pyclaw
    else:
        from clawpack import pyclaw
    import shallow_llf_2D

    riemann_solver = riemann_solver.lower()

    if solver_type == 'classic':
        if riemann_solver == 'hlle':
            solver = pyclaw.ClawSolver2D(riemann.shallow_hlle_2D)
            solver.fwave = False
        elif riemann_solver == 'llf':
            solver = pyclaw.ClawSolver2D(shallow_llf_2D)
            solver.num_eqn = 3
            solver.num_waves = 3
            solver.fwave = False
        else:
            raise Exception('Unrecognized Riemann solver') 
        solver.dimensional_split=True
        solver.limiters = pyclaw.limiters.tvd.minmod
        solver.cfl_max     = 0.46
        solver.cfl_desired = 0.45
        solver.transverse_waves = 2
        solver.order = 1
    elif solver_type == 'sharpclaw':
        raise Exception('No sharpclaw for now')

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

    xlower = -rOutflow;  xupper =  rOutflow
    ylower = -rOutflow;  yupper =  rOutflow

    mx = num_cells
    my = num_cells

    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    y = pyclaw.Dimension(ylower,yupper,my,name='y')
    domain = pyclaw.Domain([x,y])

    if bathymetry:
        state = pyclaw.State(domain,3,2)
        solver.aux_bc_lower[0] = pyclaw.BC.extrap
        solver.aux_bc_upper[0] = pyclaw.BC.extrap
        solver.aux_bc_lower[1] = pyclaw.BC.extrap
        solver.aux_bc_upper[1] = pyclaw.BC.extrap
    else:
        state = pyclaw.State(domain,3,1)
        solver.aux_bc_lower[0] = pyclaw.BC.extrap
        solver.aux_bc_upper[0] = pyclaw.BC.extrap
        solver.aux_bc_lower[1] = pyclaw.BC.extrap
        solver.aux_bc_upper[1] = pyclaw.BC.extrap

        state.aux[0,:,:] = 0
    old_state = pyclaw.State(domain,3)
        
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
    state.problem_data['FrOutflow'] = FrOutflow
    state.problem_data['rOutflow'] = rOutflow
    state.problem_data['cf'] = friction_coeff

    # INITIAL CONDITIONS #
    nx=xc/(rc+1.e-7); ny=yc/(rc+1.e-7)
    if initialConditionType==0:
        state.q[0,:,:] = (rc<r0)*h0 + (rc>=r0)*0.15
        state.q[1,:,:] = (rc<r0)*h0*u0*nx
        state.q[2,:,:] = (rc<r0)*h0*u0*ny
    elif initialConditionType==1:
        rStab,hStab,uStab,hL = initial_and_boundary_data(r_jump = rJump,
                                                         r_inner = r0,
                                                         r_outer = rOutflow,
                                                         num_cells = 10001,
                                                         g=g,
                                                         h_in=h0,
                                                         u_in=u0)
        from scipy.interpolate import interp1d
        hInterp = interp1d(rStab,hStab)
        rcInterp = np.minimum(np.maximum(rc,r0),rOutflow)
        state.q[0,:,:] = (rc<r0)*h0 + (rc>=r0)*hInterp(rcInterp)
        beta=r0*h0*u0
        state.q[1,:,:] = (rc<r0)*h0*u0*nx + (rc>=r0)*beta/(rc+1.e-7)*nx
        state.q[2,:,:] = (rc<r0)*h0*u0*ny + (rc>=r0)*beta/(rc+1.e-7)*ny
        # compute Froude number at outflow boundary
        velOutflow = uStab[-1] # = beta/(rOutflow*hStab[-1])
        state.problem_data['FrOutflow'] = velOutflow/np.sqrt(g*hStab[-1])
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

    claw.old_solution = pyclaw.Solution(old_state,domain)
    return claw

if __name__ == "__main__":
    from clawpack.pyclaw.util import run_app_from_main
    # ***** Physical parameters ***** #
    tfinal = 50
    nDTOut = 10
    grav=1.0
    
    # ***** numerics ***** #
    num_cells = 500
    solver_type = 'classic'
    #riemann_solver = 'hlle'
    riemann_solver = 'llf'

    # ****************************** #
    # ***** David's parameters ***** #
    # ****************************** #
    # initial condition
    initialConditionType = 0
    rJump = 0.25
    # inflow
    hInflow = 0.5
    velInflow = 0.75
    rInflow = 0.1
    # outflow
    rOutflow = 1.0
    FrOutflow = 0.1

    # **************************** #
    # ***** mql's test cases ***** #
    # **************************** #
    #velInflow = 0.75; hOutflow = 0.2961; 
    #velInflow = 1.50; hOutflow = 0.5067; 
    #velInflow = 5.00; hOutflow = 1.5827; 
    #velInflow = 10.0; hOutflow = 3.1564;
    #velInflow = 100.0, hOutflow = 31.6102;

    initialConditionType=1
    rJump = 1.0
    hInflow = 0.5
    velInflow = 5.0
    rInflow = 0.1
    rOutflow = 4.0
    hOutflow = 1.5827
    FrOutflow = get_FrOutflow(hInflow,rInflow,velInflow,rOutflow,hOutflow,grav)
    
    claw = setup(
        # about initial condition
        initialConditionType=initialConditionType,
        rJump=rJump,
        # inflow parametrs
        h0=hInflow,
        u0=velInflow,
        r0=rInflow,
        # outflow parameters
        rOutflow=rOutflow,
        FrOutflow=FrOutflow,
        # physical parameters
        g=grav,
        tfinal=tfinal,
        num_output_times = int(nDTOut*tfinal),    
        # numerical parameters
        num_cells=num_cells,
        solver_type = solver_type,
        riemann_solver=riemann_solver,        
        use_petsc=True)
    #
    
    claw.run()

