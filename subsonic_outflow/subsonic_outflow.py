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
from scipy import integrate
from scipy.interpolate import interp1d
from scipy.integrate import ode

def step_radial_src(solver,state,dt):
    """
    Geometric source terms for SW equations with cylindrical symmetry.
    Integrated using a 2-stage, 2nd-order Runge-Kutta method.
    This is a Clawpack-style source term routine, which approximates
    the integral of the source terms over a step.
    """
    dt2 = dt/2.
    q = state.q
    rad = state.grid.p_centers[0]

    h = q[0,:,:]
    u = q[1,:,:]/h

    qstar = np.empty(q.shape)
    qstar[0,:,:] = q[0,:,:] - dt2/rad * h*u
    qstar[1,:,:] = q[1,:,:] - dt2/rad * h*u*u

    h = qstar[0,:]
    u = qstar[1,:]/h
    
    q[0,:,:] = q[0,:,:] - dt/rad * h*u
    q[1,:,:] = q[1,:,:] - dt/rad * h*u*u
#

def dq_swe_radial(solver,state,dt):
    q   = state.q
    rad = state.grid.p_centers[0]

    h = q[0,:,:]
    u = q[1,:,:]/h
    dq = np.empty(q.shape)
    dq[0,:,:] = - 1.0/rad * h*u
    dq[1,:,:] = - 1.0/rad * h*u*u
    dq[2,:,:] = 0

    return dq
#

def b4step(self, state, old_state=None):
    pass
#

def boundary_jet(state, dim, _, qbc, __, num_ghost):
    x = state.grid.x.centers
    dx = x[1]-x[0]
        
    r0 =  state.problem_data['r0']
    h0 =  state.problem_data['h0']
    u0 =  state.problem_data['u0']
    beta = r0*h0*u0

    const_extrap = False
    if const_extrap:
        qbc[0,:num_ghost,:] = h0
        qbc[1,:num_ghost,:] = h0*u0
        qbc[2,:num_ghost,:] = 0.0        
    else:
        rm1 = r0-dx/2.0
        rm2 = r0-3*dx/2.0
        
        # hm1
        hm1 = state.problem_data['hm1'] 
        um1 = beta/(rm1*hm1)
        qbc[0,1,:] = hm1
        qbc[1,1,:] = hm1*um1
        qbc[2,1,:] = 0.0

        # hm2
        hm2 = state.problem_data['hm2'] 
        um2 = beta/(rm2*hm2)
        qbc[0,0,:] = hm2
        qbc[1,0,:] = hm2*um2
        qbc[2,0,:] = 0.0
    #
#

def boundary_outflow(state, dim, _, qbc, __, num_ghost):
    x = state.grid.x.centers
    dx = x[1]-x[0]
    
    r0 = state.problem_data['r0']
    h0 =  state.problem_data['h0']
    u0 =  state.problem_data['u0']
    rOutflow = state.problem_data['rOutflow']    
    beta = r0*h0*u0

    const_extrap=False
    if const_extrap:
        # outflow via const extrap
        qbc[0,-num_ghost:,:] = qbc[0,-2*num_ghost:-num_ghost,:]
        qbc[1,-num_ghost:,:] = qbc[1,-2*num_ghost:-num_ghost,:]
        qbc[2,-num_ghost:,:] = qbc[2,-2*num_ghost:-num_ghost,:]
    else:
        rp1 = rOutflow + dx/2.0
        rp2 = rOutflow + 3*dx/2.0
        
        # hp1
        hp1 = state.problem_data['hp1']
        up1 = beta/(rp1*hp1)
        qbc[0,-2,:] = hp1
        qbc[1,-2,:] = hp1*up1
        qbc[2,-2,:] = 0.0

        # hp2
        hp2 = state.problem_data['hp2']
        up2 = beta/(rp2*hp2)
        qbc[0,-1,:] = hp2
        qbc[1,-1,:] = hp1*up2
        qbc[2,-1,:] = 0.0
    #
#

def steady_rhs(h,r,beta,g=1.):
    return h/(g/beta**2 * r**3 * h**3 - r)
#

def steady_rhs_backward(h,r,beta,g=1.):
    a=g/beta**2 * r**3
    b=r
    return h/(a * h**3 - b)
#

def steady_rhs_backward_jacobian(r,h,beta,g=1.0):
    a=g/beta**2 * r**3
    b=r
    return (-2*a*h**3-b)/(a*h**3-b)**2
#

def setup(
        # inflow parameters
        h0=0.5,
        u0=0.75,
        r0=0.1,
        # outflow parameters
        rOutflow=1.0,
        # physical parameters
        g=1.,
        tfinal=1,
        num_output_times=10,
        # numerical parameters
        num_cells=100, 
        solver_type='classic',
        riemann_solver='hlle',
        outdir='./_output',
        use_petsc=True):
    
    from clawpack import riemann
    if use_petsc:
        from clawpack import petclaw as pyclaw
    else:
        from clawpack import pyclaw

    import shallow_llf_dk_2D
    
    riemann_solver = riemann_solver.lower()

    if solver_type == 'classic':
        if riemann_solver == 'hlle':
            solver = pyclaw.ClawSolver2D(riemann.shallow_hlle_2D)
            solver.num_eqn = 3
            solver.num_waves = 2
            solver.fwave = False            
        elif riemann_solver == 'llf_dk':
            solver = pyclaw.ClawSolver2D(shallow_llf_dk_2D)
            solver.num_eqn = 3
            solver.num_waves = 2
            solver.fwave = False            
        else:
            raise Exception('Unrecognized Riemann solver') 
        solver.dimensional_split=True
        solver.limiters = pyclaw.limiters.tvd.minmod
        solver.cfl_max     = 0.46
        solver.cfl_desired = 0.45
        solver.transverse_waves = 2
        solver.order = 2

        # geometric source term
        solver.step_source = step_radial_src
        solver.source_split = 2
        solver.before_step = b4step
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver2D(shallow_llf_dk_2D)
        solver.num_eqn = 3
        solver.num_waves = 2
        solver.fwave = False
        
        solver.dq_src = dq_swe_radial
        solver.weno_order = 5
        solver.lim_type   = 2
        
        solver.cfl_max     = 0.46
        solver.cfl_desired = 0.45
    #
    solver.user_bc_lower = boundary_jet
    solver.bc_lower[0] = pyclaw.BC.custom
    solver.user_bc_upper = boundary_outflow
    solver.bc_upper[0] = pyclaw.BC.custom
    solver.bc_lower[1] = pyclaw.BC.wall
    solver.bc_upper[1] = pyclaw.BC.wall
    
    xlower = r0;  xupper =  rOutflow
    ylower = 0;  yupper =  0.1

    mx = num_cells
    my = int(num_cells/10.0)

    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    y = pyclaw.Dimension(ylower,yupper,my,name='y')
    domain = pyclaw.Domain([x,y])

    state = pyclaw.State(domain,3,1)

    solver.aux_bc_lower[0] = pyclaw.BC.extrap
    solver.aux_bc_upper[0] = pyclaw.BC.extrap
    solver.aux_bc_lower[1] = pyclaw.BC.wall
    solver.aux_bc_upper[1] = pyclaw.BC.wall
    state.aux[0,:,:] = 0
    
    xc, yc = state.p_centers

    state.problem_data['r0'] = r0
    state.problem_data['h0'] = h0
    state.problem_data['u0'] = u0
    state.problem_data['grav'] = g   # Gravitational force
    state.problem_data['rOutflow'] = rOutflow

    # ***** COMPUTE EXACT SOLUTION TO USE AT BOUNDARIES ***** #    
    dx=state.grid.x.centers[1]-state.grid.x.centers[0]
    beta = r0*h0*u0
    # inner boundary
    dh=ode(steady_rhs_backward,
           steady_rhs_backward_jacobian).set_integrator('vode',
                                                        method='bdf',
                                                        nsteps=1E5)
    dh.set_initial_value(h0,r0).set_f_params(beta).set_jac_params(beta)
    rEnd=r0-3*dx/2.0 # for 2 ghost cells
    dr=1E-6*dx
    r_backward=np.linspace(r0-dr,rEnd,1000)
    hh_backward = r_backward*0
    for i in range(len(r_backward)):
        hh_backward[i] = dh.integrate(r_backward[i])[0]
    #    
    hInterp = interp1d(r_backward,hh_backward)
    state.problem_data['hm1'] = hInterp(r0-dx/2.0)
    state.problem_data['hm2'] = hInterp(r0-3*dx/2.0)

    # outer boundary
    r_forward = np.linspace(r0,rOutflow+2*dx,1000)
    hh_forward = integrate.odeint(steady_rhs,h0,r_forward,args=(beta,g))[:,0]
    hInterp = interp1d(r_forward,hh_forward[:])
    state.problem_data['hp1'] = hInterp(rOutflow+dx/2.0)
    state.problem_data['hp2'] = hInterp(rOutflow+3*dx/2.0)
    # END OF COMPUTATION OF EXACT SOLUTION AT BOUNDARIES #    

    # ***** INITIAL CONDITIONS ***** #
    hInit = hInterp(state.p_centers[0])
    uInit = beta/(hInit*state.p_centers[0])
    state.q[0,:,:] = hInit
    state.q[1,:,:] = hInit*uInit
    state.q[2,:,:] = 0.0
    #==============================================
    # Set up controller and controller parameters
    #==============================================
    claw = pyclaw.Controller()
    claw.tfinal = tfinal
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.num_output_times = num_output_times
    claw.outdir = outdir

    return claw

if __name__ == "__main__":
    from clawpack.pyclaw.util import run_app_from_main
    # ***** Physical parameters ***** #
    tfinal = 10.0
    nDTOut = 10
    grav = 1.0
    
    # ***** numerics ***** #
    use_petsc = False
    num_cells = 50
    solver_type = 'classic'
    #solver_type = 'sharpclaw'
    #riemann_solver = 'llf_dk'
    riemann_solver = 'hlle'
    
    # inflow 
    hInflow = 0.3
    velInflow = 0.05
    rInflow = 0.1
    # outflow
    rOutflow = 1.0
    
    claw = setup(
        # inflow parametrs
        h0=hInflow,
        u0=velInflow,
        r0=rInflow,
        # outflow parameters
        rOutflow=rOutflow,
        # physical parameters
        g=grav,
        tfinal=tfinal,
        num_output_times = int(nDTOut*tfinal),
        # numerical parameters
        num_cells=num_cells,
        solver_type = solver_type,
        riemann_solver=riemann_solver,        
        use_petsc=use_petsc)
    #
    
    claw.run()

