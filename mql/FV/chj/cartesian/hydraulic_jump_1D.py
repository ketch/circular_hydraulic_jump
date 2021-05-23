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
from ent_residual import *
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

def get_FrOutflow(hInflow,rInflow,velInflow,rOutflow,hOutflow,grav):
    beta = hInflow * rInflow * velInflow
    velOutflow = beta/(rOutflow*hOutflow)
    return velOutflow/np.sqrt(grav*hOutflow)
#

def b4step(self, state, old_state=None):
    """
    Set solution value inside the centeral jet region.
    This should really be a boundary condition, but it's
    simpler to impose this way.
    """
    state.get_qbc_from_q(2,self.qbc)
    Ri = get_ent_residual(state,self.qbc)

    set_Ri = state.problem_data['set_Ri']
    if set_Ri is None:
        state.aux[state.num_aux-1,:,:] = Ri[2:-2,2:-2]
    else:
        state.aux[state.num_aux-1,:,:] = set_Ri
    #
#

def boundary_jet(state, dim, _, qbc, __, num_ghost):
    "Subsonic outflow BC at fixed Froude number.  Handles bottom/left boundaries."
    x = state.grid.x.centers
    dx = x[1]-x[0]
        
    r0 =  state.problem_data['r0']
    h0 =  state.problem_data['h0']
    u0 =  state.problem_data['u0']
    g =  state.problem_data['grav']
    rOutflow =  state.problem_data['rOutflow']
    beta = r0*h0*u0

    const_extrap = True
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
#

def boundary_outflow(state, dim, _, qbc, __, num_ghost):
    "Subsonic outflow BC at fixed Froude number.  Handles bottom/left boundaries."
    r0 = state.problem_data['r0']
    h0 =  state.problem_data['h0']
    u0 =  state.problem_data['u0']
    rOutflow = state.problem_data['rOutflow']    
    g = state.problem_data['grav']

    # get beta
    beta = r0*h0*u0

    const_extrap=True
    if const_extrap:
        # outflow via const extrap
        qbc[0,-num_ghost:,:] = qbc[0,-2*num_ghost:-num_ghost,:]
        qbc[1,-num_ghost:,:] = qbc[1,-2*num_ghost:-num_ghost,:]
        qbc[2,-num_ghost:,:] = qbc[2,-2*num_ghost:-num_ghost,:]
    else:
        x = state.grid.x.centers
        dx = x[1]-x[0]
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

    #plt.clf()
    #y = state.grid.y.centers
    #my = len(y)
    #xbc = state.grid.p_centers_with_ghost(2)[0][:,my/2]
    ##import pdb; pdb.set_trace()
    #plt.plot(xbc,qbc[0,:,my/2],'-r')
    #plt.plot(x,state.q[0,:,my/2])
    #plt.savefig("plot_with_BCs.png")
    ##input("stop")
    
#

def subsonic_outflow(state, dim, _, qbc, __, num_ghost):
    "Subsonic outflow BC at fixed Froude number.  Handles bottom/left boundaries."
    r0 = state.problem_data['r0']
    h0 =  state.problem_data['h0']
    u0 =  state.problem_data['u0']
    rOutflow = state.problem_data['rOutflow']    
    g = state.problem_data['grav']
    Fr = state.problem_data['FrOutflow'] # Froude number at boundary

    # get beta
    beta = r0*h0*u0

    # get h and u at outflow
    h = (beta/(rOutflow*Fr*np.sqrt(g)))**(2./3)
    u = beta / (rOutflow*h + 1E-15)

    #print (h, Fr)
    #input("asd")
    # set outflow boundary to generate a jump
    qbc[0,-num_ghost:,:] = h
    qbc[1,-num_ghost:,:] = h*u
    qbc[2,-num_ghost:,:] = 0.0
#

def steady_rhs_v2(r,h,beta):
    g=1.0
    a=g/beta**2 * r**3
    b=r
    return h/(a * h**3 - b)
#

def steady_rhs_v2_jacobian(r,h,beta):
    g=1.0
    a=g/beta**2 * r**3
    b=r
    return (-2*a*h**3-b)/(a*h**3-b)**2
#

def steady_rhs(h,r,beta,g=1.):
    return h/(g/beta**2 * r**3 * h**3 - r)
#

def steady_rhs_backward(h,r,beta,g=1.):
    return h/(g/beta**2 * (-r)**3 * h**3 - (-r))
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
        use_dmin_blended = 1.0,
        set_Ri = None,
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

    import shallow_es_2D
    import shallow_roe_with_efix_2D
    import shallow_llf_dk_2D
    
    riemann_solver = riemann_solver.lower()

    if solver_type == 'classic':
        if riemann_solver == 'es':
            solver = pyclaw.ClawSolver2D(shallow_es_2D)
            solver.num_eqn = 3
            solver.num_waves = 3
            solver.fwave = False
        elif riemann_solver == 'hlle':
            solver = pyclaw.ClawSolver2D(riemann.shallow_hlle_2D)
            solver.num_eqn = 3
            solver.num_waves = 2
            solver.fwave = False            
        elif riemann_solver == 'llf_dk':
            solver = pyclaw.ClawSolver2D(shallow_llf_dk_2D)
            solver.num_eqn = 3
            solver.num_waves = 2
            solver.fwave = False            
        elif riemann_solver == 'hlle':
            solver = pyclaw.ClawSolver2D(riemann.shallow_hlle_2D)
            solver.num_eqn = 3
            solver.num_waves = 2
            solver.fwave = False
        elif riemann_solver == 'roe':
            solver = pyclaw.ClawSolver2D(shallow_roe_with_efix_2D)
            solver.fwave = False
        else:
            raise Exception('Unrecognized Riemann solver') 
        solver.dimensional_split=True
        solver.limiters = pyclaw.limiters.tvd.minmod
        #solver.limiters = pyclaw.limiters.tvd.MC
        #solver.limiters = 0
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
    if boundary == 'outflow':
        solver.user_bc_upper = boundary_outflow
        solver.bc_upper[0] = pyclaw.BC.custom
        #solver.bc_upper[0] = pyclaw.BC.extrap
    else:
        solver.user_bc_upper = subsonic_outflow
        solver.bc_upper[0] = pyclaw.BC.custom
    #
    solver.bc_lower[1] = pyclaw.BC.wall
    solver.bc_upper[1] = pyclaw.BC.wall

    #solver.before_step = b4step
    
    xlower = 0.1;  xupper =  1.0
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
    state.problem_data['kalpha'] = kalpha   # Kemm's alpha
    state.problem_data['kbeta'] = kbeta   # Kemm's beta
    state.problem_data['kepsilon'] = kepsilon   # Kemm's epsilon
    state.problem_data['FrOutflow'] = FrOutflow
    state.problem_data['rOutflow'] = rOutflow
    state.problem_data['use_dmin_blended'] = use_dmin_blended
    state.problem_data['set_Ri'] = set_Ri

    # ***** COMPUTE EXACT SOLUTION TO USE AT BOUNDARIES ***** #    
    dx=state.grid.x.centers[1]-state.grid.x.centers[0]
    beta = r0*h0*u0
    # inner boundary
    dh=ode(steady_rhs_v2,
           steady_rhs_v2_jacobian).set_integrator('vode',
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

    #plt.plot(r_backward,hh_backward,'-r')
    #
    # outer boundary
    r_forward = np.linspace(r0,1.1*rOutflow,1000)
    hh_forward = integrate.odeint(steady_rhs,h0,r_forward,args=(beta,g))[:,0]
    hInterp = interp1d(r_forward,hh_forward[:])
    state.problem_data['hp1'] = hInterp(rOutflow+dx/2.0)
    state.problem_data['hp2'] = hInterp(rOutflow+3*dx/2.0)
    # END OF COMPUTATION OF EXACT SOLUTION AT BOUNDARIES #    

    #plt.plot(r_forward,hh_forward,'-b')
    #plt.plot(r0+dx/2.0,state.problem_data['h0'],'go',lw=10)
    #plt.plot(r0-dx/2.0,state.problem_data['hm1'],'go',lw=10)
    #plt.plot(r0-3*dx/2.0,state.problem_data['hm2'],'go',lw=10)
    #plt.savefig('plot.png')
    
    # INITIAL CONDITIONS #
    if initialConditionType==0:
        state.q[0,:,:] = 0.1
        state.q[1,:,:] = 0.0
        state.q[2,:,:] = 0.0

        #hInit = hInterp(state.p_centers[0])
        #uInit = beta/(hInit*state.p_centers[0])
        #state.q[0,:,:] = hInit
        #state.q[1,:,:] = hInit*uInit
        #state.q[2,:,:] = 0.0
    elif initialConditionType==1:
        rStab,hStab,uStab,hL = initial_and_boundary_data(r_jump = rJump,
                                                         r_inner = r0,
                                                         r_outer = rOutflow,
                                                         num_cells = 10001,
                                                         g=g,
                                                         h_in=h0,
                                                         u_in=u0)
        hInterp = interp1d(rStab,hStab)
        rcInterp = np.minimum(np.maximum(xc,r0),rOutflow)
        state.q[0,:,:] = (xc<r0)*h0 + (xc>=r0)*hInterp(rcInterp)
        beta=r0*h0*u0
        state.q[1,:,:] = beta/(xc+1.e-7)
        state.q[2,:,:] = 0 
        # compute Froude number at outflow boundary
        velOutflow = uStab[-1] 
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

    return claw

if __name__ == "__main__":
    from clawpack.pyclaw.util import run_app_from_main
    # ***** Physical parameters ***** #
    tfinal = 5.0
    nDTOut = 10
    grav = 1.0
    
    # ***** numerics ***** #
    num_cells = 50
    solver_type = 'classic'
    #solver_type = 'sharpclaw'
    #riemann_solver = 'roe'
    riemann_solver = 'es'
    #riemann_solver = 'llf_dk'
    riemann_solver = 'hlle'
    
    use_dmin_blended = 0.0
    set_Ri = 1.0
    
    # initial condition
    with_jump = False
    initialConditionType=0
    rJump = 0.3
    # inflow
    hInflow = 0.3
    velInflow = 0.75
    rInflow = 0.1
    # outflow
    rOutflow = 1.0
    hOutflow = 0.37387387318873766

    #velInflow = 10.0
    #hOutflow = 4.4502243956698901
    
    #velInflow = 15.0
    #hOutflow = 6.6845019298155357
        
    if with_jump:
        tfinal=10
        boundary='subsonic'
    else:
        tfinal=5
        boundary='outflow'
    #

    #velInflow=0.05 ####
    #tfinal=20.0
    FrOutflow = get_FrOutflow(hInflow,rInflow,velInflow,rOutflow,hOutflow,grav)
    
    #print velInflow/np.sqrt(grav*hInflow) ######
    #print FrOutflow ######
    
    claw = setup(
        set_Ri=set_Ri,
        use_dmin_blended = use_dmin_blended,
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
        boundary=boundary,
        # physical parameters
        g=grav,
        tfinal=tfinal,
        num_output_times = int(nDTOut*tfinal),
        #num_output_times = nDTOut,    
        # numerical parameters
        num_cells=num_cells,
        solver_type = solver_type,
        riemann_solver=riemann_solver,        
        use_petsc=True)
    #
    
    claw.run()

