"""
Options:
    - Riemann solver: GeoClaw, Roe, HLLE, HLLEMCC
"""

from matplotlib import animation
import matplotlib.pyplot as plt
from clawpack import pyclaw
import numpy as np
from ent_residual import *

tol = 1E-7

# If plotting in a Jupyter notebook, do this:
#from clawpack.visclaw.JSAnimation import IPython_display

# *********************************************************** #
# ********** FUNCTIONS TO COMPUTE ENTROPY RESIDUAL ********** #
# *********************************************************** #
def b4step(self, state, old_state=None):
    """
    Set zero normal boundary inside the cylinder
    """
    # get entropy residual #
    if use_petsc:
        state.get_qbc_from_q(2,self.qbc)
        Ri = get_ent_residual(state,self.qbc)
    else:
        Ri = get_ent_residual(state)
    #
    
    set_Ri = state.problem_data['set_Ri']
    if set_Ri is None:
        if use_petsc:
            state.aux[state.num_aux-1,:,:] = Ri[2:-2,2:-2]
        else:
            state.aux[state.num_aux-1,:,:] = Ri
    else:
        state.aux[state.num_aux-1,:,:] = set_Ri
    #

    # set values at the cylinder #
    xc, yc = state.p_centers
    rc = np.sqrt((xc-40)**2 + (yc-50)**2)
    in_solid = np.where((rc <= 20.0))

    nx = (xc-40)/(rc+1e-7)
    ny = (yc-50)/(rc+1e-7)
    hu=state.q[1,:,:]
    hv=state.q[2,:,:]
    
    hu_tan = hu*ny - hv*nx
    state.q[1,in_solid[0],in_solid[1]] =  ny[in_solid[0],in_solid[1]] * hu_tan[in_solid[0],in_solid[1]]
    state.q[2,in_solid[0],in_solid[1]] = -nx[in_solid[0],in_solid[1]] * hu_tan[in_solid[0],in_solid[1]]
#

def inflow_bcs(state, dim, _, qbc, __, num_ghost):
    qbc[0,:num_ghost,:] = state.problem_data['h0']
    qbc[1,:num_ghost,:] = state.problem_data['h0']*state.problem_data['u0']
    qbc[2,:num_ghost,:] = 0.0
#

def setup(
        use_dmin_blended = 1.0,
        set_Ri = None,
        h0=2.0,
        u0=20.0,
        g=1.,
        x_num_cells=160,
        y_num_cells=400,
        tfinal=10,
        num_output_times=100,
        solver_type='classic',
        riemann_solver='roe',
        outdir='./_output',
        use_petsc=True,
        restart_from_frame=None,
        file_prefix='claw',
        kalpha=1./3,
        kbeta=1.3,
        kepsilon=1.e-3):
    
    from clawpack import riemann
    if use_petsc:
        from clawpack import petclaw as pyclaw
    else:
        from clawpack import pyclaw
    import shallow_es_2D
    import shallow_hllemcc_2D
    
    riemann_solver = riemann_solver.lower()

    if solver_type == 'classic':
        if riemann_solver == 'es':
            solver = pyclaw.ClawSolver2D(shallow_es_2D)
            solver.num_eqn = 3
            solver.num_waves = 3
            solver.fwave = False
            solver.dimensional_split=True
            solver.transverse_waves = 2
        elif riemann_solver == 'hllemcc':
            solver = pyclaw.ClawSolver2D(shallow_hllemcc_2D)
            solver.num_eqn = 3
            solver.num_waves = 3
            solver.fwave = False
            solver.dimensional_split=False
            solver.transverse_waves = 1
        else:
            raise Exception('Unrecognized Riemann solver')
        solver.limiters = pyclaw.limiters.tvd.minmod
        solver.cfl_max     = 0.46
        solver.cfl_desired = 0.45
        solver.order = 2
    elif solver_type == 'sharpclaw':
        raise Exception('use solver_type classic for now')

    solver.bc_lower[0] = pyclaw.BC.custom
    solver.bc_upper[0] = pyclaw.BC.extrap
    solver.bc_lower[1] = pyclaw.BC.extrap
    solver.bc_upper[1] = pyclaw.BC.extrap
    solver.user_bc_lower = inflow_bcs
    solver.before_step = b4step

    #===========================================================================
    # Set up controller and controller parameters
    #===========================================================================
    claw = pyclaw.Controller()
    claw.tfinal = tfinal
    claw.solver = solver
    claw.num_output_times = num_output_times
    claw.outdir = outdir
    claw.keep_copy = False

    # mesh resolution
    mx = x_num_cells
    my = y_num_cells

    if restart_from_frame is not None:
        claw.solution = pyclaw.Solution(restart_from_frame, file_format='petsc',read_aux=False,file_prefix=file_prefix)
        grid = claw.solution.domain.grid
        claw.num_output_times = claw.num_output_times - restart_from_frame
        claw.start_frame = restart_from_frame
    else:
        # domain
        xlower = 0.0;  xupper =  40.0
        ylower = 0.0;  yupper =  100.0

        x = pyclaw.Dimension(xlower,xupper,mx,name='x')
        y = pyclaw.Dimension(ylower,yupper,my,name='y')
        domain = pyclaw.Domain([x,y])
        
        state = pyclaw.State(domain,3,1)
        solver.aux_bc_lower[0] = pyclaw.BC.extrap
        solver.aux_bc_upper[0] = pyclaw.BC.extrap
        solver.aux_bc_lower[1] = pyclaw.BC.extrap
        solver.aux_bc_upper[1] = pyclaw.BC.extrap
        state.aux[0,:,:] = 0
        
        old_state = pyclaw.State(domain,3,1)
        xc, yc = state.p_centers    
        rc = np.sqrt(xc**2 + yc**2)
        # solid

        state.problem_data['use_dmin_blended'] = use_dmin_blended
        state.problem_data['set_Ri'] = set_Ri
        state.problem_data['use_petsc'] = use_petsc
        state.problem_data['grav'] = g   # Gravitational force
        state.problem_data['h0'] = h0
        state.problem_data['u0'] = u0
        state.problem_data['kalpha'] = kalpha   # Kemm's alpha
        state.problem_data['kbeta'] = kbeta   # Kemm's beta
        state.problem_data['kepsilon'] = kepsilon   # Kemm's epsilon
        
        # INITIAL CONDITIONS #
        nx=xc/(rc+1.e-7); ny=yc/(rc+1.e-7)
        state.q[0,:,:] = state.problem_data['h0']
        state.q[1,:,:] = state.problem_data['h0']*state.problem_data['u0']
        state.q[2,:,:] = 0.0

        # create clawpack solution
        claw.solution = pyclaw.Solution(state,domain)
        claw.old_solution = pyclaw.Solution(old_state,domain)

    return claw

if __name__ == "__main__":
    from clawpack.pyclaw.util import run_app_from_main

    # ********************** #
    # ***** parameters ***** #
    # ********************** #
    tfinal = 100
    nDTOut = 1
    h0 = 1.0
    u0 = 5.0
    g = 1.0
    
    # ******************** #
    # ***** numerics ***** #
    # ******************** #
    # riemann solver
    riemann_solver = 'es'
    #riemann_solver = 'hllemcc'
    
    use_dmin_blended = 1.0
    set_Ri = None
    use_petsc = True
    
    # mesh
    refn = 0
    mx=int((2**refn)*160)
    my=int((2**refn)*400)
    
    claw = setup(
        use_petsc=use_petsc,
        use_dmin_blended = use_dmin_blended,
        set_Ri = set_Ri,
        # parameters
        tfinal=tfinal,
        num_output_times=int(nDTOut*tfinal),
        h0=h0,
        u0=u0,
        g=g,
        restart_from_frame=None,
        # numerics
        riemann_solver=riemann_solver,
        x_num_cells=mx,
        y_num_cells=my)
    claw.run()

