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
from ent_residual import *

tol = 1E-20

def b4step(self, state,old_state=None):
    state.get_qbc_from_q(2,self.qbc)
    Ri = get_ent_residual(state,self.qbc)

    set_Ri = state.problem_data['set_Ri']
    if set_Ri is None:
        state.aux[state.num_aux-1,:,:] = Ri[2:-2,2:-2]
    else:
        state.aux[state.num_aux-1,:,:] = set_Ri
    #
#
def setup(
        use_dmin_blended = 1.0,
        set_Ri = None,
        # physical parameters
        g=1.,
        tfinal=1,
        num_output_times=10,
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
    #import shallow_hllc_2D
    import shallow_hllemcc_2D
    import shallow_hllemccRoEF_2D
    import shallow_roe_exp_2D
    import shallow_llf_2D
    import shallow_ev_2D
    import shallow_es_2D        
    
    
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
        elif riemann_solver == 'ev':
            solver = pyclaw.ClawSolver2D(shallow_ev_2D)
            solver.num_eqn = 3
            solver.num_waves = 3
            solver.fwave = False
        elif riemann_solver == 'es':
            solver = pyclaw.ClawSolver2D(shallow_es_2D)
            solver.num_eqn = 3
            solver.num_waves = 3
            solver.fwave = False
        elif riemann_solver == 'roe_exp':
            solver = pyclaw.ClawSolver2D(shallow_roe_exp_2D)
            solver.num_eqn = 3
            solver.num_waves = 3
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
        solver.dimensional_split=True
        solver.limiters = pyclaw.limiters.tvd.minmod
        solver.cfl_max     = 0.46   #0.9
        solver.cfl_desired = 0.45   #0.8
        solver.transverse_waves = 2
        solver.order = 2
    elif solver_type == 'sharpclaw':
        raise Exception('Do not use sharpclaw')

    bathymetry = False
    solver.bc_lower[0] = pyclaw.BC.wall
    solver.bc_upper[0] = pyclaw.BC.wall
    solver.bc_lower[1] = pyclaw.BC.wall
    solver.bc_upper[1] = pyclaw.BC.wall
    
    solver.before_step = b4step

    xlower = 0.0;  xupper =  10.0
    ylower = 0.0;  yupper =  1.0

    my = int(num_cells/10.0)
    mx = num_cells

    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    y = pyclaw.Dimension(ylower,yupper,my,name='y')
    domain = pyclaw.Domain([x,y])

    state = pyclaw.State(domain,3,1)
    solver.aux_bc_lower[0] = pyclaw.BC.extrap
    solver.aux_bc_upper[0] = pyclaw.BC.extrap
    solver.aux_bc_lower[1] = pyclaw.BC.extrap
    solver.aux_bc_upper[1] = pyclaw.BC.extrap
    state.aux[0,:,:] = 0

    xc, yc = state.p_centers

    state.problem_data['grav'] = g   # Gravitational force
    state.problem_data['kalpha'] = kalpha   # Kemm's alpha
    state.problem_data['kbeta'] = kbeta   # Kemm's beta
    state.problem_data['kepsilon'] = kepsilon   # Kemm's epsilon
    state.problem_data['use_dmin_blended'] = use_dmin_blended
    state.problem_data['set_Ri'] = set_Ri
    
    # INITIAL CONDITIONS #
    x0=5.0
    hl=0.005
    hr=0.001    
    #hr=1E-15

    t=1.0; x0=5.0
    xA=x0-t*np.sqrt(g*hl)
    xB=x0+2*t*np.sqrt(g*hl)
    h = (xc<=xA) * hl + (xc>xA)*(xc<xB) * 4.0/(9*g)*(np.sqrt(g*hl)-(xc-x0)/(2*t))**2
    u = (xc>xA)*(xc<xB) * 2.0/3*((xc-x0)/t+np.sqrt(g*hl))
        
    state.q[0,:,:] = (xc<=x0)*hl + (xc>x0)*hr
    state.q[1,:,:] = 0.0
    state.q[2,:,:] = 0.0

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
    # ***** Physical parameters ***** #
    tfinal = 10
    nDTOut = 1
    grav=1.0
    
    # ***** numerics ***** #
    num_cells = 400
    solver_type = 'classic'
    #riemann_solver = 'hlle'
    #riemann_solver = 'roe'
    #riemann_solver = 'roe_exp'
    #riemann_solver = 'hllemcc'
    #riemann_solver = 'llf'
    #riemann_solver = 'ev'
    riemann_solver = 'es'

    use_dmin_blended = 1.0
    set_Ri = None
        
    claw = setup(
        set_Ri=set_Ri,
        use_dmin_blended = use_dmin_blended,
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
