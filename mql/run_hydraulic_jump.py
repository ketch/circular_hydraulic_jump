import sys
sys.path.append('./../')
from hydraulic_jump_2D import *

# David's Inflow parameters 
velInflow = 0.75
hInflow = 0.5
hOutflow = 0.15
rInflow = 0.1
FrOutflow = 0.1

# test cases
#velInflow = 0.75; hOutflow = 0.2961; 
#velInflow = 1.50; hOutflow = 0.5067; 
#velInflow = 5.00; hOutflow = 1.5827; 
#velInflow = 10.0; hOutflow = 3.1564;
#velInflow = 100.0, hOutflow = 31.6102;

# other parameters
tfinal = 20.0
nDTOut = 10
grav=1.0

# numerics
num_cells = 1000
solver_type = 'classic'
#riemann_solver = 'hlle'
riemann_solver = 'roe'

claw = setup(
    initialConditionType=1,
    rJump=0.25,
    # inflow parametrs
    h0=hInflow,
    u0=velInflow,
    r0=rInflow,
    # outflow parameters
    F_bdy=FrOutflow,
    h_inf=hOutflow,
    # other parameters
    g=grav,
    tfinal=tfinal,
    num_output_times = int(nDTOut*tfinal),    
    # numerics
    use_petsc=True,
    num_cells=num_cells,
    solver_type = solver_type,
    riemann_solver=riemann_solver)

claw.run()
