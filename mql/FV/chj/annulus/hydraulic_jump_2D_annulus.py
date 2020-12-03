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
theta_lower = 0.0
theta_upper = 2*np.pi
tol=1E-15

# ENTROPY AND ITS DERIVATIVES
def get_entropy(g,h,hu,hv):
    return 0.5*g*h**2 + 0.5*(hu**2+hv**2)/(h+tol)
def dEntropy_dh(g,h,hu,hv):
    return g*h - 0.5*(hu**2+hv**2)/(h**2+tol)
def dEntropy_dhu(g,h,hu,hv):
    return hu/(h+tol)
def dEntropy_dhv(g,h,hu,hv):
    return hv/(h+tol)
# FLUXES OF THE SHALLOW WATER EQUATIONS
def xFlux1(g,h,hu,hv):
    return hu
def xFlux2(g,h,hu,hv):
    return hu**2/(h+tol) + 0.5*g*h**2
def xFlux3(g,h,hu,hv):
    return hu*hv/(h+tol)
def yFlux1(g,h,hu,hv):
    return hv
def yFlux2(g,h,hu,hv):
    return hu*hv/(h+tol)
def yFlux3(g,h,hu,hv):
    return hv**2/(h+tol) + 0.5*g*h**2
# ENTROPY FLUX 
def xEntFlux(g,h,hu,hv):
    return (0.5*(hu**2+hv**2)/(h+tol) + g*h**2)*hu/(h+tol) 
def yEntFlux(g,h,hu,hv):
    return (0.5*(hu**2+hv**2)/(h+tol) + g*h**2)*hv/(h+tol)
def sumIntOverPatch(I):
    pass
#I[1:-1,1:-1] = (I[1:-1,1:-1] +
#                    I[:-2,1:-1] + I[2:,1:-1] +
#                    I[1:-1,:-2] + I[1:-1,2:])
#

def fix_parallel_BCs(state,num_ghost,qbc):
    num_eqn=1
    q_da = state._create_DA(num_eqn,num_ghost)
    gqVec = q_da.createGlobalVector()
    _q_local_vector = q_da.createLocalVector()

    gqVec.setArray(qbc[num_ghost:-num_ghost,num_ghost:-num_ghost].reshape([-1], order = 'F'))    
    shape = [n + 2*num_ghost for n in state.grid.num_cells]
    q_da.globalToLocal(gqVec, _q_local_vector)

    qbc[:] = _q_local_vector.getArray().reshape(shape, order = 'F')
#

def get_entropy_residual(state,qbc=None):
    # ************************************ #
    # ***** compute entropy residual ***** #
    # ************************************ #
    dx, dy = state.grid.delta
    g=state.problem_data['grav']

    if qbc is not None:
        hij  = qbc[0,:,:]
        huij = qbc[1,:,:]
        hvij = qbc[2,:,:]
    else:
        hij  = state.q[0,:,:]
        huij = state.q[1,:,:]
        hvij = state.q[2,:,:]
    #

    hip1j = np.zeros_like(hij); huip1j = np.zeros_like(huij); hvip1j = np.zeros_like(hvij)
    him1j = np.zeros_like(hij); huim1j = np.zeros_like(huij); hvim1j = np.zeros_like(hvij)
    hijp1 = np.zeros_like(hij); huijp1 = np.zeros_like(huij); hvijp1 = np.zeros_like(hvij)
    hijm1 = np.zeros_like(hij); huijm1 = np.zeros_like(huij); hvijm1 = np.zeros_like(hvij)

    # get the solution on the faces
    hip1j[:-1,:] = hij[1:,:]; huip1j[:-1,:] = huij[1:,:]; hvip1j[:-1,:] = hvij[1:,:]
    him1j[1:,:] = hij[:-1,:]; huim1j[1:,:] = huij[:-1,:]; hvim1j[1:,:] = hvij[:-1,:]    
    hijp1[:,:-1] = hij[:,1:]; huijp1[:,:-1] = huij[:,1:]; hvijp1[:,:-1] = hvij[:,1:]
    hijm1[:,1:] = hij[:,:-1]; huijm1[:,1:] = huij[:,:-1]; hvijm1[:,1:] = hvij[:,:-1]
    
    # fix boundary conditions
    if qbc is not None:
        fix_parallel_BCs(state,2,hip1j)
        fix_parallel_BCs(state,2,huip1j)
        fix_parallel_BCs(state,2,hvip1j)
        #
        fix_parallel_BCs(state,2,him1j)
        fix_parallel_BCs(state,2,huim1j)
        fix_parallel_BCs(state,2,hvim1j)
        #
        fix_parallel_BCs(state,2,hijp1)
        fix_parallel_BCs(state,2,huijp1)
        fix_parallel_BCs(state,2,hvijp1)
        #
        fix_parallel_BCs(state,2,hijm1)
        fix_parallel_BCs(state,2,huijm1)
        fix_parallel_BCs(state,2,hvijm1)
    else:
        hip1j[-1,:] = hij[0,:];  huip1j[-1,:] = huij[0,:];  hvip1j[-1,:] = hvij[0,:]
        him1j[0,:]  = hij[-1,:]; huim1j[0,:]  = huij[-1,:]; hvim1j[0,:]  = hvij[-1,:]
        hijp1[:,-1] = hij[:,0];  huijp1[:,-1] = huij[:,0];  hvijp1[:,-1] = hvij[:,0]
        hijm1[:,0]  = hij[:,-1]; huijm1[:,0]  = huij[:,-1]; hvijm1[:,0]  = hvij[:,-1]
    #
    
    # get the entropy fluxes
    xFij   = xEntFlux(g,hij,huij,hvij)
    xFip1j = xEntFlux(g,hip1j,huip1j,hvip1j)
    xFim1j = xEntFlux(g,him1j,huim1j,hvim1j)
    xFijp1 = xEntFlux(g,hijp1,huijp1,hvijp1)
    xFijm1 = xEntFlux(g,hijm1,huijm1,hvijm1)
    
    yFij   = yEntFlux(g,hij,huij,hvij)
    yFip1j = yEntFlux(g,hip1j,huip1j,hvip1j)
    yFim1j = yEntFlux(g,him1j,huim1j,hvim1j)
    yFijp1 = yEntFlux(g,hijp1,huijp1,hvijp1)
    yFijm1 = yEntFlux(g,hijm1,huijm1,hvijm1)

    # physical fluxes
    # term 1
    xf1ij   = xFlux1(g,hij,huij,hvij)
    xf1ip1j = xFlux1(g,hip1j,huip1j,hvip1j)
    xf1im1j = xFlux1(g,him1j,huim1j,hvim1j)
    xf1ijp1 = xFlux1(g,hijp1,huijp1,hvijp1)
    xf1ijm1 = xFlux1(g,hijm1,huijm1,hvijm1)

    yf1ij   = yFlux1(g,hij,huij,hvij)
    yf1ip1j = yFlux1(g,hip1j,huip1j,hvip1j)
    yf1im1j = yFlux1(g,him1j,huim1j,hvim1j)
    yf1ijp1 = yFlux1(g,hijp1,huijp1,hvijp1)
    yf1ijm1 = yFlux1(g,hijm1,huijm1,hvijm1)

    # term 2
    xf2ij   = xFlux2(g,hij,huij,hvij)
    xf2ip1j = xFlux2(g,hip1j,huip1j,hvip1j)
    xf2im1j = xFlux2(g,him1j,huim1j,hvim1j)
    xf2ijp1 = xFlux2(g,hijp1,huijp1,hvijp1)
    xf2ijm1 = xFlux2(g,hijm1,huijm1,hvijm1)

    yf2ij   = yFlux2(g,hij,huij,hvij)
    yf2ip1j = yFlux2(g,hip1j,huip1j,hvip1j)
    yf2im1j = yFlux2(g,him1j,huim1j,hvim1j)
    yf2ijp1 = yFlux2(g,hijp1,huijp1,hvijp1)
    yf2ijm1 = yFlux2(g,hijm1,huijm1,hvijm1)

    # term 3
    xf3ij   = xFlux3(g,hij,huij,hvij)
    xf3ip1j = xFlux3(g,hip1j,huip1j,hvip1j)
    xf3im1j = xFlux3(g,him1j,huim1j,hvim1j)
    xf3ijp1 = xFlux3(g,hijp1,huijp1,hvijp1)
    xf3ijm1 = xFlux3(g,hijm1,huijm1,hvijm1)

    yf3ij   = yFlux3(g,hij,huij,hvij)
    yf3ip1j = yFlux3(g,hip1j,huip1j,hvip1j)
    yf3im1j = yFlux3(g,him1j,huim1j,hvim1j)
    yf3ijp1 = yFlux3(g,hijp1,huijp1,hvijp1)
    yf3ijm1 = yFlux3(g,hijm1,huijm1,hvijm1)
    
    # ***** right face ***** #
    nx = np.zeros_like(huij); ny = np.zeros_like(huij); ds = np.zeros_like(huij)
    nx_tmp = np.zeros_like(huij[2:-2,2:-2]); ny_tmp = np.zeros_like(huij[2:-2,2:-2]); ds_tmp = np.zeros_like(huij[2:-2,2:-2]);
    if qbc is None:
        nx[:-1,:] = state.aux[0,1:,:]; ny[:-1,:] = state.aux[1,1:,:]; ds[:-1,:] = state.aux[2,1:,:]*dy
        nx[-1,:] = state.aux[0,0,:]; ny[-1,:] = state.aux[1,0,:]; ds[-1,:] = state.aux[2,0,:]*dy
    else:
        nx_tmp[:-1,:] = state.aux[0,1:,:]; ny_tmp[:-1,:] = state.aux[1,1:,:]; ds_tmp[:-1,:] = state.aux[2,1:,:]*dy
        # the following line is a hack to avoid extra comm in parallel
        nx_tmp[-1,:] = nx_tmp[-2,:]; ny_tmp[-1,:] = ny_tmp[-2,:]; ds_tmp[-1,:] = ds_tmp[-2,:]
        nx[2:-2,2:-2] = nx_tmp; ny[2:-2,2:-2] = ny_tmp; ds[2:-2,2:-2] = ds_tmp;  
        fix_parallel_BCs(state,2,nx)
        fix_parallel_BCs(state,2,ny)
        fix_parallel_BCs(state,2,ds)
    #
    intEntFluxRight = ds*(nx*0.5*(xFij  + xFip1j)  + ny*0.5*(yFij + yFip1j))
    intFlux1Right   = ds*(nx*0.5*(xf1ij + xf1ip1j) + ny*0.5*(yf1ij + yf1ip1j))
    intFlux2Right   = ds*(nx*0.5*(xf2ij + xf2ip1j) + ny*0.5*(yf2ij + yf2ip1j))
    intFlux3Right   = ds*(nx*0.5*(xf3ij + xf3ip1j) + ny*0.5*(yf3ij + yf3ip1j))

    # ***** left face ***** #
    if qbc is None:
        nx[:,:] = -state.aux[0,:,:]; ny[:,:] = -state.aux[1,:,:]; ds[:,:] = state.aux[2,:,:]*dy
    else:
        nx_tmp[:,:] = -state.aux[0,:,:]; ny_tmp[:,:] = -state.aux[1,:,:]; ds_tmp[:,:] = state.aux[2,:,:]*dy
        nx[2:-2,2:-2] = nx_tmp; ny[2:-2,2:-2] = ny_tmp; ds[2:-2,2:-2] = ds_tmp;  
        fix_parallel_BCs(state,2,nx)
        fix_parallel_BCs(state,2,ny)
        fix_parallel_BCs(state,2,ds)
    #
    intEntFluxLeft = ds*(nx*0.5*(xFij  + xFim1j)  + ny*0.5*(yFij  + yFim1j))
    intFlux1Left   = ds*(nx*0.5*(xf1ij + xf1im1j) + ny*0.5*(yf1ij + yf1im1j))
    intFlux2Left   = ds*(nx*0.5*(xf2ij + xf2im1j) + ny*0.5*(yf2ij + yf2im1j))
    intFlux3Left   = ds*(nx*0.5*(xf3ij + xf3im1j) + ny*0.5*(yf3ij + yf3im1j))
    
    # ***** upper face ***** #
    if qbc is None:
        nx[:,:-1] = state.aux[3,:,1:]; ny[:,:-1] = state.aux[4,:,1:]; ds[:,:-1] = state.aux[5,:,1:]*dx
        nx[:,-1] = state.aux[3,:,0]; ny[:,-1] = state.aux[4,:,0]; ds[:,-1] = state.aux[5,:,0]*dx
    else:
        nx_tmp[:,:-1] = state.aux[3,:,1:]; ny_tmp[:,:-1] = state.aux[4,:,1:]; ds_tmp[:,:-1] = state.aux[5,:,1:]*dx
        # the following line is a hack to avoid extra comm in parallel runs
        nx_tmp[:,-1] = nx_tmp[:,-2]; ny_tmp[:,-1] = ny_tmp[:,-2]; ds_tmp[:,-1] = ds_tmp[:,-2]
        nx[2:-2,2:-2] = nx_tmp; ny[2:-2,2:-2] = ny_tmp; ds[2:-2,2:-2] = ds_tmp;  
        fix_parallel_BCs(state,2,nx)
        fix_parallel_BCs(state,2,ny)
        fix_parallel_BCs(state,2,ds)
    #
    intEntFluxUpper = ds*(nx*0.5*(xFij + xFijp1) + ny*0.5*(yFij + yFijp1))
    intFlux1Upper = ds*(nx*0.5*(xf1ij + xf1ijp1) + ny*0.5*(yf1ij + yf1ijp1))
    intFlux2Upper = ds*(nx*0.5*(xf2ij + xf2ijp1) + ny*0.5*(yf2ij + yf2ijp1))
    intFlux3Upper = ds*(nx*0.5*(xf3ij + xf3ijp1) + ny*0.5*(yf3ij + yf3ijp1))
    
    # ***** bottom face ***** #
    if qbc is None:
        nx[:,:] = -state.aux[3,:,:]; ny[:,:] = -state.aux[4,:,:]; ds[:,:] = state.aux[5,:,:]*dx
    else:
        nx_tmp[:,:] = -state.aux[3,:,:]; ny_tmp[:,:] = -state.aux[4,:,:]; ds_tmp[:,:] = state.aux[5,:,:]*dx
        nx[2:-2,2:-2] = nx_tmp; ny[2:-2,2:-2] = ny_tmp; ds[2:-2,2:-2] = ds_tmp;  
        fix_parallel_BCs(state,2,nx)
        fix_parallel_BCs(state,2,ny)
        fix_parallel_BCs(state,2,ds)
    #
    intEntFluxBottom = ds*(nx*0.5*(xFij + xFijm1) + ny*0.5*(yFij + yFijm1))
    intFlux1Bottom = ds*(nx*0.5*(xf1ij + xf1ijm1) + ny*0.5*(yf1ij + yf1ijm1))
    intFlux2Bottom = ds*(nx*0.5*(xf2ij + xf2ijm1) + ny*0.5*(yf2ij + yf2ijm1))
    intFlux3Bottom = ds*(nx*0.5*(xf3ij + xf3ijm1) + ny*0.5*(yf3ij + yf3ijm1))
    
    # ***** integral of the divergence of the flux ***** #
    int_div_flux_term1 = intFlux1Right + intFlux1Upper + intFlux1Left + intFlux1Bottom
    int_div_flux_term2 = intFlux2Right + intFlux2Upper + intFlux2Left + intFlux2Bottom
    int_div_flux_term3 = intFlux3Right + intFlux3Upper + intFlux3Left + intFlux3Bottom

    # ***** integral of the divergence of the entropy flux ***** #
    int_div_ent_flux = intEntFluxRight + intEntFluxUpper + intEntFluxLeft + intEntFluxBottom
        
    sumIntOverPatch(int_div_flux_term1)
    sumIntOverPatch(int_div_flux_term2)
    sumIntOverPatch(int_div_flux_term3)
    sumIntOverPatch(int_div_ent_flux)

    if qbc is not None:
        fix_parallel_BCs(state,2,int_div_flux_term1)
        fix_parallel_BCs(state,2,int_div_flux_term2)
        fix_parallel_BCs(state,2,int_div_flux_term3)
        fix_parallel_BCs(state,2,int_div_ent_flux)
    #

    # derivative of the entropy
    eta_prime_term1 = dEntropy_dh (g,hij,huij,hvij)
    eta_prime_term2 = dEntropy_dhu(g,hij,huij,hvij)
    eta_prime_term3 = dEntropy_dhv(g,hij,huij,hvij)
    
    # get dot product of the derivative of entropy and the int of the div of the flux
    eta_prime_times_int_div_flux = (eta_prime_term1 * int_div_flux_term1 +
                                    eta_prime_term2 * int_div_flux_term2 +
                                    eta_prime_term3 * int_div_flux_term3)
    
    # get l2 norm of the derivative of the entropy
    abs_eta_prime = np.sqrt(eta_prime_term1**2 + eta_prime_term2**2 + eta_prime_term3**2) 

    # get the l2 norm of the integral of the divergence of the flux
    abs_int_div_flux = np.sqrt(int_div_flux_term1**2 + int_div_flux_term2**2 + int_div_flux_term3**2)

    aux = (  abs(eta_prime_term1) * abs(int_div_flux_term1)
             + abs(eta_prime_term2) * abs(int_div_flux_term2)
             + abs(eta_prime_term3) * abs(int_div_flux_term3) )
    
    # entropy residual 
    Num = np.abs(int_div_ent_flux - eta_prime_times_int_div_flux) 
    #Den = np.abs(int_div_ent_flux) + abs_eta_prime * abs_int_div_flux + 1E-4
    Den = np.abs(int_div_ent_flux) + aux + 1E-14
    #Den = np.abs(int_div_ent_flux) + np.abs(eta_prime_times_int_div_flux) + 1E-4

    set_Ri = state.problem_data['set_Ri']
    if set_Ri is None:
        Ri = (Num/Den)
    else:
        Ri = Num*0 + set_Ri
    #
    
    fix_parallel_BCs(state,2,Ri)
    return Ri 
#

def b4step(self, state, old_state=None):
    """
    Set solution value inside the centeral jet region.
    This should really be a boundary condition, but it's
    simpler to impose this way.
    """
    state.get_qbc_from_q(2,self.qbc)
    Ri = get_entropy_residual(state,self.qbc)
    state.aux[state.num_aux-1,:,:] = Ri[2:-2,2:-2]
#

def get_FrOutflow(hInflow,rInflow,velInflow,rOutflow,hOutflow,grav):
    beta = hInflow * rInflow * velInflow
    velOutflow = beta/(rOutflow*hOutflow)
    return velOutflow/np.sqrt(grav*hOutflow)
#

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

def jet(state, dim, _, qbc, __, num_ghost):
    "Jet inflow BC at inner boundary."
    r0 = state.problem_data['r0']
    h0 =  state.problem_data['h0']
    u0 =  state.problem_data['u0']
    
    xc, yc = state.grid.p_centers_with_ghost(num_ghost)
    rc = np.sqrt(xc**2 + yc**2)

    rand_inflow = state.problem_data['rand_inflow']
    eps = (np.random.rand(qbc[0,:num_ghost,:].shape[0],qbc[0,:num_ghost,:].shape[1])
           -
           np.random.rand(qbc[0,:num_ghost,:].shape[0],qbc[0,:num_ghost,:].shape[1]))/100.0 * (rand_inflow)    
    u = u0*(1.0+eps)
    h = h0/(1.0+eps)
            
    if dim.name == 'r':
        qbc[0,:num_ghost,:] = h
        qbc[1,:num_ghost,:] = h*u*xc[:num_ghost,:]/r0
        qbc[2,:num_ghost,:] = h*u*yc[:num_ghost,:]/r0

def subsonic_boundary_upper(state, dim, _, qbc, __, num_ghost):
    "Subsonic outflow BC at fixed Froude number."
    xc, yc = state.grid.p_centers_with_ghost(num_ghost)
    rc = np.sqrt(xc**2 + yc**2)
    #rc, thetac = state.grid.c_centers_with_ghost(num_ghost)
    r0 = state.problem_data['r0']
    h0 =  state.problem_data['h0']
    u0 =  state.problem_data['u0']
    beta = r0*h0*u0
    F = state.problem_data['FrOutflow'] # Froude number at boundary
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
        rand_inflow = 0.0,
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
        tfinal=10,
        num_output_times=100,
        friction=False,
        friction_coeff=0.01,
        boundary='subcritical',
        # numerical parameters
        num_cells_r=100,
        num_cells_theta=100, 
        solver_type='classic', 
        riemann_solver='es',
        outdir='./_output', 
        use_petsc=True, 
        kalpha=1./3, kbeta=1./3, kepsilon=1.e-3):
    
    import shallow_annular_hllemcc_2D
    import shallow_annular_es_2D

    if use_petsc:
        from clawpack import petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type == 'classic':
        if riemann_solver == 'es':
            solver = pyclaw.ClawSolver2D(shallow_annular_es_2D)
            solver.dimensional_split = True
            solver.transverse_waves = 2
        elif riemann_solver == 'hllemcc':
            solver = pyclaw.ClawSolver2D(shallow_annular_hllemcc_2D)
            solver.dimensional_split = False
            solver.transverse_waves  = 1
        else:
            raise Exception('Unrecognized Riemann solver')
        #
    #
    else:
        raise Exception('set solver_type to classic')
    #
    solver.limiters = pyclaw.limiters.tvd.minmod
    solver.cfl_max     = 0.46   #0.9
    solver.cfl_desired = 0.45   #0.8
    solver.order = 2
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
        solver.bc_upper[0] = pyclaw.BC.extrap
    elif boundary == 'subcritical':  # subcritical boundary condition
        solver.bc_upper[0] = pyclaw.BC.custom
        solver.user_bc_upper = subsonic_boundary_upper
    elif boundary == 'wall':
        solver.bc_upper[0] = pyclaw.BC.wall

    solver.before_step = b4step

    if friction:
        solver.step_source = step_friction
        solver.source_split = 1

    #===========================================================================
    # Set up controller and controller parameters
    #===========================================================================
    claw = pyclaw.Controller()
    claw.tfinal = tfinal
    #
    claw.solver = solver
    claw.num_output_times = num_output_times
    claw.outdir = outdir
    if num_cells_r < 400:
        claw.keep_copy = True
        #claw.output_format = None
    else:
        claw.keep_copy = False
    #

    restart_from_frame = None
    if restart_from_frame is not None:
        claw.solution = pyclaw.Solution(restart_from_frame, file_format='petsc',read_aux=False)
        claw.num_output_times = claw.num_output_times - restart_from_frame
        claw.start_frame = restart_from_frame
        grid = claw.solution.domain.grid
        # set aux
        a_x, a_y, length_left, b_x, b_y, length_bottom, area = compute_geometry(grid)        
        claw.solution.state.aux[0,:,:] = a_x
        claw.solution.state.aux[1,:,:] = a_y
        claw.solution.state.aux[2,:,:] = length_left
        claw.solution.state.aux[3,:,:] = b_x
        claw.solution.state.aux[4,:,:] = b_y
        claw.solution.state.aux[5,:,:] = length_bottom
        claw.solution.state.aux[6,:,:] = area
        # set start frame
        claw.solution._start_frame = restart_from_frame
    else:
        # domain
        r = pyclaw.Dimension(r0,rOutflow,num_cells_r,name='r')
        theta = pyclaw.Dimension(theta_lower,theta_upper,num_cells_theta,name='theta')
        domain = pyclaw.Domain([r,theta])
        domain.grid.mapc2p = mapc2p_annulus
        
        # state
        state = pyclaw.State(domain,3,8) # 3 waves and 8 auxiliary variables
        
        xc, yc = state.p_centers
        rc = np.sqrt(xc**2 + yc**2)
        nx=xc/(rc+1.e-7); ny=yc/(rc+1.e-7)
        
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
        state.problem_data['use_dmin_blended'] = use_dmin_blended
        state.problem_data['set_Ri'] = set_Ri
        state.problem_data['rand_inflow'] = rand_inflow
        
        a_x, a_y, length_left, b_x, b_y, length_bottom, area = compute_geometry(state.grid)        
        state.aux[0,:,:] = a_x
        state.aux[1,:,:] = a_y
        state.aux[2,:,:] = length_left
        state.aux[3,:,:] = b_x
        state.aux[4,:,:] = b_y
        state.aux[5,:,:] = length_bottom
        state.aux[6,:,:] = area
        state.index_capa = 6 # aux[6,:,:] holds the capacity function
        
        # INITIAL CONDITIONS #
        if initialConditionType==0:
            state.q[0,:,:] = 0.1
            state.q[1,:,:] = 0.
            state.q[2,:,:] = 0.
        else:
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
            state.q[0,:,:] = hInterp(rcInterp)
            beta=r0*h0*u0
            #hunorm = beta/rc
            state.q[1,:,:] = beta/(rc+1.E-7)*nx
            state.q[2,:,:] = beta/(rc+1.E-7)*ny
        #        
        claw.solution = pyclaw.Solution(state,domain)
    #
    return claw

if __name__ == "__main__":

    # ***** Physical parameters ***** #
    tfinal = 20 #tmp
    nDTOut = 10
    grav=1.0

    # ***** numerics ***** #
    num_cells_r = 200
    num_cells_theta = 200
    solver_type = 'classic'

    # ********************** #
    # ***** PARAMETERS ***** #
    # ********************** #
    # rand perturbation in the boundaries (in percentage)
    rand_inflow = 0.0
    test_case = 5
    riemann_solver = 'es'
    #riemann_solver = 'hllemcc'
    use_dmin_blended = 1.0
    set_Ri = None
    
    # **************************** #
    # ***** mql's test cases ***** #
    # **************************** #
    # initial condition
    initialConditionType = 1
    rJump = 0.3
    # inflow
    hInflow = 0.3
    rInflow = 0.1
    # outflow
    rOutflow = 1.0

    if test_case == 1:
        # Test case 1
        velInflow = 0.75
        hOutflow = 0.37387387318873766
    elif test_case == 2:
        # Test case 2
        velInflow = 1.50
        hOutflow = 0.67946669652810898
    elif test_case == 3:
        # Test case 3
        velInflow = 5.0
        hOutflow = 2.2195774944760704
    elif test_case == 4:
        # Test case 4
        velInflow = 10.0
        hOutflow = 4.4502243956698901
    else:
        # test case 5
        velInflow = 15.0
        hOutflow = 6.6845019298155357
    #
    FrOutflow = get_FrOutflow(hInflow,rInflow,velInflow,rOutflow,hOutflow,grav)
    
    claw = setup(
        rand_inflow = rand_inflow,
        use_dmin_blended = use_dmin_blended,
        set_Ri = set_Ri,
        # about initial condition
        initialConditionType=initialConditionType,
        rJump=rJump,
        # inflow parametrs
        h0=hInflow,
        u0=velInflow,
        r0=rInflow,
        # outflow parameters
        FrOutflow=FrOutflow,
        rOutflow=rOutflow,
        #boundary='outflow',
        # physical parameters
        g=grav,
        tfinal=tfinal,
        num_output_times = int(nDTOut*tfinal),    
        # numerical parameters
        num_cells_r=num_cells_r,
        num_cells_theta=num_cells_theta,
        solver_type = solver_type,
        riemann_solver=riemann_solver,        
        use_petsc=True)
    
    claw.run()

        
