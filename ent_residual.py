import numpy as np
tol = 1E-15
#tol = 1E-7

# *********************************************************** #
# ********** FUNCTIONS TO COMPUTE ENTROPY RESIDUAL ********** #
# *********************************************************** #
def reconstruct_solution(h):
    x_hL = np.zeros_like(h)
    x_hR = np.zeros_like(h)
    x_hL[2:-2,:] = (-3.*h[:-4,:] + 27*h[1:-3,:] + 47*h[2:-2,:] - 13*h[3:-1,:] + 2*h[4:,:])/60.
    x_hR[2:-2,:] = (2.*h[:-4,:] - 13*h[1:-3,:] + 47*h[2:-2,:] + 27*h[3:-1,:] - 3*h[4:,:])/60.

    y_hL = np.zeros_like(h)
    y_hR = np.zeros_like(h)
    y_hL[:,2:-2] = (-3.*h[:,:-4] + 27*h[:,1:-3] + 47*h[:,2:-2] - 13*h[:,3:-1] + 2*h[:,4:])/60.
    y_hR[:,2:-2] = (2.*h[:,:-4] - 13*h[:,1:-3] + 47*h[:,2:-2] + 27*h[:,3:-1] - 3*h[:,4:])/60.

    return x_hL, x_hR, y_hL, y_hR
#
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
# DERIVATIVES OF ENTROPY FLUX 
def xEntFlux(g,h,hu,hv):
    return (0.5*(hu**2+hv**2)/(h+tol) + g*h**2)*hu/(h+tol) 
def yEntFlux(g,h,hu,hv):
    return (0.5*(hu**2+hv**2)/(h+tol) + g*h**2)*hv/(h+tol)
# SUM INTEGRALS OVER PATCH
def sumIntOverPatch(I):
    pass
    #I[1:-1,1:-1] = (I[1:-1,1:-1] +
    #                I[:-2,1:-1] + I[2:,1:-1] +
    #                I[1:-1,:-2] + I[1:-1,2:])
# *********************************************************** #

def fix_parallel_BCs(state,num_ghost,qbc):
    num_eqn=1
    q_da = state._create_DA(num_eqn,num_ghost)
    gqVec = q_da.createGlobalVector()
    _q_local_vector = q_da.createLocalVector()
    
    gqVec.setArray(qbc[num_ghost:-num_ghost,num_ghost:-num_ghost].reshape([-1], order = 'F'))
    shape = [n + 2*num_ghost for n in state.grid.num_cells]
    q_da.globalToLocal(gqVec, _q_local_vector)
    return _q_local_vector.getArray().reshape(shape, order = 'F')
#
    
# ****************************************************************************** #
# ********** GET ENTROPY RESIDUAL BASED ON HIGH-ORDER RECONSTRUCTIONS ********** #
# ****************************************************************************** #
def get_ent_residual_based_on_reconstructions(state,qbc=None):
    # ************************************ #
    # ***** compute entropy residual ***** #
    # ************************************ #
    g  = state.problem_data['grav']
    if qbc is not None:
        h  = qbc[0,:,:]
        q1 = qbc[1,:,:]
        q2 = qbc[2,:,:]
    else:
        h  = state.q[0,:,:]
        q1 = state.q[1,:,:]
        q2 = state.q[2,:,:]
    #
    dx = state.delta[0]
    dy = state.delta[1]
        
    # reconstruct solution at the interfaces
    x_hL, x_hR, y_hL, y_hR = reconstruct_solution(h)
    x_huL, x_huR, y_huL, y_huR = reconstruct_solution(q1)
    x_hvL, x_hvR, y_hvL, y_hvR = reconstruct_solution(q2)

    
    ###############
    if qbc is not None:
        # fix h
        x_hL = fix_parallel_BCs(state,2,x_hL)
        x_hR = fix_parallel_BCs(state,2,x_hR)
        y_hL = fix_parallel_BCs(state,2,y_hL)
        y_hR = fix_parallel_BCs(state,2,y_hR)
        # fix hu
        x_huL = fix_parallel_BCs(state,2,x_huL)
        x_huR = fix_parallel_BCs(state,2,x_huR)
        y_huL = fix_parallel_BCs(state,2,y_huL)
        y_huR = fix_parallel_BCs(state,2,y_huR)
        # fix hv
        x_hvL = fix_parallel_BCs(state,2,x_hvL)
        x_hvR = fix_parallel_BCs(state,2,x_hvR)
        y_hvL = fix_parallel_BCs(state,2,y_hvL)
        y_hvR = fix_parallel_BCs(state,2,y_hvR)
    ###################

    # integral of the divergence of the entropy flux
    int_div_ent_flux = ( (xEntFlux(g,x_hR,x_huR,x_hvR) - xEntFlux(g,x_hL,x_huL,x_hvL))*dy +
                         (yEntFlux(g,y_hR,y_huR,y_hvR) - yEntFlux(g,y_hL,y_huL,y_hvL))*dx )
    sumIntOverPatch(int_div_ent_flux)
    # fix parallel BCs
    if qbc is not None:
        int_div_ent_flux = fix_parallel_BCs(state,2,int_div_ent_flux) 
    #
    # derivative of the entropy 
    eta_prime_term1 = dEntropy_dh(g,h,q1,q2)
    eta_prime_term2 = dEntropy_dhu(g,h,q1,q2)
    eta_prime_term3 = dEntropy_dhv(g,h,q1,q2)

    # integral of the divergence of the flux 
    int_div_flux_term1 = ( (xFlux1(g,x_hR,x_huR,x_hvR) - xFlux1(g,x_hL,x_huL,x_hvL))*dy +
                           (yFlux1(g,y_hR,y_huR,y_hvR) - yFlux1(g,y_hL,y_huL,y_hvL))*dx )
    int_div_flux_term2 = ( (xFlux2(g,x_hR,x_huR,x_hvR) - xFlux2(g,x_hL,x_huL,x_hvL))*dy +
                           (yFlux2(g,y_hR,y_huR,y_hvR) - yFlux2(g,y_hL,y_huL,y_hvL))*dx )
    int_div_flux_term3 = ( (xFlux3(g,x_hR,x_huR,x_hvR) - xFlux3(g,x_hL,x_huL,x_hvL))*dy +
                           (yFlux3(g,y_hR,y_huR,y_hvR) - yFlux3(g,y_hL,y_huL,y_hvL))*dx )
    sumIntOverPatch(int_div_flux_term1)
    sumIntOverPatch(int_div_flux_term2)
    sumIntOverPatch(int_div_flux_term3)
    if qbc is not None:
        # fix parallel BCs
        int_div_flux_term1 = fix_parallel_BCs(state,2,int_div_flux_term1)
        int_div_flux_term2 = fix_parallel_BCs(state,2,int_div_flux_term2)
        int_div_flux_term3 = fix_parallel_BCs(state,2,int_div_flux_term3)
    #
    # get dot product of the derivative of entropy and the int of the div of the flux
    eta_prime_times_int_div_flux = (eta_prime_term1 * int_div_flux_term1 +
                                    eta_prime_term2 * int_div_flux_term2 +
                                    eta_prime_term3 * int_div_flux_term3)

    # get l2 norm of the derivative of the entropy
    abs_eta_prime = np.sqrt(eta_prime_term1**2 + eta_prime_term2**2 + eta_prime_term3**2) 

    # get the l2 norm of the integral of the divergence of the flux
    abs_int_div_flux = np.sqrt(int_div_flux_term1**2 + int_div_flux_term2**2 + int_div_flux_term3**2)
    
    # entropy residual
    Num = np.abs(int_div_ent_flux - eta_prime_times_int_div_flux)
    aux = np.abs(eta_prime_term1*int_div_flux_term1) + np.abs(eta_prime_term2*int_div_flux_term2) + np.abs(eta_prime_term3*int_div_flux_term3)
        
    #Den = np.abs(int_div_ent_flux) + np.abs(eta_prime_times_int_div_flux) + tol
    Den  = np.abs(int_div_ent_flux) + aux + tol
    #Den = np.abs(int_div_ent_flux) + abs_eta_prime * abs_int_div_flux + tol

    Ri = (Num/Den)

    if True:
        # set Ri at left boundary to zero
        if qbc is not None:
            xc = state.grid.c_centers_with_ghost(2)[0]
        else:
            xc = state.grid.c_centers[0]
        #
        Ri[:,:] = (xc >= 0.1 + dx) * Ri
    ##################################
    
    return Ri
#

# ********************************************************************** #
# ********** GET ENTROPY RESIDUAL BASED ON DIVERGENCE THEOREM ********** #
# ********************************************************************** #
def get_ent_residual_based_on_div_thm(state,qbc=None):
    # ************************************ #
    # ***** compute entropy residual ***** #
    # ************************************ #
    dx = state.delta[0]
    dy = state.delta[1]
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
    hip1j[-1,:] = hij[0,:];  huip1j[-1,:] = huij[0,:];  hvip1j[-1,:] = hvij[0,:]
    him1j[0,:]  = hij[-1,:]; huim1j[0,:]  = huij[-1,:]; hvim1j[0,:]  = hvij[-1,:]
    hijp1[:,-1] = hij[:,0];  huijp1[:,-1] = huij[:,0];  hvijp1[:,-1] = hvij[:,0]
    hijm1[:,0]  = hij[:,-1]; huijm1[:,0]  = huij[:,-1]; hvijm1[:,0]  = hvij[:,-1]

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
    
    # right face
    nx = 1.0; ny = 0.0; ds = dy
    intEntFluxRight = ds*(nx*0.5*(xFij  + xFip1j)  + ny*0.5*(yFij + yFip1j))
    intFlux1Right   = ds*(nx*0.5*(xf1ij + xf1ip1j) + ny*0.5*(yf1ij + yf1ip1j))
    intFlux2Right   = ds*(nx*0.5*(xf2ij + xf2ip1j) + ny*0.5*(yf2ij + yf2ip1j))
    intFlux3Right   = ds*(nx*0.5*(xf3ij + xf3ip1j) + ny*0.5*(yf3ij + yf3ip1j))
    
    # left face
    nx = -1.0; ny = 0.0; ds = dy
    intEntFluxLeft = ds*(nx*0.5*(xFij  + xFim1j)  + ny*0.5*(yFij  + yFim1j))
    intFlux1Left   = ds*(nx*0.5*(xf1ij + xf1im1j) + ny*0.5*(yf1ij + yf1im1j))
    intFlux2Left   = ds*(nx*0.5*(xf2ij + xf2im1j) + ny*0.5*(yf2ij + yf2im1j))
    intFlux3Left   = ds*(nx*0.5*(xf3ij + xf3im1j) + ny*0.5*(yf3ij + yf3im1j))

    # upper face
    nx = 0.0; ny = 1.0; ds = dx
    intEntFluxUpper = ds*(nx*0.5*(xFij + xFijp1) + ny*0.5*(yFij + yFijp1))
    intFlux1Upper = ds*(nx*0.5*(xf1ij + xf1ijp1) + ny*0.5*(yf1ij + yf1ijp1))
    intFlux2Upper = ds*(nx*0.5*(xf2ij + xf2ijp1) + ny*0.5*(yf2ij + yf2ijp1))
    intFlux3Upper = ds*(nx*0.5*(xf3ij + xf3ijp1) + ny*0.5*(yf3ij + yf3ijp1))

    # bottom face
    nx = 0.0; ny = -1.0; ds = dx
    intEntFluxBottom = ds*(nx*0.5*(xFij + xFijm1) + ny*0.5*(yFij + yFijm1))
    intFlux1Bottom = ds*(nx*0.5*(xf1ij + xf1ijm1) + ny*0.5*(yf1ij + yf1ijm1))
    intFlux2Bottom = ds*(nx*0.5*(xf2ij + xf2ijm1) + ny*0.5*(yf2ij + yf2ijm1))
    intFlux3Bottom = ds*(nx*0.5*(xf3ij + xf3ijm1) + ny*0.5*(yf3ij + yf3ijm1))

    # integral of the divergence of the flux 
    int_div_flux_term1 = intFlux1Right + intFlux1Upper + intFlux1Left + intFlux1Bottom
    int_div_flux_term2 = intFlux2Right + intFlux2Upper + intFlux2Left + intFlux2Bottom
    int_div_flux_term3 = intFlux3Right + intFlux3Upper + intFlux3Left + intFlux3Bottom

    # integral of the divergence of the entropy flux
    int_div_ent_flux = intEntFluxRight + intEntFluxUpper + intEntFluxLeft + intEntFluxBottom

    sumIntOverPatch(int_div_flux_term1)
    sumIntOverPatch(int_div_flux_term2)
    sumIntOverPatch(int_div_flux_term3)
    sumIntOverPatch(int_div_ent_flux)
    
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
   
    # entropy residual 
    Num = np.abs(int_div_ent_flux - eta_prime_times_int_div_flux)
    aux = np.abs(eta_prime_term1*int_div_flux_term1) + np.abs(eta_prime_term2*int_div_flux_term2) + np.abs(eta_prime_term3*int_div_flux_term3)
    
    
    #Den = np.abs(int_div_ent_flux) + np.abs(eta_prime_times_int_div_flux) + tol
    Den = np.abs(int_div_ent_flux) + aux + tol
    #Den = np.abs(int_div_ent_flux) + abs_eta_prime * abs_int_div_flux + tol
    
    Ri = (Num/Den)
    return Ri
#

get_ent_residual = get_ent_residual_based_on_reconstructions
#get_ent_residual = get_ent_residual_based_on_div_thm
