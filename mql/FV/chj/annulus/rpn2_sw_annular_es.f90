subroutine get_ent_variables(grav,dh,dhu,dw1,dw2)
  double precision grav,dh,dhu,dw1,dw2,du
  du = dhu/dh
  dw1 = grav*dh-.5d0*du**2.d0
  dw2 = du
end subroutine get_ent_variables

subroutine get_ent_flux(grav,dh,dhu,dEntFlux)
  double precision grav,dh,dhu,dEntFlux,du
  du = dhu/dh
  dEntFlux = (.5d0*dh*du*du+grav*dh*dh)*du
end subroutine get_ent_flux

subroutine get_ent_pot(dw1,dw2,df1,df2,dEntFlux,dEntPot)
  double precision dw1,dw2,df1,df2,dEntFlux,dEntPot
  dEntPot = dw1*df1+dw2*df2-dEntFlux
end subroutine get_ent_pot

subroutine get_fZ(grav,dh,dhZ,fZ)
  double precision grav,dh,dhZ,fZ
  if (dh<= dhZ) then
     fZ = 2.d0*(dsqrt(grav*dh)-dsqrt(grav*dhZ))
  else
     fZ = (dh-dhZ)*dsqrt(0.5d0*grav*(dh+dhZ)/dh/dhZ)
  endif
  return
end subroutine get_fZ

subroutine get_f(grav,dh,dhL,dhR,velL,velR,f)
  double precision grav,dh,dhL,dhR,velL,velR,dfL,dfR,f
  call get_fZ(grav,dh,dhL,dfL)
  call get_fZ(grav,dh,dhR,dfR) 
  f = dfL + dfR + velR - velL
  return
end subroutine get_f

subroutine get_lambda_llf(grav,dhL,dhR,velL,velR,dmax)
  double precision grav,dhL,dhR,velL,velR,dmax
  double precision dhMin,dhMax,x0,dfMin,dfMax,dhStar,dlambda1,dlambda3
  dhMin = min(dhL,dhR)
  dhMax = max(dhL,dhR)            
  x0 = (2.d0*dsqrt(2.d0)-1.d0)**2.d0
  call get_f(grav,x0*dhMin,dhL,dhR,velL,velR,dfMin)
  call get_f(grav,x0*dhMax,dhL,dhR,velL,velR,dfMax)
  if (0 <= dfMin) then
     dhStar = (max(0.d0,velL-velR+2.d0*dsqrt(grav)*(dsqrt(dhL)+dsqrt(dhR))))**2.d0/(16.d0*grav)
  else
     if (0 <= dfMax) then
        dhStar = (-dsqrt(2.d0*dhMin)+dsqrt(3.d0*dhMin+2.d0*dsqrt(2.d0*dhMin*dhMax)+dsqrt(2.d0/grav)*(velL-velR)*dsqrt(dhMin)))**2.d0
     else ! fMax < 0
        dhStar = dsqrt(dhMin*dhMax)*(1.d0+(dsqrt(2.d0)*(velL-velR))/(dsqrt(grav*dhMin)+dsqrt(grav*dhMax)))
     endif
  endif
  ! Compute max wave speed based on hStar0
  dlambda1 = velL - dsqrt(grav*dhL)*dsqrt( (1.d0+max((dhStar-dhL)/(2.d0*dhL),0.d0)) * (1.d0+max((dhStar-dhL)/dhL,0.d0)))
  dlambda3 = velR + dsqrt(grav*dhR)*dsqrt( (1.d0+max((dhStar-dhR)/(2.d0*dhR),0.d0)) * (1.d0+max((dhStar-dhR)/dhR,0.d0)))
  dmax = max(dabs(dlambda1),dabs(dlambda3))
end subroutine get_lambda_llf


subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)

! Roe-solver for the 2D shallow water equations on a mapped 
! (logically quadrilateral) grid

! waves: 3
! equations: 3

! Conserved quantities:
!       1 depth
!       2 x-momentum
!       3 y-momentum

    implicit none

    integer, intent(in) :: ixy, maxm, meqn, mwaves, maux, mbc, mx
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: ql, qr
    double precision, dimension(maux,1-mbc:maxm+mbc), intent(in) :: auxl, auxr
    double precision, dimension(meqn,mwaves,1-mbc:maxm+mbc), intent(out) :: wave
    double precision, dimension(mwaves,1-mbc:maxm+mbc), intent(out) :: s
    double precision, dimension(meqn,1-mbc:maxm+mbc), intent(out) :: amdq, apdq

    ! local arrays -- common block comroe is passed to rpt2sh
    ! ------------

    integer :: depth
    integer :: i, m, mw
    integer :: inx, iny, ilenrat
    double precision :: a1, a2, a3

    integer, parameter :: maxm2 = 1602  ! assumes at most 1000x1000 grid with mbc=2
    double precision, dimension(3) :: delta, dlBlended
    double precision, dimension(-1:maxm2) :: u, v, a, h, hunor_r, hunor_l, hutan_r, hutan_l, alpha, beta
    double precision :: hsqrt_l, hsqrt_r, hsq2, dmax
    double precision :: grav, use_dmin_blended
    double precision wave_hlle(meqn,   mwaves, 1-mbc:maxm+mbc)
    double precision w11,w12,w13,w21,w22,w23,w31,w32,w33,s1,s2,Q1m,Q2m,Q3m
    double precision Rimh,c_l,c_r,h_hat,u_hat,c_hat

    double precision hi,hui,him1,huim1,ui,uim1,w1im1,w2im1,w1i,w2i
    double precision EntFluxi,EntFluxim1
    
    double precision dmin_roe,f1i,f1im1,f2i,f2im1,EntPoti,EntPotim1
    double precision delta_w1,delta_w2,delta_u1,delta_u2,deltaEntPot,delta_w_dot_flux,delta_w_dot_delta_u
    double precision delta_w_dot_sum_s_times_wave

    double precision s_times_wave_roe(2), dlBlended_times_wave_roe(2)
    double precision delta_w_dot_sum_blended_times_wave,dmin_blended
    common /cparam/ grav, use_dmin_blended
    common /comroe/ u, v, a, h

    if (-1.gt.1-mbc .or. maxm2 .lt. maxm+mbc) then
        write(6,*) 'need to increase maxm2 in rpA'
        stop
    endif

    ! Rotate velocity components to be aligned with grid normal.
    ! The normal vector for the face at the i'th Riemann problem
    ! is stored in aux(1,...) and aux(2,...) if ixy = 1
    ! and in aux(4,...) and aux(5,...) if ixy = 2.
    ! The ratio of the length of the cell side to the area of
    ! the cell is stored in aux(3,...) or aux(6,...), respectively.

    depth = 1
    if (ixy.eq.1) then
        inx = 1
        iny = 2
        ilenrat = 3
    else
        inx = 4
        iny = 5
        ilenrat = 6
    endif

    ! Compute the rotation matrix:
    ! [ alpha beta ]
    ! [-beta alpha ]
    ! Then compute the normal and tangential momentum components to the right and left.
    do i=2-mbc, mx+mbc
        alpha(i) = auxl(inx,i)
        beta(i) = auxl(iny,i)
        hunor_r(i)   =  alpha(i) * ql(2,i)   + beta(i)  * ql(3,i)
        hunor_l(i-1) =  alpha(i) * qr(2,i-1) + beta(i)  * qr(3,i-1)
        hutan_r(i)   = -beta(i)  * ql(2,i)   + alpha(i) * ql(3,i)
        hutan_l(i-1) = -beta(i)  * qr(2,i-1) + alpha(i) * qr(3,i-1)
    enddo

    ! compute the Roe-averaged variables needed in the Roe solver.
    ! These are stored in the common block comroe since they are
    ! later used in routine rpt2 to do the transverse wave splitting.
    do i = 2-mbc, mx+mbc
        hsqrt_l = dsqrt(qr(depth,i-1))
        hsqrt_r = dsqrt(ql(depth,i))
        hsq2 = hsqrt_l + hsqrt_r
        
        ! Roe averages
        h(i) = (qr(depth,i-1)+ql(depth,i))*0.50d0 ! = h_hat
        u(i) = (hunor_l(i-1)/hsqrt_l + hunor_r(i)/hsqrt_r) / hsq2
        v(i) = (hutan_l(i-1)/hsqrt_l + hutan_r(i)/hsqrt_r) / hsq2
        a(i) = dsqrt(grav*h(i)) ! = c_hat
    enddo

    ! main loop to compute the waves and speeds
    do i = 2-mbc, mx+mbc
       ! left and right states in the normal direction
       hi = ql(1,i)
       hui = hunor_r(i)
       him1 = qr(1,i-1)
       huim1 = hunor_l(i-1)
       uim1 = hunor_l(i-1) / him1 
       ui = hunor_r(i) / hi
       
       ! *************************************** !
       ! ********** Entropy stability ********** !
       ! *************************************** !
       ! compute entropy variables
       call get_ent_variables(grav,him1,huim1,w1im1,w2im1)
       call get_ent_variables(grav,hi,hui,w1i,w2i)
       
       ! compute fluxes in the normal direction
       f1i = hunor_r(i)
       f2i = hunor_r(i)*ui + 0.5*grav*hi**2
       f1im1 = hunor_l(i-1)
       f2im1 = hunor_l(i-1)*uim1 + 0.5*grav*him1**2

       ! compute entropy fluxes in the normal direction
       call get_ent_flux(grav,him1,huim1,EntFluxim1)
       call get_ent_flux(grav,hi,hui,EntFluxi)

       ! compute entropy potential in the normal direction
       call get_ent_pot(w1im1,w2im1,f1im1,f2im1,EntFluxim1,EntPotim1)
       call get_ent_pot(w1i,w2i,f1i,f2i,EntFluxi,EntPoti)

       ! compute some auxiliary variables
       delta_w1 = w1i-w1im1
       delta_w2 = w2i-w2im1
       delta_u1 = hi-him1
       delta_u2 = hui-huim1       
       deltaEntPot = EntPoti-EntPotim1
       delta_w_dot_flux = delta_w1*.5d0*(f1i+f1im1) + delta_w2*.5d0*(f2i+f2im1)
       delta_w_dot_delta_u = delta_w1*delta_u1 + delta_w2*delta_u2
       
       ! ****************************************************** !
       ! ********** Compute the Roe waves and speeds ********** !
       ! ****************************************************** !
       delta(1) = ql(depth,i) - qr(depth,i-1)
       delta(2) = hunor_r(i) - hunor_l(i-1)
       delta(3) = hutan_r(i) - hutan_l(i-1)
       a1 = ((u(i)+a(i))*delta(1) - delta(2))*(0.50d0/a(i))
       a2 = -v(i)*delta(1) + delta(3)
       a3 = (-(u(i)-a(i))*delta(1) + delta(2))*(0.50d0/a(i))

       ! Compute the waves and speed
       w11 = a1
       w12 = a1*(u(i)-a(i))
       w13 = a1*v(i)

       w21 = 0.d0
       w22 = 0.d0
       w23 = a2

       w31 = a3
       w32 = a3*(u(i)+a(i)) 
       w33 = a3*v(i)

       s(1,i) = (u(i)-a(i))
       s(2,i) = u(i) 
       s(3,i) = (u(i)+a(i)) 
     
       ! Rotate back the Roe waves 
       wave(1,1,i) = w11
       wave(2,1,i) = alpha(i)*w12 - beta(i)*w13
       wave(3,1,i) = beta(i)*w12 + alpha(i)*w13
              
       wave(1,2,i) = w21 
       wave(2,2,i) = alpha(i)*w22 - beta(i)*w23 
       wave(3,2,i) = beta(i)*w22 + alpha(i)*w23 
       
       wave(1,3,i) = w31 
       wave(2,3,i) = alpha(i)*w32 - beta(i)*w33 
       wave(3,3,i) = beta(i)*w32 + alpha(i)*w33

       ! scale speeds and viscosity due to rotation
       s(1,i) = s(1,i) * auxl(ilenrat,i)
       s(2,i) = s(2,i) * auxl(ilenrat,i)
       s(3,i) = s(3,i) * auxl(ilenrat,i)
       
       ! ***************************************************************** !
       ! ********** Compute Rusanov's (or LLF) waves and speeds ********** !
       ! ***************************************************************** !
       call get_lambda_llf(grav,him1,hi,uim1,ui,dmax)
       s1 = -dmax
       s2 = dmax

       ! compute middle state based on rotated velocity fields
       Q1m = (hunor_r(i) - hunor_l(i-1) - s2*hi + s1*him1)/(s1-s2)
       Q2m = (hunor_r(i)*(ui-s2) - hunor_l(i-1)*(uim1-s1) + 0.5*grav*(hi**2-him1**2))/(s1-s2)
       Q3m = (hutan_r(i)*ui - hutan_l(i-1)*uim1 - s2*hutan_r(i) + s1*hutan_l(i-1))/(s1-s2)

       ! compute the 1- and 3-waves for the hlle solver
       w11 = Q1m - him1
       w12 = Q2m - hunor_l(i-1)
       w13 = Q3m - hutan_l(i-1)

       w31 = hi - Q1m
       w32 = hunor_r(i) - Q2m 
       w33 = hutan_r(i) - Q3m

       ! rotate back the 1- and 3-waves
       wave_hlle(1, 1, i) = w11
       wave_hlle(2, 1, i) = alpha(i) * w12 - beta(i) * w13
       wave_hlle(3, 1, i) = beta(i) * w12 + alpha(i) * w13

       wave_hlle(1, 2, i) = 0.d0
       wave_hlle(2, 2, i) = 0.d0
       wave_hlle(3, 2, i) = 0.d0
       
       wave_hlle(1, 3, i) = w31
       wave_hlle(2, 3, i) = alpha(i) * w32 - beta(i) * w33
       wave_hlle(3, 3, i) = beta(i) * w32 + alpha(i) * w33

       ! rescale the wave speeds
       s1 = s1 * auxl(ilenrat,i)
       s2 = s2 * auxl(ilenrat,i)

       ! ********************************************** !
       ! ********** Compute waves and speeds ********** !
       ! ********************************************** !
       !wave(1, 1, i) = wave_hlle(1, 1, i)
       !wave(2, 1, i) = wave_hlle(2, 1, i)
       !wave(3, 1, i) = wave_hlle(3, 1, i)

       !wave(1, 2, i) = 0.d0
       !wave(2, 2, i) = 0.d0
       !wave(3, 2, i) = 0.d0
       
       !wave(1, 3, i) = wave_hlle(1, 3, i)
       !wave(2, 3, i) = wave_hlle(2, 3, i)
       !wave(3, 3, i) = wave_hlle(3, 3, i)

       !s(1,i) = s1
       !s(2,i) = 0
       !s(3,i) = s2

       ! *************************************** !
       ! ********** Entropy viscosity ********** !
       ! *************************************** !
       ! entropy residual at interface imh
       Rimh = max(auxl(8,i),auxl(8,i-1))
       dmax = dmax * auxl(ilenrat,i)
       
       ! Blended speeds for stabilization
       dlBlended(1) = Rimh * dmax + (1-Rimh) * dabs(s(1,i))
       dlBlended(2) = Rimh * dmax + (1-Rimh) * dabs(s(2,i))
       dlBlended(3) = Rimh * dmax + (1-Rimh) * dabs(s(3,i))

       ! compute dmin_roe operator
       s_times_wave_roe(1) = dabs(s(1,i))*w11  + dabs(s(3,i))*w31 
       s_times_wave_roe(2) = dabs(s(1,i))*w12 + dabs(s(3,i))*w32
       delta_w_dot_sum_s_times_wave = delta_w1*s_times_wave_roe(1) + delta_w2*s_times_wave_roe(2)
       dmin_roe = max((delta_w_dot_flux - deltaEntPot - 0.5*delta_w_dot_sum_s_times_wave)/(delta_w_dot_delta_u + 1E-7),0.d0)
       dmin_roe = min(dmin_roe,dmax) ! this is for safety
       dmin_roe = 0.d0
       dmin_roe = dmin_roe * auxl(ilenrat,i)
       
       ! compute dmin_blended operator
       dlBlended_times_wave_roe(1) = dlBlended(1)*w11 + dlBlended(3)*w31
       dlBlended_times_wave_roe(2) = dlBlended(1)*w12 + dlBlended(3)*w32
       delta_w_dot_sum_blended_times_wave = delta_w1*dlBlended_times_wave_roe(1) + delta_w2*dlBlended_times_wave_roe(2)
       dmin_blended=max((delta_w_dot_flux-deltaEntPot-.5d0*delta_w_dot_sum_blended_times_wave)/(.5d0*delta_w_dot_delta_u+1E-7),0.d0)
       dmin_blended=min(dmin_blended,dmax)
       dmin_blended = use_dmin_blended * dmin_blended
       dmin_blended = dmin_blended * auxl(ilenrat,i)

       do m=1,3
          amdq(m,i) = 0.d0
          apdq(m,i) = 0.d0
          do mw=1,mwaves
             amdq(m,i) = amdq(m,i) + .5d0*( s(mw,i) - dlBlended(mw) - dmin_blended)*wave(m,mw,i)
             apdq(m,i) = apdq(m,i) + .5d0*( s(mw,i) + dlBlended(mw) + dmin_blended)*wave(m,mw,i)
          enddo
       enddo

       if (s(1,i) .le. 0) then 
          s(1,i) = - (dlBlended(1) + dmin_blended)
       else
          s(1,i) =   (dlBlended(1) + dmin_blended)
       end if

       if (s(2,i) .le. 0) then 
          s(2,i) = - (dlBlended(2) + dmin_blended)
       else
          s(2,i) =   (dlBlended(2) + dmin_blended)
       end if

       if (s(3,i) .le. 0) then 
          s(3,i) = - (dlBlended(3) + dmin_blended)
       else
          s(3,i) =   (dlBlended(3) + dmin_blended)
       end if

    enddo
    
end subroutine rpn2
