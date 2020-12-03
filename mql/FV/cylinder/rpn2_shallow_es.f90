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

! =====================================================
subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =====================================================

! Roe-solver for the 2D shallow water equations
! solve Riemann problems along one slice of data.

! waves: 3
! equations: 3

! Conserved quantities:
!       1 depth
!       2 x_momentum
!       3 y_momentum

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell

! This data is along a slice in the x-direction if ixy=1
!                            or the y-direction if ixy=2.

! On output, wave contains the waves, s the speeds,
! and amdq, apdq the decomposition of the flux difference
!   f(qr(i-1)) - f(ql(i))
! into leftgoing and rightgoing parts respectively.
! With the Roe solver we have
!    amdq  =  A^- \Delta q    and    apdq  =  A^+ \Delta q
! where A is the Roe matrix.  An entropy fix can also be incorporated
! into the flux differences.

! Note that the i'th Riemann problem has left state qr(:,i-1)
!                                    and right state ql(:,i)
! From the basic clawpack routines, this routine is called with ql = qr


! This Riemann solver differs from the original clawpack Riemann solver
! for the interleaved indices

  implicit none
  !implicit double precision (a-h,o-z)

  integer, intent(in) :: ixy, maxm, meqn, mwaves, maux, mbc, mx
  double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: ql, qr
  double precision, dimension(maux,1-mbc:maxm+mbc), intent(in) :: auxl, auxr
  double precision, dimension(mwaves,1-mbc:maxm+mbc), intent(out) :: s
  double precision, dimension(meqn,mwaves,1-mbc:maxm+mbc), intent(out) :: wave
  double precision, dimension(meqn,1-mbc:maxm+mbc), intent(out) :: apdq, amdq

  double precision grav, use_dmin_blended
  integer i, m, mw, mu, mv
  double precision a1, a2, a3  
  double precision u(-6:maxm+7),v(-6:maxm+7),a(-6:maxm+7), h(-6:maxm+7)
  double precision hsqrtl, hsqrtr, hsq2  

  double precision dlBlended(mwaves), delta(3), s_roe(mwaves), s_hlle(mwaves)
  double precision wave_roe(meqn, mwaves, 1-mbc:maxm+mbc), wave_hlle(meqn, mwaves, 1-mbc:maxm+mbc)  

  double precision dhL,dhR,x_velL,x_velR,y_velL,y_velR,velL,velR,dmax
  double precision dRimh,s1,s2,h_l,h_r,u_l,u_r,h_m,hu_m,hv_m
  double precision q1rml,q2rml,q3rml,f1rml,f2rml,f3rml

  double precision dmin_roe,dmin_blended,dhuim1,dhui,dhim1,dhi,df2im1,df2i,df1im1,df1i
  double precision dw2im1, dw2i,dw1im1,dw1i
  double precision dEntPotim1,dEntPoti,dEntFluxim1,dEntFluxi,deltaEntPot
  double precision s_times_wave_roe(2), dlBlended_times_wave_roe(2)
  double precision delta_w2,delta_w1,delta_u2,delta_u1
  double precision delta_w_dot_sum_blended_times_wave,delta_w_dot_sum_s_times_wave,delta_w_dot_flux,delta_w_dot_delta_u

  double precision c_hat,c_l,c_r,u_hat,hsqrt_l,hsqrt_r,h_hat
  !   # Gravity constant set in the shallow1D.py or setprob.f90 file
  common /cparam/ grav, use_dmin_blended

!   # Set mu to point to  the component of the system that corresponds
  !   # to momentum in the direction of this slice, mv to the orthogonal
  !   # momentum:
  
  if (ixy == 1) then
     mu = 2
     mv = 3
  else
     mu = 3
     mv = 2
  endif

!   # Note that notation for u and v reflects assumption that the
!   # Riemann problems are in the x-direction with u in the normal
!   # direciton and v in the orthogonal direcion, but with the above
!   # definitions of mu and mv the routine also works with ixy=2
!   # and returns, for example, f0 as the Godunov flux g0 for the
!   # Riemann problems u_t + g(u)_y = 0 in the y-direction.

!   # Compute the Roe-averaged variables needed in the Roe solver.

  do 10 i = 2-mbc, mx+mbc
     h(i) = (qr(1,i-1)+ql(1,i))*0.50d0
     hsqrtl = dsqrt(qr(1,i-1))
     hsqrtr = dsqrt(ql(1,i))
     hsq2 = hsqrtl + hsqrtr
     u(i) = (qr(mu,i-1)/hsqrtl + ql(mu,i)/hsqrtr) / hsq2
     v(i) = (qr(mv,i-1)/hsqrtl + ql(mv,i)/hsqrtr) / hsq2
     a(i) = dsqrt(grav*h(i))       
10 END DO
    
  !   # now split the jump in q at each interface into waves
  do i = 2-mbc, mx+mbc
     ! ****************************************************** !
     ! ********** Compute the Roe waves and speeds ********** !
     ! ****************************************************** !       
     delta(1) = ql(1,i) - qr(1,i-1)
     delta(2) = ql(mu,i) - qr(mu,i-1)
     delta(3) = ql(mv,i) - qr(mv,i-1)
     a1 = ((u(i)+a(i))*delta(1) - delta(2))*(0.50d0/a(i))
     a2 = -v(i)*delta(1) + delta(3)
     a3 = (-(u(i)-a(i))*delta(1) + delta(2))*(0.50d0/a(i))
     
     ! waves and speed
     wave_roe(1,1,i) = a1
     wave_roe(mu,1,i) = a1*(u(i)-a(i))
     wave_roe(mv,1,i) = a1*v(i)
     s_roe(1) = u(i)-a(i)
     
     wave_roe(1,2,i) = 0.0d0
     wave_roe(mu,2,i) = 0.0d0
     wave_roe(mv,2,i) = a2
     s_roe(2) = u(i)
     
     wave_roe(1,3,i) = a3
     wave_roe(mu,3,i) = a3*(u(i)+a(i))
     wave_roe(mv,3,i) = a3*v(i)
     s_roe(3) = u(i)+a(i)
     
     ! ******************************************************************** !
     ! ********** maximum wave speed based on Guermond and Popov ********** !
     ! ******************************************************************** !
     dhL = qr(1,i-1)
     dhR = ql(1,i)
     x_velL = qr(mu,i-1) / qr(1,i-1)
     x_velR = ql(mu,i) / ql(1,i)
     y_velL = qr(mv,i-1) / qr(1,i-1)
     y_velR = ql(mv,i) / ql(1,i)
     
     call get_lambda_llf(grav,dhL,dhR,x_velL,x_velR,dmax)
     s1 = -dmax
     s2 = dmax
          
     ! ***** middle state (taken from HLLE solver) ***** !
     h_l = qr(1,i-1)
     h_r = ql(1,i)
     ! Velocity
     u_l = qr(mu,i-1) / qr(1,i-1)
     u_r = ql(mu,i  ) / ql(1,i  )

     ! speed based on hlle
     !c_l = dsqrt(grav*h_l)
     !c_r = dsqrt(grav*h_r)
     !hsqrt_l = dsqrt(qr(1,i-1))
     !hsqrt_r = dsqrt(ql(1,i))
     !hsq2 = hsqrt_l + hsqrt_r
     !h_hat = 0.5*(h_l + h_r)     
     !u_hat = (u_l*hsqrt_l + u_r*hsqrt_r) / hsq2
     !c_hat = dsqrt(grav*h_hat)     
     ! Speeds of non-shear waves
     !s1 = min(u_l - c_l, u_hat - c_hat)
     !s2 = max(u_r + c_r, u_hat + c_hat)
     !dmax = max(dabs(s1),dabs(s2))
     !s1 = -dmax
     !s2 = dmax
     
     h_m = (ql(mu,i) - qr(mu,i-1) - s2*ql(1,i) + s1*qr(1,i-1))/(s1-s2)
     hu_m = (ql(mu,i)*(u_r-s2) - qr(mu,i-1)*(u_l-s1) + 0.5*grav*(h_r**2 - h_l**2) ) / (s1-s2)
     hv_m = (ql(mv,i)*u_r - qr(mv,i-1)*u_l - s2*ql(mv,i) + s1*qr(mv,i-1))/(s1-s2)
     
     ! ***** left-going wave ***** !
     wave_hlle(1, 1,i) = h_m - h_l
     wave_hlle(mu,1,i) = hu_m - qr(mu,i-1)
     wave_hlle(mv,1,i) = hv_m - qr(mv,i-1)
     
     ! ***** dummy wave ***** !
     wave_hlle(1, 2,i) = 0
     wave_hlle(mu,2,i) = 0
     wave_hlle(mv,2,i) = 0
     
     ! ***** right-going wave ***** !
     wave_hlle(1, 3,i) = h_r - h_m
     wave_hlle(mu,3,i) = ql(mu,i) - hu_m
     wave_hlle(mv,3,i) = ql(mv,i) - hv_m
     
     s_hlle(1) = s1
     s_hlle(2) = 0.
     s_hlle(3) = s2
     
     ! ************************************** !
     ! ********** Waves and speeds ********** !
     ! ************************************** !
     s(1,i) = s_roe(1)
     s(2,i) = s_roe(2)
     s(3,i) = s_roe(3)
     
     wave(1,1,i) = wave_roe(1,1,i)
     wave(2,1,i) = wave_roe(2,1,i)
     wave(3,1,i) = wave_roe(3,1,i)
     
     wave(1,2,i) = wave_roe(1,2,i)
     wave(2,2,i) = wave_roe(2,2,i)
     wave(3,2,i) = wave_roe(3,2,i)
     
     wave(1,3,i) = wave_roe(1,3,i)
     wave(2,3,i) = wave_roe(2,3,i)
     wave(3,3,i) = wave_roe(3,3,i)

     !s(1,i) = s_hlle(1)
     !s(2,i) = s_hlle(2)
     !s(3,i) = s_hlle(3)
     
     !wave(1,1,i) = wave_hlle(1,1,i)
     !wave(2,1,i) = wave_hlle(2,1,i)
     !wave(3,1,i) = wave_hlle(3,1,i)
     
     !wave(1,2,i) = wave_hlle(1,2,i)
     !wave(2,2,i) = wave_hlle(2,2,i)
     !wave(3,2,i) = wave_hlle(3,2,i)
     
     !wave(1,3,i) = wave_hlle(1,3,i)
     !wave(2,3,i) = wave_hlle(2,3,i)
     !wave(3,3,i) = wave_hlle(3,3,i)
     
     ! *************************************** !
     ! ********** Entropy viscosity ********** !
     ! *************************************** !
     ! entropy residual at interface imh
     dRimh = max(auxl(1,i),auxl(1,i-1))
     !dRimh = 0.d0

     ! Blended speeds for stabilization
     dlBlended(1) = dRimh * dmax + (1-dRimh) * dabs(s_roe(1))
     dlBlended(2) = dRimh * dmax + (1-dRimh) * dabs(s_roe(2))
     dlBlended(3) = dRimh * dmax + (1-dRimh) * dabs(s_roe(3))

     ! *************************************** !
     ! ********** Entropy stability ********** !
     ! *************************************** !
     ! compute entropy variables
     dhi = ql(1,i)
     dhui = ql(mu,i)
     dhim1 = qr(1,i-1)
     dhuim1 = qr(mu,i-1)
     call get_ent_variables(grav,dhim1,dhuim1,dw1im1,dw2im1)
     call get_ent_variables(grav,dhi,dhui,dw1i,dw2i)
     
     ! compute fluxes
     df1i = ql(mu,i)
     df1im1 = qr(mu,i-1)       
     df2i = ql(mu,i)*(u_r) + 0.5*grav*h_r**2
     df2im1 = qr(mu,i-1)*(u_l) + 0.5*grav*h_l**2
     call get_ent_flux(grav,dhim1,dhuim1,dEntFluxim1)
     call get_ent_flux(grav,dhi,dhui,dEntFluxi)
     
     ! compute entropy potential
     call get_ent_pot(dw1im1,dw2im1,df1im1,df2im1,dEntFluxim1,dEntPotim1)
     call get_ent_pot(dw1i,dw2i,df1i,df2i,dEntFluxi,dEntPoti)
     
     ! compute some auxiliary variables
     delta_w1 = dw1i-dw1im1
     delta_w2 = dw2i-dw2im1
     delta_u1 = dhi-dhim1
     delta_u2 = dhui-dhuim1       
     deltaEntPot = dEntPoti-dEntPotim1
     delta_w_dot_flux = delta_w1*.5d0*(df1i+df1im1) + delta_w2*.5d0*(df2i+df2im1)
     delta_w_dot_delta_u = delta_w1*delta_u1 + delta_w2*delta_u2

     
     ! compute dmin_roe operator
     s_times_wave_roe(1) = dabs(s_roe(1))*wave_roe(1, 1,i) + dabs(s_roe(3))*wave_roe(1, 3,i)
     s_times_wave_roe(2) = dabs(s_roe(1))*wave_roe(mu,1,i) + dabs(s_roe(3))*wave_roe(mu,3,i)
     delta_w_dot_sum_s_times_wave = delta_w1*s_times_wave_roe(1) + delta_w2*s_times_wave_roe(2)
     dmin_roe=max((delta_w_dot_flux-deltaEntPot-.5*delta_w_dot_sum_s_times_wave)/(.5d0*delta_w_dot_delta_u+1E-7),0.d0)
     dmin_roe=min(dmin_roe,dmax) ! this is for safety
       
     ! compute dmin_blended operator
     dlBlended_times_wave_roe(1) = dlBlended(1)*wave_roe(1, 1,i) + dlBlended(3)*wave_roe(1, 3,i)
     dlBlended_times_wave_roe(2) = dlBlended(1)*wave_roe(mu,1,i) + dlBlended(3)*wave_roe(mu,3,i)
     delta_w_dot_sum_blended_times_wave = delta_w1*dlBlended_times_wave_roe(1) + delta_w2*dlBlended_times_wave_roe(2)
     dmin_blended=max((delta_w_dot_flux-deltaEntPot-.5d0*delta_w_dot_sum_blended_times_wave)/(.5d0*delta_w_dot_delta_u+1E-7),0.d0)
     dmin_blended=min(dmin_blended,dmax)
     dmin_blended = use_dmin_blended * dmin_blended

     !if (dmin_blended.gt.1E-20) then
     !   print *, dmin_blended
     !end if
     ! ****************************************** !
     ! ********** Compute fluctuations ********** !
     ! ****************************************** !
     do m=1,3
        amdq(m,i) = 0
        apdq(m,i) = 0
        do mw=1,mwaves
           ! Purely Roe's
           !amdq(m,i) = amdq(m,i) + .5d0*(s_roe(mw) - abs(s_roe(mw)))*wave_roe(m,mw,i)
           !apdq(m,i) = apdq(m,i) + .5d0*(s_roe(mw) + abs(s_roe(mw)))*wave_roe(m,mw,i)

           ! Blended
           !amdq(m,i) = amdq(m,i) + .5d0*(s_roe(mw) - dlBlended(mw))*wave_roe(m,mw,i)
           !apdq(m,i) = apdq(m,i) + .5d0*(s_roe(mw) + dlBlended(mw))*wave_roe(m,mw,i)
           
           ! Blended with dmin_blended
           amdq(m,i) = amdq(m,i) + .5d0*(s_roe(mw) - dlBlended(mw) - dmin_blended)*wave_roe(m,mw,i)
           apdq(m,i) = apdq(m,i) + .5d0*(s_roe(mw) + dlBlended(mw) + dmin_blended)*wave_roe(m,mw,i)

           !amdq(m,i) = amdq(m,i) + .5d0*(s_hlle(mw) - dabs(s_hlle(mw)))*wave_hlle(m,mw,i)
           !apdq(m,i) = apdq(m,i) + .5d0*(s_hlle(mw) + dabs(s_hlle(mw)))*wave_hlle(m,mw,i)

           ! Blended with dmin_roe
           !amdq(m,i) = amdq(m,i) + .5d0*(s_roe(mw) - dlBlended(mw) - dmin_roe)*wave_roe(m,mw,i)
           !apdq(m,i) = apdq(m,i) + .5d0*(s_roe(mw) + dlBlended(mw) + dmin_roe)*wave_roe(m,mw,i)
        end do
     end do

     if (s_roe(1) .le. 0) then 
        s(1,i) = - (dlBlended(1) + dmin_blended)
     else
        s(1,i) =   (dlBlended(1) + dmin_blended)
     end if
     
     if (s_roe(2) .le. 0) then 
        s(2,i) = - (dlBlended(2) + dmin_blended)
     else
        s(2,i) =   (dlBlended(2) + dmin_blended)
     end if
     
     if (s_roe(3) .le. 0) then 
        s(3,i) = - (dlBlended(3) + dmin_blended)
     else
        s(3,i) =   (dlBlended(3) + dmin_blended)
     end if
     
  END DO
  return
end subroutine rpn2
