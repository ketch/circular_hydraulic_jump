subroutine get_f(grav,dh,dhZ,f)
  double precision grav,dh,dhZ,f
  if (dh<= dhZ) then
     f = 2.0d0*(dsqrt(grav*dh)-dsqrt(grav*dhZ))
  else
     f = (dh-dhZ)*sqrt(0.5d0*grav*(dh+dhZ)/dh/dhZ)
  endif
  return
end subroutine get_f

subroutine get_phi(grav,dh,dhL,dhR,velL,velR,phi)
  double precision grav,dh,dhL,dhR,velL,velR,dfL,dfR,phi
  call get_f(grav,dh,dhL,dfL)
  call get_f(grav,dh,dhR,dfR) 
  phi = dfL + dfR + velR - velL
  return
end subroutine get_phi

subroutine get_lambda_llf(grav,dhL,dhR,velL,velR,dlambda_llf)
  double precision grav, dhL, dhR, velL, velR, dlambda_llf
  double precision dhMin, dhMax, x0, dlambda1, dlambda3
  double precision dfMin, dfMax
  dhMin = min(dhL,dhR)
  dhMax = max(dhL,dhR)            
  x0 = (2.0d0*dsqrt(2.0d0)-1.0d0)**2.0d0
  call get_phi(grav,x0*dhMin,dhL,dhR,velL,velR,dfMin)
  call get_phi(grav,x0*dhMax,dhL,dhR,velL,velR,dfMax)
  if (0 <= dfMin) then
     dhStar = (max(0.0d0,velL-velR+2.0d0*sqrt(grav)*(sqrt(hL)+sqrt(hR))))**2.0d0/16.0d0/grav
  else
     if (0 <= dfMax) then
        dhStar = (-sqrt(2.*dhMin)+sqrt(3.*dhMin+2.*sqrt(2.*dhMin*dhMax)+sqrt(2./grav)*(velL-velR)*sqrt(dhMin)))**2.d0
     else ! fMax < 0
        dhStar = sqrt(dhMin*dhMax)*(1.+(sqrt(2.d0)*(velL-velR))/(sqrt(grav*dhMin)+sqrt(grav*dhMax)))
     endif
  endif
  ! Compute max wave speed based on hStar0
  dlambda1 = ( velL - sqrt(grav*dhL)*sqrt( (1.+max((dhStar-dhL)/2.d0/dhL,0.)) * (1.+max((dhStar-dhL)/dhL,0.))))
  dlambda3 = ( velR + sqrt(grav*dhR)*sqrt( (1.+max((dhStar-dhR)/2.d0/dhR,0.)) * (1.+max((dhStar-dhR)/dhR,0.))))  
  dlambda_llf = max(abs(dlambda1),abs(dlambda2))
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
  !  implicit double precision (a-h,o-z)

  double precision ql, qr, auxl, auxr, wave, s, amdq, apdq
  integer, intent(in) :: ixy, maxm, meqn, mwaves, maux, mbc, mx
  double precision grav
  double precision a1, a2, a3
  double precision delta(3)
  double precision dlev, dlllf, drimh, hsq2, hsqrtl, hsqrtr
  integer i, m, mu, mv, mw, maxm2

  dimension   ql(meqn,           1-mbc:maxm+mbc)
    dimension   qr(meqn,           1-mbc:maxm+mbc)
    dimension   auxl(maux,           1-mbc:maxm+mbc)
    dimension   auxr(maux,           1-mbc:maxm+mbc)    
    dimension   s(mwaves,         1-mbc:maxm+mbc)
    dimension   wave(meqn,   mwaves, 1-mbc:maxm+mbc)
    dimension   apdq(meqn,          1-mbc:maxm+mbc)
    dimension   amdq(meqn,          1-mbc:maxm+mbc)

!   # Gravity constant set in the shallow1D.py or setprob.f90 file
    common /cparam/ grav

!   # Roe averages quantities of each interface
    parameter (maxm2 = 1800)
    double precision u(-6:maxm2+7),v(-6:maxm2+7),a(-6:maxm2+7), h(-6:maxm2+7)
    double precision dlambda(mwaves), dlRoe(mwaves), dlBlended(mwaves)
    double precision wave_roe(meqn, mwaves, 1-mbc:maxm+mbc), wave_llf(meqn, 2, 1-mbc:maxm+mbc)
    double precision s_roe(mwaves, 1-mbc:maxm+mbc), s_llf(2, 1-mbc:maxm+mbc)
    
!   local arrays
!   ------------
    !dimension delta(3)
    logical :: efix

    data efix /.false./    !# Use entropy fix for transonic rarefactions

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
!   # find a1 thru a3, the coefficients of the 3 eigenvectors:
    do 20 i = 2-mbc, mx+mbc
       delta(1) = ql(1,i) - qr(1,i-1)
       delta(2) = ql(mu,i) - qr(mu,i-1)
       delta(3) = ql(mv,i) - qr(mv,i-1)
       a1 = ((u(i)+a(i))*delta(1) - delta(2))*(0.50d0/a(i))
       a2 = -v(i)*delta(1) + delta(3)
       a3 = (-(u(i)-a(i))*delta(1) + delta(2))*(0.50d0/a(i))

       ! ****************************************************** !
       ! ********** Compute the Roe waves and speeds ********** !
       ! ****************************************************** !
       wave(1,1,i) = a1
       wave(mu,1,i) = a1*(u(i)-a(i))
       wave(mv,1,i) = a1*v(i)
       s(1,i) = u(i)-a(i)
       
       wave(1,2,i) = 0.0d0
       wave(mu,2,i) = 0.0d0
       wave(mv,2,i) = a2
       s(2,i) = u(i)
       
       wave(1,3,i) = a3
       wave(mu,3,i) = a3*(u(i)+a(i))
       wave(mv,3,i) = a3*v(i)
       s(3,i) = u(i)+a(i)

    20 END DO


!    # compute flux differences amdq and apdq.
!    ---------------------------------------
    !if (.false.) go to 110
!     # based on entropy residual
!     ----------------
!     # amdq = SUM s*wave   over left-going waves
!     # apdq = SUM s*wave   over right-going waves

!    do m=1,3
!       do i=2-mbc, mx+mbc
!          amdq(m,i) = 0.d0
!          apdq(m,i) = 0.d0
!
!          ! contribution of Roe
!          do mw=1,mwaves
!             !amdq(m,i) = amdq(m,i) + 0.5d0*(s_roe(mw,i) - dlRoe(mw)) * wave_roe(m,mw,i)
!             !apdq(m,i) = apdq(m,i) + 0.5d0*(s_roe(mw,i) + dlRoe(mw)) * wave_roe(m,mw,i)
!
!             !amdq(m,i) = amdq(m,i) + 0.5d0*(s_roe(mw,i) - dlRoe(mw)) * wave_roe(m,mw,i)
!             !apdq(m,i) = apdq(m,i) + 0.5d0*(s_roe(mw,i) + dlRoe(mw)) * wave_roe(m,mw,i)
             
             
!               !if (s(mw,i) < 0.d0) then
!               !amdq(m,i) = amdq(m,i) + 2*0.5*MIN(0.d0,s(mw,i)-1.0d0*dlambda(mw))*wave(m,mw,i)
!               !amdq(m,i) = amdq(m,i) + 0.5*MIN(0.d0,s(mw,i)-1.0d0*dlmax)*wave(m,mw,i)
!               !else
!               !apdq(m,i) = apdq(m,i) + 2*0.5*MAX(0.d0,s(mw,i)+1.0d0*dlambda(mw))*wave(m,mw,i)
!               !apdq(m,i) = apdq(m,i) + 0.5*MAX(0.d0,s(mw,i)+1.0d0*dlmax)*wave(m,mw,i)
!               !endif
               
!             amdq(m,i) = amdq(m,i) + MIN(0.d0,s(mw,i))*wave(m,mw,i)
!             apdq(m,i) = apdq(m,i) + MAX(0.d0,s(mw,i))*wave(m,mw,i)

               !if (m .NE. 2) then 
               !amdq(m,i) = amdq(m,i) + 0.5d0*(s(mw,i) - dlLLF) * wave(m,mw,i)
               !apdq(m,i) = apdq(m,i) + 0.5d0*(s(mw,i) + dlLLF) * wave(m,mw,i)

               !amdq(m,i) = amdq(m,i) + 0.5d0*(s(mw,i) - dlRoe(mw) - dlEV) *wave(m,mw,i)
               !apdq(m,i) = apdq(m,i) + 0.5d0*(s(mw,i) + dlRoe(mw) + dlEV) *wave(m,mw,i)

               !amdq(m,i) = amdq(m,i) + 0.5d0*(s(mw,i) - dlBlended(mw)) *wave(m,mw,i)
               !apdq(m,i) = apdq(m,i) + 0.5d0*(s(mw,i) + dlBlended(mw)) *wave(m,mw,i)               
               !endif
!          END DO
!       END DO
!    END DO
  !       go to 900


      
!     # based on entropy residual
!     ----------------
!     # amdq = SUM s*wave   over left-going waves
!     # apdq = SUM s*wave   over right-going waves

    do 100 m=1,3
        do 100 i=2-mbc, mx+mbc
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0

            !dlambda(1) = max(abs(dlambda1im1), abs(dlambda1i))
            !dlambda(2) = max(abs(dlambda2im1), abs(dlambda2i))
            !dlambda(3) = max(abs(dlambda3im1), abs(dlambda3i))

            dlLLF = max(abs(qr(mu,i-1)/qr(1,i-1) - dsqrt(grav*qr(1,i-1))), abs(ql(mu,i)/ql(1,i)     + dsqrt(grav*ql(1,i))))
            dRimh = max(auxl(1,i),auxl(1,i-1))
            dlEV = dlLLF * dRimh

            dlRoe(1) = abs(s(1,i))
            dlRoe(2) = abs(s(2,i))
            dlRoe(3) = abs(s(3,i))

            dlBlended(1) = dRimh * dlLLF + (1-dRimh) * dlRoe(1)
            dlBlended(2) = dRimh * dlLLF + (1-dRimh) * dlRoe(2)
            dlBlended(3) = dRimh * dlLLF + (1-dRimh) * dlRoe(3)

            !dlLLF = max(abs(s(1,i)),abs(s(3,i)))
            do 90 mw=1,mwaves
               !if (s(mw,i) < 0.d0) then
               !amdq(m,i) = amdq(m,i) + 2*0.5*MIN(0.d0,s(mw,i)-1.0d0*dlambda(mw))*wave(m,mw,i)
               !amdq(m,i) = amdq(m,i) + 0.5*MIN(0.d0,s(mw,i)-1.0d0*dlmax)*wave(m,mw,i)
               !else
               !apdq(m,i) = apdq(m,i) + 2*0.5*MAX(0.d0,s(mw,i)+1.0d0*dlambda(mw))*wave(m,mw,i)
               !apdq(m,i) = apdq(m,i) + 0.5*MAX(0.d0,s(mw,i)+1.0d0*dlmax)*wave(m,mw,i)
               !endif

               
               !amdq(m,i) = amdq(m,i) + MIN(0.d0,s(mw,i))*wave(m,mw,i)
               !apdq(m,i) = apdq(m,i) + MAX(0.d0,s(mw,i))*wave(m,mw,i)

               !if (m .NE. 2) then 
               !amdq(m,i) = amdq(m,i) + 0.5d0*(s(mw,i) - dlLLF) * wave(m,mw,i)
               !apdq(m,i) = apdq(m,i) + 0.5d0*(s(mw,i) + dlLLF) * wave(m,mw,i)

!               amdq(m,i) = amdq(m,i) + 0.5d0*(s(mw,i) - dlRoe(mw) - dlEV) *wave(m,mw,i)
!               apdq(m,i) = apdq(m,i) + 0.5d0*(s(mw,i) + dlRoe(mw) + dlEV) *wave(m,mw,i)
               amdq(m,i) = amdq(m,i) + 0.5d0*(s(mw,i) - dlBlended(mw)) *wave(m,mw,i)
               apdq(m,i) = apdq(m,i) + 0.5d0*(s(mw,i) + dlBlended(mw)) *wave(m,mw,i)               
               !endif
            90 END DO
    100 END DO
!    go to 900

!-----------------------------------------------------

!    110 continue

!     # With entropy fix
!     ------------------

!    # compute flux differences amdq and apdq.
!    # First compute amdq as sum of s*wave for left going waves.
!    # Incorporate entropy fix by adding a modified fraction of wave
!    # if s should change sign.

!    do 200 i=2-mbc,mx+mbc
!       # check 1-wave
!        him1 = qr(1,i-1)
!        s0 =  qr(mu,i-1)/him1 - dsqrt(grav*him1)
!!       # check for fully supersonic case :
!        if (s0 > 0.0d0 .AND. s(1,i) > 0.0d0) then
!            do 60 m=1,3
!                amdq(m,i)=0.0d0
!            60 END DO
!            goto 200
!        endif
    
!        h1 = qr(1,i-1)+wave(1,1,i)
!        hu1= qr(mu,i-1)+wave(mu,1,i)
!        s1 = hu1/h1 - dsqrt(grav*h1) !speed just to right of 1-wave
!        if (s0 < 0.0d0 .AND. s1 > 0.0d0) then
!!           # transonic rarefaction in 1-wave
!            sfract = s0*((s1-s(1,i))/(s1-s0))
!        else if (s(1,i) < 0.0d0) then
!!           # 1-wave is leftgoing
!            sfract = s(1,i)
!        else
!!           # 1-wave is rightgoing
!            sfract = 0.0d0
!        endif
!        do 120 m=1,3
!            amdq(m,i) = sfract*wave(m,1,i)
!        120 END DO
!
!!       # check 2-wave
!        if (s(2,i) > 0.0d0) then
!!           # 2 and 3 waves are right-going
!            go to 200
!        endif

!        do 140 m=1,3
!            amdq(m,i) = amdq(m,i) + s(2,i)*wave(m,2,i)
!        140 END DO
    
!!       check 3-wave
    
!        hi = ql(1,i)
!        s03 = ql(mu,i)/hi + dsqrt(grav*hi)
!        h3=ql(1,i)-wave(1,3,i)
!        hu3=ql(mu,i)-wave(mu,3,i)
!        s3=hu3/h3 + dsqrt(grav*h3)
!        if (s3 < 0.0d0 .AND. s03 > 0.0d0) then
!!           # transonic rarefaction in 3-wave
!            sfract = s3*((s03-s(3,i))/(s03-s3))
!        else if (s(3,i) < 0.0d0) then
!!           # 3-wave is leftgoing
!            sfract = s(3,i)
!        else
!!           # 3-wave is rightgoing
!            goto 200
!        endif
!        do 160 m=1,3
!            amdq(m,i) = amdq(m,i) + sfract*wave(m,3,i)
!        160 END DO
!    200 END DO

!!   compute rightgoing flux differences :

!    do 220 m=1,3
!        do 220 i = 2-mbc,mx+mbc
!            df = 0.0d0
!            do 210 mw=1,mwaves
!                df = df + s(mw,i)*wave(m,mw,i)
!            210 END DO
!            apdq(m,i)=df-amdq(m,i)
!    220 END DO


!    900 continue
    return
  end subroutine rpn2


