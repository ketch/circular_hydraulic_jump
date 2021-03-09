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

    implicit double precision (a-h,o-z)

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
    logical :: use_llf

    data use_llf /.true./    !# Use LLF solver

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

!   # now split the jump in q at each interface into waves
    do 20 i = 2-mbc, mx+mbc

       ! ***** maximum wave speed ***** !
       h_l = qr(1,i-1)
       h_r = ql(1,i)

       ! Velocity
       u_l = qr(mu,i-1) / qr(1,i-1)
       u_r = ql(mu,i  ) / ql(1,i  )
       v_l = qr(mv,i-1) / qr(1,i-1)
       v_r = ql(mv,i  ) / ql(1,i  )

       ! Sound speed
       c_l = dsqrt(grav*h_l)
       c_r = dsqrt(grav*h_r)
       
       hsqrt_l = dsqrt(qr(1,i-1))
       hsqrt_r = dsqrt(ql(1,i))
       hsq2 = hsqrt_l + hsqrt_r
       h_hat = 0.5*(h_l + h_r)
       u_hat = (u_l*hsqrt_l + u_r*hsqrt_r) / hsq2
       v_hat = (v_l*hsqrt_l + v_r*hsqrt_r) / hsq2
       c_hat = dsqrt(grav*h_hat)

       ! speed for HLLE solver
       s1 = min(u_l - c_l, u_hat - c_hat)
       s2 = max(u_r + c_r, u_hat + c_hat)

       ! speed for LLF solver
       if (use_llf) then
          dlambda_llf = max(abs(min(u_l - c_l, u_hat - c_hat)),abs(max(u_r + c_r, u_hat + c_hat)))
          s1 = -dlambda_llf
          s2 = dlambda_llf
       end if
       
       ! ***** middle state (taken from HLLE solver) ***** !
       h_m = (ql(mu,i) - qr(mu,i-1) - s2*ql(1,i) + s1*qr(1,i-1))/(s1-s2)
       hu_m = (ql(mu,i)*(u_r-s2) - qr(mu,i-1)*(u_l-s1) + 0.5*grav*(h_r**2 - h_l**2) ) / (s1-s2)
       hv_m = (ql(mv,i)*u_r - qr(mv,i-1)*u_l - s2*ql(mv,i) + s1*qr(mv,i-1))/(s1-s2)

       ! ***** left-going wave ***** !
       wave(1, 1,i) = h_m - h_l
       wave(mu,1,i) = hu_m - qr(mu,i-1)
       wave(mv,1,i) = hv_m - qr(mv,i-1)
       s(1,i) = s1

       ! ***** dummy wave ***** !
       wave(1, 2,i) = 0
       wave(mu,2,i) = 0
       wave(mv,2,i) = 0
       s(2,i) = 0.
       
       ! ***** right-going wave ***** !
       wave(1, 3,i) = h_r - h_m
       wave(mu,3,i) = ql(mu,i) - hu_m
       wave(mv,3,i) = ql(mv,i) - hv_m
       s(3,i) = s2
       
    20 END DO

!    # compute flux differences amdq and apdq.
!    ---------------------------------------
!     ----------------
!     # amdq = SUM s*wave   over left-going waves
!     # apdq = SUM s*wave   over right-going waves

    do m=1,3
       do i=2-mbc, mx+mbc
          amdq(m,i) = 0.d0
          apdq(m,i) = 0.d0

          do mw=1,mwaves
             amdq(m,i) = amdq(m,i) + MIN(0.d0,s(mw,i))*wave(m,mw,i)
             apdq(m,i) = apdq(m,i) + MAX(0.d0,s(mw,i))*wave(m,mw,i)
          END DO
       END DO
    END DO

    return
  end subroutine rpn2
