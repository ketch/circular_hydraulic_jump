! =====================================================
subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =====================================================

! HLLEM solver for the 2D shallow water equations.

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
    
    double precision, dimension(meqn, mwaves, 1-mbc:maxm+mbc), intent(out) :: wave
    double precision, dimension(mwaves, 1-mbc:maxm+mbc), intent(out) :: s
    double precision, dimension(meqn, 1-mbc:maxm+mbc), intent(out) :: amdq, apdq
    
    double precision, dimension(2,1) :: q_hat
    double precision :: u_l, u_r, c_l, c_r, u_hat, c_hat, v_l, v_r, v_star, v_hat, h_hat
    double precision :: h_l, h_r, hsqrt_l, hsqrt_r, hsq2
    double precision :: grav, a, b, c, aprime, bprime, root0, root
    double precision :: h_ml, hu_ml, hv_ml, h_mr, hu_mr, hv_mr, s1, s2, s_star, u_ml, u_mr
    integer :: depth, mu, mv
    integer :: i, m, mw, j, maxiter
    

    common /cparam/  grav

!     # set mu to point to  the component of the system that corresponds
!     # to momentum in the direction of this slice, mv to the orthogonal
!     # momentum:
!
    maxiter = 20
    depth = 1
    if (ixy.eq.1) then
        mu = 2
        mv = 3
    else
        mu = 3
        mv = 2
    endif

    do i=2-mbc,mx+mbc
        h_l = qr(depth,i-1)
        h_r = ql(depth,i)
        ! Velocity
        u_l = qr(mu,i-1) / qr(depth,i-1)
        u_r = ql(mu,i  ) / ql(depth,i  )
        v_l = qr(mv,i-1) / qr(depth,i-1)
        v_r = ql(mv,i  ) / ql(depth,i  )
        ! Sound speed
        c_l = dsqrt(grav*h_l)
        c_r = dsqrt(grav*h_r)
	
        hsqrt_l = dsqrt(qr(depth,i-1))
        hsqrt_r = dsqrt(ql(depth,i))
        hsq2 = hsqrt_l + hsqrt_r
        ! Roe averages
        h_hat = 0.5*(h_l + h_r)
        u_hat = (u_l*hsqrt_l + u_r*hsqrt_r) / hsq2
        v_hat = (v_l*hsqrt_l + v_r*hsqrt_r) / hsq2
        c_hat = dsqrt(grav*h_hat)

        ! Speeds of non-shear waves
        s1 = min(u_l - c_l, u_hat - c_hat)
        s2 = max(u_r + c_r, u_hat + c_hat)
        ! Speeds of shear wave
        s_star = (h_r*v_r*u_r-h_l*v_l*u_l)/(h_r*v_r-h_l*v_l)
        
        if (s_star>0) then
                v_star = v_l
        else
                v_star = v_r
        endif

        a = grav
        ! left middle state
        q_hat(1,1) = s1*h_l-h_l*u_l
        q_hat(2,1) = s1*h_l*u_l-h_l*u_l*u_l-0.5*grav*h_l*h_l
        b = -2*(s1*q_hat(1,1)-q_hat(2,1))
        c = 2*q_hat(1,1)*q_hat(1,1)
        aprime = 3*a
        bprime = b
        
        root0 = 0.0
        ! Newton-Raphson root finding
        do j=1, maxiter
             root = root0-((a*root0*root0*root0+b*root0+c)/(aprime*root0*root0+bprime))
             h_ml = root
             if (abs(root-root0) .lt. 0.00001d0) then                  
                  exit
             else
                  root0 = root
             endif
        end do
        u_ml = (q_hat(2,1)+0.5*grav*h_ml*h_ml)/q_hat(1,1)
        hu_ml = h_ml*u_ml
        hv_ml = h_ml*v_star
        
        ! right middle state
        q_hat(1,1) = s2*h_r-h_r*u_r
        q_hat(2,1) = s2*h_r*u_r-h_r*u_r*u_r-0.5*grav*h_r*h_r
        b = -2*(s2*q_hat(1,1)-q_hat(2,1))
        c = 2*q_hat(1,1)*q_hat(1,1)
        aprime = 3*a
        bprime = b
        
        root0 = 0.0
        ! Newton-Raphson root finding
        do j=1, maxiter
             root = root0-((a*root0*root0*root0+b*root0+c)/(aprime*root0*root0+bprime))
             h_mr = root
             if (abs(root-root0) .lt. 0.00001d0) then                  
                  exit
             else
                  root0 = root
             endif
        end do
        u_mr = (q_hat(2,1)+0.5*grav*h_mr*h_mr)/q_hat(1,1)
        hu_mr = h_mr*u_mr
        hv_mr = h_mr*v_star
        
	! 1-wave and its speed
        wave(depth,1,i) = h_ml - h_l
        wave(mu,1,i) = hu_ml - qr(mu,i-1)
        wave(mv,1,i) = hv_ml - qr(mv,i-1)
        s(1,i) = s1

	! 2-wave and its speed
        wave(depth,2,i) = h_mr - h_ml
        wave(mu,2,i) = hu_mr - hu_ml
        wave(mv,2,i) = hv_mr - hv_ml
        s(2,i) = s_star

	! 3-wave and its speed
        wave(depth,3,i) = h_r - h_mr
        wave(mu,3,i) = ql(mu,i) - hu_mr
        wave(mv,3,i) = ql(mv,i) - hv_mr
        s(3,i) = s2
    end do


    do m=1, meqn
        do i=2-mbc, mx+mbc
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do mw=1, mwaves
                if (s(mw,i) .lt. 0.d0) then
                    amdq(m,i) = amdq(m,i) + s(mw,i)*wave(m,mw,i)
                 else
                   apdq(m,i) = apdq(m,i) + s(mw,i)*wave(m,mw,i)
                 endif
            end do
        end do
    end do

end subroutine rpn2
