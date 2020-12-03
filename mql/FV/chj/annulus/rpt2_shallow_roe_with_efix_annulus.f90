! =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! =====================================================
    implicit none !double precision (a-h,o-z)

!     # Riemann solver in the transverse direction for the shallow water
!     equations .
!     # Split asdq (= A^* \Delta q, where * = + or -)
!     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
!     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)

    integer, intent(in) :: ixy, imp, maxm, meqn, mwaves, maux, mbc, mx
    double precision, dimension(meqn, 1-mbc:maxm+mbc), intent(in) :: ql,qr
    double precision, dimension(maux,1-mbc:maxm+mbc), intent(in) :: aux1, aux2, aux3
    double precision, dimension(meqn, 1-mbc:maxm+mbc), intent(in) :: asdq
    double precision, dimension(meqn, 1-mbc:maxm+mbc), intent(out) :: bmasdq, bpasdq

    double precision, dimension (3,3) :: wave
    double precision, dimension (3) :: s, delta
    double precision :: grav

    integer :: mw,m,inx,iny,ilenrat,i
    double precision :: hsqrt_r,hsqrt_l,hsq2,a1,a2,a3
    double precision :: w11,w12,w13,w21,w22,w23,w31,w32,w33,d1,d2,d3
    
!   # grav must be set elsewhere
    common /cparam/ grav

!   # Roe averages quantities of each interface
    integer, parameter :: maxm2 = 1800
    double precision, dimension (-6:maxm2+7) :: u, v, a, h, hunor_r, hunor_l, hutan_r, hutan_l, alpha, beta

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
       alpha(i) = aux2(inx,i)
       beta(i) = aux2(iny,i)
       hunor_r(i)   =  alpha(i) * ql(2,i)   + beta(i)  * ql(3,i)
       hunor_l(i-1) =  alpha(i) * qr(2,i-1) + beta(i)  * qr(3,i-1)
       hutan_r(i)   = -beta(i)  * ql(2,i)   + alpha(i) * ql(3,i)
       hutan_l(i-1) = -beta(i)  * qr(2,i-1) + alpha(i) * qr(3,i-1)
    enddo

    ! compute the Roe-averaged variables needed in the Roe solver.
    do i = 2-mbc, mx+mbc
        hsqrt_l = dsqrt(qr(1,i-1))
        hsqrt_r = dsqrt(ql(1,i))
        hsq2 = hsqrt_l + hsqrt_r
        
        ! Roe averages
        h(i) = (qr(1,i-1)+ql(1,i))*0.50d0 ! = h_hat
        u(i) = (hunor_l(i-1)/hsqrt_l + hunor_r(i)/hsqrt_r) / hsq2
        v(i) = (hutan_l(i-1)/hsqrt_l + hutan_r(i)/hsqrt_r) / hsq2
        a(i) = dsqrt(grav*h(i)) ! = c_hat
     enddo
     
     do 20 i = 2-mbc, mx+mbc
        ! compute delta and then rotate it
        d1 = asdq(1,i)
        d2 = asdq(2,i)
        d3 = asdq(3,i)
        delta(1) = d1
        delta(2) = alpha(i)*d2 + beta(i)*d3
        delta(3) = -beta(i)*d2 + alpha(i)*d3

        ! compute a coefficients
        a1 = ((u(i)+a(i))*delta(1) - delta(2))*(0.50d0/a(i))
        a2 = -v(i)*delta(1) + delta(3)
        a3 = (-(u(i)-a(i))*delta(1) + delta(2))*(0.50d0/a(i))

        ! compute the Roe waves
        w11 = a1
        w12 = a1*(u(i)-a(i))
        w13 = a1*v(i)

        w21 = 0.d0
        w22 = 0.d0
        w23 = a2

        w31 = a3
        w32 = a3*(u(i)+a(i)) 
        w33 = a3*v(i)

        ! compute speed of Roe waves
        s(1) = u(i) - a(i)
        s(2) = u(i)
        s(3) = u(i) + a(i)
        
        ! rotate back the Roe waves
        wave(1,1) = w11
        wave(2,1) = alpha(i)*w12 - beta(i)*w13
        wave(3,1) = beta(i)*w12 + alpha(i)*w13
        
        wave(1,2) = w21 
        wave(2,2) = alpha(i)*w22 - beta(i)*w23 
        wave(3,2) = beta(i)*w22 + alpha(i)*w23 
        
        wave(1,3) = w31 
        wave(2,3) = alpha(i)*w32 - beta(i)*w33 
        wave(3,3) = beta(i)*w32 + alpha(i)*w33

        ! scale speeds due to rotation
        s(1) = s(1) * aux2(ilenrat,i)
        s(2) = s(2) * aux2(ilenrat,i)
        s(3) = s(3) * aux2(ilenrat,i)
        
!       # compute the flux differences bmasdq and bpasdq    
        do 30 m=1,meqn
            bmasdq(m,i) = 0.d0
            bpasdq(m,i) = 0.d0
            do 30 mw=1,mwaves
               bmasdq(m,i) = bmasdq(m,i) + dmin1(s(mw), 0.d0) * wave(m,mw)
               bpasdq(m,i) = bpasdq(m,i) + dmax1(s(mw), 0.d0) * wave(m,mw)
        30 END DO
    
    20 END DO

    return
    end subroutine rpt2
