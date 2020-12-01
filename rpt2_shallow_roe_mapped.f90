subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)

! Transverse Roe solver for the 2D shallow water equations
! on a mapped grid.

! waves: 3
! equations: 3

! Conserved quantities:
!   1 depth
!   2 x-momentum
!   3 y-momentum

    implicit none

    double precision, dimension(meqn, 1-mbc:maxm+mbc), intent(in) :: ql, qr, asdq
    integer, intent(in) :: ixy, imp, maxm, meqn, mwaves, maux, mbc, mx
    double precision, dimension(maux, 1-mbc:maxm+mbc), intent(in) :: aux1, aux2, aux3
    double precision, dimension(meqn, 1-mbc:maxm+mbc), intent(out) :: bmasdq, bpasdq

    ! Local variables
    double precision, dimension(3,3) :: waveb
    double precision, dimension(3) :: sb
    integer :: mu, mv, i, m, mw, inx, iny, ilenrat
    double precision :: hsqrt_l, hsqrt_r, hsq2, a1, a2, a3, grav, u, v, a, h, alpha, beta
    double precision :: asdqh, asdqn, asdqt

!   # grav must be set elsewhere
    common /cparam/ grav

!   # Roe average quantities of each interface
    integer, parameter :: maxm2 = 1800
    double precision, dimension(-6:maxm2+7) :: unorl, unorr, utanl, utanr

    if (ixy == 1) then
        inx = 4
        iny = 5
        ilenrat = 6
    else
        inx = 1
        iny = 2
        ilenrat = 3
    endif
            
    do i=2-mbc, mx+mbc
        alpha = aux2(inx,i)
        beta = aux2(iny,i)
        unorl = alpha*ql(2,i) + beta*ql(3,i)
        unorr = alpha*qr(2,i-1) + beta*qr(3,i-1)
        utanl = -beta*ql(2,i) + alpha*ql(3,i)
        utanr = -beta*qr(2,i-1) + alpha*qr(3,i-1)
    enddo

    do i = 2-mbc, mx+mbc
        alpha = aux2(inx,i)
        beta = aux2(iny,i)
        h = (qr(1,i-1)+ql(1,i))*0.50d0
        hsqrt_l = dsqrt(qr(1,i-1))
        hsqrt_r = dsqrt(ql(1,i))
        hsq2 = hsqrt_l + hsqrt_r
        u = (unorr(i-1)/hsqrt_l + unorl(i)/hsqrt_r) / hsq2
        v = (utanr(i-1)/hsqrt_l + utanl(i)/hsqrt_r) / hsq2
        a = dsqrt(grav*h)

        asdqh = asdq(1,i)
        asdqt = -beta*asdq(2,i) + alpha*asdq(3,i)
        asdqn = alpha*asdq(2,i) + beta*asdq(3,i)

        a1 = (0.50d0/a)*((v+a)*asdqh-asdqt)
        a2 = asdqn - u*asdqh
        a3 = (0.50d0/a)*(-(v-a)*asdqh+asdqt)
    
        waveb(1,1) = a1
        waveb(2,1) = alpha*a1*u - beta*a1*(v-a)
        waveb(3,1) = beta*a1*u + alpha*a1*(v-a)
        sb(1) = (v - a) * aux2(ilenrat,i)

        waveb(1,2) = 0.0d0
        waveb(2,2) = alpha*a2
        waveb(3,2) = beta*a2
        sb(2) = v * aux2(ilenrat,i)
    
        waveb(1,3) = a3
        waveb(2,3) = alpha*a3*u - beta*a3*(v+a)
        waveb(3,3) = beta*a3*u + alpha*a3*(v+a)
        sb(3) = (v + a) * aux2(ilenrat,i)
    
!       # compute the flux differences bmasdq and bpasdq
    
        do m=1,meqn
            bmasdq(m,i) = 0.d0
            bpasdq(m,i) = 0.d0
            do mw=1,mwaves
                bmasdq(m,i) = bmasdq(m,i) + dmin1(sb(mw), 0.d0) * waveb(m,mw)
                bpasdq(m,i) = bpasdq(m,i) + dmax1(sb(mw), 0.d0) * waveb(m,mw)
            end do
        end do
    
    end do

    return
    end subroutine rpt2
