subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)

! Transverse Roe solver for the 2D shallow water equations on a mapped grid.
! This is just the solver https://depts.washington.edu/clawpack/clawpack-4.3/book/chap23/shallow/cylinder/rpt2swq.f
! but updated to Fortran 90 and Clawpack 5.

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
    double precision, dimension(meqn,mwaves,1-mbc:maxm+mbc) :: wave
    double precision, dimension(mwaves,1-mbc:maxm+mbc) :: s

    ! Local variables
    integer :: i, m, mw, inx, iny, ilenrat, ix1, ixm1, i1
    double precision :: a1, a2, a3, grav
    double precision, dimension(3) :: delta

    common /cparam/ grav

    ! Roe average quantities of each interface
    integer, parameter :: maxm2 = 1800
    double precision, dimension(-6:maxm2+7) :: alpha, beta, u, v, a, h

    if (-6.gt.1-mbc .or. maxm2 .lt. maxm+mbc) then
        write(6,*) 'need to increase maxm2 in rpB'
        stop
    endif

    if (ixy == 1) then
         inx = 4
         iny = 5
         ilenrat = 6
    else
         inx = 1
         iny = 2
         ilenrat = 3
    endif

    ! imp is used to flag whether wave is going to left or right,
    ! since states and grid are different on each side

      if (imp == 1) then
    ! asdq = amdq, moving to left
         ix1 = 2-mbc
         ixm1 = mx+mbc
      else
    ! asdq = apdq, moving to right
         ix1 = 1-mbc
         ixm1 = mx+mbc
      endif

    ! --------------
    ! up-going:
    ! --------------

    ! determine rotation matrix for interface above cell, using aux3
    ! [ alpha  beta ]
    ! [-beta  alpha ]

    do i=ix1,ixm1

        if (imp == 1) then
            i1 = i-1
        else
            i1 = i
        endif

         alpha(i) = aux3(inx,i1)
         beta(i) = aux3(iny,i1)
         h(i) = ql(1,i1)
         u(i) = (alpha(i)*ql(2,i1) + beta(i)*ql(3,i1)) / h(i)
         v(i) = (-beta(i)*ql(2,i1) + alpha(i)*ql(3,i1)) / h(i)
         a(i) = dsqrt(grav*h(i))
    enddo

    ! now split asdq into waves:

    ! find a1 thru a3, the coefficients of the 3 eigenvectors:
    do i = ix1,ixm1
         delta(1) = asdq(1,i) 
         delta(2) = alpha(i)*asdq(2,i) + beta(i)*asdq(3,i)
         delta(3) = -beta(i)*asdq(2,i) + alpha(i)*asdq(3,i)
         a1 = ((u(i)+a(i))*delta(1) - delta(2))*(0.50d0/a(i))
         a2 = -v(i)*delta(1) + delta(3)
         a3 = (-(u(i)-a(i))*delta(1) + delta(2))*(0.50d0/a(i))

    ! Compute the waves

         wave(1,1,i) = a1
         wave(2,1,i) = alpha(i)*a1*(u(i)-a(i)) - beta(i)*a1*v(i)
         wave(3,1,i) = beta(i)*a1*(u(i)-a(i)) + alpha(i)*a1*v(i)
         s(1,i) = (u(i)-a(i)) * aux3(ilenrat,i1)

         wave(1,2,i) = 0.0d0
         wave(2,2,i) = -beta(i)*a2
         wave(3,2,i) = alpha(i)*a2
         s(2,i) = u(i) * aux3(ilenrat,i1)

         wave(1,3,i) = a3
         wave(2,3,i) = alpha(i)*a3*(u(i)+a(i)) - beta(i)*a3*v(i)
         wave(3,3,i) = beta(i)*a3*(u(i)+a(i)) + alpha(i)*a3*v(i)
         s(3,i) = (u(i)+a(i)) * aux3(ilenrat,i1)
    end do

    ! compute flux difference bpasdq
    ! --------------------------------
    do m=1,3
        do i=ix1,ixm1
            bpasdq(m,i) = 0.d0
            do mw=1,mwaves
               bpasdq(m,i) = bpasdq(m,i) + dmax1(s(mw,i),0.d0)*wave(m,mw,i)
            end do
        end do
    end do

    ! ------------
    ! down-going:
    ! ------------
      
    ! # determine rotation matrix for interface below cell, using aux2
    ! [ alf  beta ]
    ! [-beta  alf ]

    do i=ix1,ixm1
        if (imp == 1) then
            i1 = i-1
        else
            i1 = i
        endif

        alpha(i) = aux2(inx,i1)
        beta(i) = aux2(iny,i1)
        u(i) = (alpha(i)*ql(2,i1) + beta(i)*ql(3,i1)) / h(i)
        v(i) = (-beta(i)*ql(2,i1) + alpha(i)*ql(3,i1)) / h(i)
    end do

    ! now split asdq into waves:

    ! find a1 thru a3, the coefficients of the 3 eigenvectors:
    do i = ix1,ixm1
        delta(1) = asdq(1,i) 
        delta(2) = alpha(i)*asdq(2,i) + beta(i)*asdq(3,i)
        delta(3) = -beta(i)*asdq(2,i) + alpha(i)*asdq(3,i)
        a1 = ((u(i)+a(i))*delta(1) - delta(2))*(0.50d0/a(i))
        a2 = -v(i)*delta(1) + delta(3)
        a3 = (-(u(i)-a(i))*delta(1) + delta(2))*(0.50d0/a(i))

        ! Compute the waves

        wave(1,1,i) = a1
        wave(2,1,i) = alpha(i)*a1*(u(i)-a(i)) - beta(i)*a1*v(i)
        wave(3,1,i) = beta(i)*a1*(u(i)-a(i)) + alpha(i)*a1*v(i)
        s(1,i) = (u(i)-a(i)) * aux2(ilenrat,i1)

        wave(1,2,i) = 0.0d0
        wave(2,2,i) = -beta(i)*a2
        wave(3,2,i) = alpha(i)*a2
        s(2,i) = u(i) * aux2(ilenrat,i1)

        wave(1,3,i) = a3
        wave(2,3,i) = alpha(i)*a3*(u(i)+a(i)) - beta(i)*a3*v(i)
        wave(3,3,i) = beta(i)*a3*(u(i)+a(i)) + alpha(i)*a3*v(i)
        s(3,i) = (u(i)+a(i)) * aux2(ilenrat,i1)
    end do

    ! compute flux difference bmasdq
    ! -------------------------------
    do m=1,3
        do i=ix1,ixm1
            bmasdq(m,i) = 0.d0
            do mw=1,mwaves
                bmasdq(m,i) = bmasdq(m,i) + dmin1(s(mw,i), 0.d0)*wave(m,mw,i)
            end do
        end do
    end do

    return
end
