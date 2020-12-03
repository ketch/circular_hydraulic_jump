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

    double precision, dimension(mwaves) :: splus, sminus
    integer :: depth
    integer :: i, m, mw
    integer :: inx, iny, ilenrat
    double precision :: ffroude, a1, a2, a3

    integer, parameter :: maxm2 = 1602  ! assumes at most 1000x1000 grid with mbc=2
    double precision, dimension(3) :: delta
    logical :: efix, cfix
    double precision, dimension(-1:maxm2) :: u, v, a, h, hunorl, hunorr, hutanl, hutanr, alpha, beta
    double precision :: h_l, h_r, u_l, u_r, v_l, v_r, hsqrt_l, hsqrt_r, hsq2, c_l, c_r
    double precision :: grav
    double precision :: s1l, s3r, saux, sm1fix, sm1roe, sp3fix, sp3roe, zfroude
    double precision :: renore, rere1, reremu, reremv, rhind
    double precision :: kalpha, kbeta, kepsilon  ! Kemm's parameters for carbuncle fix

    common /cparam/ grav, kalpha, kbeta, kepsilon
    common /comroe/ u, v, a, h

    data efix /.true./        !# use entropy fix for transonic rarefactions
    data cfix /.true./        !# additionally use carbuncle fix

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
    ! The variable names here are potentially confusing and should be changed.

    do i=2-mbc, mx+mbc
        alpha(i) = auxl(inx,i)
        beta(i) = auxl(iny,i)
        hunorl(i) = alpha(i)*ql(2,i) + beta(i)*ql(3,i)
        hunorr(i-1) = alpha(i)*qr(2,i-1) + beta(i)*qr(3,i-1)
        hutanl(i) = -beta(i)*ql(2,i) + alpha(i)*ql(3,i)
        hutanr(i-1) = -beta(i)*qr(2,i-1) + alpha(i)*qr(3,i-1)
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
        u(i) = (hunorr(i-1)/hsqrt_l + hunorl(i)/hsqrt_r) / hsq2
        v(i) = (hutanr(i-1)/hsqrt_l + hutanl(i)/hsqrt_r) / hsq2
        a(i) = dsqrt(grav*h(i)) ! = c_hat
    enddo

    ! now split the jump in q at each interface into waves

    ! find a1 thru a3, the coefficients of the 3 eigenvectors:
    do i = 2-mbc, mx+mbc
        delta(1) = ql(depth,i) - qr(depth,i-1)
        delta(2) = hunorl(i) - hunorr(i-1)
        delta(3) = hutanl(i) - hutanr(i-1)
        a1 = ((u(i)+a(i))*delta(1) - delta(2))*(0.50d0/a(i))
        a2 = -v(i)*delta(1) + delta(3)
        a3 = (-(u(i)-a(i))*delta(1) + delta(2))*(0.50d0/a(i))

       ! Compute the waves.

        wave(1,1,i) = a1
        wave(2,1,i) = alpha(i)*a1*(u(i)-a(i)) - beta(i)*a1*v(i)
        wave(3,1,i) = beta(i)*a1*(u(i)-a(i)) + alpha(i)*a1*v(i)
        s(1,i) = (u(i)-a(i)) * auxl(ilenrat,i)

        wave(1,2,i) = 0.0d0
        wave(2,2,i) = -beta(i)*a2
        wave(3,2,i) = alpha(i)*a2
        s(2,i) = u(i) * auxl(ilenrat,i)

        wave(1,3,i) = a3
        wave(2,3,i) = alpha(i)*a3*(u(i)+a(i)) - beta(i)*a3*v(i)
        wave(3,3,i) = beta(i)*a3*(u(i)+a(i)) + alpha(i)*a3*v(i)
        s(3,i) = (u(i)+a(i)) * auxl(ilenrat,i)
    enddo

    ! compute flux differences amdq and apdq.

    if (efix) then

    ! With entropy fix

    ! compute flux differences amdq and apdq.
    ! First compute amdq as sum of s*wave for left going waves.
    ! Incorporate entropy fix by adding a modified fraction of wave
    ! if s should change sign.

        do i=2-mbc,mx+mbc

            ! modified wave speeds to use HLLEM-entropy fix

            h_l = qr(depth,i-1)
            h_r = ql(depth,i)
            u_l = hunorr(i-1)/h_l
            u_r = hunorl(i)/h_l
            c_l = dsqrt(grav*h_l)
            c_r = dsqrt(grav*h_r)

            s1l  = (u_l - c_l) * auxl(ilenrat,i) ! u-c in left state (cell i-1)
            s3r  = (u_r + c_r) * auxl(ilenrat,i) ! u+c in right state  (cell i)

            sm1fix     = dmin1(s(1,i),0.0,s1l)
            sm1roe     = dmin1(s(1,i),0.0)
            sminus(1)  = 0.5d0*(sm1fix+sm1roe)

            sp3fix     = dmax1(s(3,i),0.0,s3r) 
            sp3roe     = dmax1(s(3,i),0.0)
            splus(3)   = 0.5d0*(sp3fix+sp3roe)

            splus(1)  = s(1,i) - sminus(1)
            sminus(3) = s(3,i) - splus(3) 


            !  carbuncle cure:
            ! -----------------

            if(cfix) then
                ! residual relative to acoustic speed
                rere1  = wave(1,3,i)-wave(1,1,i)
                reremu = wave(2,3,i)-wave(2,1,i)
                reremv = wave(3,3,i)-wave(3,1,i)
                ! norm of relative residual
                renore = dsqrt(rere1*rere1 + reremu*reremu + reremv*reremv)
                ! rhind = a(i) * max(kepsilon*renore,1.)
                zfroude = dsqrt((u(i)*u(i))/(a(i)*a(i)))
                ffroude  = dmax1(0.0d0,(1.0d0 - zfroude**kalpha))
                ! Alternative if cbrt for cubic root is available:
                ! ffroude  = dmax1(0.0d0,(1.0d0 - dcbrt(zfroude)))
               
                rhind = a(i) * auxl(ilenrat,i)
                rhind = rhind*dmin1(kepsilon*renore*ffroude,1.0d0)**kbeta 
                ! Alternative if cbrt for cubic root is available:
                ! rhind = rhind*dcbrt(dmin1(kepsilon*renore*ffroude,1.0d0))
               
                saux       = .5*(dabs(s(2,i)-rhind)+dabs(s(2,i)+rhind))

                sminus(2)  = .5*(s(2,i) - saux)
                splus(2)   = s(2,i) - sminus(2)
            else
                sminus(2) = dmin1(s(2,i),0.d0)
                splus(2)  = dmax1(s(2,i),0.d0)
            endif

            do m=1,meqn
                amdq(m,i) = 0.d0
                apdq(m,i) = 0.d0
                do mw=1,mwaves
                    amdq(m,i) = amdq(m,i) + sminus(mw)*wave(m,mw,i)
                    apdq(m,i) = apdq(m,i) + splus(mw)*wave(m,mw,i)
                enddo
            enddo
        enddo

    else ! no efix

        ! amdq = SUM s*wave   over left-going waves
        ! apdq = SUM s*wave   over right-going waves

        do m=1,3
            do i=2-mbc, mx+mbc
                amdq(m,i) = 0.d0
                apdq(m,i) = 0.d0
                do mw=1,mwaves
                    if (s(mw,i) .lt. 0.d0) then
                        amdq(m,i) = amdq(m,i) + s(mw,i)*wave(m,mw,i)
                    else
                        apdq(m,i) = apdq(m,i) + s(mw,i)*wave(m,mw,i)
                   endif
                enddo
            enddo
        enddo
    endif
    
end subroutine rpn2
