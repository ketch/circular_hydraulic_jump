subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! Modification of Kemm's code.
! The entropy fix for the 1- and 3-wave is performed as in rpn2_shallow_roe_with_efix.f90, 
! whereas the one for the 2-wave is the same as in Kemm's code.


! Roe-solver for the 2D shallow water equations on a quadrilateral grid
! 
! Solve Riemann problems along one slice of data.
! 
! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell
!     
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
! 
! Note that the i'th Riemann problem has left state qr(i-1,:)
!                                     and right state ql(i,:)
! From the basic clawpack routines, this routine is called with ql = qr

    implicit double precision (a-h,o-z)

    dimension wave(meqn, mwaves, 1-mbc:maxm+mbc)
    dimension    s(mwaves, 1-mbc:maxm+mbc)
    dimension   ql(meqn, 1-mbc:maxm+mbc)
    dimension   qr(meqn, 1-mbc:maxm+mbc)
    dimension  apdq(meqn, 1-mbc:maxm+mbc)
    dimension  amdq(meqn, 1-mbc:maxm+mbc)
    dimension auxl(maux, 1-mbc:maxm+mbc)
    dimension auxr(maux, 1-mbc:maxm+mbc)

    ! local arrays -- common block comroe is passed to rpt2sh
    ! ------------

    dimension splus(mwaves)
    dimension sminus(mwaves)

    parameter (maxm2 = 1602)  !# assumes at most 1000x1000 grid with mbc=2
    dimension delta(3)
    logical efix, cfix
    dimension unorl(-1:maxm2), unorr(-1:maxm2)
    dimension utanl(-1:maxm2), utanr(-1:maxm2)
    dimension alf(-1:maxm2)
    dimension beta(-1:maxm2)

    common /cparam/  grav
    common /comroe/ u(-1:maxm2),v(-1:maxm2),a(-1:maxm2),h(-1:maxm2)

    data efix /.true./        !# use entropy fix for transonic rarefactions
    !data efix /.false./        !# use entropy fix for transonic rarefactions
    data cfix /.true./        !# additionally use carbuncle fix

    if (-1.gt.1-mbc .or. maxm2 .lt. maxm+mbc) then
        write(6,*) 'need to increase maxm2 in rpA'
        stop
    endif

    ! rotate the velocities q(2) and q(3) so that it is aligned with grid
    ! normal.  The normal vector for the face at the i'th Riemann problem
    ! is stored in the aux array
    ! in locations (1,2) if ixy=1 or (4,5) if ixy=2.  The ratio of the
    ! length of the cell side to the length of the computational cell
    ! is stored in aux(3) or aux(6) respectively.

    if (ixy.eq.1) then
        inx = 1
        iny = 2
        ilenrat = 3
    else
        inx = 4
        iny = 5
        ilenrat = 6
    endif
    
    if (ixy == 1) then
        mu = 2
        mv = 3
    else
        mu = 3
        mv = 2
    endif

    ! determine rotation matrix
    !     [ alf  beta ]
    !     [-beta  alf ]

    ! note that this reduces to identity on standard cartesian grid

    do i=2-mbc,mx+mbc
        alf(i) = 1.d0 !auxl(inx,i)
        beta(i) = 0.d0 !auxl(iny,i)

        unorl(i) = alf(i)*ql(mu,i) + beta(i)*ql(mv,i)
        unorr(i-1) = alf(i)*qr(mu,i-1) + beta(i)*qr(mv,i-1)
        utanl(i) = -beta(i)*ql(mu,i) + alf(i)*ql(mv,i)
        utanr(i-1) = -beta(i)*qr(mu,i-1) + alf(i)*qr(mv,i-1)
    enddo

    ! compute the Roe-averaged variables needed in the Roe solver.
    ! These are stored in the common block comroe since they are
    ! later used in routine rpt2 to do the transverse wave splitting.

    do i = 2-mbc, mx+mbc
        h(i) = (qr(1,i-1)+ql(1,i))*0.50d0
        hsqrtl = dsqrt(qr(1,i-1))
        hsqrtr = dsqrt(ql(1,i))
        hsq2 = hsqrtl + hsqrtr
        ! unorr(i-1)/hsqrtl is fine because unorr = momentum 
        u(i) = (unorr(i-1)/hsqrtl + unorl(i)/hsqrtr) / hsq2
        v(i) = (utanr(i-1)/hsqrtl + utanl(i)/hsqrtr) / hsq2
        a(i) = dsqrt(grav*h(i))
    enddo

    ! now split the jump in q at each interface into waves

    ! find a1 thru a3, the coefficients of the 3 eigenvectors:
    do i = 2-mbc, mx+mbc
        delta(1) = ql(1,i) - qr(1,i-1)
        delta(2) = unorl(i) - unorr(i-1)
        delta(3) = utanl(i) - utanr(i-1)
        a1 = ((u(i)+a(i))*delta(1) - delta(2))*(0.50d0/a(i))
        a2 = -v(i)*delta(1) + delta(3)
        a3 = (-(u(i)-a(i))*delta(1) + delta(2))*(0.50d0/a(i))

       ! Compute the waves.

        wave(1,1,i) = a1
        wave(mu,1,i) = alf(i)*a1*(u(i)-a(i)) - beta(i)*a1*v(i)
        wave(mv,1,i) = beta(i)*a1*(u(i)-a(i)) + alf(i)*a1*v(i)
        s(1,i) = (u(i)-a(i)) !* auxl(ilenrat,i)

        wave(1,2,i) = 0.0d0
        wave(mu,2,i) = -beta(i)*a2
        wave(mv,2,i) = alf(i)*a2
        s(2,i) = u(i) !* auxl(ilenrat,i)

        wave(1,3,i) = a3
        wave(mu,3,i) = alf(i)*a3*(u(i)+a(i)) - beta(i)*a3*v(i)
        wave(mv,3,i) = beta(i)*a3*(u(i)+a(i)) + alf(i)*a3*v(i)
        s(3,i) = (u(i)+a(i)) !* auxl(ilenrat,i)
    enddo

    ! compute flux differences amdq and apdq.

    if (efix) then

    ! With entropy fix

    ! compute flux differences amdq and apdq.
    ! First compute amdq as sum of s*wave for left going waves.
    ! Incorporate entropy fix by adding a modified fraction of wave
    ! if s should change sign.

        do 200 i=2-mbc,mx+mbc
!           # check 1-wave
            him1 = qr(1,i-1)
            s0 =  qr(mu,i-1)/him1 - dsqrt(grav*him1)
!           # check for fully supersonic case :
            if (s0 > 0.0d0 .AND. s(1,i) > 0.0d0) then
                do 60 m=1,3
                    amdq(m,i)=0.0d0
                60 END DO
                goto 200
            endif
    
            h1 = qr(1,i-1)+wave(1,1,i)
            hu1= qr(mu,i-1)+wave(mu,1,i)
            s1 = hu1/h1 - dsqrt(grav*h1) !speed just to right of 1-wave
            if (s0 < 0.0d0 .AND. s1 > 0.0d0) then
!               # transonic rarefaction in 1-wave
                sfract = s0*((s1-s(1,i))/(s1-s0))
            else if (s(1,i) < 0.0d0) then
!               # 1-wave is leftgoing
                sfract = s(1,i)
            else
!               # 1-wave is rightgoing
                sfract = 0.0d0
            endif
            do 120 m=1,3
                amdq(m,i) = sfract*wave(m,1,i)
            120 END DO

!           check 3-wave
    
            hi = ql(1,i)
            s03 = ql(mu,i)/hi + dsqrt(grav*hi)
            h3=ql(1,i)-wave(1,3,i)
            hu3=ql(mu,i)-wave(mu,3,i)
            s3=hu3/h3 + dsqrt(grav*h3)
            if (s3 < 0.0d0 .AND. s03 > 0.0d0) then
!               # transonic rarefaction in 3-wave
                sfract = s3*((s03-s(3,i))/(s03-s3))
            else if (s(3,i) < 0.0d0) then
!               # 3-wave is leftgoing
                sfract = s(3,i)
            else
!               # 3-wave is rightgoing
                goto 200
            endif
            do 160 m=1,3
                amdq(m,i) = amdq(m,i) + sfract*wave(m,3,i)
            160 END DO
        200 END DO

!       compute rightgoing flux differences :

        do m=1,3
            do i = 2-mbc,mx+mbc
                df = 0.0d0
                df = df + s(1,i)*wave(m,1,i)
                df = df + s(3,i)*wave(m,3,i)
                apdq(m,i)=df-amdq(m,i)
            END DO
        END DO


        do i=2-mbc,mx+mbc

            !  carbuncle cure:
            ! -----------------

            if(cfix) then
                ! residual relative to acoustic speed
                rere1  = wave(1,3,i)-wave(1,1,i)
                reremu = wave(2,3,i)-wave(2,1,i)
                reremv = wave(3,3,i)-wave(3,1,i)
                ! norm of relative residual
                renore = dsqrt(rere1*rere1 + reremu*reremu + reremv*reremv)
                ! Indicator for Rankine-Hugoniot condition, original:0.01d0 
                rhparam = 0.001d0
                ! rhind = a(i) * max(rhparam*renore,1.)
                zfroude = dsqrt((u(i)*u(i))/(a(i)*a(i)))
                ffroude  = dmax1(0.0d0,(1.0d0 - zfroude**(1.0d0/3.0d0)))
                ! Alternative if cbrt for cubic root is available:
                ! ffroude  = dmax1(0.0d0,(1.0d0 - dcbrt(zfroude)))
               
                rhind = a(i) * auxl(ilenrat,i)
                rhind = rhind*dmin1(rhparam*renore*ffroude,1.0d0)**(1.0d0/3.0d0) 
                ! Alternative if cbrt for cubic root is available:
                ! rhind = rhind*dcbrt(dmin1(rhparam*renore*ffroude,1.0d0))
               
                saux       = .5*(dabs(s(2,i)-rhind)+dabs(s(2,i)+rhind))

                sminus(2)  = .5*(s(2,i) - saux)
                splus(2)   = s(2,i) - sminus(2)
            endif

            do m=1,meqn
                amdq(m,i) = amdq(m,i) + sminus(2)*wave(m,2,i)
                apdq(m,i) = apdq(m,i) + splus(2)*wave(m,2,i)
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
