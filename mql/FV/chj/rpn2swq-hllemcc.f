
c
c
c     =====================================================
      subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &     wave,s,amdq,apdq)
c     =====================================================
c     
c     # Roe-solver for the 2D shallow water equations
c     #  on a quadrilateral grid
c     
c     # solve Riemann problems along one slice of data.
c     
c     # On input, ql contains the state vector at the left edge of each
c     cell
c     #           qr contains the state vector at the right edge of each
c     cell
c     
c     # This data is along a slice in the x-direction if ixy=1
c     #                            or the y-direction if ixy=2.
c     # On output, wave contains the waves, s the speeds,
c     # and amdq, apdq the decomposition of the flux difference
c     #   f(qr(i-1)) - f(ql(i))
c     # into leftgoing and rightgoing parts respectively.
c     # With the Roe solver we have
c     #    amdq  =  A^- \Delta q    and    apdq  =  A^+ \Delta q
c     # where A is the Roe matrix.  An entropy fix can also be
c     incorporated
c     # into the flux differences.
c     
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routines, this routine is called with ql
c     = qr
c     
c     
      implicit double precision (a-h,o-z)
c     
      dimension wave(1-mbc:maxm+mbc, meqn, mwaves)
      dimension    s(1-mbc:maxm+mbc, mwaves)
      dimension   ql(1-mbc:maxm+mbc, meqn)
      dimension   qr(1-mbc:maxm+mbc, meqn)
      dimension  apdq(1-mbc:maxm+mbc, meqn)
      dimension  amdq(1-mbc:maxm+mbc, meqn)
      dimension auxl(1-mbc:maxm+mbc, 7)
      dimension auxr(1-mbc:maxm+mbc, 7)
c     
c     local arrays -- common block comroe is passed to rpt2sh
c     ------------
c     
      dimension splus(mwaves)
      dimension sminus(mwaves)
c     
      parameter (maxm2 = 1602)  !# assumes at most 1000x1000 grid with mbc=2
      dimension delta(3)
      logical efix, cfix
      dimension unorl(-1:maxm2), unorr(-1:maxm2)
      dimension utanl(-1:maxm2), utanr(-1:maxm2)
      dimension alf(-1:maxm2)
      dimension beta(-1:maxm2)

      common /sw/  g
      common /comroe/ u(-1:maxm2),v(-1:maxm2),a(-1:maxm2),h(-1:maxm2)
c     
      data efix /.true./        !# use entropy fix for transonic rarefactions
      data cfix /.true./        !# aditionally use carbuncle fix
c     
      if (-1.gt.1-mbc .or. maxm2 .lt. maxm+mbc) then
         write(6,*) 'need to increase maxm2 in rpA'
         stop
      endif
c     
c     
c     # rotate the velocities q(2) and q(3) so that it is aligned with
c     grid
c     # normal.  The normal vector for the face at the i'th Riemann
c     problem
c     # is stored in the aux array
c     # in locations (1,2) if ixy=1 or (4,5) if ixy=2.  The ratio of the
c     # length of the cell side to the length of the computational cell
c     # is stored in aux(3) or aux(6) respectively.
c     
c     
      if (ixy.eq.1) then
         inx = 1
         iny = 2
         ilenrat = 3
      else
         inx = 4
         iny = 5
         ilenrat = 6
      endif
c     
c     # determine rotation matrix
c     [ alf  beta ]
c     [-beta  alf ]
c     
c     # note that this reduces to identity on standard cartesian grid
c     
      do i=2-mbc,mx+mbc
         alf(i) = auxl(i,inx)
         beta(i) = auxl(i,iny)
         unorl(i) = alf(i)*ql(i,2) + beta(i)*ql(i,3)
         unorr(i-1) = alf(i)*qr(i-1,2) + beta(i)*qr(i-1,3)
         utanl(i) = -beta(i)*ql(i,2) + alf(i)*ql(i,3)
         utanr(i-1) = -beta(i)*qr(i-1,2) + alf(i)*qr(i-1,3)
      enddo
c     
c     
c     # compute the Roe-averaged variables needed in the Roe solver.
c     # These are stored in the common block comroe since they are
c     # later used in routine rpt2 to do the transverse wave splitting.
c     
      do 10 i = 2-mbc, mx+mbc
         h(i) = (qr(i-1,1)+ql(i,1))*0.50d0
         hsqrtl = dsqrt(qr(i-1,1))
         hsqrtr = dsqrt(ql(i,1))
         hsq2 = hsqrtl + hsqrtr
         u(i) = (unorr(i-1)/hsqrtl + unorl(i)/hsqrtr) / hsq2
         v(i) = (utanr(i-1)/hsqrtl + utanl(i)/hsqrtr) / hsq2
         a(i) = dsqrt(g*h(i))
 10   continue
c     
c     
c     # now split the jump in q at each interface into waves
c     
c     # find a1 thru a3, the coefficients of the 3 eigenvectors:
      do 20 i = 2-mbc, mx+mbc
         delta(1) = ql(i,1) - qr(i-1,1)
         delta(2) = unorl(i) - unorr(i-1)
         delta(3) = utanl(i) - utanr(i-1)
         a1 = ((u(i)+a(i))*delta(1) - delta(2))*(0.50d0/a(i))
         a2 = -v(i)*delta(1) + delta(3)
         a3 = (-(u(i)-a(i))*delta(1) + delta(2))*(0.50d0/a(i))
c     
c     # Compute the waves.
c     
         wave(i,1,1) = a1
         wave(i,2,1) = alf(i)*a1*(u(i)-a(i)) - beta(i)*a1*v(i)
         wave(i,3,1) = beta(i)*a1*(u(i)-a(i)) + alf(i)*a1*v(i)
         s(i,1) = (u(i)-a(i)) * auxl(i,ilenrat)
c     
         wave(i,1,2) = 0.0d0
         wave(i,2,2) = -beta(i)*a2
         wave(i,3,2) = alf(i)*a2
         s(i,2) = u(i) * auxl(i,ilenrat)
c     
         wave(i,1,3) = a3
         wave(i,2,3) = alf(i)*a3*(u(i)+a(i)) - beta(i)*a3*v(i)
         wave(i,3,3) = beta(i)*a3*(u(i)+a(i)) + alf(i)*a3*v(i)
         s(i,3) = (u(i)+a(i)) * auxl(i,ilenrat)
 20   continue
c     
c     
c     # compute flux differences amdq and apdq.
c     ---------------------------------------
c     
      if (efix) go to 110
c     
c     # no entropy fix
c     ----------------
c     
c     # amdq = SUM s*wave   over left-going waves
c     # apdq = SUM s*wave   over right-going waves
c     
      do 100 m=1,3
         do 100 i=2-mbc, mx+mbc
            amdq(i,m) = 0.d0
            apdq(i,m) = 0.d0
            do 90 mw=1,mwaves
               if (s(i,mw) .lt. 0.d0) then
                  amdq(i,m) = amdq(i,m) + s(i,mw)*wave(i,m,mw)
               else
                  apdq(i,m) = apdq(i,m) + s(i,mw)*wave(i,m,mw)
               endif
 90         continue
 100     continue
         go to 900
c     
c-----------------------------------------------------
c     
 110     continue
c     
c     # With entropy fix
c     ------------------
c     
c     # compute flux differences amdq and apdq.
c     # First compute amdq as sum of s*wave for left going waves.
c     # Incorporate entropy fix by adding a modified fraction of wave
c     # if s should change sign.
c     
         do 200 i=2-mbc,mx+mbc

c     modified wave speeds to use HLLEM-entropy fix

            him1 = qr(i-1,1)
            uim1 = unorr(i-1)
            cim1 = dsqrt(g*him1)
            s1l  = (uim1 - cim1) * auxl(i,ilenrat) !# u-c in left state (cell i-1)

            hi   = ql(i,1)
            ui   = unorl(i)
            ci   = dsqrt(g*hi)
            s3r  = (ui + ci) * auxl(i,ilenrat) !# u+c in right state  (cell i)

            sm1fix     = dmin1(s(i,1),0.0,s1l)
            sm1roe     = dmin1(s(i,1),0.0)
            sminus(1)  = 0.5d0*(sm1fix+sm1roe)

            sp3fix     = dmax1(s(i,3),0.0,s3r) 
            sp3roe     = dmax1(s(i,3),0.0)
            splus(3)   = 0.5d0*(sp3fix+sp3roe)


            splus(1)  = s(i,1)-sminus(1)
            sminus(3) = s(i,3) - splus(3) 


           
c     
c     # carbuncle cure:
c     -----------------
c     
            if(cfix) then
c     residual relative to acoustic speed
               rere1  = wave(i,1,3)-wave(i,1,1)
               reremu = wave(i,2,3)-wave(i,2,1)
               reremv = wave(i,3,3)-wave(i,3,1)
c     norm of relative residual
               renore = dsqrt(rere1*rere1 + reremu*reremu +
     $              reremv*reremv)
c     Indicator for Rankine-Hugoniot condition, original:0.01d0 
               rhparam = 0.001d0
c     rhind = a(i) * max(rhparam*renore,1.)
               zfroude = dsqrt((u(i)*u(i))/(a(i)*a(i)))
               ffroude  = dmax1(0.0d0,(1.0d0 - zfroude**(1.0d0/3.0d0))
c     Alternative if cbrt for cubic root is available:
c     ffroude  = dmax1(0.0d0,(1.0d0 - dcbrt(zfroude)))
               
               rhind = a(i) * auxl(i,ilenrat)
               rhind = rhind*dmin1(rhparam*renore*ffroude,1.0d0)**(1.0d0
     $              /3.0d0) 
c     Alternative if cbrt for cubic root is available:
c     rhind = rhind*dcbrt(dmin1(rhparam*renore*ffroude,1.0d0))
               
               saux       = .5*(dabs(s(i,2)-rhind)+dabs(s(i,2)+rhind))

               sminus(2)  = .5*(s(i,2) - saux)
               splus(2)   = s(i,2) - sminus(2)
            endif


            


            do 120 m=1,meqn
               amdq(i,m) = 0.d0
               apdq(i,m) = 0.d0
               do 95 mw=1,mwaves
                  amdq(i,m) = amdq(i,m) + sminus(mw)*wave(i,m,mw)
                  apdq(i,m) = apdq(i,m) + splus(mw)*wave(i,m,mw)
 95            continue
 120        continue
 200     continue
c     

         
c     
c     
 900     continue
         return
         end


