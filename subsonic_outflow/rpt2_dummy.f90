! =====================================================
subroutine rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)
! =====================================================
  implicit none !implicit double precision (a-h,o-z)

  integer, intent(in) :: ixy, imp, maxm, meqn, mwaves, maux, mbc, mx
  double precision, dimension(meqn,1-mbc:maxm+mbc), intent(in) :: ql, qr
  double precision, dimension(maux,1-mbc:maxm+mbc), intent(in) :: aux1, aux2, aux3
  double precision, dimension(meqn,1-mbc:maxm+mbc), intent(out) :: asdq, bmasdq, bpasdq

! # Dummy transverse Riemann solver, for use in dimensionally-split algorithm.

  write(*,*) 'Error: Dummy transverse Riemann solver called!'
  return
end subroutine rpt2
