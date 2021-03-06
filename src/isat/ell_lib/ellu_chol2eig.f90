!---------------------------------------------------------------------!
! OWNER: Ithaca Combustion Enterprise, LLC                            !
! COPYRIGHT: © 2012, Ithaca Combustion Enterprise, LLC                !
! LICENSE: BSD 3-Clause License (The complete text of the license can !
!  be found in the `LICENSE-ICE.txt' file included in the ISAT-CK7    !
!  source directory.)                                                 !
!---------------------------------------------------------------------!

subroutine ellu_chol2eig( n, g, u, lam )

!  g contains the lower Cholesky factor G of the
!  n x n PSD matrix A = G * G^T.  The eigen-decomposition of A is: 
!  A = u * lam^2 * u^T.  This routine returns u and lam (which are U
!  and S in the SVD of G = U S V^T).

implicit none
integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n
real(k_dp), intent(in)  :: g(n,n)
real(k_dp), intent(out) :: u(n,n), lam(n)

real(k_dp) :: vt(n,n), work(10*n*n+20*n)
integer    :: lwork, info

lwork = 10*n*n+20*n

! copy g into vt

vt = g

call dgesvd( 'A', 'N', n, n, vt, n, lam, u, n, vt, n, work, lwork, info)

if( info /= 0 ) then
   write(0,*)'ellu_chol2eig: info, lwork, work(1) = ', info, lwork, work(1)
   stop
endif

return
end subroutine ellu_chol2eig
