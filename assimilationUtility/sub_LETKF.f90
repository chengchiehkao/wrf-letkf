

!include 'mod_math.f'

subroutine LETKF( xb_mean , xb_pert , &
                  xa_mean , xa_pert , &
                  yo , yb , R , &
                  m , n , k , inflationFactor )

use math
implicit none
integer,parameter :: realKind=8
integer m,n,k  ! m=model state, n=observation number, k=ensemble size
real(realKind),dimension(m,1) :: xb_mean,xa_mean
real(realKind),dimension(m,k) :: xb_pert,xa_pert
real(realKind),dimension(n,1) :: yo 
real(realKind),dimension(n,k) :: yb
real(realKind),dimension(n)   :: R

real(realKind),dimension(n,1) :: yb_mean
real(realKind),dimension(n,k) :: yb_pert
real(realKind),dimension(m,m) :: pa
real(realKind),dimension(k,k) :: pa_tilde

real(realKind) :: inflationFactor , rho

integer i  ! loop index
!================================================

rho = 1d0/inflationFactor

do i=1,n
    yb_mean(i,1)=sum(yb(i,:))/real(k,8)
    yb_pert(i,:)=yb(i,:)-yb_mean(i,1)
enddo

pa_tilde = inv( dble(k-1)*identityMat(k)*rho + &
                matmul( matmul_latterDiag(transpose(yb_pert(:,:)),1d0/R(:),k,n),yb_pert(:,:) ) &
                ,k )

xa_mean = xb_mean + &
          matmul( &
          matmul_latterDiag( &
          matmul( &
          matmul( &
            xb_pert(:,:),pa_tilde(:,:) &
          ), &
            transpose(yb_pert(:,:)) &
          ), &
            1d0/R(:) , m , n &
          ), &
            (yo(:,:)-yb_mean(:,:)) &
          )

xa_pert = matmul( xb_pert(:,:) , matSqrt(dble(k-1)*pa_tilde(:,:),k) )

!================================================
return
stop
end subroutine LETKF
