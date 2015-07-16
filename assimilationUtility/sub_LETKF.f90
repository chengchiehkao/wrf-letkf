

!include 'mod_math.f'

subroutine LETKF( xb_mean , xb_pert , &
                  xa_mean , xa_pert , &
                  yo , yb , R , &
                  m , n , k )

use math
implicit none
integer,parameter :: realKind=8
integer m,n,k  ! m=model state, n=observation number, k=ensemble size
real(realKind),dimension(m,1) :: xb_mean,xa_mean
real(realKind),dimension(m,k) :: xb_pert,xa_pert
real(realKind),dimension(n,1) :: yo 
real(realKind),dimension(n,k) :: yb
real(realKind),dimension(n,n) :: R

real(realKind),dimension(n,1) :: yb_mean
real(realKind),dimension(n,k) :: yb_pert
real(realKind),dimension(m,m) :: pa
real(realKind),dimension(k,k) :: pa_tilde
!real(realKind),dimension(k,k) :: identityMat

integer i  ! loop index
!================================================
!print*,shape(xb_mean(:,:)),shape(xb_pert(:,:)),shape(xa_mean(:,:)),shape(xa_pert(:,:)),shape(yo(:,:)),shape(yb(:,:)),shape(R(:,:)),m,n,k

!identityMat(:,:)=0d0
!do i=1,k
!    identityMat(i,i)=1d0
!enddo
do i=1,n
    yb_mean(i,1)=sum(yb(i,:))/real(k,8)
    yb_pert(i,:)=yb(i,:)-yb_mean(i,1)
enddo

pa_tilde = inv( dble(k-1)*identityMat(k) + &
                matmul( matmul_latterDiag(transpose(yb_pert(:,:)),inv_diagMat(R(:,:),n),k,n),yb_pert(:,:) ) &
                ,k )

! pa = 1d0/dble(k-1) *matmul( xa_pert(:,:) , transpose(xa_pert(:,:)) )

xa_mean = xb_mean + &
          matmul( &
          matmul_latterDiag( &
          matmul( &
          matmul( &
            xb_pert(:,:),pa_tilde(:,:) &
          ), &
            transpose(yb_pert(:,:)) &
          ), &
            inv_diagMat(R(:,:),n) , m , n &
          ), &
            (yo(:,:)-yb_mean(:,:)) &
          )

!read(*,*)
!print*,xb_pert(:,:)
!print*,repeat('=',20)
!print*,shape(xb_mean(:,:)),shape(xb_pert(:,:)),shape(xa_mean(:,:)),shape(xa_pert(:,:)),shape(yo(:,:)),shape(yb(:,:)),shape(R(:,:)),m,n,k
!print*,matmul(xb_pert(:,:),pa_tilde(:,:))
!print*,matmul(xb_pert(1:2,1:2),pa_tilde(1:2,1:2))

xa_pert = matmul( xb_pert(:,:) , matSqrt(dble(k-1)*pa_tilde(:,:),k) )

!================================================
return
stop
end subroutine LETKF
