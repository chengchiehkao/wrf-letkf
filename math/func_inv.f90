
function inv(inArray,n) result(outArray)

implicit none

integer n  ! array dimension
real(kind=8),dimension(n,n),intent(in)  :: inArray
real(kind=8),dimension(n,n) :: outArray
integer  ipiv(n),info,lwork
real(kind=8),allocatable :: work(:)
!================================================

outArray=inArray

call dgetrf(n,n,outArray,n,ipiv,info)
if ( info.lt.0 ) print*,'inv failed, the i-th argument had an illegal value.'
if ( info.gt.0 ) print*,'inv failed, U(i,i) is exactly zero.'

lwork=-1
allocate(work(1))
call dgetri(n,outArray,n,ipiv,work,lwork,info)

lwork=work(1)
deallocate(work)
allocate(work(lwork))
call dgetri(n,outArray,n,ipiv,work,lwork,info)
if ( info.lt.0 ) print*,'inv failed, the i-th argument had an illegal value.'
if ( info.gt.0 ) print*,'inv failed, U(i,i) is exactly zero.'

!if (info.ne.0) print*,'inv failed!'
!================================================
return
stop
end function inv
