
function identityMat(n) result(outArray)

implicit none

integer n  ! array dimension
real(kind=8),dimension(n,n) :: outArray
integer i  ! loop counter
!================================================
outArray(:,:) = 0.d0

do i = 1,n
    outArray(i,i) = 1.d0
enddo
!================================================
return
stop
end function identityMat
