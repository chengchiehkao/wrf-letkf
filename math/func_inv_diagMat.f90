
function inv_diagMat(inArray,n) result(outArray)

implicit none

integer n  ! array dimension
real(kind=8),dimension(n,n),intent(in)  :: inArray
real(kind=8),dimension(n,n) :: outArray
integer iDim
!================================================
outArray(:,:)=0.d0

do iDim = 1,n
    outArray(iDim,iDim) = 1.d0 / inArray(iDim,iDim)
enddo
!================================================
return
stop
end function inv_diagMat
