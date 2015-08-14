
function matmul_formerDiag(matrixFormer,matrixLatter,dim1,dim2) result(outArray)

implicit none

integer dim1,dim2  ! array dimension
real(kind=8),dimension(dim1)     ,intent(in)  :: matrixFormer
real(kind=8),dimension(dim1,dim2),intent(in)  :: matrixLatter
real(kind=8),dimension(dim1,dim2) :: outArray
real(kind=8),dimension(dim1)      :: diagonal
integer id1,id2
!================================================

do id1 = 1,dim1
    diagonal(id1) = matrixFormer(id1)
enddo

do id2 = 1,dim2
    outArray(:,id2) = diagonal(:) * matrixLatter(:,id2)
enddo
!================================================
return
stop
end function matmul_formerDiag
