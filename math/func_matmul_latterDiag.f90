
function matmul_latterDiag(matrixFormer,matrixLatter,dim1,dim2) result(outArray)

implicit none

integer dim1,dim2  ! array dimension
real(kind=8),dimension(dim1,dim2),intent(in)  :: matrixFormer
real(kind=8),dimension(dim2)     ,intent(in)  :: matrixLatter
real(kind=8),dimension(dim1,dim2) :: outArray
real(kind=8),dimension(dim2)      :: diagonal
integer id1,id2
!================================================

do id2 = 1,dim2
    diagonal(id2) = matrixLatter(id2)
enddo

do id2 = 1,dim2
    outArray(:,id2) = diagonal(id2) * matrixFormer(:,id2)
enddo
!================================================
return
stop
end function matmul_latterDiag
