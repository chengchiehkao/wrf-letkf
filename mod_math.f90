
include 'math/func_inv.f90'
include 'math/func_matSqrt.f90'
include 'math/func_identityMat.f90'
include 'math/func_inv_diagMat.f90'
include 'math/func_matmul_formerDiag.f90'
include 'math/func_matmul_latterDiag.f90'


module math


interface

    function inv(inArray,n) result(outArray)
      implicit none
      integer n  ! array dimension
      real(kind=8),dimension(n,n),intent(in)  :: inArray
      real(kind=8),dimension(n,n) :: outArray
    end function inv

    function matSqrt(inArray,k) result(outArray)
      implicit none
      integer,parameter :: realKind=8
      integer k  ! i/o array dimension
      real(realKind),dimension(k,k),intent(in) :: inArray
      real(realKind),dimension(k,k)            :: outArray
    end function matSqrt

    function identityMat(n) result(outArray)
      implicit none
      integer n  ! array dimension
      real(kind=8),dimension(n,n) :: outArray
    end function identityMat

    function inv_diagMat(inArray,n) result(outArray)
      implicit none
      integer n  ! array dimension
      real(kind=8),dimension(n,n),intent(in)  :: inArray
      real(kind=8),dimension(n,n) :: outArray
    end function inv_diagMat

    function matmul_formerDiag(matrixFormer,matrixLatter,dim1,dim2) result(outArray)
      implicit none
      integer dim1,dim2  ! array dimension
      real(kind=8),dimension(dim1)     ,intent(in)  :: matrixFormer
      real(kind=8),dimension(dim1,dim2),intent(in)  :: matrixLatter
      real(kind=8),dimension(dim1,dim2) :: outArray
    end function matmul_formerDiag

    function matmul_latterDiag(matrixFormer,matrixLatter,dim1,dim2) result(outArray)
      implicit none
      integer dim1,dim2  ! array dimension
      real(kind=8),dimension(dim1,dim2),intent(in)  :: matrixFormer
      real(kind=8),dimension(dim2)     ,intent(in)  :: matrixLatter
      real(kind=8),dimension(dim1,dim2) :: outArray
    end function matmul_latterDiag

end interface

end module math
