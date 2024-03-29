
function matSqrt(inArray,k) result(outArray)

implicit none
!================================================
!==declaration for matSqrt itself================
integer,parameter :: realKind=8
integer k  ! i/o array dimension
real(realKind),dimension(k,k),intent(in) :: inArray
real(realKind),dimension(k,k)            :: outArray
real(realKind),dimension(k)              :: S
integer i  ! loop index
interface
    function matmul_latterDiag(matrixFormer,matrixLatter,dim1,dim2) result(outArray)
      implicit none
      integer dim1,dim2  ! array dimension
      real(kind=8),dimension(dim1,dim2),intent(in)  :: matrixFormer
      real(kind=8),dimension(dim2)     ,intent(in)  :: matrixLatter
      real(kind=8),dimension(dim1,dim2) :: outArray
    end function matmul_latterDiag

    function matmul_ext(A,B,m,k,n) result(C)
      implicit none
      integer m,k,n
      real(kind=8) :: A(m,k),B(k,n),C(m,n)
    end function matmul_ext
endinterface
!================================================
!==declaration for subprogram of LAPACK==========
character(len=1) JOBZ,RANGE,UPLO
integer N,LDA
real(8),dimension(k,k) :: A
real(8) :: VL,VU
integer :: IL,IU
real(8) :: ABSTOL
integer :: M
real(8),dimension(k)   :: W
real(8),dimension(k,k) :: Z
integer :: LDZ
integer :: ISUPPZ(2*k)
real(8),allocatable :: WORK(:)
integer :: LWORK
integer,allocatable :: IWORK(:)
integer :: LIWORK
integer :: INFO
!================================================

JOBZ='V'
RANGE='A'
UPLO='U'
N=k
A(:,:)=inArray(:,:)
LDA=k
VL=-1d30
VU=+1d30
IL=1
IU=k
ABSTOL=0d0
LDZ=k
LWORK=3*k-1
LIWORK=LWORK

allocate(WORK(LWORK))
allocate(IWORK(LWORK))
INFO=-1

call dsyevr( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, -1, IWORK, -1, INFO )

if ( INFO.eq.0 ) then
    LWORK=nint(WORK(1))
    LIWORK=IWORK(1)
    deallocate(WORK,IWORK)
    allocate(WORK(LWORK))
    allocate(IWORK(LWORK))
    A(:,:)=inArray(:,:)
    call dsyevr( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO )
else
    print*,'func_matSqrt failed, INFO=',INFO
    print*,'program stopped.'
    stop
endif
deallocate(WORK)
deallocate(IWORK)

do i=1,k
    S(i)=dsqrt(W(i))
enddo

outArray = matmul_ext(matmul_latterDiag(Z,S,k,k),transpose(Z),k,k,k)

!================================================
return
stop
end function matSqrt
