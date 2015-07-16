
function matSqrt(inArray,k) result(outArray)

implicit none
!================================================
!==declaration for matSqrt itself================
integer,parameter :: realKind=8
integer k  ! i/o array dimension
real(realKind),dimension(k,k),intent(in) :: inArray
real(realKind),dimension(k,k)            :: outArray
real(realKind),dimension(k,k)            :: S
integer i  ! loop index
!real(realKind),external :: matmul_latterDiag
interface
    function matmul_latterDiag(matrixFormer,matrixLatter,dim1,dim2) result(outArray)
      implicit none
      integer dim1,dim2  ! array dimension
      real(kind=8),dimension(dim1,dim2),intent(in)  :: matrixFormer
      real(kind=8),dimension(dim2,dim2),intent(in)  :: matrixLatter
      real(kind=8),dimension(dim1,dim2) :: outArray
    end function matmul_latterDiag
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
!goto 0001  ! shall delete after efficiency test
LWORK=3*k-1
LIWORK=LWORK

allocate(WORK(LWORK))
allocate(IWORK(LWORK))
INFO=-1

call dsyevr( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, &
             ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, -1, &
             IWORK, -1, INFO )
if ( INFO.eq.0 ) then
    LWORK=nint(WORK(1))
    LIWORK=IWORK(1)
    deallocate(WORK,IWORK)
!0001 continue  ! shall delete after efficiency test
    !LWORK=k
    !LIWORK=12
    allocate(WORK(LWORK))
    allocate(IWORK(LWORK))
    A(:,:)=inArray(:,:)
    call dsyevr( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, &
                 ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK, &
                 IWORK, LIWORK, INFO )
    !print*,real(LWORK)/real(k)
    !print*,real(LIWORK)/real(k)
    !read(*,*)
else
    if ( INFO .ne. 0 ) print*,'func_matSqrt failed! INFO=',INFO
endif
deallocate(WORK)
deallocate(IWORK)

S(:,:)=0d0
do i=1,k
    S(i,i)=dsqrt(W(i))
    !if( S(i,i).lt.1d-10 ) then
    !    print*,'!!!Warning for rank deficiency!!!'
    !    S(i,i)=0.0d0
    !    INFO=-1
    !    exit
    !endif
enddo

!outArray = matmul(matmul(Z,S),transpose(Z))
outArray = matmul(matmul_latterDiag(Z,S,k,k),transpose(Z))

!================================================
return
stop
end function matSqrt
