
!-----Quick Sort (recursive Ver.)

recursive subroutine quickSortWithIndex(a,n,arrayIndex)

implicit none
integer,intent(in) :: n  ! array length
real,intent(inout) :: a(n)
integer,intent(inout) :: arrayIndex(n)
real aMin,aMax  ! input & output array, minimum/maximum of a
real,allocatable :: aTemp(:)
integer i  ! array index
real axis
integer UNum,LNum
integer axisIndex
integer,allocatable :: arrayIndexTemp(:)
!================================================
aMin=minval(a)
aMax=maxval(a)
if ((aMax-aMin).eq.0.) return

allocate(aTemp(n),arrayIndexTemp(n))

axis=a(1)
axisIndex=arrayIndex(1)

LNum=0
UNum=n+1
do i=2,n
  if (a(i).lt.axis) then
    LNum=LNum+1
    aTemp(LNum)=a(i)
    arrayIndexTemp(LNum)=arrayIndex(i)
  else
    UNum=UNum-1
    aTemp(UNum)=a(i)
    arrayIndexTemp(UNum)=arrayIndex(i)
  endif
enddo

if (LNum.ne.(n-1)) then
  UNum=n-(LNum+1)
else
  UNum=0
endif
aTemp(LNum+1)=axis
arrayIndexTemp(LNum+1)=axisIndex

a=aTemp
arrayIndex=arrayIndexTemp

deallocate(aTemp,arrayIndexTemp)

if (LNum.gt.1) call quickSortWithIndex(a(1:LNum),LNum,arrayIndex(1:LNum))
if (UNum.gt.1) call quickSortWithIndex(a(LNum+2:n),UNum,arrayIndex(LNum+2:n))

!================================================
return
stop
end subroutine quickSortWithIndex
