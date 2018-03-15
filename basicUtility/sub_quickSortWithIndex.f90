
!-----Quick Sort (recursive Ver.)

recursive subroutine quickSortWithIndex(dataToBeSorted,sizeOfData,IndexOfData)

implicit none

integer,intent(in)    :: sizeOfData  ! size of data array
real,intent(inout)    :: dataToBeSorted(sizeOfData)
integer,intent(inout) :: IndexOfData(sizeOfData)

real                   dataMin,dataMax  ! input & output array, minimum/maximum of a
real,allocatable ::    data_temp(:)
integer                iD  ! loop counter
real                   axis
integer                UpperPartCount,lowerPartCount
integer                indexOfAxis
integer,allocatable :: IndexOfData_temp(:)
!================================================

dataMin = minval(dataToBeSorted)
dataMax = maxval(dataToBeSorted)

if ( dataMax-dataMin .eq. 0. ) return

allocate( data_temp(sizeOfData) , IndexOfData_temp(sizeOfData) )

axis        = dataToBeSorted(1)
indexOfAxis = IndexOfData(1)

lowerPartCount = 0            ! initial value
UpperPartCount = sizeOfData+1 ! initial value
do iD = 2,sizeOfData
    if ( dataToBeSorted(iD) .lt. axis ) then
        lowerPartCount = lowerPartCount+1
        data_temp(lowerPartCount)        = dataToBeSorted(iD)
        IndexOfData_temp(lowerPartCount) = IndexOfData(iD)
    else
        UpperPartCount = UpperPartCount-1
        data_temp(UpperPartCount)        = dataToBeSorted(iD)
        IndexOfData_temp(UpperPartCount) = IndexOfData(iD)
    endif
enddo

if ( lowerPartCount .ne. (sizeOfData-1) ) then
    UpperPartCount = sizeOfData-(lowerPartCount+1)
else
    UpperPartCount = 0
endif
data_temp(lowerPartCount+1)        = axis
IndexOfData_temp(lowerPartCount+1) = indexOfAxis

dataToBeSorted = data_temp
IndexOfData    = IndexOfData_temp

deallocate(data_temp,IndexOfData_temp)

if (lowerPartCount.gt.1) call quickSortWithIndex( dataToBeSorted(1:lowerPartCount) , lowerPartCount , IndexOfData(1:lowerPartCount) )
if (UpperPartCount.gt.1) call quickSortWithIndex( dataToBeSorted(lowerPartCount+2:sizeOfData) , UpperPartCount , IndexOfData(lowerPartCount+2:sizeOfData) )

!================================================
return
stop
end subroutine quickSortWithIndex
