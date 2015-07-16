
subroutine locateAsIndex1d(valueRef,valueRefSize,valueTarget,indexReturned)
!====Description of this subroutine====
! Purpose:
!     Input a target value and return the index of which interval of 
!     referenced values it lies.
!     Basically, valueRef(indexReturned) <= valueTarget < valueRef(indexReturned+1),
!     but indexReturned=valueRefSize-1  when  valueTarget=valueRef(valueRefSize).
!     This subroutine is base on binary search algorithm.
! Input:
!     valueRef(rank-1, default real, input only), should be strictly increasing.
!     valueRefSize(rank-0, default integer, input only), should be greater than or equal to 1.
!     valueTarget(rank-0, default real, input only).
! output:
!     indexReturned(rank-0, default integer).
!     indexReturned=0 when the target value is smaller than all referenced values.
!     indexReturned=valuRefLen+1 when the target value is greater than all referenced values.

! Wrote by Cheng-Chieh Kao. (2013-12-18)
!====End of the description====

implicit none
integer,intent(in)                      :: valueRefSize  ! size of referenced values.
real,intent(in),dimension(valueRefSize) :: valueRef      ! referenced value array.
real,intent(in)                         :: valueTarget   ! target value.
integer,intent(out)                     :: indexReturned

integer                                 :: searchBlockLowerBound
integer                                 :: searchBlockUpperBound
integer                                 :: searchBlockMedian
!=================================================

! Return specific index if target value is less than or greater than all referenced values.
if ( valueTarget .lt. valueRef(1) ) then
    indexReturned = 0
    return
endif
if ( valueTarget .gt. valueRef(valueRefSize) ) then
    indexReturned = valueRefSize+1
    return
endif


searchBlockLowerBound = 1
searchBlockUpperBound = valueRefSize

! Locate the target value and return the index of interval which it lies.
do  ! End looping when index is located.
    searchBlockMedian = (searchBlockLowerBound+searchBlockUpperBound)/2
    if ( valueTarget .ge. valueRef(searchBlockMedian) ) then
        searchBlockLowerBound = searchBlockMedian
    else
        searchBlockUpperBound = searchBlockMedian
    endif

    if ( (searchBlockUpperBound-searchBlockLowerBound).eq.1 ) then
        indexReturned = searchBlockLowerBound
        return
    endif

enddo

!=================================================
return
stop
end subroutine locateAsIndex1d
