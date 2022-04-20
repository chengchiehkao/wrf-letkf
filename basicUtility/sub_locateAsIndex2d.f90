
!  include 'func_inPolygon.f90'

subroutine locateAsIndex2d(xRef,yRef,refRankOneSize,refRankTwoSize, &
                           xTarget,yTarget, &
                           rankOneIndexReturned,rankTwoIndexReturned)
!====Description of this subroutine====
! Purpose:
!     Given a set of rank-2 grid points and return the index of 1st rank and 2nd rank of the grids
!     which the target point lies.
!     The target point is in the grid of
!     xRef(rankOneIndexReturned:rankOneIndexReturned+1,rankTwoIndexReturned:rankTwoIndexReturned+1) &
!     yRef(rankOneIndexReturned:rankOneIndexReturned+1,rankTwoIndexReturned:rankTwoIndexReturned+1).
!     but indexReturned=valueRefSize-1  when  valueTarget=valueRef(valueRefSize).
!     This subroutine is base on binary search algorithm in 2 way(1st & 2nd rank).
! Input:
!     xRef(rank-2, default real, input only).
!     yRef(rank-2, default real, input only).
!     refRankOneSize(rank-0, default integer, input only), should be greater than or equal to 1.
!     refRankTwoSize(rank-0, default integer, input only), should be greater than or equal to 1.
!     xTarget(rank-0, default real, input only).
!     yTarget(rank-0, default real, input only).
! output:
!     rankOneIndexReturned(rank-0, default integer, output only).
!     rankTwoIndexReturned(rank-0, default integer, output inly).
!     Return zeros means the target point does not lie in the grids.

! Wrote by Cheng-Chieh Kao. (2014-01-03)
! Bug Fixed @ 2014-08-15: (1) Exception handling for minimum side length of 2.
!                         (2) Indices of vertor subscript shall not be constant.

!====End of the description====

implicit none
integer,intent(in)                      :: refRankOneSize,refRankTwoSize  ! size of 1st and 2nd rank of the grids.
real,intent(in),dimension(refRankOneSize,refRankTwoSize) :: xRef,yRef     ! x- and y-value of the grids.
real,intent(in)                         :: xTarget,yTarget  ! x- and y-value of the target point.
integer,intent(out)                     :: rankOneIndexReturned,rankTwoIndexReturned

integer                                 :: blockRankOneLowerBound
integer                                 :: blockRankOneUpperBound
integer                                 :: blockRankOneMedian
integer                                 :: blockRankTwoLowerBound
integer                                 :: blockRankTwoUpperBound
integer                                 :: blockRankTwoMedian

integer                                 :: ir  ! loop counter

logical,external                        :: inPolygon
!=================================================

! Return zeros if the target point does not lie in the grids.
if ( .not. inPolygon( (/xRef(1:refRankOneSize-1 ,1                  ), &
                        xRef(refRankOneSize     ,1:refRankTwoSize-1 ), &
                        xRef(refRankOneSize:2:-1,refRankTwoSize     ), &
                        xRef(1                  ,refRankTwoSize:2:-1)/) , &
                      (/yRef(1:refRankOneSize-1 ,1                  ), &
                        yRef(refRankOneSize     ,1:refRankTwoSize-1 ), &
                        yRef(refRankOneSize:2:-1,refRankTwoSize     ), &
                        yRef(1                  ,refRankTwoSize:2:-1)/) , &
                        2*(refRankOneSize+refRankTwoSize)-4, &
                        xTarget,yTarget) ) then
    rankOneIndexReturned = 0
    rankTwoIndexReturned = 0
    return
else
    if ( refRankOneSize.eq.2 .and. refRankTwoSize.eq.2 ) then
        rankOneIndexReturned = 1
        rankTwoIndexReturned = 1
        return
    endif
endif



blockRankOneLowerBound = 1
blockRankOneUpperBound = refRankOneSize

blockRankTwoLowerBound = 1
blockRankTwoUpperBound = refRankTwoSize

! Locate the target point and return the index of grids which it lies.
do  ! End looping when indices were retrived.

    if ( (blockRankOneUpperBound-blockRankOneLowerBound).gt.1 ) then
        blockRankOneMedian = (blockRankOneLowerBound+blockRankOneUpperBound)/2
    else
        blockRankOneMedian = blockRankOneUpperBound
    endif

    if ( size((/xRef(blockRankOneLowerBound:blockRankOneMedian-1   ,blockRankTwoLowerBound                            ), &
                      xRef(blockRankOneMedian                            ,blockRankTwoLowerBound:blockRankTwoUpperBound-1   ), &
                      xRef(blockRankOneMedian:blockRankOneLowerBound+1:-1,blockRankTwoUpperBound                            ), &
                      xRef(blockRankOneLowerBound                        ,blockRankTwoUpperBound:blockRankTwoLowerBound+1:-1)/)) .ne. &
         2*((blockRankOneMedian-blockRankOneLowerBound+1)+(blockRankTwoUpperBound-blockRankTwoLowerBound+1))-4 ) then
         print*,'Got an inequality in locateAsIndex2d.'
         stop
    endif

    if ( inPolygon( (/xRef(blockRankOneLowerBound:blockRankOneMedian-1   ,blockRankTwoLowerBound                            ), &
                      xRef(blockRankOneMedian                            ,blockRankTwoLowerBound:blockRankTwoUpperBound-1   ), &
                      xRef(blockRankOneMedian:blockRankOneLowerBound+1:-1,blockRankTwoUpperBound                            ), &
                      xRef(blockRankOneLowerBound                        ,blockRankTwoUpperBound:blockRankTwoLowerBound+1:-1)/) , &
                    (/yRef(blockRankOneLowerBound:blockRankOneMedian-1   ,blockRankTwoLowerBound                            ), &
                      yRef(blockRankOneMedian                            ,blockRankTwoLowerBound:blockRankTwoUpperBound-1   ), &
                      yRef(blockRankOneMedian:blockRankOneLowerBound+1:-1,blockRankTwoUpperBound                            ), &
                      yRef(blockRankOneLowerBound                        ,blockRankTwoUpperBound:blockRankTwoLowerBound+1:-1)/) , &
                      2*((blockRankOneMedian-blockRankOneLowerBound+1)+(blockRankTwoUpperBound-blockRankTwoLowerBound+1))-4, &
                      xTarget,yTarget) ) then

        if ( (blockRankOneUpperBound-blockRankOneLowerBound).gt.1 ) then
            blockRankOneUpperBound = blockRankOneMedian
        endif
    else
        if ( (blockRankOneUpperBound-blockRankOneLowerBound).gt.1 ) then
            blockRankOneLowerBound = blockRankOneMedian
        endif
    endif

    if ( (blockRankOneUpperBound-blockRankOneLowerBound).eq.1 .and. &
         (blockRankTwoUpperBound-blockRankTwoLowerBound).eq.1 ) then
        rankOneIndexReturned = blockRankOneLowerBound
        rankTwoIndexReturned = blockRankTwoLowerBound
        return
    endif

!=====

    if ( (blockRankTwoUpperBound-blockRankTwoLowerBound).gt.1 )  then
        blockRankTwoMedian = (blockRankTwoLowerBound+blockRankTwoUpperBound)/2
    else
        blockRankTwoMedian = blockRankTwoUpperBound
    endif

    if ( size((/xRef(blockRankOneLowerBound:blockRankOneUpperBound-1   ,blockRankTwoLowerBound                        ), &
                      xRef(blockRankOneUpperBound                            ,blockRankTwoLowerBound:blockRankTwoMedian-1   ), &
                      xRef(blockRankOneUpperBound:blockRankOneLowerBound+1:-1,blockRankTwoMedian                            ), &
                      xRef(blockRankOneLowerBound                            ,blockRankTwoMedian:blockRankTwoLowerBound+1:-1)/)) .ne. &
         2*((blockRankOneUpperBound-blockRankOneLowerBound+1)+(blockRankTwoMedian-blockRankTwoLowerBound+1))-4 ) then
         print*,'Got an inequality in locateAsIndex2d.'
         stop
    endif

    if ( inPolygon( (/xRef(blockRankOneLowerBound:blockRankOneUpperBound-1   ,blockRankTwoLowerBound                        ), &
                      xRef(blockRankOneUpperBound                            ,blockRankTwoLowerBound:blockRankTwoMedian-1   ), &
                      xRef(blockRankOneUpperBound:blockRankOneLowerBound+1:-1,blockRankTwoMedian                            ), &
                      xRef(blockRankOneLowerBound                            ,blockRankTwoMedian:blockRankTwoLowerBound+1:-1)/) , &
                    (/yRef(blockRankOneLowerBound:blockRankOneUpperBound-1   ,blockRankTwoLowerBound                        ), &
                      yRef(blockRankOneUpperBound                            ,blockRankTwoLowerBound:blockRankTwoMedian-1   ), &
                      yRef(blockRankOneUpperBound:blockRankOneLowerBound+1:-1,blockRankTwoMedian                            ), &
                      yRef(blockRankOneLowerBound                            ,blockRankTwoMedian:blockRankTwoLowerBound+1:-1)/) , &
                      2*((blockRankOneUpperBound-blockRankOneLowerBound+1)+(blockRankTwoMedian-blockRankTwoLowerBound+1))-4, &
                      xTarget,yTarget) ) then

        if ( (blockRankTwoUpperBound-blockRankTwoLowerBound).gt.1 )  then
            blockRankTwoUpperBound = blockRankTwoMedian
        endif
    else
        if ( (blockRankTwoUpperBound-blockRankTwoLowerBound).gt.1 )  then
            blockRankTwoLowerBound = blockRankTwoMedian
        endif
    endif

    if ( (blockRankOneUpperBound-blockRankOneLowerBound).eq.1 .and. &
         (blockRankTwoUpperBound-blockRankTwoLowerBound).eq.1 ) then
        rankOneIndexReturned = blockRankOneLowerBound
        rankTwoIndexReturned = blockRankTwoLowerBound
        return
    endif

enddo

!=================================================
return
stop
end subroutine locateAsIndex2d
