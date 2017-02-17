
!  include 'func_isCollinear.f90'
!  include 'sub_interp1d.f90'

logical function inPolygon(xVertex,yVertex,vertexSize,xTarget,yTarget)
!====Description of this subroutine====
! Purpose:
!     Given a set of vertices of a simple polygon and
!     determine whether a target point is inside the polygon.
!     This function is base on "Ray Method".
! Input:
!     xVertex(rank-1, default real, input only).
!     yVertex(rank-1, default real, input only).
!     vertexSize(rank-0, default integer, input only), should be greater than or equal to 1.
!     xTarget(rank-0, default real, input only).
!     yTarget(rank-0, default real, input only).
! output:
!     inPolygon(rank-0, logical), as function itself.

! Wrote by Cheng-Chieh Kao. (2014-01-03)
!====End of the description====

implicit none
integer,intent(in)                    :: vertexSize     ! size of vertices.
real,intent(in),dimension(vertexSize) :: xVertex,yVertex  ! x-value and y-value of vertices.
real,intent(in)                       :: xTarget,yTarget  ! x-value and y-value of target point.

real,dimension(1)                     :: yTargetDummy,xProjected
integer                               :: intersectionCounter

real,allocatable,dimension(:)         :: yDiffOfTargetAndVertex
integer                               :: iv      ! loop counter
real,parameter                        :: invalidValue=-9.e16  ! for interp1d, meaningless in inPolygon.
logical,external                      :: isCollinear
!=================================================

intersectionCounter = 0

!================
!=====STEP 1=====
!================

! check if the target point is one of vertices of the polygon.
!do iv=1,vertexSize
!    if ( xTarget .eq. xVertex(iv) ) then
!        if ( yTarget .eq. yVertex(iv) ) then
!            inPolygon = .true.
!            return
!        endif
!    endif
!enddo

! check if the target point is on one of the sides of the polygon (Except head-tail side).
do iv=1,vertexSize-1
    if ( isCollinear( xVertex(iv:iv+1) , yVertex(iv:iv+1) , xTarget , yTarget ) ) then
        inPolygon = .true.
        return
    endif
enddo

! check if the target point is on the side between end points,
! basically the same as the previous section.
! Isolated this part to save memory and also prevent from illegal access of memory address.
if ( isCollinear( (/xVertex(1),xVertex(vertexSize)/) , (/yVertex(1),yVertex(vertexSize)/) , xTarget , yTarget ) ) then
    inPolygon = .true.
    return
endif

!================
!=====STEP 2=====
!================

! check if the eastward ray from the target point cross the sides.
do iv=1,vertexSize-1
    if ( xVertex(iv) .gt. xTarget .or. xVertex(iv+1) .gt. xTarget ) then
        if ( (yVertex(iv)-yTarget)*(yVertex(iv+1)-yTarget) .lt. 0. ) then
            yTargetDummy(1:1) = yTarget
            if ( yVertex(iv+1) .gt. yVertex(iv) ) then
                call interp1d( yVertex(iv:iv+1)  , xVertex(iv:iv+1) , 2 , &
                               yTargetDummy(1:1) , xProjected(1:1)  , 1 , &
                               2 , .false. , invalidValue )
            else
                call interp1d( (/yVertex(iv+1),yVertex(iv)/)  , (/xVertex(iv+1),xVertex(iv)/) , 2 , &
                               yTargetDummy(1:1) , xProjected(1:1)  , 1 , &
                               2 , .false. , invalidValue )
            endif

            if ( xProjected(1) .gt. xTarget ) then
                intersectionCounter = intersectionCounter + 1
            endif
        endif
    endif
enddo

! check if the eastward ray from the target point cross the sides between endpoints,
! basically the same as the previous section.
! Isolated this part to save memory and also prevent from illegal access of memory address.
if ( xVertex(vertexSize) .gt. xTarget .or. xVertex(1) .gt. xTarget ) then
    if ( (yVertex(vertexSize)-yTarget)*(yVertex(1)-yTarget) .lt. 0. ) then
        yTargetDummy(1:1) = yTarget
        if ( yVertex(vertexSize) .gt. yVertex(1) ) then
            call interp1d( (/yVertex(1),yVertex(vertexSize)/) , (/xVertex(1),xVertex(vertexSize)/) , 2 , &
                           yTargetDummy(1:1)         , xProjected(1:1)           , 1 , &
                           2 , .false. , invalidValue )
        else
            call interp1d( (/yVertex(vertexSize),yVertex(1)/) , (/xVertex(vertexSize),xVertex(1)/) , 2 , &
                           yTargetDummy(1:1)         , xProjected(1:1)           , 1 , &
                           2 , .false. , invalidValue )
        endif

        if ( xProjected(1) .gt. xTarget ) then
            intersectionCounter = intersectionCounter + 1
        endif
    endif
endif

!================
!=====STEP 3=====
!================

! check if the eastward ray from the target point tangentially cross the vertiecs.

allocate(yDiffOfTargetAndVertex(vertexSize))
yDiffOfTargetAndVertex(:) = yVertex(:) - yTarget

! check if the eastward ray from the target point tangentially cross the indexically 1st vertex.
if ( yDiffOfTargetAndVertex(1) .eq. 0. .and. xTarget.lt.xVertex(1) ) then
    if ( yDiffOfTargetAndVertex(vertexSize)*yDiffOfTargetAndVertex(2) .gt. 0. ) then
        intersectionCounter = intersectionCounter + 2
    elseif ( yDiffOfTargetAndVertex(vertexSize)*yDiffOfTargetAndVertex(2) .lt. 0. ) then
        intersectionCounter = intersectionCounter + 1
    else
        if ( yDiffOfTargetAndVertex(vertexSize)+yDiffOfTargetAndVertex(2) .gt. 0. ) then
            intersectionCounter = intersectionCounter + 2
        elseif ( yDiffOfTargetAndVertex(vertexSize)+yDiffOfTargetAndVertex(2) .lt. 0. ) then
            intersectionCounter = intersectionCounter - 1
        endif
    endif
endif

! check if the eastward ray from the target point tangentially cross the vertices between
!  2nd to (vertexSize-1)th vertex.
do iv=2,vertexSize-1
    if ( yDiffOfTargetAndVertex(iv) .eq. 0. .and. xTarget.lt.xVertex(iv) ) then
        if ( yDiffOfTargetAndVertex(iv-1)*yDiffOfTargetAndVertex(iv+1) .gt. 0. ) then
            intersectionCounter = intersectionCounter + 2
        elseif ( yDiffOfTargetAndVertex(iv-1)*yDiffOfTargetAndVertex(iv+1) .lt. 0. ) then
            intersectionCounter = intersectionCounter + 1
        else
            if ( yDiffOfTargetAndVertex(iv-1)+yDiffOfTargetAndVertex(iv+1) .gt. 0. ) then
                intersectionCounter = intersectionCounter + 2
            elseif ( yDiffOfTargetAndVertex(iv-1)+yDiffOfTargetAndVertex(iv+1) .lt. 0. ) then
                intersectionCounter = intersectionCounter - 1
            endif
        endif
    endif
enddo

! check if the eastward ray from the target point tangentially cross the indexically last vertex.
if ( yDiffOfTargetAndVertex(vertexSize) .eq. 0. .and. xTarget.lt.xVertex(vertexSize) ) then
    if ( yDiffOfTargetAndVertex(vertexSize-1)*yDiffOfTargetAndVertex(1) .gt. 0. ) then
        intersectionCounter = intersectionCounter + 2
    elseif ( yDiffOfTargetAndVertex(vertexSize-1)*yDiffOfTargetAndVertex(1) .lt. 0. ) then
        intersectionCounter = intersectionCounter + 1
    else
        if ( yDiffOfTargetAndVertex(vertexSize-1)+yDiffOfTargetAndVertex(1) .gt. 0. ) then
            intersectionCounter = intersectionCounter + 2
        elseif ( yDiffOfTargetAndVertex(vertexSize-1)+yDiffOfTargetAndVertex(1) .lt. 0. ) then
            intersectionCounter = intersectionCounter - 1
        endif
    endif
endif

deallocate(yDiffOfTargetAndVertex)


if ( mod(abs(intersectionCounter),2) .eq. 1 ) then
    inPolygon = .true.
else
    inPolygon = .false.
endif
!=================================================
return
stop
end function inPolygon
