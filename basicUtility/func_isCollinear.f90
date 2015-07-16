
logical function isCollinear(xEndpoint,yEndpoint,xTarget,yTarget)
!====Description of this function====
! Purpose:
!     Determine whether the target point is collinear with the endpoints of the line segment.
!     Return .true. if it is, otherwise .false..
! Input:
!     xEndpoint(rank-1, default real, input only), size of 2.
!     yEndpoint(rank-1, default real, input only), size of 2.
!     xTarget(rank-0, default real, input only).
!     yTarget(rank-0, default real, input only).
! output:
!     isCollinear(rank-0, logical), as function itself.

! Wrote by Cheng-Chieh Kao. (2013-12-23)
!====End of the description====

implicit none
real,intent(in),dimension(2) :: xEndpoint,yEndpoint  ! Endpoints of the line segment.
real,intent(in)              :: xTarget  ,yTarget    ! target point.
!=================================================

! Check whether 3 points are lying on the same ray.
if ( (yEndpoint(2)-yEndpoint(1))*(xTarget-xEndpoint(1)) .eq. &
     (xEndpoint(2)-xEndpoint(1))*(yTarget-yEndpoint(1)) ) then
    ! Check whether the target point is lying between the endpoints.
    if ( (yEndpoint(2)-yTarget)*(yTarget-yEndpoint(1)) .ge. 0. ) then
        if ( (xEndpoint(2)-xTarget)*(xTarget-xEndpoint(1)) .ge. 0. ) then
            isCollinear = .true.
            return
        endif
    endif
endif

isCollinear = .false.
!=================================================
return
stop
end function isCollinear
