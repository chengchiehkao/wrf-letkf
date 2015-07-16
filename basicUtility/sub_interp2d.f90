
!  include 'sub_locateAsIndex2d.f90'
!  include 'sub_interp1d.f90'

subroutine interp2d( xRef           , yRef           , zRef , &
                     refRankOneSize , refRankTwoSize , &
                     xTarget , yTarget , zTarget , targetSize, &
                     method , invalidValue )
!====Description of this subroutine====
! Purpose:
!     Interpolate target z on target x & y
!     by the surface consist of referenced x, y & z which z=f(x,y).
! Input:
!     refRankOneSize(rank-0, default integer, input only).
!     refRankTwoSize(rank-0, default integer, input only).
!     xRef(rank-2, default real, input only), shall be strictly increasing along rank 1.
!     yRef(rank-2, default real, input only), shall be strictly increasing along rank 2.
!     targetSize(rank-0, default integer, input only).
!     xTarget(rank-1, default real, input only).
!     yTarget(rank-1, default real, input only).
!     method(rank-0, default integer, input only), 1 = bilinear.
!     invalidValue(rank-0, default real, input only).
! output:
!     zTarget(rank-1, default real, output only).

! Wrote by Cheng-Chieh Kao. (2014-01-10)
! Buf Fixed @ 2014-08-27: Arguments with wrong rank for interp1d.
!====End of the description====

implicit none

integer,intent(in)                                       :: refRankOneSize,refRankTwoSize  ! size of referenced values.
real,intent(in),dimension(refRankOneSize,refRankTwoSize) :: xRef,yRef,zRef  ! referenced x & y.

integer,intent(in)                     :: targetSize  ! size of target values.
real,intent(in) ,dimension(targetSize) :: xTarget,yTarget  ! target x & y.
real,intent(out),dimension(targetSize) :: zTarget  ! target z.

integer,intent(in) :: method  ! 1=bilinear
real,intent(in)    :: invalidValue

real,dimension(1) :: zProjectedOnLowerIndSide,zProjectedOnUpperIndSide
real,dimension(1) :: yProjectedOnLowerIndSide,yProjectedOnUpperIndSide
integer rankOneIndexReturned,rankTwoIndexReturned  ! index of interval of referenced data.

integer it  ! loop counter
!=================================================

select case (method)
  case(1)  ! 1=bilinear

      do it =1,targetSize
          call locateAsIndex2d(xRef(:,:),yRef(:,:),refRankOneSize,refRankTwoSize, &
                               xTarget(it),yTarget(it),rankOneIndexReturned,rankTwoIndexReturned)
          if ( rankOneIndexReturned .eq. 0 .or. rankTwoIndexReturned .eq. 0 ) then
              zTarget(it) = invalidValue
              cycle
          endif

          call interp1d( xRef(rankOneIndexReturned:rankOneIndexReturned+1,rankTwoIndexReturned) , &
                         yRef(rankOneIndexReturned:rankOneIndexReturned+1,rankTwoIndexReturned) , &
                         2 , &
                         xTarget(it:it) , yProjectedOnLowerIndSide(1:1) , 1 , &
                         2 , .true. , invalidValue )
          call interp1d( xRef(rankOneIndexReturned:rankOneIndexReturned+1,rankTwoIndexReturned) , &
                         zRef(rankOneIndexReturned:rankOneIndexReturned+1,rankTwoIndexReturned) , &
                         2 , &
                         xTarget(it:it) , zProjectedOnLowerIndSide(1:1) , 1 , &
                         2 , .true. , invalidValue )

          call interp1d( xRef(rankOneIndexReturned:rankOneIndexReturned+1,rankTwoIndexReturned+1) , &
                         yRef(rankOneIndexReturned:rankOneIndexReturned+1,rankTwoIndexReturned+1) , &
                         2 , &
                         xTarget(it:it) , yProjectedOnUpperIndSide(1:1) , 1 , &
                         2 , .true. , invalidValue )
          call interp1d( xRef(rankOneIndexReturned:rankOneIndexReturned+1,rankTwoIndexReturned+1) , &
                         zRef(rankOneIndexReturned:rankOneIndexReturned+1,rankTwoIndexReturned+1) , &
                         2 , &
                         xTarget(it:it) , zProjectedOnUpperIndSide(1:1) , 1 , &
                         2 , .true. , invalidValue )

          call interp1d( (/yProjectedOnLowerIndSide,yProjectedOnUpperIndSide/) , &
                         (/zProjectedOnLowerIndSide,zProjectedOnUpperIndSide/) , &
                         2 , &
                         yTarget(it:it) , zTarget(it:it) , 1 , &
                         2 , .true. , invalidValue )

      enddo

end select


!=================================================
return
stop
end subroutine interp2d
