
!  include 'sub_locateAsIndex1d.f90'

subroutine interp1d( xRef    , yRef    , refSize   , &
                     xTarget , yTarget , targetSize, &
                     method , extrap , invalidValue )
!====Description of this subroutine====
! Purpose:
!     Interpolate target y on target x by the curve consist of referenced x & y.
! Input:
!     refSize(rank-0, default integer, input only).
!     xRef(rank-1, default real, input only), shall be strictly increasing.
!     yRef(rank-1, default real, input only).
!     targetSize(rank-0, default integer, input only).
!     xTarget(rank-1, default real, input only).
!     method(rank-0, default integer, input only).
!     extrap(rank-0, logical, input only).
!     invalidValue(rank-0, default real, input only).
! output:
!     yRef(rank-1, default real, output only).

! Wrote by Cheng-Chieh Kao. (2013-12-23)
!====End of the description====

implicit none

integer,intent(in)                    :: refSize    ! size of referenced values.
real,intent(in),dimension(refSize)    :: xRef,yRef  ! referenced x & y.

integer,intent(in)                     :: targetSize  ! size of target values.
real,intent(in) ,dimension(targetSize) :: xTarget     ! target x.
real,intent(out),dimension(targetSize) :: yTarget     ! target y.

integer,intent(in) :: method
logical,intent(in) :: extrap
real,intent(in)    :: invalidValue

integer indOfRef  ! index of interval of referenced data.
integer it  ! loop counter
!=================================================

select case (method)
  case(1)  ! 1=nearest

      if ( .not. extrap) then
          do it =1,targetSize
              call locateAsIndex1d(xRef(:),refSize,xTarget(it),indOfRef)
              if ( indOfRef .lt. 1 .or. indOfRef .gt. refSize ) then
                  yTarget(it) = invalidValue
                  cycle
              endif
              if ( abs(xTarget(it)-xRef(indOfRef)) .lt. abs(xTarget(it)-xRef(indOfRef+1)) ) then
                  yTarget(it) = yRef(indOfRef)
              else
                  yTarget(it) = yRef(indOfRef+1)
              endif
          enddo
      else
          do it =1,targetSize
              call locateAsIndex1d(xRef(:),refSize,xTarget(it),indOfRef)
              if ( indOfRef .lt. 1 ) then
                  indOfRef = 1
              elseif ( indOfRef .gt. refSize ) then
                  indOfRef = refSize-1
              endif

              if ( abs(xTarget(it)-xRef(indOfRef)) .lt. abs(xTarget(it)-xRef(indOfRef+1)) ) then
                  yTarget(it) = yRef(indOfRef)
              else
                  yTarget(it) = yRef(indOfRef+1)
              endif
          enddo
      endif

  case(2)  ! 2=linear

      if ( .not. extrap) then
          do it =1,targetSize
              call locateAsIndex1d(xRef(:),refSize,xTarget(it),indOfRef)
              if ( indOfRef .lt. 1 .or. indOfRef .gt. refSize ) then
                  yTarget(it) = invalidValue
                  cycle
              endif
              yTarget(it) = yRef(indOfRef) + &
                           (xTarget(it)-xRef(indOfRef)) * &
                           (yRef(indOfRef+1)-yRef(indOfRef)) / &
                           (xRef(indOfRef+1)-xRef(indOfRef))
          enddo
      else
          do it =1,targetSize
              call locateAsIndex1d(xRef(:),refSize,xTarget(it),indOfRef)
              if ( indOfRef .lt. 1 ) then
                  indOfRef = 1
              elseif ( indOfRef .gt. refSize ) then
                  indOfRef = refSize-1
              endif
              yTarget(it) = yRef(indOfRef) + &
                           (xTarget(it)-xRef(indOfRef)) * &
                           (yRef(indOfRef+1)-yRef(indOfRef)) / &
                           (xRef(indOfRef+1)-xRef(indOfRef))
          enddo
      endif

end select


!=================================================
return
stop
end subroutine interp1d
