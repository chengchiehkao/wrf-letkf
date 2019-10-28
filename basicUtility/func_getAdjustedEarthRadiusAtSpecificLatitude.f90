
real function getAdjustedEarthRadiusAtSpecificLatitude(latitude)

! Comments preserved:
! "
! Ad Hoc Radius
! Smithsonian Meteorological Tables p.218 Equation 6
! Ginned up earth radius to compensate for centrifugal force
! variation with latitide
! Note that this is not the Earth ellipsoid radius!
! "

implicit none
real, intent(in) :: latitude  ! unit: degree
real, external :: getGravityAtSpecificLatitude , getVerticalGravityGradientAtSpecificLatitude
!================================================
getAdjustedEarthRadiusAtSpecificLatitude = -2. * getGravityAtSpecificLatitude(latitude) / getVerticalGravityGradientAtSpecificLatitude(latitude)
!================================================
return
end function getAdjustedEarthRadiusAtSpecificLatitude
