
real function convertGPHToGMHAtSpecificLatitude(gph,latitude)

! Comments preserved:
! "
! Convert geometric altitude to geopotential height
! Zg ........ Geometric altitude (km)
! Latitude .. Degrees
! Return .... Geopotential height
!
! Verified against Smithsonian Meteorological Tables (page 220-221)
! if 9.80665 set to 9.8 as in tables
! "

implicit none
real, intent(in) :: gph  ! unit: meter
real, intent(in) :: latitude  ! unit: degree
real :: adjustedEarthRadius
real, external :: getAdjustedEarthRadiusAtSpecificLatitude , getGravityAtSpecificLatitude
!================================================
adjustedEarthRadius = getAdjustedEarthRadiusAtSpecificLatitude(latitude)
convertGPHToGMHAtSpecificLatitude = adjustedEarthRadius * (getGravityAtSpecificLatitude(latitude) / 9.80665) * (gph/(gph + adjustedEarthRadius))
!================================================
return
end function convertGPHToGMHAtSpecificLatitude
