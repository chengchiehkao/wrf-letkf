    
real function getGravityAtSpecificLatitude(latitude)

! Comments preserved:
! "
! Normal Gravity
! Normal Gravity vs latitude from Smithsonian Meteorological Tables (page 488)
! "

implicit none
real,intent(in) :: latitude  ! unit: degree
real            :: cos2L
real,parameter :: radian = 57.29577951308232
!================================================
cos2L = cos( 2. * latitude * radian )
getGravityAtSpecificlatitude = 9.80616 * (1.00000 - 0.0026373*cos2L + 0.0000059*(cos2L ** 2.) )
!================================================
return
end function getGravityAtSpecificLatitude

