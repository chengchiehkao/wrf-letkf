
real function getVerticalGravityGradientAtSpecificLatitude(latitude)

! Comments preserved:
! "
! Vertical Variation in Normal Gravity
!
! Rate of change of gravity with altitude
! From the Smithsonian Meteorological Tables p.218 Equation 7
! "

implicit none
real, intent(in) :: latitude  ! unit: degree
real             :: cos2L, cos4L
real,parameter :: radian = 57.29577951308232
!================================================
cos2L = cos(2. * radian * latitude)
cos4L = cos(4. * radian * latitude)
getVerticalGravityGradientAtSpecificlatitude = -( 3.085462 * 10.**(-6.) + 2.27*10.**(-9.)*cos2L - 2*10.**(-12.)*cos4L )
!================================================
return
end function getVerticalGravityGradientAtSpecificLatitude
