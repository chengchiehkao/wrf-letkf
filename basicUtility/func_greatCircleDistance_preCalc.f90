real function greatCircleDistance_preCalc( point1_lon , point2_lon , &
                                       sin_point1_lat , sin_point2_lat , &
                                       cos_point1_lat , cos_point2_lat )

implicit none
real,intent(in) :: point1_lon,point2_lon
real,intent(in) :: sin_point1_lat,sin_point2_lat
real,intent(in) :: cos_point1_lat,cos_point2_lat
real,parameter :: Re = 6371000. ! unit: meter
!real,parameter :: radian = 57.29577951308232
real temp
!================================================

temp = sin_point1_lat*sin_point2_lat + &
       cos_point1_lat*cos_point2_lat*cos(abs(point1_lon-point2_lon))

if ( abs(temp) .gt. 1. ) temp = sign(temp,1.)

greatCircleDistance_preCalc = Re * acos(temp)

!================================================
return
stop
end function greatCircleDistance_preCalc

