subroutine greatCircleDistance_preCalc_vec( point1_lon , point2_lon , &
                                            sin_point1_lat , sin_point2_lat , &
                                            cos_point1_lat , cos_point2_lat , size1 , distance )

implicit none
integer,intent(in) :: size1
real,intent(in) :: point1_lon(size1),point2_lon
real,intent(in) :: sin_point1_lat(size1),sin_point2_lat
real,intent(in) :: cos_point1_lat(size1),cos_point2_lat
real,parameter :: Re = 6371000. ! unit: meter
!real,parameter :: radian = 57.29577951308232
real temp(size1)
!================================================

temp(:) = sin_point1_lat(:)*sin_point2_lat + &
          cos_point1_lat(:)*cos_point2_lat*cos(abs(point1_lon(:)-point2_lon))

where ( abs(temp) .gt. 1. ) temp = sign(temp,1.)

greatCircleDistance_preCalc = Re * acos(temp)

!================================================
return
stop
end subroutine greatCircleDistance_preCalc_vec

