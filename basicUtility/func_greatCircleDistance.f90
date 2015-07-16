
real function greatCircleDistance(point1,point2)

implicit none
real,intent(in) :: point1(2),point2(2)  ! 1 = longitude ; 2 = latitude
real point1InRadian(2),point2InRadian(2)
real,parameter :: Re = 6371000. ! unit: meter
real,parameter :: radian = 57.29577951308232
real temp
!================================================
point1InRadian(:) = point1(:)/radian
point2InRadian(:) = point2(:)/radian

temp = sin(point1InRadian(2))*sin(point2InRadian(2)) + &
       cos(point1InRadian(2))*cos(point2InRadian(2))*cos(abs(point2InRadian(1)-point1InRadian(1)))

if ( abs(temp) .gt. 1. ) temp = sign(temp,1.)

greatCircleDistance = Re * acos(temp)

!================================================
return
stop
end function greatCircleDistance
