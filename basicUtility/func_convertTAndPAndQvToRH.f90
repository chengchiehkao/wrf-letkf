
real function convertTAndPAndQvToRH( t_k , p , qv )

!====Description of this function====
! Purpose:
!     Get relative humidity from temperature(t_k) and pressure(p) and water vapor mixing ratio.
! Input:
!     t_k(rank-0, default real, input only).
!     p(rank-0, default real, input only).
!     qv(rank-0, default real, input only).
! output:
!     convertTAndPAndQvToRH(function itself)(rank-0, default real, output only).
! Wrote by Cheng-Chieh Kao. (2016-09-21)
!====End of the description====

implicit none

real,intent(in) :: t_k  ! normal temperature (unit: K)
real,intent(in) :: p    ! pressure (unit: Pa)
real,intent(in) :: qv   ! water vapor mixing ratio (unit: kg/kg)

real :: t_c  ! normal temperature (unit: degree-C)
real :: es   ! saturation vapor pressure in hPa
real :: e    ! vapor pressure in hPa
real :: RH   ! Relative Humidity
!================================================

t_c = t_k-273.15
es  = 6.1094*exp(17.625*t_c/(t_c+243.04))
e   = (p/100.)*qv/(0.622+qv)

RH  = (e/es)*100.

convertTAndPAndQvToRH = RH

!================================================
return
stop
end function convertTAndPAndQvToRH
