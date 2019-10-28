
real function convertTAndPAndQvToRefractivity( t_k , p , qv )

!====Description of this function====
! Purpose:
!     Get refractivity from temperature(t_k) and pressure(p) and water vapor mixing ratio.
! Input:
!     t_k(rank-0, default real, input only).
!     p(rank-0, default real, input only).
!     qv(rank-0, default real, input only).
! output:
!     convertTAndPAndQvToRefractivity(function itself)(rank-0, default real, output only).
! Wrote by Cheng-Chieh Kao. (2018-03-19)
!====End of the description====

implicit none

real,intent(in) :: t_k  ! normal temperature (unit: K)
real,intent(in) :: p    ! pressure (unit: Pa)
real,intent(in) :: qv   ! water vapor mixing ratio (unit: kg/kg)

real :: pw  ! vapor pressure (unit: Pa)
real,parameter :: PaPerhPa=100.
!================================================

pw = p * qv/(0.622+qv)

convertTAndPAndQvToRefractivity = 77.6 * (p/PaPerhPa)/t_k + 3.73e+5* (pw/PaPerhPa)/(t_k**2.)

!================================================
return
end function convertTAndPAndQvToRefractivity
