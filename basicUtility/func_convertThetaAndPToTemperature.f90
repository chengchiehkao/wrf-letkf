
real function convertThetaAndPToTemperature(theta,p)

!====Description of this function====
! Purpose:
!     Get normal temperature from potential temperature(theta) and pressure(p).
! Input:
!     theta(rank-0, default real, input only).
!     p(rank-0, default real, input only).
! output:
!     convertThetaAndPToTemperature(function itself)(rank-0, default real, output only).
! Wrote by Cheng-Chieh Kao. (2014-08-11)
!====End of the description====

implicit none

real,intent(in) :: theta  ! potential temperature (unit: K)
real,intent(in) :: p      ! pressure (unit: Pa)

real,parameter :: p0 = 100000. ! standard pressure (unit: Pa)
real,parameter :: R  = 287.058  ! specific gas constant for dry air (unit: J*kg^-1*K^-1)
real,parameter :: cp = 1003.5  ! isobaric mass heat capacity (unit: J*kg^-1*K^-1)
real,parameter :: kappa = R/cp
!================================================

convertThetaAndPToTemperature = theta * (p/p0)**kappa

!================================================
return
stop
end function convertThetaAndPToTemperature
