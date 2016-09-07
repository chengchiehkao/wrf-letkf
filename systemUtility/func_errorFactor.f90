
real(kind=8) function errorFactor(rc_h,rc_z,dh,dz)

implicit none

real(kind=8) rc_h,rc_z
real(kind=8) dh,dz
!================================================

errorFactor = dsqrt( dexp( 0.5d0*( (dh/rc_h)**2.d0 + (dz/rc_z)**2.d0 ) ) )

!================================================
return
end function errorFactor
