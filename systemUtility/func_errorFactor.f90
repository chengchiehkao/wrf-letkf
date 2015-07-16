
real(kind=8) function errorFactor(rc,rlev,dh,dz)

implicit none

real(kind=8) rc,rlev
real(kind=8) dh,dz
!================================================

errorFactor = dsqrt( dexp( 0.5d0*( (dh/rc)**2.d0 + (dz/rlev)**2.d0 ) ) )

!================================================
return
end function errorFactor
