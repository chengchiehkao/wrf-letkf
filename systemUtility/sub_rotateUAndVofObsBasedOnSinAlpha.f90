
subroutine rotateUAndVofObsBasedOnSinAlpha( u_org , v_org , sinalpha , u_rotated , v_rotated )

implicit none

real(kind=8),intent(in)  :: u_org , v_org , sinalpha
real(kind=8),intent(out) :: u_rotated , v_rotated

real(kind=8) :: cosalpha
!================================================

cosalpha = dsqrt( 1d0 - sinalpha**2.d0 )

u_rotated =  (u_org * cosalpha) + (v_org * sinalpha)
v_rotated =  (v_org * cosalpha) - (u_org * sinalpha)

!================================================
return
end subroutine rotateUAndVofObsBasedOnSinAlpha
