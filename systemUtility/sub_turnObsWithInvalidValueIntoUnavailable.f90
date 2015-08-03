
subroutine turnObsWithInvalidValueIntoUnavailable(obs)

use derivedType
implicit none
type(obsParent),intent(inout) :: obs
real(kind=8),parameter :: invalidValue = -888888.d0
!================================================

where ( dabs(obs%obs(:)%lon - invalidValue) .lt. dabs(invalidValue*epsilon(1.d0)) .or. &
        dabs(obs%obs(:)%lat - invalidValue) .lt. dabs(invalidValue*epsilon(1.d0)) .or. &
        dabs(obs%obs(:)%z   - invalidValue) .lt. dabs(invalidValue*epsilon(1.d0)) .or. &
        dabs(obs%obs(:)%value - invalidValue) .lt. dabs(invalidValue*epsilon(1.d0)) )

    obs%obs(:)%available = .false.

endwhere

!================================================
return
stop
end subroutine turnObsWithInvalidValueIntoUnavailable
