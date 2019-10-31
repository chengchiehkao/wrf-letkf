
!include 'sub_interp1d.f90'

subroutine setCYGNSSError(cygnss)

use derivedType
use basicUtility

implicit none

type(obsParent),intent(inout)   :: cygnss
!real(kind=8),parameter :: thresholdFactor = 3.d0

integer io  ! loop counter
!================================================

#ifndef PGI
!$omp parallel do default(private) shared(cygnss) schedule(dynamic,100)
#endif
do io = 1 , cygnss%obsNum

    if ( cygnss%obs(io)%available ) then

        cygnss%obs(io)%error = cygnss%obs(io+1)%value  ! Assumed latter field is error.

        if ( cygnss%obs(io)%value .le. 10. .and. dabs(cygnss%obs(io)%innov) .gt. 1. * cygnss%obs(io)%error ) then
            cygnss%obs(io)%available = .false.
            deallocate( cygnss%obs(io)%background )
        endif

        if ( cygnss%obs(io)%value .gt. 10. .and. dabs(cygnss%obs(io)%innov) .gt. 3. * cygnss%obs(io)%error ) then
            cygnss%obs(io)%available = .false.
            deallocate( cygnss%obs(io)%background )
        endif

    endif

enddo
#ifndef PGI
!$omp end parallel do
#endif


!================================================
return
stop
end subroutine setCYGNSSError

