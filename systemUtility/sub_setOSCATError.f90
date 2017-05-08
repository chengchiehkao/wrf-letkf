
!include 'sub_interp1d.f90'

subroutine setOSCATError(oscat)

use derivedType
use basicUtility

implicit none

type(obsParent),intent(inout)   :: oscat
real(kind=8),parameter :: thresholdFactor = 5.d0

integer io  ! loop counter
!================================================


!$omp parallel do default(private) shared(oscat) schedule(dynamic,100)
do io = 1 , oscat%obsNum

    if ( oscat%obs(io)%available ) then

        oscat%obs(io)%error = 2.0d0

        if ( dabs(oscat%obs(io)%innov) .gt. thresholdFactor * oscat%obs(io)%error ) then
            oscat%obs(io)%available = .false.
            deallocate( oscat%obs(io)%background )
        endif

    endif

enddo
!$omp end parallel do


!================================================
return
stop
end subroutine setOSCATError

