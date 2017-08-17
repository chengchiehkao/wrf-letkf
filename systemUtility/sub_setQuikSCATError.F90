
!include 'sub_interp1d.f90'

subroutine setQuikSCATError(quikscat)

use derivedType
use basicUtility

implicit none

type(obsParent),intent(inout)   :: quikscat
real(kind=8),parameter :: thresholdFactor = 5.d0

integer io  ! loop counter
!================================================


#ifndef PGI
!$omp parallel do default(private) shared(quikscat) schedule(dynamic,100)
#endif
do io = 1 , quikscat%obsNum

    if ( quikscat%obs(io)%available ) then

        if ( quikscat%obs(io)%value .lt. 10.d0 ) then
            quikscat%obs(io)%error = 3.9d0
        elseif ( quikscat%obs(io)%value .ge. 10.d0 .and. quikscat%obs(io)%value .le. 17.2d0 ) then
            quikscat%obs(io)%error = 3.7d0
        elseif ( quikscat%obs(io)%value .gt. 17.2d0 ) then
            quikscat%obs(io)%error = 4.5d0
        endif

        if ( dabs(quikscat%obs(io)%innov) .gt. thresholdFactor * quikscat%obs(io)%error ) then
            quikscat%obs(io)%available = .false.
            deallocate( quikscat%obs(io)%background )
        endif

    endif

enddo
#ifndef PGI
!$omp end parallel do
#endif


!================================================
return
stop
end subroutine setQuikSCATError

