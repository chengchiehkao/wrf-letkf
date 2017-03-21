
!include 'sub_interp1d.f90'

subroutine setASCATError(ascat)

use derivedType
use basicUtility

implicit none

type(obsParent),intent(inout)   :: ascat
real(kind=8),parameter :: thresholdFactor = 5.d0

integer io  ! loop counter
!================================================


!$omp parallel do default(private) shared(ascat) schedule(dynamic,100)
do io = 1 , ascat%obsNum

    if ( ascat%obs(io)%available ) then

        if ( ascat%obs(io)%value .lt. 12.d0 ) then
            ascat%obs(io)%error = 3.8d0
        elseif ( ascat%obs(io)%value .ge. 12.d0 .and. ascat%obs(io)%value .le. 18.d0 ) then
            ascat%obs(io)%error = 2.1d0
        elseif ( ascat%obs(io)%value .gt. 18.d0 ) then
            ascat%obs(io)%error = 8.7d0
        endif

        if ( dabs(ascat%obs(io)%innov) .gt. thresholdFactor * ascat%obs(io)%error ) then
            ascat%obs(io)%available = .false.
            deallocate( ascat%obs(io)%background )
        endif

    endif

enddo
!$omp end parallel do


!================================================
return
stop
end subroutine setASCATError

