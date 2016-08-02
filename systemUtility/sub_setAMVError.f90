
!include 'sub_interp1d.f90'

subroutine setAMVError(amv)

use derivedType
use basicUtility

implicit none

type(obsParent),intent(inout)   :: amv
real(kind=8),parameter :: thresholdFactor = 5.d0

integer io  ! loop counter
!================================================


!$omp parallel do default(private) shared(amv) schedule(dynamic,100)
do io = 1 , amv%obsNum

    if ( amv%obs(io)%available ) then

        if ( amv%obs(io)%z .ge. 70000.d0 ) then
            amv%obs(io)%error = 4.d0
        elseif ( amv%obs(io)%z .lt. 70000.d0 .and. amv%obs(io)%z .ge. 40000.d0 ) then
            amv%obs(io)%error = 5.d0
        else
            amv%obs(io)%error = 7.d0
        endif

        if ( dabs(amv%obs(io)%innov) .gt. thresholdFactor * amv%obs(io)%error ) then
            amv%obs(io)%available = .false.
            deallocate( amv%obs(io)%background )
        endif

    endif

enddo
!$omp end parallel do


!================================================
return
stop
end subroutine setAMVError

