
!include 'sub_interp1d.f90'

subroutine setWindSatError(windsat)

use derivedType
use basicUtility

implicit none

type(obsParent),intent(inout)   :: windsat
real(kind=8),parameter :: thresholdFactor = 5.d0

integer io  ! loop counter
!================================================


#ifndef PGI
!$omp parallel do default(private) shared(windsat) schedule(dynamic,100)
#endif
do io = 1 , windsat%obsNum

    if ( windsat%obs(io)%available ) then

        windsat%obs(io)%error = 0.91d0  ! Statistics value in tropical pacific ocean for median refquency product. Refers to Zhang Lei et al. Acta Oceanol. Sin., 2016, Vol. 35, No. 1, P. 67-78.

        if ( dabs(windsat%obs(io)%innov) .gt. thresholdFactor * windsat%obs(io)%error ) then
            windsat%obs(io)%available = .false.
            deallocate( windsat%obs(io)%background )
        endif

    endif

enddo
#ifndef PGI
!$omp end parallel do
#endif


!================================================
return
stop
end subroutine setWindSatError

