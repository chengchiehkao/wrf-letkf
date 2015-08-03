
subroutine setSoundingError(sounding)

use derivedType

implicit none

type(obsParent),intent(inout)   :: sounding

integer io  ! loop counter
!================================================

!$omp parallel do default(none) &
!$omp private(io) &
!$omp shared(sounding) &
!$omp schedule(dynamic,100)
do io = 1 , sounding%obsNum

    if ( sounding%obs(io)%available ) then  ! SHALL AWARE OF PSFC

        select case ( trim(adjustl(sounding%obs(io)%varName)) )
        !select case (sounding%obs(io)%varName(1:1))
        case ( 'U' )
            sounding%obs(io)%error = 1.d0
        case ( 'V' )
            sounding%obs(io)%error = 1.d0
        case ( 'T' )
            sounding%obs(io)%error = 0.5d0
        case ( 'QVAPOR' )
        !case ( 'Q' )
            sounding%obs(io)%error = 0.2d-3
        end select

    endif

enddo
!$omp end parallel do


!================================================
return
stop
end subroutine setSoundingError

