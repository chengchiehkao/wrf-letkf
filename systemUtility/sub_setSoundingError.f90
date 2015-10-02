
!include 'sub_interp1d.f90'

subroutine setSoundingError(sounding)

use derivedType
use basicUtility

implicit none

type(obsParent),intent(inout)   :: sounding
integer,parameter :: windErrorProfileLevelNum = 28
real(kind=8),dimension(windErrorProfileLevelNum) :: windErrorProfile_value   =(/1.1d0,1.1d0,1.1d0,1.1d0,1.1d0,1.1d0,1.1d0, &
                                                                                1.1d0,1.1d0,1.3d0,1.4d0,1.6d0,1.8d0,2.0d0, &
                                                                                2.3d0,2.5d0,2.8d0,3.0d0,3.3d0,3.3d0,3.3d0, &
                                                                                3.0d0,2.7d0,2.7d0,2.1d0,2.7d0,2.7d0,2.7d0/)

real(kind=8),dimension(windErrorProfileLevelNum) :: windErrorProfile_pressure = (/115000.d0,111000.d0,110000.d0,105000.d0,100000.d0,95000.d0,90000.d0, &
                                                                                  85000.d0, 80000.d0, 75000.d0, 70000.d0, 65000.d0,60000.d0,55000.d0, &
                                                                                  50000.d0, 45000.d0, 40000.d0, 35000.d0, 30000.d0,25000.d0,20000.d0, &
                                                                                  15000.d0, 10000.d0,  5000.d0,  4000.d0,  3000.d0, 2000.d0, 1000.d0/)

real(kind=8) :: dummyArg_soundingError(1)
real(kind=8),parameter :: dummyInvalidValue = -888888.d0
real(kind=8),parameter :: thresholdFactor = 5.d0

integer io  ! loop counter
!================================================


!$omp parallel do default(private) shared(sounding,windErrorProfile_value,windErrorProfile_pressure) schedule(dynamic,100)
do io = 1 , sounding%obsNum

    if ( sounding%obs(io)%available ) then

        if ( sounding%obs(io)%z .gt. windErrorProfile_pressure(1) .or. &
             sounding%obs(io)%z .lt. windErrorProfile_pressure(windErrorProfileLevelNum) ) then
            sounding%obs(io)%available = .false.
            cycle
        endif

        select case ( trim(adjustl(sounding%obs(io)%varName)) )
        case ( 'U' )
            call interp1d( windErrorProfile_pressure(windErrorProfileLevelNum:1:-1) , windErrorProfile_value(windErrorProfileLevelNum:1:-1) , windErrorProfileLevelNum , &
                           (/sounding%obs(io)%z/) , dummyArg_soundingError(1:1) , 1 , &
                           2 , .false. , dummyInvalidValue )
            sounding%obs(io)%error = dummyArg_soundingError(1)
        case ( 'V' )
            call interp1d( windErrorProfile_pressure(windErrorProfileLevelNum:1:-1) , windErrorProfile_value(windErrorProfileLevelNum:1:-1) , windErrorProfileLevelNum , &
                           (/sounding%obs(io)%z/) , dummyArg_soundingError(1:1) , 1 , &
                           2 , .false. , dummyInvalidValue )
            sounding%obs(io)%error = dummyArg_soundingError(1)
        case ( 'T' )
            sounding%obs(io)%error = 1.d0
        case ( 'QVAPOR' )
            sounding%obs(io)%error = 2d-4
        case ( 'PSFC' )
            sounding%obs(io)%available = .false.  ! Disable PSFC by default before figure out how to assimilate it.
        end select


        if ( dabs(sounding%obs(io)%innov) .gt. thresholdFactor * sounding%obs(io)%error ) then
            sounding%obs(io)%available = .false.
            deallocate( sounding%obs(io)%background )
        endif

    endif

enddo
!$omp end parallel do


!================================================
return
stop
end subroutine setSoundingError

