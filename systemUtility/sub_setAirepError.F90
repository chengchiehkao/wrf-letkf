
!include 'sub_interp1d.f90'

subroutine setAirepError(airep)

use derivedType
use basicUtility

implicit none

type(obsParent),intent(inout)   :: airep

integer,parameter :: windErrorProfileLevelNum = 2
real(kind=8),dimension(windErrorProfileLevelNum) :: windErrorProfile_value    = (/3.6d0,3.6d0/)
real(kind=8),dimension(windErrorProfileLevelNum) :: windErrorProfile_pressure = (/100000.d0,20000.d0/)

integer,parameter :: tErrorProfileLevelNum = 2
real(kind=8),dimension(tErrorProfileLevelNum) :: tErrorProfile_value    = (/1.0d0,1.0d0/)
real(kind=8),dimension(tErrorProfileLevelNum) :: tErrorProfile_pressure = (/100000d0,20000d0/)

real(kind=8) :: dummyArg_airepError(1)
real(kind=8),parameter :: dummyInvalidValue = -888888.d0
real(kind=8),parameter :: thresholdFactor = 5.d0

integer io  ! loop counter
!================================================


#ifndef PGI
!$omp parallel do default(private) shared(airep,windErrorProfile_value,windErrorProfile_pressure,tErrorProfile_value,tErrorProfile_pressure) schedule(dynamic,100)
#endif
do io = 1 , airep%obsNum

    if ( airep%obs(io)%available ) then

        select case ( trim(adjustl(airep%obs(io)%varName)) )
        case ( 'U' )
            if ( airep%obs(io)%z .gt. windErrorProfile_pressure(1) .or. &
                 airep%obs(io)%z .lt. windErrorProfile_pressure(windErrorProfileLevelNum) ) then
                airep%obs(io)%available = .false.
                cycle
            endif
            call interp1d( dlog(windErrorProfile_pressure(windErrorProfileLevelNum:1:-1)) , windErrorProfile_value(windErrorProfileLevelNum:1:-1) , windErrorProfileLevelNum , &
                           dlog((/airep%obs(io)%z/)) , dummyArg_airepError(1:1) , 1 , &
                           2 , .false. , dummyInvalidValue )
            airep%obs(io)%error = dummyArg_airepError(1)
        case ( 'V' )
            if ( airep%obs(io)%z .gt. windErrorProfile_pressure(1) .or. &
                 airep%obs(io)%z .lt. windErrorProfile_pressure(windErrorProfileLevelNum) ) then
                airep%obs(io)%available = .false.
                cycle
            endif
            call interp1d( dlog(windErrorProfile_pressure(windErrorProfileLevelNum:1:-1)) , windErrorProfile_value(windErrorProfileLevelNum:1:-1) , windErrorProfileLevelNum , &
                           dlog((/airep%obs(io)%z/)) , dummyArg_airepError(1:1) , 1 , &
                           2 , .false. , dummyInvalidValue )
            airep%obs(io)%error = dummyArg_airepError(1)
        case ( 'T' )
            if ( airep%obs(io)%z .gt. tErrorProfile_pressure(1) .or. &
                 airep%obs(io)%z .lt. tErrorProfile_pressure(tErrorProfileLevelNum) ) then
                airep%obs(io)%available = .false.
                cycle
            endif
            call interp1d( dlog(tErrorProfile_pressure(tErrorProfileLevelNum:1:-1)) , tErrorProfile_value(tErrorProfileLevelNum:1:-1) , tErrorProfileLevelNum , &
                           dlog((/airep%obs(io)%z/)) , dummyArg_airepError(1:1) , 1 , &
                           2 , .false. , dummyInvalidValue )
            airep%obs(io)%error = dummyArg_airepError(1)
        case ( 'QVAPOR' )
            airep%obs(io)%available = .false.  ! No QVAPOR for AIREP.
            cycle
        end select


        if ( dabs(airep%obs(io)%innov) .gt. thresholdFactor * airep%obs(io)%error ) then
            airep%obs(io)%available = .false.
            deallocate( airep%obs(io)%background )
        endif

    endif

enddo
#ifndef PGI
!$omp end parallel do
#endif


!================================================
return
stop
end subroutine setAirepError

