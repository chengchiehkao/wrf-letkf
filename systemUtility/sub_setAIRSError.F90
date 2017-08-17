
!include 'sub_interp1d.f90'

subroutine setAIRSError(airs)

use derivedType
use basicUtility

implicit none

type(obsParent),intent(inout)   :: airs
integer,parameter :: tErrorProfileLevelNum = 3
real(kind=8),dimension(tErrorProfileLevelNum) :: tErrorProfile_value    = (/2.0d0,2.0d0,4.0d0/)
real(kind=8),dimension(tErrorProfileLevelNum) :: tErrorProfile_pressure = (/100000.d0,75000.d0,1000.d0/)

integer,parameter :: qvErrorProfileLevelNum = 2
real(kind=8),dimension(qvErrorProfileLevelNum) :: qvErrorProfile_value    = (/0.0025d0,0.00001d0/)
real(kind=8),dimension(qvErrorProfileLevelNum) :: qvErrorProfile_pressure = (/100000.d0,10000.d0/)

real(kind=8) :: dummyArg_airsError(1)
real(kind=8),parameter :: dummyInvalidValue = -888888.d0
real(kind=8),parameter :: thresholdFactor = 5.d0

integer io  ! loop counter
!================================================


#ifndef PGI
!$omp parallel do default(private) shared(airs,tErrorProfile_value,tErrorProfile_pressure,qvErrorProfile_value,qvErrorProfile_pressure) schedule(dynamic,100)
#endif
do io = 1 , airs%obsNum

    if ( airs%obs(io)%available ) then

        select case ( trim(adjustl(airs%obs(io)%varName)) )
        case ( 'U' )
            airs%obs(io)%available = .false.  ! No retrieval of U for AIRS.
            cycle
        case ( 'V' )
            airs%obs(io)%available = .false.  ! No retrieval of V for AIRS.
            cycle
        case ( 'T' )
            if ( airs%obs(io)%z .gt. tErrorProfile_pressure(1) .or. &
                 airs%obs(io)%z .lt. tErrorProfile_pressure(tErrorProfileLevelNum) ) then
                airs%obs(io)%available = .false.
                cycle
            endif
            call interp1d( dlog(tErrorProfile_pressure(tErrorProfileLevelNum:1:-1)) , tErrorProfile_value(tErrorProfileLevelNum:1:-1) , tErrorProfileLevelNum , &
                           dlog((/airs%obs(io)%z/)) , dummyArg_airsError(1:1) , 1 , &
                           2 , .false. , dummyInvalidValue )
            airs%obs(io)%error = dummyArg_airsError(1)
        case ( 'QVAPOR' )
            if ( airs%obs(io)%z .gt. qvErrorProfile_pressure(1) .or. &
                 airs%obs(io)%z .lt. qvErrorProfile_pressure(qvErrorProfileLevelNum) ) then
                airs%obs(io)%available = .false.
                cycle
            endif
            call interp1d( dlog(qvErrorProfile_pressure(qvErrorProfileLevelNum:1:-1)) , qvErrorProfile_value(qvErrorProfileLevelNum:1:-1) , qvErrorProfileLevelNum , &
                           dlog((/airs%obs(io)%z/)) , dummyArg_airsError(1:1) , 1 , &
                           2 , .false. , dummyInvalidValue )
            airs%obs(io)%error = dummyArg_airsError(1)
        end select


        if ( dabs(airs%obs(io)%innov) .gt. thresholdFactor * airs%obs(io)%error ) then
            airs%obs(io)%available = .false.
            deallocate( airs%obs(io)%background )
        endif

    endif

enddo
#ifndef PGI
!$omp end parallel do
#endif


!================================================
return
stop
end subroutine setAIRSError

