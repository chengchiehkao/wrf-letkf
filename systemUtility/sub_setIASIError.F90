
!include 'sub_interp1d.f90'

subroutine setIASIError(iasi)

use derivedType
use basicUtility

implicit none

type(obsParent),intent(inout)   :: iasi
integer,parameter :: tErrorProfileLevelNum = 14
real(kind=8),dimension(tErrorProfileLevelNum) :: tErrorProfile_value    = (/1.8d0,0.4d0,0.6d0,0.3d0,0.5d0,0.3d0,0.5d0,0.2d0,0.4d0,1.6d0,1.1d0,1.2d0,0.6d0,0.9d0/)
real(kind=8),dimension(tErrorProfileLevelNum) :: tErrorProfile_pressure = (/98000.d0,90000.d0,85000.d0,75000.d0,55000.d0,50000.d0,37500.d0,35000.d0,30000.d0,25000.d0,20000.d0,15000.d0,12500.d0,10000.d0/)

integer,parameter :: rhErrorProfileLevelNum = 19
real(kind=8),dimension(rhErrorProfileLevelNum) :: rhErrorProfile_value    = (/11.1d0,11.1d0,11.2d0,9.5d0,7.d0,7.3d0,1.5d0,4.5d0,5.d0,4.d0,4.d0,5.d0,8.5d0,6.d0,5.7d0,8.d0,2.2d0,2.2d0,1.d0/)
real(kind=8),dimension(rhErrorProfileLevelNum) :: rhErrorProfile_pressure = (/98000d0,90000d0,82500d0,76000d0,65000d0,59000d0,55000d0,50000d0,45000d0,42000d0,38000d0,34000d0,30000d0,27000d0,23500d0,20000d0,16000d0,13000d0,10000d0/)

real(kind=8) :: dummyArg_iasiError(1)
real(kind=8),parameter :: dummyInvalidValue = -888888.d0
real(kind=8),parameter :: thresholdFactor = 5.d0

integer io  ! loop counter
!================================================


#ifndef PGI
!$omp parallel do default(private) shared(iasi,tErrorProfile_value,tErrorProfile_pressure,rhErrorProfile_value,rhErrorProfile_pressure) schedule(dynamic,100)
#endif
do io = 1 , iasi%obsNum

    if ( iasi%obs(io)%available ) then

        select case ( trim(adjustl(iasi%obs(io)%varName)) )
        case ( 'T' )
            if ( iasi%obs(io)%z .gt. tErrorProfile_pressure(1) .or. &
                 iasi%obs(io)%z .lt. tErrorProfile_pressure(tErrorProfileLevelNum) ) then
                iasi%obs(io)%available = .false.
                cycle
            endif
            call interp1d( dlog(tErrorProfile_pressure(tErrorProfileLevelNum:1:-1)) , tErrorProfile_value(tErrorProfileLevelNum:1:-1) , tErrorProfileLevelNum , &
                           dlog((/iasi%obs(io)%z/)) , dummyArg_iasiError(1:1) , 1 , &
                           2 , .false. , dummyInvalidValue )
            iasi%obs(io)%error = dummyArg_iasiError(1)
        case ( 'RH' )
            if ( iasi%obs(io)%z .gt. rhErrorProfile_pressure(1) .or. &
                 iasi%obs(io)%z .lt. rhErrorProfile_pressure(rhErrorProfileLevelNum) ) then
                iasi%obs(io)%available = .false.
                cycle
            endif
            call interp1d( dlog(rhErrorProfile_pressure(rhErrorProfileLevelNum:1:-1)) , rhErrorProfile_value(rhErrorProfileLevelNum:1:-1) , rhErrorProfileLevelNum , &
                           dlog((/iasi%obs(io)%z/)) , dummyArg_iasiError(1:1) , 1 , &
                           2 , .false. , dummyInvalidValue )
            iasi%obs(io)%error = dummyArg_iasiError(1)
        case ( 'QVAPOR' )
            if ( iasi%obs(io)%z .gt. 85000.d0 .or. &
                 iasi%obs(io)%z .lt. rhErrorProfile_pressure(rhErrorProfileLevelNum) ) then
                iasi%obs(io)%available = .false.
                cycle
            endif
            if ( iasi%obs(io)%value .eq. 0.d0 ) then
                iasi%obs(io)%error = 1d-5
                cycle
            endif
            call interp1d( dlog(rhErrorProfile_pressure(rhErrorProfileLevelNum:1:-1)) , rhErrorProfile_value(rhErrorProfileLevelNum:1:-1) , rhErrorProfileLevelNum , &
                           dlog((/iasi%obs(io)%z/)) , dummyArg_iasiError(1:1) , 1 , &
                           2 , .false. , dummyInvalidValue )
            iasi%obs(io)%error = dummyArg_iasiError(1)
            iasi%obs(io)%error = iasi%obs(io)%value * dabs( iasi%obs(io)%error / convertTAndPAndQvToRH( iasi%obs(io-1)%value , iasi%obs(io)%z , iasi%obs(io)%value ) )  ! Assumed iasi%obs(io-1) is T.
        end select


        if ( dabs(iasi%obs(io)%innov) .gt. thresholdFactor * iasi%obs(io)%error ) then
            iasi%obs(io)%available = .false.
            deallocate( iasi%obs(io)%background )
        endif

    endif

enddo
#ifndef PGI
!$omp end parallel do
#endif


!================================================
return
stop
end subroutine setIASIError

