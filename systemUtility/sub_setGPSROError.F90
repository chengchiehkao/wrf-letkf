
!include 'sub_interp1d.f90'

subroutine setGPSROError(gpsro)

use derivedType
use basicUtility

implicit none

type(obsParent),intent(inout)   :: gpsro

integer,parameter :: refractivityErrorProfileLevelNum = 2
real(kind=8),dimension(refractivityErrorProfileLevelNum) :: refractivityErrorProfile_value = (/2.d0,2.d0/)
real(kind=8),dimension(refractivityErrorProfileLevelNum) :: refractivityErrorProfile_gmh   = (/0.d0,100000d0/)

real(kind=8) :: dummyArg_gpsroError(1)
real(kind=8),parameter :: dummyInvalidValue = -888888.d0
real(kind=8),parameter :: thresholdFactor = 5.d0

integer io  ! loop counter
!================================================


#ifndef PGI
!$omp parallel do default(private) shared(gpsro,refractivityErrorProfile_value,refractivityErrorProfile_gmh) schedule(dynamic,100)
#endif
do io = 1 , gpsro%obsNum

    if ( gpsro%obs(io)%available ) then

        select case ( trim(adjustl(gpsro%obs(io)%varName)) )
        case ( 'REF' )
            if ( gpsro%obs(io)%z .lt. refractivityErrorProfile_gmh(1) .or. &
                 gpsro%obs(io)%z .gt. refractivityErrorProfile_gmh(refractivityErrorProfileLevelNum) ) then
                gpsro%obs(io)%available = .false.
                cycle
            endif
            call interp1d( refractivityErrorProfile_gmh(:) , refractivityErrorProfile_value(:) , refractivityErrorProfileLevelNum , &
                           (/gpsro%obs(io)%z/) , dummyArg_gpsroError(1:1) , 1 , &
                           2 , .false. , dummyInvalidValue )
            gpsro%obs(io)%error = dummyArg_gpsroError(1)
        end select


        if ( dabs(gpsro%obs(io)%innov) .gt. thresholdFactor * gpsro%obs(io)%error ) then
            gpsro%obs(io)%available = .false.
            deallocate( gpsro%obs(io)%background )
        endif

    endif

enddo
#ifndef PGI
!$omp end parallel do
#endif


!================================================
return
stop
end subroutine setGPSROError

