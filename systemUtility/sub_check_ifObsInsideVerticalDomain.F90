
!include 'sub_locateAsIndex2d.f90'
!include 'sub_interp2d.f90'

subroutine check_ifObsInsideVerticalDomain(domain,domainSize,obs)

use derivedType
use basicUtility

implicit none
integer,intent(in) :: domainSize
type(domainInfo),intent(in) :: domain(domainSize)
type(obsParent),intent(inout) :: obs

real(kind=8) zBottom(1),zTop(1)
real(kind=8),parameter :: invalidValue = -9.d6

integer obsIndexRankOne , obsIndexRankTwo

integer io , id  ! loop counter
!================================================


! Shall aware of unit conversion.

#ifndef PGI
!$omp parallel do default(private) shared(domainSize,domain,obs) schedule(dynamic,100)
#endif
do io=1,obs%obsNum

    if ( .not. obs%obs(io)%available ) cycle

    call locateAsIndex2d( domain(1)%lon(:,:)        , domain(1)%lat(:,:) , &
                          domain(1)%size_westToEast , domain(1)%size_southToNorth , &
                          obs%obs(io)%lon , obs%obs(io)%lat , &
                          obsIndexRankOne , obsIndexRankTwo )

    if ( obsIndexRankOne .lt. 1 .or. obsIndexRankTwo .lt. 1 ) cycle

    do id=1,domainSize

        if ( obs%obs(io)%available .and. &
             ( trim(obs%obs(io)%zName) .eq. 'P' .or. trim(obs%obs(io)%zName) .eq. 'GPH' ) ) then

            if ( trim(obs%obs(io)%zName) .eq. 'P' )  then 
                call interp2d( domain(id)%lon(obsIndexRankOne:obsIndexRankOne+1,obsIndexRankTwo:obsIndexRankTwo+1) , &
                               domain(id)%lat(obsIndexRankOne:obsIndexRankOne+1,obsIndexRankTwo:obsIndexRankTwo+1) , &
                               domain(id)%pressure(obsIndexRankOne:obsIndexRankOne+1,obsIndexRankTwo:obsIndexRankTwo+1,1) , &
                               2 , 2 , &
                               obs%obs(io:io)%lon , obs%obs(io:io)%lat , zBottom(1:1) , 1 , &
                               1 , invalidValue )

                call interp2d( domain(id)%lon(obsIndexRankOne:obsIndexRankOne+1,obsIndexRankTwo:obsIndexRankTwo+1) , &
                               domain(id)%lat(obsIndexRankOne:obsIndexRankOne+1,obsIndexRankTwo:obsIndexRankTwo+1) , &
                               domain(id)%pressure(obsIndexRankOne:obsIndexRankOne+1,obsIndexRankTwo:obsIndexRankTwo+1,domain(id)%size_bottomToTop) , &
                               2 , 2 , &
                               obs%obs(io:io)%lon , obs%obs(io:io)%lat , zTop(1:1) , 1 , &
                               1 , invalidValue )

                if ( obs%obs(io)%z > zBottom(1) .or. obs%obs(io)%z < zTop(1) ) then
                    obs%obs(io)%available = .false.
                endif
            elseif ( trim(obs%obs(io)%zName) .eq. 'GPH' )  then 
                call interp2d( domain(id)%lon(obsIndexRankOne:obsIndexRankOne+1,obsIndexRankTwo:obsIndexRankTwo+1) , &
                               domain(id)%lat(obsIndexRankOne:obsIndexRankOne+1,obsIndexRankTwo:obsIndexRankTwo+1) , &
                               domain(id)%GPH(obsIndexRankOne:obsIndexRankOne+1,obsIndexRankTwo:obsIndexRankTwo+1,1) , &
                               2 , 2 , &
                               obs%obs(io:io)%lon , obs%obs(io:io)%lat , zBottom(1:1) , 1 , &
                               1 , invalidValue )

                call interp2d( domain(id)%lon(obsIndexRankOne:obsIndexRankOne+1,obsIndexRankTwo:obsIndexRankTwo+1) , &
                               domain(id)%lat(obsIndexRankOne:obsIndexRankOne+1,obsIndexRankTwo:obsIndexRankTwo+1) , &
                               domain(id)%GPH(obsIndexRankOne:obsIndexRankOne+1,obsIndexRankTwo:obsIndexRankTwo+1,domain(id)%size_bottomToTop_stag) , &
                               2 , 2 , &
                               obs%obs(io:io)%lon , obs%obs(io:io)%lat , zTop(1:1) , 1 , &
                               1 , invalidValue )

                if ( obs%obs(io)%z < zBottom(1) .or. obs%obs(io)%z > zTop(1) ) then
                    obs%obs(io)%available = .false.
                endif
            endif

        endif

    enddo
enddo
#ifndef PGI
!$omp end parallel do
#endif

!================================================
return
stop
end subroutine check_ifObsInsideVerticalDomain
