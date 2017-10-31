
!include 'sub_locateAsIndex2d.f90'
!include 'sub_interp1d.f90'
!include 'sub_interp2d.f90'

subroutine convertBackgroundToWindSat(background,ensembleSize,domain,domain_mean,windsat)

use derivedType
use basicUtility

implicit none

integer,intent(in)              :: ensembleSize
type(backgroundInfo),intent(in) :: background(ensembleSize)
type(domainInfo),intent(in)     :: domain(ensembleSize)
type(domainInfo),intent(in)     :: domain_mean
type(obsParent),intent(inout)   :: windsat

integer obsVarIndexRankOne , obsVarIndexRankTwo
integer obsZIndexRankOne   , obsZIndexRankTwo

real(kind=8),dimension(1) :: obsVarBuffer , obsZBuffer  ! SHALL AWARE OF DIFFERENCES BETWEEN P & GPH

real(kind=8),parameter :: invalidValue = -9.d6

integer io,iens,iz  ! loop counter
!================================================

#ifndef PGI
!$omp parallel do default(none) &
!$omp private(io,iens,iz,obsVarIndexRankOne,obsVarIndexRankTwo,obsZIndexRankOne,obsZIndexRankTwo,obsVarBuffer,obsZBuffer) &
!$omp shared(windsat,domain,domain_mean,background,ensembleSize) &
!$omp schedule(dynamic,100)
#endif
do io = 1 , windsat%obsNum


    if ( windsat%obs(io)%available ) then

        allocate( windsat%obs(io)%background(ensembleSize) )
        windsat%obs(io)%background(:) = 0.d0

        select case ( trim(adjustl(windsat%obs(io)%varName)) )
        case ( 'U10' )
            call locateAsIndex2d( domain(1)%lon(:,:)        , domain(1)%lat(:,:) , &
                                  domain(1)%size_westToEast , domain(1)%size_southToNorth , &
                                  windsat%obs(io)%lon , windsat%obs(io)%lat , &
                                  obsVarIndexRankOne , obsVarIndexRankTwo )
            if ( obsVarIndexRankOne.eq.0 .or. obsVarIndexRankTwo.eq.0 ) then
                deallocate( windsat%obs(io)%background )
                windsat%obs(io)%available = .false.
                cycle
            endif
            obsZIndexRankOne = obsVarIndexRankOne
            obsZIndexRankTwo = obsVarIndexRankTwo
        case ( 'V10' )
            call locateAsIndex2d( domain(1)%lon(:,:)        , domain(1)%lat(:,:) , &
                                  domain(1)%size_westToEast , domain(1)%size_southToNorth , &
                                  windsat%obs(io)%lon , windsat%obs(io)%lat , &
                                  obsVarIndexRankOne , obsVarIndexRankTwo )
            if ( obsVarIndexRankOne.eq.0 .or. obsVarIndexRankTwo.eq.0 ) then
                deallocate( windsat%obs(io)%background )
                windsat%obs(io)%available = .false.
                cycle
            endif
            obsZIndexRankOne = obsVarIndexRankOne
            obsZIndexRankTwo = obsVarIndexRankTwo
        end select


        do iens = 1 , ensembleSize

            select case ( trim(adjustl(windsat%obs(io)%varName)) )
            case ( 'U10' )
                call interp2d( domain(1)%lon(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                               domain(1)%lat(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                               background(iens)%u10(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                               2 , 2 , &
                               (/windsat%obs(io)%lon/) , (/windsat%obs(io)%lat/) , obsVarBuffer(1:1) , 1, &
                               1 , invalidValue )
                call interp2d( domain_mean%lon(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                               domain_mean%lat(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                               domain_mean%pressure_w(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1,1) , &
                               2 , 2 , &
                               (/windsat%obs(io)%lon/) , (/windsat%obs(io)%lat/) , obsZBuffer(1:1) , 1, &
                               1 , invalidValue )
            case ( 'V10' )
                call interp2d( domain(1)%lon(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                               domain(1)%lat(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                               background(iens)%v10(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                               2 , 2 , &
                               (/windsat%obs(io)%lon/) , (/windsat%obs(io)%lat/) , obsVarBuffer(1:1) , 1, &
                               1 , invalidValue )
                call interp2d( domain_mean%lon(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                               domain_mean%lat(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                               domain_mean%pressure_w(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1,1) , &
                               2 , 2 , &
                               (/windsat%obs(io)%lon/) , (/windsat%obs(io)%lat/) , obsZBuffer(1:1) , 1, &
                               1 , invalidValue )
            end select


            select case ( trim(adjustl(windsat%obs(io)%zName)) )
            case ( 'P' )
                windsat%obs(io)%z = obsZBuffer(1) - 100.d0  ! Assumed 10-m high means 1-hPa less than PSFC of mean domain.
                windsat%obs(io)%background(iens) = obsVarBuffer(1)
            end select

        enddo

        windsat%obs(io)%innov = windsat%obs(io)%value - sum(windsat%obs(io)%background(:))/real(ensembleSize,8)

    endif

enddo
#ifndef PGI
!$omp end parallel do
#endif



!================================================
return
stop
end subroutine convertBackgroundToWindSat

