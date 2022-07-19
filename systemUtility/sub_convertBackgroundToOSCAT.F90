
!include 'sub_locateAsIndex2d.f90'
!include 'sub_interp1d.f90'
!include 'sub_interp2d.f90'

subroutine convertBackgroundToOSCAT(background,ensembleSize,domain,domain_mean,oscat,systemParameters,domainID)

use derivedType
use basicUtility

implicit none

integer,intent(in)              :: ensembleSize
type(backgroundInfo),intent(in) :: background(ensembleSize)
type(domainInfo),intent(in)     :: domain(ensembleSize)
type(domainInfo),intent(in)     :: domain_mean
type(obsParent),intent(inout)   :: oscat
type(systemParameter),intent(in) :: systemParameters
integer,intent(in)              :: domainID

integer obsVarIndexRankOne , obsVarIndexRankTwo
integer obsZIndexRankOne   , obsZIndexRankTwo

real(kind=8),dimension(1) :: obsVarBuffer , obsZBuffer  ! SHALL AWARE OF DIFFERENCES BETWEEN P & GPH
real(kind=8) :: u_rotated , v_rotated
real(kind=8),dimension(1) :: sinalpha_interpolated

real(kind=8),parameter :: invalidValue = -9.d6

integer io,iens,iz  ! loop counter
!================================================


if ( systemParameters%rotateUAndVOfObsBasedOnWRFMapProjection ) then

    do io = 1 , oscat%obsNum

        if ( oscat%obs(io)%available .and. count(oscat%obs(io)%insideHorizontalDomain(:)).eq.domainID ) then

            call locateAsIndex2d( domain(1)%lon(:,:)        , domain(1)%lat(:,:) , &
                                  domain(1)%size_westToEast , domain(1)%size_southToNorth , &
                                  oscat%obs(io)%lon , oscat%obs(io)%lat , &
                                  obsZIndexRankOne , obsZIndexRankTwo )

            select case ( trim(adjustl(oscat%obs(io)%varName)) )
            case ( 'V10' )
                call interp2d( domain(1)%lon(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                               domain(1)%lat(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                               domain(1)%sinalpha(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                               2 , 2 , &
                               (/oscat%obs(io)%lon/) , (/oscat%obs(io)%lat/) , sinalpha_interpolated(1:1) , 1, &
                               1 , invalidValue )

                call rotateUAndVofObsBasedOnSinAlpha( oscat%obs(io-1)%value , oscat%obs(io)%value , sinalpha_interpolated(1) , u_rotated , v_rotated )
                oscat%obs(io-1)%value = u_rotated
                oscat%obs(io  )%value = v_rotated
            end select

        endif

    enddo

endif


#ifndef PGI
!$omp parallel do default(none) &
!$omp private(io,iens,iz,obsVarIndexRankOne,obsVarIndexRankTwo,obsZIndexRankOne,obsZIndexRankTwo,obsVarBuffer,obsZBuffer) &
!$omp shared(oscat,domain,domain_mean,background,ensembleSize,domainID) &
!$omp schedule(dynamic,100)
#endif
do io = 1 , oscat%obsNum


    if ( oscat%obs(io)%available .and. count(oscat%obs(io)%insideHorizontalDomain(:)).eq.domainID ) then

        allocate( oscat%obs(io)%background(ensembleSize) )
        oscat%obs(io)%background(:) = 0.d0

        select case ( trim(adjustl(oscat%obs(io)%varName)) )
        case ( 'U10' )
            call locateAsIndex2d( domain(1)%lon(:,:)        , domain(1)%lat(:,:) , &
                                  domain(1)%size_westToEast , domain(1)%size_southToNorth , &
                                  oscat%obs(io)%lon , oscat%obs(io)%lat , &
                                  obsVarIndexRankOne , obsVarIndexRankTwo )
            if ( obsVarIndexRankOne.eq.0 .or. obsVarIndexRankTwo.eq.0 ) then
                deallocate( oscat%obs(io)%background )
                oscat%obs(io)%available = .false.
                cycle
            endif
            obsZIndexRankOne = obsVarIndexRankOne
            obsZIndexRankTwo = obsVarIndexRankTwo
        case ( 'V10' )
            call locateAsIndex2d( domain(1)%lon(:,:)        , domain(1)%lat(:,:) , &
                                  domain(1)%size_westToEast , domain(1)%size_southToNorth , &
                                  oscat%obs(io)%lon , oscat%obs(io)%lat , &
                                  obsVarIndexRankOne , obsVarIndexRankTwo )
            if ( obsVarIndexRankOne.eq.0 .or. obsVarIndexRankTwo.eq.0 ) then
                deallocate( oscat%obs(io)%background )
                oscat%obs(io)%available = .false.
                cycle
            endif
            obsZIndexRankOne = obsVarIndexRankOne
            obsZIndexRankTwo = obsVarIndexRankTwo
        end select


        do iens = 1 , ensembleSize

            select case ( trim(adjustl(oscat%obs(io)%varName)) )
            case ( 'U10' )
                call interp2d( domain(1)%lon(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                               domain(1)%lat(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                               background(iens)%u10(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                               2 , 2 , &
                               (/oscat%obs(io)%lon/) , (/oscat%obs(io)%lat/) , obsVarBuffer(1:1) , 1, &
                               1 , invalidValue )
                call interp2d( domain_mean%lon(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                               domain_mean%lat(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                               domain_mean%pressure_w(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1,1) , &
                               2 , 2 , &
                               (/oscat%obs(io)%lon/) , (/oscat%obs(io)%lat/) , obsZBuffer(1:1) , 1, &
                               1 , invalidValue )
            case ( 'V10' )
                call interp2d( domain(1)%lon(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                               domain(1)%lat(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                               background(iens)%v10(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                               2 , 2 , &
                               (/oscat%obs(io)%lon/) , (/oscat%obs(io)%lat/) , obsVarBuffer(1:1) , 1, &
                               1 , invalidValue )
                call interp2d( domain_mean%lon(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                               domain_mean%lat(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                               domain_mean%pressure_w(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1,1) , &
                               2 , 2 , &
                               (/oscat%obs(io)%lon/) , (/oscat%obs(io)%lat/) , obsZBuffer(1:1) , 1, &
                               1 , invalidValue )
            end select


            select case ( trim(adjustl(oscat%obs(io)%zName)) )
            case ( 'P' )
                oscat%obs(io)%z = obsZBuffer(1) - 100.d0  ! Assumed 10-m high means 1-hPa less than PSFC of mean domain.
                oscat%obs(io)%background(iens) = obsVarBuffer(1)
            end select

        enddo

        oscat%obs(io)%innov = oscat%obs(io)%value - sum(oscat%obs(io)%background(:))/real(ensembleSize,8)

    endif

enddo
#ifndef PGI
!$omp end parallel do
#endif



!================================================
return
stop
end subroutine convertBackgroundToOSCAT

