
!include 'sub_locateAsIndex2d.f90'
!include 'sub_interp1d.f90'
!include 'sub_interp2d.f90'

subroutine mapGPSROFromGMHToP(domain_mean,gpsro,domainID)

use derivedType
use basicUtility

implicit none

type(domainInfo),intent(in)     :: domain_mean
type(obsParent),intent(inout)   :: gpsro
integer,intent(in) :: domainID

integer obsVarIndexRankOne , obsVarIndexRankTwo
integer obsZIndexRankOne   , obsZIndexRankTwo

real(kind=8),dimension(domain_mean%size_bottomToTop) :: obsVarBuffer , obsZBuffer  ! SHALL AWARE OF DIFFERENCES BETWEEN P & GPH

real(kind=8) :: dummy_output(1)

real(kind=8),parameter :: invalidValue = -9.d6

integer io,iens,iz  ! loop counter
!================================================

#ifndef PGI
!$omp parallel do default(none) &
!$omp private(io,iens,iz,obsVarIndexRankOne,obsVarIndexRankTwo,obsZIndexRankOne,obsZIndexRankTwo,obsVarBuffer,obsZBuffer,dummy_output) &
!$omp shared(gpsro,domain_mean,domainID) &
!$omp schedule(dynamic,100)
#endif
do io = 1 , gpsro%obsNum


    if ( gpsro%obs(io)%available  .and. count(gpsro%obs(io)%insideHorizontalDomain(:)).eq.domainID ) then  ! SHALL AWARE OF PSFC


        call locateAsIndex2d( domain_mean%lon(:,:)        , domain_mean%lat(:,:) , &
                              domain_mean%size_westToEast , domain_mean%size_southToNorth , &
                              gpsro%obs(io)%lon , gpsro%obs(io)%lat , &
                              obsVarIndexRankOne , obsVarIndexRankTwo )
        if ( obsVarIndexRankOne.eq.0 .or. obsVarIndexRankTwo.eq.0 ) then
            deallocate( gpsro%obs(io)%background )
            gpsro%obs(io)%available = .false.
            cycle
        endif
        obsZIndexRankOne = obsVarIndexRankOne
        obsZIndexRankTwo = obsVarIndexRankTwo 



        do iz = 1 , domain_mean%size_bottomToTop

            call interp2d( domain_mean%lon(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                           domain_mean%lat(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                           domain_mean%GMH_unstag(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1,iz) , &
                           2 , 2 , &
                           (/gpsro%obs(io)%lon/) , (/gpsro%obs(io)%lat/) , obsZBuffer(iz:iz) , 1, &
                           1 , invalidValue )

            call interp2d( domain_mean%lon(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                           domain_mean%lat(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                           domain_mean%pressure(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1,iz) , &
                           2 , 2 , &
                           (/gpsro%obs(io)%lon/) , (/gpsro%obs(io)%lat/) , obsVarBuffer(iz:iz) , 1, &
                           1 , invalidValue )

        enddo

        select case ( trim(adjustl(gpsro%obs(io)%zName)) )
        case ( 'GMH' )
            call interp1d( obsZBuffer(:) , obsVarBuffer(:) , &
                           domain_mean%size_bottomToTop , &
                           (/gpsro%obs(io)%z/) , dummy_output(1:1) , &
                           1 , &
                           2 , .false. , invalidValue )
            gpsro%obs(io)%zName = 'P         '
            gpsro%obs(io)%z     = dummy_output(1)
        end select

    endif

enddo
#ifndef PGI
!$omp end parallel do
#endif



!================================================
return
end subroutine mapGPSROFromGMHToP

