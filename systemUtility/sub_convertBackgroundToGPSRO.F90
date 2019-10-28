
!include 'sub_locateAsIndex2d.f90'
!include 'sub_interp1d.f90'
!include 'sub_interp2d.f90'

subroutine convertBackgroundToGPSRO(background,ensembleSize,domain,gpsro)

use derivedType
use basicUtility

implicit none

integer,intent(in)              :: ensembleSize
type(backgroundInfo),intent(in) :: background(ensembleSize)
type(domainInfo),intent(in)     :: domain(ensembleSize)
type(obsParent),intent(inout)   :: gpsro

integer obsVarIndexRankOne , obsVarIndexRankTwo
integer obsZIndexRankOne   , obsZIndexRankTwo

real(kind=8),dimension(domain(1)%size_bottomToTop) :: obsVarBuffer , obsZBuffer  ! SHALL AWARE OF DIFFERENCES BETWEEN P & GPH

real(kind=8),parameter :: invalidValue = -9.d6

integer io,iens,iz  ! loop counter
!================================================

#ifndef PGI
!$omp parallel do default(none) &
!$omp private(io,iens,iz,obsVarIndexRankOne,obsVarIndexRankTwo,obsZIndexRankOne,obsZIndexRankTwo,obsVarBuffer,obsZBuffer) &
!$omp shared(gpsro,domain,background,ensembleSize) &
!$omp schedule(dynamic,100)
#endif
do io = 1 , gpsro%obsNum


    if ( gpsro%obs(io)%available ) then  ! SHALL AWARE OF PSFC

        allocate( gpsro%obs(io)%background(ensembleSize) )
        gpsro%obs(io)%background(:) = 0.d0

        select case ( trim(adjustl(gpsro%obs(io)%varName)) )
        case ( 'REF' )
            call locateAsIndex2d( domain(1)%lon(:,:)        , domain(1)%lat(:,:) , &
                                  domain(1)%size_westToEast , domain(1)%size_southToNorth , &
                                  gpsro%obs(io)%lon , gpsro%obs(io)%lat , &
                                  obsVarIndexRankOne , obsVarIndexRankTwo )
            if ( obsVarIndexRankOne.eq.0 .or. obsVarIndexRankTwo.eq.0 ) then
                deallocate( gpsro%obs(io)%background )
                gpsro%obs(io)%available = .false.
                cycle
            endif
            obsZIndexRankOne = obsVarIndexRankOne
            obsZIndexRankTwo = obsVarIndexRankTwo 
        end select


        do iens = 1 , ensembleSize
            do iz = 1 , domain(1)%size_bottomToTop  ! if obs%zName is P or GMH_unstag; shall be domain(1)%size_bottomToTop_stag for GPH

                select case ( trim(adjustl(gpsro%obs(io)%varName)) )
                case ( 'REF' )
                    call interp2d( domain(1)%lon(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                                   domain(1)%lat(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                                   background(iens)%refractivity(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1,iz) , &
                                   2 , 2 , &
                                   (/gpsro%obs(io)%lon/) , (/gpsro%obs(io)%lat/) , obsVarBuffer(iz:iz) , 1, &
                                   1 , invalidValue )
                    call interp2d( domain(1)%lon(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                                   domain(1)%lat(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                                   domain(iens)%GMH_unstag(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1,iz) , &
                                   2 , 2 , &
                                   (/gpsro%obs(io)%lon/) , (/gpsro%obs(io)%lat/) , obsZBuffer(iz:iz) , 1, &
                                   1 , invalidValue )
                end select

            enddo

            select case ( trim(adjustl(gpsro%obs(io)%zName)) )
            case ( 'GMH' )
                call interp1d( obsZBuffer(:) , obsVarBuffer(:) , &
                               domain(1)%size_bottomToTop , &
                               (/gpsro%obs(io)%z/) , gpsro%obs(io)%background(iens:iens) , &
                               1 , &
                               2 , .false. , invalidValue )
            end select

        enddo

        gpsro%obs(io)%innov = gpsro%obs(io)%value - sum(gpsro%obs(io)%background(:))/real(ensembleSize,8)

    endif

enddo
#ifndef PGI
!$omp end parallel do
#endif



!================================================
return
stop
end subroutine convertBackgroundToGPSRO

