
!include 'sub_locateAsIndex2d.f90'
!include 'sub_interp1d.f90'
!include 'sub_interp2d.f90'

subroutine convertBackgroundToAMV(background,ensembleSize,domain,amv,systemParameters,domainID)

use derivedType
use basicUtility

implicit none

integer,intent(in)              :: ensembleSize
type(backgroundInfo),intent(in) :: background(ensembleSize)
type(domainInfo),intent(in)     :: domain(ensembleSize)
type(obsParent),intent(inout)   :: amv
type(systemParameter),intent(in) :: systemParameters
integer,intent(in)              :: domainID

integer obsVarIndexRankOne , obsVarIndexRankTwo
integer obsZIndexRankOne   , obsZIndexRankTwo

real(kind=8),dimension(domain(1)%size_bottomToTop) :: obsVarBuffer , obsZBuffer  ! SHALL AWARE OF DIFFERENCES BETWEEN P & GPH
real(kind=8) :: u_rotated , v_rotated
real(kind=8),dimension(1) :: sinalpha_interpolated

real(kind=8),parameter :: invalidValue = -9.d6

integer io,iens,iz  ! loop counter
!================================================


if ( systemParameters%rotateUAndVOfObsBasedOnWRFMapProjection ) then

    do io = 1 , amv%obsNum

        if ( amv%obs(io)%available .and. count(amv%obs(io)%insideHorizontalDomain(:)).eq.domainID ) then

            call locateAsIndex2d( domain(1)%lon(:,:)        , domain(1)%lat(:,:) , &
                                  domain(1)%size_westToEast , domain(1)%size_southToNorth , &
                                  amv%obs(io)%lon , amv%obs(io)%lat , &
                                  obsZIndexRankOne , obsZIndexRankTwo )

            select case ( trim(adjustl(amv%obs(io)%varName)) )
            case ( 'V' )
                call interp2d( domain(1)%lon(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                               domain(1)%lat(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                               domain(1)%sinalpha(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                               2 , 2 , &
                               (/amv%obs(io)%lon/) , (/amv%obs(io)%lat/) , sinalpha_interpolated(1:1) , 1, &
                               1 , invalidValue )

                call rotateUAndVofObsBasedOnSinAlpha( amv%obs(io-1)%value , amv%obs(io)%value , sinalpha_interpolated(1) , u_rotated , v_rotated )
                amv%obs(io-1)%value = u_rotated
                amv%obs(io  )%value = v_rotated
            end select

        endif

    enddo

endif


#ifndef PGI
!$omp parallel do default(none) &
!$omp private(io,iens,iz,obsVarIndexRankOne,obsVarIndexRankTwo,obsZIndexRankOne,obsZIndexRankTwo,obsVarBuffer,obsZBuffer) &
!$omp shared(amv,domain,background,ensembleSize,domainID) &
!$omp schedule(dynamic,100)
#endif
do io = 1 , amv%obsNum


    if ( amv%obs(io)%available .and. count(amv%obs(io)%insideHorizontalDomain(:)).eq.domainID ) then  ! SHALL AWARE OF PSFC

        allocate( amv%obs(io)%background(ensembleSize) )
        amv%obs(io)%background(:) = 0.d0

        select case ( trim(adjustl(amv%obs(io)%varName)) )
        case ( 'U' )
            call locateAsIndex2d( domain(1)%lon_u(:,:)        , domain(1)%lat_u(:,:) , &
                                  domain(1)%size_westToEast_stag , domain(1)%size_southToNorth , &
                                  amv%obs(io)%lon , amv%obs(io)%lat , &
                                  obsVarIndexRankOne , obsVarIndexRankTwo )
            if ( obsVarIndexRankOne.eq.0 .or. obsVarIndexRankTwo.eq.0 ) then
                deallocate( amv%obs(io)%background )
                amv%obs(io)%available = .false.
                cycle
            endif
            call locateAsIndex2d( domain(1)%lon(:,:)        , domain(1)%lat(:,:) , &
                                  domain(1)%size_westToEast , domain(1)%size_southToNorth , &
                                  amv%obs(io)%lon , amv%obs(io)%lat , &
                                  obsZIndexRankOne , obsZIndexRankTwo )
        case ( 'V' )
            call locateAsIndex2d( domain(1)%lon_v(:,:)        , domain(1)%lat_v(:,:) , &
                                  domain(1)%size_westToEast , domain(1)%size_southToNorth_stag , &
                                  amv%obs(io)%lon , amv%obs(io)%lat , &
                                  obsVarIndexRankOne , obsVarIndexRankTwo )
            if ( obsVarIndexRankOne.eq.0 .or. obsVarIndexRankTwo.eq.0 ) then
                deallocate( amv%obs(io)%background )
                amv%obs(io)%available = .false.
                cycle
            endif
            call locateAsIndex2d( domain(1)%lon(:,:)        , domain(1)%lat(:,:) , &
                                  domain(1)%size_westToEast , domain(1)%size_southToNorth , &
                                  amv%obs(io)%lon , amv%obs(io)%lat , &
                                  obsZIndexRankOne , obsZIndexRankTwo )
        end select


        do iens = 1 , ensembleSize
            do iz = 1 , domain(1)%size_bottomToTop  ! if obs%zName is P; shall be domain(1)%size_bottomToTop_stag for GPH

                select case ( trim(adjustl(amv%obs(io)%varName)) )
                case ( 'U' )
                    call interp2d( domain(1)%lon_u(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                                   domain(1)%lat_u(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                                   background(iens)%u(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1,iz) , &
                                   2 , 2 , &
                                   (/amv%obs(io)%lon/) , (/amv%obs(io)%lat/) , obsVarBuffer(iz:iz) , 1, &
                                   1 , invalidValue )
                    call interp2d( domain(1)%lon(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                                   domain(1)%lat(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                                   domain(iens)%pressure(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1,iz) , &
                                   2 , 2 , &
                                   (/amv%obs(io)%lon/) , (/amv%obs(io)%lat/) , obsZBuffer(iz:iz) , 1, &
                                   1 , invalidValue )
                case ( 'V' )
                    call interp2d( domain(1)%lon_v(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                                   domain(1)%lat_v(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                                   background(iens)%v(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1,iz) , &
                                   2 , 2 , &
                                   (/amv%obs(io)%lon/) , (/amv%obs(io)%lat/) , obsVarBuffer(iz:iz) , 1, &
                                   1 , invalidValue )
                    call interp2d( domain(1)%lon(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                                   domain(1)%lat(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                                   domain(iens)%pressure(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1,iz) , &
                                   2 , 2 , &
                                   (/amv%obs(io)%lon/) , (/amv%obs(io)%lat/) , obsZBuffer(iz:iz) , 1, &
                                   1 , invalidValue )
                end select

            enddo

            select case ( trim(adjustl(amv%obs(io)%zName)) )
            case ( 'P' )
                call interp1d( dlog(obsZBuffer(domain(1)%size_bottomToTop:1:-1)) , obsVarBuffer(domain(1)%size_bottomToTop:1:-1) , &
                               domain(1)%size_bottomToTop , &
                               dlog((/amv%obs(io)%z/)) , amv%obs(io)%background(iens:iens) , &
                               1 , &
                               2 , .false. , invalidValue )
            end select

        enddo

        amv%obs(io)%innov = amv%obs(io)%value - sum(amv%obs(io)%background(:))/real(ensembleSize,8)

    endif

enddo
#ifndef PGI
!$omp end parallel do
#endif



!================================================
return
stop
end subroutine convertBackgroundToAMV

