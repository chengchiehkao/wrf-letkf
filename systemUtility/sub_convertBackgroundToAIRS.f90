
!include 'sub_locateAsIndex2d.f90'
!include 'sub_interp1d.f90'
!include 'sub_interp2d.f90'

subroutine convertBackgroundToAIRS(background,ensembleSize,domain,airs)

use derivedType
use basicUtility

implicit none

integer,intent(in)              :: ensembleSize
type(backgroundInfo),intent(in) :: background(ensembleSize)
type(domainInfo),intent(in)     :: domain(ensembleSize)
type(obsParent),intent(inout)   :: airs

integer obsVarIndexRankOne , obsVarIndexRankTwo
integer obsZIndexRankOne   , obsZIndexRankTwo

real(kind=8),dimension(domain(1)%size_bottomToTop) :: obsVarBuffer , obsZBuffer  ! SHALL AWARE OF DIFFERENCES BETWEEN P & GPH

real(kind=8),parameter :: invalidValue = -9.d6

integer io,iens,iz  ! loop counter
!================================================

!$omp parallel do default(none) &
!$omp private(io,iens,iz,obsVarIndexRankOne,obsVarIndexRankTwo,obsZIndexRankOne,obsZIndexRankTwo,obsVarBuffer,obsZBuffer) &
!$omp shared(airs,domain,background,ensembleSize) &
!$omp schedule(dynamic,100)
do io = 1 , airs%obsNum


    if ( airs%obs(io)%available ) then  ! SHALL AWARE OF PSFC

        allocate( airs%obs(io)%background(ensembleSize) )
        airs%obs(io)%background(:) = 0.d0  ! MAY REMOVED AFTER BUG SOLVED

        if ( .not. associated( airs%obs(io)%background ) ) then
            print*,'found abnormal alllocation status.'
        endif

        !select case ( trim(adjustl(airs%obs(io)%varName)) )
        select case ( airs%obs(io)%varName(1:1) )
        case ( 'U' )
            call locateAsIndex2d( domain(1)%lon_u(:,:)        , domain(1)%lat_u(:,:) , &
                                  domain(1)%size_westToEast_stag , domain(1)%size_southToNorth , &
                                  airs%obs(io)%lon , airs%obs(io)%lat , &
                                  obsVarIndexRankOne , obsVarIndexRankTwo )
            if ( obsVarIndexRankOne.eq.0 .or. obsVarIndexRankTwo.eq.0 ) then
                deallocate( airs%obs(io)%background )
                airs%obs(io)%available = .false.
                cycle
            endif
            call locateAsIndex2d( domain(1)%lon(:,:)        , domain(1)%lat(:,:) , &
                                  domain(1)%size_westToEast , domain(1)%size_southToNorth , &
                                  airs%obs(io)%lon , airs%obs(io)%lat , &
                                  obsZIndexRankOne , obsZIndexRankTwo )
        case ( 'V' )
            call locateAsIndex2d( domain(1)%lon_v(:,:)        , domain(1)%lat_v(:,:) , &
                                  domain(1)%size_westToEast , domain(1)%size_southToNorth_stag , &
                                  airs%obs(io)%lon , airs%obs(io)%lat , &
                                  obsVarIndexRankOne , obsVarIndexRankTwo )
            if ( obsVarIndexRankOne.eq.0 .or. obsVarIndexRankTwo.eq.0 ) then
                deallocate( airs%obs(io)%background )
                airs%obs(io)%available = .false.
                cycle
            endif
            call locateAsIndex2d( domain(1)%lon(:,:)        , domain(1)%lat(:,:) , &
                                  domain(1)%size_westToEast , domain(1)%size_southToNorth , &
                                  airs%obs(io)%lon , airs%obs(io)%lat , &
                                  obsZIndexRankOne , obsZIndexRankTwo )
        case ( 'T' )
            call locateAsIndex2d( domain(1)%lon(:,:)        , domain(1)%lat(:,:) , &
                                  domain(1)%size_westToEast , domain(1)%size_southToNorth , &
                                  airs%obs(io)%lon , airs%obs(io)%lat , &
                                  obsVarIndexRankOne , obsVarIndexRankTwo )
            if ( obsVarIndexRankOne.eq.0 .or. obsVarIndexRankTwo.eq.0 ) then
                deallocate( airs%obs(io)%background )
                airs%obs(io)%available = .false.
                cycle
            endif
            obsZIndexRankOne = obsVarIndexRankOne
            obsZIndexRankTwo = obsVarIndexRankTwo
        case ( 'Q' )
        !case ( 'QVAPOR' )
            call locateAsIndex2d( domain(1)%lon(:,:)        , domain(1)%lat(:,:) , &
                                  domain(1)%size_westToEast , domain(1)%size_southToNorth , &
                                  airs%obs(io)%lon , airs%obs(io)%lat , &
                                  obsVarIndexRankOne , obsVarIndexRankTwo )
            if ( obsVarIndexRankOne.eq.0 .or. obsVarIndexRankTwo.eq.0 ) then
                deallocate( airs%obs(io)%background )
                airs%obs(io)%available = .false.
               cycle
            endif
            obsZIndexRankOne = obsVarIndexRankOne
            obsZIndexRankTwo = obsVarIndexRankTwo 
        end select


        do iens = 1 , ensembleSize
            do iz = 1 , domain(1)%size_bottomToTop  ! if obs%zName is P; shall be domain(1)%size_bottomToTop_stag for GPH

                select case ( airs%obs(io)%varName(1:1) )
                !select case ( trim(adjustl(airs%obs(io)%varName)) )
                case ( 'U' )
                    call interp2d( domain(1)%lon_u(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                                   domain(1)%lat_u(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                                   background(iens)%u(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1,iz) , &
                                   2 , 2 , &
                                   (/airs%obs(io)%lon/) , (/airs%obs(io)%lat/) , obsVarBuffer(iz:iz) , 1, &
                                   1 , invalidValue )
                    call interp2d( domain(1)%lon(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                                   domain(1)%lat(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                                   domain(iens)%pressure(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1,iz) , &
                                   2 , 2 , &
                                   (/airs%obs(io)%lon/) , (/airs%obs(io)%lat/) , obsZBuffer(iz:iz) , 1, &
                                   1 , invalidValue )
                case ( 'V' )
                    call interp2d( domain(1)%lon_v(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                                   domain(1)%lat_v(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                                   background(iens)%v(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1,iz) , &
                                   2 , 2 , &
                                   (/airs%obs(io)%lon/) , (/airs%obs(io)%lat/) , obsVarBuffer(iz:iz) , 1, &
                                   1 , invalidValue )
                    call interp2d( domain(1)%lon(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                                   domain(1)%lat(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                                   domain(iens)%pressure(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1,iz) , &
                                   2 , 2 , &
                                   (/airs%obs(io)%lon/) , (/airs%obs(io)%lat/) , obsZBuffer(iz:iz) , 1, &
                                   1 , invalidValue )
                case ( 'T' )
                    call interp2d( domain(1)%lon(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                                   domain(1)%lat(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                                   background(iens)%normalT(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1,iz) , &
                                   2 , 2 , &
                                   (/airs%obs(io)%lon/) , (/airs%obs(io)%lat/) , obsVarBuffer(iz:iz) , 1, &
                                   1 , invalidValue )
                    call interp2d( domain(1)%lon(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                                   domain(1)%lat(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                                   domain(iens)%pressure(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1,iz) , &
                                   2 , 2 , &
                                   (/airs%obs(io)%lon/) , (/airs%obs(io)%lat/) , obsZBuffer(iz:iz) , 1, &
                                   1 , invalidValue )
                case ( 'Q' )
                !case ( 'QVAPOR' )
                    call interp2d( domain(1)%lon(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                                   domain(1)%lat(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                                   background(iens)%qvapor(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1,iz) , &
                                   2 , 2 , &
                                   (/airs%obs(io)%lon/) , (/airs%obs(io)%lat/) , obsVarBuffer(iz:iz) , 1, &
                                   1 , invalidValue )
                    call interp2d( domain(1)%lon(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                                   domain(1)%lat(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                                   domain(iens)%pressure(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1,iz) , &
                                   2 , 2 , &
                                   (/airs%obs(io)%lon/) , (/airs%obs(io)%lat/) , obsZBuffer(iz:iz) , 1, &
                                   1 , invalidValue )
                end select

            enddo

            !select case ( trim(adjustl(airs%obs(io)%zName)) )
            select case ( airs%obs(io)%zName(1:1) )
            case ( 'P' )
                call interp1d( dlog(obsZBuffer(domain(1)%size_bottomToTop:1:-1)) , obsVarBuffer(domain(1)%size_bottomToTop:1:-1) , &
                               domain(1)%size_bottomToTop , &
                               dlog((/airs%obs(io)%z/)) , airs%obs(io)%background(iens:iens) , &
                               1 , &
                               2 , .false. , invalidValue )
            end select

        enddo

        airs%obs(io)%innov = airs%obs(io)%value - sum(airs%obs(io)%background(:))/real(ensembleSize,8)

    endif

enddo
!$omp end parallel do



!================================================
return
stop
end subroutine convertBackgroundToAIRS

