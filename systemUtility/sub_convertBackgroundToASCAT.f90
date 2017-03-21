
!include 'sub_locateAsIndex2d.f90'
!include 'sub_interp1d.f90'
!include 'sub_interp2d.f90'

subroutine convertBackgroundToASCAT(background,ensembleSize,domain,domain_mean,ascat)

use derivedType
use basicUtility

implicit none

integer,intent(in)              :: ensembleSize
type(backgroundInfo),intent(in) :: background(ensembleSize)
type(domainInfo),intent(in)     :: domain(ensembleSize)
type(domainInfo),intent(in)     :: domain_mean
type(obsParent),intent(inout)   :: ascat

integer obsVarIndexRankOne , obsVarIndexRankTwo
integer obsZIndexRankOne   , obsZIndexRankTwo

real(kind=8),dimension(1) :: obsVarBuffer , obsZBuffer  ! SHALL AWARE OF DIFFERENCES BETWEEN P & GPH

real(kind=8),parameter :: invalidValue = -9.d6

integer io,iens,iz  ! loop counter
!================================================

!$omp parallel do default(none) &
!$omp private(io,iens,iz,obsVarIndexRankOne,obsVarIndexRankTwo,obsZIndexRankOne,obsZIndexRankTwo,obsVarBuffer,obsZBuffer) &
!$omp shared(ascat,domain,domain_mean,background,ensembleSize) &
!$omp schedule(dynamic,100)
do io = 1 , ascat%obsNum


    if ( ascat%obs(io)%available ) then

        allocate( ascat%obs(io)%background(ensembleSize) )
        ascat%obs(io)%background(:) = 0.d0  ! MAY REMOVED AFTER BUG SOLVED

        if ( .not. associated( ascat%obs(io)%background ) ) then
            print*,'found abnormal alllocation status.'
        endif

        !select case ( trim(adjustl(ascat%obs(io)%varName)) )
        select case ( ascat%obs(io)%varName(1:3) )
        case ( 'U10' )
            call locateAsIndex2d( domain(1)%lon(:,:)        , domain(1)%lat(:,:) , &
                                  domain(1)%size_westToEast , domain(1)%size_southToNorth , &
                                  ascat%obs(io)%lon , ascat%obs(io)%lat , &
                                  obsVarIndexRankOne , obsVarIndexRankTwo )
            if ( obsVarIndexRankOne.eq.0 .or. obsVarIndexRankTwo.eq.0 ) then
                deallocate( ascat%obs(io)%background )
                ascat%obs(io)%available = .false.
                cycle
            endif
            obsZIndexRankOne = obsVarIndexRankOne
            obsZIndexRankTwo = obsVarIndexRankTwo
        case ( 'V10' )
            call locateAsIndex2d( domain(1)%lon(:,:)        , domain(1)%lat(:,:) , &
                                  domain(1)%size_westToEast , domain(1)%size_southToNorth , &
                                  ascat%obs(io)%lon , ascat%obs(io)%lat , &
                                  obsVarIndexRankOne , obsVarIndexRankTwo )
            if ( obsVarIndexRankOne.eq.0 .or. obsVarIndexRankTwo.eq.0 ) then
                deallocate( ascat%obs(io)%background )
                ascat%obs(io)%available = .false.
                cycle
            endif
            obsZIndexRankOne = obsVarIndexRankOne
            obsZIndexRankTwo = obsVarIndexRankTwo
        end select


        do iens = 1 , ensembleSize

            select case ( ascat%obs(io)%varName(1:3) )
            !select case ( trim(adjustl(ascat%obs(io)%varName)) )
            case ( 'U10' )
                call interp2d( domain(1)%lon(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                               domain(1)%lat(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                               background(iens)%u10(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                               2 , 2 , &
                               (/ascat%obs(io)%lon/) , (/ascat%obs(io)%lat/) , obsVarBuffer(1:1) , 1, &
                               1 , invalidValue )
                call interp2d( domain_mean%lon(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                               domain_mean%lat(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                               domain_mean%pressure_w(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1,1) , &
                               2 , 2 , &
                               (/ascat%obs(io)%lon/) , (/ascat%obs(io)%lat/) , obsZBuffer(1:1) , 1, &
                               1 , invalidValue )
            case ( 'V10' )
                call interp2d( domain(1)%lon(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                               domain(1)%lat(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                               background(iens)%v10(obsVarIndexRankOne:obsVarIndexRankOne+1,obsVarIndexRankTwo:obsVarIndexRankTwo+1) , &
                               2 , 2 , &
                               (/ascat%obs(io)%lon/) , (/ascat%obs(io)%lat/) , obsVarBuffer(1:1) , 1, &
                               1 , invalidValue )
                call interp2d( domain_mean%lon(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                               domain_mean%lat(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1) , &
                               domain_mean%pressure_w(obsZIndexRankOne:obsZIndexRankOne+1,obsZIndexRankTwo:obsZIndexRankTwo+1,1) , &
                               2 , 2 , &
                               (/ascat%obs(io)%lon/) , (/ascat%obs(io)%lat/) , obsZBuffer(1:1) , 1, &
                               1 , invalidValue )
            end select


            !select case ( trim(adjustl(ascat%obs(io)%zName)) )
            select case ( ascat%obs(io)%zName(1:1) )
            case ( 'P' )
                ascat%obs(io)%z = obsZBuffer(1) - 100.d0  ! Assumed 10-m high means 1-hPa less than PSFC of mean domain.
                ascat%obs(io)%background(iens) = obsVarBuffer(1)
            end select

        enddo

        ascat%obs(io)%innov = ascat%obs(io)%value - sum(ascat%obs(io)%background(:))/real(ensembleSize,8)

    endif

enddo
!$omp end parallel do



!================================================
return
stop
end subroutine convertBackgroundToASCAT

