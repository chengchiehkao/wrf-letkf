
!include 'sub_locateAsIndex2d.f90'
!include 'sub_interp2d.f90'

subroutine build_pressureU(domain,domainSize)

use derivedType
use basicUtility

implicit none

integer,intent(in)             :: domainSize
type(domainInfo),intent(inout) :: domain(domainSize)

integer targetRankOneIndex,targetRankTwoIndex
real(kind=8),dimension(1) :: dummy_pressureU
real(kind=8),parameter :: invalidValue = -888888.d0
integer,parameter :: boundaryOffset = 1  ! Assumed outmost grids are outside mass grid space.

integer id,iwe,isn,iz  ! loop conuter
!================================================


!
!  Assume horizontal coordinates of domain(1) & others are the same.
!
#ifndef PGI
!$omp parallel do default(private) shared(domain,domainSize) schedule(dynamic,5)
#endif
do isn = 1+boundaryOffset,domain(1)%size_southToNorth-boundaryOffset
do iwe = 1+boundaryOffset,domain(1)%size_westToEast_stag-boundaryOffset

    call locateAsIndex2d( domain(1)%lon(:,:) , domain(1)%lat(:,:) , domain(1)%size_westToEast , domain(1)%size_southToNorth , &
                          domain(1)%lon_u(iwe,isn) , domain(1)%lat_u(iwe,isn) , &
                          targetRankOneIndex , targetRankTwoIndex  )
    if ( targetRankOneIndex .eq. 0 .or. targetRankTwoIndex .eq. 0 ) then
       print*,'Out of boundary when locating (iwe,isn)=',iwe,isn
       stop
    endif

    do id = 1,domainSize

        do iz = 1,domain(id)%size_bottomToTop
            call interp2d( domain(id)%lon(targetRankOneIndex:targetRankOneIndex+1,targetRankTwoIndex:targetRankTwoIndex+1) , &
                           domain(id)%lat(targetRankOneIndex:targetRankOneIndex+1,targetRankTwoIndex:targetRankTwoIndex+1) , &
                           domain(id)%pressure(targetRankOneIndex:targetRankOneIndex+1,targetRankTwoIndex:targetRankTwoIndex+1,iz) , &
                           2 , 2 , &
                           (/domain(id)%lon_u(iwe,isn)/) , (/domain(id)%lat_u(iwe,isn)/) , dummy_pressureU(:) , &
                           1 , &
                           1 , invalidValue )
            domain(id)%pressure_u(iwe,isn,iz) = dummy_pressureU(1)
        enddo

    enddo

enddo
enddo
#ifndef PGI
!$omp end parallel do
#endif


!================================================
return
stop
end subroutine build_pressureU
