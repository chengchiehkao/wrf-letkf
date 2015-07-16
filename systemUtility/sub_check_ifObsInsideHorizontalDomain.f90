
!include 'func_inPolygon.f90'

subroutine check_ifObsInsideHorizontalDomain(domain,obs)

use derivedType
use basicUtility
implicit none
type(domainInfo),intent(in) :: domain
type(obsParent),intent(inout) :: obs

integer io

!logical,external :: inPolygon
!================================================

!domain%size_westToEast = size(domain%lon(:,:),1)
!domain%size_southToNorth = size(domain%lon(:,:),2)

do io=1,obs%obsNum

    obs%obs(io)%available = inPolygon( (/domain%lon(1:domain%size_westToEast-1,1), &
                                         domain%lon(domain%size_westToEast,1:domain%size_southToNorth-1), &
                                         domain%lon(domain%size_westToEast:2:-1,domain%size_southToNorth), &
                                         domain%lon(1,domain%size_southToNorth:2:-1)/) , &
                                       (/domain%lat(1:domain%size_westToEast-1,1), &
                                         domain%lat(domain%size_westToEast,1:domain%size_southToNorth-1 ), &
                                         domain%lat(domain%size_westToEast:2:-1,domain%size_southToNorth), &
                                         domain%lat(1,domain%size_southToNorth:2:-1)/) , &
                                         2*(domain%size_westToEast+domain%size_southToNorth)-4, &
                                         obs%obs(io)%lon,obs%obs(io)%lat )

enddo

!================================================
return
stop
end subroutine check_ifObsInsideHorizontalDomain
