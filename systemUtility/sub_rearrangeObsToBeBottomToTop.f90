
subroutine rearrangeObsToBeBottomToTop(obs)

use derivedType
use basicUtility

implicit none
type(obsParent),intent(inout) :: obs

real(kind=8),allocatable,dimension(:) :: dummy_z
integer,allocatable,dimension(:) :: indexOfObs

integer io
!================================================

if ( obs%obsNum .gt. 0 ) then

    allocate( dummy_z( obs%obsNum ) )
    allocate( indexOfObs( obs%obsNum ) )

    dummy_z(:) = obs%obs(:)%z
    forall(io=1:size(indexOfObs)) indexOfObs(io) = io


    call quickSortWithIndex( dummy_z(:) , size(dummy_z(:)) , indexOfObs(:) )

    obs%obs(:) = obs%obs( indexOfObs(obs%obsNum:1:-1) )  ! reverse for pressure coordinate

endif

!================================================
return
end subroutine rearrangeObsToBeBottomToTop
