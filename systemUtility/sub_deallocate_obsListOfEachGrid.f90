

subroutine deallocate_obsListOfEachGrid(obsListOfEachGrid)

use derivedType

implicit none
type(integerVector),pointer,dimension(:,:,:) :: obsListOfEachGrid

integer rank1UpperBound , rank2UpperBound , rank3UpperBound
integer iwe,isn,iz  ! loop counter
!================================================

if ( .not. associated(obsListOfEachGrid) ) then
    print*,'obsListOfEachGrid is disassociated before deallocated, deallocation failed.'
    return
endif

rank1UpperBound = ubound( obsListOfEachGrid , 1 )
rank2UpperBound = ubound( obsListOfEachGrid , 2 )
rank3UpperBound = ubound( obsListOfEachGrid , 3 )

do iz = 1,rank3UpperBound
do isn = 1,rank2UpperBound
do iwe = 1,rank1UpperBound
    if ( obsListOfEachGrid(iwe,isn,iz)%vectorSize .ne. 0 ) then
        deallocate( obsListOfEachGrid(iwe,isn,iz)%vector )
    endif
enddo
enddo
enddo


deallocate( obsListOfEachGrid )

!================================================
return
stop
end subroutine deallocate_obsListOfEachGrid 
