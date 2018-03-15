
subroutine mergeObs(allObs,obsToBeMerged)

use derivedType

implicit none
type(obsParent),intent(inout) :: allObs
type(obsParent),intent(inout) :: obsToBeMerged

type(obsChild),allocatable :: buffer(:)
integer bufferSize

integer io,ib  ! loop counter
!================================================

bufferSize = allObs%obsNum + obsToBeMerged%obsNum
allocate( buffer(bufferSize) )

do io = 1,allObs%obsNum
    buffer(io)%available  = allObs%obs(io)%available
    buffer(io)%lon        = allObs%obs(io)%lon
    buffer(io)%lat        = allObs%obs(io)%lat
    buffer(io)%z          = allObs%obs(io)%z
    buffer(io)%value      = allObs%obs(io)%value
    buffer(io)%error      = allObs%obs(io)%error
    buffer(io)%innov      = allObs%obs(io)%innov
    buffer(io)%rfict      = allObs%obs(io)%rfict
    buffer(io)%type       = allObs%obs(io)%type
    buffer(io)%instrument = allObs%obs(io)%instrument
    buffer(io)%zName      = allObs%obs(io)%zName
    buffer(io)%varName    = allObs%obs(io)%varName
    if ( associated(allObs%obs(io)%background) ) then
        allocate( buffer(io)%background( size(allObs%obs(io)%background(:)) ) )
        buffer(io)%background(:) = allObs%obs(io)%background(:)
        deallocate(allObs%obs(io)%background)
    endif
enddo

do io = 1,obsToBeMerged%obsNum
    ib = io + allObs%obsNum
    buffer(ib)%available  = obsToBeMerged%obs(io)%available
    buffer(ib)%lon        = obsToBeMerged%obs(io)%lon
    buffer(ib)%lat        = obsToBeMerged%obs(io)%lat
    buffer(ib)%z          = obsToBeMerged%obs(io)%z
    buffer(ib)%value      = obsToBeMerged%obs(io)%value
    buffer(ib)%error      = obsToBeMerged%obs(io)%error
    buffer(ib)%innov      = obsToBeMerged%obs(io)%innov
    buffer(ib)%rfict      = obsToBeMerged%obs(io)%rfict
    buffer(ib)%type       = obsToBeMerged%obs(io)%type
    buffer(ib)%instrument = obsToBeMerged%obs(io)%instrument
    buffer(ib)%zName      = obsToBeMerged%obs(io)%zName
    buffer(ib)%varName    = obsToBeMerged%obs(io)%varName
    if ( associated(obsToBeMerged%obs(io)%background) ) then
        allocate( buffer(ib)%background( size(obsToBeMerged%obs(io)%background(:)) ) )
        buffer(ib)%background(:) = obsToBeMerged%obs(io)%background(:)
    endif
enddo

if ( associated(allObs%obs) ) then
    deallocate(allObs%obs)
endif
allObs%obsNum = 0
allObs%obsNum = bufferSize
allocate( allObs%obs(allObs%obsNum) )

do io = 1,bufferSize
    allObs%obs(io)%available  = buffer(io)%available
    allObs%obs(io)%lon        = buffer(io)%lon
    allObs%obs(io)%lat        = buffer(io)%lat
    allObs%obs(io)%z          = buffer(io)%z
    allObs%obs(io)%value      = buffer(io)%value
    allObs%obs(io)%error      = buffer(io)%error
    allObs%obs(io)%innov      = buffer(io)%innov
    allObs%obs(io)%rfict      = buffer(io)%rfict
    allObs%obs(io)%type       = buffer(io)%type
    allObs%obs(io)%instrument = buffer(io)%instrument
    allObs%obs(io)%zName      = buffer(io)%zName
    allObs%obs(io)%varName    = buffer(io)%varName
    if ( associated(buffer(io)%background) ) then
        allocate( allObs%obs(io)%background( size(buffer(io)%background(:)) ) )
        allObs%obs(io)%background(:) = buffer(io)%background(:)
        deallocate( buffer(io)%background )
    endif
enddo

deallocate( buffer )

!================================================
return
end subroutine mergeObs
