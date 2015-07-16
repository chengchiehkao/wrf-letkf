
!  include 'mod_derivedType.f90'
!  include 'func_availableFileID.f90'

subroutine getSynop(synop,varList,varListSize)

use derivedType

implicit none
type(obsParent),intent(out) :: synop
integer,intent(in)                                   :: varListSize
character(len=10),dimension(varListSize),intent(in)  :: varList

real(kind=8) lon,lat
integer levNum,serialNum
real :: dummyArg1,dummyArg2

integer ioStatus,fileID
integer i  ! loop counter

integer,external :: availableFileID
!================================================

Synop%obsNum = 0

fileID = availableFileID()  ! get an available file id.

! First, find out the number of observation.
print*,'reading observation: synop'
open(fileID,file='input/obs_input.synop',status='old')
do 
    read(fileID,*,iostat=ioStatus) dummyArg1,dummyArg2 , levNum
    if ( ioStatus .lt. 0 ) then
        exit
    endif
    synop%obsNum = synop%obsNum + varListSize*levNum
    do i = 1,levNum
        read(fileID,*)  ! just read through lines which are not headers.
    enddo
enddo

! Then, allocate all necessary components.
allocate( synop%obs(synop%obsNum) )

! Initilize the components just allocated.
synop%obs(:)%lon   = 0.0d0
synop%obs(:)%lat   = 0.0d0
synop%obs(:)%z     = 0.0d0
synop%obs(:)%var   = 0.0d0
synop%obs(:)%error = 0.0d0
synop%obs(:)%instrument = 'SYNOP'
synop%obs(:)%zName      = 'TerrainHgt'
synop%obs(:)%varName    = ''
synop%obs(:)%available  = .true.

rewind(fileID)

! Finally, read all observations.
serialNum = 0
do 
    read(fileID,*,iostat=ioStatus) lon,lat,levNum
    if ( ioStatus .lt. 0 ) exit
    do i = 1,levNum
        serialNum = serialNum + 1
        read(fileID,*) synop%obs(serialNum:serialNum+varListSize-1)%var,synop%obs(serialNum)%z

        synop%obs(serialNum:serialNum+varListSize-1)%lon = lon
        synop%obs(serialNum:serialNum+varListSize-1)%lat = lat
        synop%obs(serialNum:serialNum+varListSize-1)%z = synop%obs(serialNum)%z
        synop%obs(serialNum:serialNum+varListSize-1)%varName = varList(1:varListSize)
        serialNum = serialNum + (varListSize-1)
    enddo
enddo

close(fileID)
!================================================
return
stop
end subroutine getSynop
